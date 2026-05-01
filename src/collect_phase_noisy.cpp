/**
 * @file collect_phase_noisy.cpp
 * @brief Step 4: iterative noisy-region MSA variant calling.
 *
 * Ports longcallD step 4 from `collect_var.c` / `align.c`:
 *   `sort_noisy_regs`, `collect_noisy_vars1`, and the outer while-loop in
 *   `collect_var_main` that gates k-means re-runs on `kCandGermlineVarCate`.
 *
 * Alignment-heavy functions (`collect_noisy_reg_aln_strs` and helpers) live in
 * `align.cpp`; this file owns the MSA-to-candidate conversion, read profile
 * update, old/new candidate merge, and iterative k-means trigger.
 */

#include "collect_phase_noisy.hpp"
#include "collect_phase.hpp" // assign_hap_based_on_germline_het_vars_kmeans, kCandGermlineVarCate
#include "collect_var.hpp"   // exact_comp_cand_var, exact_comp_var_site

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <numeric>

namespace pgphase_collect {

namespace {

constexpr uint8_t kGapBase = 5;
constexpr const char* kNt4Bases = "ACGTN";

char var_type_code_local(VariantType type) {
    switch (type) {
        case VariantType::Snp:
            return 'X';
        case VariantType::Insertion:
            return 'I';
        case VariantType::Deletion:
            return 'D';
    }
    return '?';
}

const char* category_name_local(VariantCategory category) {
    switch (category) {
        case VariantCategory::LowCoverage:         return "LOW_COV";
        case VariantCategory::LowAlleleFraction:   return "LOW_AF";
        case VariantCategory::StrandBias:          return "STRAND_BIAS";
        case VariantCategory::CleanHetSnp:         return "CLEAN_HET_SNP";
        case VariantCategory::CleanHetIndel:       return "CLEAN_HET_INDEL";
        case VariantCategory::CleanHom:            return "CLEAN_HOM";
        case VariantCategory::NoisyCandHet:        return "NOISY_CAND_HET";
        case VariantCategory::NoisyCandHom:        return "NOISY_CAND_HOM";
        case VariantCategory::NoisyResolved:       return "NOISY_RESOLVED";
        case VariantCategory::RepeatHetIndel:      return "REP_HET_INDEL";
        case VariantCategory::NonVariant:          return "NON_VAR";
    }
    return "UNKNOWN";
}

char nt4_to_base(uint8_t b) {
    return b < 4 ? kNt4Bases[b] : 'N';
}

uint8_t base_to_nt4(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
        case 'U':
        case 'u':
            return 3;
        default:
            return 4;
    }
}

int alt_len_for_key(const VariantKey& key) {
    return key.type == VariantType::Deletion ? 0 : static_cast<int>(key.alt.size());
}

void update_variant_depth_fields(CandidateVariant& var) {
    if (var.counts.alle_covs.empty())
        var.counts.alle_covs.assign(static_cast<size_t>(std::max(2, var.counts.n_uniq_alles)), 0);
    var.counts.n_uniq_alles = static_cast<int>(var.counts.alle_covs.size());
    var.counts.ref_cov = var.counts.alle_covs.empty() ? 0 : var.counts.alle_covs[0];
    var.counts.alt_cov = var.counts.alle_covs.size() > 1 ? var.counts.alle_covs[1] : 0;
    int alle_cov_sum = 0;
    for (int cov : var.counts.alle_covs) alle_cov_sum += cov;
    if (var.counts.total_cov <= 0) var.counts.total_cov = alle_cov_sum;
    var.counts.allele_fraction =
        var.counts.total_cov > 0 ? static_cast<double>(var.counts.alt_cov) / var.counts.total_cov : 0.0;
}

void set_noisy_category(CandidateVariant& var, VariantCategory category) {
    var.counts.category = category;
    var.lcd_var_i_to_cate = category_to_flag(category);
}

static void restore_stashed_initial_if_any(BamChunk& chunk, CandidateVariant& var) {
    auto& stash = chunk.erased_clean_signal_initial;
    for (size_t i = 0; i < stash.size(); ++i) {
        if (exact_comp_var_site(&stash[i].first, &var.key) != 0) continue;
        var.counts.candvarcate_initial = stash[i].second;
        if (i + 1 < stash.size()) stash[i] = std::move(stash.back());
        stash.pop_back();
        return;
    }
}

void update_read_var_profile_with_allele(int var_idx, int allele, int alt_qi, ReadVariantProfile& profile) {
    if (var_idx < 0) return;
    if (profile.start_var_idx < 0) {
        profile.start_var_idx = var_idx;
        profile.end_var_idx = var_idx;
        profile.alleles.assign(1, allele);
        profile.alt_qi.assign(1, alt_qi);
        return;
    }
    if (var_idx < profile.start_var_idx) {
        const int prefix = profile.start_var_idx - var_idx;
        profile.alleles.insert(profile.alleles.begin(), static_cast<size_t>(prefix), -1);
        profile.alt_qi.insert(profile.alt_qi.begin(), static_cast<size_t>(prefix), -1);
        profile.start_var_idx = var_idx;
    }
    if (var_idx > profile.end_var_idx) {
        const int gap = var_idx - profile.end_var_idx - 1;
        profile.alleles.insert(profile.alleles.end(), static_cast<size_t>(gap), -1);
        profile.alt_qi.insert(profile.alt_qi.end(), static_cast<size_t>(gap), -1);
        profile.end_var_idx = var_idx;
        profile.alleles.push_back(allele);
        profile.alt_qi.push_back(alt_qi);
        return;
    }
    const int offset = var_idx - profile.start_var_idx;
    profile.alleles[static_cast<size_t>(offset)] = allele;
    profile.alt_qi[static_cast<size_t>(offset)] = alt_qi;
}

std::vector<ReadVariantProfile> init_read_profiles(size_t n_reads) {
    std::vector<ReadVariantProfile> profiles(n_reads);
    for (size_t i = 0; i < n_reads; ++i)
        profiles[i].read_id = static_cast<int>(i);
    return profiles;
}

bool var_is_homopolymer_indel(const BamChunk& chunk,
                              hts_pos_t ref_pos,
                              VariantType type,
                              int ref_len,
                              const std::string& alt) {
    // Literal port of longcallD collect_var.c:1720 `var_is_homopolymer_indel`.
    // INS: ref loop compares raw FASTA bytes to nt4 alt (same as LCD); DEL: raw bytes throughout.
    if (type == VariantType::Snp) return false;
    if (type == VariantType::Insertion) {
        if (alt.empty()) return false;
        const uint8_t ins_base0 = base_to_nt4(alt[0]);
        for (size_t i = 1; i < alt.size(); ++i) {
            if (base_to_nt4(alt[i]) != ins_base0) return false;
        }
        for (int i = 0; i < 5; ++i) {
            const size_t idx = static_cast<size_t>(ref_pos - chunk.ref_beg + static_cast<hts_pos_t>(i));
            if (static_cast<uint8_t>(static_cast<unsigned char>(chunk.ref_seq[idx])) != ins_base0) return false;
        }
        return true;
    }
    const size_t idx0 = static_cast<size_t>(ref_pos - chunk.ref_beg);
    const uint8_t ref_base0 = static_cast<uint8_t>(static_cast<unsigned char>(chunk.ref_seq[idx0]));
    for (int i = 1; i < ref_len; ++i) {
        const size_t idx = idx0 + static_cast<size_t>(i);
        if (static_cast<uint8_t>(static_cast<unsigned char>(chunk.ref_seq[idx])) != ref_base0) return false;
    }
    for (int i = 0; i < 5; ++i) {
        const size_t idx = idx0 + static_cast<size_t>(i);
        if (static_cast<uint8_t>(static_cast<unsigned char>(chunk.ref_seq[idx])) != ref_base0) return false;
    }
    return true;
}

CandidateVariant make_noisy_candidate(const BamChunk& chunk,
                                      hts_pos_t ref_pos,
                                      VariantType type,
                                      int ref_len,
                                      uint8_t ref_base,
                                      std::string alt,
                                      uint8_t alt_ref_base,
                                      bool is_homopolymer_indel) {
    CandidateVariant var;
    var.key.tid = chunk.region.tid;
    var.key.pos = ref_pos;
    var.key.type = type;
    var.key.ref_len = ref_len;
    var.key.alt = std::move(alt);
    var.ref_base = type == VariantType::Snp && ref_base < 4 ? ref_base : 4;
    // INS/DEL: keep consensus anchor encoding (e.g. abPOA gap == 5) for longcallD VCF parity;
    // SNP path does not use alt_ref_base (longcallD fills separately).
    var.alt_ref_base = type == VariantType::Snp ? static_cast<uint8_t>(4) : alt_ref_base;
    var.is_homopolymer_indel = is_homopolymer_indel;
    var.counts.n_uniq_alles = 2;
    var.counts.alle_covs.assign(2, 0);
    var.counts.category = VariantCategory::NoisyCandHom;
    var.lcd_var_i_to_cate = category_to_flag(VariantCategory::NoisyCandHom);
    return var;
}

std::vector<CandidateVariant> make_cand_vars_from_baln0(const Options& opts,
                                                        const BamChunk& chunk,
                                                        hts_pos_t noisy_reg_beg,
                                                        const std::vector<uint8_t>& ref_msa_seq,
                                                        const std::vector<uint8_t>& cons_msa_seq,
                                                        bool no_end_var) {
    std::vector<CandidateVariant> vars;
    vars.reserve(ref_msa_seq.size());

    hts_pos_t ref_pos = noisy_reg_beg;
    int i = 0;
    const int msa_len = static_cast<int>(ref_msa_seq.size());
    while (i < msa_len) {
        const uint8_t ref_base = ref_msa_seq[static_cast<size_t>(i)];
        const uint8_t cons_base = cons_msa_seq[static_cast<size_t>(i)];
        if (ref_base == kGapBase && cons_base == kGapBase) {
            ++i;
            continue;
        }
        if (ref_base == cons_base) {
            ++i;
            ++ref_pos;
            continue;
        }
        if (ref_base != kGapBase && cons_base != kGapBase) {
            const bool next_ref_non_gap = (i + 1 == msa_len) || ref_msa_seq[static_cast<size_t>(i + 1)] != kGapBase;
            const bool next_cons_non_gap = (i + 1 == msa_len) || cons_msa_seq[static_cast<size_t>(i + 1)] != kGapBase;
            if (next_ref_non_gap && next_cons_non_gap) {
                vars.push_back(make_noisy_candidate(chunk, ref_pos, VariantType::Snp, 1,
                                                    ref_base, std::string(1, nt4_to_base(cons_base)),
                                                    4, false));
            }
            ++i;
            ++ref_pos;
            continue;
        }
        if (ref_base == kGapBase) {
            int gap_len = 1;
            while (i + gap_len < msa_len &&
                   ref_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase &&
                   cons_msa_seq[static_cast<size_t>(i + gap_len)] != kGapBase) {
                ++gap_len;
            }
            if (no_end_var && (i - 1 < 0 || i + gap_len >= msa_len ||
                               ref_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                               ref_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase ||
                               cons_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                               cons_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase)) {
                i += gap_len;
                continue;
            }
            std::string alt;
            alt.reserve(static_cast<size_t>(gap_len));
            for (int k = 0; k < gap_len; ++k)
                alt.push_back(nt4_to_base(cons_msa_seq[static_cast<size_t>(i + k)]));

            const bool hp = gap_len < opts.min_sv_len &&
                            var_is_homopolymer_indel(chunk, ref_pos, VariantType::Insertion, 0, alt);
            const uint8_t alt_ref_base = i - 1 >= 0 ? cons_msa_seq[static_cast<size_t>(i - 1)] : 4;
            vars.push_back(make_noisy_candidate(chunk, ref_pos, VariantType::Insertion, 0,
                                                4, alt, alt_ref_base, hp));
            i += gap_len;
            continue;
        }

        int gap_len = 1;
        while (i + gap_len < msa_len &&
               ref_msa_seq[static_cast<size_t>(i + gap_len)] != kGapBase &&
               cons_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase) {
            ++gap_len;
        }
        if (no_end_var && (i - 1 < 0 || i + gap_len >= msa_len ||
                           ref_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                           ref_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase ||
                           cons_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                           cons_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase)) {
            i += gap_len;
            ref_pos += gap_len;
            continue;
        }
        const bool hp = gap_len < opts.min_sv_len &&
                        var_is_homopolymer_indel(chunk, ref_pos, VariantType::Deletion, gap_len, {});
        const uint8_t alt_ref_base = i - 1 >= 0 ? cons_msa_seq[static_cast<size_t>(i - 1)] : 4;
        vars.push_back(make_noisy_candidate(chunk, ref_pos, VariantType::Deletion, gap_len,
                                            4, {}, alt_ref_base, hp));
        i += gap_len;
        ref_pos += gap_len;
    }
    return vars;
}

std::vector<CandidateVariant> make_cand_vars_from_msa(const Options& opts,
                                                      const BamChunk& chunk,
                                                      hts_pos_t noisy_reg_beg,
                                                      const std::vector<uint8_t>& ref_msa_seq,
                                                      const std::vector<uint8_t>& cons_msa_seq,
                                                      int msa_len,
                                                      bool no_end_var) {
    std::vector<uint8_t> ref_filtered;
    std::vector<uint8_t> cons_filtered;
    ref_filtered.reserve(static_cast<size_t>(msa_len));
    cons_filtered.reserve(static_cast<size_t>(msa_len));
    for (int i = 0; i < msa_len; ++i) {
        if (ref_msa_seq[static_cast<size_t>(i)] != kGapBase ||
            cons_msa_seq[static_cast<size_t>(i)] != kGapBase) {
            ref_filtered.push_back(ref_msa_seq[static_cast<size_t>(i)]);
            cons_filtered.push_back(cons_msa_seq[static_cast<size_t>(i)]);
        }
    }
    return make_cand_vars_from_baln0(opts, chunk, noisy_reg_beg, ref_filtered, cons_filtered, no_end_var);
}

int is_match_aln_str(const AlnStr& aln_str,
                     int target_pos,
                     int len,
                     float cons_sim_thres,
                     int* full_cover) {
    int cur_pos = -1;
    int n_eq = 0;
    int n_xid = 0;
    int cover_start = 0;
    int cover_end = 0;
    const int start_pos = target_pos < 0 ? 0 : target_pos;
    const int end_pos = target_pos < 0 ? len - 1 : target_pos + len - 1;

    for (int i = 0; i < aln_str.aln_len; ++i) {
        if (aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_pos;
        if (cur_pos == target_pos + len) break;
        if (i < aln_str.query_beg || i < aln_str.target_beg) continue;
        if (i > aln_str.query_end || i > aln_str.target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cur_pos >= target_pos) {
            if (aln_str.query_aln[static_cast<size_t>(i)] == aln_str.target_aln[static_cast<size_t>(i)]) ++n_eq;
            else ++n_xid;
        }
    }
    *full_cover = cover_start && cover_end;
    if (len >= 10) {
        if (n_eq >= len * cons_sim_thres) return 1;
        return *full_cover ? 0 : -1;
    }
    if (n_eq == len && n_xid == 0) return 1;
    return *full_cover ? 0 : -1;
}

int is_match_aln_str_del(const AlnStr& aln_str,
                         int target_del_left,
                         int target_del_right,
                         int* full_cover) {
    int cur_pos = -1;
    int started_check_del = 0;
    int n_non_del = 0;
    int cover_start = 0;
    int cover_end = 0;
    const int start_pos = target_del_left < 0 ? 0 : target_del_left;
    const int end_pos = target_del_left < 0 ? target_del_right : target_del_right;

    for (int i = 0; i < aln_str.aln_len; ++i) {
        if (aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_pos;
        if (cur_pos > target_del_right) break;
        if (i < aln_str.query_beg || i < aln_str.target_beg) continue;
        if (i > aln_str.query_end || i > aln_str.target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cur_pos >= target_del_left && cur_pos < target_del_right) {
            if (started_check_del == 0) {
                started_check_del = 1;
            } else if (aln_str.query_aln[static_cast<size_t>(i)] != kGapBase) {
                ++n_non_del;
            }
        }
    }
    if (cover_start && cover_end) {
        *full_cover = 1;
        return n_non_del == 0 ? 1 : 0;
    }
    *full_cover = 0;
    return -1;
}

int get_var_allele_i_from_cons_aln_str(const AlnStr& cons_aln_str,
                                       const CandidateVariant& var,
                                       int alt_pos,
                                       float cons_sim_thres,
                                       int* full_cover) {
    *full_cover = 0;
    if (var.key.type == VariantType::Snp) {
        assert(alt_len_for_key(var.key) == 1);
        return is_match_aln_str(cons_aln_str, alt_pos, 1, cons_sim_thres, full_cover);
    }
    if (var.key.type == VariantType::Insertion) {
        return is_match_aln_str(cons_aln_str, alt_pos, alt_len_for_key(var.key),
                                cons_sim_thres, full_cover);
    }
    return is_match_aln_str_del(cons_aln_str, alt_pos - 1, alt_pos, full_cover);
}

int is_cover_aln_str(const AlnStr& aln_str, int target_pos, int len) {
    int cur_pos = -1;
    int cover_start = 0;
    int cover_end = 0;
    const int start_pos = target_pos < 0 ? 0 : target_pos;
    const int end_pos = target_pos < 0 ? len - 1 : target_pos + len - 1;

    for (int i = 0; i < aln_str.aln_len; ++i) {
        if (aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_pos;
        if (i < aln_str.query_beg || i < aln_str.target_beg) continue;
        if (i > aln_str.query_end || i > aln_str.target_end) break;
        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cover_start && cover_end) return 1;
    }
    return 0;
}

int get_full_cover_from_cons_aln_str(const AlnStr& cons_aln_str,
                                     const CandidateVariant& var,
                                     int alt_pos) {
    if (var.key.type == VariantType::Snp) {
        assert(var.key.ref_len == 1);
        return is_cover_aln_str(cons_aln_str, alt_pos, 1);
    }
    if (var.key.type == VariantType::Insertion) {
        return is_cover_aln_str(cons_aln_str, alt_pos, var.key.ref_len + 1);
    }
    return is_cover_aln_str(cons_aln_str, alt_pos - 1, var.key.ref_len + 1);
}

int get_full_cover_from_ref_cons_aln_str(const AlnStr& cons_aln_str,
                                         const AlnStr& ref_cons_aln_str,
                                         int beg_in_ref,
                                         int end_in_ref) {
    int cur_ref_pos = -1;
    int cur_cons_pos = -1;
    int beg_in_cons = -1;
    int end_in_cons = -1;
    int reach_end = 0;
    for (int i = 0; i < ref_cons_aln_str.aln_len; ++i) {
        if (ref_cons_aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_ref_pos;
        if (ref_cons_aln_str.query_aln[static_cast<size_t>(i)] != kGapBase) ++cur_cons_pos;
        if (i < ref_cons_aln_str.query_beg || i < ref_cons_aln_str.target_beg) continue;
        if (i > ref_cons_aln_str.query_end || i > ref_cons_aln_str.target_end) break;

        if (cur_ref_pos == beg_in_ref && beg_in_cons == -1) beg_in_cons = cur_cons_pos;
        if (cur_ref_pos == end_in_ref) reach_end = 1;
        if (reach_end && ref_cons_aln_str.query_aln[static_cast<size_t>(i)] != kGapBase) {
            end_in_cons = cur_cons_pos;
            break;
        }
    }
    return is_cover_aln_str(cons_aln_str, beg_in_cons, end_in_cons - beg_in_cons + 1);
}

void update_cand_var_profile_from_cons_aln_str(const AlnStr& cons_aln_str,
                                               hts_pos_t ref_pos_beg,
                                               std::vector<CandidateVariant>& vars,
                                               ReadVariantProfile& profile) {
    constexpr float cons_sim_thres = 0.9f;
    int delta_ref_alt = 0;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        CandidateVariant& var = vars[static_cast<size_t>(i)];
        const int var_ref_pos = static_cast<int>(var.key.pos - ref_pos_beg);
        int full_cover = 0;
        const int allele_i = get_var_allele_i_from_cons_aln_str(
            cons_aln_str, var, var_ref_pos - delta_ref_alt, cons_sim_thres, &full_cover);
        if (full_cover) {
            ++var.counts.total_cov;
            if (var.counts.alle_covs.empty()) var.counts.alle_covs.assign(2, 0);
            if (allele_i >= 0 && allele_i < static_cast<int>(var.counts.alle_covs.size()))
                ++var.counts.alle_covs[static_cast<size_t>(allele_i)];
            update_read_var_profile_with_allele(i, allele_i, -1, profile);
        }
        if (var.key.type == VariantType::Insertion) delta_ref_alt -= alt_len_for_key(var.key);
        else if (var.key.type == VariantType::Deletion) delta_ref_alt += var.key.ref_len;
    }
}

void update_cand_var_profile_from_cons_aln_str1(int clu_n_seqs,
                                                const std::vector<int>& clu_read_ids,
                                                const std::vector<AlnStr>& clu_aln_strs,
                                                hts_pos_t ref_pos_beg,
                                                std::vector<CandidateVariant>& noisy_vars,
                                                std::vector<ReadVariantProfile>& profiles) {
    for (int i = 0; i < clu_n_seqs; ++i) {
        const size_t cons_read_slot = static_cast<size_t>((i + 1) * 2 - 1);
        if (cons_read_slot >= clu_aln_strs.size()) continue;
        const int read_id = clu_read_ids[static_cast<size_t>(i)];
        if (read_id < 0 || static_cast<size_t>(read_id) >= profiles.size()) continue;
        update_cand_var_profile_from_cons_aln_str(
            clu_aln_strs[cons_read_slot], ref_pos_beg, noisy_vars, profiles[static_cast<size_t>(read_id)]);
    }
}

void update_cand_var_profile_from_cons_aln_str21(int clu_idx,
                                                 const AlnStr& cons_aln_str,
                                                 const AlnStr& ref_cons_aln_str,
                                                 hts_pos_t ref_pos_beg,
                                                 std::vector<CandidateVariant>& vars,
                                                 const std::vector<int>& var_from_cons_idx,
                                                 ReadVariantProfile& profile) {
    constexpr float cons_sim_thres = 0.9f;
    int delta_ref_alt = 0;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        CandidateVariant& var = vars[static_cast<size_t>(i)];
        const int var_beg_in_ref_str = static_cast<int>(var.key.pos - ref_pos_beg);
        const int var_end_in_ref_str = var.key.type == VariantType::Insertion
                                           ? var_beg_in_ref_str
                                           : var_beg_in_ref_str + var.key.ref_len - 1;
        int full_cover = 0;
        int allele_i = -1;
        if (var_from_cons_idx[static_cast<size_t>(i)] & clu_idx) {
            allele_i = get_var_allele_i_from_cons_aln_str(
                cons_aln_str, var, var_beg_in_ref_str - delta_ref_alt,
                cons_sim_thres, &full_cover);
        } else {
            if (var.key.type != VariantType::Deletion) {
                full_cover = get_full_cover_from_cons_aln_str(
                    cons_aln_str, var, var_beg_in_ref_str - delta_ref_alt);
            } else {
                full_cover = get_full_cover_from_ref_cons_aln_str(
                    cons_aln_str, ref_cons_aln_str,
                    var_beg_in_ref_str - 1, var_end_in_ref_str + 1);
            }
            allele_i = 0;
        }
        if (full_cover) {
            ++var.counts.total_cov;
            if (var.counts.alle_covs.empty()) var.counts.alle_covs.assign(2, 0);
            if (allele_i >= 0 && allele_i < static_cast<int>(var.counts.alle_covs.size()))
                ++var.counts.alle_covs[static_cast<size_t>(allele_i)];
            update_read_var_profile_with_allele(i, allele_i, -1, profile);
        }
        if (var_from_cons_idx[static_cast<size_t>(i)] & clu_idx) {
            if (var.key.type == VariantType::Insertion) delta_ref_alt -= alt_len_for_key(var.key);
            else if (var.key.type == VariantType::Deletion) delta_ref_alt += var.key.ref_len;
        }
    }
}

int update_cand_var_profile_from_cons_aln_str2(const Options& opts,
                                               const BamChunk& chunk,
                                               const std::array<int, 2>& clu_n_seqs,
                                               const std::array<std::vector<int>, 2>& clu_read_ids,
                                               const std::array<std::vector<AlnStr>, 2>& aln_strs,
                                               hts_pos_t noisy_reg_beg,
                                               std::vector<CandidateVariant>& hap1_vars,
                                               std::vector<CandidateVariant>& hap2_vars,
                                               std::vector<CandidateVariant>& noisy_vars,
                                               std::vector<VariantCategory>& noisy_var_cate,
                                               std::vector<ReadVariantProfile>& profiles) {
    if (hap1_vars.empty() && hap2_vars.empty()) return 0;

    std::vector<int> var_from_cons_idx;
    noisy_vars.clear();
    noisy_var_cate.clear();
    var_from_cons_idx.reserve(hap1_vars.size() + hap2_vars.size());

    size_t i1 = 0, i2 = 0;
    while (i1 < hap1_vars.size() && i2 < hap2_vars.size()) {
        const int ret = exact_comp_var_site(&hap1_vars[i1].key, &hap2_vars[i2].key);
        if (ret < 0) {
            set_noisy_category(hap1_vars[i1], VariantCategory::NoisyCandHet);
            noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
            var_from_cons_idx.push_back(1);
            noisy_vars.push_back(hap1_vars[i1++]);
        } else if (ret > 0) {
            set_noisy_category(hap2_vars[i2], VariantCategory::NoisyCandHet);
            noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
            var_from_cons_idx.push_back(2);
            noisy_vars.push_back(hap2_vars[i2++]);
        } else {
            set_noisy_category(hap1_vars[i1], VariantCategory::NoisyCandHom);
            noisy_var_cate.push_back(VariantCategory::NoisyCandHom);
            var_from_cons_idx.push_back(3);
            noisy_vars.push_back(hap1_vars[i1++]);
            ++i2;
        }
    }
    while (i1 < hap1_vars.size()) {
        set_noisy_category(hap1_vars[i1], VariantCategory::NoisyCandHet);
        noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
        var_from_cons_idx.push_back(1);
        noisy_vars.push_back(hap1_vars[i1++]);
    }
    while (i2 < hap2_vars.size()) {
        set_noisy_category(hap2_vars[i2], VariantCategory::NoisyCandHet);
        noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
        var_from_cons_idx.push_back(2);
        noisy_vars.push_back(hap2_vars[i2++]);
    }

    profiles = init_read_profiles(chunk.reads.size());
    for (int ci = 0; ci < 2; ++ci) {
        const std::vector<AlnStr>& clu_aln_strs = aln_strs[static_cast<size_t>(ci)];
        if (clu_aln_strs.empty()) continue;
        const AlnStr& ref_cons_aln_str = clu_aln_strs[0];
        for (int j = 0; j < clu_n_seqs[static_cast<size_t>(ci)]; ++j) {
            const size_t cons_read_slot = static_cast<size_t>((j + 1) * 2 - 1);
            if (cons_read_slot >= clu_aln_strs.size()) continue;
            const int read_id = clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(j)];
            if (read_id < 0 || static_cast<size_t>(read_id) >= profiles.size()) continue;
            update_cand_var_profile_from_cons_aln_str21(
                ci + 1, clu_aln_strs[cons_read_slot], ref_cons_aln_str,
                noisy_reg_beg, noisy_vars, var_from_cons_idx,
                profiles[static_cast<size_t>(read_id)]);
        }
    }
    for (CandidateVariant& var : noisy_vars) update_variant_depth_fields(var);

    if (opts.verbose >= 2 && !noisy_vars.empty()) {
        for (size_t i = 0; i < noisy_vars.size(); ++i) {
            const CandidateVariant& v = noisy_vars[i];
            const VariantCounts& c = v.counts;
            std::fprintf(stderr, "Var: %" PRId64 ", %d-%c-%d %d(%d,%d) %s\t%s\n",
                         static_cast<int64_t>(v.key.pos),
                         v.key.ref_len,
                         var_type_code_local(v.key.type),
                         alt_len_for_key(v.key),
                         c.total_cov,
                         c.ref_cov,
                         c.alt_cov,
                         category_name_local(c.category),
                         v.key.alt.c_str());
        }
        for (int hap = 1; hap <= 2; ++hap) {
            const int ci = hap - 1;
            for (int j = 0; j < clu_n_seqs[static_cast<size_t>(ci)]; ++j) {
                const int read_id = clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(j)];
                if (read_id < 0 || static_cast<size_t>(read_id) >= chunk.reads.size()) continue;
                const ReadVariantProfile& p = profiles[static_cast<size_t>(read_id)];
                std::fprintf(stderr, "Hap%d-Read: %s start_var_i: %d, end_var_i: %d\n",
                             hap,
                             chunk.reads[static_cast<size_t>(read_id)].qname.c_str(),
                             p.start_var_idx,
                             p.end_var_idx);
                if (p.start_var_idx < 0) continue;
                for (int k = 0; k <= (p.end_var_idx - p.start_var_idx); ++k) {
                    const int var_i = p.start_var_idx + k;
                    if (var_i < 0 || static_cast<size_t>(var_i) >= noisy_vars.size()) continue;
                    const CandidateVariant& v = noisy_vars[static_cast<size_t>(var_i)];
                    const int allele = p.alleles[static_cast<size_t>(k)];
                    std::fprintf(stderr, "P\tVar: (%d) %" PRId64 " %d-%c-%d, allele: %d\n",
                                 k,
                                 static_cast<int64_t>(v.key.pos),
                                 v.key.ref_len,
                                 var_type_code_local(v.key.type),
                                 alt_len_for_key(v.key),
                                 allele);
                }
            }
        }
    }
    return static_cast<int>(noisy_vars.size());
}

void merge_read_var_profile_entries(const ReadVariantProfile* old_profile,
                                    const std::vector<int>& old_to_merged,
                                    const ReadVariantProfile* new_profile,
                                    const std::vector<int>& new_to_merged,
                                    ReadVariantProfile& merged_profile,
                                    int merged_var_limit) {
    int old_var_i = old_profile != nullptr ? old_profile->start_var_idx : 1;
    int old_end_var_i = old_profile != nullptr ? old_profile->end_var_idx : 0;
    int new_var_i = new_profile != nullptr ? new_profile->start_var_idx : 1;
    int new_end_var_i = new_profile != nullptr ? new_profile->end_var_idx : 0;

    if (old_profile == nullptr || old_to_merged.empty() || old_var_i < 0 || old_end_var_i < old_var_i) {
        old_var_i = 1;
        old_end_var_i = 0;
    }
    if (new_profile == nullptr || new_to_merged.empty() || new_var_i < 0 || new_end_var_i < new_var_i) {
        new_var_i = 1;
        new_end_var_i = 0;
    }

    while (true) {
        while (old_var_i <= old_end_var_i &&
               (old_var_i >= static_cast<int>(old_to_merged.size()) || old_to_merged[static_cast<size_t>(old_var_i)] < 0)) {
            ++old_var_i;
        }
        while (new_var_i <= new_end_var_i &&
               (new_var_i >= static_cast<int>(new_to_merged.size()) || new_to_merged[static_cast<size_t>(new_var_i)] < 0)) {
            ++new_var_i;
        }

        const int old_merged_i = old_var_i <= old_end_var_i
                                     ? old_to_merged[static_cast<size_t>(old_var_i)]
                                     : merged_var_limit;
        const int new_merged_i = new_var_i <= new_end_var_i
                                     ? new_to_merged[static_cast<size_t>(new_var_i)]
                                     : merged_var_limit;
        if (old_merged_i == merged_var_limit && new_merged_i == merged_var_limit) break;

        if (old_merged_i <= new_merged_i) {
            const int old_profile_i = old_var_i - old_profile->start_var_idx;
            update_read_var_profile_with_allele(
                old_merged_i,
                old_profile->alleles[static_cast<size_t>(old_profile_i)],
                old_profile->alt_qi[static_cast<size_t>(old_profile_i)],
                merged_profile);
            ++old_var_i;
        } else {
            const int new_profile_i = new_var_i - new_profile->start_var_idx;
            update_read_var_profile_with_allele(
                new_merged_i,
                new_profile->alleles[static_cast<size_t>(new_profile_i)],
                new_profile->alt_qi[static_cast<size_t>(new_profile_i)],
                merged_profile);
            ++new_var_i;
        }
    }
}

} // namespace

// ════════════════════════════════════════════════════════════════════════════
// sort_noisy_regs
// ════════════════════════════════════════════════════════════════════════════

std::vector<int> sort_noisy_regs(const BamChunk& chunk) {
    const int n = static_cast<int>(chunk.noisy_regions.size());
    std::vector<int> idx(static_cast<size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);

    // mirrors longcallD bubble sort: primary key = noisy_reg_var_sizes[i] (= cr_label = Interval::label)
    // ascending; secondary key = noisy_reg_lens[i] (= cr_end - cr_start = end - beg) ascending.
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const Interval& a = chunk.noisy_regions[static_cast<size_t>(idx[i])];
            const Interval& b = chunk.noisy_regions[static_cast<size_t>(idx[j])];
            bool do_swap = false;
            if (a.label > b.label) {
                do_swap = true;
            } else if (a.label == b.label) {
                if ((a.end - a.beg) > (b.end - b.beg)) do_swap = true;
            }
            if (do_swap) std::swap(idx[i], idx[j]);
        }
    }
    return idx;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_reg_reads
// ════════════════════════════════════════════════════════════════════════════

std::vector<int> collect_noisy_reg_reads(const BamChunk& chunk,
                                         hts_pos_t beg, hts_pos_t end) {
    std::vector<int> result;
    result.reserve(chunk.reads.size());

    // mirrors longcallD collect_noisy_reg_reads1: iterate ordered_read_ids, skip
    // skipped reads and reads with no overlap (beg > end || end <= beg).
    // ReadRecord::beg/end ↔ longcallD chunk->digars[read_i].beg/end.
    const auto visit = [&](int ri) {
        if (ri < 0 || static_cast<size_t>(ri) >= chunk.reads.size()) return;
        const ReadRecord& r = chunk.reads[static_cast<size_t>(ri)];
        if (r.is_skipped) return;
        if (r.beg > end || r.end <= beg) return;
        result.push_back(ri);
    };
    if (chunk.ordered_read_ids.empty()) {
        for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i)
            visit(static_cast<int>(read_i));
    } else {
        for (int ri : chunk.ordered_read_ids) visit(ri);
    }
    return result;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_reg_ref_bseq
// ════════════════════════════════════════════════════════════════════════════

// 2-bit encoding matching longcallD nst_nt4_table (A=0, C=1, G=2, T/U=3, other=4).
static const uint8_t kNstNt4Table[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x00-0x0F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x10-0x1F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x20-0x2F (space, !"#...)
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x30-0x3F (0-9, :;<=>?)
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4, // 0x40-0x4F (@A-O): A=0,C=1,G=2
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4, // 0x50-0x5F (P-_): T=3,U=3
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4, // 0x60-0x6F (`a-o): a=0,c=1,g=2
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4, // 0x70-0x7F (p-DEL): t=3,u=3
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x80-0x8F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x90-0x9F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xA0-0xAF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xB0-0xBF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xC0-0xCF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xD0-0xDF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xE0-0xEF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xF0-0xFF
};

std::vector<uint8_t> collect_reg_ref_bseq(const BamChunk& chunk,
                                           hts_pos_t& beg, hts_pos_t& end) {
    // mirrors longcallD collect_reg_ref_bseq: clip beg/end to ref slice boundaries.
    if (beg < chunk.ref_beg) beg = chunk.ref_beg;
    if (end > chunk.ref_end) end = chunk.ref_end;

    const hts_pos_t len = end - beg + 1;
    if (len <= 0 || chunk.ref_seq.empty()) return {};

    std::vector<uint8_t> out(static_cast<size_t>(len));
    for (hts_pos_t i = beg; i <= end; ++i) {
        out[static_cast<size_t>(i - beg)] =
            kNstNt4Table[static_cast<unsigned char>(
                chunk.ref_seq[static_cast<size_t>(i - chunk.ref_beg)])];
    }
    return out;
}

// ════════════════════════════════════════════════════════════════════════════
// make_vars_from_msa_cons_aln
// ════════════════════════════════════════════════════════════════════════════

int make_vars_from_msa_cons_aln(
    const Options& opts, BamChunk& chunk,
    int /*n_noisy_reads*/, const std::vector<int>& /*read_ids*/,
    hts_pos_t noisy_reg_beg,
    int n_cons,
    const std::array<int, 2>& clu_n_seqs,
    const std::array<std::vector<int>, 2>& clu_read_ids,
    const std::array<std::vector<AlnStr>, 2>& aln_strs,
    std::vector<CandidateVariant>& noisy_vars,
    std::vector<VariantCategory>& noisy_var_cate,
    std::vector<ReadVariantProfile>& noisy_rvp) {
    noisy_vars.clear();
    noisy_var_cate.clear();
    noisy_rvp.clear();

    // longcallD `make_vars_from_msa_cons_aln` (`collect_var.c`): branch on the `n_cons`
    // returned by `collect_noisy_reg_aln_strs`, not on whether auxiliary vectors are empty.
    if (n_cons == 0) return 0;

    if (n_cons == 1) {
        const std::vector<AlnStr>& clu_aln_strs = aln_strs[0];
        if (clu_aln_strs.empty()) return 0;
        const AlnStr& ref_cons_aln_str = clu_aln_strs[0];
        noisy_vars = make_cand_vars_from_msa(opts, chunk, noisy_reg_beg,
                                             ref_cons_aln_str.target_aln,
                                             ref_cons_aln_str.query_aln,
                                             ref_cons_aln_str.aln_len,
                                             false);
        if (noisy_vars.empty()) return 0;
        noisy_var_cate.assign(noisy_vars.size(), VariantCategory::NoisyCandHom);
        for (CandidateVariant& var : noisy_vars) var.counts.category = VariantCategory::NoisyCandHom;

        noisy_rvp = init_read_profiles(chunk.reads.size());
        update_cand_var_profile_from_cons_aln_str1(
            clu_n_seqs[0], clu_read_ids[0], clu_aln_strs,
            noisy_reg_beg, noisy_vars, noisy_rvp);
        for (CandidateVariant& var : noisy_vars) update_variant_depth_fields(var);
        return static_cast<int>(noisy_vars.size());
    }

    std::vector<CandidateVariant> hap1_vars;
    std::vector<CandidateVariant> hap2_vars;
    if (!aln_strs[0].empty()) {
        const AlnStr& ref_cons_aln_str = aln_strs[0][0];
        hap1_vars = make_cand_vars_from_msa(opts, chunk, noisy_reg_beg,
                                            ref_cons_aln_str.target_aln,
                                            ref_cons_aln_str.query_aln,
                                            ref_cons_aln_str.aln_len,
                                            false);
    }
    if (!aln_strs[1].empty()) {
        const AlnStr& ref_cons_aln_str = aln_strs[1][0];
        hap2_vars = make_cand_vars_from_msa(opts, chunk, noisy_reg_beg,
                                            ref_cons_aln_str.target_aln,
                                            ref_cons_aln_str.query_aln,
                                            ref_cons_aln_str.aln_len,
                                            false);
    }
    return update_cand_var_profile_from_cons_aln_str2(
        opts, chunk, clu_n_seqs, clu_read_ids, aln_strs, noisy_reg_beg,
        hap1_vars, hap2_vars, noisy_vars, noisy_var_cate, noisy_rvp);
}

// ════════════════════════════════════════════════════════════════════════════
// merge_var_profile
// ════════════════════════════════════════════════════════════════════════════

void merge_var_profile(BamChunk& chunk,
                       const std::vector<CandidateVariant>& noisy_vars,
                       const std::vector<VariantCategory>& noisy_var_cate,
                       const std::vector<ReadVariantProfile>& noisy_rvp) {
    if (noisy_vars.empty()) return;

    std::vector<CandidateVariant> new_vars = noisy_vars;
    std::vector<VariantCategory> new_cats = noisy_var_cate;
    if (new_cats.size() != new_vars.size()) {
        new_cats.resize(new_vars.size());
        for (size_t i = 0; i < new_vars.size(); ++i)
            new_cats[i] = new_vars[i].counts.category;
    }

    std::vector<size_t> order(new_vars.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return exact_comp_cand_var(&new_vars[a], &new_vars[b]) < 0;
    });
    std::vector<CandidateVariant> sorted_new;
    std::vector<VariantCategory> sorted_cats;
    sorted_new.reserve(new_vars.size());
    sorted_cats.reserve(new_cats.size());
    for (size_t idx : order) {
        sorted_new.push_back(std::move(new_vars[idx]));
        sorted_cats.push_back(new_cats[idx]);
    }
    new_vars = std::move(sorted_new);
    new_cats = std::move(sorted_cats);

    const CandidateTable old_vars = chunk.candidates;
    const std::vector<ReadVariantProfile> old_profiles = chunk.read_var_profile.empty()
                                                             ? init_read_profiles(chunk.reads.size())
                                                             : chunk.read_var_profile;

    CandidateTable merged_vars;
    merged_vars.reserve(old_vars.size() + new_vars.size());
    std::vector<int> old_to_merged(old_vars.size(), -1);
    std::vector<int> new_to_merged(new_vars.size(), -1);

    size_t old_i = 0;
    size_t new_i = 0;
    while (old_i < old_vars.size() && new_i < new_vars.size()) {
        const int ret = exact_comp_var_site(&old_vars[old_i].key, &new_vars[new_i].key);
        if (ret < 0) {
            old_to_merged[old_i] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(old_vars[old_i++]);
        } else if (ret > 0) {
            set_noisy_category(new_vars[new_i], new_cats[new_i]);
            restore_stashed_initial_if_any(chunk, new_vars[new_i]);
            new_to_merged[new_i] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(new_vars[new_i++]);
        } else {
            // longcallD merge_var_profile: exact-key collision always keeps old variant.
            old_to_merged[old_i] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(old_vars[old_i++]);
            ++new_i;
        }
    }
    while (old_i < old_vars.size()) {
        old_to_merged[old_i] = static_cast<int>(merged_vars.size());
        merged_vars.push_back(old_vars[old_i++]);
    }
    while (new_i < new_vars.size()) {
        set_noisy_category(new_vars[new_i], new_cats[new_i]);
        restore_stashed_initial_if_any(chunk, new_vars[new_i]);
        new_to_merged[new_i] = static_cast<int>(merged_vars.size());
        merged_vars.push_back(new_vars[new_i++]);
    }

    std::vector<ReadVariantProfile> merged_profiles = init_read_profiles(chunk.reads.size());
    const std::vector<ReadVariantProfile> empty_new_profiles =
        noisy_rvp.empty() ? init_read_profiles(chunk.reads.size()) : std::vector<ReadVariantProfile>{};
    const std::vector<ReadVariantProfile>& new_profiles = noisy_rvp.empty() ? empty_new_profiles : noisy_rvp;

    const int n_reads = static_cast<int>(chunk.reads.size());
    for (int ord = 0; ord < (chunk.ordered_read_ids.empty() ? n_reads : static_cast<int>(chunk.ordered_read_ids.size())); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || read_i >= n_reads) continue;
        if (chunk.reads[static_cast<size_t>(read_i)].is_skipped) continue;
        const ReadVariantProfile* old_p = static_cast<size_t>(read_i) < old_profiles.size()
                                              ? &old_profiles[static_cast<size_t>(read_i)]
                                              : nullptr;
        const ReadVariantProfile* new_p = static_cast<size_t>(read_i) < new_profiles.size()
                                              ? &new_profiles[static_cast<size_t>(read_i)]
                                              : nullptr;
        merge_read_var_profile_entries(old_p, old_to_merged,
                                       new_p, new_to_merged,
                                       merged_profiles[static_cast<size_t>(read_i)],
                                       static_cast<int>(merged_vars.size()));
    }

    cgranges_t* merged_cr = cr_init();
    for (int ord = 0; ord < (chunk.ordered_read_ids.empty() ? n_reads : static_cast<int>(chunk.ordered_read_ids.size())); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || read_i >= n_reads) continue;
        if (chunk.reads[static_cast<size_t>(read_i)].is_skipped) continue;
        const ReadVariantProfile& p = merged_profiles[static_cast<size_t>(read_i)];
        if (p.start_var_idx < 0 || p.end_var_idx < p.start_var_idx) continue;
        cr_add(merged_cr, "cr", p.start_var_idx, p.end_var_idx + 1,
               static_cast<int32_t>(read_i));
    }
    cr_index(merged_cr);

    chunk.candidates = std::move(merged_vars);
    chunk.read_var_profile = std::move(merged_profiles);
    chunk.read_var_cr.reset(merged_cr);
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_vars1
// ════════════════════════════════════════════════════════════════════════════

int collect_noisy_vars1(BamChunk& chunk, const Options& opts, int noisy_reg_i) {
    const Interval& reg = chunk.noisy_regions[static_cast<size_t>(noisy_reg_i)];
    // longcallD `collect_noisy_vars1` enters with `noisy_reg_beg = cr_start(chunk_noisy_regs, i)`.
    hts_pos_t noisy_reg_beg = reg.beg;
    hts_pos_t noisy_reg_end = reg.end;

    // mirrors longcallD: skip regions longer than max_noisy_reg_len (return 0 = done, no vars).
    if (noisy_reg_end - noisy_reg_beg + 1 > static_cast<hts_pos_t>(opts.max_noisy_reg_len))
        return 0;

    // longcallD first clips beg/end via collect_reg_ref_bseq, then collects reads on the
    // finalized region bounds.
    std::vector<uint8_t> ref_seq =
        collect_reg_ref_bseq(chunk, noisy_reg_beg, noisy_reg_end);
    if (ref_seq.empty())
        return 0;

    const std::vector<int> read_ids =
        collect_noisy_reg_reads(chunk, noisy_reg_beg, noisy_reg_end);
    const int n_noisy_reads = static_cast<int>(read_ids.size());

    // mirrors longcallD: skip regions with more than max_noisy_reg_cov reads.
    if (n_noisy_reads > opts.max_noisy_reg_cov)
        return 0;

    // MSA + consensus via WFA2-lib / abPOA (align.cpp).
    std::array<int, 2> clu_n_seqs{};
    std::array<std::vector<int>, 2> clu_read_ids{};
    std::array<std::vector<AlnStr>, 2> aln_strs{};
    const int n_cons = collect_noisy_reg_aln_strs(
        opts, chunk, noisy_reg_beg, noisy_reg_end,
        read_ids, ref_seq,
        clu_n_seqs, clu_read_ids, aln_strs);

    // mirrors longcallD: n_cons == 0 → MSA could not resolve; return -1 so the
    // outer loop leaves this region undone and may retry if another region makes progress.
    if (n_cons == 0)
        return -1;

    std::vector<CandidateVariant>    noisy_vars;
    std::vector<VariantCategory>     noisy_var_cate;
    std::vector<ReadVariantProfile>  noisy_rvp;

    const int n_noisy_vars = make_vars_from_msa_cons_aln(
        opts, chunk,
        n_noisy_reads, read_ids, noisy_reg_beg,
        n_cons, clu_n_seqs, clu_read_ids, aln_strs,
        noisy_vars, noisy_var_cate, noisy_rvp);

    merge_var_profile(chunk, noisy_vars, noisy_var_cate, noisy_rvp);

    return n_noisy_vars;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_vars_step4
// ════════════════════════════════════════════════════════════════════════════

void collect_noisy_vars_step4(BamChunk& chunk, const Options& opts) {
    if (chunk.noisy_regions.empty()) return;

    // mirrors longcallD step 4 in collect_var_main.
    const std::vector<int> sorted = sort_noisy_regs(chunk);
    const int n_regs = static_cast<int>(chunk.noisy_regions.size());
    std::vector<bool> done(static_cast<size_t>(n_regs), false);

    while (true) {
        bool any_done = false, any_new_var = false;
        for (int reg_idx : sorted) {
            if (done[static_cast<size_t>(reg_idx)]) continue;
            const int ret = collect_noisy_vars1(chunk, opts, reg_idx);
            // ret >= 0: region attempted (done regardless of whether new vars found).
            // ret < 0:  MSA failed; leave undone so it can be retried if another region
            //           makes progress and changes the phasing context.
            if (ret >= 0) {
                done[static_cast<size_t>(reg_idx)] = true;
                any_done = true;
                if (ret > 0) any_new_var = true;
            }
        }
        // mirrors longcallD: re-run k-means with kCandGermlineVarCate whenever new
        // noisy variants were merged, incorporating NOISY_CAND_HET / NOISY_CAND_HOM.
        if (any_new_var)
            assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
        // stop when no region made progress in this pass.
        if (!any_done) break;
    }
}

} // namespace pgphase_collect
