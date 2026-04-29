/**
 * @file collect_phase.cpp
 * @brief k-means read-haplotype clustering, porting longcallD assign_hap.c.
 */

#include "collect_phase.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cstdlib>
#include <iterator>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

/** longcallD `LONGCALLD_DEF_PLOID` (assign_hap.c flip loop bound). */
constexpr int kLongcalldDefPloid = 2;

/** Read visit count / index: longcallD always uses `ordered_read_ids[i]` with `n_reads` slots. */
static inline int read_visit_count(const BamChunk& chunk) {
    return chunk.ordered_read_ids.empty() ? static_cast<int>(chunk.reads.size())
                                         : static_cast<int>(chunk.ordered_read_ids.size());
}

static inline int read_at_visit_ord(const BamChunk& chunk, int ord) {
    return chunk.ordered_read_ids.empty() ? ord : chunk.ordered_read_ids[static_cast<size_t>(ord)];
}

static bool parse_debug_site_pos(const std::string& site, hts_pos_t& pos_out) {
    const size_t colon = site.find(':');
    if (colon == std::string::npos) return false;
    std::string pos_s = site.substr(colon + 1);
    pos_s.erase(std::remove(pos_s.begin(), pos_s.end(), ','), pos_s.end());
    if (pos_s.empty()) return false;
    try {
        const long long p = std::stoll(pos_s);
        if (p < 1) return false;
        pos_out = static_cast<hts_pos_t>(p);
        return true;
    } catch (...) {
        return false;
    }
}

static uint32_t aux_read_u32_le(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

static int extract_hs_set_ids(const ReadRecord& read, std::vector<uint32_t>& out_ids) {
    out_ids.clear();
    if (!read.alignment) return 0;
    uint8_t* hs = bam_aux_get(read.alignment.get(), "hs");
    if (hs == nullptr || hs[0] != 'B') return 0;
    const uint8_t subtype = hs[1];
    const uint32_t n = aux_read_u32_le(hs + 2);
    const uint8_t* data = hs + 6;
    out_ids.reserve(static_cast<size_t>(n));

    if (subtype == 'I' || subtype == 'i') {
        for (uint32_t i = 0; i < n; ++i) out_ids.push_back(aux_read_u32_le(data + i * 4));
        return 1;
    }
    if (subtype == 'S' || subtype == 's') {
        for (uint32_t i = 0; i < n; ++i) {
            const uint8_t* p = data + i * 2;
            out_ids.push_back(static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8));
        }
        return 1;
    }
    if (subtype == 'C' || subtype == 'c') {
        for (uint32_t i = 0; i < n; ++i) out_ids.push_back(static_cast<uint32_t>(data[i]));
        return 1;
    }
    return 0;
}

static bool collect_thread_set_for_group(const BamChunk& chunk,
                                         int hap,
                                         hts_pos_t phase_set,
                                         const PgbamSidecarData& sidecar,
                                         std::vector<uint64_t>& threads) {
    threads.clear();
    std::vector<uint32_t> hs_ids;
    bool saw_any_valid_read = false;

    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped || chunk.haps[i] != hap || chunk.phase_sets[i] != phase_set) continue;
        const int hs_ok = extract_hs_set_ids(chunk.reads[i], hs_ids);
        if (!hs_ok || hs_ids.empty()) continue;
        bool read_resolved = false;
        for (uint32_t set_id : hs_ids) {
            auto it = sidecar.set_to_threads.find(set_id);
            if (it == sidecar.set_to_threads.end()) continue;
            for (uint64_t tid : it->second) threads.push_back(tid);
            read_resolved = true;
        }
        if (read_resolved) saw_any_valid_read = true;
    }
    if (!threads.empty()) {
        std::sort(threads.begin(), threads.end());
        threads.erase(std::unique(threads.begin(), threads.end()), threads.end());
    }
    return saw_any_valid_read;
}

static int intersection_size(const std::vector<uint64_t>& lhs, const std::vector<uint64_t>& rhs) {
    int n = 0;
    auto li = lhs.begin();
    auto ri = rhs.begin();
    while (li != lhs.end() && ri != rhs.end()) {
        if (*li < *ri) ++li;
        else if (*ri < *li) ++ri;
        else {
            ++n;
            ++li;
            ++ri;
        }
    }
    return n;
}

static void sort_unique_threads(std::vector<uint64_t>& threads) {
    if (threads.empty()) return;
    std::sort(threads.begin(), threads.end());
    threads.erase(std::unique(threads.begin(), threads.end()), threads.end());
}

static std::vector<uint64_t> set_difference_threads(const std::vector<uint64_t>& lhs,
                                                    const std::vector<uint64_t>& rhs) {
    std::vector<uint64_t> out;
    out.reserve(lhs.size());
    auto li = lhs.begin();
    auto ri = rhs.begin();
    while (li != lhs.end()) {
        while (ri != rhs.end() && *ri < *li) ++ri;
        if (ri == rhs.end() || *li < *ri) out.push_back(*li);
        ++li;
    }
    return out;
}

static std::vector<uint64_t> union_threads(const std::vector<uint64_t>& lhs,
                                           const std::vector<uint64_t>& rhs) {
    std::vector<uint64_t> out;
    out.reserve(lhs.size() + rhs.size());
    std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(out));
    return out;
}

struct PhaseBlockThreadState {
    hts_pos_t anchor = -1;
    std::array<std::vector<uint64_t>, 3> hap_threads;
    std::array<std::vector<uint64_t>, 3> hap_unique_threads;
};

static void refresh_phase_block_unique_threads(PhaseBlockThreadState& state) {
    state.hap_unique_threads[1] = set_difference_threads(state.hap_threads[1], state.hap_threads[2]);
    state.hap_unique_threads[2] = set_difference_threads(state.hap_threads[2], state.hap_threads[1]);
}

static void build_read_thread_cache(const BamChunk& chunk,
                                    const PgbamSidecarData& sidecar,
                                    std::vector<std::vector<uint64_t>>& read_threads) {
    read_threads.clear();
    read_threads.resize(chunk.reads.size());
    std::vector<uint32_t> hs_ids;

    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped) continue;
        if (!extract_hs_set_ids(chunk.reads[i], hs_ids) || hs_ids.empty()) continue;
        std::vector<uint64_t>& threads = read_threads[i];
        for (uint32_t set_id : hs_ids) {
            auto it = sidecar.set_to_threads.find(set_id);
            if (it == sidecar.set_to_threads.end()) continue;
            for (uint64_t tid : it->second) threads.push_back(tid);
        }
        sort_unique_threads(threads);
    }
}

static std::vector<PhaseBlockThreadState>
build_phase_block_thread_states(const BamChunk& chunk,
                                const std::vector<hts_pos_t>& psets,
                                const std::vector<std::vector<uint64_t>>& read_threads) {
    std::vector<PhaseBlockThreadState> states(psets.size());
    std::unordered_map<hts_pos_t, size_t> pset_to_state;
    pset_to_state.reserve(psets.size());
    for (size_t i = 0; i < psets.size(); ++i) {
        states[i].anchor = psets[i];
        pset_to_state[psets[i]] = i;
    }

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        if (chunk.reads[read_i].is_skipped) continue;
        const int hap = chunk.haps[read_i];
        if (hap != 1 && hap != 2) continue;
        auto it = pset_to_state.find(chunk.phase_sets[read_i]);
        if (it == pset_to_state.end()) continue;
        const std::vector<uint64_t>& threads = read_threads[read_i];
        states[it->second].hap_threads[hap].insert(states[it->second].hap_threads[hap].end(),
                                                   threads.begin(), threads.end());
    }

    for (PhaseBlockThreadState& state : states) {
        sort_unique_threads(state.hap_threads[1]);
        sort_unique_threads(state.hap_threads[2]);
        refresh_phase_block_unique_threads(state);
    }
    return states;
}

static PhaseBlockThreadState collect_phase_block_thread_state(const BamChunk& chunk,
                                                              hts_pos_t phase_set,
                                                              const PgbamSidecarData& sidecar) {
    PhaseBlockThreadState state;
    state.anchor = phase_set;
    collect_thread_set_for_group(chunk, 1, phase_set, sidecar, state.hap_threads[1]);
    collect_thread_set_for_group(chunk, 2, phase_set, sidecar, state.hap_threads[2]);
    refresh_phase_block_unique_threads(state);
    return state;
}

static bool decide_phase_block_concordance(const PhaseBlockThreadState& left,
                                           const PhaseBlockThreadState& right,
                                           bool& do_flip) {
    const int s11 = intersection_size(left.hap_unique_threads[1], right.hap_unique_threads[1]);
    const int s12 = intersection_size(left.hap_unique_threads[1], right.hap_unique_threads[2]);
    const int s21 = intersection_size(left.hap_unique_threads[2], right.hap_unique_threads[1]);
    const int s22 = intersection_size(left.hap_unique_threads[2], right.hap_unique_threads[2]);
    const int score_nonflip = s11 + s22;
    const int score_flip = s12 + s21;

    if (score_nonflip == score_flip) return false;
    if (std::max(score_nonflip, score_flip) < 2) return false;
    do_flip = score_flip > score_nonflip;
    return true;
}

static PhaseBlockThreadState merge_phase_block_thread_states(const PhaseBlockThreadState& left,
                                                             const PhaseBlockThreadState& right,
                                                             bool do_flip) {
    PhaseBlockThreadState merged;
    merged.anchor = left.anchor;
    const int right_hap_for_1 = do_flip ? 2 : 1;
    const int right_hap_for_2 = do_flip ? 1 : 2;
    merged.hap_threads[1] = union_threads(left.hap_threads[1], right.hap_threads[right_hap_for_1]);
    merged.hap_threads[2] = union_threads(left.hap_threads[2], right.hap_threads[right_hap_for_2]);
    refresh_phase_block_unique_threads(merged);
    return merged;
}

// ════════════════════════════════════════════════════════════════════════════
// Category flag mapping
// ════════════════════════════════════════════════════════════════════════════

uint32_t category_to_flag(VariantCategory c) {
    switch (c) {
        case VariantCategory::CleanHetSnp:   return kCandCleanHetSnp;
        case VariantCategory::CleanHetIndel: return kCandCleanHetIndel;
        case VariantCategory::CleanHom:      return kCandCleanHom;
        case VariantCategory::NoisyCandHet:  return kCandNoisyCandHet;
        case VariantCategory::NoisyCandHom:  return kCandNoisyCandHom;
        default: return 0;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Helper functions (static, matching longcallD assign_hap.c internals)
// ════════════════════════════════════════════════════════════════════════════

// Mirrors longcallD read_init_hap_phase_set.
static void read_init_hap_phase_set(BamChunk& chunk) {
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        chunk.haps[i] = 0;
        chunk.phase_sets[i] = -1;
    }
}

/** Width of `hap_to_alle_profile[1/2]` (longcallD `n_uniq_alles` with optional `alle_covs`). */
static int variant_allele_slots(const CandidateVariant& var) {
    const VariantCounts& c = var.counts;
    if (c.n_uniq_alles > 0) return c.n_uniq_alles;
    if (!c.alle_covs.empty()) return static_cast<int>(c.alle_covs.size());
    return 2;
}

// Returns majority allele index or -1 if inconclusive (longcallD `get_var_init_max_cov_allele`).
// Strict `>` on coverage so ties keep the lower index (ref-first), matching longcallD.
static int get_var_init_max_cov_allele(bool is_ont, const CandidateVariant& var) {
    if (is_ont && var.is_homopolymer_indel) return -1;
    const VariantCounts& c = var.counts;
    int max_cov = 0;
    int max_cov_alle_i = -1;
    if (!c.alle_covs.empty()) {
        const int n = c.n_uniq_alles > 0
                          ? c.n_uniq_alles
                          : static_cast<int>(c.alle_covs.size());
        for (int i = 0; i < n; ++i) {
            if (c.alle_covs[i] > max_cov) {
                max_cov = c.alle_covs[i];
                max_cov_alle_i = i;
            }
        }
        return max_cov_alle_i;
    }
    const int n = (c.n_uniq_alles > 0) ? c.n_uniq_alles : 2;
    for (int i = 0; i < n; ++i) {
        const int cov = (i == 0) ? c.ref_cov : (i == 1 ? c.alt_cov : 0);
        if (cov > max_cov) {
            max_cov = cov;
            max_cov_alle_i = i;
        }
    }
    return max_cov_alle_i;
}

// Initialize hap_to_alle_profile (zeroed) and hap_to_cons_alle for all valid vars.
// Hom vars get cons_alle[1]=cons_alle[2]=1; het vars get -1/-1.
// Mirrors longcallD var_init_hap_profile_cons_allele.
static void var_init_hap_profile_cons_allele(bool is_ont,
                                              CandidateTable& variants,
                                              const std::vector<int>& valid_var_idx) {
    for (int vi : valid_var_idx) {
        CandidateVariant& var = variants[vi];
        const int na = variant_allele_slots(var);
        const bool first_init = var.hap_to_alle_profile[1].empty() && var.hap_to_alle_profile[2].empty();
        if (first_init) {
            var.hap_to_alle_profile[0].assign(na, 0);
        }
        var.hap_to_alle_profile[1].assign(na, 0);
        var.hap_to_alle_profile[2].assign(na, 0);
        var.hap_to_cons_alle[0] = get_var_init_max_cov_allele(is_ont, var);
        const uint32_t flag = category_to_flag(var.counts.category);
        if (flag == kCandCleanHom || flag == kCandNoisyCandHom) {
            var.hap_to_cons_alle[1] = 1;
            var.hap_to_cons_alle[2] = 1;
        } else {
            var.hap_to_cons_alle[1] = -1;
            var.hap_to_cons_alle[2] = -1;
        }
    }
}

// Zero all hap allele profiles (0, 1, 2) for the valid vars; preserves cons_alle.
// Mirrors longcallD var_init_hap_to_alle_profile (j = 0..PLOID).
static void var_init_hap_to_alle_profile(CandidateTable& variants,
                                          const std::vector<int>& valid_var_idx) {
    for (int vi : valid_var_idx) {
        CandidateVariant& v = variants[vi];
        const int na = variant_allele_slots(v);
        v.hap_to_alle_profile[0].assign(na, 0);
        v.hap_to_alle_profile[1].assign(na, 0);
        v.hap_to_alle_profile[2].assign(na, 0);
    }
}

// Pick pivot: deepest CleanHetSnp > CleanHetIndel > NoisyHetSnp > NoisyHetIndel.
// Returns index into valid_var_idx, or -1.
// Mirrors longcallD select_init_var.
static int select_init_var(const CandidateTable& variants,
                            const std::vector<int>& valid_var_idx) {
    int snp_i = -1, indel_i = -1, noisy_snp_i = -1, noisy_indel_i = -1;
    int snp_dp = 0, indel_dp = 0, noisy_snp_dp = 0, noisy_indel_dp = 0;
    for (int _vi = 0; _vi < (int)valid_var_idx.size(); ++_vi) {
        const CandidateVariant& var = variants[valid_var_idx[_vi]];
        const uint32_t flag = category_to_flag(var.counts.category);
        if (flag == kCandCleanHetSnp) {
            if (snp_i == -1 || var.counts.total_cov > snp_dp) {
                snp_i = _vi; snp_dp = var.counts.total_cov;
            }
        } else if (flag == kCandCleanHetIndel) {
            if (indel_i == -1 || var.counts.total_cov > indel_dp) {
                indel_i = _vi; indel_dp = var.counts.total_cov;
            }
        } else if (flag == kCandNoisyCandHet) {
            if (var.key.type == VariantType::Snp) {
                if (noisy_snp_i == -1 || var.counts.total_cov > noisy_snp_dp) {
                    noisy_snp_i = _vi; noisy_snp_dp = var.counts.total_cov;
                }
            } else if (!var.is_homopolymer_indel) {
                if (noisy_indel_i == -1 || var.counts.total_cov > noisy_indel_dp) {
                    noisy_indel_i = _vi; noisy_indel_dp = var.counts.total_cov;
                }
            }
        }
    }
    if (snp_i != -1) return snp_i;
    if (indel_i != -1) return indel_i;
    if (noisy_snp_i != -1) return noisy_snp_i;
    return noisy_indel_i;
}

// Update hap_to_cons_alle[hap] via argmax of the allele profile, with ONT 67% guard.
// Mirrors longcallD update_var_hap_to_cons_alle.
static void update_var_hap_to_cons_alle(bool is_ont, CandidateVariant& var, int hap) {
    if (hap == 0) return;
    const auto& prof = var.hap_to_alle_profile[hap];
    int max_cov = 0, max_alle = -1, total = 0;
    for (size_t a = 0; a < prof.size(); ++a) {
        total += prof[a];
        if (prof[a] > max_cov) {
            max_cov = prof[a];
            max_alle = static_cast<int>(a);
        }
    }
    if (is_ont && var.is_homopolymer_indel && max_cov < total * 0.67) max_alle = -1;
    var.hap_to_cons_alle[hap] = max_alle;
}

// Score a read allele against the hap consensus; infers complement when one hap is unknown.
// Returns +var_score (match), −var_score (mismatch), or 0 (no info).
// Side effect: may set hap_to_cons_alle[hap] or [3-hap] when one is -1.
// Mirrors longcallD read_to_cons_allele_score.
static int read_to_cons_allele_score(CandidateVariant& var, int hap, uint32_t vflag, int allele_i) {
    const int var_score = (vflag == kCandCleanHetSnp || vflag == kCandCleanHetIndel) ? 2 : 1;
    if (var.hap_to_cons_alle[hap] == -1 && var.hap_to_cons_alle[3 - hap] == -1) return 0;
    if (var.hap_to_cons_alle[hap] == -1) var.hap_to_cons_alle[hap] = 1 - var.hap_to_cons_alle[3 - hap];
    if (var.hap_to_cons_alle[3 - hap] == -1) var.hap_to_cons_alle[3 - hap] = 1 - var.hap_to_cons_alle[hap];
    if (var.hap_to_cons_alle[hap] == allele_i) return var_score;
    if (var.hap_to_cons_alle[hap] == -1) return 0;
    return -var_score;
}

// Assign a read to hap 1, 2, 0 (tied), or -1 (no informative variants).
// Updates n_clean_agree_snps / n_clean_conflict_snps on the read (for max_hap only).
// CleanHom variants contribute to agree/conflict stats but not to hap_scores.
// Mirrors longcallD init_assign_read_hap_based_on_cons_alle.
static int init_assign_read_hap(BamChunk& chunk, int read_i, uint32_t flags) {
    ReadRecord& read = chunk.reads[read_i];
    read.n_clean_agree_snps = 0;
    read.n_clean_conflict_snps = 0;

    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return -1;

    int hap_scores[3] = {0, 0, 0};
    int n_vars_used[3] = {0, 0, 0};
    int n_clean_agree[3] = {0, 0, 0};
    int n_clean_conflict[3] = {0, 0, 0};

    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        CandidateVariant& var = chunk.candidates[vi];
        const uint32_t vflag = category_to_flag(var.counts.category);
        if ((vflag & flags) == 0) continue;
        // Homopolymer indels and NoisyCandHom are skipped entirely.
        if (var.is_homopolymer_indel || vflag == kCandNoisyCandHom) continue;

        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;

        for (int hap = 1; hap <= 2; ++hap) {
            const int score = read_to_cons_allele_score(var, hap, vflag, aidx);
            // CleanHom does not contribute to hap_scores or n_vars_used.
            if (vflag != kCandCleanHom) {
                if (score != 0) n_vars_used[hap]++;
                hap_scores[hap] += score;
            }
            // Count agree/conflict for all clean SNPs (including CleanHom SNPs).
            if ((vflag & kCandGermlineClean) && var.key.type == VariantType::Snp && score != 0) {
                if (score > 0) n_clean_agree[hap]++;
                else n_clean_conflict[hap]++;
            }
        }
    }

    int max_hap = 0, max_score = 0, min_hap = 0, min_score = 0;
    for (int hap = 1; hap <= 2; ++hap) {
        if (hap_scores[hap] > max_score) { max_hap = hap; max_score = hap_scores[hap]; }
        else if (hap_scores[hap] < min_score) { min_hap = hap; min_score = hap_scores[hap]; }
    }

    if (n_vars_used[1] == 0 && n_vars_used[2] == 0) return -1;
    if (max_score == 0 && min_score == 0) return 0;
    if (max_score > 0) {
        read.n_clean_agree_snps = n_clean_agree[max_hap];
        read.n_clean_conflict_snps = n_clean_conflict[max_hap];
        return max_hap;
    }
    return 3 - min_hap;
}

// Update allele profile + cons_alle for all vars a read covers (used in Phase 1).
// hap==0 means unassigned — both hap 1 and 2 are updated identically.
// Mirrors longcallD update_var_hap_profile_cons_alle_based_on_read_hap.
static void update_var_hap_profile_cons_alle(BamChunk& chunk, bool is_ont,
                                              int read_i, int hap, uint32_t flags) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return;
    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        const uint32_t vflag = category_to_flag(chunk.candidates[vi].counts.category);
        if ((vflag & flags) == 0) continue;
        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;
        CandidateVariant& var = chunk.candidates[vi];
        if (hap == 0) {
            for (int h = 1; h <= 2; ++h) {
                if (aidx >= static_cast<int>(var.hap_to_alle_profile[h].size())) continue;
                var.hap_to_alle_profile[h][aidx]++;
                update_var_hap_to_cons_alle(is_ont, var, h);
            }
        } else {
            if (aidx >= static_cast<int>(var.hap_to_alle_profile[hap].size())) continue;
            var.hap_to_alle_profile[hap][aidx]++;
            update_var_hap_to_cons_alle(is_ont, var, hap);
        }
    }
}

// Update allele profile only for all vars a read covers (used in Phase 2).
// Mirrors longcallD update_var_hap_profile_based_on_read_hap.
static void update_var_hap_profile(BamChunk& chunk, int read_i, int hap, uint32_t flags) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return;
    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        const uint32_t vflag = category_to_flag(chunk.candidates[vi].counts.category);
        if ((vflag & flags) == 0) continue;
        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;
        CandidateVariant& cand = chunk.candidates[vi];
        if (hap == 0) {
            if (aidx < static_cast<int>(cand.hap_to_alle_profile[1].size()))
                cand.hap_to_alle_profile[1][aidx]++;
            if (aidx < static_cast<int>(cand.hap_to_alle_profile[2].size()))
                cand.hap_to_alle_profile[2][aidx]++;
        } else {
            if (aidx < static_cast<int>(cand.hap_to_alle_profile[hap].size()))
                cand.hap_to_alle_profile[hap][aidx]++;
        }
    }
}

// Returns 1=agree, 0=conflict, -1=uninformative.
// Mirrors longcallD check_agree_haps.
static int check_agree_haps(const BamChunk& chunk, int read_i, int hap, int var1, int var2) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (var1 < prof.start_var_idx || var2 > prof.end_var_idx) return -1;
    if (hap == 0) return -1;
    const int a1 = prof.alleles[var1 - prof.start_var_idx];
    const int a2 = prof.alleles[var2 - prof.start_var_idx];
    if (a1 < 0 || a2 < 0) return -1;
    const CandidateVariant& v1 = chunk.candidates[var1];
    const CandidateVariant& v2 = chunk.candidates[var2];
    const bool agree   = (v1.hap_to_cons_alle[hap] == a1 && v2.hap_to_cons_alle[hap] == a2);
    const bool conflict = (v1.hap_to_cons_alle[hap] == a1 && v2.hap_to_cons_alle[3 - hap] == a2);
    if (agree) return 1;
    if (conflict) return 0;
    return -1;
}

// Phase-set assignment + flip for one k-means iteration.
// Returns 1 if any flip occurred (changed), 0 if converged.
// Mirrors longcallD iter_update_var_hap_cons_phase_set.
static int iter_update_var_hap_cons_phase_set(BamChunk& chunk,
                                               const std::vector<int>& valid_var_idx,
                                               const Options& opts) {
    const int n = (int)valid_var_idx.size();

    // Collect het var positions (indices into valid_var_idx).
    std::vector<int> het_var_idx;
    std::vector<bool> is_het(n, false);
    for (int _vi = 0; _vi < n; ++_vi) {
        const CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
        if (var.hap_to_cons_alle[1] != -1 && var.hap_to_cons_alle[2] != -1 &&
            var.hap_to_cons_alle[1] != var.hap_to_cons_alle[2] && !var.is_homopolymer_indel) {
            is_het[_vi] = true;
            het_var_idx.push_back(_vi);
        }
    }

    std::vector<int> n_agree(n, 0), n_conflict(n, 0);
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    cgranges_t* cr = chunk.read_var_cr.get();

    // For each adjacent pair of het vars, count spanning read agree/conflict.
    for (int hi = 1; hi < (int)het_var_idx.size(); ++hi) {
        const int _vi = het_var_idx[hi];
        const int vi = valid_var_idx[_vi];
        const int prev_vi = valid_var_idx[het_var_idx[hi - 1]];

        const int64_t ovlp_n = cr_overlap(cr, "cr", prev_vi, vi + 1, &ovlp_b, &max_b);
        for (int64_t oi = 0; oi < ovlp_n; ++oi) {
            const int read_i = (int)cr_label(cr, ovlp_b[oi]);
            if (chunk.reads[read_i].is_skipped) continue;
            const int agree = check_agree_haps(chunk, read_i, chunk.haps[read_i], prev_vi, vi);
            if (agree > 0) n_agree[_vi]++;
            else if (agree == 0) n_conflict[_vi]++;
        }
    }
    free(ovlp_b);

    int changed = 0;
    int flip = 0;
    hts_pos_t phase_set = -1;
    for (int _vi = 0; _vi < n; ++_vi) {
        const int vi = valid_var_idx[_vi];
        CandidateVariant& var = chunk.candidates[vi];
        if (_vi == 0) {
            phase_set = var.key.sort_pos();
            var.phase_set = phase_set;
            continue;
        }
        if (opts.verbose >= 2) {
            std::fprintf(stderr, "%" PRIi64 " %d %d\n",
                         static_cast<int64_t>(var.key.pos),
                         n_agree[_vi], n_conflict[_vi]);
        }
        if (is_het[_vi]) {
            if (n_agree[_vi] < 2 && n_conflict[_vi] < 2) {
                // Insufficient spanning reads — start a new phase set.
                phase_set = var.key.sort_pos();
            } else if (n_conflict[_vi] > n_agree[_vi]) {
                flip ^= 1;
            }
            if (flip == 1) {
                changed = 1;
                // longcallD: swap hap_to_cons_alle[hap] with [3-hap] for hap=1..LONGCALLD_DEF_PLOID
                // (with ploidy 2 this is two swaps → net identity; matches released assign_hap.c).
                for (int hap = 1; hap <= kLongcalldDefPloid; ++hap) {
                    const int tmp = var.hap_to_cons_alle[hap];
                    var.hap_to_cons_alle[hap] = var.hap_to_cons_alle[3 - hap];
                    var.hap_to_cons_alle[3 - hap] = tmp;
                }
            }
        }
        var.phase_set = phase_set;
    }
    return changed;
}

// Rebuild allele profiles and consensus; returns 1 if any cons_alle changed.
// Mirrors longcallD iter_update_var_hap_to_cons_alle.
static int iter_update_var_hap_to_cons_alle(BamChunk& chunk, bool is_ont,
                                             const std::vector<int>& valid_var_idx,
                                             uint32_t flags) {
    const int n = (int)valid_var_idx.size();

    // Save current consensus for convergence check.
    std::vector<std::array<int, 3>> saved(n);
    for (int _vi = 0; _vi < n; ++_vi)
        saved[_vi] = chunk.candidates[valid_var_idx[_vi]].hap_to_cons_alle;

    var_init_hap_to_alle_profile(chunk.candidates, valid_var_idx);

    for (int oi = 0; oi < read_visit_count(chunk); ++oi) {
        const int read_i = read_at_visit_ord(chunk, oi);
        if (read_i < 0 || static_cast<size_t>(read_i) >= chunk.reads.size()) continue;
        if (chunk.reads[read_i].is_skipped) continue;
        int hap = init_assign_read_hap(chunk, read_i, flags);
        if (hap == -1) hap = 0;
        chunk.haps[read_i] = hap;
        update_var_hap_profile(chunk, read_i, hap, flags);
    }

    for (int _vi = 0; _vi < n; ++_vi) {
        CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
        for (int hap = 1; hap <= 2; ++hap)
            update_var_hap_to_cons_alle(is_ont, var, hap);
    }

    int changed = 0;
    for (int _vi = 0; _vi < n && !changed; ++_vi) {
        const auto& cur = chunk.candidates[valid_var_idx[_vi]].hap_to_cons_alle;
        for (int hap = 1; hap <= 2; ++hap)
            if (cur[hap] != saved[_vi][hap]) { changed = 1; break; }
    }
    return changed;
}

// Assign per-read phase_sets from the first phased het var each read covers.
// Mirrors longcallD update_read_phase_set.
static void update_read_phase_set(BamChunk& chunk, const std::vector<bool>& var_is_valid) {
    for (int oi = 0; oi < read_visit_count(chunk); ++oi) {
        const int read_i = read_at_visit_ord(chunk, oi);
        if (read_i < 0 || static_cast<size_t>(read_i) >= chunk.reads.size()) continue;
        if (chunk.reads[read_i].is_skipped) continue;
        const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
        if (prof.start_var_idx < 0) continue;
        hts_pos_t ps = -1;
        for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
            if (!var_is_valid[vi]) continue;
            const CandidateVariant& var = chunk.candidates[vi];
            if (var.hap_to_cons_alle[1] != -1 && var.hap_to_cons_alle[2] != -1 &&
                var.hap_to_cons_alle[1] != var.hap_to_cons_alle[2]) {
                ps = var.phase_set;
            }
            if (ps != -1) break;
        }
        chunk.phase_sets[read_i] = ps;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Public entry point
// ════════════════════════════════════════════════════════════════════════════

void assign_hap_based_on_germline_het_vars_kmeans(BamChunk& chunk,
                                                   const Options& opts,
                                                   uint32_t flags) {
    const int n_cands = (int)chunk.candidates.size();
    std::vector<int> valid_var_idx;
    valid_var_idx.reserve(n_cands);
    std::vector<bool> var_is_valid(n_cands, false);
    for (int i = 0; i < n_cands; ++i) {
        if ((category_to_flag(chunk.candidates[i].counts.category) & flags) == 0) continue;
        valid_var_idx.push_back(i);
        var_is_valid[i] = true;
    }
    if (valid_var_idx.empty()) return;

    const bool is_ont = opts.is_ont();
    const size_t n_reads = chunk.reads.size();

    chunk.haps.assign(n_reads, 0);
    chunk.phase_sets.assign(n_reads, -1);

    read_init_hap_phase_set(chunk);
    var_init_hap_profile_cons_allele(is_ont, chunk.candidates, valid_var_idx);

    // Phase 1: initial sweep from highest-confidence pivot variant outward.
    const int init_vi = select_init_var(chunk.candidates, valid_var_idx);
    if (init_vi != -1) {
        const int nv = (int)valid_var_idx.size();

        // Build sweep order: [pivot, pivot-1, ..., 0, pivot+1, ..., nv-1].
        std::vector<int> var_ii(nv);
        var_ii[0] = init_vi;
        for (int vi = init_vi - 1; vi >= 0; --vi) var_ii[init_vi - vi] = vi;
        for (int vi = init_vi + 1; vi < nv; ++vi) var_ii[vi] = vi;

        int64_t* ovlp_b = nullptr;
        int64_t max_b = 0;
        cgranges_t* cr = chunk.read_var_cr.get();

        for (int idx = 0; idx < nv; ++idx) {
            const int vi = valid_var_idx[var_ii[idx]];
            const uint32_t vflag = category_to_flag(chunk.candidates[vi].counts.category);
            // Hom variants are not used in Phase 1.
            if (vflag == kCandCleanHom || vflag == kCandNoisyCandHom) continue;

            const int64_t ovlp_n = cr_overlap(cr, "cr", vi, vi + 1, &ovlp_b, &max_b);
            for (int64_t oi = 0; oi < ovlp_n; ++oi) {
                const int read_i = (int)cr_label(cr, ovlp_b[oi]);
                if (chunk.reads[read_i].is_skipped || chunk.haps[read_i] != 0) continue;
                int hap = init_assign_read_hap(chunk, read_i, flags);
                if (hap == -1) hap = 1; // no informative vars yet — seed new phase set as hap1
                chunk.haps[read_i] = hap;
                update_var_hap_profile_cons_alle(chunk, is_ont, read_i, hap, flags);
            }
        }
        free(ovlp_b);
    }

    // Phase 2: iterative k-means (up to 10 rounds, stop on convergence).
    for (int iter = 0; iter < 10; ++iter) {
        const int c1 = iter_update_var_hap_cons_phase_set(chunk, valid_var_idx, opts);
        const int c2 = iter_update_var_hap_to_cons_alle(chunk, is_ont, valid_var_idx, flags);
        if (c1 == 0 && c2 == 0) break;
    }

    // Phase 3: finalize per-read phase_sets.
    update_read_phase_set(chunk, var_is_valid);

    // Phase 4: fill hap_alt / hap_ref from finalized hap_to_cons_alle.
    // Mirrors longcallD make_variants() allele-resolution logic (collect_var.c):
    //   both -1 -> fall back to hap_to_cons_alle[0] (hom_idx);
    //   one -1  -> treat as ref (0), matching LONGCALLD_REF_ALLELE;
    //   any non-zero allele index is treated as ALT for GT projection.
    for (auto& var : chunk.candidates) {
        int c1 = var.hap_to_cons_alle[1];
        int c2 = var.hap_to_cons_alle[2];
        if (c1 == -1 && c2 == -1) {
            c1 = c2 = var.hap_to_cons_alle[0]; // hom_idx fallback
        }
        if (c1 == -1) c1 = 0; // unknown hap → ref
        if (c2 == -1) c2 = 0;
        const bool h1_alt = (c1 != 0);
        const bool h2_alt = (c2 != 0);
        if (h1_alt && h2_alt) {
            var.hap_alt = 3; var.hap_ref = 0;
        } else if (h1_alt && !h2_alt) {
            var.hap_alt = 1; var.hap_ref = 2;
        } else if (!h1_alt && h2_alt) {
            var.hap_alt = 2; var.hap_ref = 1;
        }
        // both ref/unresolved: leave hap_alt/hap_ref as 0/0

        if (opts.verbose >= 2 && !opts.debug_site.empty()) {
            hts_pos_t dbg = 0;
            if (parse_debug_site_pos(opts.debug_site, dbg)) {
                hts_pos_t vcf_pos = var.key.pos;
                if (var.key.type == VariantType::Insertion || var.key.type == VariantType::Deletion) {
                    vcf_pos = std::max<hts_pos_t>(1, var.key.pos - 1);
                }
                if (vcf_pos == dbg) {
                    const char t = (var.key.type == VariantType::Snp ? 'X'
                                 : (var.key.type == VariantType::Insertion ? 'I' : 'D'));
                    std::fprintf(stderr,
                                 "DebugPhase\tpos=%" PRId64 "\tkey_pos=%" PRId64 "\ttype=%c\talt=%s\tc1=%d\tc2=%d\thap_alt=%d\thap_ref=%d\thp=%d\n",
                                 static_cast<int64_t>(vcf_pos),
                                 static_cast<int64_t>(var.key.pos),
                                 t,
                                 var.key.alt.c_str(),
                                 var.hap_to_cons_alle[1],
                                 var.hap_to_cons_alle[2],
                                 var.hap_alt,
                                 var.hap_ref,
                                 var.is_homopolymer_indel ? 1 : 0);
                }
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Chunk-boundary stitching (mirrors longcallD flip_variant_hap / stitch_var_main)
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Stitch one pair of adjacent chunks (strict longcallD flip_variant_hap).
 *
 * Mirrors longcallD `flip_variant_hap` + `update_chunk_var_hap_phase_set1`
 * + `update_chunk_read_hap_phase_set1` from collect_var.c lines 1593–1695.
 */
static void apply_chunk_flip_and_merge(BamChunk& cur, bool do_flip, hts_pos_t max_pre_ps, hts_pos_t min_cur_ps) {
    if (do_flip && min_cur_ps != INT64_MAX) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set != min_cur_ps) continue;
            std::swap(v.hap_to_cons_alle[1], v.hap_to_cons_alle[2]);
        }
        for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
            if (cur.reads[read_i].is_skipped || cur.haps[read_i] == 0) continue;
            if (cur.phase_sets[read_i] == min_cur_ps) cur.haps[read_i] = 3 - cur.haps[read_i];
        }
    }
    if (max_pre_ps != -1 && min_cur_ps != INT64_MAX) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set == -1) continue;
            if (v.phase_set == min_cur_ps) v.phase_set = max_pre_ps;
        }
        for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
            if (cur.phase_sets[read_i] == -1) continue;
            if (cur.phase_sets[read_i] == min_cur_ps) cur.phase_sets[read_i] = max_pre_ps;
        }
    }
}

static void stitch_phase_blocks_within_chunk(BamChunk& chunk, const PgbamSidecarData& sidecar) {
    std::vector<hts_pos_t> psets;
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped || chunk.haps[i] == 0 || chunk.phase_sets[i] < 0) continue;
        psets.push_back(chunk.phase_sets[i]);
    }
    std::sort(psets.begin(), psets.end());
    psets.erase(std::unique(psets.begin(), psets.end()), psets.end());
    if (psets.size() < 2) return;

    std::vector<std::vector<uint64_t>> read_threads;
    build_read_thread_cache(chunk, sidecar, read_threads);
    std::vector<PhaseBlockThreadState> states =
        build_phase_block_thread_states(chunk, psets, read_threads);

    for (size_t i = 1; i < psets.size(); ++i) {
        const hts_pos_t left_ps  = psets[i - 1];
        const hts_pos_t right_ps = psets[i];

        bool do_flip = false;
        if (!decide_phase_block_concordance(states[i - 1], states[i], do_flip)) continue;
        apply_chunk_flip_and_merge(chunk, do_flip, left_ps, right_ps);
        states[i] = merge_phase_block_thread_states(states[i - 1], states[i], do_flip);
        psets[i] = left_ps;  // propagate merged identity so next iteration uses the correct left anchor
    }
}

static bool flip_chunk_hap(BamChunk& pre,
                           BamChunk& cur,
                           const Options* opts,
                           const PgbamSidecarData* pgbam_sidecar) {
    if (pre.region.tid != cur.region.tid) return false;
    int n_cur_ovlp_reads = 0;
    int n_pre_ovlp_reads = 0;
    const size_t n_bams = cur.up_ovlp_read_i.size();
    for (size_t bi = 0; bi < n_bams; ++bi) {
        n_cur_ovlp_reads += static_cast<int>(cur.up_ovlp_read_i[bi].size());
        if (bi < pre.down_ovlp_read_i.size()) {
            n_pre_ovlp_reads += static_cast<int>(pre.down_ovlp_read_i[bi].size());
        }
    }
    if (n_cur_ovlp_reads != n_pre_ovlp_reads) {
        throw std::runtime_error("overlap read count mismatch between adjacent chunks");
    }
    int flip_score = 0;
    hts_pos_t max_pre_ps = -1;
    hts_pos_t min_cur_ps = INT64_MAX;

    if (n_cur_ovlp_reads > 0) {
        if (pre.candidates.empty() || cur.candidates.empty()) return false;
        for (size_t bi = 0; bi < n_bams; ++bi) {
            const std::vector<int>& cur_list = cur.up_ovlp_read_i[bi];
            const std::vector<int>& pre_list = pre.down_ovlp_read_i[bi];
            for (size_t j = 0; j < cur_list.size(); ++j) {
                if (j >= pre_list.size()) {
                    throw std::runtime_error("overlap read pairing mismatch between adjacent chunks");
                }
                const int cur_read_i = cur_list[j];
                const int pre_read_i = pre_list[j];
                if (pre_read_i < 0 || static_cast<size_t>(pre_read_i) >= pre.reads.size() ||
                    cur_read_i < 0 || static_cast<size_t>(cur_read_i) >= cur.reads.size()) {
                    throw std::runtime_error("overlap read index out of bounds during chunk stitching");
                }
                if (pre.reads[static_cast<size_t>(pre_read_i)].is_skipped || pre.haps[static_cast<size_t>(pre_read_i)] == 0 ||
                    cur.reads[static_cast<size_t>(cur_read_i)].is_skipped || cur.haps[static_cast<size_t>(cur_read_i)] == 0) continue;
                const int pre_hap = pre.haps[static_cast<size_t>(pre_read_i)];
                const hts_pos_t pre_ps = pre.phase_sets[static_cast<size_t>(pre_read_i)];
                const int cur_hap = cur.haps[static_cast<size_t>(cur_read_i)];
                const hts_pos_t cur_ps = cur.phase_sets[static_cast<size_t>(cur_read_i)];
                if (pre_hap == cur_hap) --flip_score;
                else ++flip_score;
                if (max_pre_ps < pre_ps) max_pre_ps = pre_ps;
                if (min_cur_ps > cur_ps) min_cur_ps = cur_ps;
            }
        }
        if (flip_score != 0) {
            apply_chunk_flip_and_merge(cur, flip_score > 0, max_pre_ps, min_cur_ps);
            return true;
        }
    }

    if (opts == nullptr || pgbam_sidecar == nullptr || opts->pgbam_file.empty()) return false;

    max_pre_ps = -1;
    min_cur_ps = INT64_MAX;
    for (size_t i = 0; i < pre.reads.size(); ++i) {
        if (pre.reads[i].is_skipped || pre.haps[i] == 0 || pre.phase_sets[i] < 0) continue;
        if (pre.phase_sets[i] > max_pre_ps) max_pre_ps = pre.phase_sets[i];
    }
    for (size_t i = 0; i < cur.reads.size(); ++i) {
        if (cur.reads[i].is_skipped || cur.haps[i] == 0 || cur.phase_sets[i] < 0) continue;
        if (cur.phase_sets[i] < min_cur_ps) min_cur_ps = cur.phase_sets[i];
    }
    if (max_pre_ps < 0 || min_cur_ps == INT64_MAX) return false;

    const PhaseBlockThreadState pre_state =
        collect_phase_block_thread_state(pre, max_pre_ps, *pgbam_sidecar);
    const PhaseBlockThreadState cur_state =
        collect_phase_block_thread_state(cur, min_cur_ps, *pgbam_sidecar);
    bool do_flip = false;
    if (!decide_phase_block_concordance(pre_state, cur_state, do_flip)) return false;
    apply_chunk_flip_and_merge(cur, do_flip, max_pre_ps, min_cur_ps);
    return true;
}

void stitch_chunk_haps(std::vector<BamChunk>& chunks,
                       const Options* opts,
                       const PgbamSidecarData* pgbam_sidecar) {
    if (pgbam_sidecar != nullptr && opts != nullptr && !opts->pgbam_file.empty()) {
        for (BamChunk& chunk : chunks)
            stitch_phase_blocks_within_chunk(chunk, *pgbam_sidecar);
    }
    for (size_t ii = 1; ii < chunks.size(); ++ii)
        flip_chunk_hap(chunks[ii - 1], chunks[ii], opts, pgbam_sidecar);
}

} // namespace pgphase_collect
