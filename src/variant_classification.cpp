#include "variant_classification.hpp"

#include "noisy_regions.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace pgphase_collect {

/** longcallD `log_hypergeometric` / `fisher_exact_test` (math_utils.c), using `lgamma` only. */
static double log_hypergeom_pg(int a, int b, int c, int d) {
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N = n1 + n2;
    if (N <= 0) {
        return -std::numeric_limits<double>::infinity();
    }
    if (n1 > n2) {
        return log_hypergeom_pg(c, d, a, b);
    }
    if (m1 > m2) {
        return log_hypergeom_pg(b, a, d, c);
    }
    return std::lgamma(static_cast<double>(n1 + 1)) + std::lgamma(static_cast<double>(n2 + 1)) +
           std::lgamma(static_cast<double>(m1 + 1)) + std::lgamma(static_cast<double>(m2 + 1)) -
           (std::lgamma(static_cast<double>(a + 1)) + std::lgamma(static_cast<double>(b + 1)) +
            std::lgamma(static_cast<double>(c + 1)) + std::lgamma(static_cast<double>(d + 1)) +
            std::lgamma(static_cast<double>(N + 1)));
}

static double fisher_exact_two_tail(int a, int b, int c, int d) {
    if (a + b + c + d <= 0) {
        return 1.0;
    }
    const double log_p_observed = log_hypergeom_pg(a, b, c, d);
    const double p_observed = std::exp(log_p_observed);
    double total_p = 0.0;
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    const int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    const int denom = a + b + c + d;
    const int mode_a =
        denom > 0 ? static_cast<int>((static_cast<double>(a + b) * static_cast<double>(a + c)) / static_cast<double>(denom)) : 0;

    for (int delta = 0; delta <= max_a - min_a; ++delta) {
        int current_a = mode_a + delta;
        if (current_a <= max_a) {
            int current_b = (a + b) - current_a;
            int current_c = (a + c) - current_a;
            int current_d = (b + d) - current_b;
            if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                const double log_p = log_hypergeom_pg(current_a, current_b, current_c, current_d);
                const double p = std::exp(log_p);
                if (p <= p_observed + DBL_EPSILON) {
                    total_p += p;
                }
            }
        }
        if (delta > 0) {
            current_a = mode_a - delta;
            if (current_a >= min_a) {
                int current_b = (a + b) - current_a;
                int current_c = (a + c) - current_a;
                int current_d = (b + d) - current_b;
                if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                    const double log_p = log_hypergeom_pg(current_a, current_b, current_c, current_d);
                    const double p = std::exp(log_p);
                    if (p <= p_observed + DBL_EPSILON) {
                        total_p += p;
                    }
                }
            }
        }
    }
    return total_p;
}

static int nt4_from_ref_char(char ch) {
    switch (std::toupper(static_cast<unsigned char>(ch))) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4;
    }
}

static int ref_nt4_at(const std::string& ref_seq, hts_pos_t ref_beg, hts_pos_t abs_pos) {
    if (ref_seq.empty() || abs_pos < ref_beg) return 4;
    const hts_pos_t rel = abs_pos - ref_beg;
    if (rel < 0 || static_cast<size_t>(rel) >= ref_seq.size()) return 4;
    return nt4_from_ref_char(ref_seq[static_cast<size_t>(rel)]);
}

/** longcallD `var_is_homopolymer` (collect_var.c). */
static bool var_is_homopolymer_pg(const VariantKey& var,
                                  const std::string& ref_seq,
                                  hts_pos_t ref_beg,
                                  hts_pos_t ref_end,
                                  int xid) {
    if (ref_seq.empty()) return false;
    hts_pos_t start_pos = 0;
    hts_pos_t end_pos = 0;
    if (var.type == VariantType::Snp) {
        start_pos = var.pos - 1;
        end_pos = var.pos + 1;
    } else if (var.type == VariantType::Insertion) {
        const int alt_len = static_cast<int>(var.alt.size());
        if (alt_len > xid) return false;
        start_pos = var.pos - 1;
        end_pos = var.pos;
    } else {
        if (var.ref_len > xid) return false;
        start_pos = var.pos + var.ref_len - 1;
        end_pos = var.pos;
    }
    (void)ref_end;
    int is_homopolymer = 1;
    constexpr int max_unit_len = 6;
    constexpr int n_check_copy_num = 3;
    uint8_t ref_bases[6];
    for (int i = 0; i < 6; ++i) {
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, end_pos + i));
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, end_pos + i * rep_unit_len + j) != ref_bases[j]) {
                    is_homopolymer = 0;
                    break;
                }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) break;
    }
    if (is_homopolymer) return true;
    for (int i = 0; i < 6; ++i) {
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, start_pos - i));
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, start_pos - i * rep_unit_len - j) != ref_bases[j]) {
                    is_homopolymer = 0;
                    break;
                }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) break;
    }
    return is_homopolymer != 0;
}

/** longcallD `var_is_repeat_region` (collect_var.c), short indels only. */
static bool var_is_repeat_region_pg(const VariantKey& var,
                                    const std::string& ref_seq,
                                    hts_pos_t ref_beg,
                                    hts_pos_t ref_end,
                                    int xid) {
    if (ref_seq.empty()) return false;
    const hts_pos_t pos = var.pos;
    if (var.type == VariantType::Deletion) {
        const int del_len = var.ref_len;
        if (del_len > xid) return false;
        const int len = del_len * 3;
        if (pos < ref_beg || pos + del_len + len > ref_end) return false;
        const size_t off = static_cast<size_t>(pos - ref_beg);
        const size_t off2 = static_cast<size_t>(pos + del_len - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size() || off2 + static_cast<size_t>(len) > ref_seq.size()) {
            return false;
        }
        return std::memcmp(ref_seq.data() + off, ref_seq.data() + off2, static_cast<size_t>(len)) == 0;
    }
    if (var.type == VariantType::Insertion) {
        const int ins_len = static_cast<int>(var.alt.size());
        if (ins_len > xid) return false;
        const int len = ins_len * 3;
        if (pos < ref_beg || pos + len > ref_end) return false;
        const size_t off = static_cast<size_t>(pos - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size()) return false;
        std::string ref_b = ref_seq.substr(off, static_cast<size_t>(len));
        std::string alt_b = ref_b;
        for (int j = ins_len; j < len; ++j) {
            alt_b[static_cast<size_t>(j)] = alt_b[static_cast<size_t>(j - ins_len)];
        }
        for (int j = 0; j < ins_len; ++j) {
            alt_b[static_cast<size_t>(j)] = var.alt[static_cast<size_t>(j)];
        }
        return ref_b == alt_b;
    }
    return false;
}

static double noisy_reads_ratio_in_span(const BamChunk& chunk, hts_pos_t var_start, hts_pos_t var_end) {
    int total = 0;
    int noisy = 0;
    for (const ReadRecord& r : chunk.reads) {
        if (r.is_skipped) continue;
        if (r.end < var_start || r.beg > var_end) continue;
        ++total;
        bool hit = false;
        for (const Interval& iv : r.noisy_regions) {
            if (iv.end < var_start || iv.beg > var_end) continue;
            hit = true;
            break;
        }
        if (hit) ++noisy;
    }
    return total > 0 ? static_cast<double>(noisy) / static_cast<double>(total) : 0.0;
}

/** longcallD `cr_add_var_cr` (collect_var.c). */
static void cr_add_var_to_noisy_cr(cgranges_t* var_cr,
                                   cgranges_t* low_comp_cr,
                                   const VariantKey& var,
                                   const BamChunk& chunk,
                                   bool check_noisy_reads_ratio,
                                   const Options& opts) {
    hts_pos_t var_start = 0;
    hts_pos_t var_end = 0;
    variant_genomic_span(var, var_start, var_end);
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        int64_t* low_comp_b = nullptr;
        int64_t max_low_comp_n = 0;
        const int64_t low_comp_n = cr_overlap(
            low_comp_cr,
            "cr",
            static_cast<int32_t>(var_start - 1),
            static_cast<int32_t>(var_end),
            &low_comp_b,
            &max_low_comp_n);
        for (int64_t j = 0; j < low_comp_n; ++j) {
            const int32_t start = cr_start(low_comp_cr, low_comp_b[j]) + 1;
            const int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
            if (start < var_start) var_start = start;
            if (end > var_end) var_end = end;
        }
        std::free(low_comp_b);
    }
    if (!check_noisy_reads_ratio || noisy_reads_ratio_in_span(chunk, var_start, var_end) >= opts.min_af) {
        cr_add(var_cr,
               "cr",
               static_cast<int32_t>(var_start - 1),
               static_cast<int32_t>(var_end),
               1);
    }
}

static bool has_alt_strand_bias(const VariantCounts& counts) {
    const int major = std::max(counts.forward_alt, counts.reverse_alt);
    const int minor = std::min(counts.forward_alt, counts.reverse_alt);
    return counts.alt_cov >= 4 && minor == 0 && major >= 4;
}

/** longcallD `classify_var_cate` first stage (collect_var.c). */
static VariantCategory classify_variant_initial(const VariantKey& key,
                                                VariantCounts& counts,
                                                const std::string& ref_slice,
                                                hts_pos_t ref_beg,
                                                hts_pos_t ref_end,
                                                const Options& opts) {
    const int depth_with_low_quality = counts.total_cov + counts.low_qual_cov;
    counts.allele_fraction =
        counts.total_cov == 0 ? 0.0 : static_cast<double>(counts.alt_cov) / static_cast<double>(counts.total_cov);

    if (depth_with_low_quality < opts.min_depth) {
        return VariantCategory::LowCoverage;
    }
    if (counts.alt_cov < opts.min_alt_depth) {
        return VariantCategory::LowCoverage;
    }
    if (opts.is_ont) {
        const int fa = counts.forward_alt;
        const int ra = counts.reverse_alt;
        const int expected = (fa + ra) / 2;
        if (expected > 0) {
            const double p = fisher_exact_two_tail(fa, ra, expected, expected);
            if (p < opts.strand_bias_pval) {
                return VariantCategory::StrandBias;
            }
        }
    } else if (has_alt_strand_bias(counts)) {
        return VariantCategory::StrandBias;
    }
    if (counts.allele_fraction < opts.min_af) {
        return VariantCategory::LowAlleleFraction;
    }
    if (counts.allele_fraction > opts.max_af) {
        return VariantCategory::CleanHom;
    }
    if ((key.type == VariantType::Insertion || key.type == VariantType::Deletion) &&
        (var_is_homopolymer_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps) ||
         var_is_repeat_region_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps))) {
        return VariantCategory::RepeatHetIndel;
    }
    return key.type == VariantType::Snp ? VariantCategory::CleanHetSnp : VariantCategory::CleanHetIndel;
}

/** longcallD `classify_cand_vars` core: noisy overlap, REP/dense loci → `noisy_var_cr`, `cr_merge2`. */
void classify_cand_vars_pgphase(BamChunk& chunk, const Options& opts) {
    CandidateTable& variants = chunk.candidates;
    if (variants.empty()) return;

    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(chunk.low_complexity_regions));
    cgranges_t* low_comp = low_own.cr;
    if (low_comp != nullptr) {
        cr_index(low_comp);
    }

    std::vector<VariantCategory> cats;
    cats.reserve(variants.size());
    for (auto& cv : variants) {
        cats.push_back(classify_variant_initial(
            cv.key, cv.counts, chunk.ref_seq, chunk.ref_beg, chunk.ref_end, opts));
    }

    cgranges_t* var_pos_cr = cr_init();
    for (size_t i = 0; i < variants.size(); ++i) {
        const VariantCategory c = cats[i];
        if (c == VariantCategory::LowCoverage) continue;
        if (opts.is_ont && c == VariantCategory::StrandBias) continue;
        const VariantKey& k = variants[i].key;
        if (k.type == VariantType::Insertion) {
            cr_add(var_pos_cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), 1);
        } else {
            cr_add(var_pos_cr,
                   "cr",
                   static_cast<int32_t>(k.pos - 1),
                   static_cast<int32_t>(k.pos + k.ref_len - 1),
                   1);
        }
    }
    cr_index(var_pos_cr);

    cgranges_t* chunk_noisy = intervals_to_cr(chunk.noisy_regions);
    if (chunk_noisy == nullptr) {
        chunk_noisy = cr_init();
    } else {
        cr_index(chunk_noisy);
    }

    cgranges_t* noisy_var_cr = cr_init();
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;

    const hts_pos_t reg_beg = chunk.region.beg;
    const hts_pos_t reg_end = chunk.region.end;

    for (size_t i = 0; i < variants.size(); ++i) {
        VariantCategory c = cats[i];
        if (c == VariantCategory::StrandBias) continue;

        if (!opts.is_ont && chunk_noisy->n_r > 0) {
            int64_t n = 0;
            if (variants[i].key.type == VariantType::Insertion) {
                n = cr_overlap(chunk_noisy,
                               "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos),
                               &ovlp_b,
                               &max_b);
            } else {
                n = cr_overlap(chunk_noisy,
                               "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                               &ovlp_b,
                               &max_b);
            }
            if (n > 0) {
                cats[i] = VariantCategory::NonVariant;
                continue;
            }
        }

        if (c == VariantCategory::LowCoverage) continue;

        if (c == VariantCategory::RepeatHetIndel) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, false, opts);
            }
            continue;
        }

        int64_t n_ov = 0;
        if (variants[i].key.type == VariantType::Insertion) {
            n_ov = cr_overlap(var_pos_cr,
                              "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos),
                              &ovlp_b,
                              &max_b);
        } else {
            n_ov = cr_overlap(var_pos_cr,
                              "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                              &ovlp_b,
                              &max_b);
        }
        if (n_ov > 1) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, true, opts);
            }
        }

        if (c == VariantCategory::LowAlleleFraction) {
            cats[i] = VariantCategory::LowCoverage;
        }
    }

    std::free(ovlp_b);
    ovlp_b = nullptr;
    cr_destroy(var_pos_cr);

    if (noisy_var_cr->n_r > 0) {
        cr_index(noisy_var_cr);
        cgranges_t* merged =
            cr_merge2(chunk_noisy, noisy_var_cr, -1, opts.noisy_reg_merge_dis, opts.min_sv_len);
        cr_destroy(chunk_noisy);
        cr_destroy(noisy_var_cr);
        chunk_noisy = merged;
    } else {
        cr_destroy(noisy_var_cr);
    }

    intervals_from_cr(chunk_noisy, chunk.noisy_regions);
    cr_destroy(chunk_noisy);

    for (size_t i = 0; i < variants.size(); ++i) {
        variants[i].counts.category = cats[i];
    }
}

static void resolve_simple_noisy_candidates(CandidateTable& variants) {
    for (CandidateVariant& candidate : variants) {
        if (candidate.counts.category != VariantCategory::NoisyCandidate &&
            candidate.counts.category != VariantCategory::RepeatHetIndel) {
            continue;
        }
        if (candidate.key.alt.size() >= static_cast<size_t>(kMinSvLen) || candidate.key.ref_len >= kMinSvLen) {
            candidate.counts.category = VariantCategory::NoisyResolved;
        }
    }
}

void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    (void)header;
    classify_cand_vars_pgphase(chunk, opts);
    resolve_simple_noisy_candidates(chunk.candidates);
}

} // namespace pgphase_collect
