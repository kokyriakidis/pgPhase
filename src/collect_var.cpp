#include "collect_var.hpp"

#include "sdust.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Variant site comparison
// ════════════════════════════════════════════════════════════════════════════

/** 1-based alt length; matches longcallD var_site_t (SNP/INS: sequence; DEL: 0). */
static int var_site_alt_len(const VariantKey& v) {
    if (v.type == VariantType::Deletion) return 0;
    return static_cast<int>(v.alt.size());
}

/**
 * longcallD `exact_comp_var_site`: total order for qsort of var_site_t from digars.
 * Return <0 if a<b, 0 if equal, >0 if a>b.
 */
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2) {
    if (var1->tid != var2->tid) return var1->tid < var2->tid ? -1 : 1;
    const hts_pos_t p1 = var1->type == VariantType::Snp ? var1->pos : var1->pos - 1;
    const hts_pos_t p2 = var2->type == VariantType::Snp ? var2->pos : var2->pos - 1;
    if (p1 < p2) return -1;
    if (p1 > p2) return 1;
    const int t1 = static_cast<int>(var1->type);
    const int t2 = static_cast<int>(var2->type);
    if (t1 < t2) return -1;
    if (t1 > t2) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    const int a1 = var_site_alt_len(*var1);
    const int a2 = var_site_alt_len(*var2);
    if (a1 < a2) return -1;
    if (a1 > a2) return 1;
    if (var1->type == VariantType::Snp || var1->type == VariantType::Insertion) {
        return std::memcmp(var1->alt.data(), var2->alt.data(), static_cast<size_t>(a1));
    }
    return 0;
}

/**
 * longcallD `exact_comp_var_site_ins`: same sort keys, but large insertions (>= min_sv_len)
 * merge when min(alt_len) >= 0.8 * max(alt_len); small insertions still require exact sequence.
 */
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len) {
    if (var1->tid != var2->tid) return var1->tid < var2->tid ? -1 : 1;
    const hts_pos_t p1 = var1->type == VariantType::Snp ? var1->pos : var1->pos - 1;
    const hts_pos_t p2 = var2->type == VariantType::Snp ? var2->pos : var2->pos - 1;
    if (p1 < p2) return -1;
    if (p1 > p2) return 1;
    const int t1 = static_cast<int>(var1->type);
    const int t2 = static_cast<int>(var2->type);
    if (t1 < t2) return -1;
    if (t1 > t2) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->type == VariantType::Snp) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < a2) return -1;
        if (a1 > a2) return 1;
        return std::memcmp(var1->alt.data(), var2->alt.data(), static_cast<size_t>(a1));
    }
    if (var1->type == VariantType::Insertion) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < min_sv_len) {
            if (a1 < a2) return -1;
            if (a1 > a2) return 1;
            return std::memcmp(var1->alt.data(), var2->alt.data(), static_cast<size_t>(a1));
        }
        const int min_len = a1 < a2 ? a1 : a2;
        const int max_len = a1 > a2 ? a1 : a2;
        if (static_cast<double>(min_len) >= static_cast<double>(max_len) * 0.8) return 0;
        return a1 - a2;
    }
    return 0;
}

/** STL comparator wrapping exact_comp_var_site; used for std::map / std::set keyed on VariantKey. */
bool VariantKeyLess::operator()(const VariantKey& lhs, const VariantKey& rhs) const {
    return exact_comp_var_site(&lhs, &rhs) < 0;
}

// ════════════════════════════════════════════════════════════════════════════
// Candidate site collection
// ════════════════════════════════════════════════════════════════════════════

/** Convert one digar op to a VariantKey suitable for the CandidateTable.
 *  SNP: ref_len=1; INS: ref_len=0 (insertion has no consumed reference bases); DEL: ref_len=op.len. */
VariantKey variant_key_from_digar(int tid, const DigarOp& op) {
    if (op.type == DigarType::Snp) {
        return VariantKey{tid, op.pos, VariantType::Snp, 1, op.alt};
    }
    if (op.type == DigarType::Insertion) {
        return VariantKey{tid, op.pos, VariantType::Insertion, 0, op.alt};
    }
    return VariantKey{tid, op.pos, VariantType::Deletion, op.len, ""};
}

/** longcallD `is_collectible_var_digar` (reg_{beg,end} == -1 disables that side). */
static bool is_collectible_var_digar(const DigarOp& digar, hts_pos_t reg_beg, hts_pos_t reg_end) {
    const hts_pos_t digar_pos = digar.pos;
    if (reg_beg != -1 && digar_pos < reg_beg) return false;
    if (reg_end != -1 && digar_pos > reg_end) return false;
    if (digar.low_quality) return false;
    return digar.type == DigarType::Snp || digar.type == DigarType::Insertion ||
           digar.type == DigarType::Deletion;
}

/** True for the three variant types that participate in the allele-counting sweep; clips/equal skip. */
static bool is_variant_digar_for_cand_sweep(const DigarOp& d) {
    return d.type == DigarType::Snp || d.type == DigarType::Insertion || d.type == DigarType::Deletion;
}

/** longcallD `get_var_site_start` (bam_utils.c) on sorted candidate `VariantKey`s. */
static int get_var_site_start(const CandidateTable& v, hts_pos_t start, int n_total) {
    hts_pos_t target = start > 0 ? start - 1 : start;
    int left = 0;
    int right = n_total;
    while (left < right) {
        const int mid = left + (right - left) / 2;
        const hts_pos_t mid_pos = v[static_cast<size_t>(mid)].key.sort_pos();
        if (mid_pos < target) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    while (left < n_total && v[static_cast<size_t>(left)].key.pos < start) {
        ++left;
    }
    return left;
}

/** longcallD `collect_all_cand_var_sites` merge pass: qsort order + `exact_comp_var_site_ins`. */
void collapse_fuzzy_large_insertions(CandidateTable& variants) {
    std::sort(variants.begin(), variants.end(), [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_var_site(&a.key, &b.key) < 0;
    });
    if (variants.empty()) return;
    size_t write_i = 1;
    for (size_t read_i = 1; read_i < variants.size(); ++read_i) {
        if (exact_comp_var_site_ins(
                &variants[write_i - 1].key, &variants[read_i].key, kLongcalldMinSvLen) == 0) {
            continue;
        }
        if (write_i != read_i) variants[write_i] = std::move(variants[read_i]);
        ++write_i;
    }
    variants.resize(write_i);
}

/** Gather one raw CandidateVariant per digar variant op, then deduplicate via
 *  collapse_fuzzy_large_insertions. Corresponds to the first half of longcallD
 *  `collect_all_cand_var_sites` (collect_var.c). */
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants) {
    variants.clear();
    size_t max_var_sites = 0;
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            ++max_var_sites;
        }
    }
    if (max_var_sites == 0) return;

    variants.reserve(max_var_sites);
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            variants.push_back(CandidateVariant{variant_key_from_digar(read.tid, digar), VariantCounts{}});
        }
    }
    collapse_fuzzy_large_insertions(variants);
}

// ════════════════════════════════════════════════════════════════════════════
// Allele counting
// ════════════════════════════════════════════════════════════════════════════

/** longcallD `get_digar_ave_qual` (bam_utils.c) — average BQ over bases supporting the call. */
static int get_digar_ave_qual(const DigarOp& d, const std::vector<uint8_t>& qual) {
    if (d.low_quality) return 0;
    if (d.qi < 0) return 0;
    if (qual.empty()) return 0;
    int q_start = 0;
    int q_end = 0;
    if (d.type == DigarType::Deletion) {
        if (d.qi == 0) {
            q_start = 0;
            q_end = 0;
        } else {
            q_start = d.qi - 1;
            q_end = d.qi;
        }
    } else {
        q_start = d.qi;
        q_end = d.qi + d.len - 1;
    }
    if (q_start < 0) return 0;
    if (q_end >= static_cast<int>(qual.size())) return 0;
    int sum = 0;
    for (int i = q_start; i <= q_end; ++i) sum += static_cast<int>(qual[i]);
    return sum / (q_end - q_start + 1);
}

/** Increment the appropriate coverage counters for one read observation at a candidate site.
 *  Palindromic ONT reads count toward total/ref/alt depth but NOT toward strand counts,
 *  because the fold-back artefact makes strand information unreliable for those reads. */
static void add_coverage(VariantCounts& counts, bool reverse, bool alt, bool low_quality,
                         bool is_ont_palindrome) {
    if (low_quality) {
        counts.low_qual_cov++;
        return;
    }
    counts.total_cov++;
    if (alt) {
        counts.alt_cov++;
        if (!is_ont_palindrome) reverse ? counts.reverse_alt++ : counts.forward_alt++;
    } else {
        counts.ref_cov++;
        if (!is_ont_palindrome) reverse ? counts.reverse_ref++ : counts.forward_ref++;
    }
}

/**
 * longcallD `init_cand_vars_based_on_sites` + `update_cand_vars_from_digar` (bam_utils.c):
 * one pass over each read's digars, merged with the sorted `CandidateTable` using
 * `exact_comp_var_site_ins`; alt BQ = max(digar flag, mean qual < min_bq) like
 * is_low_qual || ave_qual < opt->min_bq; trailing sites with pos <= read end get
 * ref unless past pos_end.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq) {
    const int n = static_cast<int>(variants.size());
    if (n == 0) return;

    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        if (read.tid < 0) continue;
        const hts_pos_t pos_start = read.beg;
        const hts_pos_t pos_end = read.end;
        int site_i = get_var_site_start(variants, pos_start, n);
        size_t digar_i = 0;
        const std::vector<DigarOp>& dig = read.digars;

        while (site_i < n && digar_i < dig.size()) {
            if (const int stid = variants[static_cast<size_t>(site_i)].key.tid; stid != read.tid) {
                if (stid < read.tid) {
                    ++site_i;
                } else {
                    break;
                }
                continue;
            }
            if (!is_variant_digar_for_cand_sweep(dig[digar_i])) {
                ++digar_i;
                continue;
            }
            const DigarOp& d = dig[digar_i];
            VariantKey dkey = variant_key_from_digar(read.tid, d);
            const int ave = get_digar_ave_qual(d, read.qual);
            const int ret = exact_comp_var_site_ins(
                &variants[static_cast<size_t>(site_i)].key, &dkey, kLongcalldMinSvLen);
            if (ret < 0) {
                VariantCounts& c = variants[static_cast<size_t>(site_i)].counts;
                add_coverage(c, read.reverse, false, false, read.is_ont_palindrome);
                if (read_support_out != nullptr && chunk_region != nullptr) {
                    const auto& k = variants[static_cast<size_t>(site_i)].key;
                    read_support_out->push_back(ReadSupportRow{k.tid,
                                                               k.pos,
                                                               k.type,
                                                               k.ref_len,
                                                               k.alt,
                                                               read.qname,
                                                               0,
                                                               0,
                                                               read.reverse,
                                                               read.mapq,
                                                               chunk_region->beg,
                                                               chunk_region->end});
                }
                ++site_i;
            } else if (ret == 0) {
                const bool lq = d.low_quality || ave < min_bq;
                add_coverage(
                    variants[static_cast<size_t>(site_i)].counts, read.reverse, true, lq,
                    read.is_ont_palindrome);
                if (read_support_out != nullptr && chunk_region != nullptr) {
                    const auto& k = variants[static_cast<size_t>(site_i)].key;
                    read_support_out->push_back(
                        ReadSupportRow{k.tid,
                                       k.pos,
                                       k.type,
                                       k.ref_len,
                                       k.alt,
                                       read.qname,
                                       1,
                                       lq ? 1 : 0,
                                       read.reverse,
                                       read.mapq,
                                       chunk_region->beg,
                                       chunk_region->end});
                }
                ++site_i;
            } else {
                ++digar_i;
            }
        }

        for (; site_i < n; ++site_i) {
            if (const int stid = variants[static_cast<size_t>(site_i)].key.tid; stid != read.tid) {
                if (stid < read.tid) continue;
                break;
            }
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) break;
            add_coverage(variants[static_cast<size_t>(site_i)].counts, read.reverse, false, false,
                         read.is_ont_palindrome);
            if (read_support_out != nullptr && chunk_region != nullptr) {
                const auto& k = variants[static_cast<size_t>(site_i)].key;
                read_support_out->push_back(ReadSupportRow{k.tid,
                                                           k.pos,
                                                           k.type,
                                                           k.ref_len,
                                                           k.alt,
                                                           read.qname,
                                                           0,
                                                           0,
                                                           read.reverse,
                                                           read.mapq,
                                                           chunk_region->beg,
                                                           chunk_region->end});
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Interval utilities
// ════════════════════════════════════════════════════════════════════════════

/** Sort intervals by start, then merge overlapping or adjacent ones in place.
 *  The label of the merged interval is the max of the inputs (used to preserve noisy-region size). */
void merge_intervals(std::vector<Interval>& intervals) {
    if (intervals.empty()) return;
    std::sort(intervals.begin(), intervals.end(), [](const Interval& lhs, const Interval& rhs) {
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        return lhs.end < rhs.end;
    });
    size_t write_i = 0;
    for (size_t read_i = 1; read_i < intervals.size(); ++read_i) {
        if (intervals[read_i].beg <= intervals[write_i].end + 1) {
            intervals[write_i].end = std::max(intervals[write_i].end, intervals[read_i].end);
            intervals[write_i].label = std::max(intervals[write_i].label, intervals[read_i].label);
        } else {
            intervals[++write_i] = intervals[read_i];
        }
    }
    intervals.resize(write_i + 1);
}

/** Convert a vector of Intervals to a cgranges tree (0-based half-open internally).
 *  Returns nullptr if the vector is empty or all intervals are degenerate (beg > end).
 *  The returned tree is NOT indexed; callers must call cr_index before overlap queries. */
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals) {
    if (intervals.empty()) return nullptr;
    cgranges_t* cr = cr_init();
    int added = 0;
    for (const Interval& iv : intervals) {
        if (iv.beg > iv.end) continue;
        const int32_t st = static_cast<int32_t>(iv.beg - 1);
        const int32_t en = static_cast<int32_t>(iv.end);
        const int32_t label = iv.label > 0 ? iv.label : 1;
        cr_add(cr, "cr", st, en, label);
        ++added;
    }
    if (added == 0) {
        cr_destroy(cr);
        return nullptr;
    }
    // Leave unindexed: callers run cr_index once before overlap/merge.
    return cr;
}

/** Extract all intervals stored under the "cr" contig back into a vector (inverse of intervals_to_cr). */
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out) {
    out.clear();
    if (cr == nullptr || cr->n_r == 0) return;
    const int32_t ctg = cr_get_ctg(cr, "cr");
    if (ctg < 0) return;
    const cr_ctg_t* c = &cr->ctg[ctg];
    for (int64_t i = c->off; i < c->off + c->n; ++i) {
        const hts_pos_t beg = static_cast<hts_pos_t>(cr_start(cr, i) + 1);
        const hts_pos_t end = static_cast<hts_pos_t>(cr_end(cr, i));
        out.push_back(Interval{beg, end, cr_label(cr, i)});
    }
}

/** Build an indexed cgranges tree from a single read's noisy_regions; used in pre_process_noisy_regs
 *  to test whether the read has any errors overlapping each chunk noisy region. */
cgranges_t* build_read_noisy_cr(const ReadRecord& read) {
    if (read.noisy_regions.empty()) return nullptr;
    cgranges_t* cr = intervals_to_cr(read.noisy_regions);
    if (cr != nullptr) cr_index(cr);
    return cr;
}

static bool read_has_noisy_overlap(const ReadRecord& read, hts_pos_t start, hts_pos_t end) {
    for (const Interval& iv : read.noisy_regions) {
        if (iv.end < start) continue;
        if (iv.beg > end) return false;
        return true;
    }
    return false;
}

/** Return the 1-based inclusive [var_start, var_end] reference span for a variant.
 *  For insertions, var_end < var_start (zero-length reference span) to match longcallD's convention. */
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end) {
    if (key.type == VariantType::Snp) {
        var_start = var_end = key.pos;
        return;
    }
    if (key.type == VariantType::Deletion) {
        var_start = key.pos;
        var_end = key.pos + key.ref_len - 1;
        return;
    }
    // Insertions have ref_len=0; longcallD uses pos..pos+ref_len-1.
    var_start = key.pos;
    var_end = key.pos - 1;
}

// ════════════════════════════════════════════════════════════════════════════
// Noisy region processing
// ════════════════════════════════════════════════════════════════════════════

/** Expand [start, end] to include any overlapping low-complexity intervals.
 *  Used by cr_extend_noisy_regs_with_low_comp to anchor noisy regions to sdust windows. */
static void low_comp_cr_start_end(cgranges_t* low_comp_cr,
                                  int32_t start, int32_t end,
                                  int32_t* new_start, int32_t* new_end) {
    *new_start = start;
    *new_end = end;
    if (low_comp_cr == nullptr || low_comp_cr->n_r == 0) return;
    int64_t* low_comp_b = nullptr;
    int64_t max_low_comp_n = 0;
    const int64_t low_comp_n =
        cr_overlap(low_comp_cr, "cr", start - 1, end, &low_comp_b, &max_low_comp_n);
    for (int64_t j = 0; j < low_comp_n; ++j) {
        const int32_t s = cr_start(low_comp_cr, low_comp_b[j]) + 1;
        const int32_t e = cr_end(low_comp_cr, low_comp_b[j]);
        if (s < *new_start) *new_start = s;
        if (e > *new_end) *new_end = e;
    }
    std::free(low_comp_b);
}

/**
 * longcallD `cr_extend_noisy_regs_with_low_comp` (collect_var.c): extend into
 * low-complexity intervals, then always `cr_merge(…, -1, merge_dis, min_sv_len)`.
 */
static cgranges_t* cr_extend_noisy_regs_with_low_comp(cgranges_t* chunk_noisy_regs,
                                                      cgranges_t* low_comp_cr,
                                                      int merge_dis,
                                                      int min_sv_len) {
    if (chunk_noisy_regs == nullptr) return nullptr;
    cgranges_t* cur = chunk_noisy_regs;
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        cgranges_t* new_noisy_regs = cr_init();
        for (int i = 0; i < cur->n_r; ++i) {
            const int32_t start = cr_start(cur, i) + 1;
            const int32_t end = cr_end(cur, i);
            int32_t new_start = start;
            int32_t new_end = end;
            low_comp_cr_start_end(low_comp_cr, start, end, &new_start, &new_end);
            cr_add(new_noisy_regs, "cr", new_start - 1, new_end, cr_label(cur, i));
        }
        cr_index(new_noisy_regs);
        cr_destroy(cur);
        cur = new_noisy_regs;
    }
    return cr_merge(cur, -1, merge_dis, min_sv_len);
}

/** Categories that should not anchor a noisy-region flank in collect_noisy_reg_start_end. */
static bool category_skipped_for_noisy_flank(VariantCategory c) {
    return c == VariantCategory::LowCoverage || c == VariantCategory::StrandBias ||
           c == VariantCategory::NonVariant;
}

/** longcallD `collect_noisy_reg_start_end` (collect_var.c): for each noisy region, find the nearest
 *  clean candidate variants on either side and extend the region boundary to include them plus
 *  kNoisyRegFlankLen bp. This widens noisy windows to capture the flanking het variants that define
 *  the haplotype context for MSA realignment. */
static void collect_noisy_reg_start_end_pgphase(const cgranges_t* chunk_noisy_regs,
                                                const CandidateTable& cand,
                                                const std::vector<VariantCategory>& cand_cate,
                                                std::vector<int32_t>& start_out,
                                                std::vector<int32_t>& end_out) {
    const int n_noisy = chunk_noisy_regs->n_r;
    start_out.assign(n_noisy, 0);
    end_out.assign(n_noisy, 0);
    if (cand.empty()) {
        for (int reg_i = 0; reg_i < n_noisy; ++reg_i) {
            const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
            const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
            start_out[static_cast<size_t>(reg_i)] = ori_reg_start - kNoisyRegFlankLen;
            end_out[static_cast<size_t>(reg_i)] = ori_reg_end + kNoisyRegFlankLen;
        }
        return;
    }
    std::vector<int> max_left_var_i(n_noisy, -1);
    std::vector<int> min_right_var_i(n_noisy, -1);

    int reg_i = 0;
    int var_i = 0;
    while (reg_i < n_noisy && var_i < static_cast<int>(cand.size())) {
        if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(var_i)])) {
            ++var_i;
            continue;
        }
        hts_pos_t var_start = 0;
        hts_pos_t var_end = 0;
        variant_genomic_span(cand[static_cast<size_t>(var_i)].key, var_start, var_end);
        const int32_t reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t reg_end = cr_end(chunk_noisy_regs, reg_i);
        if (static_cast<int32_t>(var_start) > reg_end) {
            if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
                min_right_var_i[static_cast<size_t>(reg_i)] = var_i;
            }
            ++reg_i;
        } else if (static_cast<int32_t>(var_end) < reg_start) {
            max_left_var_i[static_cast<size_t>(reg_i)] = var_i;
            ++var_i;
        } else {
            ++var_i;
        }
    }

    for (reg_i = 0; reg_i < n_noisy; ++reg_i) {
        if (max_left_var_i[static_cast<size_t>(reg_i)] == -1) {
            max_left_var_i[static_cast<size_t>(reg_i)] =
                std::min(std::max(0, static_cast<int>(cand.size()) - 1), 0);
        }
        if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
            min_right_var_i[static_cast<size_t>(reg_i)] =
                std::max(0, static_cast<int>(cand.size()) - 1);
        }
        const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
        int32_t cur_start = ori_reg_start - kNoisyRegFlankLen;
        for (int v = max_left_var_i[static_cast<size_t>(reg_i)]; v >= 0; --v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            variant_genomic_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(ve) < cur_start - 1) break;
            if (static_cast<int32_t>(vs) - kNoisyRegFlankLen < cur_start) {
                cur_start = static_cast<int32_t>(vs) - kNoisyRegFlankLen;
            }
        }
        start_out[static_cast<size_t>(reg_i)] = cur_start;
        int32_t cur_end = ori_reg_end + kNoisyRegFlankLen;
        for (int v = min_right_var_i[static_cast<size_t>(reg_i)]; v < static_cast<int>(cand.size()); ++v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            variant_genomic_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(vs) > cur_end + 1) break;
            if (static_cast<int32_t>(ve) + kNoisyRegFlankLen > cur_end) {
                cur_end = static_cast<int32_t>(ve) + kNoisyRegFlankLen;
            }
        }
        end_out[static_cast<size_t>(reg_i)] = cur_end;
    }
}

/** longcallD `pre_process_noisy_regs` (collect_var.c): merge raw per-read noisy intervals into
 *  chunk-level regions, extend into low-complexity windows, then filter out regions where fewer than
 *  min_alt_depth reads or a ratio below min_af are actually noisy (likely noise rather than true SVs). */
void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts) {
    if (chunk.noisy_regions.empty()) return;

    CrangesOwner noisy_own;
    noisy_own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (noisy_own.cr == nullptr) {
        chunk.noisy_regions.clear();
        return;
    }
    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(chunk.low_complexity_regions));

    cr_index(noisy_own.cr);
    if (low_own.cr != nullptr) cr_index(low_own.cr);

    const int merge_dis = opts.noisy_reg_merge_dis;
    const int msv = opts.min_sv_len;

    // longcallD: cr_extend_noisy_regs_with_low_comp then a second cr_merge.
    cgranges_t* noisy = noisy_own.release();
    noisy = cr_extend_noisy_regs_with_low_comp(noisy, low_own.cr, merge_dis, msv);
    noisy = cr_merge(noisy, -1, merge_dis, msv);
    noisy_own.adopt(noisy);

    cgranges_t* noisy_regs = noisy_own.cr;
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    std::vector<uint8_t> skip_noisy_reg(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_total(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_noisy(static_cast<size_t>(noisy_regs->n_r), 0);

    const std::vector<int>& read_order = chunk.ordered_read_ids;
    const auto visit_read = [&](const ReadRecord& read) {
        if (read.is_skipped) return;
        const int64_t beg = read.beg;
        const int64_t end = read.end;
        const int64_t ovlp_n = cr_overlap(
            noisy_regs, "cr", static_cast<int32_t>(beg - 1), static_cast<int32_t>(end), &ovlp_b, &max_b);
        for (int64_t ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            const int r_i = static_cast<int>(ovlp_b[ovlp_i]);
            noisy_reg_to_total[static_cast<size_t>(r_i)]++;
            const int32_t noisy_reg_start = cr_start(noisy_regs, r_i) + 1;
            const int32_t noisy_reg_end = cr_end(noisy_regs, r_i);
            if (read_has_noisy_overlap(read, noisy_reg_start, noisy_reg_end)) {
                noisy_reg_to_noisy[static_cast<size_t>(r_i)]++;
            }
        }
    };
    if (read_order.empty() && !chunk.reads.empty()) {
        for (const ReadRecord& r : chunk.reads) visit_read(r);
    } else {
        for (int ord : read_order) {
            if (ord < 0 || static_cast<size_t>(ord) >= chunk.reads.size()) continue;
            visit_read(chunk.reads[static_cast<size_t>(ord)]);
        }
    }
    std::free(ovlp_b);

    const int min_noisy_reg_reads = opts.min_alt_depth;
    const float min_noisy_reg_ratio = static_cast<float>(opts.min_af);
    const int min_overlap_reads = opts.min_noisy_reg_total_depth;
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        const int n_noisy = noisy_reg_to_noisy[static_cast<size_t>(i)];
        const int n_total = noisy_reg_to_total[static_cast<size_t>(i)];
        if ((min_overlap_reads > 0 && n_total < min_overlap_reads) ||
            n_noisy < min_noisy_reg_reads ||
            (n_total > 0 &&
             static_cast<float>(n_noisy) / static_cast<float>(n_total) < min_noisy_reg_ratio)) {
            skip_noisy_reg[static_cast<size_t>(i)] = 1;
            ++n_skipped;
        }
    }

    if (n_skipped > 0) {
        cgranges_t* new_noisy_regs = cr_init();
        for (int i = 0; i < noisy_regs->n_r; ++i) {
            if (skip_noisy_reg[static_cast<size_t>(i)]) continue;
            cr_add(new_noisy_regs, "cr", cr_start(noisy_regs, i), cr_end(noisy_regs, i),
                   cr_label(noisy_regs, i));
        }
        cr_index(new_noisy_regs);
        noisy_own.adopt(new_noisy_regs);
    }

    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

/** longcallD `post_process_noisy_regs` (collect_var.c): after classification, expand each noisy
 *  region boundary outward to the nearest flanking clean variant (plus kNoisyRegFlankLen), then
 *  re-merge. This ensures the MSA window in Step 2 contains all variants needed for phasing. */
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand) {
    if (chunk.noisy_regions.empty()) return;
    CrangesOwner noisy_own;
    noisy_own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (noisy_own.cr == nullptr) {
        chunk.noisy_regions.clear();
        return;
    }
    cr_index(noisy_own.cr);
    cgranges_t* chunk_noisy_regs = noisy_own.release();
    const int n_noisy_regs = chunk_noisy_regs->n_r;
    std::vector<int32_t> noisy_reg_start(static_cast<size_t>(n_noisy_regs));
    std::vector<int32_t> noisy_reg_end(static_cast<size_t>(n_noisy_regs));
    std::vector<VariantCategory> cand_cate;
    cand_cate.reserve(cand.size());
    for (const auto& cv : cand) cand_cate.push_back(cv.counts.category);
    collect_noisy_reg_start_end_pgphase(
        chunk_noisy_regs, cand, cand_cate, noisy_reg_start, noisy_reg_end);

    cgranges_t* noisy_regs = cr_init();
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        cr_add(noisy_regs, "cr",
               noisy_reg_start[static_cast<size_t>(reg_i)],
               noisy_reg_end[static_cast<size_t>(reg_i)],
               cr_label(chunk_noisy_regs, reg_i));
    }
    cr_index(noisy_regs);
    cr_destroy(chunk_noisy_regs);
    cgranges_t* merged = cr_merge(noisy_regs, 0, -1, -1);
    noisy_own.adopt(merged);
    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

/** After `post_process_noisy_regs_pgphase`, longcallD `cr_is_contained` sweep. */
void apply_noisy_containment_filter(BamChunk& chunk) {
    if (chunk.noisy_regions.empty()) return;
    CrangesOwner own;
    own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (own.cr == nullptr || own.cr->n_r == 0) return;
    cr_index(own.cr);
    int64_t* b = nullptr;
    int64_t m = 0;
    for (auto& cv : chunk.candidates) {
        if (cv.counts.category == VariantCategory::NonVariant) continue;
        const VariantKey& k = cv.key;
        int64_t n = 0;
        if (k.type == VariantType::Insertion) {
            n = cr_is_contained(
                own.cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), &b, &m);
        } else {
            n = cr_is_contained(own.cr, "cr",
                                static_cast<int32_t>(k.pos - 1),
                                static_cast<int32_t>(k.pos + k.ref_len),
                                &b, &m);
        }
        if (n > 0) cv.counts.category = VariantCategory::NonVariant;
    }
    std::free(b);
}

/** Run sdust on the chunk's reference slice to identify low-complexity regions (tandem repeats,
 *  homopolymers). These are used to anchor and extend noisy regions in pre_process_noisy_regs. */
void populate_low_complexity_intervals(BamChunk& chunk) {
    chunk.low_complexity_regions.clear();
    if (chunk.ref_seq.empty()) return;

    const hts_pos_t beg = std::max(chunk.region.beg, chunk.ref_beg);
    const hts_pos_t end = std::min(chunk.region.end, chunk.ref_end);
    if (beg > end) return;

    const size_t offset = static_cast<size_t>(beg - chunk.ref_beg);
    const int len = static_cast<int>(end - beg + 1);
    int n = 0;
    uint64_t* intervals = sdust(nullptr,
                                reinterpret_cast<const uint8_t*>(chunk.ref_seq.data() + offset),
                                len, kSdustThreshold, kSdustWindow, &n);
    for (int i = 0; i < n; ++i) {
        const hts_pos_t rel_beg = static_cast<hts_pos_t>(intervals[i] >> 32);
        const hts_pos_t rel_end = static_cast<hts_pos_t>(static_cast<uint32_t>(intervals[i]));
        if (rel_end <= rel_beg) continue;
        chunk.low_complexity_regions.push_back(
            Interval{beg + rel_beg, beg + rel_end - 1, static_cast<int>(rel_end - rel_beg)});
    }
    std::free(intervals);
}

/** Build chunk.ordered_read_ids (index list for noisy-region sweeps) and collect all per-read
 *  noisy_regions into chunk.noisy_regions, merging overlapping intervals. */
void populate_chunk_read_indexes(BamChunk& chunk) {
    chunk.ordered_read_ids.clear();
    chunk.noisy_regions.clear();
    chunk.ordered_read_ids.reserve(chunk.reads.size());
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        chunk.ordered_read_ids.push_back(static_cast<int>(read_i));
        if (read.is_skipped) continue;
        chunk.noisy_regions.insert(
            chunk.noisy_regions.end(), read.noisy_regions.begin(), read.noisy_regions.end());
    }
    merge_intervals(chunk.noisy_regions);
}

// ════════════════════════════════════════════════════════════════════════════
// Variant classification — math utilities
// ════════════════════════════════════════════════════════════════════════════

/** longcallD `log_hypergeometric` / `fisher_exact_test` (math_utils.c), using `lgamma` only. */
static double log_hypergeom_pg(int a, int b, int c, int d) {
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N = n1 + n2;
    if (N <= 0) return -std::numeric_limits<double>::infinity();
    if (n1 > n2) return log_hypergeom_pg(c, d, a, b);
    if (m1 > m2) return log_hypergeom_pg(b, a, d, c);
    return std::lgamma(static_cast<double>(n1 + 1)) + std::lgamma(static_cast<double>(n2 + 1)) +
           std::lgamma(static_cast<double>(m1 + 1)) + std::lgamma(static_cast<double>(m2 + 1)) -
           (std::lgamma(static_cast<double>(a + 1)) + std::lgamma(static_cast<double>(b + 1)) +
            std::lgamma(static_cast<double>(c + 1)) + std::lgamma(static_cast<double>(d + 1)) +
            std::lgamma(static_cast<double>(N + 1)));
}

static double fisher_exact_two_tail(int a, int b, int c, int d) {
    if (a + b + c + d <= 0) return 1.0;
    const double log_p_observed = log_hypergeom_pg(a, b, c, d);
    const double p_observed = std::exp(log_p_observed);
    double total_p = 0.0;
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    const int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    const int denom = a + b + c + d;
    const int mode_a =
        denom > 0
            ? static_cast<int>((static_cast<double>(a + b) * static_cast<double>(a + c)) /
                               static_cast<double>(denom))
            : 0;
    for (int delta = 0; delta <= max_a - min_a; ++delta) {
        int current_a = mode_a + delta;
        if (current_a <= max_a) {
            int current_b = (a + b) - current_a;
            int current_c = (a + c) - current_a;
            int current_d = (b + d) - current_b;
            if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                const double p = std::exp(log_hypergeom_pg(current_a, current_b, current_c, current_d));
                if (p <= p_observed + DBL_EPSILON) total_p += p;
            }
        }
        if (delta > 0) {
            current_a = mode_a - delta;
            if (current_a >= min_a) {
                int current_b = (a + b) - current_a;
                int current_c = (a + c) - current_a;
                int current_d = (b + d) - current_b;
                if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                    const double p =
                        std::exp(log_hypergeom_pg(current_a, current_b, current_c, current_d));
                    if (p <= p_observed + DBL_EPSILON) total_p += p;
                }
            }
        }
    }
    return total_p;
}

static int nt4_from_ref_char(char ch) {
    switch (std::toupper(static_cast<unsigned char>(ch))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
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
        if (off + static_cast<size_t>(len) > ref_seq.size() ||
            off2 + static_cast<size_t>(len) > ref_seq.size()) return false;
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

// ════════════════════════════════════════════════════════════════════════════
// Noisy reads ratio cache  (longcallD `build_var_noisy_reads_ratio_cache`)
// ════════════════════════════════════════════════════════════════════════════

/** Build cgranges interval trees for fast var_noisy_reads_ratio queries (lazy, called once). */
static void build_var_noisy_read_cache(BamChunk& chunk) {
    if (chunk.var_noisy_read_cov_cr) return;

    chunk.var_noisy_read_cov_cr.reset(cr_init());
    chunk.var_noisy_read_err_cr.reset(cr_init());
    chunk.var_noisy_read_marks.assign(chunk.reads.size(), 0);
    chunk.var_noisy_read_mark_id = 0;

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& r = chunk.reads[read_i];
        if (r.is_skipped || r.beg > r.end) continue;

        cr_add(chunk.var_noisy_read_cov_cr.get(), "cr",
               static_cast<int32_t>(r.beg - 1), static_cast<int32_t>(r.end),
               static_cast<int32_t>(read_i));

        // Merge adjacent XID events per read → noisy intervals.
        hts_pos_t noisy_start = -1;
        hts_pos_t noisy_end = -1;
        bool has_noisy = false;
        for (const DigarOp& d : r.digars) {
            if (d.type != DigarType::Snp && d.type != DigarType::Insertion &&
                d.type != DigarType::Deletion) continue;
            const hts_pos_t cur_start = d.pos - 1;
            hts_pos_t cur_end = d.pos;
            if (d.type == DigarType::Snp || d.type == DigarType::Deletion) cur_end += d.len - 1;
            if (!has_noisy) {
                noisy_start = cur_start;
                noisy_end = cur_end;
                has_noisy = true;
                continue;
            }
            if (cur_start < noisy_end) {
                if (cur_end > noisy_end) noisy_end = cur_end;
                continue;
            }
            cr_add(chunk.var_noisy_read_err_cr.get(), "cr",
                   static_cast<int32_t>(noisy_start), static_cast<int32_t>(noisy_end),
                   static_cast<int32_t>(read_i));
            noisy_start = cur_start;
            noisy_end = cur_end;
        }
        if (has_noisy) {
            cr_add(chunk.var_noisy_read_err_cr.get(), "cr",
                   static_cast<int32_t>(noisy_start), static_cast<int32_t>(noisy_end),
                   static_cast<int32_t>(read_i));
        }
    }
    cr_index(chunk.var_noisy_read_cov_cr.get());
    cr_index(chunk.var_noisy_read_err_cr.get());
}

/** longcallD `var_noisy_reads_ratio` (collect_var.c). */
static double var_noisy_reads_ratio(BamChunk& chunk, hts_pos_t var_start, hts_pos_t var_end) {
    build_var_noisy_read_cache(chunk);

    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    const int64_t total_n = cr_overlap(chunk.var_noisy_read_cov_cr.get(), "cr",
                                       static_cast<int32_t>(var_start - 1),
                                       static_cast<int32_t>(var_end), &ovlp_b, &max_b);
    const int total_reads = static_cast<int>(total_n);
    int noisy_reads = 0;

    if (total_reads > 0) {
        ++chunk.var_noisy_read_mark_id;
        if (chunk.var_noisy_read_mark_id == INT_MAX) {
            std::fill(chunk.var_noisy_read_marks.begin(), chunk.var_noisy_read_marks.end(), 0);
            chunk.var_noisy_read_mark_id = 1;
        }
        int64_t* noisy_b = nullptr;
        int64_t noisy_max_b = 0;
        const int64_t noisy_n = cr_overlap(chunk.var_noisy_read_err_cr.get(), "cr",
                                           static_cast<int32_t>(var_start - 1),
                                           static_cast<int32_t>(var_end), &noisy_b, &noisy_max_b);
        for (int64_t i = 0; i < noisy_n; ++i) {
            const int ri =
                static_cast<int>(cr_label(chunk.var_noisy_read_err_cr.get(), noisy_b[i]));
            if (ri < 0 || static_cast<size_t>(ri) >= chunk.var_noisy_read_marks.size()) continue;
            if (chunk.var_noisy_read_marks[static_cast<size_t>(ri)] == chunk.var_noisy_read_mark_id)
                continue;
            chunk.var_noisy_read_marks[static_cast<size_t>(ri)] = chunk.var_noisy_read_mark_id;
            ++noisy_reads;
        }
        std::free(noisy_b);
    }
    std::free(ovlp_b);
    return total_reads > 0 ? static_cast<double>(noisy_reads) / static_cast<double>(total_reads) : 0.0;
}

/** longcallD `cr_add_var_cr` (collect_var.c). */
static void cr_add_var_to_noisy_cr(cgranges_t* var_cr,
                                   cgranges_t* low_comp_cr,
                                   const VariantKey& var,
                                   BamChunk& chunk,
                                   bool check_noisy_reads_ratio,
                                   const Options& opts) {
    hts_pos_t var_start = 0;
    hts_pos_t var_end = 0;
    variant_genomic_span(var, var_start, var_end);
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        int64_t* low_comp_b = nullptr;
        int64_t max_low_comp_n = 0;
        const int64_t low_comp_n = cr_overlap(low_comp_cr, "cr",
                                              static_cast<int32_t>(var_start - 1),
                                              static_cast<int32_t>(var_end),
                                              &low_comp_b, &max_low_comp_n);
        for (int64_t j = 0; j < low_comp_n; ++j) {
            const int32_t start = cr_start(low_comp_cr, low_comp_b[j]) + 1;
            const int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
            if (start < var_start) var_start = start;
            if (end > var_end) var_end = end;
        }
        std::free(low_comp_b);
    }
    if (!check_noisy_reads_ratio || var_noisy_reads_ratio(chunk, var_start, var_end) >= opts.min_af) {
        cr_add(var_cr, "cr",
               static_cast<int32_t>(var_start - 1), static_cast<int32_t>(var_end), 1);
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Variant classification
// ════════════════════════════════════════════════════════════════════════════

/** longcallD `classify_var_cate` first stage (strand bias only if `is_ont`). */
static VariantCategory classify_variant_initial(const VariantKey& key,
                                                VariantCounts& counts,
                                                const std::string& ref_slice,
                                                hts_pos_t ref_beg,
                                                hts_pos_t ref_end,
                                                const Options& opts) {
    const int depth_with_low_quality = counts.total_cov + counts.low_qual_cov;
    counts.allele_fraction =
        counts.total_cov == 0
            ? 0.0
            : static_cast<double>(counts.alt_cov) / static_cast<double>(counts.total_cov);

    if (depth_with_low_quality < opts.min_depth) return VariantCategory::LowCoverage;
    if (counts.alt_cov < opts.min_alt_depth) return VariantCategory::LowCoverage;

    if (opts.is_ont()) {
        const int fa = counts.forward_alt;
        const int ra = counts.reverse_alt;
        const int expected = (fa + ra) / 2;
        if (expected > 0) {
            const double p = fisher_exact_two_tail(fa, ra, expected, expected);
            if (p < opts.strand_bias_pval) return VariantCategory::StrandBias;
        }
    }

    if (counts.allele_fraction < opts.min_af) return VariantCategory::LowAlleleFraction;
    if (counts.allele_fraction > opts.max_af) return VariantCategory::CleanHom;

    if ((key.type == VariantType::Insertion || key.type == VariantType::Deletion) &&
        (var_is_homopolymer_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps) ||
         var_is_repeat_region_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps))) {
        return VariantCategory::RepeatHetIndel;
    }
    return key.type == VariantType::Snp ? VariantCategory::CleanHetSnp : VariantCategory::CleanHetIndel;
}

/** longcallD `classify_cand_vars` core: noisy overlap, REP/dense loci → noisy_var_cr. */
static void classify_cand_vars_pgphase(BamChunk& chunk, const Options& opts) {
    CandidateTable& variants = chunk.candidates;
    if (variants.empty()) return;

    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(chunk.low_complexity_regions));
    cgranges_t* low_comp = low_own.cr;
    if (low_comp != nullptr) cr_index(low_comp);

    // Populate ref_base (SNPs) and alt_ref_base (INS/DEL) from reference sequence.
    for (auto& cv : variants) {
        if (cv.key.type == VariantType::Snp) {
            cv.ref_base = static_cast<uint8_t>(ref_nt4_at(chunk.ref_seq, chunk.ref_beg, cv.key.pos));
        } else {
            cv.alt_ref_base =
                static_cast<uint8_t>(ref_nt4_at(chunk.ref_seq, chunk.ref_beg, cv.key.pos - 1));
        }
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
        if (opts.is_ont() && c == VariantCategory::StrandBias) continue;
        const VariantKey& k = variants[i].key;
        if (k.type == VariantType::Insertion) {
            cr_add(var_pos_cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), 1);
        } else {
            cr_add(var_pos_cr, "cr",
                   static_cast<int32_t>(k.pos - 1),
                   static_cast<int32_t>(k.pos + k.ref_len - 1), 1);
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

        if (!opts.is_ont() && chunk_noisy->n_r > 0) {
            int64_t n = 0;
            if (variants[i].key.type == VariantType::Insertion) {
                n = cr_overlap(chunk_noisy, "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos),
                               &ovlp_b, &max_b);
            } else {
                n = cr_overlap(chunk_noisy, "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                               &ovlp_b, &max_b);
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
            n_ov = cr_overlap(var_pos_cr, "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos),
                              &ovlp_b, &max_b);
        } else {
            n_ov = cr_overlap(var_pos_cr, "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                              &ovlp_b, &max_b);
        }
        if (n_ov > 1) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, true, opts);
            }
        }

        if (c == VariantCategory::LowAlleleFraction) cats[i] = VariantCategory::LowCoverage;
    }

    std::free(ovlp_b);
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

/** Promote NoisyCandidate / RepeatHetIndel variants whose alt or ref span is >= kMinSvLen to
 *  NoisyResolved, so they are retained in the output TSV as large SV evidence. */
static void resolve_simple_noisy_candidates(CandidateTable& variants) {
    for (CandidateVariant& candidate : variants) {
        if (candidate.counts.category != VariantCategory::NoisyCandidate &&
            candidate.counts.category != VariantCategory::RepeatHetIndel) continue;
        if (candidate.key.alt.size() >= static_cast<size_t>(kMinSvLen) ||
            candidate.key.ref_len >= kMinSvLen) {
            candidate.counts.category = VariantCategory::NoisyResolved;
        }
    }
}

/** Top-level classification entry point for one chunk: classify all candidates, then resolve
 *  large noisy SVs. Called once per chunk after allele counting, before post_process_noisy_regs. */
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    (void)header;
    classify_cand_vars_pgphase(chunk, opts);
    resolve_simple_noisy_candidates(chunk.candidates);
}

} // namespace pgphase_collect
