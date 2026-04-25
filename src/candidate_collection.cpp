#include "candidate_collection.hpp"

#include <algorithm>
#include <cstring>
#include <iterator>
#include <utility>

namespace pgphase_collect {

VariantKey variant_key_from_digar(int tid, const DigarOp& op) {
    if (op.type == DigarType::Snp) {
        return VariantKey{tid, op.pos, VariantType::Snp, 1, op.alt};
    }
    if (op.type == DigarType::Insertion) {
        return VariantKey{tid, op.pos, VariantType::Insertion, 0, op.alt};
    }
    return VariantKey{tid, op.pos, VariantType::Deletion, op.len, ""};
}

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

bool VariantKeyLess::operator()(const VariantKey& lhs, const VariantKey& rhs) const {
    return exact_comp_var_site(&lhs, &rhs) < 0;
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

static bool is_variant_digar_for_cand_sweep(const DigarOp& d) {
    return d.type == DigarType::Snp || d.type == DigarType::Insertion || d.type == DigarType::Deletion;
}

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
    CandidateTable collapsed;
    collapsed.reserve(variants.size());
    for (auto& candidate : variants) {
        if (!collapsed.empty() &&
            exact_comp_var_site_ins(&collapsed.back().key, &candidate.key, kLongcalldMinSvLen) == 0) {
            continue;
        }
        collapsed.push_back(std::move(candidate));
    }
    variants.swap(collapsed);
}
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants) {
    variants.clear();
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            variants.push_back(CandidateVariant{variant_key_from_digar(read.tid, digar), VariantCounts{}});
        }
    }

    collapse_fuzzy_large_insertions(variants);
}

void add_coverage(VariantCounts& counts, bool reverse, bool alt, bool low_quality) {
    if (low_quality) {
        counts.low_qual_cov++;
        return;
    }

    counts.total_cov++;
    if (alt) {
        counts.alt_cov++;
        reverse ? counts.reverse_alt++ : counts.forward_alt++;
    } else {
        counts.ref_cov++;
        reverse ? counts.reverse_ref++ : counts.forward_ref++;
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
                add_coverage(c, read.reverse, false, false);
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
                    variants[static_cast<size_t>(site_i)].counts, read.reverse, true, lq);
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
                if (stid < read.tid) {
                    continue;
                }
                break;
            }
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) {
                break;
            }
            add_coverage(variants[static_cast<size_t>(site_i)].counts, read.reverse, false, false);
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

} // namespace pgphase_collect
