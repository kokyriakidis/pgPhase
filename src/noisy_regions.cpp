#include "noisy_regions.hpp"

#include "sdust.h"

#include <algorithm>
#include <cstdlib>

namespace pgphase_collect {

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
    // Leave unindexed (longcallD-style): callers run `cr_index` once before overlap/merge.
    return cr;
}

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

void low_comp_cr_start_end(cgranges_t* low_comp_cr, int32_t start, int32_t end, int32_t* new_start, int32_t* new_end) {
    *new_start = start;
    *new_end = end;
    if (low_comp_cr == nullptr || low_comp_cr->n_r == 0) return;
    int64_t* low_comp_b = nullptr;
    int64_t low_comp_n = 0;
    int64_t max_low_comp_n = 0;
    low_comp_n = cr_overlap(low_comp_cr, "cr", start - 1, end, &low_comp_b, &max_low_comp_n);
    for (int64_t j = 0; j < low_comp_n; ++j) {
        const int32_t s = cr_start(low_comp_cr, low_comp_b[j]) + 1;
        const int32_t e = cr_end(low_comp_cr, low_comp_b[j]);
        if (s < *new_start) *new_start = s;
        if (e > *new_end) *new_end = e;
    }
    std::free(low_comp_b);
}

/**
 * longcallD `cr_extend_noisy_regs_with_low_comp` (collect_var.c): optional extension into
 * low-complexity intervals, then **always** `cr_merge(…, -1, merge_dis, min_sv_len)`.
 */
cgranges_t* cr_extend_noisy_regs_with_low_comp(cgranges_t* chunk_noisy_regs,
                                               cgranges_t* low_comp_cr,
                                               int merge_dis,
                                               int min_sv_len) {
    if (chunk_noisy_regs == nullptr) {
        return nullptr;
    }
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

cgranges_t* build_read_noisy_cr(const ReadRecord& read) {
    if (read.noisy_regions.empty()) return nullptr;
    cgranges_t* cr = intervals_to_cr(read.noisy_regions);
    if (cr != nullptr) cr_index(cr);
    return cr;
}

bool category_skipped_for_noisy_flank(VariantCategory c) {
    return c == VariantCategory::LowCoverage || c == VariantCategory::StrandBias ||
           c == VariantCategory::NonVariant;
}

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
    // longcallD uses `var->pos .. var->pos + ref_len - 1` here; insertions have ref_len=0.
    var_start = key.pos;
    var_end = key.pos - 1;
}
/** After `post_process_noisy_regs_pgphase`, longcallD `cr_is_contained` sweep (collect_var.c:1007–1017). */
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
            n = cr_is_contained(own.cr,
                                "cr",
                                static_cast<int32_t>(k.pos - 1),
                                static_cast<int32_t>(k.pos + k.ref_len),
                                &b,
                                &m);
        }
        if (n > 0) {
            cv.counts.category = VariantCategory::NonVariant;
        }
    }
    std::free(b);
}

void collect_noisy_reg_start_end_pgphase(const cgranges_t* chunk_noisy_regs,
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
            max_left_var_i[static_cast<size_t>(reg_i)] = std::min(std::max(0, static_cast<int>(cand.size()) - 1), 0);
        }
        if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
            min_right_var_i[static_cast<size_t>(reg_i)] = std::max(0, static_cast<int>(cand.size()) - 1);
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
    if (low_own.cr != nullptr) {
        cr_index(low_own.cr);
    }

    const int merge_dis = opts.noisy_reg_merge_dis;
    const int msv = opts.min_sv_len;

    cgranges_t* noisy = noisy_own.release();
    // longcallD: `cr_extend_noisy_regs_with_low_comp` then a second `cr_merge` (collect_var.c:567–568).
    noisy = cr_extend_noisy_regs_with_low_comp(noisy, low_own.cr, merge_dis, msv);
    noisy = cr_merge(noisy, -1, merge_dis, msv);
    noisy_own.adopt(noisy);

    cgranges_t* noisy_regs = noisy_own.cr;
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    std::vector<uint8_t> skip_noisy_reg(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_total(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_noisy(static_cast<size_t>(noisy_regs->n_r), 0);

    // longcallD: `ordered_read_ids` and skip `is_skipped[read_i]`.
    const std::vector<int>& read_order = chunk.ordered_read_ids;
    const auto visit_read = [&](const ReadRecord& read) {
        if (read.is_skipped) return;
        CrangesOwner read_noisy_own;
        read_noisy_own.adopt(build_read_noisy_cr(read));
        const int64_t beg = read.beg;
        const int64_t end = read.end;
        const int64_t ovlp_n = cr_overlap(
            noisy_regs, "cr", static_cast<int32_t>(beg - 1), static_cast<int32_t>(end), &ovlp_b, &max_b);
        for (int64_t ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            const int r_i = static_cast<int>(ovlp_b[ovlp_i]);
            noisy_reg_to_total[static_cast<size_t>(r_i)]++;
            const int32_t noisy_reg_start = cr_start(noisy_regs, r_i) + 1;
            const int32_t noisy_reg_end = cr_end(noisy_regs, r_i);
            int64_t* noisy_digar_b = nullptr;
            int64_t noisy_digar_max_b = 0;
            int64_t noisy_digar_ovlp_n = 0;
            if (read_noisy_own.cr != nullptr) {
                noisy_digar_ovlp_n = cr_overlap(
                    read_noisy_own.cr,
                    "cr",
                    noisy_reg_start - 1,
                    noisy_reg_end,
                    &noisy_digar_b,
                    &noisy_digar_max_b);
            }
            if (noisy_digar_ovlp_n > 0) {
                noisy_reg_to_noisy[static_cast<size_t>(r_i)]++;
            }
            std::free(noisy_digar_b);
        }
    };
    if (read_order.empty() && !chunk.reads.empty()) {
        for (const ReadRecord& r : chunk.reads) {
            visit_read(r);
        }
    } else {
        for (int ord : read_order) {
            if (ord < 0 || static_cast<size_t>(ord) >= chunk.reads.size()) continue;
            visit_read(chunk.reads[static_cast<size_t>(ord)]);
        }
    }
    std::free(ovlp_b);
    ovlp_b = nullptr;

    const int min_noisy_reg_reads = opts.min_alt_depth;
    const float min_noisy_reg_ratio = static_cast<float>(opts.min_af);
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        const int n_noisy = noisy_reg_to_noisy[static_cast<size_t>(i)];
        const int n_total = noisy_reg_to_total[static_cast<size_t>(i)];
        if (n_noisy < min_noisy_reg_reads ||
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
            cr_add(new_noisy_regs,
                   "cr",
                   cr_start(noisy_regs, i),
                   cr_end(noisy_regs, i),
                   cr_label(noisy_regs, i));
        }
        cr_index(new_noisy_regs);
        // Let `noisy_own.adopt` destroy the previous cgranges via `reset()`; do not
        // `cr_destroy(noisy_regs)` here or the owner still holds a dangling pointer.
        noisy_own.adopt(new_noisy_regs);
        noisy_regs = noisy_own.cr;
    }

    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

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
    collect_noisy_reg_start_end_pgphase(chunk_noisy_regs, cand, cand_cate, noisy_reg_start, noisy_reg_end);

    cgranges_t* noisy_regs = cr_init();
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        cr_add(noisy_regs,
               "cr",
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
void populate_low_complexity_intervals(BamChunk& chunk) {
    chunk.low_complexity_regions.clear();
    if (chunk.ref_seq.empty()) return;

    const hts_pos_t beg = std::max(chunk.region.beg, chunk.ref_beg);
    const hts_pos_t end = std::min(chunk.region.end, chunk.ref_end);
    if (beg > end) return;

    const size_t offset = static_cast<size_t>(beg - chunk.ref_beg);
    const int len = static_cast<int>(end - beg + 1);
    int n = 0;
    uint64_t* intervals = sdust(
        nullptr,
        reinterpret_cast<const uint8_t*>(chunk.ref_seq.data() + offset),
        len,
        kSdustThreshold,
        kSdustWindow,
        &n);
    for (int i = 0; i < n; ++i) {
        const hts_pos_t rel_beg = static_cast<hts_pos_t>(intervals[i] >> 32);
        const hts_pos_t rel_end = static_cast<hts_pos_t>(static_cast<uint32_t>(intervals[i]));
        if (rel_end <= rel_beg) continue;
        chunk.low_complexity_regions.push_back(Interval{
            beg + rel_beg,
            beg + rel_end - 1,
            static_cast<int>(rel_end - rel_beg)});
    }
    std::free(intervals);
}

void populate_chunk_read_indexes(BamChunk& chunk) {
    chunk.ordered_read_ids.clear();
    chunk.noisy_regions.clear();

    chunk.ordered_read_ids.reserve(chunk.reads.size());
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        if (read.is_skipped) continue;
        chunk.ordered_read_ids.push_back(static_cast<int>(read_i));
        chunk.noisy_regions.insert(chunk.noisy_regions.end(), read.noisy_regions.begin(), read.noisy_regions.end());
    }
    merge_intervals(chunk.noisy_regions);
}

} // namespace pgphase_collect
