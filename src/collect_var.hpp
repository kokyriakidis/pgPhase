#ifndef PGPHASE_COLLECT_VAR_HPP
#define PGPHASE_COLLECT_VAR_HPP

#include "collect_types.hpp"

#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

/**
 * @brief RAII wrapper around `cgranges_t` (move-only, destroys on reset).
 */
struct CrangesOwner {
    cgranges_t* cr = nullptr;
    CrangesOwner() = default;
    CrangesOwner(const CrangesOwner&) = delete;
    CrangesOwner& operator=(const CrangesOwner&) = delete;
    CrangesOwner(CrangesOwner&& other) noexcept : cr(other.cr) { other.cr = nullptr; }
    CrangesOwner& operator=(CrangesOwner&& other) noexcept {
        if (this != &other) { reset(); cr = other.cr; other.cr = nullptr; }
        return *this;
    }
    ~CrangesOwner() { reset(); }
    void reset() { if (cr) { cr_destroy(cr); cr = nullptr; } }
    void adopt(cgranges_t* p) { reset(); if (p) cr = p; }
    cgranges_t* release() { cgranges_t* t = cr; cr = nullptr; return t; }
};

/** @brief Total order on `VariantKey` (longcallD `exact_comp_var_site`). */
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2);
/** @brief Same order with large-insertion fuzzy merge (longcallD `exact_comp_var_site_ins`). */
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len);

/**
 * @brief Strict-weak ordering functor delegating to `exact_comp_var_site`.
 */
struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const;
};

/**
 * @brief Maps one `DigarOp` to a `VariantKey` for the candidate table.
 */
VariantKey variant_key_from_digar(int tid, const DigarOp& op);

/**
 * @brief Deduplicates adjacent-equivalent candidates after sorting (large INS fuzzy rule).
 */
void collapse_fuzzy_large_insertions(CandidateTable& variants);

/**
 * @brief Collects candidate sites from reads overlapping a chunk, then fuzzy-deduplicates.
 */
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants);

/**
 * @brief Allele and strand depth counts plus optional read-support rows.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq);

/**
 * @brief Merges overlapping/adjacent intervals in place (max label wins).
 */
void merge_intervals(std::vector<Interval>& intervals);

/**
 * @brief Converts intervals to an unindexed `cgranges_t` on synthetic contig `"cr"`.
 */
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals);

/**
 * @brief Copies intervals from `"cr"` back into a vector.
 */
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out);

/**
 * @brief Indexed interval tree of one read's noisy subregions.
 */
cgranges_t* build_read_noisy_cr(const ReadRecord& read);

/**
 * @brief Reference span for a variant key (1-based; insertion is zero-width).
 */
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end);

/**
 * @brief Merge and filter read-level noisy intervals into `chunk.noisy_regions`.
 */
void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts);

/**
 * @brief Expand noisy bounds using classified candidates; re-merge.
 */
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand);

/**
 * @brief Marks candidates fully inside noisy spans as `NonVariant` (non-ONT path).
 */
void apply_noisy_containment_filter(BamChunk& chunk);

/**
 * @brief Runs sdust on chunk reference to fill `low_complexity_regions`.
 */
void populate_low_complexity_intervals(BamChunk& chunk);

/**
 * @brief Builds `ordered_read_ids` and unions read noisy intervals into the chunk.
 */
void populate_chunk_read_indexes(BamChunk& chunk);

/**
 * @brief Classifies all candidates in `chunk` (strand bias, depth, noisy overlap, containment).
 */
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
