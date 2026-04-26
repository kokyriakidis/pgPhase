#ifndef PGPHASE_COLLECT_VAR_HPP
#define PGPHASE_COLLECT_VAR_HPP

/**
 * @file collect_var.hpp
 * @brief Candidate collection, intervals, noisy regions, and classification API.
 */

#include "collect_types.hpp"

#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

/**
 * @brief Move-only RAII owner for heap-allocated `cgranges_t`.
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
    /** @brief Destroys held tree and nulls pointer. */
    void reset() { if (cr) { cr_destroy(cr); cr = nullptr; } }
    /** @brief Replaces ownership with \a p (previous tree destroyed). */
    void adopt(cgranges_t* p) { reset(); if (p) cr = p; }
    /** @brief Returns raw pointer and clears ownership without destroy. */
    cgranges_t* release() { cgranges_t* t = cr; cr = nullptr; return t; }
};

/**
 * @brief Total order on `VariantKey` (longcallD `exact_comp_var_site`).
 * @return Negative if `*var1 < *var2`, zero if equal, positive if greater.
 */
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2);

/**
 * @brief Same ordering with large-insertion fuzzy collapse (longcallD `exact_comp_var_site_ins`).
 * @param min_sv_len Length threshold for insertion fuzzy merge rule.
 * @return Negative / zero / positive as `exact_comp_var_site`.
 */
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len);

/**
 * @brief Strict-weak ordering functor delegating to `exact_comp_var_site`.
 */
struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const;
};

/**
 * @brief Maps one `DigarOp` to a `VariantKey` (longcallD `make_var_site_from_digar`).
 * @param tid Read contig index.
 * @param op SNP/indel digar operation.
 */
VariantKey variant_key_from_digar(int tid, const DigarOp& op);

/**
 * @brief Sorts and deduplicates candidates using fuzzy large-insertion equivalence.
 * @param variants In/out table.
 */
void collapse_fuzzy_large_insertions(CandidateTable& variants);

/**
 * @brief Gathers candidate keys from digars overlapping `chunk`, then collapses fuzzy INS.
 * @param chunk 1-based inclusive region bounds.
 * @param reads Parsed reads for the chunk.
 * @param variants Output candidate rows (counts still zero).
 */
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants);

/**
 * @brief Fills allele/strand depth from digars vs sorted candidates; optional read-support log.
 * @param reads Chunk reads.
 * @param variants Sorted candidate table.
 * @param chunk_region If non-null with `read_support_out`, stamped into support rows.
 * @param read_support_out Optional observation list.
 * @param min_bq Minimum base quality for high-quality alt tally.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq);

/**
 * @brief Sorts and merges overlapping/adjacent intervals; merged label is max of inputs.
 * @param intervals In/out vector.
 */
void merge_intervals(std::vector<Interval>& intervals);

/**
 * @brief Builds unindexed `cgranges_t` on synthetic contig `"cr"` (0-based half-open storage).
 * @return New tree or nullptr if nothing to add.
 */
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals);

/**
 * @brief Reads intervals from contig `"cr"` into \a out.
 * @param cr Source tree (may be null).
 * @param out Cleared and filled with 1-based inclusive intervals.
 */
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out);

/**
 * @brief Indexed `cgranges_t` for one read's `noisy_regions`.
 */
cgranges_t* build_read_noisy_cr(const ReadRecord& read);

/**
 * @brief Reference span used for overlap and noisy logic (insertion is zero-width).
 * @param key Variant locus.
 * @param var_start Set to first reference base of span (1-based).
 * @param var_end Set to last reference base (may be `key.pos - 1` for INS).
 */
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end);

/**
 * @brief Merges and filters read-level noisy intervals into `chunk.noisy_regions`.
 * @param chunk Chunk state with reads and optional low-complexity intervals.
 * @param opts Depth/ratio/merge thresholds.
 */
void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts);

/**
 * @brief Expands noisy intervals using classified candidates; re-merges (longcallD post-process).
 * @param chunk Chunk whose `noisy_regions` are updated.
 * @param cand Current candidate table (categories consulted).
 */
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand);

/**
 * @brief Sets category to `NonVariant` for sites fully contained in noisy spans (non-ONT).
 * @param chunk Chunk with candidates and `noisy_regions`.
 */
void apply_noisy_containment_filter(BamChunk& chunk);

/**
 * @brief Runs sdust on `chunk.ref_seq` slice to fill `low_complexity_regions`.
 * @param chunk Must have reference sequence loaded for the window.
 */
void populate_low_complexity_intervals(BamChunk& chunk);

/**
 * @brief Fills `ordered_read_ids` and unions per-read noisy intervals into the chunk list.
 * @param chunk Chunk whose `reads` are already populated.
 */
void populate_chunk_read_indexes(BamChunk& chunk);

/**
 * @brief Runs full candidate classification for one chunk (two-pass longcallD-shaped logic).
 * @param chunk Chunk with counts populated.
 * @param opts Technology and thresholds.
 * @param header BAM header (reserved for future use; currently unused).
 */
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
