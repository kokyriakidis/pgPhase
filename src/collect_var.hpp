#ifndef PGPHASE_COLLECT_VAR_HPP
#define PGPHASE_COLLECT_VAR_HPP

#include "collect_types.hpp"

#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

// ── RAII wrapper for cgranges_t ──────────────────────────────────────────────
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

// ── Variant site comparison ──────────────────────────────────────────────────
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2);
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len);

struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const;
};

// ── Candidate site collection ────────────────────────────────────────────────
/**
 * @brief Isolates variant metrics out from a parsed `DigarOp` node into a searchable key.
 */
VariantKey variant_key_from_digar(int tid, const DigarOp& op);

/** Collapses overlapping/adjacent insertion elements bridging fragmented representation mappings inside `longcallD`. */
void collapse_fuzzy_large_insertions(CandidateTable& variants);

/** Core parser evaluating `ReadRecord` groups pulling mapped components into the statistical candidate set. */
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants);

// ── Allele counting ──────────────────────────────────────────────────────────
/**
 * @brief Resolves depth coverage arrays spanning overlapping candidate SNPs/Indels.
 * 
 * Functions identically to `update_cand_var_depth`, feeding specific read inclusion constraints
 * based on minimum qualities and explicit base mismatches. Populates statistical weights.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq);

// ── Interval utilities ───────────────────────────────────────────────────────

/** Merges overlapping structural intervals inside a vector. (Maps to `cr_merge`/cgranges). */
void merge_intervals(std::vector<Interval>& intervals);

/** Converts a vector of coordinate spans into an optimized C-based interval tree interval tree (`cgranges_t`). */
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals);

/** Deserializes coordinate spans out from `cgranges_t` back into vectors. */
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out);

/** Builds an interval tree encapsulating all noisy clipped stretches inside a specific read. */
cgranges_t* build_read_noisy_cr(const ReadRecord& read);

/** Gets 1-based start and end boundary positions out of a Variant token. */
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end);

// ── Noisy region processing ──────────────────────────────────────────────────

/** 
 * @brief Expands structurally noisy alignments via standard heuristics mapping `longcallD`'s `pre_process_noisy_regs`. 
 */
void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts);

/** 
 * @brief Identifies whether a previously clean candidate overlaps with the final bounds of a noisy block. 
 * Connects to `apply_noisy_containment_filter()`.
 */
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand);

/** 
 * @brief Core filtering phase applying `cgranges` overlap scans against `var_t` arrays masking false SNVs heavily contained inside noisy spans. 
 */
void apply_noisy_containment_filter(BamChunk& chunk);

/** Maps intervals containing long repeating simple subsequences using `sdust` library loops. */
void populate_low_complexity_intervals(BamChunk& chunk);

/** Connects index values map to efficiently query raw alignments during grouping. */
void populate_chunk_read_indexes(BamChunk& chunk);

// ── Variant classification ───────────────────────────────────────────────────

/** Aggregates variant categories passing through filtering (AF, Depth, Content) mapping tags inside `CandidateTable`. */
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
