#ifndef PGPHASE_COLLECT_VAR_HPP
#define PGPHASE_COLLECT_VAR_HPP

// Candidate collection, intervals, noisy regions, and classification API.
//
// The public entry point for the per-chunk biological workflow is
// collect_var_main.  Higher-level code owns BAM/FASTA I/O, chunk loading,
// threading, and output streaming.

#include "collect_types.hpp"

#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

// Move-only RAII owner for heap-allocated cgranges_t.
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

// Total order on VariantKey: negative if *var1 < *var2, zero if equal, positive if greater.
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2);

// Total order on CandidateVariant via canonical VariantKey encoding.
int exact_comp_cand_var(const CandidateVariant* var1, const CandidateVariant* var2);

// Same ordering with large-insertion fuzzy collapse: insertions within
// min_sv_len of each other at the same position are treated as equivalent.
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len);

// Strict-weak ordering functor delegating to exact_comp_var_site.
struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const;
};

// Maps one DigarOp (SNP/indel alignment event) to a VariantKey.
VariantKey variant_key_from_digar(int tid, const DigarOp& op);

// Sort and deduplicate candidates using fuzzy large-insertion equivalence.
// SNP/INS ALT tie-break uses nt4-encoded bytes, not raw ASCII order.
void collapse_fuzzy_large_insertions(CandidateTable& variants, int min_sv_len);

// Sort and drop rows whose VariantKey compares equal under exact_comp_var_site.
// Used when concatenating per-window candidate tables — cross-chunk dedupe
// must not fuzzy-merge large insertions (only per-chunk collapse does that).
void collapse_exact_duplicate_variants(CandidateTable& variants);

// Gather candidate keys from digars overlapping the chunk region, then
// collapse fuzzy large insertions.  Output counts are still zero.
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants,
                                          int min_sv_len);

// Fill allele/strand depth from digars vs sorted candidates.
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        int min_bq,
                                        int min_sv_len);

// Build an unindexed cgranges_t on synthetic contig "cr" (0-based half-open).
// Returns nullptr if intervals is empty.
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals);

// Read intervals from contig "cr" into out (1-based inclusive).
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out);

// Build an indexed cgranges_t for one read's noisy_regions.
cgranges_t* build_read_noisy_cr(const ReadRecord& read);

// Reference span for overlap and noisy logic (insertion is zero-width).
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end);

// Merge and filter read-level noisy intervals into chunk.noisy_regions.
// Called after candidate sites and allele counts, before classification.
void pre_process_noisy_regs_pgphase(PhasingChunk& chunk, const Options& opts);

// Expand noisy intervals using classified candidates and re-merge.
// Called after classification to capture additional noisy spans.
void post_process_noisy_regs_pgphase(PhasingChunk& chunk, const CandidateTable& cand);

// Set category to NonVariant for sites fully contained in noisy spans (non-ONT).
void apply_noisy_containment_filter(PhasingChunk& chunk);

// Run sdust on chunk.ref_seq to fill low_complexity_regions.
void populate_low_complexity_intervals(PhasingChunk& chunk);

// Fill ordered_read_ids and union per-read noisy intervals into the chunk list.
void populate_chunk_read_indexes(PhasingChunk& chunk);

// Two-pass candidate classification: first pass assigns initial categories
// (clean het/hom, low-cov, strand-bias, repeat-het), second pass applies
// noisy-region adjustments and AF→LowCov rewrites.
void classify_chunk_candidates(PhasingChunk& chunk, const Options& opts, const bam_hdr_t* header);

// Steps 1-2: candidate discovery + classification.  Collects sites from
// digars, counts alleles, classifies candidates, and prunes NON_VAR.
// After this call, chunk.candidates is sorted and classified but no
// read profiles or k-means phasing has been done.
void collect_var_classify(PhasingChunk& chunk,
                          const Options& opts,
                          const bam_hdr_t* header);

// Step 3.1: build per-read variant profiles from digars.
// Populates chunk.read_var_profile and chunk.read_var_cr.
void collect_var_build_profiles(PhasingChunk& chunk, const Options& opts);

// Steps 3.2-4: k-means phasing and noisy-region MSA recall.
// Expects read profiles to be populated (by collect_var_build_profiles
// and/or manual injection for hybrid reads).
void collect_var_run_phasing(PhasingChunk& chunk, const Options& opts);

// Steps 3-4 combined: build profiles + run phasing.
void collect_var_phase(PhasingChunk& chunk,
                       const Options& opts);

// Full pipeline: collect_var_classify + collect_var_phase.
// Expects chunk.reads, chunk.ref_seq, read indexes, low-complexity intervals,
// and initial noisy regions to already be populated by the BAM/FASTA layer.
void collect_var_main(PhasingChunk& chunk,
                      const Options& opts,
                      const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
