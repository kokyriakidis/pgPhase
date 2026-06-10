#ifndef PGPHASE_HYBRID_INJECT_HPP
#define PGPHASE_HYBRID_INJECT_HPP

/// @file hybrid_inject.hpp
/// @brief Augment a BAM-derived PhasingChunk with graph snarl site
///        observations.  Used by the hybrid BAM+graph phasing mode.

#include "graph_query.hpp"
#include "graph_sites.hpp"
#include "phasing_types.hpp"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pgphase_collect {

/// Maps graph site keys to candidate indices in the augmented table.
using SiteToCandidateMap = std::unordered_map<std::string, int>;

/// Result of augmenting a BAM chunk with graph sites and reads.
struct HybridInjectionResult {
    int sites_added = 0;        // graph-only sites added to candidate table
    int sites_bridged = 0;      // existing BAM sites matched to graph sites
    int reads_injected = 0;     // graph-only reads added
};

/// Phase A: Augment the candidate table with graph-only sites.
///
/// For each snarl site in the chunk region, find matching BAM candidates
/// (exact position + REF/ALT match).  Sites with no BAM match are added
/// as new candidates.  Returns a map from graph site key to candidate index.
///
/// Call BEFORE collect_var_build_profiles.
SiteToCandidateMap inject_graph_sites(
    PhasingChunk& chunk,
    const GraphSiteCatalogView& graph_sites,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts,
    int* sites_bridged_out,
    int* sites_added_out,
    std::unordered_set<int>* graph_only_candidates_out = nullptr);

/// Phase B: Inject graph-only reads and extend doubly-mapped read profiles.
///
/// For each graph read NOT already in chunk.reads, build a synthetic
/// ReadRecord and ReadVariantProfile from its allele observations at
/// mapped candidate sites.  Appends to chunk.reads, chunk.read_var_profile,
/// and rebuilds chunk.read_var_cr.
///
/// For reads present in both BAM and GAF, extend their existing BAM
/// profiles with graph observations at graph-only candidate sites
/// (positions where the BAM profile has no informative call).
///
/// Call AFTER collect_var_build_profiles (so BAM profiles are already built).
int inject_graph_reads(
    PhasingChunk& chunk,
    const std::vector<GraphReadAllele>& graph_rows,
    const SiteToCandidateMap& site_to_candidate,
    const std::unordered_set<int>& graph_only_candidates,
    const Options& opts,
    int* reads_extended_out = nullptr);

/// Apply noise filter to graph-only candidates.
///
/// Reclassifies CleanHetIndel candidates added by inject_graph_sites
/// that sit in homopolymer/repeat/low-complexity reference contexts.
///
/// @param chunk                  Chunk with augmented candidates.
/// @param ref_seq                Reference sequence slice for the chunk region.
/// @param ref_beg                1-based start of ref_seq.
/// @param ref_end                1-based end of ref_seq (inclusive).
/// @param graph_only_candidates  Set of candidate indices added by inject_graph_sites.
/// @param max_xgaps              Maximum indel span to check.
void apply_hybrid_noise_filter(
    PhasingChunk& chunk,
    const std::string& ref_seq,
    hts_pos_t ref_beg,
    hts_pos_t ref_end,
    const std::unordered_set<int>& graph_only_candidates,
    int max_xgaps);

/// Backfill allele counts on graph-only candidates from BAM read profiles.
///
/// After collect_var_build_profiles, BAM reads have allele observations at
/// graph-only candidate positions (typically ref=0) but the candidate's
/// ref_cov/alt_cov/total_cov were never updated.  This function scans
/// existing BAM profiles and accumulates the missing counts.
///
/// Call AFTER collect_var_build_profiles, BEFORE inject_graph_reads.
void backfill_graph_candidate_counts(
    PhasingChunk& chunk,
    const std::unordered_set<int>& graph_only_candidates);

}  // namespace pgphase_collect

#endif  // PGPHASE_HYBRID_INJECT_HPP
