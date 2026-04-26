#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

/**
 * Public API for region chunking and collect-bam-variation orchestration.
 * Implementation and CLI parsing: collect_pipeline.cpp.
 */

// ── Region chunking ──────────────────────────────────────────────────────────
RegionFilter parse_region(const std::string& region);
std::vector<RegionFilter> load_bed_regions(const std::string& path);
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai);
std::vector<RegionChunk> load_region_chunks(const Options& opts);

// ── Pipeline ─────────────────────────────────────────────────────────────────
CandidateTable collect_chunks_parallel(
    const Options& opts,
    const std::vector<RegionChunk>& chunks,
    std::vector<std::vector<ReadSupportRow>>* read_support_batches);

void run_collect_bam_variation(const Options& opts);

} // namespace pgphase_collect

// ── CLI entry point ──────────────────────────────────────────────────────────
int collect_bam_variation(int argc, char* argv[]);

#endif
