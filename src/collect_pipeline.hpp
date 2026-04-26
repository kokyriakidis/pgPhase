#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

/**
 * @brief Public API for region chunking and parallel collect-bam-variation orchestration.
 * 
 * Maps to the upper-level flow control of `longcallD`, breaking down genomic queries
 * into thread-safe discrete blocks (`RegionChunk`), distributing them to parser workers, 
 * and funneling the extracted variants into unified result models.
 */

// ── Region chunking ──────────────────────────────────────────────────────────

/** Parses a chromosome region string (e.g., "chr1:100-200") into a filter struct. */
RegionFilter parse_region(const std::string& region);

/** Loads standard 3-column BED regions into inclusion filters. */
std::vector<RegionFilter> load_bed_regions(const std::string& path);

/** Divides genomic contigs into evenly sized chunks for parallel read digestion mapping to `longcallD` thread dispatchers. */
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai);

/** Identifies genomic ranges to scan, pulling from opts.region_file or full indices. */
std::vector<RegionChunk> load_region_chunks(const Options& opts);

// ── Pipeline ─────────────────────────────────────────────────────────────────

/**
 * @brief Core parallelization dispatcher replacing `longcallD`'s thread pools.
 * 
 * Dispatches chunks to `process_chunk`, waits for all variants, and correctly merges 
 * the results back into a globally-sorted, duplicate-free `CandidateTable`.
 */
CandidateTable collect_chunks_parallel(
    const Options& opts,
    const std::vector<RegionChunk>& chunks,
    std::vector<std::vector<ReadSupportRow>>* read_support_batches);

/** Main driver method that boots the full pipeline step-by-step. */
void run_collect_bam_variation(const Options& opts);

} // namespace pgphase_collect

// ── CLI entry point ──────────────────────────────────────────────────────────
int collect_bam_variation(int argc, char* argv[]);

#endif
