#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

/**
 * @brief Public API for region chunking and parallel collect-bam-variation orchestration.
 *
 * Mirrors longcallD-style region partitioning into `RegionChunk` units, worker dispatch,
 * and merged candidate tables.
 */

/**
 * @brief Parses `chr`, `chr:pos`, or `chr:start-end` (commas allowed) into `RegionFilter`.
 */
RegionFilter parse_region(const std::string& region);

/**
 * @brief Loads 3-column BED regions as inclusion filters (0-based BED → 1-based inclusive).
 */
std::vector<RegionFilter> load_bed_regions(const std::string& path);

/**
 * @brief Tiles the genome (or user filters) into chunks of `opts.chunk_size` and annotates neighbors.
 */
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai);

/**
 * @brief Opens primary BAM + reference index and returns chunks from `build_region_chunks`.
 */
std::vector<RegionChunk> load_region_chunks(const Options& opts);

/**
 * @brief Runs all chunks through a worker pool and merges results into one `CandidateTable`.
 *
 * @param opts Thread count and paths.
 * @param chunks Region list from `load_region_chunks`.
 * @param read_support_batches If non-null, receives one batch of rows per chunk when requested by caller.
 */
CandidateTable collect_chunks_parallel(
    const Options& opts,
    const std::vector<RegionChunk>& chunks,
    std::vector<std::vector<ReadSupportRow>>* read_support_batches);

/**
 * @brief Streaming driver: batch by `reg_chunk_i`, write TSV/VCF/read-support incrementally.
 */
void run_collect_bam_variation(const Options& opts);

} // namespace pgphase_collect

/**
 * @brief CLI entry for `collect-bam-variation` (argv without subcommand name).
 */
int collect_bam_variation(int argc, char* argv[]);

#endif
