#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

/**
 * @file collect_pipeline.hpp
 * @brief Public declarations for region chunking and collect-bam-variation orchestration.
 */

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

/**
 * @brief Public API for region chunking and parallel collect-bam-variation orchestration.
 *
 * @details Mirrors longcallD-style region partitioning into `RegionChunk` units, worker dispatch,
 * and merged candidate tables.
 */

/**
 * @brief Parses `chr`, `chr:pos`, or `chr:start-end` (commas allowed) into `RegionFilter`.
 * @param region Region literal from CLI.
 * @return Filter; empty string yields `enabled == false`.
 */
RegionFilter parse_region(const std::string& region);

/**
 * @brief Loads 3-column BED regions as inclusion filters (0-based BED → 1-based inclusive).
 * @param path BED file path.
 * @return List of filters; throws on I/O or parse errors.
 */
std::vector<RegionFilter> load_bed_regions(const std::string& path);

/**
 * @brief Tiles the genome (or user filters) into chunks of `opts.chunk_size` and annotates neighbors.
 * @param opts Chunk size and region inputs.
 * @param header BAM header.
 * @param fai Reference index (contig presence checked).
 * @return Annotated chunk list.
 */
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai);

/**
 * @brief Opens primary BAM + reference index and returns chunks from `build_region_chunks`.
 * @param opts Must set `primary_bam_file()`, `ref_fasta`, and region fields.
 * @return Chunk list; throws if BAM/index/FAI cannot be opened.
 */
std::vector<RegionChunk> load_region_chunks(const Options& opts);

/**
 * @brief Runs all chunks through a worker pool and merges results into one `CandidateTable`.
 *
 * @param opts Thread count and paths.
 * @param chunks Region list from `load_region_chunks`.
 * @param read_support_batches If non-null, receives per-chunk read-support batches when the caller collects them.
 * @return Merged candidates in chunk processing order.
 */
CandidateTable collect_chunks_parallel(
    const Options& opts,
    const std::vector<RegionChunk>& chunks,
    std::vector<std::vector<ReadSupportRow>>* read_support_batches);

/**
 * @brief Streaming driver: batch by `reg_chunk_i`, write TSV/VCF/read-support incrementally.
 * @param opts Output paths, reference, and BAM list.
 */
void run_collect_bam_variation(const Options& opts);

} // namespace pgphase_collect

/**
 * @brief CLI entry for `collect-bam-variation` (argv without subcommand name).
 * @param argc Argument count.
 * @param argv Arguments after subcommand removal.
 * @return Exit code (0 success, 1 error).
 */
int collect_bam_variation(int argc, char* argv[]);

#endif
