#ifndef PGPHASE_COLLECT_PHASE_PGBAM_HPP
#define PGPHASE_COLLECT_PHASE_PGBAM_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

/**
 * @brief Merge phase blocks within one chunk using annotated-BAM `hs` tags and `.pgbam` threads.
 */
void stitch_phase_blocks_with_pgbam(BamChunk& chunk, const PgbamSidecarData& sidecar);

/**
 * @brief Merge phase blocks within one chunk with explicit `.pgbam` evidence thresholds.
 */
void stitch_phase_blocks_with_pgbam(BamChunk& chunk,
                                    const PgbamSidecarData& sidecar,
                                    int min_winning_threads,
                                    int polarity_margin);

/**
 * @brief Cleanup pass over final contig phase blocks using annotated-BAM `hs` tags and `.pgbam` threads.
 */
void stitch_phase_blocks_with_pgbam(std::vector<BamChunk>& chunks,
                                    const PgbamSidecarData& sidecar,
                                    int min_winning_threads);

/**
 * @brief Cleanup pass over final contig phase blocks with an explicit thread-polarity margin.
 */
void stitch_phase_blocks_with_pgbam(std::vector<BamChunk>& chunks,
                                    const PgbamSidecarData& sidecar,
                                    int min_winning_threads,
                                    int polarity_margin);

/**
 * @brief Fallback stitch for adjacent chunks when common-read overlap has no decisive signal.
 */
bool stitch_adjacent_chunks_with_pgbam(BamChunk& pre, BamChunk& cur, const PgbamSidecarData& sidecar);

/**
 * @brief Fallback stitch for adjacent chunks with explicit `.pgbam` evidence thresholds.
 */
bool stitch_adjacent_chunks_with_pgbam(BamChunk& pre,
                                       BamChunk& cur,
                                       const PgbamSidecarData& sidecar,
                                       int min_winning_threads,
                                       int polarity_margin);

} // namespace pgphase_collect

#endif
