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
 * @brief Fallback stitch for adjacent chunks when common-read overlap has no decisive signal.
 */
bool stitch_adjacent_chunks_with_pgbam(BamChunk& pre, BamChunk& cur, const PgbamSidecarData& sidecar);

} // namespace pgphase_collect

#endif
