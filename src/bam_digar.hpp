#ifndef PGPHASE_BAM_DIGAR_HPP
#define PGPHASE_BAM_DIGAR_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

/** Calculates bounds logic to drop alignments intersecting previously parsed map regions, avoiding processing overlapping duplicates. */
bool read_overlaps_prev_region(const RegionChunk& chunk, int read_tid, hts_pos_t read_beg, hts_pos_t read_end);

/** Identifies if a read bleeds into the boundary space of the chunk physically next in the alignment timeline. */
bool read_overlaps_next_region(const RegionChunk& chunk, int read_tid, hts_pos_t read_beg, hts_pos_t read_end);

/**
 * @brief Memory-centric driver orchestrating parsing alignments mapped within an isolation block.
 * 
 * Maps directly to `longcallD` sequence loading routines mapping alignments across threading targets using HTSlib queries. 
 * Performs CIGAR mapping into variant states via `build_digars_*`.
 */
std::vector<ReadRecord> load_read_records_for_chunk(const Options& opts,
                                                    const RegionChunk& chunk,
                                                    int input_index,
                                                    SamFile& bam,
                                                    bam_hdr_t* header,
                                                    const hts_idx_t* index,
                                                    ReferenceCache& ref,
                                                    OverlapSkipCounts* overlap_skip_counts = nullptr);

/** Executes cleanups closing off parsing state. Triggers FASTA cache bounds constraints mapped against contig lengths. */
void finalize_bam_chunk(BamChunk& chunk, ReferenceCache& ref, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
