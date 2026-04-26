#ifndef PGPHASE_BAM_DIGAR_HPP
#define PGPHASE_BAM_DIGAR_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

bool read_overlaps_prev_region(const RegionChunk& chunk, int read_tid, hts_pos_t read_beg, hts_pos_t read_end);
bool read_overlaps_next_region(const RegionChunk& chunk, int read_tid, hts_pos_t read_beg, hts_pos_t read_end);
std::vector<ReadRecord> load_read_records_for_chunk(const Options& opts,
                                                    const RegionChunk& chunk,
                                                    int input_index,
                                                    SamFile& bam,
                                                    bam_hdr_t* header,
                                                    const hts_idx_t* index,
                                                    ReferenceCache& ref,
                                                    OverlapSkipCounts* overlap_skip_counts = nullptr);
void finalize_bam_chunk(BamChunk& chunk, ReferenceCache& ref, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
