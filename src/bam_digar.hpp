#ifndef PGPHASE_BAM_DIGAR_HPP
#define PGPHASE_BAM_DIGAR_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

std::vector<ReadRecord> load_read_records_for_chunk(const Options& opts,
                                                    const RegionChunk& chunk,
                                                    SamFile& bam,
                                                    bam_hdr_t* header,
                                                    const hts_idx_t* index,
                                                    ReferenceCache& ref);
void finalize_bam_chunk(BamChunk& chunk, ReferenceCache& ref, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
