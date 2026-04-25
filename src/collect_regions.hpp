#ifndef PGPHASE_COLLECT_REGIONS_HPP
#define PGPHASE_COLLECT_REGIONS_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

RegionFilter parse_region(const std::string& region);
std::vector<RegionFilter> load_bed_regions(const std::string& path);
std::vector<RegionChunk> build_region_chunks(const Options& opts, const bam_hdr_t* header, const faidx_t* fai);
std::vector<RegionChunk> load_region_chunks(const Options& opts);

} // namespace pgphase_collect

#endif
