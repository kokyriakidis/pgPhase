#include "region_excludes.hpp"

#include <algorithm>

namespace pgphase_collect {

std::vector<std::pair<hts_pos_t, hts_pos_t>> subtract_excludes(
    const std::string& chrom, hts_pos_t beg, hts_pos_t end,
    const std::vector<RegionFilter>& excludes) {
    std::vector<std::pair<hts_pos_t, hts_pos_t>> clipped;
    for (const RegionFilter& e : excludes) {
        if (e.chrom != chrom) continue;
        const hts_pos_t eb = std::max(beg, e.beg);
        const hts_pos_t ee = std::min(end, e.end);
        if (eb <= ee) clipped.push_back({eb, ee});
    }
    std::sort(clipped.begin(), clipped.end());

    std::vector<std::pair<hts_pos_t, hts_pos_t>> kept;
    hts_pos_t cursor = beg;
    for (const auto& ex : clipped) {
        if (ex.first > cursor) kept.push_back({cursor, ex.first - 1});
        cursor = std::max(cursor, ex.second + 1);
        if (cursor > end) break;
    }
    if (cursor <= end) kept.push_back({cursor, end});
    return kept;
}

} // namespace pgphase_collect
