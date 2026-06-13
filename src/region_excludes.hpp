#ifndef PGPHASE_REGION_EXCLUDES_HPP
#define PGPHASE_REGION_EXCLUDES_HPP

// Pure interval arithmetic for --exclude-bed support, factored out so it can be
// unit tested without linking the full pipeline.

#include "phasing_types.hpp"

#include <string>
#include <utility>
#include <vector>

#include <htslib/hts.h>

namespace pgphase_collect {

/**
 * @brief Subtracts exclude intervals from one 1-based inclusive include range.
 *
 * `excludes` may be unsorted and overlapping; only entries on @p chrom that
 * intersect `[beg, end]` are applied.
 *
 * @return The portions of `[beg, end]` that remain after removing every exclude
 *         interval (empty if the range is fully excluded).
 */
std::vector<std::pair<hts_pos_t, hts_pos_t>> subtract_excludes(
    const std::string& chrom, hts_pos_t beg, hts_pos_t end,
    const std::vector<RegionFilter>& excludes);

} // namespace pgphase_collect

#endif
