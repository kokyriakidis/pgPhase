#include "region_excludes.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace pgphase_collect;

namespace {

bool check(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

using Range = std::pair<hts_pos_t, hts_pos_t>;

// Build an exclude filter list from {chrom, beg, end} triples (1-based inclusive).
std::vector<RegionFilter> ex(std::vector<RegionFilter> v) { return v; }

bool eq(const std::vector<Range>& got, const std::vector<Range>& want) {
    if (got.size() != want.size()) return false;
    for (size_t i = 0; i < got.size(); ++i)
        if (got[i] != want[i]) return false;
    return true;
}

} // namespace

int main() {
    bool ok = true;

    // No excludes: whole range survives.
    ok &= check(eq(subtract_excludes("chr1", 1, 100, {}),
                   {{1, 100}}),
                "no excludes keeps full range");

    // Exclude on a different contig is ignored.
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr2", 10, 20}})),
                   {{1, 100}}),
                "exclude on other contig ignored");

    // Middle cut splits into two ranges.
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr1", 30, 40}})),
                   {{1, 29}, {41, 100}}),
                "middle exclude splits range");

    // Left-edge cut advances begin.
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr1", 1, 10}})),
                   {{11, 100}}),
                "left-edge exclude advances begin");

    // Right-edge cut shortens end.
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr1", 90, 100}})),
                   {{1, 89}}),
                "right-edge exclude shortens end");

    // Fully covering exclude removes everything.
    ok &= check(eq(subtract_excludes("chr1", 10, 50,
                       ex({{true, "chr1", 1, 100}})),
                   {}),
                "covering exclude removes range");

    // Exclude exactly equal to the range removes everything.
    ok &= check(eq(subtract_excludes("chr1", 10, 50,
                       ex({{true, "chr1", 10, 50}})),
                   {}),
                "exact exclude removes range");

    // Overlapping excludes are merged (no negative/empty fragments).
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr1", 30, 50}, {true, "chr1", 40, 60}})),
                   {{1, 29}, {61, 100}}),
                "overlapping excludes merge");

    // Unsorted excludes handled (same expected result as sorted).
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr1", 70, 80}, {true, "chr1", 20, 30}})),
                   {{1, 19}, {31, 69}, {81, 100}}),
                "unsorted excludes handled");

    // Adjacent (touching) excludes leave a single internal gap removed.
    ok &= check(eq(subtract_excludes("chr1", 1, 100,
                       ex({{true, "chr1", 30, 39}, {true, "chr1", 40, 49}})),
                   {{1, 29}, {50, 100}}),
                "adjacent excludes leave one gap");

    // Exclude entirely outside the range is a no-op.
    ok &= check(eq(subtract_excludes("chr1", 50, 100,
                       ex({{true, "chr1", 1, 20}})),
                   {{50, 100}}),
                "exclude outside range is no-op");

    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
