#include "graph_query.hpp"
#include "gbz_ffi.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

#include <htslib/hts.h>

namespace pgphase_collect {

namespace {

// Compact handle encoding: pack (node_id, orientation) into a uint32_t.
// Bit layout: [31:1] = node_id, [0] = reverse flag.
// This matches the GBWT handle encoding (2*node_id + orientation), so
// FFI-provided GBWT handles can be truncated directly to CompactHandle.
using CompactHandle = uint32_t;
constexpr uint64_t kMaxCompactNodeId =
    static_cast<uint64_t>(std::numeric_limits<CompactHandle>::max() >> 1);

CompactHandle encode_handle(uint64_t node, bool reverse) {
    return static_cast<CompactHandle>((node << 1) |
                                      static_cast<uint64_t>(reverse ? 1 : 0));
}

CompactHandle reverse_handle(CompactHandle handle) {
    return handle ^ 1ULL;
}

bool parse_unsigned_node(std::string_view value, uint64_t& out) {
    if (value.empty()) return false;
    uint64_t v = 0;
    for (char c : value) {
        if (c < '0' || c > '9') return false;
        const uint64_t digit = static_cast<uint64_t>(c - '0');
        if (v > (std::numeric_limits<uint64_t>::max() - digit) / 10ULL) return false;
        v = v * 10ULL + digit;
    }
    out = v;
    return true;
}

// Parse a GAF path field (e.g. ">12>34<56>78") into a vector of CompactHandles.
// Allele walks stored in a single contiguous buffer for cache locality.
// Each allele is described by an offset+length into the shared buffer.
struct FlatAllelePack {
    std::vector<CompactHandle> data;       // contiguous walk data for all alleles
    std::vector<uint32_t>      offsets;    // offsets[i] = start of allele i in data
    std::vector<uint16_t>      lengths;    // lengths[i] = number of handles in allele i
    size_t n_alleles = 0;

    const CompactHandle* allele_ptr(size_t i) const { return data.data() + offsets[i]; }
    size_t allele_len(size_t i) const { return lengths[i]; }
};

// A snarl site ready for numeric matching.  Stores the ref walk's boundary
// handles (left/right and their reverse complements) plus all allele walks
// packed into a FlatAllelePack.
struct CompactGraphSite {
    size_t original_index = 0;   // index into GraphSiteCatalog::sites
    std::string site_id;
    std::string chrom;
    hts_pos_t pos = 0;
    CompactHandle left = 0;      // first handle of ref walk (forward)
    CompactHandle right = 0;     // last handle of ref walk (forward)
    CompactHandle left_rev = 0;  // reverse complement of left
    CompactHandle right_rev = 0; // reverse complement of right
    FlatAllelePack alleles;
    size_t max_walk_len = 0;     // longest allele walk (for early-break)
};

// Maps a boundary handle to the compact site it belongs to.
// All entries are sorted by handle for binary-search lookup.
struct BoundarySiteEntry {
    CompactHandle handle;
    uint32_t site_index;
};

struct CompactGraphSiteIndex {
    std::vector<CompactGraphSite> sites;
    std::vector<BoundarySiteEntry> boundary_entries;  // sorted by handle

    // Find all site indices for a given boundary handle via binary search.
    // Returns a pair of pointers into boundary_entries [begin, end).
    std::pair<const BoundarySiteEntry*, const BoundarySiteEntry*>
    find_boundary(CompactHandle handle) const {
        auto cmp = [](const BoundarySiteEntry& a, const BoundarySiteEntry& b) {
            return a.handle < b.handle;
        };
        BoundarySiteEntry key{handle, 0};
        auto beg = std::lower_bound(boundary_entries.begin(), boundary_entries.end(), key, cmp);
        if (beg == boundary_entries.end() || beg->handle != handle)
            return {nullptr, nullptr};
        auto end = std::upper_bound(beg, boundary_entries.end(), key, cmp);
        const BoundarySiteEntry* base = boundary_entries.data();
        return {base + (beg - boundary_entries.begin()),
                base + (end - boundary_entries.begin())};
    }
};

void add_boundary_entry(std::vector<BoundarySiteEntry>& entries,
                        CompactHandle handle, uint32_t site_index) {
    entries.push_back({handle, site_index});
}

// Convert a GraphSiteCatalog into a CompactGraphSiteIndex: flatten allele
// walks into contiguous FlatAllelePacks and build a sorted boundary→site
// index for O(log n) handle lookup during read scanning.
bool build_compact_graph_site_index(const GraphSiteCatalog& catalog,
                                    CompactGraphSiteIndex& out) {
    out = {};
    out.sites.reserve(catalog.sites.size());
    // Pre-allocate boundary entries: ~4 entries per site (left, right, left_rev, right_rev).
    out.boundary_entries.reserve(catalog.sites.size() * 4);

    for (size_t i = 0; i < catalog.sites.size(); ++i) {
        const GraphSite& site = catalog.sites[i];
        if (!graph_site_is_queryable(site)) continue;
        if (site.allele_walks.empty() || site.allele_walks[0].size() < 2) continue;

        CompactGraphSite compact;
        compact.original_index = i;
        compact.site_id = graph_site_key_str(site);
        compact.chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
        compact.pos = site.pos;

        // Flatten allele walks into a contiguous buffer.
        FlatAllelePack& pack = compact.alleles;
        pack.n_alleles = site.allele_walks.size();
        pack.offsets.resize(pack.n_alleles);
        pack.lengths.resize(pack.n_alleles);
        size_t total_handles = 0;
        for (const GraphWalk& walk : site.allele_walks) total_handles += walk.size();
        pack.data.reserve(total_handles);

        for (size_t ai = 0; ai < site.allele_walks.size(); ++ai) {
            const GraphWalk& walk = site.allele_walks[ai];
            pack.offsets[ai] = static_cast<uint32_t>(pack.data.size());
            if (walk.empty()) {
                // Spanning deletion (*) — zero-length entry, never matches.
                pack.lengths[ai] = 0;
                continue;
            }
            const size_t walk_start = pack.data.size();
            for (const GraphWalkStep& step : walk) {
                uint64_t node = 0;
                if (!parse_unsigned_node(step.node, node)) return false;
                if (node > kMaxCompactNodeId) return false;
                pack.data.push_back(encode_handle(node, step.reverse));
            }
            const size_t walk_len = pack.data.size() - walk_start;
            if (walk_len == 0) return false;
            if (walk_len > std::numeric_limits<uint16_t>::max()) return false;
            pack.lengths[ai] = static_cast<uint16_t>(walk_len);
            compact.max_walk_len = std::max(compact.max_walk_len, walk_len);
        }

        // Ref walk (allele 0) boundaries.
        compact.left = pack.data[pack.offsets[0]];
        compact.right = pack.data[pack.offsets[0] + pack.lengths[0] - 1];
        compact.left_rev = reverse_handle(compact.left);
        compact.right_rev = reverse_handle(compact.right);

        const auto compact_index = static_cast<uint32_t>(out.sites.size());
        out.sites.push_back(std::move(compact));
        add_boundary_entry(out.boundary_entries, out.sites.back().left, compact_index);
        add_boundary_entry(out.boundary_entries, out.sites.back().right, compact_index);
        add_boundary_entry(out.boundary_entries, out.sites.back().left_rev, compact_index);
        add_boundary_entry(out.boundary_entries, out.sites.back().right_rev, compact_index);
    }

    // Sort boundary entries by handle for binary search.
    std::sort(out.boundary_entries.begin(), out.boundary_entries.end(),
              [](const BoundarySiteEntry& a, const BoundarySiteEntry& b) {
                  return a.handle < b.handle ||
                         (a.handle == b.handle && a.site_index < b.site_index);
              });
    // Deduplicate (same handle+site_index).
    out.boundary_entries.erase(
        std::unique(out.boundary_entries.begin(), out.boundary_entries.end(),
                    [](const BoundarySiteEntry& a, const BoundarySiteEntry& b) {
                        return a.handle == b.handle && a.site_index == b.site_index;
                    }),
        out.boundary_entries.end());

    return !out.sites.empty();
}

// Per-read boundary position index: flat sorted vector of (handle, position).
// Replaces unordered_map<CompactHandle, vector<size_t>> to avoid per-read
// heap allocation. Sorted by handle for binary-search lookup.
struct FlatBoundaryPositions {
    struct Entry { CompactHandle handle; uint32_t pos; };
    std::vector<Entry> entries;

    void clear() { entries.clear(); }
    void add(CompactHandle h, uint32_t pos) { entries.push_back({h, pos}); }
    void sort() {
        std::sort(entries.begin(), entries.end(),
                  [](const Entry& a, const Entry& b) {
                      if (a.handle != b.handle) return a.handle < b.handle;
                      return a.pos < b.pos;
                  });
    }

    // Returns [begin, end) pointers for all positions with the given handle.
    std::pair<const Entry*, const Entry*> find(CompactHandle h) const {
        auto cmp = [](const Entry& a, const Entry& b) { return a.handle < b.handle; };
        Entry key{h, 0};
        auto beg = std::lower_bound(entries.begin(), entries.end(), key, cmp);
        if (beg == entries.end() || beg->handle != h) return {nullptr, nullptr};
        auto end = std::upper_bound(beg, entries.end(), key, cmp);
        const Entry* base = entries.data();
        return {base + (beg - entries.begin()),
                base + (end - entries.begin())};
    }
};

// Compare a read sub-walk against one allele walk.  Forward comparison uses
// memcmp; reverse comparison walks the allele backwards, flipping orientation.
bool compact_span_matches_allele(const CompactHandle* read_ptr, size_t span_len,
                                 const CompactHandle* allele_ptr, size_t allele_len,
                                 bool reverse) {
    if (span_len != allele_len || allele_len == 0) return false;
    if (!reverse) {
        return std::memcmp(read_ptr, allele_ptr, span_len * sizeof(CompactHandle)) == 0;
    }
    for (size_t i = 0; i < span_len; ++i) {
        if (read_ptr[i] != reverse_handle(allele_ptr[span_len - 1 - i])) return false;
    }
    return true;
}

// Try to match a read span against all allele walks of a site.
// Returns the allele index (≥0), kGraphAlleleMissing, or kGraphAlleleAmbiguous
// if multiple alleles match (shouldn't happen with well-formed snarls).
int match_compact_span(const CompactHandle* read_ptr, size_t span_len,
                       const CompactGraphSite& site,
                       bool reverse) {
    int match = kGraphAlleleMissing;
    const FlatAllelePack& pack = site.alleles;
    for (size_t allele_i = 0; allele_i < pack.n_alleles; ++allele_i) {
        if (!compact_span_matches_allele(read_ptr, span_len,
                                         pack.allele_ptr(allele_i),
                                         pack.allele_len(allele_i), reverse)) {
            continue;
        }
        if (match != kGraphAlleleMissing) return kGraphAlleleAmbiguous;
        match = static_cast<int>(allele_i);
    }
    return match;
}

// Determine which allele a read traverses at a given site.
// Looks up the site's boundary handles in the read's position index, extracts
// the sub-walk between left and right boundaries, and compares against each
// allele walk.  Tries forward orientation first, then reverse complement.
int match_compact_site_on_read(
    const std::vector<CompactHandle>& read_walk,
    const FlatBoundaryPositions& boundary_positions,
    const CompactGraphSite& site,
    bool* reverse_out = nullptr) {
    auto [left_beg, left_end] = boundary_positions.find(site.left);
    auto [right_beg, right_end] = boundary_positions.find(site.right);
    if (left_beg && right_beg) {
        for (auto lp = left_beg; lp != left_end; ++lp) {
            for (auto rp = right_beg; rp != right_end; ++rp) {
                if (rp->pos < lp->pos) continue;
                const size_t span = rp->pos - lp->pos + 1;
                const int allele = match_compact_span(
                    read_walk.data() + lp->pos, span, site, false);
                if (allele >= 0) { if (reverse_out) *reverse_out = false; return allele; }
                if (allele == kGraphAlleleAmbiguous) return allele;
                if (site.max_walk_len > 0 && span > site.max_walk_len) break;
            }
        }
    }

    auto [rrev_beg, rrev_end] = boundary_positions.find(site.right_rev);
    auto [lrev_beg, lrev_end] = boundary_positions.find(site.left_rev);
    if (rrev_beg && lrev_beg) {
        for (auto rp = rrev_beg; rp != rrev_end; ++rp) {
            for (auto lp = lrev_beg; lp != lrev_end; ++lp) {
                if (lp->pos < rp->pos) continue;
                const size_t span = lp->pos - rp->pos + 1;
                const int allele = match_compact_span(
                    read_walk.data() + rp->pos, span, site, true);
                if (allele >= 0) { if (reverse_out) *reverse_out = true; return allele; }
                if (allele == kGraphAlleleAmbiguous) return allele;
                if (site.max_walk_len > 0 && span > site.max_walk_len) break;
            }
        }
    }

    return kGraphAlleleMissing;
}

// Callback context for the Rust FFI path: receives pre-parsed GBWT handle
// arrays instead of GAF text.  GBWT handles (2*node_id + orientation) share
// the same encoding as CompactHandle, so no conversion is needed.
struct CompactFFIContext {
    const CompactGraphSiteIndex* index;
    int min_mapq;
    std::vector<GraphReadAllele>* results;
    // Per-callback scratch buffers (avoid repeated allocation).
    std::vector<CompactHandle> read_walk;
    FlatBoundaryPositions boundary_positions;
    std::vector<size_t> candidates;
    std::vector<uint32_t> site_marks;
    uint32_t mark_stamp = 0;
};

extern "C" void structured_alignment_callback(
    const unsigned char* name, size_t name_len,
    const uint64_t* nodes, size_t node_count,
    int mapq, void* user_data)
{
    auto* ctx = static_cast<CompactFFIContext*>(user_data);
    if (mapq < ctx->min_mapq) return;
    if (node_count == 0) return;

    const CompactGraphSiteIndex& index = *ctx->index;

    // Convert GBWT handles to CompactHandles (same encoding, just truncate).
    // Skip the entire read if any node exceeds the CompactHandle range to
    // avoid corrupting the walk (matches parse_gaf_path_compact behavior).
    ctx->read_walk.clear();
    ctx->read_walk.reserve(node_count);
    for (size_t i = 0; i < node_count; ++i) {
        if (nodes[i] / 2 > kMaxCompactNodeId) {
            ctx->read_walk.clear();
            return;
        }
        ctx->read_walk.push_back(static_cast<CompactHandle>(nodes[i]));
    }
    if (ctx->read_walk.empty()) return;

    // Find candidate sites via boundary handle lookup.
    ctx->candidates.clear();
    ctx->boundary_positions.clear();
    if (++(ctx->mark_stamp) == 0) {
        std::fill(ctx->site_marks.begin(), ctx->site_marks.end(), 0);
        ctx->mark_stamp = 1;
    }

    for (size_t step_i = 0; step_i < ctx->read_walk.size(); ++step_i) {
        const CompactHandle handle = ctx->read_walk[step_i];
        auto [beg, end] = index.find_boundary(handle);
        if (!beg) continue;
        ctx->boundary_positions.add(handle, static_cast<uint32_t>(step_i));
        for (auto it = beg; it != end; ++it) {
            if (ctx->site_marks[it->site_index] == ctx->mark_stamp) continue;
            ctx->site_marks[it->site_index] = ctx->mark_stamp;
            ctx->candidates.push_back(it->site_index);
        }
    }
    if (ctx->candidates.empty()) return;
    ctx->boundary_positions.sort();

    const std::string read_name(reinterpret_cast<const char*>(name), name_len);

    for (size_t compact_site_index : ctx->candidates) {
        const CompactGraphSite& site = index.sites[compact_site_index];
        bool rev = false;
        const int allele = match_compact_site_on_read(
            ctx->read_walk, ctx->boundary_positions, site, &rev);
        if (allele < 0) continue;

        ctx->results->push_back(GraphReadAllele{
            site.site_id,
            site.chrom,
            site.pos,
            read_name,
            allele,
            mapq,
            rev
        });
    }
}

} // namespace

std::vector<GraphReadAllele>
query_gbz_interval_gaf_ffi(void* gbz_handle,
                            void* gaf_handle,
                            const std::string& sample,
                            const std::string& contig,
                            hts_pos_t beg,
                            hts_pos_t end,
                            const GraphSiteCatalog& catalog,
                            int min_mapq)
{
    if (!gbz_handle) throw std::runtime_error("null GBZ handle for FFI query");
    if (!gaf_handle) throw std::runtime_error("null GAF handle for FFI query");
    if (contig.empty() || end <= beg)
        throw std::runtime_error("invalid interval for graph query: " + contig +
                                 ":" + std::to_string(beg) + ".." + std::to_string(end));

    // Build compact numeric site index from catalog walks.
    CompactGraphSiteIndex compact_index;
    if (!build_compact_graph_site_index(catalog, compact_index)) {
        // Fallback: catalog has no compactable sites.
        return {};
    }

    std::vector<GraphReadAllele> results;
    CompactFFIContext ctx{};
    ctx.index = &compact_index;
    ctx.min_mapq = min_mapq;
    ctx.results = &results;
    ctx.site_marks.resize(compact_index.sites.size(), 0);

    char* err = nullptr;
    int rc = pgphase_gbz_query_interval_structured(
        gbz_handle, gaf_handle,
        sample.empty() ? nullptr : sample.c_str(),
        contig.c_str(),
        static_cast<uint64_t>(beg),
        static_cast<uint64_t>(end),
        structured_alignment_callback,
        &ctx,
        &err);
    if (rc != 0) {
        std::string msg = err ? std::string(err) : "unknown FFI error";
        if (err) pgphase_gbz_free_string(err);
        throw std::runtime_error("gbz_query_interval failed: " + msg);
    }

    return results;
}

PathRange
query_gbz_path_range(void* gbz_handle,
                     const std::string& sample,
                     const std::string& contig)
{
    PathRange result;
    if (!gbz_handle) return result;

    uint64_t out_start = 0, out_end = 0;
    char* err = nullptr;
    int rc = pgphase_gbz_path_range(
        gbz_handle,
        sample.empty() ? nullptr : sample.c_str(),
        contig.c_str(),
        &out_start, &out_end, &err);
    if (rc != 0) {
        if (err) pgphase_gbz_free_string(err);
        return result;
    }

    result.start = static_cast<hts_pos_t>(out_start);
    result.end   = static_cast<hts_pos_t>(out_end);
    result.valid  = true;
    return result;
}

} // namespace pgphase_collect
