#include "graph_query.hpp"
#include "gbz_ffi.h"

#include <array>
#include <cctype>
#include <cstdio>
#include <functional>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

namespace pgphase_collect {

using GraphReadAlleleThreadEmitter = std::function<void(size_t, GraphReadAllele&&)>;

namespace {

std::string graph_site_key_str(const GraphSite& site, size_t /*index*/) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
}

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

bool graph_walk_to_compact_handles(const GraphWalk& walk,
                                   std::vector<CompactHandle>& out) {
    out.clear();
    out.reserve(walk.size());
    for (const GraphWalkStep& step : walk) {
        uint64_t node = 0;
        if (!parse_unsigned_node(step.node, node)) return false;
        if (node > kMaxCompactNodeId) return false;
        out.push_back(encode_handle(node, step.reverse));
    }
    return true;
}

bool parse_gaf_path_compact(std::string_view walk,
                            std::vector<CompactHandle>& out) {
    out.clear();
    size_t i = 0;
    while (i < walk.size()) {
        const char orient = walk[i];
        if (orient != '>' && orient != '<') return false;
        const bool reverse = orient == '<';
        ++i;
        const size_t beg = i;
        uint64_t node = 0;
        while (i < walk.size() && walk[i] != '>' && walk[i] != '<') {
            const char c = walk[i];
            if (c < '0' || c > '9') return false;
            const uint64_t digit = static_cast<uint64_t>(c - '0');
            if (node > (std::numeric_limits<uint64_t>::max() - digit) / 10ULL) return false;
            node = node * 10ULL + digit;
            if (node > kMaxCompactNodeId) return false;
            ++i;
        }
        if (beg == i) return false;
        out.push_back(encode_handle(node, reverse));
    }
    return true;
}

struct GafCoreFields {
    std::string_view read_name;
    std::string_view walk;
    int mapq = 0;
};

bool parse_gaf_core_fields_from_column(std::string_view line,
                                       int first_gaf_column,
                                       GafCoreFields& out);

bool parse_gaf_core_fields_from_column(std::string_view line,
                                       int first_gaf_column,
                                       GafCoreFields& out) {
    std::array<std::pair<const char*, size_t>, 12> fld{};
    int fi = -first_gaf_column;
    const char* fs = line.data();
    const char* le = fs + line.size();
    for (const char* p = fs; p <= le && fi < 12; ++p) {
        if (p == le || *p == '\t') {
            if (fi >= 0) fld[fi] = {fs, static_cast<size_t>(p - fs)};
            ++fi;
            fs = p + 1;
        }
    }
    if (fi < 12) return false;

    int mapq = 0;
    const char* mp = fld[11].first;
    const char* me = mp + fld[11].second;
    if (mp == me) return false;
    for (; mp < me; ++mp) {
        if (*mp < '0' || *mp > '9') return false;
        mapq = mapq * 10 + (*mp - '0');
    }
    out.read_name = std::string_view(fld[0].first, fld[0].second);
    out.walk = std::string_view(fld[5].first, fld[5].second);
    out.mapq = mapq;
    return true;
}

struct CompactGraphSite {
    size_t original_index = 0;
    std::string site_id;
    std::string chrom;
    hts_pos_t pos = 0;
    CompactHandle left = 0;
    CompactHandle right = 0;
    CompactHandle left_rev = 0;
    CompactHandle right_rev = 0;
    std::vector<std::vector<CompactHandle>> allele_walks;
};

struct CompactGraphSiteIndex {
    std::vector<CompactGraphSite> sites;
    std::unordered_map<CompactHandle, std::vector<size_t>> boundary_to_sites;
};

bool add_boundary_site(std::unordered_map<CompactHandle, std::vector<size_t>>& boundary_to_sites,
                       CompactHandle handle,
                       size_t site_index) {
    std::vector<size_t>& sites = boundary_to_sites[handle];
    if (sites.empty() || sites.back() != site_index) sites.push_back(site_index);
    return true;
}

bool build_compact_graph_site_index(GraphSiteCatalog& catalog,
                                    CompactGraphSiteIndex& out,
                                    bool release_catalog_walks) {
    out = {};
    out.sites.reserve(catalog.sites.size());
    for (size_t i = 0; i < catalog.sites.size(); ++i) {
        GraphSite& site = catalog.sites[i];
        if (!graph_site_is_queryable(site)) continue;
        if (site.allele_walks.empty() || site.allele_walks[0].size() < 2) continue;

        CompactGraphSite compact;
        compact.original_index = i;
        compact.site_id = graph_site_key_str(site, i);
        compact.chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
        compact.pos = site.pos;
        compact.allele_walks.reserve(site.allele_walks.size());
        for (const GraphWalk& walk : site.allele_walks) {
            std::vector<CompactHandle> handles;
            if (!graph_walk_to_compact_handles(walk, handles)) return false;
            if (handles.empty()) return false;
            compact.allele_walks.push_back(std::move(handles));
        }
        compact.left = compact.allele_walks[0].front();
        compact.right = compact.allele_walks[0].back();
        compact.left_rev = reverse_handle(compact.left);
        compact.right_rev = reverse_handle(compact.right);

        const size_t compact_index = out.sites.size();
        out.sites.push_back(std::move(compact));
        add_boundary_site(out.boundary_to_sites, out.sites.back().left, compact_index);
        add_boundary_site(out.boundary_to_sites, out.sites.back().right, compact_index);
        add_boundary_site(out.boundary_to_sites, out.sites.back().left_rev, compact_index);
        add_boundary_site(out.boundary_to_sites, out.sites.back().right_rev, compact_index);

        if (release_catalog_walks) {
            const size_t allele_count = site.allele_walks.size();
            std::vector<GraphWalk> released(allele_count);
            site.allele_walks.swap(released);
            std::vector<std::string>().swap(site.allele_traversals);
        }
    }
    return !out.sites.empty();
}

bool compact_span_matches_allele(const std::vector<CompactHandle>& read_walk,
                                 size_t beg,
                                 size_t end,
                                 const std::vector<CompactHandle>& allele,
                                 bool reverse) {
    if (end < beg) return false;
    const size_t len = end - beg + 1;
    if (len != allele.size()) return false;
    if (!reverse) {
        for (size_t i = 0; i < len; ++i) {
            if (read_walk[beg + i] != allele[i]) return false;
        }
    } else {
        for (size_t i = 0; i < len; ++i) {
            if (read_walk[beg + i] != reverse_handle(allele[len - 1 - i])) return false;
        }
    }
    return true;
}

int match_compact_span(const std::vector<CompactHandle>& read_walk,
                       size_t beg,
                       size_t end,
                       const CompactGraphSite& site,
                       bool reverse) {
    int match = kGraphAlleleMissing;
    for (size_t allele_i = 0; allele_i < site.allele_walks.size(); ++allele_i) {
        if (!compact_span_matches_allele(read_walk, beg, end,
                                         site.allele_walks[allele_i], reverse)) {
            continue;
        }
        if (match != kGraphAlleleMissing) return kGraphAlleleAmbiguous;
        match = static_cast<int>(allele_i);
    }
    return match;
}

int match_compact_site_on_read(
    const std::vector<CompactHandle>& read_walk,
    const std::unordered_map<CompactHandle, std::vector<size_t>>& boundary_positions,
    const CompactGraphSite& site,
    bool* reverse_out = nullptr) {
    auto left_it = boundary_positions.find(site.left);
    auto right_it = boundary_positions.find(site.right);
    if (left_it != boundary_positions.end() && right_it != boundary_positions.end()) {
        for (size_t left_pos : left_it->second) {
            for (size_t right_pos : right_it->second) {
                if (right_pos < left_pos) continue;
                const int allele = match_compact_span(read_walk, left_pos, right_pos,
                                                      site, false);
                if (allele >= 0) { if (reverse_out) *reverse_out = false; return allele; }
                if (allele == kGraphAlleleAmbiguous) return allele;
                if (!site.allele_walks.empty() &&
                    right_pos - left_pos + 1 > site.allele_walks[0].size()) {
                    break;
                }
            }
        }
    }

    auto right_rev_it = boundary_positions.find(site.right_rev);
    auto left_rev_it = boundary_positions.find(site.left_rev);
    if (right_rev_it != boundary_positions.end() && left_rev_it != boundary_positions.end()) {
        for (size_t right_pos : right_rev_it->second) {
            for (size_t left_pos : left_rev_it->second) {
                if (left_pos < right_pos) continue;
                const int allele = match_compact_span(read_walk, right_pos, left_pos,
                                                      site, true);
                if (allele >= 0) { if (reverse_out) *reverse_out = true; return allele; }
                if (allele == kGraphAlleleAmbiguous) return allele;
                if (!site.allele_walks.empty() &&
                    left_pos - right_pos + 1 > site.allele_walks[0].size()) {
                    break;
                }
            }
        }
    }

    return kGraphAlleleMissing;
}

size_t scan_gaf_line_compact(std::string_view line,
                             int first_gaf_column,
                             const CompactGraphSiteIndex& compact_index,
                             int min_mapq,
                             size_t worker_id,
                             std::vector<CompactHandle>& read_walk,
                             std::unordered_map<CompactHandle, std::vector<size_t>>& boundary_positions,
                             std::vector<size_t>& candidates,
                             std::vector<uint32_t>& site_marks,
                             uint32_t& mark_stamp,
                             const GraphReadAlleleThreadEmitter& emit) {
    if (line.empty() || line[0] == '#') return 0;

    GafCoreFields fields;
    if (!parse_gaf_core_fields_from_column(line, first_gaf_column, fields)) return 0;
    if (fields.mapq < min_mapq) return 0;
    if (!parse_gaf_path_compact(fields.walk, read_walk) || read_walk.empty()) return 0;

    candidates.clear();
    boundary_positions.clear();
    if (++mark_stamp == 0) {
        std::fill(site_marks.begin(), site_marks.end(), 0);
        mark_stamp = 1;
    }

    for (size_t step_i = 0; step_i < read_walk.size(); ++step_i) {
        const CompactHandle handle = read_walk[step_i];
        auto site_it = compact_index.boundary_to_sites.find(handle);
        if (site_it == compact_index.boundary_to_sites.end()) continue;
        boundary_positions[handle].push_back(step_i);
        for (size_t site_index : site_it->second) {
            if (site_marks[site_index] == mark_stamp) continue;
            site_marks[site_index] = mark_stamp;
            candidates.push_back(site_index);
        }
    }
    if (candidates.empty()) return 0;

    size_t emitted = 0;
    const std::string read_name(fields.read_name);
    for (size_t compact_site_index : candidates) {
        const CompactGraphSite& site = compact_index.sites[compact_site_index];
        bool rev = false;
        const int allele = match_compact_site_on_read(read_walk, boundary_positions, site, &rev);
        if (allele < 0) continue;
        emit(worker_id, GraphReadAllele{
            site.site_id,
            site.chrom,
            site.pos,
            read_name,
            allele,
            "",
            fields.mapq,
            rev
        });
        ++emitted;
    }
    return emitted;
}

// Returns the sub-walk of `walk` that spans between left_bnd and right_bnd
// (inclusive), in whichever orientation the read traverses the snarl.
// Returns empty if the read does not span both boundaries.
// match_graph_allele_exact handles reverse-complement internally, so we can
// pass the raw sub-walk without reversing it first.
GraphWalk extract_subwalk_between(const GraphWalk& walk,
                                   const GraphWalkStep& left_bnd,
                                   const GraphWalkStep& right_bnd) {
    // Forward orientation: find left_bnd then right_bnd
    for (size_t i = 0; i < walk.size(); ++i) {
        if (!(walk[i] == left_bnd)) continue;
        for (size_t j = i; j < walk.size(); ++j) {
            if (walk[j] == right_bnd)
                return GraphWalk(walk.begin() + static_cast<std::ptrdiff_t>(i),
                                 walk.begin() + static_cast<std::ptrdiff_t>(j) + 1);
        }
    }
    // Reverse orientation: read traverses snarl right-to-left
    const GraphWalkStep right_rev{right_bnd.node, !right_bnd.reverse};
    const GraphWalkStep left_rev{left_bnd.node, !left_bnd.reverse};
    for (size_t i = 0; i < walk.size(); ++i) {
        if (!(walk[i] == right_rev)) continue;
        for (size_t j = i; j < walk.size(); ++j) {
            if (walk[j] == left_rev)
                return GraphWalk(walk.begin() + static_cast<std::ptrdiff_t>(i),
                                 walk.begin() + static_cast<std::ptrdiff_t>(j) + 1);
        }
    }
    return {};
}

} // namespace

// pggaf coordinate-indexed GAF tabix often names sequences by reference suffix ("chr20")
// while graph-site VCFs use pangenome paths ("CHM13#0#chr20"). Mirrors graph_sites
// chrom_suffix matching so tabix queries succeed without a separate FAI map.
static int tbx_seq_tid_with_pangenome_fallback(tbx_t* tbx, const std::string& contig) {
    int tid = tbx_name2id(tbx, contig.c_str());
    if (tid >= 0) return tid;
    const size_t h = contig.rfind('#');
    if (h == std::string::npos) return -1;
    return tbx_name2id(tbx, contig.substr(h + 1).c_str());
}

void require_indexed_gaf(const std::string& indexed_gaf_file) {
    htsFile* fp = hts_open(indexed_gaf_file.c_str(), "r");
    if (fp == nullptr) {
        throw std::runtime_error("failed to open indexed GAF: " + indexed_gaf_file);
    }
    hts_close(fp);

    tbx_t* tbx = tbx_index_load(indexed_gaf_file.c_str());
    if (tbx == nullptr) {
        throw std::runtime_error(
            "indexed GAF requires a tabix index (.tbi): " + indexed_gaf_file);
    }
    tbx_destroy(tbx);
}

std::vector<GraphReadAllele>
scan_indexed_gaf_chunk(const std::string& indexed_gaf_file,
                       const std::string& contig,
                       hts_pos_t beg,
                       hts_pos_t end,
                       const GraphSiteCatalog& catalog,
                       int min_mapq) {
    if (contig.empty() || end <= beg) return {};
    if (catalog.sites.empty()) return {};

    GraphSiteCatalog& mutable_catalog = const_cast<GraphSiteCatalog&>(catalog);
    CompactGraphSiteIndex compact_index;
    if (!build_compact_graph_site_index(mutable_catalog, compact_index, false)) {
        throw std::runtime_error(
            "indexed GAF chunk scan requires queryable graph sites with numeric node IDs");
    }

    htsFile* raw_fp = hts_open(indexed_gaf_file.c_str(), "r");
    if (raw_fp == nullptr) {
        throw std::runtime_error("failed to open indexed GAF: " + indexed_gaf_file);
    }
    struct HtsFileGuard {
        htsFile* fp = nullptr;
        ~HtsFileGuard() { if (fp != nullptr) hts_close(fp); }
    } fp_guard{raw_fp};

    tbx_t* raw_tbx = tbx_index_load(indexed_gaf_file.c_str());
    if (raw_tbx == nullptr) {
        throw std::runtime_error(
            "indexed GAF requires a tabix index (.tbi): " + indexed_gaf_file);
    }
    struct TbxGuard {
        tbx_t* tbx = nullptr;
        ~TbxGuard() { if (tbx != nullptr) tbx_destroy(tbx); }
    } tbx_guard{raw_tbx};

    const int tid = tbx_seq_tid_with_pangenome_fallback(raw_tbx, contig);
    if (tid < 0) return {};

    hts_itr_t* raw_itr = tbx_itr_queryi(raw_tbx, tid, beg, end);
    if (raw_itr == nullptr) return {};
    struct IterGuard {
        hts_itr_t* itr = nullptr;
        ~IterGuard() { if (itr != nullptr) hts_itr_destroy(itr); }
    } itr_guard{raw_itr};

    std::vector<GraphReadAllele> rows;
    uint32_t mark_stamp = 1;
    std::vector<uint32_t> site_marks(compact_index.sites.size(), 0);
    std::vector<size_t> candidates;
    std::vector<CompactHandle> read_walk;
    std::unordered_map<CompactHandle, std::vector<size_t>> boundary_positions;

    kstring_t line = KS_INITIALIZE;
    struct KStringGuard {
        kstring_t* line = nullptr;
        ~KStringGuard() { if (line != nullptr) ks_free(line); }
    } line_guard{&line};

    while (tbx_itr_next(raw_fp, raw_tbx, raw_itr, &line) >= 0) {
        scan_gaf_line_compact(
            std::string_view(line.s, line.l), 3, compact_index, min_mapq, 0,
            read_walk, boundary_positions, candidates, site_marks, mark_stamp,
            [&](size_t, GraphReadAllele&& row) { rows.push_back(std::move(row)); });
    }

    return rows;
}

// ── FFI-based interval query (structured, no GAF text) ────────────────────

namespace {

// Convert a GBWT node handle to a GraphWalkStep.
// GBWT encoding: handle = 2 * node_id + orientation (Forward=0, Reverse=1).
inline GraphWalkStep gbwt_handle_to_step(uint64_t handle) {
    return GraphWalkStep{std::to_string(handle / 2), (handle & 1) != 0};
}

// Context for the structured alignment callback.
struct StructuredQueryContext {
    const std::unordered_map<std::string, std::vector<size_t>>* node_to_sites;
    const GraphSiteCatalog* catalog;
    int min_mapq;
    std::vector<GraphReadAllele>* results;
};

extern "C" void structured_alignment_callback(
    const unsigned char* name, size_t name_len,
    const uint64_t* nodes, size_t node_count,
    int mapq, void* user_data)
{
    auto* ctx = static_cast<StructuredQueryContext*>(user_data);
    if (mapq < ctx->min_mapq) return;
    if (node_count == 0) return;

    // Build GraphWalk from GBWT handles.
    GraphWalk read_walk;
    read_walk.reserve(node_count);
    for (size_t i = 0; i < node_count; ++i)
        read_walk.push_back(gbwt_handle_to_step(nodes[i]));

    // Find candidate sites whose boundary node appears in this walk.
    std::unordered_set<size_t> candidates;
    for (const GraphWalkStep& step : read_walk) {
        auto it = ctx->node_to_sites->find(step.node);
        if (it == ctx->node_to_sites->end()) continue;
        for (size_t si : it->second) candidates.insert(si);
    }
    if (candidates.empty()) return;

    const std::string read_name(reinterpret_cast<const char*>(name), name_len);

    for (size_t si : candidates) {
        const GraphSite& site = ctx->catalog->sites[si];
        const GraphWalkStep& left_bnd  = site.allele_walks[0].front();
        const GraphWalkStep& right_bnd = site.allele_walks[0].back();

        const GraphWalk subwalk = extract_subwalk_between(read_walk, left_bnd, right_bnd);
        if (subwalk.empty()) continue;

        bool rev = false;
        const int allele = match_graph_allele_exact(subwalk, site.allele_walks, &rev);
        if (allele < 0) continue;

        ctx->results->push_back(GraphReadAllele{
            graph_site_key_str(site, si),
            site.ref_contig.empty() ? site.chrom : site.ref_contig,
            site.pos,
            read_name,
            allele,
            graph_walk_to_string(subwalk),
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

    // Build boundary-node index for catalog matching.
    std::unordered_map<std::string, std::vector<size_t>> node_to_sites;
    for (size_t i = 0; i < catalog.sites.size(); ++i) {
        const GraphSite& site = catalog.sites[i];
        if (!graph_site_is_queryable(site)) continue;
        const GraphWalk& w = site.allele_walks[0];
        if (w.size() < 2) continue;
        node_to_sites[w.front().node].push_back(i);
        if (w.back().node != w.front().node)
            node_to_sites[w.back().node].push_back(i);
    }

    std::vector<GraphReadAllele> results;
    StructuredQueryContext ctx{&node_to_sites, &catalog, min_mapq, &results};

    char* err = nullptr;
    const int rc = pgphase_gbz_query_interval_structured(
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

} // namespace pgphase_collect
