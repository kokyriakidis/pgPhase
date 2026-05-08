#include "graph_query.hpp"

#include <array>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

namespace pgphase_collect {

namespace {

std::string shell_quote(const std::string& value) {
    std::string out = "'";
    for (char c : value) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

std::string run_command_capture_stdout(const std::string& cmd) {
    std::array<char, 4096> buffer{};
    std::string output;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) throw std::runtime_error("failed to run command: " + cmd);
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
        output += buffer.data();
    }
    const int status = pclose(pipe);
    if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        throw std::runtime_error("command failed: " + cmd);
    }
    return output;
}

std::string temp_path(const std::string& prefix, const std::string& suffix) {
    static std::atomic<int> counter{0};
    std::ostringstream out;
    out << "/tmp/" << prefix << "_" << getpid() << "_" << counter++ << suffix;
    return out.str();
}

std::string ensure_gaf_db(const GraphQueryConfig& config, bool& remove_when_done) {
    remove_when_done = false;
    if (!config.gaf_db.empty()) return config.gaf_db;
    if (config.gaf_file.empty()) throw std::runtime_error("graph query requires --gaf or --gaf-db");

    const std::string db_path = temp_path("pgphase_gaf_base", ".db");
    std::string cmd = shell_quote(config.gaf2db_bin) +
                      " --ref-free " + shell_quote(config.gbz_db) +
                      " --output " + shell_quote(db_path) +
                      " --overwrite " + shell_quote(config.gaf_file);
    (void)run_command_capture_stdout(cmd);
    remove_when_done = true;
    return db_path;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream in(line);
    while (std::getline(in, cur, '\t')) out.push_back(cur);
    return out;
}

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

bool parse_gaf_core_fields(const std::string& line, GafCoreFields& out) {
    return parse_gaf_core_fields_from_column(std::string_view(line), 0, out);
}

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
    const CompactGraphSite& site) {
    auto left_it = boundary_positions.find(site.left);
    auto right_it = boundary_positions.find(site.right);
    if (left_it != boundary_positions.end() && right_it != boundary_positions.end()) {
        for (size_t left_pos : left_it->second) {
            for (size_t right_pos : right_it->second) {
                if (right_pos < left_pos) continue;
                const int allele = match_compact_span(read_walk, left_pos, right_pos,
                                                      site, false);
                if (allele >= 0) return allele;
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
                if (allele >= 0) return allele;
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

size_t file_size_or_throw(const std::string& path) {
    struct stat st {};
    if (stat(path.c_str(), &st) != 0) {
        throw std::runtime_error("failed to stat GAF file: " + path);
    }
    if (st.st_size < 0) throw std::runtime_error("invalid GAF file size: " + path);
    return static_cast<size_t>(st.st_size);
}

size_t scan_gaf_range_compact(const std::string& gaf_file,
                              const CompactGraphSiteIndex& compact_index,
                              int min_mapq,
                              size_t worker_id,
                              size_t byte_begin,
                              size_t byte_end,
                              const GraphReadAlleleThreadEmitter& emit) {
    std::ifstream in(gaf_file, std::ios::binary);
    if (!in) throw std::runtime_error("failed to open GAF file: " + gaf_file);
    in.seekg(static_cast<std::streamoff>(byte_begin));
    if (byte_begin > 0) {
        std::string discard;
        std::getline(in, discard);
    }

    size_t emitted = 0;
    uint32_t mark_stamp = 1;
    std::vector<uint32_t> site_marks(compact_index.sites.size(), 0);
    std::vector<size_t> candidates;
    std::vector<CompactHandle> read_walk;
    std::unordered_map<CompactHandle, std::vector<size_t>> boundary_positions;
    std::string line;

    while (true) {
        const std::streampos line_pos = in.tellg();
        if (line_pos < 0) break;
        if (static_cast<size_t>(line_pos) >= byte_end) break;
        if (!std::getline(in, line)) break;
        if (line.empty() || line[0] == '#') continue;

        GafCoreFields fields;
        if (!parse_gaf_core_fields(line, fields)) continue;
        if (fields.mapq < min_mapq) continue;
        if (!parse_gaf_path_compact(fields.walk, read_walk) || read_walk.empty()) continue;

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
        if (candidates.empty()) continue;

        const std::string read_name(fields.read_name);
        for (size_t compact_site_index : candidates) {
            const CompactGraphSite& site = compact_index.sites[compact_site_index];
            const int allele = match_compact_site_on_read(read_walk, boundary_positions, site);
            if (allele < 0) continue;
            emit(worker_id, GraphReadAllele{
                site.site_id,
                site.chrom,
                site.pos,
                read_name,
                allele,
                "",
                fields.mapq
            });
            ++emitted;
        }
    }
    return emitted;
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
        const int allele = match_compact_site_on_read(read_walk, boundary_positions, site);
        if (allele < 0) continue;
        emit(worker_id, GraphReadAllele{
            site.site_id,
            site.chrom,
            site.pos,
            read_name,
            allele,
            "",
            fields.mapq
        });
        ++emitted;
    }
    return emitted;
}

std::vector<GraphReadAllele> parse_gaf_read_alleles(const std::string& path,
                                                    const GraphSite& site,
                                                    int min_mapq) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open graph query GAF output: " + path);

    std::vector<GraphReadAllele> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> fields = split_tab(line);
        if (fields.size() < 12) continue;
        const int mapq = std::stoi(fields[11]);
        if (mapq < min_mapq) continue;
        const std::string& qname = fields[0];
        const std::string& walk_text = fields[5];
        GraphWalk walk;
        try {
            walk = parse_graph_walk(walk_text);
        } catch (const std::exception&) {
            continue;
        }
        const int allele = match_graph_allele_exact(walk, site.allele_walks);
        if (allele < 0) continue;
        rows.push_back(GraphReadAllele{
            site.id,
            site.ref_contig.empty() ? site.chrom : site.ref_contig,
            site.pos,
            qname,
            allele,
            walk_text,
            mapq
        });
    }
    return rows;
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

size_t scan_gaf_for_catalog_emit_string(const std::string& gaf_file,
                                        const GraphSiteCatalog& catalog,
                                        int min_mapq,
                                        const GraphReadAlleleEmitter& emit) {
    // Build boundary-node index: node_id → site indices whose first allele walk
    // has that node as either the left or right boundary anchor.
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

    std::ifstream in(gaf_file);
    if (!in) throw std::runtime_error("failed to open GAF file: " + gaf_file);

    size_t emitted = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        GafCoreFields fields;
        if (!parse_gaf_core_fields(line, fields)) continue;
        if (fields.mapq < min_mapq) continue;

        GraphWalk read_walk;
        try {
            read_walk = parse_graph_walk(std::string(fields.walk));
        } catch (const std::exception&) { continue; }
        if (read_walk.empty()) continue;

        // Collect candidate sites — any site whose boundary node appears in this walk.
        std::unordered_set<size_t> candidates;
        for (const GraphWalkStep& step : read_walk) {
            auto it = node_to_sites.find(step.node);
            if (it == node_to_sites.end()) continue;
            for (size_t si : it->second) candidates.insert(si);
        }
        if (candidates.empty()) continue;

        const std::string read_name(fields.read_name);

        for (size_t si : candidates) {
            const GraphSite& site = catalog.sites[si];
            const GraphWalkStep& left_bnd  = site.allele_walks[0].front();
            const GraphWalkStep& right_bnd = site.allele_walks[0].back();

            const GraphWalk subwalk = extract_subwalk_between(read_walk, left_bnd, right_bnd);
            if (subwalk.empty()) continue;

            const int allele = match_graph_allele_exact(subwalk, site.allele_walks);
            if (allele < 0) continue;

            emit(GraphReadAllele{
                graph_site_key_str(site, si),
                site.ref_contig.empty() ? site.chrom : site.ref_contig,
                site.pos,
                read_name,
                allele,
                graph_walk_to_string(subwalk),
                fields.mapq
            });
            ++emitted;
        }
    }
    return emitted;
}

} // namespace

std::vector<GraphReadAllele>
scan_gaf_for_catalog(const std::string& gaf_file,
                     const GraphSiteCatalog& catalog,
                     int min_mapq) {
    std::vector<GraphReadAllele> results;
    scan_gaf_for_catalog_emit(gaf_file, catalog, min_mapq,
                              [&](GraphReadAllele&& row) {
                                  results.push_back(std::move(row));
                              });
    return results;
}

size_t scan_gaf_for_catalog_emit(const std::string& gaf_file,
                                 const GraphSiteCatalog& catalog,
                                 int min_mapq,
                                 const GraphReadAlleleEmitter& emit) {
    return scan_gaf_for_catalog_emit_parallel(
        gaf_file, catalog, min_mapq, 1,
        [&](size_t, GraphReadAllele&& row) { emit(std::move(row)); });
}

size_t scan_gaf_for_catalog_emit_parallel_impl(
    const std::string& gaf_file,
    GraphSiteCatalog& catalog,
    int min_mapq,
    size_t threads,
    bool release_catalog_walks,
    const GraphReadAlleleThreadEmitter& emit) {
    CompactGraphSiteIndex compact_index;
    if (!build_compact_graph_site_index(catalog, compact_index, release_catalog_walks)) {
        size_t emitted = 0;
        emitted = scan_gaf_for_catalog_emit_string(
            gaf_file, catalog, min_mapq,
            [&](GraphReadAllele&& row) { emit(0, std::move(row)); });
        return emitted;
    }

    const size_t file_size = file_size_or_throw(gaf_file);
    if (file_size == 0) return 0;
    const size_t worker_count = std::max<size_t>(1, std::min<size_t>(threads, file_size));
    if (worker_count == 1) {
        return scan_gaf_range_compact(gaf_file, compact_index, min_mapq,
                                      0, 0, file_size, emit);
    }

    std::atomic<size_t> emitted_total{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (size_t worker = 0; worker < worker_count; ++worker) {
        const size_t beg = (file_size * worker) / worker_count;
        const size_t end = (file_size * (worker + 1)) / worker_count;
        workers.emplace_back([&, worker, beg, end]() {
            try {
                emitted_total.fetch_add(
                    scan_gaf_range_compact(gaf_file, compact_index, min_mapq,
                                           worker, beg, end, emit),
                    std::memory_order_relaxed);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);
    return emitted_total.load(std::memory_order_relaxed);
}

size_t scan_gaf_for_catalog_emit_parallel(const std::string& gaf_file,
                                          const GraphSiteCatalog& catalog,
                                          int min_mapq,
                                          size_t threads,
                                          const GraphReadAlleleThreadEmitter& emit) {
    GraphSiteCatalog& mutable_catalog = const_cast<GraphSiteCatalog&>(catalog);
    return scan_gaf_for_catalog_emit_parallel_impl(
        gaf_file, mutable_catalog, min_mapq, threads, false, emit);
}

size_t scan_gaf_for_catalog_emit_parallel_releasing_walks(
    const std::string& gaf_file,
    GraphSiteCatalog& catalog,
    int min_mapq,
    size_t threads,
    const GraphReadAlleleThreadEmitter& emit) {
    return scan_gaf_for_catalog_emit_parallel_impl(
        gaf_file, catalog, min_mapq, threads, true, emit);
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

    const int tid = tbx_name2id(raw_tbx, contig.c_str());
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

std::vector<GraphReadAllele>
query_gbz_interval_gaf(const GraphQueryConfig& config,
                       const std::string& sample,
                       const std::string& contig,
                       hts_pos_t beg,
                       hts_pos_t end,
                       const GraphSiteCatalog& catalog)
{
    if (config.gbz_db.empty()) throw std::runtime_error("--gbz-db is required for graph interval query");
    if (config.gaf_db.empty()) throw std::runtime_error("--gaf-db is required for graph interval query");
    if (contig.empty() || end <= beg)
        throw std::runtime_error("invalid interval for graph query: " + contig +
                                 ":" + std::to_string(beg) + ".." + std::to_string(end));

    const std::string gaf_out = temp_path("pgphase_chunk", ".gaf");
    // RAII: always remove temp file, even on exception.
    struct TempFileCleanup {
        const std::string& path;
        ~TempFileCleanup() { std::remove(path.c_str()); }
    } cleanup{gaf_out};

    std::ostringstream cmd;
    cmd << shell_quote(config.query_bin);
    if (!sample.empty())
        cmd << " --sample " << shell_quote(sample);
    cmd << " --contig "    << shell_quote(contig)
        << " --interval "  << beg << ".." << end
        << " --gaf-base "  << shell_quote(config.gaf_db)
        << " --gaf-output "<< shell_quote(gaf_out)
        << " --alignments overlapping "
        << shell_quote(config.gbz_db)
        << " > /dev/null"; // discard GFA subgraph written to stdout
    (void)run_command_capture_stdout(cmd.str());

    return scan_gaf_for_catalog(gaf_out, catalog, config.min_mapq);
}

std::vector<GraphReadAllele>
collect_graph_read_alleles_for_catalog(const GraphQueryConfig& config,
                                       const GraphSiteCatalog& catalog) {
    if (config.gbz_db.empty()) throw std::runtime_error("graph query requires --gbz-db");
    bool remove_gaf_db = false;
    const std::string gaf_db = ensure_gaf_db(config, remove_gaf_db);

    std::vector<GraphReadAllele> all_rows;
    try {
        for (const GraphSite& site : catalog.sites) {
            const std::string skip_reason = graph_site_validation_skip_reason(site);
            if (!site.eligible || !skip_reason.empty()) {
                std::cerr << "Skipping graph site "
                          << (site.id.empty() ? "." : site.id)
                          << ": " << (skip_reason.empty() ? site.skip_reason : skip_reason)
                          << "\n";
                continue;
            }
            const std::string between = graph_site_between_query(site);
            if (between.empty()) continue;

            const std::string gaf_out = temp_path("pgphase_graph_site", ".gaf");
            try {
                std::string cmd = shell_quote(config.query_bin) +
                                  " --between " + shell_quote(between) +
                                  " --limit " + std::to_string(config.between_limit_nodes) +
                                  " --gaf-base " + shell_quote(gaf_db) +
                                  " --gaf-output " + shell_quote(gaf_out) +
                                  " --alignments clipped " +
                                  shell_quote(config.gbz_db);
                (void)run_command_capture_stdout(cmd);

                std::vector<GraphReadAllele> rows = parse_gaf_read_alleles(gaf_out, site, config.min_mapq);
                all_rows.insert(all_rows.end(), rows.begin(), rows.end());
            } catch (const std::exception& e) {
                std::cerr << "Skipping graph site "
                          << (site.id.empty() ? "." : site.id)
                          << " after query failure: " << e.what() << "\n";
            }
            std::remove(gaf_out.c_str());
        }
    } catch (...) {
        if (remove_gaf_db) std::remove(gaf_db.c_str());
        throw;
    }
    if (remove_gaf_db) std::remove(gaf_db.c_str());
    return all_rows;
}

GraphQueryConfig prepare_graph_query_config(const GraphQueryConfig& config,
                                            std::string* temporary_gaf_db) {
    if (temporary_gaf_db != nullptr) temporary_gaf_db->clear();
    if (!config.gaf_db.empty()) return config;
    if (config.gaf_file.empty()) throw std::runtime_error("graph query requires --gaf or --gaf-db");

    bool remove_when_done = false;
    const std::string db_path = ensure_gaf_db(config, remove_when_done);
    GraphQueryConfig out = config;
    out.gaf_db = db_path;
    out.gaf_file.clear();
    if (temporary_gaf_db != nullptr && remove_when_done) *temporary_gaf_db = db_path;
    return out;
}

std::string query_gbz_interval_gfa(const GraphQueryConfig& config,
                                   const std::string& sample,
                                   const std::string& contig,
                                   hts_pos_t beg,
                                   hts_pos_t end,
                                   bool extend_snarls) {
    if (config.gbz_db.empty()) throw std::runtime_error("graph interval query requires --gbz-db");
    if (sample.empty() || contig.empty() || beg < 0 || end <= beg) {
        throw std::runtime_error("invalid graph interval query coordinates");
    }
    std::ostringstream interval;
    interval << beg << ".." << end;

    std::string cmd = shell_quote(config.query_bin) +
                      " --sample " + shell_quote(sample) +
                      " --contig " + shell_quote(contig) +
                      " --interval " + shell_quote(interval.str()) +
                      " --format gfa";
    if (extend_snarls) cmd += " --extend-snarls";
    cmd += " " + shell_quote(config.gbz_db);
    return run_command_capture_stdout(cmd);
}

void write_graph_read_alleles_tsv_header(std::ostream& out) {
    out << "SITE_ID\tCHROM\tPOS\tREAD\tMAPQ\tALLELE\tWALK\n";
}

void write_graph_read_alleles_tsv_rows(std::ostream& out,
                                       const std::vector<GraphReadAllele>& rows) {
    for (const GraphReadAllele& row : rows) {
        out << row.site_id << '\t'
            << row.chrom << '\t'
            << row.pos << '\t'
            << row.read_name << '\t'
            << row.mapq << '\t'
            << row.allele << '\t'
            << row.walk << '\n';
    }
}

void write_graph_read_alleles_tsv(std::ostream& out,
                                  const std::vector<GraphReadAllele>& rows) {
    write_graph_read_alleles_tsv_header(out);
    write_graph_read_alleles_tsv_rows(out, rows);
}

} // namespace pgphase_collect
