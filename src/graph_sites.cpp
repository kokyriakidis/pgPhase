#include "graph_sites.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

namespace pgphase_collect {

namespace {

std::vector<std::string> split_char(const std::string& value, char delim) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream in(value);
    while (std::getline(in, cur, delim)) out.push_back(cur);
    if (!value.empty() && value.back() == delim) out.emplace_back();
    return out;
}

std::unordered_map<std::string, std::string> parse_info(const std::string& field) {
    std::unordered_map<std::string, std::string> out;
    if (field.empty() || field == ".") return out;
    for (const std::string& item : split_char(field, ';')) {
        if (item.empty()) continue;
        const size_t eq = item.find('=');
        if (eq == std::string::npos) {
            out.emplace(item, "");
        } else {
            out.emplace(item.substr(0, eq), item.substr(eq + 1));
        }
    }
    return out;
}

int parse_optional_int(const std::unordered_map<std::string, std::string>& info,
                       const char* key) {
    auto it = info.find(key);
    if (it == info.end() || it->second.empty()) return -1;
    return std::stoi(it->second);
}

std::string find_first_info_value(const std::unordered_map<std::string, std::string>& info,
                                  const std::vector<const char*>& keys) {
    for (const char* key : keys) {
        auto it = info.find(key);
        if (it != info.end()) return it->second;
    }
    return "";
}

// Extract the AT (allele traversal) field from VCF INFO.  Falls back to UT
// and LVAT for compatibility with older vg deconstruct output.
std::vector<std::string> parse_allele_traversals(const std::unordered_map<std::string, std::string>& info) {
    const std::string value = find_first_info_value(info, {"AT", "UT", "LVAT"});
    if (value.empty()) return {};
    return split_char(value, ',');
}

std::vector<int> parse_optional_int_list(const std::unordered_map<std::string, std::string>& info,
                                         const std::vector<const char*>& keys) {
    const std::string value = find_first_info_value(info, keys);
    std::vector<int> out;
    if (value.empty() || value == ".") return out;
    for (const std::string& field : split_char(value, ',')) {
        if (field.empty() || field == ".") continue;
        out.push_back(std::stoi(field));
    }
    return out;
}

// Parse each AT entry into a GraphWalk.  Spanning deletions ('*') and
// missing entries ('.') produce empty walks.  Sets `malformed` if any
// entry fails to parse (the site will be marked ineligible).
std::vector<GraphWalk> parse_allele_walks(const std::vector<std::string>& traversals,
                                          bool& malformed) {
    std::vector<GraphWalk> out;
    out.reserve(traversals.size());
    for (const std::string& traversal : traversals) {
        if (traversal == "*" || traversal == ".") {
            out.emplace_back();
            continue;
        }
        try {
            out.push_back(parse_graph_walk(traversal));
        } catch (const std::exception&) {
            malformed = true;
            out.emplace_back();
        }
    }
    return out;
}

bool contains_spanning_deletion(const std::vector<std::string>& alts) {
    return std::find(alts.begin(), alts.end(), "*") != alts.end();
}

// Extract the suffix after the last '#' for pangenome path matching.
// "CHM13#0#chr20" → "chr20", "chr20" → "chr20".
std::string chrom_suffix_match_key(const std::string& s) {
    const size_t h = s.rfind('#');
    return h == std::string::npos ? s : s.substr(h + 1);
}

bool chrom_names_match_vcf(const std::string& a, const std::string& b) {
    return a == b || chrom_suffix_match_key(a) == chrom_suffix_match_key(b);
}

// Sort sites by (chrom, pos, id) and run structural validation on each.
// Parse-time skip reasons (e.g. malformed_allele_traversal) are preserved
// because the validator can't detect them from the stored empty walk.
void finalize_graph_site_catalog_inplace(GraphSiteCatalog& catalog) {
    std::stable_sort(catalog.sites.begin(), catalog.sites.end(),
                     [](const GraphSite& lhs, const GraphSite& rhs) {
                         if (lhs.chrom != rhs.chrom) return lhs.chrom < rhs.chrom;
                         if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
                         return lhs.id < rhs.id;
                     });
    for (GraphSite& site : catalog.sites) {
        // Preserve parse-time skip reasons (e.g. malformed_allele_traversal)
        // that the structural validator cannot detect because the malformed
        // walk was stored as an empty entry.
        if (!site.skip_reason.empty()) continue;
        site.skip_reason = graph_site_validation_skip_reason(site);
        site.eligible = site.skip_reason.empty();
    }
}

// Parse one VCF data line into a GraphSite and append it to the catalog.
// Extracts CHROM, POS, ID, REF, ALT, and INFO fields (AT, LV, PS, PA, END).
void append_graph_site_from_vcf_data_line(const char* buf,
                                          size_t len,
                                          GraphSiteCatalog& catalog,
                                          bool keep_allele_traversal_strings) {
    if (len == 0 || buf[0] == '#') return;
    const std::string line(buf, len);
    std::vector<std::string> fields = split_char(line, '\t');
    if (fields.size() < 8) return;

    GraphSite site;
    site.chrom = fields[0];
    const auto hash_pos = fields[0].rfind('#');
    site.ref_contig = (hash_pos != std::string::npos)
                      ? fields[0].substr(hash_pos + 1)
                      : fields[0];
    site.pos = static_cast<hts_pos_t>(std::stoll(fields[1]));
    site.id = fields[2];
    site.ref = fields[3];
    if (fields[4] != ".") site.alts = split_char(fields[4], ',');
    const auto info = parse_info(fields[7]);
    site.allele_traversals = parse_allele_traversals(info);
    bool malformed_walk = false;
    site.allele_walks = parse_allele_walks(site.allele_traversals, malformed_walk);
    site.level = parse_optional_int(info, "LV");
    site.parent = find_first_info_value(info, {"PS", "PARENT", "Parent"});
    site.root = find_first_info_value(info, {"RS", "ROOT", "Root"});
    const std::string rc_value = find_first_info_value(info, {"RC", "REF_CONTIG", "RefContig"});
    if (!rc_value.empty()) {
        const auto hp = rc_value.rfind('#');
        site.ref_contig = (hp != std::string::npos) ? rc_value.substr(hp + 1) : rc_value;
    }
    site.ref_beg = site.pos;
    site.ref_end = site.pos;
    const std::string end_value = find_first_info_value(info, {"END"});
    if (!end_value.empty()) site.ref_end = static_cast<hts_pos_t>(std::stoll(end_value));
    site.conditional_parent_alleles =
        parse_optional_int_list(info, {"PA", "PARENT_ALLELE", "PARENT_ALLELES"});
    site.has_spanning_deletion = contains_spanning_deletion(site.alts);
    if (malformed_walk) {
        site.eligible = false;
        site.skip_reason = "malformed_allele_traversal";
    } else {
        site.skip_reason = graph_site_validation_skip_reason(site);
        site.eligible = site.skip_reason.empty();
    }
    if (!keep_allele_traversal_strings) {
        std::vector<std::string>().swap(site.allele_traversals);
    }
    catalog.sites.push_back(std::move(site));
}

// Load sites from a tabix-indexed VCF, restricted to the given regions.
// Handles pangenome path name mismatches (e.g. filter asks for "chr20" but
// VCF uses "CHM13#0#chr20") via suffix matching.
void append_graph_sites_tabix_filtered(htsFile* fp,
                                       tbx_t* tbx,
                                       const std::vector<RegionFilter>& filters,
                                       GraphSiteCatalog& catalog,
                                       bool keep_allele_traversal_strings) {
    struct HtsItrDeleter {
        void operator()(hts_itr_t* p) const {
            if (p) hts_itr_destroy(p);
        }
    };

    int n_vcf_seqs = 0;
    const char** vcf_seqs = tbx_seqnames(tbx, &n_vcf_seqs);
    std::unordered_map<std::string, std::string> filter_to_vcf_chrom;
    for (const RegionFilter& filter : filters) {
        if (!filter.enabled || filter_to_vcf_chrom.count(filter.chrom)) continue;
        if (tbx_name2id(tbx, filter.chrom.c_str()) >= 0) {
            filter_to_vcf_chrom[filter.chrom] = filter.chrom;
        } else {
            for (int si = 0; si < n_vcf_seqs; ++si) {
                if (chrom_names_match_vcf(filter.chrom, std::string(vcf_seqs[si]))) {
                    filter_to_vcf_chrom[filter.chrom] = vcf_seqs[si];
                    break;
                }
            }
        }
    }
    free(vcf_seqs);

    kstring_t str = KS_INITIALIZE;
    struct KsGuard {
        kstring_t* s;
        ~KsGuard() {
            if (s) ks_free(s);
        }
    } ks_guard{&str};

    auto make_rstr = [](const std::string& chrom, hts_pos_t beg, hts_pos_t end) {
        std::string s = chrom;
        if (beg > 0) {
            s += ":" + std::to_string(beg);
            if (end > 0) s += "-" + std::to_string(end);
        }
        return s;
    };

    for (const RegionFilter& filter : filters) {
        if (!filter.enabled) continue;
        auto it = filter_to_vcf_chrom.find(filter.chrom);
        if (it == filter_to_vcf_chrom.end()) continue;
        std::unique_ptr<hts_itr_t, HtsItrDeleter> itr(
            tbx_itr_querys(tbx,
                           make_rstr(it->second, filter.beg, filter.end).c_str()));
        if (!itr) continue;
        while (tbx_itr_next(fp, tbx, itr.get(), &str) >= 0) {
            if (str.l > 0)
                append_graph_site_from_vcf_data_line(str.s, str.l, catalog,
                                                     keep_allele_traversal_strings);
        }
    }
}

struct GfaWalkRecord {
    std::string sample;
    std::string hap;
    std::string contig;
    hts_pos_t start = 0;
    hts_pos_t end = 0;
    GraphWalk walk;
};

std::string step_key(const GraphWalkStep& step) {
    std::string out = step.node;
    out.push_back(step.reverse ? '-' : '+');
    return out;
}

hts_pos_t parse_pos_or_zero(const std::string& value) {
    if (value.empty() || value == "*") return 0;
    return static_cast<hts_pos_t>(std::stoll(value));
}

std::vector<GfaWalkRecord>
parse_gfa_walk_records(const std::string& gfa_text,
                       std::unordered_map<std::string, std::string>& node_seqs) {
    std::vector<GfaWalkRecord> walks;
    std::istringstream in(gfa_text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_char(line, '\t');
        if (fields.empty()) continue;
        if (fields[0] == "S" && fields.size() >= 3) {
            node_seqs[fields[1]] = fields[2];
        } else if (fields[0] == "W" && fields.size() >= 7) {
            try {
                GfaWalkRecord rec;
                rec.sample = fields[1];
                rec.hap = fields[2];
                rec.contig = fields[3];
                rec.start = parse_pos_or_zero(fields[4]);
                rec.end = parse_pos_or_zero(fields[5]);
                rec.walk = parse_graph_walk(fields[6]);
                if (!rec.walk.empty()) walks.push_back(std::move(rec));
            } catch (const std::exception&) {
                continue;
            }
        }
    }
    return walks;
}

size_t find_step_at_or_after(const GraphWalk& walk, const GraphWalkStep& step, size_t beg) {
    for (size_t i = beg; i < walk.size(); ++i) {
        if (walk[i] == step) return i;
    }
    return std::numeric_limits<size_t>::max();
}

GraphWalk slice_walk(const GraphWalk& walk, size_t beg, size_t end_inclusive) {
    if (beg >= walk.size() || end_inclusive >= walk.size() || end_inclusive < beg) return {};
    return GraphWalk(walk.begin() + static_cast<std::ptrdiff_t>(beg),
                     walk.begin() + static_cast<std::ptrdiff_t>(end_inclusive + 1));
}

GraphWalk orient_walk_to_reference(const GraphWalk& walk, const GraphWalk& reference_walk) {
    if (walk.empty() || reference_walk.empty()) return walk;

    std::unordered_set<std::string> ref_steps;
    ref_steps.reserve(reference_walk.size() * 2);
    for (const GraphWalkStep& step : reference_walk) ref_steps.insert(step_key(step));

    const auto count_ref_steps = [&](const GraphWalk& candidate) {
        size_t n = 0;
        for (const GraphWalkStep& step : candidate) {
            if (ref_steps.count(step_key(step))) ++n;
        }
        return n;
    };

    const GraphWalk rev = reverse_graph_walk(walk);
    return count_ref_steps(rev) > count_ref_steps(walk) ? rev : walk;
}

std::vector<hts_pos_t>
reference_step_positions(const GfaWalkRecord& ref,
                         const std::unordered_map<std::string, std::string>& node_seqs) {
    std::vector<hts_pos_t> positions;
    positions.reserve(ref.walk.size());
    hts_pos_t pos = ref.start;
    for (const GraphWalkStep& step : ref.walk) {
        positions.push_back(pos);
        auto it = node_seqs.find(step.node);
        pos += (it == node_seqs.end() || it->second.empty())
                 ? 1
                 : static_cast<hts_pos_t>(it->second.size());
    }
    return positions;
}

} // namespace

GraphSiteCatalog load_graph_site_catalog_from_vcf(
    const std::string& path,
    const std::vector<RegionFilter>& filters,
    bool keep_allele_traversal_strings)
{
    GraphSiteCatalog catalog;

    bool has_filters = false;
    for (const RegionFilter& f : filters) {
        if (f.enabled) { has_filters = true; break; }
    }

    struct TbxDeleter { void operator()(tbx_t* p) const { tbx_destroy(p); } };
    struct HtsFileDeleter { void operator()(htsFile* p) const { hts_close(p); } };

    if (has_filters) {
        // Region-filtered: require tabix index for efficient random access.
        std::unique_ptr<tbx_t, TbxDeleter> tbx(tbx_index_load(path.c_str()));
        if (!tbx)
            throw std::runtime_error(
                "sites VCF must be bgzip-compressed with a tabix index: " + path);
        std::unique_ptr<htsFile, HtsFileDeleter> fp(hts_open(path.c_str(), "r"));
        if (!fp) throw std::runtime_error("failed to open graph site VCF: " + path);
        append_graph_sites_tabix_filtered(fp.get(), tbx.get(), filters, catalog,
                                           keep_allele_traversal_strings);
    } else {
        // No region filters: stream the entire file.
        std::unique_ptr<htsFile, HtsFileDeleter> fp(hts_open(path.c_str(), "r"));
        if (!fp) throw std::runtime_error("failed to open graph site VCF: " + path);
        kstring_t str = KS_INITIALIZE;
        while (hts_getline(fp.get(), '\n', &str) >= 0) {
            if (str.l == 0 || str.s[0] == '#') continue;
            append_graph_site_from_vcf_data_line(str.s, str.l, catalog,
                                                 keep_allele_traversal_strings);
        }
        ks_free(&str);
    }

    finalize_graph_site_catalog_inplace(catalog);
    return catalog;
}

bool graph_site_vcf_has_tabix_index(const std::string& path) {
    tbx_t* tbx = tbx_index_load(path.c_str());
    if (!tbx) return false;
    tbx_destroy(tbx);
    return true;
}

void require_graph_site_vcf_tabix_index(const std::string& path) {
    if (graph_site_vcf_has_tabix_index(path)) return;
    throw std::runtime_error(
        "graph-site VCF must be bgzip-compressed with a tabix index (.tbi or .csi): " +
        path);
}

GraphSiteCatalog load_graph_site_catalog_from_gfa_text(const std::string& gfa_text,
                                                       const std::string& reference_sample,
                                                       const std::string& contig) {
    std::unordered_map<std::string, std::string> node_seqs;
    const std::vector<GfaWalkRecord> walks = parse_gfa_walk_records(gfa_text, node_seqs);

    int ref_i = -1;
    for (size_t i = 0; i < walks.size(); ++i) {
        if (walks[i].contig == contig && walks[i].sample == reference_sample) {
            ref_i = static_cast<int>(i);
            break;
        }
    }
    if (ref_i < 0) {
        for (size_t i = 0; i < walks.size(); ++i) {
            if (walks[i].contig == contig) {
                ref_i = static_cast<int>(i);
                break;
            }
        }
    }
    if (ref_i < 0) {
        throw std::runtime_error("GBZ interval GFA did not contain a W line for contig: " + contig);
    }

    const GfaWalkRecord& ref = walks[static_cast<size_t>(ref_i)];
    const std::vector<hts_pos_t> ref_positions = reference_step_positions(ref, node_seqs);

    GraphSiteCatalog catalog;
    std::unordered_map<std::string, size_t> site_by_key;

    const auto add_event = [&](size_t ref_beg, size_t ref_end,
                               const GraphWalk& ref_sub, const GraphWalk& alt_sub) {
        if (ref_sub.size() < 2 || alt_sub.size() < 2 || ref_sub == alt_sub) return;
        const std::string ref_text = graph_walk_to_string(ref_sub);
        const std::string alt_text = graph_walk_to_string(alt_sub);
        const std::string key = contig + "|" + step_key(ref_sub.front()) + "|" +
                                step_key(ref_sub.back()) + "|" + ref_text;

        size_t site_i = 0;
        auto it = site_by_key.find(key);
        if (it == site_by_key.end()) {
            GraphSite site;
            site.chrom = contig;
            site.ref_contig = contig;
            site.pos = ref_beg < ref_positions.size() ? ref_positions[ref_beg] + 1 : 1;
            site.id = "auto_" + std::to_string(catalog.sites.size() + 1) + "_" +
                      step_key(ref_sub.front()) + "_" + step_key(ref_sub.back());
            site.ref = ".";
            site.level = 0;
            site.ref_beg = site.pos;
            site.ref_end = ref_end < ref_positions.size() ? ref_positions[ref_end] + 1 : site.pos;
            site.allele_traversals.push_back(ref_text);
            site.allele_walks.push_back(ref_sub);
            catalog.sites.push_back(std::move(site));
            site_i = catalog.sites.size() - 1;
            site_by_key.emplace(key, site_i);
        } else {
            site_i = it->second;
        }

        GraphSite& site = catalog.sites[site_i];
        if (std::find(site.allele_traversals.begin(),
                      site.allele_traversals.end(),
                      alt_text) == site.allele_traversals.end()) {
            site.allele_traversals.push_back(alt_text);
            site.allele_walks.push_back(alt_sub);
            site.alts.push_back(".");
        }
        (void)ref_end;
    };

    for (size_t walk_i = 0; walk_i < walks.size(); ++walk_i) {
        if (static_cast<int>(walk_i) == ref_i || walks[walk_i].contig != ref.contig) continue;
        const GraphWalk alt = orient_walk_to_reference(walks[walk_i].walk, ref.walk);
        size_t ref_pos = 0;
        size_t alt_scan = 0;
        while (ref_pos + 1 < ref.walk.size()) {
            const size_t alt_left = find_step_at_or_after(alt, ref.walk[ref_pos], alt_scan);
            if (alt_left == std::numeric_limits<size_t>::max()) {
                ++ref_pos;
                continue;
            }
            if (alt_left + 1 < alt.size() && alt[alt_left + 1] == ref.walk[ref_pos + 1]) {
                alt_scan = alt_left + 1;
                ++ref_pos;
                continue;
            }

            size_t ref_right = std::numeric_limits<size_t>::max();
            size_t alt_right = std::numeric_limits<size_t>::max();
            for (size_t j = ref_pos + 1; j < ref.walk.size(); ++j) {
                const size_t found = find_step_at_or_after(alt, ref.walk[j], alt_left + 1);
                if (found != std::numeric_limits<size_t>::max()) {
                    ref_right = j;
                    alt_right = found;
                    break;
                }
            }
            if (ref_right == std::numeric_limits<size_t>::max()) break;

            add_event(ref_pos, ref_right,
                      slice_walk(ref.walk, ref_pos, ref_right),
                      slice_walk(alt, alt_left, alt_right));
            ref_pos = ref_right;
            alt_scan = alt_right;
        }
    }

    std::stable_sort(catalog.sites.begin(), catalog.sites.end(),
                     [](const GraphSite& lhs, const GraphSite& rhs) {
                         if (lhs.chrom != rhs.chrom) return lhs.chrom < rhs.chrom;
                         if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
                         return lhs.id < rhs.id;
                     });
    return catalog;
}

GraphWalk parse_graph_walk(const std::string& walk) {
    GraphWalk out;
    size_t i = 0;
    while (i < walk.size()) {
        const char c = walk[i];
        if (c != '>' && c != '<') {
            throw std::runtime_error("invalid graph walk step in: " + walk);
        }
        const bool reverse = (c == '<');
        ++i;
        const size_t node_beg = i;
        while (i < walk.size() && walk[i] != '>' && walk[i] != '<') ++i;
        if (node_beg == i) {
            throw std::runtime_error("empty node id in graph walk: " + walk);
        }
        out.push_back(GraphWalkStep{walk.substr(node_beg, i - node_beg), reverse});
    }
    return out;
}

std::string graph_walk_to_string(const GraphWalk& walk) {
    std::string out;
    for (const GraphWalkStep& step : walk) {
        out.push_back(step.reverse ? '<' : '>');
        out += step.node;
    }
    return out;
}

GraphWalk reverse_graph_walk(const GraphWalk& walk) {
    GraphWalk out;
    out.reserve(walk.size());
    for (auto it = walk.rbegin(); it != walk.rend(); ++it) {
        out.push_back(GraphWalkStep{it->node, !it->reverse});
    }
    return out;
}

int match_graph_allele_exact(const GraphWalk& read_walk,
                             const std::vector<GraphWalk>& allele_walks,
                             bool* reverse_out) {
    int match = kGraphAlleleMissing;
    bool match_reverse = false;
    const GraphWalk reverse_walk = reverse_graph_walk(read_walk);
    for (size_t allele_i = 0; allele_i < allele_walks.size(); ++allele_i) {
        if (read_walk == allele_walks[allele_i]) {
            if (match != kGraphAlleleMissing) return kGraphAlleleAmbiguous;
            match = static_cast<int>(allele_i);
            match_reverse = false;
        } else if (reverse_walk == allele_walks[allele_i]) {
            if (match != kGraphAlleleMissing) return kGraphAlleleAmbiguous;
            match = static_cast<int>(allele_i);
            match_reverse = true;
        }
    }
    if (match >= 0 && reverse_out) *reverse_out = match_reverse;
    return match;
}

// Structural validation: checks that allele walks have consistent boundaries,
// at least 2 unique alleles, and walks of length ≥2 (left + right boundary).
std::string graph_site_validation_skip_reason(const GraphSite& site) {
    if (!site.skip_reason.empty() && !site.eligible) return site.skip_reason;
    if (!site.allele_traversals.empty() &&
        site.allele_walks.size() != site.allele_traversals.size()) {
        return "allele_walk_traversal_count_mismatch";
    }
    if (site.allele_walks.size() < 2) return "too_few_alleles";

    const GraphWalk* boundary_walk = nullptr;
    std::unordered_set<std::string> unique_walks;
    for (const GraphWalk& walk : site.allele_walks) {
        if (walk.empty()) continue;
        if (walk.size() < 2) return "allele_walk_without_boundaries";
        unique_walks.insert(graph_walk_to_string(walk));
        if (boundary_walk == nullptr) {
            boundary_walk = &walk;
            continue;
        }
        if (!(walk.front() == boundary_walk->front()) ||
            !(walk.back() == boundary_walk->back())) {
            return "incompatible_allele_boundaries";
        }
    }
    if (boundary_walk == nullptr) return "no_queryable_allele_walks";
    if (unique_walks.size() < 2) return "too_few_unique_alleles";
    return "";
}

bool graph_site_is_queryable(const GraphSite& site) {
    return site.eligible && graph_site_validation_skip_reason(site).empty();
}

// Build a "left+:right+" boundary query string for the GBZ FFI interval query.
std::string graph_site_between_query(const GraphSite& site) {
    if (!graph_site_is_queryable(site)) return "";
    const GraphWalk* query_walk = nullptr;
    for (const GraphWalk& candidate : site.allele_walks) {
        if (!candidate.empty()) {
            query_walk = &candidate;
            break;
        }
    }
    if (query_walk == nullptr) return "";
    const GraphWalk& walk = *query_walk;
    const GraphWalkStep& left = walk.front();
    const GraphWalkStep& right = walk.back();
    std::string out;
    out += left.node;
    out.push_back(left.reverse ? '-' : '+');
    out.push_back(':');
    out += right.node;
    out.push_back(right.reverse ? '-' : '+');
    return out;
}


void write_graph_site_catalog_tsv(std::ostream& out, const GraphSiteCatalog& catalog) {
    out << "CHROM\tREF_CONTIG\tPOS\tEND\tID\tLEVEL\tPARENT\tROOT\tCOND_PARENT_ALLELES\tN_ALTS\tN_TRAVERSALS\tHAS_STAR\tELIGIBLE\tSKIP_REASON\tREF\tALT\n";
    for (const GraphSite& site : catalog.sites) {
        const std::string skip_reason = graph_site_validation_skip_reason(site);
        out << site.chrom << '\t'
            << (site.ref_contig.empty() ? site.chrom : site.ref_contig) << '\t'
            << site.pos << '\t'
            << site.ref_end << '\t'
            << site.id << '\t'
            << site.level << '\t'
            << (site.parent.empty() ? "." : site.parent) << '\t'
            << (site.root.empty() ? "." : site.root) << '\t';
        for (size_t i = 0; i < site.conditional_parent_alleles.size(); ++i) {
            if (i > 0) out << ',';
            out << site.conditional_parent_alleles[i];
        }
        if (site.conditional_parent_alleles.empty()) out << '.';
        out << '\t'
            << site.alts.size() << '\t'
            << site.allele_traversals.size() << '\t'
            << (site.has_spanning_deletion ? 1 : 0) << '\t'
            << (skip_reason.empty() && site.eligible ? 1 : 0) << '\t'
            << (skip_reason.empty() ? "." : skip_reason) << '\t'
            << site.ref << '\t';
        for (size_t i = 0; i < site.alts.size(); ++i) {
            if (i > 0) out << ',';
            out << site.alts[i];
        }
        out << '\n';
    }
}

} // namespace pgphase_collect
