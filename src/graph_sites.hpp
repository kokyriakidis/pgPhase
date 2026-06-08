#ifndef PGPHASE_GRAPH_SITES_HPP
#define PGPHASE_GRAPH_SITES_HPP

// Snarl site catalog: parsing VCF sites produced by `vg deconstruct -a`,
// validating allele traversals, and providing the GraphSiteCatalog used by
// graph_query for allele matching and graph_bam_adapter for phasing.

#include "collect_types.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace pgphase_collect {

// One step in an allele traversal walk: a graph node ID and its orientation.
struct GraphWalkStep {
    std::string node;
    bool reverse = false;

    bool operator==(const GraphWalkStep& other) const {
        return node == other.node && reverse == other.reverse;
    }
};

using GraphWalk = std::vector<GraphWalkStep>;

constexpr int kGraphAlleleMissing = -1;    // read doesn't span the site
constexpr int kGraphAlleleAmbiguous = -2;  // read matches multiple alleles

// A snarl variation site parsed from a VCF record.
struct GraphSite {
    std::string chrom;          // VCF CHROM (pangenome path name)
    std::string ref_contig;     // suffix after last '#' (e.g. "chr20")
    hts_pos_t pos = 0;          // VCF 1-based POS
    std::string id;              // VCF ID field (snarl identifier)
    std::string ref;
    std::vector<std::string> alts;
    std::vector<std::string> allele_traversals;  // raw AT field entries
    std::vector<GraphWalk> allele_walks;          // parsed walks (one per allele)
    int level = -1;              // snarl nesting level (LV info field)
    std::string parent;          // parent snarl site key (PS info field)
    std::string root;            // root snarl (unused, reserved)
    hts_pos_t ref_beg = 0;      // linearized ref interval start
    hts_pos_t ref_end = 0;      // linearized ref interval end
    std::vector<int> conditional_parent_alleles;  // PA field: parent alleles that gate this child
    bool has_spanning_deletion = false;  // true if any ALT is '*'
    bool eligible = true;        // false if validation failed
    std::string skip_reason;     // why the site was marked ineligible

    hts_pos_t order_pos() const { return pos; }
};

struct GraphSiteCatalogView;

struct GraphSiteCatalog {
    std::vector<GraphSite> sites;

    // Return a view covering all sites (no copy).
    GraphSiteCatalogView view_all() const;
};

// Lightweight non-owning view into a subset of a GraphSiteCatalog's sites.
// Avoids deep-copying GraphSite objects when building per-chunk catalogs.
// The referenced catalog must outlive the view.
struct GraphSiteCatalogView {
    const std::vector<GraphSite>* source = nullptr;
    std::vector<size_t> indices;

    size_t size() const { return indices.size(); }
    bool empty() const { return indices.empty(); }
    const GraphSite& operator[](size_t i) const { return (*source)[indices[i]]; }
};

inline GraphSiteCatalogView GraphSiteCatalog::view_all() const {
    GraphSiteCatalogView v;
    v.source = &sites;
    v.indices.resize(sites.size());
    for (size_t i = 0; i < sites.size(); ++i) v.indices[i] = i;
    return v;
}

// When `filters` is non-empty only sites overlapping at least one filter are loaded.
// For bgzipped VCFs with a .tbi/.csi index, tabix is used to seek directly to each
// region; otherwise the file is streamed and non-overlapping lines are skipped.
GraphSiteCatalog load_graph_site_catalog_from_vcf(
    const std::string& path,
    const std::vector<RegionFilter>& filters = {},
    bool keep_allele_traversal_strings = true);

// True if `path` has a tabix index (.tbi/.csi) loadable by htslib.
bool graph_site_vcf_has_tabix_index(const std::string& path);

// Throws std::runtime_error if `path` has no tabix index (.tbi/.csi).
void require_graph_site_vcf_tabix_index(const std::string& path);

GraphSiteCatalog load_graph_site_catalog_from_gfa_text(const std::string& gfa_text,
                                                       const std::string& reference_sample,
                                                       const std::string& contig);

void write_graph_site_catalog_tsv(std::ostream& out, const GraphSiteCatalog& catalog);


GraphWalk parse_graph_walk(const std::string& walk);

std::string graph_walk_to_string(const GraphWalk& walk);

GraphWalk reverse_graph_walk(const GraphWalk& walk);

int match_graph_allele_exact(const GraphWalk& read_walk,
                             const std::vector<GraphWalk>& allele_walks,
                             bool* reverse_out = nullptr);

std::string graph_site_between_query(const GraphSite& site);

std::string graph_site_validation_skip_reason(const GraphSite& site);

bool graph_site_is_queryable(const GraphSite& site);

// Canonical site key: VCF ID if present, otherwise "chrom:pos:ref".
inline std::string graph_site_key_str(const GraphSite& site) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
}

} // namespace pgphase_collect

#endif
