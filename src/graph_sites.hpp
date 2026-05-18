#ifndef PGPHASE_GRAPH_SITES_HPP
#define PGPHASE_GRAPH_SITES_HPP

#include "collect_types.hpp"

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

struct GraphWalkStep {
    std::string node;
    bool reverse = false;

    bool operator==(const GraphWalkStep& other) const {
        return node == other.node && reverse == other.reverse;
    }
};

using GraphWalk = std::vector<GraphWalkStep>;

constexpr int kGraphAlleleMissing = -1;
constexpr int kGraphAlleleAmbiguous = -2;

struct GraphSite {
    std::string chrom;
    std::string ref_contig;
    hts_pos_t pos = 0; // VCF 1-based POS.
    std::string id;
    std::string ref;
    std::vector<std::string> alts;
    std::unordered_map<std::string, std::string> info;
    std::vector<std::string> allele_traversals;
    std::vector<GraphWalk> allele_walks;
    int level = -1;
    std::string parent;
    std::string root;
    hts_pos_t ref_beg = 0;
    hts_pos_t ref_end = 0;
    std::vector<int> conditional_parent_alleles;
    bool has_spanning_deletion = false;
    bool eligible = true;
    std::string skip_reason;

    hts_pos_t order_pos() const { return pos; }
};

struct GraphSiteCatalog {
    std::vector<GraphSite> sites;
};

// When `filters` is non-empty only sites overlapping at least one filter are loaded.
// For bgzipped VCFs with a .tbi/.csi index, tabix is used to seek directly to each
// region; otherwise the file is streamed and non-overlapping lines are skipped.
GraphSiteCatalog load_graph_site_catalog_from_vcf(
    const std::string& path,
    const std::vector<RegionFilter>& filters = {},
    bool keep_allele_traversal_strings = true);

/**
 * @brief True if `path` has a tabix index (.tbi/.csi) loadable by htslib.
 */
bool graph_site_vcf_has_tabix_index(const std::string& path);

/** @throws std::runtime_error if \p path has no tabix index (.tbi/.csi). */
void require_graph_site_vcf_tabix_index(const std::string& path);

/**
 * @brief Reference length (bp) from VCF ##contig meta for a logical sequence name.
 *
 * Matches CHROM to ##contig ID with the same pangenome suffix rules as site loading.
 * @throws std::runtime_error if no matching ##contig record carries length=
 */
hts_pos_t graph_site_contig_end_bp_from_vcf_header(const std::string& vcf_path,
                                                   const std::string& logical_chrom);

/**
 * @brief Incremental graph-site loader: one open VCF + tabix index, many interval queries.
 *
 * Construct one instance per worker thread. Uses half-open \p beg0..\p end0 in the same
 * coordinate convention as `phase-graph` chunk tiling (0-based, end exclusive).
 */
class GraphSitesTabixReader {
public:
    explicit GraphSitesTabixReader(const std::string& vcf_path);
    ~GraphSitesTabixReader();

    GraphSitesTabixReader(const GraphSitesTabixReader&) = delete;
    GraphSitesTabixReader& operator=(const GraphSitesTabixReader&) = delete;

    GraphSitesTabixReader(GraphSitesTabixReader&&) noexcept;
    GraphSitesTabixReader& operator=(GraphSitesTabixReader&&) noexcept;

    GraphSiteCatalog load_half_open_interval(const std::string& logical_chrom,
                                             hts_pos_t beg0,
                                             hts_pos_t end0,
                                             bool keep_allele_traversal_strings);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Returns the ordered list of distinct reference-contig names present in a VCF.
// For bgzipped+indexed VCFs, reads the tabix/csi sequence name table (fast, no
// data scan).  For plain-text VCFs, streams all data lines (linear time).
// Names are normalised: pangenome prefixes like "CHM13#0#chr20" → "chr20".
// Deduplication preserves VCF-header order for indexed files and uses lexicographic
// order for the streaming fallback.
std::vector<std::string> load_graph_site_contig_names(const std::string& vcf_path);

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

} // namespace pgphase_collect

#endif
