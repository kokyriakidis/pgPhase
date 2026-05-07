#ifndef PGPHASE_GRAPH_QUERY_HPP
#define PGPHASE_GRAPH_QUERY_HPP

#include "collect_types.hpp"
#include "graph_sites.hpp"

#include <functional>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

struct GraphReadAllele {
    std::string site_id;
    std::string chrom;
    hts_pos_t pos = 0;
    std::string read_name;
    int allele = kGraphAlleleMissing;
    std::string walk;
    int mapq = 255;
};

struct GraphQueryConfig {
    std::string gbz_db;
    std::string gaf_file;
    std::string gaf_db;
    std::string query_bin = "query";
    std::string gaf2db_bin = "gaf2db";
    int between_limit_nodes = 100000;
    int min_mapq = kDefaultMinMapq;
};

// Fast single-pass scan of a raw GAF file against a pre-loaded site catalog.
// Builds a boundary-node index from the catalog, streams the GAF, and for each
// read extracts and matches the sub-walk spanning each candidate snarl.
// No subprocess calls, no gbz-base database required.
std::vector<GraphReadAllele>
scan_gaf_for_catalog(const std::string& gaf_file,
                     const GraphSiteCatalog& catalog,
                     int min_mapq);

using GraphReadAlleleEmitter = std::function<void(GraphReadAllele&&)>;
using GraphReadAlleleThreadEmitter = std::function<void(size_t, GraphReadAllele&&)>;

size_t scan_gaf_for_catalog_emit(const std::string& gaf_file,
                                 const GraphSiteCatalog& catalog,
                                 int min_mapq,
                                 const GraphReadAlleleEmitter& emit);

// Parallel variant for raw GAF scans. The emitter receives the worker id so callers can
// write thread-local shard files without locking on every observation.
size_t scan_gaf_for_catalog_emit_parallel(const std::string& gaf_file,
                                          const GraphSiteCatalog& catalog,
                                          int min_mapq,
                                          size_t threads,
                                          const GraphReadAlleleThreadEmitter& emit);

// Same as above, but the fast numeric path may release per-site string walk storage
// from the catalog after it has converted those walks to compact handles. This is
// intended for collect-graph-variation's raw-GAF path, where downstream chunk
// phasing only needs allele counts/VCF alleles after observations are gathered.
size_t scan_gaf_for_catalog_emit_parallel_releasing_walks(
    const std::string& gaf_file,
    GraphSiteCatalog& catalog,
    int min_mapq,
    size_t threads,
    const GraphReadAlleleThreadEmitter& emit);

// Per-chunk indexed GAF query using the query binary and a prebuilt gaf.db.
// Calls: query [--sample S] --contig C --interval BEG..END
//              --gaf-base gaf_db --gaf-output TMPFILE --alignments overlapping GBZ
// then scans the temporary GAF for observations of sites in the catalog.
// Peak memory per call = one chunk's worth of reads (mirrors BAM per-chunk access).
std::vector<GraphReadAllele>
query_gbz_interval_gaf(const GraphQueryConfig& config,
                       const std::string& sample,
                       const std::string& contig,
                       hts_pos_t beg,
                       hts_pos_t end,
                       const GraphSiteCatalog& catalog);

// Legacy subprocess-based collection (used for GBZ-path and gaf-db inputs).
std::vector<GraphReadAllele>
collect_graph_read_alleles_for_catalog(const GraphQueryConfig& config,
                                       const GraphSiteCatalog& catalog);

GraphQueryConfig prepare_graph_query_config(const GraphQueryConfig& config,
                                            std::string* temporary_gaf_db);

std::string query_gbz_interval_gfa(const GraphQueryConfig& config,
                                   const std::string& sample,
                                   const std::string& contig,
                                   hts_pos_t beg,
                                   hts_pos_t end,
                                   bool extend_snarls);

void write_graph_read_alleles_tsv(std::ostream& out,
                                  const std::vector<GraphReadAllele>& rows);

} // namespace pgphase_collect

#endif
