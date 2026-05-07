#ifndef PGPHASE_GRAPH_BAM_ADAPTER_HPP
#define PGPHASE_GRAPH_BAM_ADAPTER_HPP

#include "collect_types.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace pgphase_collect {

struct FilteredGraphSite {
    std::string site_id;
    int ref_cov = 0;
    int alt_cov = 0;
    int total_cov = 0;
    double allele_fraction = 0.0;
    std::string filter_reason; // "ref_only", "low_depth", "high_af", "low_af"
};

struct GraphBamChunkBuildResult {
    BamChunk chunk;
    std::vector<std::string> site_ids;
    // site_allele_orig_idx[i][j] = original allele-walk index in GraphSite for surviving allele j
    // (0 = ref walk, 1 = first alt walk, etc.). Identity mapping before Phase-1 alt filtering.
    std::vector<std::vector<int>> site_allele_orig_idx;
    // Sites dropped by depth/AF thresholds in Phase 2 (not included in phasing).
    std::vector<FilteredGraphSite> filtered_sites;
};

GraphBamChunkBuildResult build_graph_bam_chunk(const GraphSiteCatalog& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& contig,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts);

void populate_graph_chunk_overlaps(std::vector<GraphBamChunkBuildResult>& graph_chunks);

void phase_graph_bam_chunks(std::vector<GraphBamChunkBuildResult>& graph_chunks,
                            const Options& opts);

void write_graph_bam_site_counts_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks);

void write_graph_bam_read_profiles_tsv(std::ostream& out,
                                       const std::vector<GraphBamChunkBuildResult>& graph_chunks);

void write_graph_bam_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks);

// Writes the phased graph sites in the standard variant TSV format
// (CHROM POS TYPE REF ALT DP ... PHASE_SET HAP_ALT HAP_REF).
// Uses VCF REF/ALT sequences from the catalog — no reference FASTA required.
void write_graph_bam_variants_tsv(std::ostream& out,
                                  const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                                  const GraphSiteCatalog& catalog);

void write_graph_bam_phase_reads_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                                     const std::vector<std::string>& all_read_names = {});

// Writes an unaligned BAM (name-sorted) with HP and PS aux tags per read.
// Reads with no phased assignment get no HP/PS tag. Mirrors the HP/PS tagging
// of the BAM phasing path but without requiring alignment coordinates.
void write_graph_bam_phase_bam(const std::string& out_path,
                               const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                               const std::vector<std::string>& all_read_names = {});

void write_graph_bam_filtered_sites_tsv(std::ostream& out,
                                        const std::vector<GraphBamChunkBuildResult>& graph_chunks);

} // namespace pgphase_collect

#endif
