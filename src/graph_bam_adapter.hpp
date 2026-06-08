#ifndef PGPHASE_GRAPH_BAM_ADAPTER_HPP
#define PGPHASE_GRAPH_BAM_ADAPTER_HPP

// Bridge between the graph pipeline (snarl sites + allele walk matching) and
// the shared phasing pipeline (k-means clustering, chunk stitching, VCF
// output).  build_graph_chunk converts per-read allele observations from
// graph_query into a PhasingChunk with CandidateVariant + ReadVariantProfile,
// applying parent-snarl gating, multi-allelic→biallelic decomposition, and
// three-phase depth/AF filtering along the way.

#include "collect_types.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <htslib/sam.h>

#include <iosfwd>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pgphase_collect {

// Per-read phasing result accumulated across chunks during streaming output.
struct PhaseReadOutputRow {
    std::string read_name;
    int chunk_id = -1;
    int hap = 0;                // 1 or 2 when phased, 0 when unphased
    hts_pos_t phase_set = -1;   // genomic position anchoring the phase block
    bool has_phased_assignment = false;
    int copies = 0;             // number of chunks that observed this read
    std::unordered_map<std::string, int> allele_by_site;  // site_id → allele
};

// Record for a site that was dropped during depth/AF filtering (Phase 1-2).
struct FilteredGraphSite {
    std::string site_id;
    int ref_cov = 0;
    int alt_cov = 0;
    int total_cov = 0;
    double allele_fraction = 0.0;
    std::string filter_reason;  // e.g. "ref_only", "low_depth", "high_af"
};

// VCF-level metadata for a snarl site, stored per candidate so the output
// converter doesn't need a global site lookup table.
struct GraphSiteMeta {
    std::string chrom;       // output contig (ref_contig or chrom)
    hts_pos_t pos = 0;
    std::string ref;
    std::vector<std::string> alts;
};

// Output of build_graph_chunk: a PhasingChunk ready for k-means phasing,
// plus graph-specific bookkeeping for VCF output and diagnostics.
struct GraphChunkBuildResult {
    PhasingChunk chunk;
    // Snarl site ID per candidate (parallel to chunk.candidates).
    std::vector<std::string> site_ids;
    // VCF-level metadata per candidate (parallel to chunk.candidates).
    std::vector<GraphSiteMeta> site_meta;
    // Maps final biallelic allele indices back to original allele-walk indices
    // in GraphSite (0 = ref walk, 1 = first alt walk, etc.).
    std::vector<std::vector<int>> site_allele_orig_idx;
    // Sites that survived Phase 1 but were dropped in Phase 2.
    std::vector<FilteredGraphSite> filtered_sites;
};

// Convert graph-space allele observations into a PhasingChunk for phasing.
// Applies parent-snarl gating, multi-allelic→biallelic decomposition,
// and three-phase depth/AF filtering.
GraphChunkBuildResult build_graph_chunk(const GraphSiteCatalogView& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& contig,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts);

// Find reads shared between adjacent chunks (by name merge-intersect)
// and record them in PhasingChunk overlap bookkeeping for stitching.
void populate_graph_chunk_overlaps(std::vector<GraphChunkBuildResult>& graph_chunks);

// K-means hap assignment for every chunk, then overlap detection and stitching.
void phase_graph_chunks(std::vector<GraphChunkBuildResult>& graph_chunks,
                            const Options& opts);

// Fold one stitched chunk's per-read hap/PS assignments into a running map.
void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphChunkBuildResult& gc);

// Write phased-BAM records for reads whose chunks are fully stitched.
// Reads still needed by the next chunk (next_chunk_qnames) are held back;
// pass nullptr to flush all remaining reads (end of contig).
void flush_graph_phase_bam_after_merge(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::unordered_set<std::string>* next_chunk_qnames,
    std::unordered_set<std::string>& emitted_read_names);

// Phase-sites TSV output (used by tests).
void write_graph_phase_sites_tsv_header(std::ostream& out);
void write_graph_phase_sites_tsv_rows(std::ostream& out,
                                          const GraphChunkBuildResult& gc);
void write_graph_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphChunkBuildResult>& graph_chunks);

} // namespace pgphase_collect

#endif
