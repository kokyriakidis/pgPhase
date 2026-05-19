#ifndef PGPHASE_GRAPH_BAM_ADAPTER_HPP
#define PGPHASE_GRAPH_BAM_ADAPTER_HPP

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

// Compact per-read phasing result; accumulated across chunks during streaming.
struct PhaseReadOutputRow {
    std::string read_name;
    int chunk_id = -1;
    int hap = 0;
    hts_pos_t phase_set = -1;
    bool has_phased_assignment = false;
    int copies = 0;
    std::unordered_map<std::string, int> allele_by_site;
};

struct FilteredGraphSite {
    std::string site_id;
    int ref_cov = 0;
    int alt_cov = 0;
    int total_cov = 0;
    double allele_fraction = 0.0;
    std::string filter_reason; // "ref_only", "low_depth", "high_af", "low_af"
};

/** CHROM/POS/REF/ALT snapshot per final candidate for streaming TSV (no global site map). */
struct GraphVariantEmitRow {
    std::string chrom;
    hts_pos_t pos = 0;
    std::string snarl_id;
    std::string ref;
    std::string alt;
};

struct GraphBamChunkBuildResult {
    BamChunk chunk;
    std::vector<std::string> site_ids;
    // site_allele_orig_idx[i][j] = original allele-walk index in GraphSite for surviving allele j
    // (0 = ref walk, 1 = first alt walk, etc.). Identity mapping before Phase-1 alt filtering.
    std::vector<std::vector<int>> site_allele_orig_idx;
    // Sites dropped by depth/AF thresholds in Phase 2 (not included in phasing).
    std::vector<FilteredGraphSite> filtered_sites;
    /** Parallel to chunk.candidates; filled by build_graph_bam_chunk for variant TSV emission. */
    std::vector<GraphVariantEmitRow> variant_emit_rows;
    /** Reference/pangenome contig name for phase-read TSV CHROM columns (may be empty). */
    std::string graph_phase_contig;
};

GraphBamChunkBuildResult build_graph_bam_chunk(const GraphSiteCatalog& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& contig,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts,
                                               const ReadWalkMap* read_walks = nullptr);

// ── Overlap / stitch helpers ──────────────────────────────────────────────────

// Populate overlap bookkeeping for all adjacent pairs in the full vector
// (used by the batch-at-once path in collect_pipeline.cpp).
void populate_graph_chunk_overlaps(std::vector<GraphBamChunkBuildResult>& graph_chunks);

// Populate overlap bookkeeping for exactly one adjacent pair (streaming path).
void populate_graph_chunk_pair_overlaps(GraphBamChunkBuildResult& pre,
                                        GraphBamChunkBuildResult& cur);

// Per-chunk k-means hap assignment (no cross-chunk stitching).
void assign_graph_chunk_hap(GraphBamChunkBuildResult& gc, const Options& opts);

// Stitch exactly two adjacent chunks: populate overlaps, then flip cur if needed.
void stitch_graph_chunk_pair(GraphBamChunkBuildResult& pre,
                             GraphBamChunkBuildResult& cur,
                             const Options& opts);

// Full batch: assign hap for every chunk, populate all overlaps, stitch.
void phase_graph_bam_chunks(std::vector<GraphBamChunkBuildResult>& graph_chunks,
                            const Options& opts);

// ── Streaming per-read accumulator ───────────────────────────────────────────

// Merge one finalized chunk's per-read hap/PS assignments into a persistent map.
// Call after the chunk has been fully stitched.
void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphBamChunkBuildResult& gc);

// ── Per-output-file header writers (call once before the streaming loop) ─────

void write_graph_bam_site_counts_tsv_header(std::ostream& out);
void write_graph_bam_read_profiles_tsv_header(std::ostream& out);
void write_graph_bam_filtered_sites_tsv_header(std::ostream& out);
void write_graph_bam_phase_sites_tsv_header(std::ostream& out);
void write_graph_bam_variants_tsv_header(std::ostream& out);

// ── Per-chunk row writers (call in streaming loop, after the matching header) ─

void write_graph_bam_site_counts_tsv_rows(std::ostream& out,
                                          const GraphBamChunkBuildResult& gc);
void write_graph_bam_read_profiles_tsv_rows(std::ostream& out,
                                            const GraphBamChunkBuildResult& gc);
void write_graph_bam_filtered_sites_tsv_rows(std::ostream& out,
                                             const GraphBamChunkBuildResult& gc);
void write_graph_bam_phase_sites_tsv_rows(std::ostream& out,
                                          const GraphBamChunkBuildResult& gc);

// Requires a pre-built site_id → GraphSite* map from the catalog (build once, reuse).
void write_graph_bam_variants_tsv_rows(
    std::ostream& out,
    const GraphBamChunkBuildResult& gc,
    const std::unordered_map<std::string, const GraphSite*>& site_map);

// Same as above with owned GraphSite rows (used when join maps are built externally).
void write_graph_bam_variants_tsv_rows(
    std::ostream& out,
    const GraphBamChunkBuildResult& gc,
    const std::unordered_map<std::string, GraphSite>& site_by_id);

// Uses gc.variant_emit_rows (chunk-local REF/ALT); no catalog map.
void write_graph_bam_variants_tsv_rows(std::ostream& out,
                                       const GraphBamChunkBuildResult& gc);

// ── Phase-BAM (streaming phase-graph path) ──────────────────────────────────

// After merging a stitched chunk into rows_by_read, drop reads that cannot appear
// in next_chunk_qnames (typically qnames in the next genomic chunk). Pass nullptr
// for next_chunk_qnames to flush everything (end of contig / pipeline).
void flush_graph_phase_bam_after_merge(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::unordered_set<std::string>* next_chunk_qnames,
    std::unordered_set<std::string>& emitted_read_names);

// Reads from all_read_names_sorted not present in emitted_read_names (no graph observations).
void write_graph_phase_bam_for_unobserved(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    const std::unordered_set<std::string>& emitted_read_names,
    const std::vector<std::string>& all_read_names_sorted);

// ── End-of-run writers from accumulated rows_by_read map ─────────────────────

// Write an unaligned BAM with HP/PS tags from a pre-built per-read map.
// all_read_names: additional read names (e.g. from GAF) to include as unphased.
void write_graph_bam_phase_bam_from_rows(
    const std::string& out_path,
    const std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::vector<std::string>& all_read_names);

// ── Legacy multi-chunk writers (used by collect_pipeline.cpp) ────────────────

void write_graph_bam_site_counts_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks);

void write_graph_bam_read_profiles_tsv(std::ostream& out,
                                       const std::vector<GraphBamChunkBuildResult>& graph_chunks);

void write_graph_bam_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks);

// Writes the phased graph sites in the standard variant TSV format
// (CHROM POS TYPE REF ALT DP ... PHASE_SET HAP_ALT HAP_REF).
// Each chunk must have variant_emit_rows from build_graph_bam_chunk (embedded REF/ALT).
// catalog is unused (kept for a stable call signature).
void write_graph_bam_variants_tsv(std::ostream& out,
                                  const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                                  const GraphSiteCatalog& catalog);

// Writes an unaligned BAM (name-sorted) with HP and PS aux tags per read.
// Reads with no phased assignment get no HP/PS tag.
void write_graph_bam_phase_bam(const std::string& out_path,
                               const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                               const std::vector<std::string>& all_read_names = {});

void write_graph_bam_filtered_sites_tsv(std::ostream& out,
                                        const std::vector<GraphBamChunkBuildResult>& graph_chunks);

} // namespace pgphase_collect

#endif
