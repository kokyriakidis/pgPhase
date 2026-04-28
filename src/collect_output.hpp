#ifndef PGPHASE_COLLECT_OUTPUT_HPP
#define PGPHASE_COLLECT_OUTPUT_HPP

/**
 * @file collect_output.hpp
 * @brief TSV/VCF/read-support writer declarations for collect-bam-variation.
 */

#include "collect_types.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace pgphase_collect {

/**
 * @brief Declarations for TSV/VCF/read-support serialization of collect-bam-variation results.
 *
 * @details Contract: TSV/read-support outputs remain **candidate-space** diagnostics. VCF outputs use a
 * longcallD-like final-call projection (germline/noisy called categories that pass depth/alt-depth gates)
 * rather than dumping all candidate rows. TSV includes CATEGORY (final, post-containment) and INIT_CAT
 * (longcallD first `classify_var_cate` only, for parity checks). Candidate TSV includes `PHASE_SET` /
 * `HAP_ALT` / `HAP_REF` from k-means. Optional `--phase-read-tsv` lists per-read scaffold fields
 * (`HAP`, `PHASE_SET`) for debugging phasing.
 */

/**
 * @brief Maps `VariantType` to TSV/VCF tokens (`SNP`, `INS`, `DEL`, or `UNKNOWN`).
 */
std::string type_name(VariantType type);

/**
 * @brief Maps `VariantCategory` to longcallD-style labels (e.g. `LOW_COV`, `CLEAN_HET_SNP`).
 */
std::string category_name(VariantCategory category);

/**
 * @brief Writes the column header for `--read-support` TSV.
 * @param out Stream receiving the single header line.
 */
void write_read_support_header(std::ostream& out);

/**
 * @brief Appends read×candidate observation lines (CHROM from \a header).
 * @param out Output stream.
 * @param header BAM header for `target_name[tid]`.
 * @param rows Batch of support rows for one chunk or merge step.
 */
void write_read_support_rows(std::ostream& out,
                             const bam_hdr_t* header,
                             const std::vector<ReadSupportRow>& rows);

/**
 * @brief Header line for `--phase-read-tsv` (per-read fields after `assign_hap_based_on_germline_het_vars_kmeans`).
 */
void write_phase_read_tsv_header(std::ostream& out);

/**
 * @brief Appends one line per read in \a chunk with haplotype scaffold columns.
 * @param out Output stream.
 * @param header BAM header for sequence names.
 * @param chunk Chunk after `collect_var_main` (phasing may be a no-op if there are no scaffold variants).
 */
void write_phase_read_tsv_rows(std::ostream& out, const bam_hdr_t* header, const BamChunk& chunk);

/**
 * @brief Writes the main candidate-table TSV header (one line).
 * @param out Output stream.
 */
void write_variants_tsv_header(std::ostream& out);

/**
 * @brief Writes one TSV row per `CandidateVariant` (REF/ALT from \a ref where applicable).
 * @param out Output stream.
 * @param header BAM header for contig names.
 * @param ref Reference sequence cache.
 * @param variants Rows to emit.
 */
void write_variants_tsv_records(std::ostream& out,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants);

/**
 * @brief Writes VCF v4.2 meta lines, FILTER/INFO, ##contig, and `#CHROM` header row.
 * @param out Output stream.
 * @param opts Options (reserved for future SOURCE fields; may be unused).
 * @param header BAM header for contig IDs and lengths.
 */
void write_variants_vcf_header(std::ostream& out, const Options& opts, const bam_hdr_t* header);

/**
 * @brief Writes VCF data lines for all candidates (left-normalized REF/ALT, FILTER, INFO).
 * @param out Output stream.
 * @param opts Thresholds such as `min_sv_len` for SVTYPE/SVLEN.
 * @param header BAM header for CHROM names.
 * @param ref Reference bases for VCF alleles.
 * @param variants Candidates to serialize.
 */
void write_variants_vcf_records(std::ostream& out,
                                const Options& opts,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants);

/**
 * @brief Writes VCF v4.2 header with FORMAT/GT:PS columns for phased candidate output.
 */
void write_phased_variants_vcf_header(std::ostream& out, const Options& opts, const bam_hdr_t* header);

/**
 * @brief Writes phased VCF records with sample `GT:PS` from candidate `hap_*` / `phase_set`.
 */
void write_phased_variants_vcf_records(std::ostream& out,
                                       const Options& opts,
                                       const bam_hdr_t* header,
                                       ReferenceCache& ref,
                                       const CandidateTable& variants);

/**
 * @brief Writes a standalone read-support file from precomputed per-chunk batches.
 * @param opts Must set `read_support_tsv` and BAM/ref paths.
 * @param by_chunk Observations grouped by chunk order.
 * @throws std::runtime_error On BAM/header or file open failure.
 */
void write_read_support_tsv(const Options& opts, const std::vector<std::vector<ReadSupportRow>>& by_chunk);

/**
 * @brief Writes `opts.output_tsv` (header + all rows) in one shot.
 * @param opts I/O paths and primary BAM.
 * @param fai FASTA index for `ReferenceCache`.
 * @param variants Full candidate table.
 * @throws std::runtime_error On BAM/header or output open failure.
 */
void write_variants(const Options& opts, faidx_t* fai, const CandidateTable& variants);

/**
 * @brief Writes `opts.output_vcf` if the path is non-empty; no-op otherwise.
 * @param opts Must include non-empty `output_vcf` when VCF is desired.
 * @param fai FASTA index.
 * @param variants Candidates to emit.
 * @throws std::runtime_error On BAM/header or VCF open failure when output is requested.
 */
void write_variants_vcf(const Options& opts, faidx_t* fai, const CandidateTable& variants);

} // namespace pgphase_collect

#endif
