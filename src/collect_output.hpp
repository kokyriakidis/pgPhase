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
 * Final-call projection (germline/noisy called categories that pass depth/alt-depth gates)
 * rather than dumping all candidate rows. TSV includes CATEGORY (final, post-containment) and INIT_CAT
 * (first classify pass only, for parity checks). Candidate TSV includes `PHASE_SET` /
 * `HAP_ALT` / `HAP_REF` from k-means. Optional `--phase-read-tsv` lists per-read scaffold fields
 * (`HAP`, `PHASE_SET`) for debugging phasing.
 */

/**
 * @brief Maps `VariantType` to TSV/VCF tokens (`SNP`, `INS`, `DEL`, or `UNKNOWN`).
 */
std::string type_name(VariantType type);

/**
 * Maps VariantCategory to short labels (e.g. LOW_COV, CLEAN_HET_SNP).
 */
std::string category_name(VariantCategory category);

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
