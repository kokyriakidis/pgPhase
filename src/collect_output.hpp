#ifndef PGPHASE_COLLECT_OUTPUT_HPP
#define PGPHASE_COLLECT_OUTPUT_HPP

#include "collect_types.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace pgphase_collect {

/**
 * TSV/VCF/read-support writers for collect-bam-variation.
 *
 * Contract: outputs describe **pre-phasing candidates**, not final diploid genotypes. TSV includes
 * CATEGORY (final, post-containment) and INIT_CAT (longcallD first `classify_var_cate` only, for parity
 * checks). VCF uses FILTER from final category; INFO.CAT matches final labels. PHASE_SET / hap in TSV
 * are placeholders (0). Read-support TSV is
 * optional auxiliary evidence for downstream phasing tools.
 */

/** @brief TSV/VCF label for variant type (SNP/INS/DEL). */
std::string type_name(VariantType type);
/** @brief Token for `VariantCategory` columns and VCF INFO.CAT. */
std::string category_name(VariantCategory category);
/** @brief Column header line for `--read-support` TSV. */
void write_read_support_header(std::ostream& out);
/** @brief Appends read×site rows for one batch. */
void write_read_support_rows(std::ostream& out,
                             const bam_hdr_t* header,
                             const std::vector<ReadSupportRow>& rows);
/** @brief Header line for main candidate TSV. */
void write_variants_tsv_header(std::ostream& out);
/** @brief One data row per candidate (REF/ALT from cache). */
void write_variants_tsv_records(std::ostream& out,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants);
/** @brief VCF 4.2 meta and contig lines through #CHROM. */
void write_variants_vcf_header(std::ostream& out, const Options& opts, const bam_hdr_t* header);
/** @brief VCF body records for all candidates. */
void write_variants_vcf_records(std::ostream& out,
                                const Options& opts,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants);
/** @brief Writes read-support TSV by reopening BAM for SQ names (batch mode helper). */
void write_read_support_tsv(const Options& opts, const std::vector<std::vector<ReadSupportRow>>& by_chunk);
/** @brief Writes full candidate TSV via `opts.output_tsv`. */
void write_variants(const Options& opts, faidx_t* fai, const CandidateTable& variants);
/** @brief Optional candidate VCF when `opts.output_vcf` is set. */
void write_variants_vcf(const Options& opts, faidx_t* fai, const CandidateTable& variants);

} // namespace pgphase_collect

#endif
