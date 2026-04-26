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
 * Contract: outputs describe **pre-phasing candidates**, not final diploid genotypes. TSV includes all
 * categories (filter with CATEGORY). VCF uses FILTER=LowQual/RefCall/NoCall for non-clean rows; INFO.CAT
 * holds the same labels as TSV. PHASE_SET / hap columns in TSV are placeholders (0). Read-support TSV is
 * optional auxiliary evidence for downstream phasing tools.
 */

std::string type_name(VariantType type);
std::string category_name(VariantCategory category);
void write_read_support_header(std::ostream& out);
void write_read_support_rows(std::ostream& out,
                             const bam_hdr_t* header,
                             const std::vector<ReadSupportRow>& rows);
void write_variants_tsv_header(std::ostream& out);
void write_variants_tsv_records(std::ostream& out,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants);
void write_variants_vcf_header(std::ostream& out, const Options& opts, const bam_hdr_t* header);
void write_variants_vcf_records(std::ostream& out,
                                const Options& opts,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants);
void write_read_support_tsv(const Options& opts, const std::vector<std::vector<ReadSupportRow>>& by_chunk);
void write_variants(const Options& opts, faidx_t* fai, const CandidateTable& variants);
void write_variants_vcf(const Options& opts, faidx_t* fai, const CandidateTable& variants);

} // namespace pgphase_collect

#endif
