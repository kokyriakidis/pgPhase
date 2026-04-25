#ifndef PGPHASE_COLLECT_OUTPUT_HPP
#define PGPHASE_COLLECT_OUTPUT_HPP

#include "collect_types.hpp"

#include <string>
#include <vector>

namespace pgphase_collect {

std::string type_name(VariantType type);
std::string category_name(VariantCategory category);
void write_read_support_tsv(const Options& opts, const std::vector<std::vector<ReadSupportRow>>& by_chunk);
void write_variants(const Options& opts, faidx_t* fai, const CandidateTable& variants);
void write_variants_vcf(const Options& opts, faidx_t* fai, const CandidateTable& variants);

} // namespace pgphase_collect

#endif
