#ifndef PGPHASE_CANDIDATE_COLLECTION_HPP
#define PGPHASE_CANDIDATE_COLLECTION_HPP

#include "collect_types.hpp"

#include <vector>

namespace pgphase_collect {

VariantKey variant_key_from_digar(int tid, const DigarOp& op);
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2);
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len);

struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const;
};

void collapse_fuzzy_large_insertions(CandidateTable& variants);
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants);
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq);

} // namespace pgphase_collect

#endif
