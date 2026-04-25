#ifndef PGPHASE_VARIANT_CLASSIFICATION_HPP
#define PGPHASE_VARIANT_CLASSIFICATION_HPP

#include "collect_types.hpp"

namespace pgphase_collect {

void classify_cand_vars_pgphase(BamChunk& chunk, const Options& opts);
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
