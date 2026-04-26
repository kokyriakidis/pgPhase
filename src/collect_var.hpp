#ifndef PGPHASE_COLLECT_VAR_HPP
#define PGPHASE_COLLECT_VAR_HPP

#include "collect_types.hpp"

#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

// ── RAII wrapper for cgranges_t ──────────────────────────────────────────────
struct CrangesOwner {
    cgranges_t* cr = nullptr;
    CrangesOwner() = default;
    CrangesOwner(const CrangesOwner&) = delete;
    CrangesOwner& operator=(const CrangesOwner&) = delete;
    CrangesOwner(CrangesOwner&& other) noexcept : cr(other.cr) { other.cr = nullptr; }
    CrangesOwner& operator=(CrangesOwner&& other) noexcept {
        if (this != &other) { reset(); cr = other.cr; other.cr = nullptr; }
        return *this;
    }
    ~CrangesOwner() { reset(); }
    void reset() { if (cr) { cr_destroy(cr); cr = nullptr; } }
    void adopt(cgranges_t* p) { reset(); if (p) cr = p; }
    cgranges_t* release() { cgranges_t* t = cr; cr = nullptr; return t; }
};

// ── Variant site comparison ──────────────────────────────────────────────────
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2);
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len);

struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const;
};

// ── Candidate site collection ────────────────────────────────────────────────
VariantKey variant_key_from_digar(int tid, const DigarOp& op);
void collapse_fuzzy_large_insertions(CandidateTable& variants);
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants);

// ── Allele counting ──────────────────────────────────────────────────────────
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq);

// ── Interval utilities ───────────────────────────────────────────────────────
void merge_intervals(std::vector<Interval>& intervals);
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals);
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out);
cgranges_t* build_read_noisy_cr(const ReadRecord& read);
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end);

// ── Noisy region processing ──────────────────────────────────────────────────
void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts);
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand);
void apply_noisy_containment_filter(BamChunk& chunk);
void populate_low_complexity_intervals(BamChunk& chunk);
void populate_chunk_read_indexes(BamChunk& chunk);

// ── Variant classification ───────────────────────────────────────────────────
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header);

} // namespace pgphase_collect

#endif
