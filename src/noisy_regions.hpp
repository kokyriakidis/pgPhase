#ifndef PGPHASE_NOISY_REGIONS_HPP
#define PGPHASE_NOISY_REGIONS_HPP

#include "collect_types.hpp"

#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

struct CrangesOwner {
    cgranges_t* cr = nullptr;
    CrangesOwner() = default;
    CrangesOwner(const CrangesOwner&) = delete;
    CrangesOwner& operator=(const CrangesOwner&) = delete;
    CrangesOwner(CrangesOwner&& other) noexcept : cr(other.cr) { other.cr = nullptr; }
    CrangesOwner& operator=(CrangesOwner&& other) noexcept {
        if (this != &other) {
            reset();
            cr = other.cr;
            other.cr = nullptr;
        }
        return *this;
    }
    ~CrangesOwner() { reset(); }
    void reset() {
        if (cr != nullptr) {
            cr_destroy(cr);
            cr = nullptr;
        }
    }
    void adopt(cgranges_t* p) {
        reset();
        if (p == nullptr) return;
        cr = p;
    }
    cgranges_t* release() {
        cgranges_t* t = cr;
        cr = nullptr;
        return t;
    }
};

void merge_intervals(std::vector<Interval>& intervals);
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals);
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out);
cgranges_t* build_read_noisy_cr(const ReadRecord& read);
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end);

void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts);
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand);
void apply_noisy_containment_filter(BamChunk& chunk);
void populate_low_complexity_intervals(BamChunk& chunk);
void populate_chunk_read_indexes(BamChunk& chunk);

} // namespace pgphase_collect

#endif
