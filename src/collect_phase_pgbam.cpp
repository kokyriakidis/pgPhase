/**
 * @file collect_phase_pgbam.cpp
 * @brief Annotated-BAM/.pgbam thread support for phase-block stitching.
 */

#include "collect_phase_pgbam.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <iterator>
#include <unordered_map>

#include <htslib/sam.h>

namespace pgphase_collect {

static uint32_t aux_read_u32_le(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

static int extract_hs_set_ids(const ReadRecord& read, std::vector<uint32_t>& out_ids) {
    out_ids.clear();
    if (!read.alignment) return 0;
    uint8_t* hs = bam_aux_get(read.alignment.get(), "hs");
    if (hs == nullptr || hs[0] != 'B') return 0;
    const uint8_t subtype = hs[1];
    const uint32_t n = aux_read_u32_le(hs + 2);
    const uint8_t* data = hs + 6;
    out_ids.reserve(static_cast<size_t>(n));

    if (subtype == 'I' || subtype == 'i') {
        for (uint32_t i = 0; i < n; ++i) out_ids.push_back(aux_read_u32_le(data + i * 4));
        return 1;
    }
    if (subtype == 'S' || subtype == 's') {
        for (uint32_t i = 0; i < n; ++i) {
            const uint8_t* p = data + i * 2;
            out_ids.push_back(static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8));
        }
        return 1;
    }
    if (subtype == 'C' || subtype == 'c') {
        for (uint32_t i = 0; i < n; ++i) out_ids.push_back(static_cast<uint32_t>(data[i]));
        return 1;
    }
    return 0;
}

static void sort_unique_threads(std::vector<uint64_t>& threads) {
    if (threads.empty()) return;
    std::sort(threads.begin(), threads.end());
    threads.erase(std::unique(threads.begin(), threads.end()), threads.end());
}

static void collect_thread_set_for_group(const BamChunk& chunk,
                                         int hap,
                                         hts_pos_t phase_set,
                                         const PgbamSidecarData& sidecar,
                                         std::vector<uint64_t>& threads) {
    threads.clear();
    std::vector<uint32_t> hs_ids;
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped || chunk.haps[i] != hap || chunk.phase_sets[i] != phase_set) continue;
        if (!extract_hs_set_ids(chunk.reads[i], hs_ids) || hs_ids.empty()) continue;
        for (uint32_t set_id : hs_ids) {
            auto it = sidecar.set_to_threads.find(set_id);
            if (it == sidecar.set_to_threads.end()) continue;
            threads.insert(threads.end(), it->second.begin(), it->second.end());
        }
    }
    sort_unique_threads(threads);
}

static int intersection_size(const std::vector<uint64_t>& lhs, const std::vector<uint64_t>& rhs) {
    int n = 0;
    auto li = lhs.begin();
    auto ri = rhs.begin();
    while (li != lhs.end() && ri != rhs.end()) {
        if (*li < *ri) ++li;
        else if (*ri < *li) ++ri;
        else {
            ++n;
            ++li;
            ++ri;
        }
    }
    return n;
}

static std::vector<uint64_t> set_difference_threads(const std::vector<uint64_t>& lhs,
                                                    const std::vector<uint64_t>& rhs) {
    std::vector<uint64_t> out;
    out.reserve(lhs.size());
    auto li = lhs.begin();
    auto ri = rhs.begin();
    while (li != lhs.end()) {
        while (ri != rhs.end() && *ri < *li) ++ri;
        if (ri == rhs.end() || *li < *ri) out.push_back(*li);
        ++li;
    }
    return out;
}

static std::vector<uint64_t> union_threads(const std::vector<uint64_t>& lhs,
                                           const std::vector<uint64_t>& rhs) {
    std::vector<uint64_t> out;
    out.reserve(lhs.size() + rhs.size());
    std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(out));
    return out;
}

struct PhaseBlockThreadState {
    hts_pos_t anchor = -1;
    std::array<std::vector<uint64_t>, 3> hap_threads;
    std::array<std::vector<uint64_t>, 3> hap_unique_threads;
};

static void refresh_phase_block_unique_threads(PhaseBlockThreadState& state) {
    state.hap_unique_threads[1] = set_difference_threads(state.hap_threads[1], state.hap_threads[2]);
    state.hap_unique_threads[2] = set_difference_threads(state.hap_threads[2], state.hap_threads[1]);
}

static void build_read_thread_cache(const BamChunk& chunk,
                                    const PgbamSidecarData& sidecar,
                                    std::vector<std::vector<uint64_t>>& read_threads) {
    read_threads.clear();
    read_threads.resize(chunk.reads.size());
    std::vector<uint32_t> hs_ids;

    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped) continue;
        if (!extract_hs_set_ids(chunk.reads[i], hs_ids) || hs_ids.empty()) continue;
        std::vector<uint64_t>& threads = read_threads[i];
        for (uint32_t set_id : hs_ids) {
            auto it = sidecar.set_to_threads.find(set_id);
            if (it == sidecar.set_to_threads.end()) continue;
            threads.insert(threads.end(), it->second.begin(), it->second.end());
        }
        sort_unique_threads(threads);
    }
}

static std::vector<PhaseBlockThreadState>
build_phase_block_thread_states(const BamChunk& chunk,
                                const std::vector<hts_pos_t>& psets,
                                const std::vector<std::vector<uint64_t>>& read_threads) {
    std::vector<PhaseBlockThreadState> states(psets.size());
    std::unordered_map<hts_pos_t, size_t> pset_to_state;
    pset_to_state.reserve(psets.size());
    for (size_t i = 0; i < psets.size(); ++i) {
        states[i].anchor = psets[i];
        pset_to_state[psets[i]] = i;
    }

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        if (chunk.reads[read_i].is_skipped) continue;
        const int hap = chunk.haps[read_i];
        if (hap != 1 && hap != 2) continue;
        auto it = pset_to_state.find(chunk.phase_sets[read_i]);
        if (it == pset_to_state.end()) continue;
        const std::vector<uint64_t>& threads = read_threads[read_i];
        states[it->second].hap_threads[hap].insert(states[it->second].hap_threads[hap].end(),
                                                   threads.begin(), threads.end());
    }

    for (PhaseBlockThreadState& state : states) {
        sort_unique_threads(state.hap_threads[1]);
        sort_unique_threads(state.hap_threads[2]);
        refresh_phase_block_unique_threads(state);
    }
    return states;
}

static PhaseBlockThreadState collect_phase_block_thread_state(const BamChunk& chunk,
                                                              hts_pos_t phase_set,
                                                              const PgbamSidecarData& sidecar) {
    PhaseBlockThreadState state;
    state.anchor = phase_set;
    collect_thread_set_for_group(chunk, 1, phase_set, sidecar, state.hap_threads[1]);
    collect_thread_set_for_group(chunk, 2, phase_set, sidecar, state.hap_threads[2]);
    refresh_phase_block_unique_threads(state);
    return state;
}

static bool decide_phase_block_concordance(const PhaseBlockThreadState& left,
                                           const PhaseBlockThreadState& right,
                                           bool& do_flip) {
    const int s11 = intersection_size(left.hap_unique_threads[1], right.hap_unique_threads[1]);
    const int s12 = intersection_size(left.hap_unique_threads[1], right.hap_unique_threads[2]);
    const int s21 = intersection_size(left.hap_unique_threads[2], right.hap_unique_threads[1]);
    const int s22 = intersection_size(left.hap_unique_threads[2], right.hap_unique_threads[2]);
    const int score_nonflip = s11 + s22;
    const int score_flip = s12 + s21;

    if (score_nonflip == score_flip) return false;
    if (std::max(score_nonflip, score_flip) < 2) return false;
    do_flip = score_flip > score_nonflip;
    return true;
}

static PhaseBlockThreadState merge_phase_block_thread_states(const PhaseBlockThreadState& left,
                                                             const PhaseBlockThreadState& right,
                                                             bool do_flip) {
    PhaseBlockThreadState merged;
    merged.anchor = left.anchor;
    const int right_hap_for_1 = do_flip ? 2 : 1;
    const int right_hap_for_2 = do_flip ? 1 : 2;
    merged.hap_threads[1] = union_threads(left.hap_threads[1], right.hap_threads[right_hap_for_1]);
    merged.hap_threads[2] = union_threads(left.hap_threads[2], right.hap_threads[right_hap_for_2]);
    refresh_phase_block_unique_threads(merged);
    return merged;
}

static void apply_pgbam_phase_merge(BamChunk& cur,
                                    bool do_flip,
                                    hts_pos_t target_ps,
                                    hts_pos_t merged_ps) {
    if (do_flip && target_ps != INT64_MAX && target_ps != static_cast<hts_pos_t>(-1)) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set != target_ps) continue;
            std::swap(v.hap_to_cons_alle[1], v.hap_to_cons_alle[2]);
        }
        for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
            if (cur.reads[read_i].is_skipped || cur.haps[read_i] == 0) continue;
            if (cur.phase_sets[read_i] == target_ps) cur.haps[read_i] = 3 - cur.haps[read_i];
        }
    }
    if (merged_ps != -1 && target_ps != INT64_MAX) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set == target_ps) v.phase_set = merged_ps;
        }
        for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
            if (cur.phase_sets[read_i] == target_ps) cur.phase_sets[read_i] = merged_ps;
        }
    }
}

void stitch_phase_blocks_with_pgbam(BamChunk& chunk, const PgbamSidecarData& sidecar) {
    std::vector<hts_pos_t> psets;
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped || chunk.haps[i] == 0 || chunk.phase_sets[i] < 0) continue;
        psets.push_back(chunk.phase_sets[i]);
    }
    std::sort(psets.begin(), psets.end());
    psets.erase(std::unique(psets.begin(), psets.end()), psets.end());
    if (psets.size() < 2) return;

    std::vector<std::vector<uint64_t>> read_threads;
    build_read_thread_cache(chunk, sidecar, read_threads);
    std::vector<PhaseBlockThreadState> states =
        build_phase_block_thread_states(chunk, psets, read_threads);

    for (size_t i = 1; i < psets.size(); ++i) {
        const hts_pos_t left_ps = psets[i - 1];
        const hts_pos_t right_ps = psets[i];

        bool do_flip = false;
        if (!decide_phase_block_concordance(states[i - 1], states[i], do_flip)) continue;
        apply_pgbam_phase_merge(chunk, do_flip, right_ps, left_ps);
        states[i] = merge_phase_block_thread_states(states[i - 1], states[i], do_flip);
        psets[i] = left_ps;
    }
}

bool stitch_adjacent_chunks_with_pgbam(BamChunk& pre, BamChunk& cur, const PgbamSidecarData& sidecar) {
    if (pre.region.tid != cur.region.tid) return false;

    hts_pos_t max_pre_ps = -1;
    hts_pos_t min_cur_ps = INT64_MAX;
    for (size_t i = 0; i < pre.reads.size(); ++i) {
        if (pre.reads[i].is_skipped || pre.haps[i] == 0 || pre.phase_sets[i] < 0) continue;
        if (pre.phase_sets[i] > max_pre_ps) max_pre_ps = pre.phase_sets[i];
    }
    for (size_t i = 0; i < cur.reads.size(); ++i) {
        if (cur.reads[i].is_skipped || cur.haps[i] == 0 || cur.phase_sets[i] < 0) continue;
        if (cur.phase_sets[i] < min_cur_ps) min_cur_ps = cur.phase_sets[i];
    }
    if (max_pre_ps < 0 || min_cur_ps == INT64_MAX) return false;

    const PhaseBlockThreadState pre_state = collect_phase_block_thread_state(pre, max_pre_ps, sidecar);
    const PhaseBlockThreadState cur_state = collect_phase_block_thread_state(cur, min_cur_ps, sidecar);
    bool do_flip = false;
    if (!decide_phase_block_concordance(pre_state, cur_state, do_flip)) return false;
    apply_pgbam_phase_merge(cur, do_flip, min_cur_ps, max_pre_ps);
    return true;
}

} // namespace pgphase_collect
