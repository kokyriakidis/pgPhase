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

constexpr int kPgbamThreadPolarityMargin = 2;
constexpr int kPgbamMinWinningThreads = 2;

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

static void collect_read_threads(const ReadRecord& read,
                                 const PgbamSidecarData& sidecar,
                                 std::vector<uint64_t>& threads) {
    threads.clear();
    std::vector<uint32_t> hs_ids;
    if (!extract_hs_set_ids(read, hs_ids) || hs_ids.empty()) return;
    for (uint32_t set_id : hs_ids) {
        auto it = sidecar.set_to_threads.find(set_id);
        if (it == sidecar.set_to_threads.end()) continue;
        threads.insert(threads.end(), it->second.begin(), it->second.end());
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

struct ThreadSupport {
    uint64_t thread = 0;
    int hap1_count = 0;
    int hap2_count = 0;
};

struct PhaseBlockThreadState {
    hts_pos_t anchor = -1;
    std::array<std::vector<uint64_t>, 3> hap_read_threads;
    std::vector<ThreadSupport> thread_support;
    std::array<std::vector<uint64_t>, 3> hap_polarized_threads;
};

static int count_thread_occurrences(const std::vector<uint64_t>& threads, size_t& i, uint64_t thread) {
    int count = 0;
    while (i < threads.size() && threads[i] == thread) {
        ++count;
        ++i;
    }
    return count;
}

static std::vector<ThreadSupport>
build_thread_support_from_hap_threads(std::vector<uint64_t>& hap1_threads,
                                      std::vector<uint64_t>& hap2_threads) {
    std::sort(hap1_threads.begin(), hap1_threads.end());
    std::sort(hap2_threads.begin(), hap2_threads.end());

    std::vector<ThreadSupport> out;
    out.reserve(hap1_threads.size() + hap2_threads.size());
    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < hap1_threads.size() || i2 < hap2_threads.size()) {
        uint64_t thread = 0;
        if (i2 == hap2_threads.size() || (i1 < hap1_threads.size() && hap1_threads[i1] < hap2_threads[i2])) {
            thread = hap1_threads[i1];
        } else {
            thread = hap2_threads[i2];
        }

        ThreadSupport s;
        s.thread = thread;
        if (i1 < hap1_threads.size() && hap1_threads[i1] == thread) {
            s.hap1_count = count_thread_occurrences(hap1_threads, i1, thread);
        }
        if (i2 < hap2_threads.size() && hap2_threads[i2] == thread) {
            s.hap2_count = count_thread_occurrences(hap2_threads, i2, thread);
        }
        out.push_back(s);
    }
    return out;
}

static void refresh_phase_block_polarized_threads(PhaseBlockThreadState& state, int polarity_margin) {
    state.hap_polarized_threads[1].clear();
    state.hap_polarized_threads[2].clear();
    state.hap_polarized_threads[1].reserve(state.thread_support.size());
    state.hap_polarized_threads[2].reserve(state.thread_support.size());
    for (const ThreadSupport& s : state.thread_support) {
        if (s.hap1_count - s.hap2_count >= polarity_margin) {
            state.hap_polarized_threads[1].push_back(s.thread);
        } else if (s.hap2_count - s.hap1_count >= polarity_margin) {
            state.hap_polarized_threads[2].push_back(s.thread);
        }
    }
}

static void finalize_phase_block_thread_state(PhaseBlockThreadState& state, int polarity_margin) {
    state.thread_support =
        build_thread_support_from_hap_threads(state.hap_read_threads[1], state.hap_read_threads[2]);
    state.hap_read_threads[1].clear();
    state.hap_read_threads[2].clear();
    refresh_phase_block_polarized_threads(state, polarity_margin);
}

static void build_read_thread_cache(const BamChunk& chunk,
                                    const PgbamSidecarData& sidecar,
                                    std::vector<std::vector<uint64_t>>& read_threads) {
    read_threads.clear();
    read_threads.resize(chunk.reads.size());
    std::vector<uint32_t> hs_ids;

    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped) continue;
        collect_read_threads(chunk.reads[i], sidecar, read_threads[i]);
    }
}

static std::vector<PhaseBlockThreadState>
build_phase_block_thread_states(const BamChunk& chunk,
                                const std::vector<hts_pos_t>& psets,
                                const std::vector<std::vector<uint64_t>>& read_threads,
                                int polarity_margin) {
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
        states[it->second].hap_read_threads[hap].insert(states[it->second].hap_read_threads[hap].end(),
                                                        threads.begin(), threads.end());
    }

    for (PhaseBlockThreadState& state : states) {
        finalize_phase_block_thread_state(state, polarity_margin);
    }
    return states;
}

static PhaseBlockThreadState collect_phase_block_thread_state(const BamChunk& chunk,
                                                              hts_pos_t phase_set,
                                                              const PgbamSidecarData& sidecar,
                                                              int polarity_margin) {
    PhaseBlockThreadState state;
    state.anchor = phase_set;
    std::vector<uint64_t> threads;
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        if (chunk.reads[i].is_skipped || chunk.phase_sets[i] != phase_set) continue;
        const int hap = chunk.haps[i];
        if (hap != 1 && hap != 2) continue;
        collect_read_threads(chunk.reads[i], sidecar, threads);
        state.hap_read_threads[hap].insert(state.hap_read_threads[hap].end(),
                                           threads.begin(), threads.end());
    }
    finalize_phase_block_thread_state(state, polarity_margin);
    return state;
}

static bool decide_phase_block_concordance(const PhaseBlockThreadState& left,
                                           const PhaseBlockThreadState& right,
                                           bool& do_flip,
                                           int min_winning_threads) {
    const int s11 = intersection_size(left.hap_polarized_threads[1], right.hap_polarized_threads[1]);
    const int s12 = intersection_size(left.hap_polarized_threads[1], right.hap_polarized_threads[2]);
    const int s21 = intersection_size(left.hap_polarized_threads[2], right.hap_polarized_threads[1]);
    const int s22 = intersection_size(left.hap_polarized_threads[2], right.hap_polarized_threads[2]);
    const int score_nonflip = s11 + s22;
    const int score_flip = s12 + s21;

    if (score_nonflip == score_flip) return false;
    if (std::max(score_nonflip, score_flip) < min_winning_threads) return false;
    do_flip = score_flip > score_nonflip;
    return true;
}

static ThreadSupport orient_thread_support(const ThreadSupport& support, bool flip) {
    if (!flip) return support;
    ThreadSupport out = support;
    std::swap(out.hap1_count, out.hap2_count);
    return out;
}

static PhaseBlockThreadState merge_phase_block_thread_states(const PhaseBlockThreadState& left,
                                                             const PhaseBlockThreadState& right,
                                                             bool do_flip,
                                                             int polarity_margin) {
    PhaseBlockThreadState merged;
    merged.anchor = left.anchor;
    merged.thread_support.reserve(left.thread_support.size() + right.thread_support.size());
    size_t li = 0;
    size_t ri = 0;
    while (li < left.thread_support.size() || ri < right.thread_support.size()) {
        if (ri == right.thread_support.size() ||
            (li < left.thread_support.size() &&
             left.thread_support[li].thread < right.thread_support[ri].thread)) {
            merged.thread_support.push_back(left.thread_support[li++]);
        } else if (li == left.thread_support.size() ||
                   right.thread_support[ri].thread < left.thread_support[li].thread) {
            merged.thread_support.push_back(orient_thread_support(right.thread_support[ri++], do_flip));
        } else {
            ThreadSupport out = left.thread_support[li++];
            const ThreadSupport right_oriented = orient_thread_support(right.thread_support[ri++], do_flip);
            out.hap1_count += right_oriented.hap1_count;
            out.hap2_count += right_oriented.hap2_count;
            merged.thread_support.push_back(out);
        }
    }
    refresh_phase_block_polarized_threads(merged, polarity_margin);
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

static void apply_pgbam_phase_merge(std::vector<BamChunk>& chunks,
                                    bool do_flip,
                                    hts_pos_t target_ps,
                                    hts_pos_t merged_ps) {
    for (BamChunk& chunk : chunks) {
        apply_pgbam_phase_merge(chunk, do_flip, target_ps, merged_ps);
    }
}

static void collect_phase_sets_from_chunks(const std::vector<BamChunk>& chunks,
                                           std::vector<hts_pos_t>& psets) {
    psets.clear();
    for (const BamChunk& chunk : chunks) {
        for (size_t i = 0; i < chunk.reads.size(); ++i) {
            if (chunk.reads[i].is_skipped || chunk.haps[i] == 0 || chunk.phase_sets[i] < 0) continue;
            psets.push_back(chunk.phase_sets[i]);
        }
    }
    std::sort(psets.begin(), psets.end());
    psets.erase(std::unique(psets.begin(), psets.end()), psets.end());
}

static std::vector<PhaseBlockThreadState>
build_phase_block_thread_states(const std::vector<BamChunk>& chunks,
                                const std::vector<hts_pos_t>& psets,
                                const PgbamSidecarData& sidecar,
                                int polarity_margin) {
    std::vector<PhaseBlockThreadState> states(psets.size());
    std::unordered_map<hts_pos_t, size_t> pset_to_state;
    pset_to_state.reserve(psets.size());
    for (size_t i = 0; i < psets.size(); ++i) {
        states[i].anchor = psets[i];
        pset_to_state[psets[i]] = i;
    }

    std::vector<uint64_t> threads;
    for (const BamChunk& chunk : chunks) {
        for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
            if (chunk.reads[read_i].is_skipped) continue;
            const int hap = chunk.haps[read_i];
            if (hap != 1 && hap != 2) continue;
            auto it = pset_to_state.find(chunk.phase_sets[read_i]);
            if (it == pset_to_state.end()) continue;
            collect_read_threads(chunk.reads[read_i], sidecar, threads);
            states[it->second].hap_read_threads[hap].insert(states[it->second].hap_read_threads[hap].end(),
                                                            threads.begin(), threads.end());
        }
    }

    for (PhaseBlockThreadState& state : states) {
        finalize_phase_block_thread_state(state, polarity_margin);
    }
    return states;
}

static void stitch_phase_blocks_with_pgbam(std::vector<BamChunk>& chunks,
                                           const PgbamSidecarData& sidecar,
                                           int min_winning_threads,
                                           int polarity_margin,
                                           std::vector<hts_pos_t>& psets,
                                           std::vector<PhaseBlockThreadState>& states) {
    (void)sidecar;
    for (size_t i = 1; i < psets.size(); ++i) {
        const hts_pos_t left_ps = psets[i - 1];
        const hts_pos_t right_ps = psets[i];
        if (left_ps == right_ps) continue;

        bool do_flip = false;
        if (!decide_phase_block_concordance(states[i - 1], states[i], do_flip, min_winning_threads)) continue;
        apply_pgbam_phase_merge(chunks, do_flip, right_ps, left_ps);
        states[i] = merge_phase_block_thread_states(states[i - 1], states[i], do_flip, polarity_margin);
        psets[i] = left_ps;
    }
}

void stitch_phase_blocks_with_pgbam(BamChunk& chunk, const PgbamSidecarData& sidecar) {
    stitch_phase_blocks_with_pgbam(chunk, sidecar, kPgbamMinWinningThreads, kPgbamThreadPolarityMargin);
}

void stitch_phase_blocks_with_pgbam(BamChunk& chunk,
                                    const PgbamSidecarData& sidecar,
                                    int min_winning_threads,
                                    int polarity_margin) {
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
        build_phase_block_thread_states(chunk, psets, read_threads, polarity_margin);

    for (size_t i = 1; i < psets.size(); ++i) {
        const hts_pos_t left_ps = psets[i - 1];
        const hts_pos_t right_ps = psets[i];

        bool do_flip = false;
        if (!decide_phase_block_concordance(states[i - 1], states[i], do_flip, min_winning_threads)) continue;
        apply_pgbam_phase_merge(chunk, do_flip, right_ps, left_ps);
        states[i] = merge_phase_block_thread_states(states[i - 1], states[i], do_flip, polarity_margin);
        psets[i] = left_ps;
    }
}

void stitch_phase_blocks_with_pgbam(std::vector<BamChunk>& chunks,
                                    const PgbamSidecarData& sidecar,
                                    int min_winning_threads) {
    stitch_phase_blocks_with_pgbam(chunks, sidecar, min_winning_threads, kPgbamThreadPolarityMargin);
}

void stitch_phase_blocks_with_pgbam(std::vector<BamChunk>& chunks,
                                    const PgbamSidecarData& sidecar,
                                    int min_winning_threads,
                                    int polarity_margin) {
    std::vector<hts_pos_t> psets;
    collect_phase_sets_from_chunks(chunks, psets);
    if (psets.size() < 2) return;

    std::vector<PhaseBlockThreadState> states =
        build_phase_block_thread_states(chunks, psets, sidecar, polarity_margin);
    stitch_phase_blocks_with_pgbam(chunks, sidecar, min_winning_threads, polarity_margin, psets, states);
}

bool stitch_adjacent_chunks_with_pgbam(BamChunk& pre, BamChunk& cur, const PgbamSidecarData& sidecar) {
    return stitch_adjacent_chunks_with_pgbam(pre, cur, sidecar, kPgbamMinWinningThreads, kPgbamThreadPolarityMargin);
}

bool stitch_adjacent_chunks_with_pgbam(BamChunk& pre,
                                       BamChunk& cur,
                                       const PgbamSidecarData& sidecar,
                                       int min_winning_threads,
                                       int polarity_margin) {
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

    const PhaseBlockThreadState pre_state =
        collect_phase_block_thread_state(pre, max_pre_ps, sidecar, polarity_margin);
    const PhaseBlockThreadState cur_state =
        collect_phase_block_thread_state(cur, min_cur_ps, sidecar, polarity_margin);
    bool do_flip = false;
    if (!decide_phase_block_concordance(pre_state, cur_state, do_flip, min_winning_threads)) return false;
    apply_pgbam_phase_merge(cur, do_flip, min_cur_ps, max_pre_ps);
    return true;
}

} // namespace pgphase_collect
