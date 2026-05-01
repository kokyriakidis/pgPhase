/**
 * @file test_phase_block_stitch.cpp
 * @brief Unit tests for longcallD-equivalent chunk stitching (`stitch_chunk_haps` / `flip_variant_hap`).
 *
 * Only boundary overlap read pairing is exercised — matching longcallD `stitch_var_main`, which does
 * not run pangenome-graph or within-chunk merges.
 *
 * Build (from repo root):
 *   make collect_phase.o && g++ -O0 -g -std=c++17 -Wall -Wextra \
 *       src/test_phase_block_stitch.cpp collect_phase.o -lhts -lz -lpthread \
 *       -o test_phase_block_stitch && ./test_phase_block_stitch
 */

#include "collect_phase.hpp"
#include "collect_types.hpp"

#include <cstdio>
#include <vector>

using namespace pgphase_collect;

static bool check(bool cond, const char* msg) {
    if (!cond) {
        std::printf("FAIL: %s\n", msg);
        return false;
    }
    return true;
}

static ReadRecord min_read() {
    ReadRecord r;
    r.is_skipped = false;
    return r;
}

static CandidateVariant dummy_cand(hts_pos_t ps) {
    CandidateVariant v;
    v.phase_set = ps;
    v.hap_to_cons_alle = {-1, 1, 0};
    return v;
}

// Two chunks, one overlapping pair: pre down index 0, cur up index 0.
static bool test_cross_chunk_flip_when_haps_disagree() {
    std::printf("--- test_cross_chunk_flip_when_haps_disagree ---\n");
    std::vector<BamChunk> chunks;
    chunks.resize(2);
    BamChunk& pre = chunks[0];
    BamChunk& cur = chunks[1];
    pre.region.tid = 0;
    cur.region.tid = 0;

    pre.reads.push_back(min_read());
    pre.haps.push_back(1);
    pre.phase_sets.push_back(500);
    pre.candidates.push_back(dummy_cand(500));
    pre.down_ovlp_read_i = {{0}};
    pre.up_ovlp_read_i = {{}};

    cur.reads.push_back(min_read());
    cur.haps.push_back(2);
    cur.phase_sets.push_back(900);
    cur.candidates.push_back(dummy_cand(900));
    cur.up_ovlp_read_i = {{0}};
    cur.down_ovlp_read_i = {{}};

    Options opts;
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = true;
    const CandidateVariant& v = chunks[1].candidates[0];
    ok &= check(v.phase_set == 500, "phase_set merged to pre anchor");
    ok &= check(v.hap_to_cons_alle[1] == 0 && v.hap_to_cons_alle[2] == 1, "hap consensus swapped");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// Equal votes: one agree (−1) and one disagree (+1) ⇒ flip_hap_score == 0 ⇒ no updates.
static bool test_flip_score_zero_no_merge() {
    std::printf("--- test_flip_score_zero_no_merge ---\n");
    BamChunk pre, cur;
    pre.region.tid = cur.region.tid = 0;

    for (int i = 0; i < 2; ++i) {
        pre.reads.push_back(min_read());
        pre.haps.push_back(1);
        pre.phase_sets.push_back(100);
    }
    pre.candidates.push_back(dummy_cand(100));
    pre.down_ovlp_read_i = {{0, 1}};

    for (int i = 0; i < 2; ++i) {
        cur.reads.push_back(min_read());
        cur.haps.push_back(0);
        cur.phase_sets.push_back(200);
    }
    cur.haps[0] = 1;  // agree with pre
    cur.haps[1] = 2;  // disagree
    cur.candidates.push_back(dummy_cand(200));
    cur.up_ovlp_read_i = {{0, 1}};

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(pre));
    chunks.push_back(std::move(cur));

    Options opts;
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = true;
    ok &= check(chunks[1].candidates[0].phase_set == 200, "phase_set unchanged when flip score ties");
    ok &= check(chunks[1].candidates[0].hap_to_cons_alle[1] == 1 &&
                    chunks[1].candidates[0].hap_to_cons_alle[2] == 0,
                "hap consensus unchanged when flip score ties");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_skipped_overlap_read_ignored() {
    std::printf("--- test_skipped_overlap_read_ignored ---\n");
    BamChunk pre, cur;
    pre.region.tid = cur.region.tid = 0;

    ReadRecord rskip = min_read();
    rskip.is_skipped = true;
    pre.reads.push_back(std::move(rskip));
    pre.haps.push_back(1);
    pre.phase_sets.push_back(10);
    pre.candidates.push_back(dummy_cand(10));
    pre.down_ovlp_read_i = {{0}};

    cur.reads.push_back(min_read());
    cur.haps.push_back(2);
    cur.phase_sets.push_back(20);
    cur.candidates.push_back(dummy_cand(20));
    cur.up_ovlp_read_i = {{0}};

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(pre));
    chunks.push_back(std::move(cur));

    Options opts;
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = check(chunks[1].candidates[0].phase_set == 20, "no informative overlap reads → no merge");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

int main() {
    int failures = 0;
    failures += test_cross_chunk_flip_when_haps_disagree() ? 0 : 1;
    failures += test_flip_score_zero_no_merge() ? 0 : 1;
    failures += test_skipped_overlap_read_ignored() ? 0 : 1;
    if (failures == 0)
        std::printf("ALL PASS\n");
    else
        std::printf("FAILURES: %d\n", failures);
    return failures == 0 ? 0 : 1;
}
