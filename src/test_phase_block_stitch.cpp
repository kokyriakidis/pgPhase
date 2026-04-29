/**
 * @file test_phase_block_stitch.cpp
 * @brief Unit test for within-chunk and cross-chunk phase block stitching via stitch_chunk_haps.
 *
 * Builds synthetic BamChunk objects with real bam1_t records carrying 'hs' B/I array tags,
 * then verifies that stitch_chunk_haps correctly orients adjacent phase blocks using the
 * pgbam per-read voting fallback.
 *
 * Build:
 *   g++ -O0 -g -std=c++17 -Wall -Wextra \
 *       src/test_phase_block_stitch.cpp src/collect_phase.o \
 *       -lhts -lm -lz -lpthread -o test_phase_block_stitch
 *   ./test_phase_block_stitch
 */

#include "collect_phase.hpp"
#include "collect_types.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

extern "C" {
#include <htslib/sam.h>
}

using namespace pgphase_collect;

// ─── helpers ────────────────────────────────────────────────────────────────

static void write_le32(uint8_t* dst, uint32_t v) {
    dst[0] =  v        & 0xFF;
    dst[1] = (v >>  8) & 0xFF;
    dst[2] = (v >> 16) & 0xFF;
    dst[3] = (v >> 24) & 0xFF;
}

// Create a bam1_t carrying an hs B/I array tag with the given set IDs.
static std::unique_ptr<bam1_t, AlignmentDeleter> make_bam_hs(const std::vector<uint32_t>& ids) {
    bam1_t* b = bam_init1();

    // Manually build the data payload for bam_aux_append('B'):
    //   byte 0:   subtype 'I'
    //   bytes 1-4: count (uint32 LE)
    //   bytes 5+:  elements (uint32 LE each)
    const uint32_t n = static_cast<uint32_t>(ids.size());
    std::vector<uint8_t> payload(1 + 4 + n * 4);
    payload[0] = 'I';
    write_le32(payload.data() + 1, n);
    for (uint32_t i = 0; i < n; ++i)
        write_le32(payload.data() + 5 + i * 4, ids[i]);

    bam_aux_append(b, "hs", 'B', static_cast<int>(payload.size()), payload.data());
    return std::unique_ptr<bam1_t, AlignmentDeleter>(b);
}

// Append a read to chunk.reads / haps / phase_sets.
static void push_read(BamChunk& chunk, int hap, hts_pos_t ps, const std::vector<uint32_t>& hs_ids) {
    ReadRecord r;
    r.is_skipped  = false;
    r.alignment   = make_bam_hs(hs_ids);
    chunk.reads.push_back(std::move(r));
    chunk.haps.push_back(hap);
    chunk.phase_sets.push_back(ps);
}

// Build opts with a non-empty pgbam_file so the pgbam path is enabled.
static Options make_opts() {
    Options o;
    o.pgbam_file = "test.pgbam";
    return o;
}

// Sidecar:  set_id 1→thread 10, 2→11, 3→20, 4→21
//   group A threads: {10, 11}  (set_ids {1, 2})
//   group B threads: {20, 21}  (set_ids {3, 4})
static PgbamSidecarData make_sidecar() {
    PgbamSidecarData s;
    s.set_to_threads[1] = {10};
    s.set_to_threads[2] = {11};
    s.set_to_threads[3] = {20};
    s.set_to_threads[4] = {21};
    return s;
}

static bool check(bool cond, const char* msg) {
    if (!cond) { std::printf("FAIL: %s\n", msg); return false; }
    return true;
}

// ─── test 1: within-chunk flip ──────────────────────────────────────────────
//
// One chunk, two phase blocks:
//   PS=1000  hap1→group A, hap2→group B
//   PS=2000  hap1→group B, hap2→group A   ← FLIPPED relative to PS=1000
//
// Expected after stitching:
//   PS=2000 reads flipped: former hap1→hap2, former hap2→hap1
//   All PS=2000 phase_sets renamed to 1000
static bool test_within_chunk_flip() {
    std::printf("--- test_within_chunk_flip ---\n");
    auto sidecar = make_sidecar();
    auto opts    = make_opts();

    BamChunk chunk;
    chunk.region.tid = 0;
    // no overlap reads needed (pgbam path only)

    // PS=1000: 3× hap1 (group A), 3× hap2 (group B)
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 1000, {1, 2}); // threads {10,11}
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 1000, {3, 4}); // threads {20,21}

    // PS=2000: 3× hap1 (group B) — WRONG orientation, should flip to hap2
    //          3× hap2 (group A) — WRONG orientation, should flip to hap1
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 2000, {3, 4}); // threads {20,21}
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 2000, {1, 2}); // threads {10,11}

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(chunk));
    stitch_chunk_haps(chunks, &opts, &sidecar);

    BamChunk& c = chunks[0];
    bool ok = true;

    // PS=1000 reads (idx 0-5): unchanged
    for (int i = 0; i < 3; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "PS1000 hap1 unchanged");
    for (int i = 3; i < 6; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "PS1000 hap2 unchanged");

    // PS=2000 reads (idx 6-8): were hap=1, should be hap=2 after flip, ps=1000
    for (int i = 6; i < 9; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "PS2000 former-hap1 → hap2 ps1000");

    // PS=2000 reads (idx 9-11): were hap=2, should be hap=1 after flip, ps=1000
    for (int i = 9; i < 12; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "PS2000 former-hap2 → hap1 ps1000");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── test 2: within-chunk nonflip ───────────────────────────────────────────
//
// Same setup but PS=2000 is already in the CORRECT orientation.
// Expected: haps unchanged, phase_sets renamed to 1000.
static bool test_within_chunk_nonflip() {
    std::printf("--- test_within_chunk_nonflip ---\n");
    auto sidecar = make_sidecar();
    auto opts    = make_opts();

    BamChunk chunk;
    chunk.region.tid = 0;

    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 1000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 1000, {3, 4});
    // PS=2000 same orientation: hap1→group A, hap2→group B
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 2000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 2000, {3, 4});

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(chunk));
    stitch_chunk_haps(chunks, &opts, &sidecar);

    BamChunk& c = chunks[0];
    bool ok = true;

    for (int i = 0; i < 3; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "PS1000 hap1 unchanged");
    for (int i = 3; i < 6; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "PS1000 hap2 unchanged");
    // PS=2000 → merged to 1000, haps unchanged
    for (int i = 6; i < 9; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "PS2000 hap1 stays hap1 ps1000");
    for (int i = 9; i < 12; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "PS2000 hap2 stays hap2 ps1000");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── test 3: within-chunk chain (3 phase blocks, middle flipped) ─────────────
//
// PS=1000 correct orientation.
// PS=2000 FLIPPED → after stitching with 1000, merged as corrected.
// PS=3000 correct orientation relative to PS=1000.
// The third stitch compares PS=3000 against the merged (1000+2000) block.
// PS=3000 should be merged into 1000 without flipping.
static bool test_within_chunk_chain() {
    std::printf("--- test_within_chunk_chain ---\n");
    auto sidecar = make_sidecar();
    auto opts    = make_opts();

    BamChunk chunk;
    chunk.region.tid = 0;

    // PS=1000: hap1→A, hap2→B
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 1000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 1000, {3, 4});
    // PS=2000: FLIPPED (hap1→B, hap2→A)
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 2000, {3, 4});
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 2000, {1, 2});
    // PS=3000: correct orientation (hap1→A, hap2→B)
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 3000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 3000, {3, 4});

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(chunk));
    stitch_chunk_haps(chunks, &opts, &sidecar);

    BamChunk& c = chunks[0];
    bool ok = true;

    // PS=1000 (idx 0-5): unchanged
    for (int i = 0; i < 3; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "chain PS1000 hap1");
    for (int i = 3; i < 6; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "chain PS1000 hap2");

    // PS=2000 (idx 6-11): flipped and merged
    for (int i = 6; i < 9; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "chain PS2000 former-hap1 → hap2");
    for (int i = 9; i < 12; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "chain PS2000 former-hap2 → hap1");

    // PS=3000 (idx 12-17): merged nonflip (stitched against the now-merged block which includes
    //   corrected PS=2000 reads; the dominant thread signal is still group A for hap1)
    for (int i = 12; i < 15; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "chain PS3000 hap1 nonflip");
    for (int i = 15; i < 18; ++i)
        ok &= check(c.haps[i] == 2 && c.phase_sets[i] == 1000, "chain PS3000 hap2 nonflip");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── test 4: cross-chunk flip via pgbam fallback ─────────────────────────────
//
// Two chunks, each with one phase block, no overlap reads.
// chunk0 PS=1000: hap1→A, hap2→B
// chunk1 PS=2000: hap1→B, hap2→A  ← FLIPPED
// Expected: chunk1 PS=2000 flipped and merged to 1000.
static bool test_cross_chunk_flip() {
    std::printf("--- test_cross_chunk_flip ---\n");
    auto sidecar = make_sidecar();
    auto opts    = make_opts();

    BamChunk c0, c1;
    c0.region.tid = 0;
    c1.region.tid = 0;

    for (int i = 0; i < 3; ++i) push_read(c0, 1, 1000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(c0, 2, 1000, {3, 4});

    for (int i = 0; i < 3; ++i) push_read(c1, 1, 2000, {3, 4}); // FLIPPED
    for (int i = 0; i < 3; ++i) push_read(c1, 2, 2000, {1, 2}); // FLIPPED

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(c0));
    chunks.push_back(std::move(c1));
    stitch_chunk_haps(chunks, &opts, &sidecar);

    bool ok = true;
    BamChunk& pre = chunks[0];
    BamChunk& cur = chunks[1];

    // chunk0 unchanged
    for (int i = 0; i < 3; ++i)
        ok &= check(pre.haps[i] == 1 && pre.phase_sets[i] == 1000, "cross chunk0 hap1");
    for (int i = 3; i < 6; ++i)
        ok &= check(pre.haps[i] == 2 && pre.phase_sets[i] == 1000, "cross chunk0 hap2");

    // chunk1: flipped and merged into chunk0's phase set
    for (int i = 0; i < 3; ++i)
        ok &= check(cur.haps[i] == 2 && cur.phase_sets[i] == 1000, "cross chunk1 former-hap1 → hap2");
    for (int i = 3; i < 6; ++i)
        ok &= check(cur.haps[i] == 1 && cur.phase_sets[i] == 1000, "cross chunk1 former-hap2 → hap1");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── test 5: no sidecar → no stitching ──────────────────────────────────────
static bool test_no_sidecar_no_stitch() {
    std::printf("--- test_no_sidecar_no_stitch ---\n");
    auto opts = make_opts();

    BamChunk chunk;
    chunk.region.tid = 0;
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 1000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 1000, {3, 4});
    for (int i = 0; i < 3; ++i) push_read(chunk, 1, 2000, {3, 4}); // would flip if sidecar present
    for (int i = 0; i < 3; ++i) push_read(chunk, 2, 2000, {1, 2});

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(chunk));
    stitch_chunk_haps(chunks, &opts, nullptr); // no sidecar

    BamChunk& c = chunks[0];
    bool ok = true;

    // Nothing should have changed
    for (int i = 0; i < 3; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 1000, "no-sidecar PS1000 hap1 untouched");
    for (int i = 6; i < 9; ++i)
        ok &= check(c.haps[i] == 1 && c.phase_sets[i] == 2000, "no-sidecar PS2000 hap1 untouched");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── test 6: within-chunk chain then cross-chunk ─────────────────────────────
//
// Chunk 0 has three phase blocks:
//   PS=1000  hap1→A hap2→B                  (correct)
//   PS=2000  hap1→B hap2→A  FLIPPED         (corrected by within-chunk sweep)
//   PS=3000  hap1→A hap2→B                  (correct, merged after sweep)
//
// After within-chunk stitching every read in chunk 0 carries phase_set=1000.
// The cross-chunk anchor therefore points at PS=1000 and draws thread signal
// from all 18 reads (not just the 6 that were originally in PS=3000, which is
// what max_pre_ps would have picked up had within-chunk stitching not run).
//
// Chunk 1 has one phase block:
//   PS=4000  hap1→B hap2→A  FLIPPED
//
// Expected: chunk 1 flipped and merged into PS=1000.
static bool test_chain_then_cross_chunk() {
    std::printf("--- test_chain_then_cross_chunk ---\n");
    auto sidecar = make_sidecar();
    auto opts    = make_opts();

    // ── chunk 0: three phase blocks ──────────────────────────────────────────
    BamChunk c0;
    c0.region.tid = 0;

    // PS=1000: 3× hap1 (A), 3× hap2 (B)   idx 0–5
    for (int i = 0; i < 3; ++i) push_read(c0, 1, 1000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(c0, 2, 1000, {3, 4});
    // PS=2000: FLIPPED — 3× hap1 (B), 3× hap2 (A)   idx 6–11
    for (int i = 0; i < 3; ++i) push_read(c0, 1, 2000, {3, 4});
    for (int i = 0; i < 3; ++i) push_read(c0, 2, 2000, {1, 2});
    // PS=3000: correct — 3× hap1 (A), 3× hap2 (B)   idx 12–17
    for (int i = 0; i < 3; ++i) push_read(c0, 1, 3000, {1, 2});
    for (int i = 0; i < 3; ++i) push_read(c0, 2, 3000, {3, 4});

    // ── chunk 1: one phase block, flipped ────────────────────────────────────
    BamChunk c1;
    c1.region.tid = 0;
    // PS=4000: FLIPPED — 3× hap1 (B), 3× hap2 (A)   idx 0–5
    for (int i = 0; i < 3; ++i) push_read(c1, 1, 4000, {3, 4});
    for (int i = 0; i < 3; ++i) push_read(c1, 2, 4000, {1, 2});

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(c0));
    chunks.push_back(std::move(c1));
    stitch_chunk_haps(chunks, &opts, &sidecar);

    bool ok = true;
    BamChunk& pre = chunks[0];
    BamChunk& cur = chunks[1];

    // ── verify chunk 0 ───────────────────────────────────────────────────────
    // PS=1000 reads (idx 0-5): unchanged
    for (int i = 0; i < 3; ++i)
        ok &= check(pre.haps[i] == 1 && pre.phase_sets[i] == 1000, "c0 PS1000 hap1 unchanged");
    for (int i = 3; i < 6; ++i)
        ok &= check(pre.haps[i] == 2 && pre.phase_sets[i] == 1000, "c0 PS1000 hap2 unchanged");
    // PS=2000 reads (idx 6-11): flipped and merged into PS=1000
    for (int i = 6; i < 9; ++i)
        ok &= check(pre.haps[i] == 2 && pre.phase_sets[i] == 1000, "c0 PS2000 former-hap1→hap2 ps1000");
    for (int i = 9; i < 12; ++i)
        ok &= check(pre.haps[i] == 1 && pre.phase_sets[i] == 1000, "c0 PS2000 former-hap2→hap1 ps1000");
    // PS=3000 reads (idx 12-17): nonflip, merged into PS=1000
    for (int i = 12; i < 15; ++i)
        ok &= check(pre.haps[i] == 1 && pre.phase_sets[i] == 1000, "c0 PS3000 hap1 nonflip ps1000");
    for (int i = 15; i < 18; ++i)
        ok &= check(pre.haps[i] == 2 && pre.phase_sets[i] == 1000, "c0 PS3000 hap2 nonflip ps1000");

    // ── verify chunk 1 ───────────────────────────────────────────────────────
    // Cross-chunk stitching uses max_pre_ps=1000 (all chunk0 reads after merge),
    // drawing signal from all 18 reads rather than only the original 6 PS=3000 reads.
    // chunk1 is flipped and merged into PS=1000.
    for (int i = 0; i < 3; ++i)
        ok &= check(cur.haps[i] == 2 && cur.phase_sets[i] == 1000, "c1 former-hap1→hap2 ps1000");
    for (int i = 3; i < 6; ++i)
        ok &= check(cur.haps[i] == 1 && cur.phase_sets[i] == 1000, "c1 former-hap2→hap1 ps1000");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── test 7: merged block carries all hap-thread evidence ───────────────────
//
// A-B can merge with support 2. B-C alone has no shared hap-unique threads,
// but merged AB-C has support 2 because C shares one hap1 and one hap2 thread
// with A. This catches regressions where the sequential sweep uses only the
// immediately previous original block instead of the live merged block state.
static bool test_within_chunk_merged_state_uses_all_evidence() {
    std::printf("--- test_within_chunk_merged_state_uses_all_evidence ---\n");
    PgbamSidecarData sidecar;
    for (uint32_t id = 10; id <= 17; ++id) sidecar.set_to_threads[id] = {id};
    auto opts = make_opts();

    BamChunk chunk;
    chunk.region.tid = 0;

    // PS=1000, block A: hap1 {10,11}, hap2 {14,15}
    push_read(chunk, 1, 1000, {10, 11});
    push_read(chunk, 2, 1000, {14, 15});

    // PS=2000, block B: hap1 {10,12}, hap2 {14,16}
    // A-B nonflip support is 2: thread 10 on hap1, thread 14 on hap2.
    push_read(chunk, 1, 2000, {10, 12});
    push_read(chunk, 2, 2000, {14, 16});

    // PS=3000, block C: hap1 {11,13}, hap2 {15,17}
    // B-C alone has support 0, but AB-C has support 2 through A's 11 and 15.
    push_read(chunk, 1, 3000, {11, 13});
    push_read(chunk, 2, 3000, {15, 17});

    std::vector<BamChunk> chunks;
    chunks.push_back(std::move(chunk));
    stitch_chunk_haps(chunks, &opts, &sidecar);

    BamChunk& c = chunks[0];
    bool ok = true;
    for (size_t i = 0; i < c.reads.size(); ++i)
        ok &= check(c.phase_sets[i] == 1000, "all blocks merged to PS1000");
    ok &= check(c.haps[0] == 1 && c.haps[1] == 2, "block A haps unchanged");
    ok &= check(c.haps[2] == 1 && c.haps[3] == 2, "block B haps unchanged");
    ok &= check(c.haps[4] == 1 && c.haps[5] == 2, "block C haps unchanged");

    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ─── main ────────────────────────────────────────────────────────────────────
int main() {
    int failures = 0;
    failures += test_within_chunk_flip()         ? 0 : 1;
    failures += test_within_chunk_nonflip()      ? 0 : 1;
    failures += test_within_chunk_chain()        ? 0 : 1;
    failures += test_cross_chunk_flip()          ? 0 : 1;
    failures += test_no_sidecar_no_stitch()      ? 0 : 1;
    failures += test_chain_then_cross_chunk()    ? 0 : 1;
    failures += test_within_chunk_merged_state_uses_all_evidence() ? 0 : 1;
    std::printf("\n%s  (%d failure%s)\n",
                failures == 0 ? "ALL PASS" : "SOME FAILURES",
                failures, failures == 1 ? "" : "s");
    return failures == 0 ? 0 : 1;
}
