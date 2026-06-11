/// @file test_hybrid_inject.cpp
/// @brief Unit tests for the hybrid graph-candidate quality gate.

#include "hybrid_inject.hpp"

#include "collect_phase.hpp"
#include "phasing_types.hpp"

#include <iostream>
#include <string>
#include <unordered_set>

using namespace pgphase_collect;

static bool check(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

/// Build a graph-only SNP candidate with the given ref/alt coverage, mirroring
/// the unclassified state add_graph_only_candidate leaves behind (flag 0).
static CandidateVariant make_graph_snp(hts_pos_t pos, int ref_cov, int alt_cov) {
    CandidateVariant cand;
    cand.key.tid = 0;
    cand.key.pos = pos;
    cand.key.type = VariantType::Snp;
    cand.key.ref_len = 1;
    cand.key.alt = "A";
    cand.counts.n_uniq_alles = 2;
    cand.counts.ref_cov = ref_cov;
    cand.counts.alt_cov = alt_cov;
    cand.counts.total_cov = ref_cov + alt_cov;
    cand.counts.category = VariantCategory::LowCoverage;
    cand.counts.candvarcate_initial = VariantCategory::LowCoverage;
    cand.lcd_var_i_to_cate = 0;
    cand.lcd_make_variants_region_pass = true;
    return cand;
}

int main() {
    bool ok = true;

    // Chunk with a reference slice long enough for the candidate positions.
    PhasingChunk chunk;
    chunk.ref_beg = 1;
    chunk.ref_end = 200;
    chunk.ref_seq = std::string(200, 'C');

    // Defaults: min_depth=5, min_alt_depth=2, min_af=0.20, max_af=0.80.
    Options opts;

    // 0: clean het      — 5 ref / 5 alt  (AF 0.50)  -> promoted, CleanHetSnp
    // 1: homozygous      — 0 ref / 10 alt (AF 1.00)  -> CleanHom, not het
    // 2: low total depth — 2 ref / 1 alt  (depth 3)  -> LowCoverage
    // 3: low alt depth   — 9 ref / 1 alt  (alt 1)    -> LowCoverage
    // 4: low AF          — 9 ref / 1 alt? use depth ok but AF<0.20:
    //                       17 ref / 3 alt (AF 0.15) -> LowAlleleFraction
    chunk.candidates.push_back(make_graph_snp(10, 5, 5));
    chunk.candidates.push_back(make_graph_snp(20, 0, 10));
    chunk.candidates.push_back(make_graph_snp(30, 2, 1));
    chunk.candidates.push_back(make_graph_snp(40, 9, 1));
    chunk.candidates.push_back(make_graph_snp(50, 17, 3));

    std::unordered_set<int> graph_only = {0, 1, 2, 3, 4};

    const int promoted = classify_graph_only_candidates(chunk, graph_only, opts);

    ok &= check(promoted == 1, "exactly one candidate promoted to clean het");

    ok &= check(chunk.candidates[0].lcd_var_i_to_cate == kCandCleanHetSnp,
                "clean het SNP gets CleanHetSnp flag");
    ok &= check(chunk.candidates[0].counts.category == VariantCategory::CleanHetSnp,
                "clean het SNP category is CleanHetSnp");

    ok &= check(chunk.candidates[1].counts.category == VariantCategory::CleanHom,
                "homozygous site classified CleanHom");
    ok &= check((chunk.candidates[1].lcd_var_i_to_cate & kCandCleanHetSnp) == 0,
                "homozygous site excluded from het mask");

    ok &= check(chunk.candidates[2].counts.category == VariantCategory::LowCoverage,
                "low total depth classified LowCoverage");
    ok &= check((chunk.candidates[2].lcd_var_i_to_cate & kCandCleanHetSnp) == 0,
                "low total depth excluded from het mask");

    ok &= check(chunk.candidates[3].counts.category == VariantCategory::LowCoverage,
                "low alt depth classified LowCoverage");
    ok &= check((chunk.candidates[3].lcd_var_i_to_cate & kCandCleanHetSnp) == 0,
                "low alt depth excluded from het mask");

    ok &= check(chunk.candidates[4].counts.category == VariantCategory::LowAlleleFraction,
                "low AF classified LowAlleleFraction");
    ok &= check((chunk.candidates[4].lcd_var_i_to_cate & kCandCleanHetSnp) == 0,
                "low AF excluded from het mask");

    // Empty input is a no-op.
    std::unordered_set<int> empty_set;
    ok &= check(classify_graph_only_candidates(chunk, empty_set, opts) == 0,
                "empty candidate set promotes nothing");

    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
