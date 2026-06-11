/// @file test_hybrid_inject.cpp
/// @brief Unit tests for the hybrid graph-candidate quality gate.

#include "hybrid_inject.hpp"

#include "collect_phase.hpp"
#include "collect_var.hpp"
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

    // ── vcf_to_variant_key: deletion / insertion normalization ───────────────
    // Deletions must strip the full shared prefix to match the BAM convention
    // (variant_key_from_digar): alt = "", ref_len = deleted span, pos = first
    // deleted base.  A homopolymer deletion like TAA->TA is a 1 bp deletion and
    // must NOT keep a residual alt base or an inflated ref_len, otherwise it
    // collides with a genuine BAM deletion at the same locus (exact_comp_var_site
    // ignores alt for deletions) and overwrites it during merge.
    {
        // SNP: G->T at pos 100.
        VariantKey snp = vcf_to_variant_key(0, 100, "G", "T");
        ok &= check(snp.type == VariantType::Snp && snp.pos == 100 &&
                        snp.ref_len == 1 && snp.alt == "T",
                    "SNP G->T maps to pos=100 ref_len=1 alt=T");

        // Single-base-anchor deletion: TA->T at pos 100 is a 1 bp deletion.
        VariantKey del1 = vcf_to_variant_key(0, 100, "TA", "T");
        ok &= check(del1.type == VariantType::Deletion && del1.pos == 101 &&
                        del1.ref_len == 1 && del1.alt.empty(),
                    "DEL TA->T maps to pos=101 ref_len=1 alt=\"\"");

        // Homopolymer 1 bp deletion: TAA->TA.  Old code yielded ref_len=2,
        // alt=\"A\"; normalized form is pos=102 ref_len=1 alt=\"\".
        VariantKey del_hp1 = vcf_to_variant_key(0, 100, "TAA", "TA");
        ok &= check(del_hp1.type == VariantType::Deletion && del_hp1.pos == 102 &&
                        del_hp1.ref_len == 1 && del_hp1.alt.empty(),
                    "DEL TAA->TA normalizes to pos=102 ref_len=1 alt=\"\"");

        // Homopolymer 2 bp deletion: TAA->T.  pos=101 ref_len=2 alt=\"\".
        VariantKey del_hp2 = vcf_to_variant_key(0, 100, "TAA", "T");
        ok &= check(del_hp2.type == VariantType::Deletion && del_hp2.pos == 101 &&
                        del_hp2.ref_len == 2 && del_hp2.alt.empty(),
                    "DEL TAA->T normalizes to pos=101 ref_len=2 alt=\"\"");

        // The 1 bp and 2 bp homopolymer deletions must encode distinct keys so
        // they do not alias under exact_comp_var_site (which keys deletions on
        // pos and ref_len).
        ok &= check(del_hp1.pos != del_hp2.pos || del_hp1.ref_len != del_hp2.ref_len,
                    "1bp and 2bp homopolymer deletions are distinct keys");

        // Single-base-anchor insertion: T->TA at pos 100 -> ref_len=0 alt=A at
        // pos 101 (unchanged from the pre-normalization behavior).
        VariantKey ins = vcf_to_variant_key(0, 100, "T", "TA");
        ok &= check(ins.type == VariantType::Insertion && ins.pos == 101 &&
                        ins.ref_len == 0 && ins.alt == "A",
                    "INS T->TA maps to pos=101 ref_len=0 alt=A");

        // Multi-base-anchor insertion: TA->TAAA at pos 100 is a 2 bp insertion
        // after the shared "TA" run.  Full-prefix strip yields pos=102,
        // ref_len=0, alt="AA"; the old single-base strip mis-encoded this as
        // pos=101 ref_len=1 alt="AAA".
        VariantKey ins_hp = vcf_to_variant_key(0, 100, "TA", "TAAA");
        ok &= check(ins_hp.type == VariantType::Insertion && ins_hp.pos == 102 &&
                        ins_hp.ref_len == 0 && ins_hp.alt == "AA",
                    "INS TA->TAAA normalizes to pos=102 ref_len=0 alt=AA");
    }

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
    //                       17 ref / 3 alt (AF 0.15) -> LowCoverage
    //                       (LowAlleleFraction is folded to LowCoverage so it
    //                        is pruned from output, matching the BAM pipeline)
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

    ok &= check(chunk.candidates[4].counts.category == VariantCategory::LowCoverage,
                "low AF folded to LowCoverage (so it is pruned from output)");
    ok &= check((chunk.candidates[4].lcd_var_i_to_cate & kCandCleanHetSnp) == 0,
                "low AF excluded from het mask");

    // Empty input is a no-op.
    std::unordered_set<int> empty_set;
    ok &= check(classify_graph_only_candidates(chunk, empty_set, opts) == 0,
                "empty candidate set promotes nothing");

    // ── apply_hybrid_noise_filter: indel low-complexity screening ────────────
    // Reference with a 40 bp homopolymer A-run (1-based 121..160) flanked by
    // mixed sequence.  SDUST flags the run; positions ~150 ARE low-complexity.
    // SNPs are intentionally NOT demoted: both the BAM pipeline (NOISY_CAND_HET
    // recall) and the standalone graph pipeline (CLEAN_HET_SNP) keep low-
    // complexity het SNPs as real calls, so the hybrid pipeline must too.  Only
    // indels are screened (homopolymer / repeat / low-complexity).
    const std::string mixed =
        "ACGTTGCAGATCCTGAGTACGTCAGTTGACCATGGATCAGTACTGGCATGACTTAGCATGC"
        "TGACAGTCATGCATGACTGCATGCTAGCATCGATCGATTGCATGCATGCTAGCTGATCAGT";
    PhasingChunk nchunk;
    nchunk.ref_beg = 1;
    nchunk.ref_end = static_cast<hts_pos_t>(mixed.size() * 2 + 40);
    nchunk.ref_seq = mixed + std::string(40, 'A') + mixed;

    // 0: het SNP in low-complexity run (pos 130)  -> kept (SNPs not demoted)
    // 1: het indel in low-complexity run (pos 150)-> demoted to RepeatHetIndel
    CandidateVariant snp_in_lc = make_graph_snp(130, 5, 5);
    snp_in_lc.counts.category = VariantCategory::CleanHetSnp;
    snp_in_lc.counts.candvarcate_initial = VariantCategory::CleanHetSnp;
    snp_in_lc.lcd_var_i_to_cate = kCandCleanHetSnp;

    CandidateVariant ins_noisy = make_graph_snp(150, 5, 5);
    ins_noisy.key.type = VariantType::Insertion;
    ins_noisy.key.ref_len = 0;
    ins_noisy.key.alt = "A";
    ins_noisy.counts.category = VariantCategory::CleanHetIndel;
    ins_noisy.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
    ins_noisy.lcd_var_i_to_cate = kCandCleanHetIndel;

    nchunk.candidates.push_back(snp_in_lc);
    nchunk.candidates.push_back(ins_noisy);

    std::unordered_set<int> noise_set = {0, 1};
    apply_hybrid_noise_filter(nchunk, nchunk.ref_seq, nchunk.ref_beg,
                              nchunk.ref_end, noise_set, opts.noisy_reg_max_xgaps);

    ok &= check(nchunk.candidates[0].counts.category == VariantCategory::CleanHetSnp,
                "het SNP in low-complexity kept (not demoted)");
    ok &= check(nchunk.candidates[0].lcd_var_i_to_cate == kCandCleanHetSnp,
                "het SNP in low-complexity keeps CleanHetSnp flag");

    ok &= check(nchunk.candidates[1].counts.category == VariantCategory::RepeatHetIndel,
                "noisy indel demoted to RepeatHetIndel");
    ok &= check(nchunk.candidates[1].lcd_var_i_to_cate == kLongcalldRepHetVar,
                "noisy indel gets RepHetVar flag");

    // ── apply_hybrid_noise_filter: original VCF strings override key-based ────
    // When a GraphOnlyVcfAlleles map supplies the catalog (pos, ref, alt), the
    // filter screens on those instead of reconstructing from the VariantKey.
    // This mirrors the standalone graph pipeline (apply_graph_noise_filter),
    // which screens on catalog allele strings, and prevents over-demoting graph
    // het indels whose key-reconstructed form lands in a different context.
    {
        // Put the candidate's KEY position inside the homopolymer run (would be
        // demoted via reconstruction), but supply original VCF coordinates in a
        // unique, non-low-complexity context so is_noisy_site returns false.
        PhasingChunk ochunk;
        ochunk.ref_beg = 1;
        ochunk.ref_end = static_cast<hts_pos_t>(mixed.size() * 2 + 40);
        ochunk.ref_seq = mixed + std::string(40, 'A') + mixed;

        CandidateVariant ins = make_graph_snp(150, 5, 5);  // key pos in A-run
        ins.key.type = VariantType::Insertion;
        ins.key.ref_len = 0;
        ins.key.alt = "A";
        ins.counts.category = VariantCategory::CleanHetIndel;
        ins.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
        ins.lcd_var_i_to_cate = kCandCleanHetIndel;
        ochunk.candidates.push_back(ins);

        // Original VCF: a clean 1 bp insertion at pos 10 (mixed sequence,
        // not low-complexity, not homopolymer).
        GraphOnlyVcfAlleles alleles;
        alleles[0] = GraphOnlyVcfAllele{10, "A", "AC"};

        std::unordered_set<int> oset = {0};
        apply_hybrid_noise_filter(ochunk, ochunk.ref_seq, ochunk.ref_beg,
                                  ochunk.ref_end, oset,
                                  opts.noisy_reg_max_xgaps, &alleles);

        ok &= check(ochunk.candidates[0].counts.category ==
                        VariantCategory::CleanHetIndel,
                    "indel kept when original VCF strings are non-noisy "
                    "(even though key pos is in a homopolymer)");
    }

    // ── prune_not_candidate_variants: drop non-call categories ───────────────
    // Mirrors the BAM pipeline's end-of-classification prune.  The hybrid
    // pipeline re-runs this after appending graph-only candidates so failed
    // gates (LowCoverage / NonVariant / StrandBias) do not leak into output.
    PhasingChunk pchunk;
    CandidateVariant keep_snp = make_graph_snp(10, 5, 5);
    keep_snp.counts.category = VariantCategory::CleanHetSnp;
    CandidateVariant drop_lowcov = make_graph_snp(20, 1, 1);
    drop_lowcov.counts.category = VariantCategory::LowCoverage;
    CandidateVariant drop_nonvar = make_graph_snp(30, 5, 5);
    drop_nonvar.counts.category = VariantCategory::NonVariant;
    CandidateVariant keep_hom = make_graph_snp(40, 0, 10);
    keep_hom.counts.category = VariantCategory::CleanHom;
    pchunk.candidates.push_back(keep_snp);
    pchunk.candidates.push_back(drop_lowcov);
    pchunk.candidates.push_back(drop_nonvar);
    pchunk.candidates.push_back(keep_hom);

    prune_not_candidate_variants(pchunk);

    ok &= check(pchunk.candidates.size() == 2,
                "prune keeps only the two real-call candidates");
    bool has_lowcov_or_nonvar = false;
    for (const CandidateVariant& c : pchunk.candidates) {
        if (c.counts.category == VariantCategory::LowCoverage ||
            c.counts.category == VariantCategory::NonVariant ||
            c.counts.category == VariantCategory::StrandBias)
            has_lowcov_or_nonvar = true;
    }
    ok &= check(!has_lowcov_or_nonvar,
                "prune removes LowCoverage / NonVariant / StrandBias");

    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
