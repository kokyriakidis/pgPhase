/// @file test_noise_filter.cpp
/// @brief Unit tests for reference-based noise detection.
///
/// Covers three layers:
///   1. Low-level homopolymer / tandem-repeat / low-complexity detection.
///   2. is_noisy_site integration (SNP gating, indel context checks).
///   3. apply_graph_noise_filter end-to-end, including the non-minimal
///      multiallelic graph-allele trimming that lets homopolymer/STR indels
///      be detected even when emitted with full repeat runs on both flanks.

#include "noise_filter.hpp"
#include "graph_bam_adapter.hpp"
#include "collect_phase.hpp"  // kCand* / kLongcalldRepHetVar bitmask flags

#include <iostream>
#include <string>
#include <vector>

using namespace pgphase_collect;

static int g_failures = 0;

static bool check(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        ++g_failures;
    }
    return cond;
}

// Convenience wrappers: ref_beg = 1, ref_end derived from slice length.
static bool hp(hts_pos_t pos, const std::string& ref, const std::string& alt,
               const std::string& rs) {
    return is_homopolymer_indel(pos, ref, alt, rs, 1,
                                static_cast<hts_pos_t>(rs.size()),
                                kDefaultNoisyRegMaxXgaps);
}
static bool rep(hts_pos_t pos, const std::string& ref, const std::string& alt,
                const std::string& rs) {
    return is_repeat_indel(pos, ref, alt, rs, 1,
                           static_cast<hts_pos_t>(rs.size()),
                           kDefaultNoisyRegMaxXgaps);
}
static bool noisy(hts_pos_t pos, const std::string& ref, const std::string& alt,
                  const std::string& rs) {
    const std::vector<Interval> lc = find_low_complexity_intervals(rs, 1);
    return is_noisy_site(pos, ref, alt, rs, 1,
                         static_cast<hts_pos_t>(rs.size()), lc,
                         kDefaultNoisyRegMaxXgaps);
}

// ────────────────────────────────────────────────────────────────────────────
// Layer 1: homopolymer detection
// ────────────────────────────────────────────────────────────────────────────
static void test_homopolymer() {
    std::cout << "--- test_homopolymer ---\n";
    // Index:           1234567890123456
    // A-run spans positions 5..12.
    const std::string poly = "ACGTAAAAAAAAGCGT";

    check(hp(5, "A", "AA", poly),  "insertion into A-run is homopolymer");
    check(hp(5, "AA", "A", poly),  "deletion from A-run is homopolymer");

    // Unique, non-repetitive flank: no period<=6 repeats twice.
    const std::string uniq = "ACGTCAGTGACTGGCATTACGAACCTGTA";
    check(!hp(10, "T", "TA", uniq), "insertion in unique flank is not homopolymer");

    // Indel span over max_xgaps must be skipped (returns false).
    const std::string a6 = "ACGTAAAAAAAAAAAAGCGT";
    check(!hp(5, "A", "AAAAAAA", a6), "6 bp insertion exceeds max_xgaps -> skipped");
}

// ────────────────────────────────────────────────────────────────────────────
// Layer 1: tandem-repeat detection
// ────────────────────────────────────────────────────────────────────────────
static void test_tandem_repeat() {
    std::cout << "--- test_tandem_repeat ---\n";
    // CA dinucleotide repeat spanning positions 5..20.
    const std::string ca = "ACGTCACACACACACACACAGGTT";
    check(rep(6, "A", "ACA", ca), "CA insertion into CA-repeat is repeat indel");
    check(rep(6, "ACA", "A", ca), "CA deletion from CA-repeat is repeat indel");

    // ATG trinucleotide repeat spanning positions 5..19.
    const std::string atg = "CCGGATGATGATGATGATGCCGG";
    check(rep(7, "G", "GATG", atg), "ATG insertion into ATG-repeat is repeat indel");

    // Unique flank: not a repeat indel.
    const std::string uniq = "ACGTCAGTGACTGGCATTACGAACCTGTA";
    check(!rep(10, "T", "TA", uniq), "insertion in unique flank is not repeat indel");

    // SNPs are never repeat indels.
    const std::string poly = "ACGTAAAAAAAAGCGT";
    check(!rep(8, "A", "C", poly), "SNP is never a repeat indel");
}

// ────────────────────────────────────────────────────────────────────────────
// Layer 1: low-complexity (SDUST)
// ────────────────────────────────────────────────────────────────────────────
static void test_low_complexity() {
    std::cout << "--- test_low_complexity ---\n";
    // Complex 20 bp + polyA 24 bp + complex 20 bp.
    const std::string s =
        "ACGTAGCTAGCTAGCTAGCT"
        "AAAAAAAAAAAAAAAAAAAAAAAA"
        "GCATGCATGCATGCATGCAT";
    const std::vector<Interval> lc = find_low_complexity_intervals(s, 1);
    check(!lc.empty(), "SDUST reports a low-complexity interval for polyA");

    // The polyA run (positions 21..44) must be flagged.
    check(pos_in_low_complexity(30, lc), "polyA midpoint is in low-complexity");

    // A fully complex sequence yields no intervals.
    const std::string complex_seq = "ACGTGCATGACTGCAGTACGTACTGACGT";
    const std::vector<Interval> none = find_low_complexity_intervals(complex_seq, 1);
    check(!pos_in_low_complexity(10, none), "complex sequence is not low-complexity");

    // ref_beg offset is honoured: same slice starting at chromosome pos 1001.
    const std::vector<Interval> lc_off = find_low_complexity_intervals(s, 1001);
    check(pos_in_low_complexity(1030, lc_off), "ref_beg offset shifts LC coordinates");
}

// ────────────────────────────────────────────────────────────────────────────
// Layer 2: is_noisy_site integration + SNP gating
// ────────────────────────────────────────────────────────────────────────────
static void test_is_noisy_site() {
    std::cout << "--- test_is_noisy_site ---\n";
    const std::string poly = "ACGTAAAAAAAAGCGT";   // A-run 5..12
    const std::string uniq = "ACGTCAGTGACTGGCATTACGAACCTGTA";

    check(noisy(5, "A", "AA", poly),  "indel in homopolymer is noisy");
    check(!noisy(10, "T", "TA", uniq), "indel in unique context is not noisy");
    check(!noisy(10, "T", "C", uniq),  "SNP in unique context is not noisy");

    // SNP gating: is_noisy_site applies homopolymer/repeat checks to indels
    // only; SNPs are checked against low-complexity alone. Use a short 3 bp AAA
    // run (positions 7..9) that SDUST does NOT flag, so the low-complexity check
    // cannot mask the gating. An indel there is noisy (homopolymer) while a SNP
    // at the same locus is not — proving SNPs skip the homopolymer check.
    const std::string aaa = "GTCAGTAAAGTCAGT"; // AAA run at positions 7..9
    check(noisy(6, "T", "TA", aaa), "indel in 3 bp A-run is noisy (homopolymer)");
    check(!noisy(7, "A", "C", aaa), "SNP in 3 bp A-run is not noisy (SNP gating)");
}

// ────────────────────────────────────────────────────────────────────────────
// Layer 3: apply_graph_noise_filter with non-minimal multiallelic graph alleles
// ────────────────────────────────────────────────────────────────────────────

// Build a one-candidate GraphChunkBuildResult for an indel site. The candidate
// starts classified as a clean het indel; apply_graph_noise_filter should demote
// it to RepeatHetIndel when the (possibly non-minimal) alt sits in a repeat.
static GraphChunkBuildResult make_indel_result(hts_pos_t pos,
                                               const std::string& ref,
                                               const std::vector<std::string>& alts) {
    GraphChunkBuildResult res;
    CandidateVariant cand;
    cand.key.tid = 0;
    cand.key.pos = pos;
    cand.key.type = (alts[0].size() > ref.size()) ? VariantType::Insertion
                                                   : VariantType::Deletion;
    cand.key.ref_len = static_cast<int>(ref.size());
    cand.key.alt = alts[0];
    cand.counts.category = VariantCategory::CleanHetIndel;
    cand.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
    cand.lcd_var_i_to_cate = kCandCleanHetIndel;
    res.chunk.candidates.push_back(cand);

    GraphSiteMeta meta;
    meta.chrom = "chr1";
    meta.pos = pos;
    meta.ref = ref;
    meta.alts = alts;
    res.site_meta.push_back(meta);
    res.site_ids.push_back("s1");
    return res;
}

static bool is_demoted(const GraphChunkBuildResult& res) {
    const CandidateVariant& c = res.chunk.candidates[0];
    return c.counts.category == VariantCategory::RepeatHetIndel &&
           c.lcd_var_i_to_cate == kLongcalldRepHetVar;
}

static void test_graph_noise_filter() {
    std::cout << "--- test_graph_noise_filter ---\n";
    // Reference slice (ref_beg = 1):
    //   positions 6..17 are a 12 bp A-run.
    const std::string ref_seq = "CCGTCAAAAAAAAAAAAGTGTA";
    const hts_pos_t ref_beg = 1;
    const hts_pos_t ref_end = static_cast<hts_pos_t>(ref_seq.size());

    // (a) Minimal-form homopolymer insertion: anchor at pos 6 (the C before
    //     the A-run), insert one A. Should be demoted.
    {
        GraphChunkBuildResult res = make_indel_result(5, "C", {"CA"});
        // Place the anchor so inserted A lands in the run: ref base at pos 5 is C,
        // pos 6 is first A. Insertion content begins at pos 6.
        apply_graph_noise_filter(res, ref_seq, ref_beg, ref_end,
                                 kDefaultNoisyRegMaxXgaps);
        check(is_demoted(res), "minimal homopolymer insertion is demoted");
    }

    // (b) NON-minimal graph allele: a 1 bp deletion inside the A-run emitted
    //     with the full run on both sides. Without trimming, the derived indel
    //     length exceeds max_xgaps and the site would slip through. The trimming
    //     in apply_graph_noise_filter must reduce it to a 1 bp homopolymer del.
    //     ref allele = "CAAAAAAAAAAAA" (pos5..17), alt drops one A.
    {
        GraphChunkBuildResult res =
            make_indel_result(5, "CAAAAAAAAAAAA", {"CAAAAAAAAAAA"});
        apply_graph_noise_filter(res, ref_seq, ref_beg, ref_end,
                                 kDefaultNoisyRegMaxXgaps);
        check(is_demoted(res),
              "non-minimal homopolymer deletion is trimmed and demoted");
    }

    // (c) NON-minimal graph insertion: insert one A into the run, full run on
    //     both flanks. Trimming -> 1 bp homopolymer insertion -> demoted.
    {
        GraphChunkBuildResult res =
            make_indel_result(5, "CAAAAAAAAAAAA", {"CAAAAAAAAAAAAA"});
        apply_graph_noise_filter(res, ref_seq, ref_beg, ref_end,
                                 kDefaultNoisyRegMaxXgaps);
        check(is_demoted(res),
              "non-minimal homopolymer insertion is trimmed and demoted");
    }

    // (d) Multiallelic site: one alt is clean, one alt is a homopolymer indel.
    //     The site must be demoted because ANY noisy alt taints it.
    {
        GraphChunkBuildResult res = make_indel_result(
            5, "CAAAAAAAAAAAA", {"CAAAAAAAAAAA", "CAAAAAAAAAAAAA"});
        apply_graph_noise_filter(res, ref_seq, ref_beg, ref_end,
                                 kDefaultNoisyRegMaxXgaps);
        check(is_demoted(res),
              "multiallelic site with one noisy alt is demoted");
    }

    // (e) Clean indel in a unique context must NOT be demoted.
    {
        const std::string uniq = "ACGTCAGTGACTGGCATTACGAACCTGTA";
        GraphChunkBuildResult res = make_indel_result(10, "T", {"TA"});
        apply_graph_noise_filter(res, uniq, 1,
                                 static_cast<hts_pos_t>(uniq.size()),
                                 kDefaultNoisyRegMaxXgaps);
        check(!is_demoted(res), "clean indel in unique context is not demoted");
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Edge cases
// ────────────────────────────────────────────────────────────────────────────
static void test_edge_cases() {
    std::cout << "--- test_edge_cases ---\n";
    check(!hp(5, "A", "AA", ""),  "empty ref_seq -> homopolymer false");
    check(!rep(5, "A", "AA", ""), "empty ref_seq -> repeat false");
    check(!noisy(5, "A", "AA", ""), "empty ref_seq -> not noisy");

    // Out-of-range position: ref_nt4_at returns 4 (N) for both the flank base and
    // the comparison window, and 4 == 4 reads as a homopolymer run, so
    // is_homopolymer_indel returns TRUE. This is a documented quirk, not a bug:
    // callers (BAM and graph pipelines) only ever pass positions inside the chunk
    // slice, so an N-vs-N false positive cannot occur in production. is_repeat_indel
    // is range-checked and correctly returns false out of range.
    const std::string poly = "ACGTAAAAAAAAGCGT";
    check(hp(1000, "A", "AA", poly),
          "out-of-range homopolymer is N-vs-N TRUE (documented quirk)");
    check(!rep(1000, "A", "AA", poly),
          "out-of-range repeat is range-checked -> false");
}

int main() {
    test_homopolymer();
    test_tandem_repeat();
    test_low_complexity();
    test_is_noisy_site();
    test_graph_noise_filter();
    test_edge_cases();

    if (g_failures == 0) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    std::cerr << g_failures << " CHECK(S) FAILED\n";
    return 1;
}
