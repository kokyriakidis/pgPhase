#ifndef PGPHASE_COLLECT_PHASE_NOISY_HPP
#define PGPHASE_COLLECT_PHASE_NOISY_HPP

/**
 * @file collect_phase_noisy.hpp
 * @brief Step 4: iterative noisy-region MSA variant calling and re-phasing.
 *
 * Mirrors longcallD `collect_var_main` step 4 from `collect_var.c`:
 *   sort_noisy_regs → collect_noisy_vars1 loop → assign_hap k-means (kCandGermlineVarCate).
 *
 * Alignment-heavy functions (WFA2-lib / abPOA) live in align.hpp / align.cpp.
 * This file covers the outer loop, MSA variant extraction, and merge into the chunk.
 */

#include "collect_types.hpp"
#include "align.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Step 4 top-level entry
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Iteratively call variants in noisy regions and re-run k-means phasing.
 *
 * Mirrors longcallD step 4 in `collect_var_main` (`collect_var.c`).
 * Called after `assign_hap_based_on_germline_het_vars_kmeans(kCandGermlineClean)`.
 *
 * Sorts noisy regions by `Interval::label` (longcallD `cr_label` / var_size) ascending,
 * then by length ascending (same bubble-sort as longcallD `sort_noisy_regs`).  For each
 * region in sorted order, calls `collect_noisy_vars1`; if any attempt returns ≥0 the
 * region is marked done.  When at least one region yields new variants (`ret > 0`),
 * re-runs k-means with `kCandGermlineVarCate` to incorporate noisy candidates.
 * Loops until no region makes further progress.
 *
 * @note No-ops when `chunk.noisy_regions` is empty.
 */
void collect_noisy_vars_step4(BamChunk& chunk, const Options& opts);

/**
 * @brief Process one noisy region: collect reads → extract reference → MSA →
 *        variant extraction → merge into chunk.
 *
 * Mirrors longcallD `collect_noisy_vars1` (`collect_var.c`).
 *
 * Return values (same semantics as longcallD):
 *   > 0  new variants found; caller should trigger k-means re-run.
 *   = 0  region was attempted but no new variants (too long, too deep, or MSA
 *        returned empty consensus); region is marked done by the outer loop.
 *   < 0  MSA could not resolve the region; region stays undone and may be retried.
 *
 * @param noisy_reg_i  Index into `chunk.noisy_regions`.
 */
int collect_noisy_vars1(BamChunk& chunk, const Options& opts, int noisy_reg_i);

// ════════════════════════════════════════════════════════════════════════════
// Sorting
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Return sorted noisy-region indices: label asc, then length asc.
 *
 * Mirrors longcallD `sort_noisy_regs` (bubble sort matching longcallD exactly).
 * `label` ↔ longcallD `cr_label` (var_size recorded during noisy-window detection).
 */
std::vector<int> sort_noisy_regs(const BamChunk& chunk);

// ════════════════════════════════════════════════════════════════════════════
// Read and reference collection
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Collect indices of non-skipped reads overlapping `[beg, end]`.
 *
 * Mirrors longcallD `collect_noisy_reg_reads1`.  Iterates in
 * `chunk.ordered_read_ids` order (same as longcallD `chunk->ordered_read_ids`).
 * Overlap test: `r.beg > end || r.end <= beg` → skip (longcallD `beg/end` from
 * `chunk->digars[read_i].beg/end`; pgPhase equivalent: `ReadRecord::beg/end`).
 */
std::vector<int> collect_noisy_reg_reads(const BamChunk& chunk,
                                         hts_pos_t beg, hts_pos_t end);

/**
 * @brief Extract the reference slice `[beg, end]` as 2-bit encoded bytes.
 *
 * Mirrors longcallD `collect_reg_ref_bseq` (`seq.c`).  Encoding follows
 * longcallD `nst_nt4_table`: A=0, C=1, G=2, T/U=3, other=4.  `beg` and `end`
 * are clipped to `[chunk.ref_beg, chunk.ref_end]` in-place (same as longcallD).
 *
 * @param beg in/out — clipped to chunk.ref_beg if smaller.
 * @param end in/out — clipped to chunk.ref_end if larger.
 * @return Vector of 2-bit bases; empty when the clipped interval is degenerate
 *         or `chunk.ref_seq` is empty.
 */
std::vector<uint8_t> collect_reg_ref_bseq(const BamChunk& chunk,
                                           hts_pos_t& beg, hts_pos_t& end);

// ════════════════════════════════════════════════════════════════════════════
// Variant extraction from MSA
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Extract `NOISY_CAND_HET` / `NOISY_CAND_HOM` candidates from MSA strings.
 *
 * Mirrors longcallD `make_vars_from_msa_cons_aln` (`collect_var.c`).  Walks
 * consensus-vs-reference and read-vs-consensus alignment strings to locate
 * positions where the consensus differs from the reference, producing new
 * `CandidateVariant` objects and per-read allele profiles for the noisy region.
 *
 * @param noisy_vars     [out] Newly discovered variants.
 * @param noisy_var_cate [out] Category for each new variant.
 * @param noisy_rvp      [out] Per-read allele profiles covering the new variants.
 * @return Number of new variants produced.
 */
int make_vars_from_msa_cons_aln(
    const Options& opts, BamChunk& chunk,
    int n_noisy_reads, const std::vector<int>& read_ids,
    hts_pos_t noisy_reg_beg,
    int n_cons,
    const std::array<int, 2>& clu_n_seqs,
    const std::array<std::vector<int>, 2>& clu_read_ids,
    const std::array<std::vector<AlnStr>, 2>& aln_strs,
    std::vector<CandidateVariant>& noisy_vars,
    std::vector<VariantCategory>& noisy_var_cate,
    std::vector<ReadVariantProfile>& noisy_rvp);

// ════════════════════════════════════════════════════════════════════════════
// Merging into chunk
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Merge newly discovered noisy variants into `chunk.candidates` and
 *        `chunk.read_var_profile`, then rebuild `chunk.read_var_cr`.
 *
 * Mirrors longcallD `merge_var_profile` (`collect_var.c`).  Inserts new
 * `CandidateVariant` objects into the sorted candidate table, re-indexes all
 * existing `ReadVariantProfile` entries to account for shifted variant indices,
 * merges `noisy_rvp` allele observations into the affected reads' profiles, and
 * rebuilds the `chunk.read_var_cr` interval tree so that the subsequent
 * k-means call (`kCandGermlineVarCate`) sees the complete, merged profile.
 */
void merge_var_profile(BamChunk& chunk,
                       const std::vector<CandidateVariant>& noisy_vars,
                       const std::vector<VariantCategory>& noisy_var_cate,
                       const std::vector<ReadVariantProfile>& noisy_rvp);

} // namespace pgphase_collect

#endif
