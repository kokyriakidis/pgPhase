#ifndef PGPHASE_COLLECT_PHASE_HPP
#define PGPHASE_COLLECT_PHASE_HPP

/**
 * @file collect_phase.hpp
 * @brief k-means read-haplotype clustering declarations for collect-bam-variation.
 *
 * Implements the germline phasing scaffold from longcallD's `assign_hap.c`:
 * iterative read→haplotype assignment driven by per-variant consensus allele
 * profiles (`hap_to_cons_alle`), up to 10 refinement rounds.
 */

#include "collect_types.hpp"

#include <cstdint>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Candidate-category bitmask flags
// ════════════════════════════════════════════════════════════════════════════

/**
 * @name Category flags
 * @brief Bitmask constants for `assign_hap_based_on_germline_het_vars_kmeans`.
 *
 * Each flag selects which `VariantCategory` values participate in phasing.
 * Combine with `|` to build the `flags` argument (longcallD `target_var_cate`). Mirrors
 * `LONGCALLD_*` masks in longcallD `collect_var.h`.
 */
///@{
constexpr uint32_t kCandCleanHetSnp   = 1u << 0; ///< `LONGCALLD_CLEAN_HET_SNP`
constexpr uint32_t kCandCleanHetIndel = 1u << 1; ///< `LONGCALLD_CLEAN_HET_INDEL`
constexpr uint32_t kCandCleanHom      = 1u << 2; ///< `LONGCALLD_CLEAN_HOM_VAR`
constexpr uint32_t kCandNoisyCandHet  = 1u << 3; ///< `LONGCALLD_NOISY_CAND_HET_VAR`
constexpr uint32_t kCandNoisyCandHom  = 1u << 4; ///< `LONGCALLD_NOISY_CAND_HOM_VAR`

/** `LONGCALLD_CAND_HET_VAR_CATE`: phased het sites (clean + noisy het), no hom. */
constexpr uint32_t kCandHetVarCate =
    kCandCleanHetSnp | kCandCleanHetIndel | kCandNoisyCandHet;

/** `LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE` — sole k-means call in longcallD `collect_var_main` phase 3. */
constexpr uint32_t kCandGermlineClean = kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom;

/** `LONGCALLD_CAND_GERMLINE_VAR_CATE` (used after longcallD noisy MSA step 4, not in phase 3 block). */
constexpr uint32_t kCandGermlineVarCate =
    kCandGermlineClean | kCandNoisyCandHet | kCandNoisyCandHom;
///@}

// ════════════════════════════════════════════════════════════════════════════
// Public interface
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Maps a `VariantCategory` to its bitmask flag (0 for unmapped categories).
 */
uint32_t category_to_flag(VariantCategory c);

/**
 * @brief Stitch haplotype assignments across chunk boundaries.
 *
 * Mirrors longcallD `stitch_var_main` / `flip_variant_hap` from `collect_var.c`.
 *
 * For each pair of adjacent chunks on the same contig, reads that span the boundary
 * (already catalogued in `down_ovlp_read_i` / `up_ovlp_read_i`) are matched by query
 * name. If the majority have inconsistent haplotype labels (hap1 in one chunk → hap2
 * in the other), the current chunk's haplotype assignments are flipped. In either case,
 * if both chunks have a valid phase block touching the boundary, the current chunk's
 * earliest phase-set anchor is rewritten to the previous chunk's latest anchor, merging
 * the two blocks into one.
 *
 * @param chunks Ordered, adjacent `BamChunk` objects from `collect_chunk_batch_parallel`,
 *               each already phased by `assign_hap_based_on_germline_het_vars_kmeans`.
 *               Chunks on different contigs are silently skipped at the boundary.
 */
void stitch_chunk_haps(std::vector<BamChunk>& chunks,
                       const Options* opts = nullptr,
                       const PgbamSidecarData* pgbam_sidecar = nullptr);

/**
 * @brief Assign haplotypes and phase sets to reads via iterative k-means clustering.
 *
 * Mirrors longcallD `assign_hap_based_on_germline_het_vars_kmeans` from `assign_hap.c`.
 *
 * @note Parity is with **assign_hap.c only**: category masks use pgPhase `kCand*` bit positions
 *       (not longcallD `0x004`/`0x008`/… hex values); they are self-consistent via
 *       `category_to_flag`. This function does **not** set VCF-style output; pgPhase additionally
 *       fills `hap_alt`/`hap_ref` after the k-means loop using logic shaped like longcallD
 *       `make_variants`, which downstream longcallD does in a separate pass.
 *
 * **Phase 1** — initial sweep from the highest-confidence pivot variant outward,
 * assigning each read to hap 1 or 2 based on the current per-variant consensus
 * allele profile (`hap_to_cons_alle`), and updating the profile after each assignment.
 *
 * **Phase 2** — up to 10 k-means iterations:
 *  - `iter_update_var_hap_cons_phase_set`: detect phase-set breaks (too few spanning reads)
 *    and flip `hap_to_cons_alle[1]/[2]` when conflict reads outnumber agreement reads.
 *  - `iter_update_var_hap_to_cons_alle`: re-assign all reads from scratch and rebuild
 *    allele profiles; stop early if consensus alleles converge.
 *
 * **Phase 3** — assign `chunk.phase_sets[read_i]` from the first phased het variant
 * each read overlaps.
 *
 * **Phase 4** — fill `hap_alt` / `hap_ref` on every `CandidateVariant` from the
 * finalized `hap_to_cons_alle[1..2]` (0=ref, 1=alt, −1=unresolved).
 *
 * @param chunk    Working chunk; `read_var_profile` and `read_var_cr` must be populated
 *                 by `collect_read_var_profile` before calling this function.
 * @param opts     Thresholds — `is_ont()` controls the ONT homopolymer 67% guard.
 * @param flags    Bitmask of `kCand*` constants selecting which categories participate.
 */
void assign_hap_based_on_germline_het_vars_kmeans(BamChunk& chunk, const Options& opts, uint32_t flags);

} // namespace pgphase_collect

#endif
