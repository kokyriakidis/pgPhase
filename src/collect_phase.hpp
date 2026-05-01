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
 * Values and composites are copied verbatim from longcallD `collect_var.h` (lines 11–27): the same
 * integers used for `var_i_to_cate[]` and `target_var_cate` in `assign_hap.c`.
 */
///@{
/** `LONGCALLD_LOW_COV_VAR` */
constexpr uint32_t kLongcalldLowCovVar = 0x001u;
/** `LONGCALLD_STRAND_BIAS_VAR` */
constexpr uint32_t kLongcalldStrandBiasVar = 0x002u;
/** `LONGCALLD_CLEAN_HET_SNP` */
constexpr uint32_t kCandCleanHetSnp = 0x004u;
/** `LONGCALLD_CLEAN_HET_INDEL` */
constexpr uint32_t kCandCleanHetIndel = 0x008u;
/** `LONGCALLD_REP_HET_VAR` */
constexpr uint32_t kLongcalldRepHetVar = 0x010u;
/** `LONGCALLD_CAND_SOMATIC_VAR` */
constexpr uint32_t kLongcalldCandSomaticVar = 0x040u;
/** `LONGCALLD_CLEAN_HOM_VAR` */
constexpr uint32_t kCandCleanHom = 0x080u;
/** `LONGCALLD_NOISY_CAND_HET_VAR` */
constexpr uint32_t kCandNoisyCandHet = 0x100u;
/** `LONGCALLD_NOISY_CAND_HOM_VAR` */
constexpr uint32_t kCandNoisyCandHom = 0x200u;
/** `LONGCALLD_LOW_AF_VAR` */
constexpr uint32_t kLongcalldLowAfVar = 0x400u;
/** `LONGCALLD_NON_VAR` */
constexpr uint32_t kLongcalldNonVar = 0x800u;

/** `LONGCALLD_NOT_CAND_VAR_CATE` — skipped before `cr_is_contained` in classify compaction. */
constexpr uint32_t kLongcalldNotCandVarCate =
    kLongcalldNonVar | kLongcalldLowCovVar | kLongcalldStrandBiasVar;

/** `LONGCALLD_CAND_HET_VAR_CATE` */
constexpr uint32_t kCandHetVarCate =
    kCandCleanHetSnp | kCandCleanHetIndel | kCandNoisyCandHet;

/** `LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE` */
constexpr uint32_t kCandGermlineClean = kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom;

/** `LONGCALLD_CAND_GERMLINE_VAR_CATE` */
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
 * Matches longcallD `stitch_var_main` exactly: for each adjacent pair on the same contig,
 * calls `flip_variant_hap` logic only — overlapping boundary reads (`down_ovlp_read_i` /
 * `up_ovlp_read_i`, same pairing order as longcallD load) vote on orientation; when the
 * vote is decisive (`flip_hap_score != 0`), applies `update_chunk_var_hap_phase_set1` and,
 * only if phased alignment output is requested (`opts != nullptr &&
 * !opts->output_aln.empty()`, analogous to `opt->out_aln_fp != NULL`), read-level hap/PS updates.
 *
 * The `pgbam_sidecar` argument is ignored (kept for call-site compatibility); longcallD does
 * not perform pangenome-graph stitching or within-chunk phase-block merges here.
 *
 * @param chunks Ordered, adjacent `BamChunk` objects, each already phased by k-means.
 */
void stitch_chunk_haps(std::vector<BamChunk>& chunks,
                       const Options* opts = nullptr,
                       const PgbamSidecarData* pgbam_sidecar = nullptr);

/**
 * @brief Assign haplotypes and phase sets to reads via iterative k-means clustering.
 *
 * Mirrors longcallD `assign_hap_based_on_germline_het_vars_kmeans` from `assign_hap.c`.
 *
 * @note Parity is with **assign_hap.c**: `flags` is a `LONGCALLD_*` composite (`target_var_cate`);
 *       masking uses `CandidateVariant::lcd_var_i_to_cate` (longcallD `var_i_to_cate[i]`), not
 *       `counts.category` alone — matching `assign_hap_based_on_germline_het_vars_kmeans` /
 *       `init_assign_read_hap_based_on_cons_alle` / `select_init_var` / profile updates.
 *       This function does **not** set VCF-style output; pgPhase additionally
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
