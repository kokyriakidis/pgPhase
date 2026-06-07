#ifndef PGPHASE_COLLECT_PHASE_HPP
#define PGPHASE_COLLECT_PHASE_HPP

// k-means read-haplotype clustering: iterative read→haplotype assignment
// driven by per-variant consensus allele profiles (hap_to_cons_alle),
// up to 10 refinement rounds, plus cross-chunk stitching.

#include "collect_types.hpp"

#include <cstdint>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Candidate-category bitmask flags
// ════════════════════════════════════════════════════════════════════════════

// Candidate-category bitmask flags for lcd_var_i_to_cate and k-means target selection.
// Individual category bits:
constexpr uint32_t kLongcalldLowCovVar      = 0x001u;  // low coverage / low alt depth
constexpr uint32_t kLongcalldStrandBiasVar   = 0x002u;  // strand bias
constexpr uint32_t kCandCleanHetSnp          = 0x004u;  // clean het SNP
constexpr uint32_t kCandCleanHetIndel        = 0x008u;  // clean het indel
constexpr uint32_t kLongcalldRepHetVar       = 0x010u;  // repeat/homopolymer het indel
constexpr uint32_t kLongcalldCandSomaticVar  = 0x040u;  // candidate somatic variant
constexpr uint32_t kCandCleanHom             = 0x080u;  // clean homozygous
constexpr uint32_t kCandNoisyCandHet         = 0x100u;  // noisy-region MSA het
constexpr uint32_t kCandNoisyCandHom         = 0x200u;  // noisy-region MSA hom
constexpr uint32_t kLongcalldLowAfVar        = 0x400u;  // low allele fraction
constexpr uint32_t kLongcalldNonVar          = 0x800u;  // non-variant placeholder

// Composite masks:
// Categories excluded from noisy-region containment checks.
constexpr uint32_t kLongcalldNotCandVarCate =
    kLongcalldNonVar | kLongcalldLowCovVar | kLongcalldStrandBiasVar;
// All het candidate categories.
constexpr uint32_t kCandHetVarCate =
    kCandCleanHetSnp | kCandCleanHetIndel | kCandNoisyCandHet;
// Clean germline categories (used for VCF INFO CLEAN flag).
constexpr uint32_t kCandGermlineClean = kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom;
// All germline categories (clean + noisy-recalled).
constexpr uint32_t kCandGermlineVarCate =
    kCandGermlineClean | kCandNoisyCandHet | kCandNoisyCandHom;

// ════════════════════════════════════════════════════════════════════════════
// Public interface
// ════════════════════════════════════════════════════════════════════════════

// Maps a VariantCategory enum to its bitmask flag (0 for unmapped categories).
uint32_t category_to_flag(VariantCategory c);

// Stitch haplotype assignments across chunk boundaries.  For each adjacent
// pair on the same contig, overlapping boundary reads vote on whether the
// downstream chunk's hap labels should be flipped.  When the vote is decisive,
// candidate and read-level hap/PS assignments are updated.
//
// When pgbam_sidecar is provided, annotated-BAM thread IDs are used to merge
// phase blocks that lack decisive common-read overlap.
void stitch_chunk_haps(std::vector<PhasingChunk>& chunks,
                       const Options* opts = nullptr,
                       const PgbamSidecarData* pgbam_sidecar = nullptr);

// Assign haplotypes and phase sets to reads via iterative k-means clustering.
//
// Phase 1: initial sweep from the highest-confidence pivot variant outward,
//   assigning each read to hap 1 or 2 based on per-variant consensus allele
//   profiles (hap_to_cons_alle), updating profiles after each assignment.
//
// Phase 2: up to 10 k-means iterations — detect phase-set breaks, flip
//   consensus alleles when conflict reads outnumber agreement, re-assign all
//   reads and rebuild profiles.  Stops early on convergence.
//
// Phase 3: assign chunk.phase_sets[read_i] from the first phased het variant
//   each read overlaps.
//
// Phase 4: fill hap_alt/hap_ref on every CandidateVariant from the finalized
//   consensus alleles (0=ref, 1=alt, -1=unresolved).
//
// `flags` is a bitmask of kCand* constants selecting which candidate categories
// participate.  Masking uses lcd_var_i_to_cate, not counts.category.
// read_var_profile and read_var_cr must be populated before calling.
void assign_hap_based_on_germline_het_vars_kmeans(PhasingChunk& chunk, const Options& opts, uint32_t flags);

} // namespace pgphase_collect

#endif
