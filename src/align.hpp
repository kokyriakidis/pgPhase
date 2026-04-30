#ifndef PGPHASE_ALIGN_HPP
#define PGPHASE_ALIGN_HPP

/**
 * @file align.hpp
 * @brief WFA2-lib / abPOA alignment wrappers and noisy-region read-info helpers.
 *
 * Mirrors longcallD align.h / align.c.  All functions that directly call WFA2-lib
 * or abPOA live here; collect_phase_noisy.cpp calls into this file.
 */

#include "collect_types.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Cover-flag constants  (mirrors longcallD align.h LONGCALLD_NOISY_* macros)
// ════════════════════════════════════════════════════════════════════════════

constexpr int kNoisyNoCover    = 0x0000;
constexpr int kNoisyRightGap   = 0x0001;
constexpr int kNoisyLeftGap    = 0x0002;
constexpr int kNoisyRightCover = 0x0004; // 4, 6
constexpr int kNoisyLeftCover  = 0x0008; // 8, 9
constexpr int kNoisyBothCover  = 0x000C; // 12

inline bool noisyIsBothCover(int c)  { return (c & kNoisyLeftCover) && (c & kNoisyRightCover); }
inline bool noisyIsLeftCover(int c)  { return (c & kNoisyLeftCover) && !(c & kNoisyRightCover); }
inline bool noisyIsRightCover(int c) { return !(c & kNoisyLeftCover) && (c & kNoisyRightCover); }
inline bool noisyIsLeftGap(int c)    { return (c & kNoisyLeftGap) != 0; }
inline bool noisyIsRightGap(int c)   { return (c & kNoisyRightGap) != 0; }
inline bool noisyIsNotCover(int c)   { return !(c & kNoisyLeftCover) && !(c & kNoisyRightCover); }

// ════════════════════════════════════════════════════════════════════════════
// Alignment scoring  (mirrors longcallD align.h / call_var_main.h constants)
// ════════════════════════════════════════════════════════════════════════════

constexpr int kMatchScore     = 2;
constexpr int kMismatchScore  = 6;
constexpr int kGapOpen1Score  = 6;
constexpr int kGapExt1Score   = 2;
constexpr int kGapOpen2Score  = 24;
constexpr int kGapExt2Score   = 1;

constexpr int kGapLeftAln  = 1; // default: gaps placed at left-most position
constexpr int kGapRightAln = 2;

constexpr int kExtAlnLeftToRight = 1;
constexpr int kExtAlnRightToLeft = 2;

constexpr int kWfaNoHeuristic = 0;
constexpr int kWfaAdaptive    = 1;
constexpr int kWfaZdrop       = 2;

constexpr int kWfaAffine1p = 0;
constexpr int kWfaAffine2p = 1;

constexpr double kPartialAlnRatio = 1.1; // max longer/shorter ratio for partial alignment

// ════════════════════════════════════════════════════════════════════════════
// AlnStr  (mirrors longcallD aln_str_t in collect_var.h)
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Pairwise alignment string for one (target, query) pair from the MSA.
 *
 * target_aln / query_aln are parallel byte arrays; each byte is a 2-bit base
 * (0–3 = ACGT) or 5 for a gap character, matching longcallD's convention.
 */
struct AlnStr {
    std::vector<uint8_t> target_aln; ///< reference / consensus bases (5 = gap)
    std::vector<uint8_t> query_aln;  ///< read / consensus bases    (5 = gap)
    int aln_len    = 0;
    int target_beg = 0, target_end = 0;
    int query_beg  = 0, query_end  = 0;
};

// ════════════════════════════════════════════════════════════════════════════
// Per-region read information
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Per-read data extracted from a noisy region for MSA input.
 *
 * Mirrors what longcallD collect_noisy_read_info allocates and fills.
 * All arrays are parallel and indexed [0, n_reads).
 */
struct NoisyReadInfo {
    int n_reads = 0;
    std::vector<int>                noisy_read_ids; ///< indices into chunk.reads
    std::vector<int>                lens;
    std::vector<std::vector<uint8_t>> seqs;
    std::vector<std::vector<uint8_t>> quals;
    std::vector<uint8_t>            strands;
    std::vector<int>                fully_covers;
    std::vector<int>                haps;
    std::vector<hts_pos_t>          phase_sets;
    // per-read (read_id → cover/beg/end) maps indexed by read_id for digar update
    std::vector<int>                read_id_to_full_covers; // size = chunk.reads.size()
    std::vector<int>                read_reg_beg;           // size = chunk.reads.size()
    std::vector<int>                read_reg_end;           // size = chunk.reads.size()
};

// ════════════════════════════════════════════════════════════════════════════
// Function declarations
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Compute per-base expected error rate from Phred quality bytes.
 * Mirrors longcallD calc_read_error_rate (seq.c).
 */
double calc_read_error_rate(int len, const uint8_t* qual);

/**
 * @brief Compare two cover flags: BOTH_COVER > LEFT/RIGHT_COVER > NO_COVER.
 * Mirrors longcallD full_cover_cmp (align.c).
 */
int full_cover_cmp(int cover1, int cover2);

/**
 * @brief Extract per-read sequences and cover metadata for a noisy region.
 *
 * Mirrors longcallD collect_noisy_read_info (align.c lines 1377–1461).
 * Walks digars to locate query positions spanning [reg_beg, reg_end], extracts
 * bases via bam_get_seq, computes the fully_covers bitmask, and records haps/phase_sets.
 */
NoisyReadInfo collect_noisy_read_info(const Options& opts, const BamChunk& chunk,
                                      hts_pos_t reg_beg, hts_pos_t reg_end,
                                      const std::vector<int>& noisy_read_ids);

/**
 * @brief Sort reads in a noisy region: BOTH_COVER first, then by error rate (opt.), then by length desc.
 * Mirrors longcallD sort_noisy_region_reads (align.c).
 * Mutates info in place (parallel arrays).
 */
void sort_noisy_region_reads(NoisyReadInfo& info, bool use_error_rate);

/**
 * @brief Find the best phase set with >= min_hap_full_reads per haplotype.
 *
 * Mirrors longcallD collect_phase_set_with_both_haps (align.c).
 * Returns -1 when no valid phase set exists.
 */
hts_pos_t collect_phase_set_with_both_haps(const NoisyReadInfo& info,
                                           int min_hap_full_reads,
                                           int min_hap_all_reads);

/**
 * @brief Trim AlnStr for partially covering reads (left- or right-cover only).
 * Mirrors longcallD wfa_trim_aln_str (align.c).
 */
void wfa_trim_aln_str(int full_cover, AlnStr& aln_str);

/**
 * @brief WFA end-to-end alignment of pattern vs text, returns alignment strings.
 *
 * Mirrors longcallD wfa_end2end_aln (align.c).
 * pattern_alg / text_alg are filled with 2-bit bases (5=gap); caller owns memory.
 */
int wfa_end2end_aln(const uint8_t* pattern, int plen, const uint8_t* text, int tlen,
                    int gap_aln, int b, int q, int e, int q2, int e2,
                    int heuristic, int affine_gap,
                    std::vector<uint8_t>& pattern_alg, std::vector<uint8_t>& text_alg);

/**
 * @brief Align target vs query into an AlnStr, handling both full- and partial-cover reads.
 * Mirrors longcallD wfa_collect_aln_str (align.c).
 */
int wfa_collect_aln_str(const Options& opts, const uint8_t* target, int tlen,
                        const uint8_t* query, int qlen, int full_cover,
                        int heuristic, int affine_gap, AlnStr& aln_str);

/**
 * @brief Build cons-vs-read AlnStr from abPOA MSA result rows.
 * Mirrors longcallD make_cons_read_aln_str (align.c).
 */
int make_cons_read_aln_str(const uint8_t* cons_msa_row, const uint8_t* read_msa_row,
                           int msa_len, int full_cover, AlnStr& out);

/**
 * @brief Build ref-vs-read AlnStr by merging ref-cons and cons-read AlnStrs.
 * Mirrors longcallD make_ref_read_aln_str (align.c).
 */
int make_ref_read_aln_str(const Options& opts, const AlnStr& ref_cons,
                          const AlnStr& cons_read, AlnStr& ref_read);

/**
 * @brief Run abPOA de-novo MSA and produce ≤2 consensus sequences (haplotype-unaware path).
 *
 * Mirrors longcallD abpoa_aln_msa_cons (align.c).
 * Returns number of consensus sequences (0, 1, or 2).
 * clu_n_seqs[i] = # reads assigned to cluster i; clu_read_ids[i] = their original read indices.
 * msa_seqs[i][j] = MSA row (length msa_len[i]) for read j within cluster i (j == clu_n_seqs[i] → consensus row).
 */
int abpoa_aln_msa_cons(const Options& opts, int n_reads,
                       const std::vector<int>& read_ids,
                       const std::vector<uint8_t*>& read_seqs,
                       const std::vector<int>& read_lens,
                       int max_n_cons,
                       std::vector<int>& cons_lens,
                       std::vector<std::vector<uint8_t>>& cons_seqs,
                       std::array<int, 2>& clu_n_seqs,
                       std::array<std::vector<int>, 2>& clu_read_ids,
                       std::vector<int>& msa_seq_lens,
                       std::array<std::vector<std::vector<uint8_t>>, 2>& msa_seqs);

/**
 * @brief Run abPOA partial MSA for a single haplotype and produce one consensus.
 *
 * Mirrors longcallD abpoa_partial_aln_msa_cons (align.c).
 * Returns 1 on success, 0 if no consensus was produced.
 */
int abpoa_partial_aln_msa_cons(const Options& opts, int sampling_reads,
                               int n_reads,
                               const std::vector<int>& read_ids,
                               const std::vector<uint8_t*>& read_seqs,
                               const std::vector<uint8_t*>& read_quals,
                               const std::vector<int>& read_lens,
                               const std::vector<int>& read_full_covers,
                               int max_n_cons,
                               int& cons_len,
                               std::vector<uint8_t>& cons_seq,
                               int& out_clu_n_seqs,
                               std::vector<int>& out_clu_read_ids,
                               int& msa_seq_len,
                               std::vector<std::vector<uint8_t>>& msa_seqs);

/**
 * @brief Build MSA alignment strings for a noisy region (haplotype-unaware path).
 *
 * Mirrors longcallD wfa_collect_noisy_aln_str_no_ps_hap (align.c).
 * Returns number of consensus sequences (0, 1, or 2).
 */
int wfa_collect_noisy_aln_str_no_ps_hap(const Options& opts, NoisyReadInfo& info,
                                         const uint8_t* ref_seq, int ref_seq_len,
                                         bool collect_ref_read_aln_str,
                                         std::array<int, 2>& clu_n_seqs,
                                         std::array<std::vector<int>, 2>& clu_read_ids,
                                         std::array<std::vector<AlnStr>, 2>& aln_strs);

/**
 * @brief Build MSA alignment strings for a noisy region (haplotype-aware path).
 *
 * Mirrors longcallD wfa_collect_noisy_aln_str_with_ps_hap (align.c).
 * Returns 2 on success (one consensus per haplotype), 0 otherwise.
 */
int wfa_collect_noisy_aln_str_with_ps_hap(const Options& opts, bool sampling_reads,
                                           NoisyReadInfo& info,
                                           hts_pos_t ps,
                                           int min_hap_full_reads, int min_hap_all_reads,
                                           const uint8_t* ref_seq, int ref_seq_len,
                                           bool collect_ref_read_aln_str,
                                           std::array<int, 2>& clu_n_seqs,
                                           std::array<std::vector<int>, 2>& clu_read_ids,
                                           std::array<std::vector<AlnStr>, 2>& aln_strs);

/**
 * @brief Top-level entry: collect reads → sort → find phase set → MSA → fill aln_strs.
 *
 * Mirrors longcallD collect_noisy_reg_aln_strs (align.c / align.h).
 * Returns number of consensus sequences produced (0, 1, or 2).
 */
int collect_noisy_reg_aln_strs(const Options& opts, BamChunk& chunk,
                                hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                                const std::vector<int>& noisy_read_ids,
                                const std::vector<uint8_t>& ref_seq,
                                std::array<int, 2>& clu_n_seqs,
                                std::array<std::vector<int>, 2>& clu_read_ids,
                                std::array<std::vector<AlnStr>, 2>& aln_strs);

} // namespace pgphase_collect

#endif
