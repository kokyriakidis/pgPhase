/**
 * @file collect_var.cpp
 * @brief Candidate sites, allele counting, noisy regions, and longcallD-shaped classification.
 */

#include "collect_var.hpp"
#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"

#include "sdust.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace pgphase_collect {

namespace {

// longcallD `seq.c` `nst_nt4_table` — encodes SNP/INS `alt_seq` bytes used in `memcmp` by
// `exact_comp_var_site` / `exact_comp_var_site_ins` (`collect_var.c`).
alignas(64) static constexpr uint8_t kLcdNstNt4Table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

/** `memcmp(var1->alt_seq, var2->alt_seq, len)` on longcallD-encoded ALT bytes (LCD table per char). */
int cmp_alt_seq_lcd_nt4(const std::string& a, const std::string& b, int len) {
    for (int i = 0; i < len; ++i) {
        const uint8_t ca = kLcdNstNt4Table[static_cast<unsigned char>(a[static_cast<size_t>(i)])];
        const uint8_t cb = kLcdNstNt4Table[static_cast<unsigned char>(b[static_cast<size_t>(i)])];
        if (ca != cb) return (ca < cb) ? -1 : 1;
    }
    return 0;
}

static int cmp_alt_bytes_lcd_nt4(const uint8_t* a, const uint8_t* b, int len) {
    for (int i = 0; i < len; ++i) {
        const uint8_t ca = kLcdNstNt4Table[static_cast<unsigned char>(a[i])];
        const uint8_t cb = kLcdNstNt4Table[static_cast<unsigned char>(b[i])];
        if (ca != cb) return (ca < cb) ? -1 : 1;
    }
    return 0;
}

/** Shape of longcallD `var_site_t` for `update_cand_vars_from_digar` merge only (`collect_var.c`). */
struct LcdVarSiteMerge {
    hts_pos_t pos = 0;
    int var_type = 0;
    int ref_len = 0;
    int alt_len = 0;
    const uint8_t* alt_seq = nullptr;
};

static LcdVarSiteMerge lcd_merge_site_from_candidate_key(const VariantKey& k) {
    LcdVarSiteMerge s{};
    s.pos = k.pos;
    if (k.type == VariantType::Snp) {
        s.var_type = BAM_CDIFF;
        s.ref_len = 1;
        s.alt_len = static_cast<int>(k.alt.size());
        s.alt_seq = reinterpret_cast<const uint8_t*>(k.alt.data());
    } else if (k.type == VariantType::Insertion) {
        s.var_type = BAM_CINS;
        s.ref_len = 0;
        s.alt_len = static_cast<int>(k.alt.size());
        s.alt_seq = reinterpret_cast<const uint8_t*>(k.alt.data());
    } else {
        s.var_type = BAM_CDEL;
        s.ref_len = k.ref_len;
        s.alt_len = 0;
        s.alt_seq = nullptr;
    }
    return s;
}

static int digar_to_bam_var_type(DigarType t) {
    switch (t) {
        case DigarType::Snp:
            return BAM_CDIFF;
        case DigarType::Insertion:
            return BAM_CINS;
        case DigarType::Deletion:
            return BAM_CDEL;
        case DigarType::SoftClip:
            return BAM_CSOFT_CLIP;
        case DigarType::HardClip:
            return BAM_CHARD_CLIP;
        case DigarType::RefSkip:
            return BAM_CREF_SKIP;
        default:
            return BAM_CEQUAL;
    }
}

/** longcallD `make_var_site_from_digar` (`collect_var.c`) for merge-site fields only. */
static LcdVarSiteMerge lcd_merge_site_from_digar_op(const DigarOp& op) {
    LcdVarSiteMerge s{};
    s.pos = op.pos;
    switch (op.type) {
        case DigarType::Snp:
            s.var_type = BAM_CDIFF;
            s.ref_len = 1;
            s.alt_len = static_cast<int>(op.alt.size());
            s.alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return s;
        case DigarType::Insertion:
            s.var_type = BAM_CINS;
            s.ref_len = 0;
            s.alt_len = static_cast<int>(op.alt.size());
            s.alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return s;
        case DigarType::Deletion:
            s.var_type = BAM_CDEL;
            s.ref_len = op.len;
            s.alt_len = 0;
            s.alt_seq = nullptr;
            return s;
        default:
            s.var_type = digar_to_bam_var_type(op.type);
            s.ref_len = 1;
            s.alt_len = op.len;
            s.alt_seq = nullptr;
            return s;
    }
}

/** longcallD `exact_comp_var_site_ins(opt, …)` (`collect_var.c`) on merge-shaped sites (no `tid`). */
static int lcd_exact_comp_var_site_ins_merge(const LcdVarSiteMerge* var1, const LcdVarSiteMerge* var2, int min_sv_len) {
    const hts_pos_t var1_pos = var1->var_type == BAM_CDIFF ? var1->pos : var1->pos - 1;
    const hts_pos_t var2_pos = var2->var_type == BAM_CDIFF ? var2->pos : var2->pos - 1;
    if (var1_pos < var2_pos) return -1;
    if (var1_pos > var2_pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->var_type == BAM_CDIFF) {
        if (var1->alt_len < var2->alt_len) return -1;
        if (var1->alt_len > var2->alt_len) return 1;
        return cmp_alt_bytes_lcd_nt4(var1->alt_seq, var2->alt_seq, var1->alt_len);
    }
    if (var1->var_type == BAM_CINS) {
        if (var1->alt_len < min_sv_len) {
            if (var1->alt_len < var2->alt_len) return -1;
            if (var1->alt_len > var2->alt_len) return 1;
            return cmp_alt_bytes_lcd_nt4(var1->alt_seq, var2->alt_seq, var1->alt_len);
        }
        const int min_len = var1->alt_len < var2->alt_len ? var1->alt_len : var2->alt_len;
        const int max_len = var1->alt_len > var2->alt_len ? var1->alt_len : var2->alt_len;
        if (static_cast<double>(min_len) >= static_cast<double>(max_len) * 0.8) return 0;
        return var1->alt_len - var2->alt_len;
    }
    return 0;
}

} // namespace

// ════════════════════════════════════════════════════════════════════════════
// Variant site comparison
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Effective alternate length for sorting and comparison (longcallD `var_site_t` convention).
 *
 * Deletions use length 0; SNPs and insertions use `alt.size()`.
 * @param v Candidate key.
 * @return Alternate length used in comparators.
 */
static int var_site_alt_len(const VariantKey& v) {
    if (v.type == VariantType::Deletion) return 0;
    return static_cast<int>(v.alt.size());
}

/**
 * @brief Total order for candidate keys (longcallD `exact_comp_var_site`).
 *
 * Compares tid, sort position, type, `ref_len`, and alternate sequence length / bytes.
 * Return \<0 if \a var1 \< \a var2, 0 if equal, \>0 if \a var1 \> \a var2.
 *
 * @param var1 First variant key.
 * @param var2 Second variant key.
 * @return \<0, 0, or \>0 according to total order.
 */
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2) {
    if (var1->tid != var2->tid) return var1->tid < var2->tid ? -1 : 1;
    const hts_pos_t p1 = var1->type == VariantType::Snp ? var1->pos : var1->pos - 1;
    const hts_pos_t p2 = var2->type == VariantType::Snp ? var2->pos : var2->pos - 1;
    if (p1 < p2) return -1;
    if (p1 > p2) return 1;
    const int t1 = static_cast<int>(var1->type);
    const int t2 = static_cast<int>(var2->type);
    if (t1 < t2) return -1;
    if (t1 > t2) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    const int a1 = var_site_alt_len(*var1);
    const int a2 = var_site_alt_len(*var2);
    if (a1 < a2) return -1;
    if (a1 > a2) return 1;
    if (var1->type == VariantType::Snp || var1->type == VariantType::Insertion) {
        return cmp_alt_seq_lcd_nt4(var1->alt, var2->alt, a1);
    }
    return 0;
}

int exact_comp_cand_var(const CandidateVariant* var1, const CandidateVariant* var2) {
    return exact_comp_var_site(&var1->key, &var2->key);
}

/**
 * @brief Variant comparison with large-insertion fuzzy merge (longcallD `exact_comp_var_site_ins`).
 *
 * Same ordering keys as `exact_comp_var_site` for SNPs and deletions. For insertions with
 * `alt` length ≥ \a min_sv_len, treats two alleles as equal when
 * `min(len)/max(len) ≥ 0.8`; shorter insertions still require exact sequence match.
 *
 * @param var1 First variant key.
 * @param var2 Second variant key.
 * @param min_sv_len Minimum insertion length to apply the length-ratio rule.
 * @return \<0, 0, or \>0; zero means merge-equivalent under fuzzy INS rule.
 */
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len) {
    if (var1->tid != var2->tid) return var1->tid < var2->tid ? -1 : 1;
    const hts_pos_t p1 = var1->type == VariantType::Snp ? var1->pos : var1->pos - 1;
    const hts_pos_t p2 = var2->type == VariantType::Snp ? var2->pos : var2->pos - 1;
    if (p1 < p2) return -1;
    if (p1 > p2) return 1;
    const int t1 = static_cast<int>(var1->type);
    const int t2 = static_cast<int>(var2->type);
    if (t1 < t2) return -1;
    if (t1 > t2) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->type == VariantType::Snp) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < a2) return -1;
        if (a1 > a2) return 1;
        return cmp_alt_seq_lcd_nt4(var1->alt, var2->alt, a1);
    }
    if (var1->type == VariantType::Insertion) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < min_sv_len) {
            if (a1 < a2) return -1;
            if (a1 > a2) return 1;
            return cmp_alt_seq_lcd_nt4(var1->alt, var2->alt, a1);
        }
        const int min_len = a1 < a2 ? a1 : a2;
        const int max_len = a1 > a2 ? a1 : a2;
        if (static_cast<double>(min_len) >= static_cast<double>(max_len) * 0.8) return 0;
        return a1 - a2;
    }
    return 0;
}

/**
 * @brief Strict-weak ordering for `VariantKey` using `exact_comp_var_site`.
 * @return True if `lhs` is ordered before `rhs`.
 */
bool VariantKeyLess::operator()(const VariantKey& lhs, const VariantKey& rhs) const {
    return exact_comp_var_site(&lhs, &rhs) < 0;
}

// ════════════════════════════════════════════════════════════════════════════
// Candidate site collection
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Builds a `VariantKey` from one digar op (longcallD `make_var_site_from_digar`).
 *
 * SNP: `ref_len` 1; insertion: `ref_len` 0; deletion: `ref_len` = `op.len`, empty `alt`.
 *
 * @param tid Contig index for the read.
 * @param op Parsed SNP/indel operation.
 * @return Normalized `VariantKey` for the candidate table.
 */
VariantKey variant_key_from_digar(int tid, const DigarOp& op) {
    if (op.type == DigarType::Snp) {
        return VariantKey{tid, op.pos, VariantType::Snp, 1, op.alt};
    }
    if (op.type == DigarType::Insertion) {
        return VariantKey{tid, op.pos, VariantType::Insertion, 0, op.alt};
    }
    return VariantKey{tid, op.pos, VariantType::Deletion, op.len, ""};
}

/**
 * @brief Whether a digar op is a non–low-quality SNP/indel inside the region bounds.
 *
 * Matches longcallD `is_collectible_var_digar`; use `reg_beg` or `reg_end` == -1 to disable that bound.
 *
 * @param digar Parsed operation.
 * @param reg_beg Minimum 1-based position, or -1 for no lower bound.
 * @param reg_end Maximum 1-based position, or -1 for no upper bound.
 * @return True if SNP/indel, not low-quality, and inside bounds.
 */
static bool is_collectible_var_digar(const DigarOp& digar, hts_pos_t reg_beg, hts_pos_t reg_end) {
    const hts_pos_t digar_pos = digar.pos;
    if (reg_beg != -1 && digar_pos < reg_beg) return false;
    if (reg_end != -1 && digar_pos > reg_end) return false;
    if (digar.low_quality) return false;
    return digar.type == DigarType::Snp || digar.type == DigarType::Insertion ||
           digar.type == DigarType::Deletion;
}

/**
 * @brief Lower-bound index into sorted candidates for a read span (longcallD `get_var_site_start`).
 * @param v Sorted candidates by `sort_pos` / `pos`.
 * @param start Read start (1-based).
 * @param n_total `static_cast<int>(v.size())`.
 * @return First index to consider in the merge walk.
 */
static int get_var_site_start(const CandidateTable& v, hts_pos_t start, int n_total) {
    hts_pos_t target = start > 0 ? start - 1 : start;
    int left = 0;
    int right = n_total;
    while (left < right) {
        const int mid = left + (right - left) / 2;
        const hts_pos_t mid_pos = v[static_cast<size_t>(mid)].key.sort_pos();
        if (mid_pos < target) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    while (left < n_total && v[static_cast<size_t>(left)].key.pos < start) {
        ++left;
    }
    return left;
}

/**
 * @brief Deduplicates sorted candidates using `exact_comp_var_site_ins` (longcallD merge pass).
 *
 * Sorts by `exact_comp_var_site`, then keeps one row per equivalence class where fuzzy large
 * insertions may collapse. Equivalent to the merge step in longcallD `collect_all_cand_var_sites`.
 *
 * @param variants In/out table; replaced with unique keys in sort order.
 */
void collapse_fuzzy_large_insertions(CandidateTable& variants, int min_sv_len) {
    std::sort(variants.begin(), variants.end(), [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_var_site(&a.key, &b.key) < 0;
    });
    if (variants.empty()) return;
    size_t write_i = 1;
    for (size_t read_i = 1; read_i < variants.size(); ++read_i) {
        if (exact_comp_var_site_ins(
                &variants[write_i - 1].key, &variants[read_i].key, min_sv_len) == 0) {
            continue;
        }
        if (write_i != read_i) variants[write_i] = std::move(variants[read_i]);
        ++write_i;
    }
    variants.resize(write_i);
}

void collapse_exact_duplicate_variants(CandidateTable& variants) {
    std::sort(variants.begin(), variants.end(), [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_var_site(&a.key, &b.key) < 0;
    });
    if (variants.empty()) return;
    size_t write_i = 1;
    for (size_t read_i = 1; read_i < variants.size(); ++read_i) {
        if (exact_comp_var_site(&variants[write_i - 1].key, &variants[read_i].key) == 0) {
            variants[write_i - 1].lcd_make_variants_region_pass |= variants[read_i].lcd_make_variants_region_pass;
            continue;
        }
        if (write_i != read_i) variants[write_i] = std::move(variants[read_i]);
        ++write_i;
    }
    variants.resize(write_i);
}

/**
 * @brief Builds the raw candidate table from per-read digars for one chunk.
 *
 * Emits one `CandidateVariant` per collectible SNP/indel overlapping `chunk`, then runs
 * `collapse_fuzzy_large_insertions`. Matches the site-gathering half of longcallD
 * `collect_all_cand_var_sites`.
 *
 * @param chunk Genomic interval (1-based inclusive bounds).
 * @param reads Parsed alignments with digar ops.
 * @param variants Cleared and filled with deduplicated candidate keys.
 */
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants,
                                          int min_sv_len) {
    variants.clear();
    size_t max_var_sites = 0;
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            ++max_var_sites;
        }
    }
    if (max_var_sites == 0) return;

    variants.reserve(max_var_sites);
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            variants.push_back(CandidateVariant{variant_key_from_digar(read.tid, digar), VariantCounts{}});
        }
    }
    collapse_fuzzy_large_insertions(variants, min_sv_len);
}

// ════════════════════════════════════════════════════════════════════════════
// Allele counting
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Mean base quality over bases implicated by a digar op (longcallD `get_digar_ave_qual`).
 * @param d Digar (uses `qi`, `len`, deletion semantics).
 * @param qual Per-read QUAL array.
 * @return Rounded mean Phred, or 0 if not applicable.
 */
static int get_digar_ave_qual(const DigarOp& d, const std::vector<uint8_t>& qual) {
    if (d.low_quality) return 0;
    if (d.qi < 0) return 0;
    if (qual.empty()) return 0;
    int q_start = 0;
    int q_end = 0;
    if (d.type == DigarType::Deletion) {
        if (d.qi == 0) {
            q_start = 0;
            q_end = 0;
        } else {
            q_start = d.qi - 1;
            q_end = d.qi;
        }
    } else {
        q_start = d.qi;
        q_end = d.qi + d.len - 1;
    }
    if (q_start < 0) return 0;
    if (q_end >= static_cast<int>(qual.size())) return 0;
    int sum = 0;
    for (int i = q_start; i <= q_end; ++i) sum += static_cast<int>(qual[i]);
    return sum / (q_end - q_start + 1);
}

/**
 * @brief Updates strand-aware depth counters for one read observation at a site.
 *
 * Low-quality observations increment only `low_qual_cov`. Palindromic ONT reads contribute to
 * depth and ref/alt counts but not to forward/reverse strand tallies.
 *
 * @param counts Counters to update.
 * @param reverse Read reverse flag.
 * @param alt True if observation supports alternate allele.
 * @param low_quality True if counted toward `low_qual_cov` only.
 * @param is_ont_palindrome If true, strand subcounts are skipped.
 */
static void add_coverage(VariantCounts& counts, bool reverse, bool alt, bool low_quality,
                         bool is_ont_palindrome) {
    if (low_quality) {
        counts.low_qual_cov++;
        return;
    }
    counts.total_cov++;
    if (alt) {
        counts.alt_cov++;
        if (!is_ont_palindrome) reverse ? counts.reverse_alt++ : counts.forward_alt++;
    } else {
        counts.ref_cov++;
        if (!is_ont_palindrome) reverse ? counts.reverse_ref++ : counts.forward_ref++;
    }
}

/**
 * @brief Sweeps reads against sorted candidates to fill allele and strand counts.
 *
 * Implements longcallD `init_cand_vars_based_on_sites` + `update_cand_vars_from_digar`
 * (`bam_utils.c`): merge walk skipping only `BAM_CEQUAL` digars, comparing candidate vs digar
 * with `exact_comp_var_site_ins` on merge-shaped `var_site_t` (clips participate in ordering).
 * Low-quality alt when `digar.low_quality` or mean BQ \< \a min_bq; implicit ref for trailing
 * sites before `read.end`. Optionally appends `ReadSupportRow` rows when both
 * \a read_support_out and \a chunk_region are non-null.
 *
 * @param reads Parsed chunk reads.
 * @param variants Sorted candidate table (same order as site collection).
 * @param chunk_region If non-null with \a read_support_out, written into support rows.
 * @param read_support_out Optional per-read x site observation log.
 * @param min_bq Minimum base quality threshold for counting an alt as high-quality.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq,
                                        int min_sv_len) {
    const int n = static_cast<int>(variants.size());
    if (n == 0) return;

    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        if (read.tid < 0) continue;
        const hts_pos_t pos_start = read.beg;
        const hts_pos_t pos_end = read.end;
        int site_i = get_var_site_start(variants, pos_start, n);
        size_t digar_i = 0;
        const std::vector<DigarOp>& dig = read.digars;

        while (site_i < n && digar_i < dig.size()) {
            if (const int stid = variants[static_cast<size_t>(site_i)].key.tid; stid != read.tid) {
                if (stid < read.tid) {
                    ++site_i;
                } else {
                    break;
                }
                continue;
            }
            if (dig[digar_i].type == DigarType::Equal) {
                ++digar_i;
                continue;
            }
            const DigarOp& d = dig[digar_i];
            const int ave = get_digar_ave_qual(d, read.qual);
            const LcdVarSiteMerge cand_site = lcd_merge_site_from_candidate_key(variants[static_cast<size_t>(site_i)].key);
            const LcdVarSiteMerge dig_site = lcd_merge_site_from_digar_op(d);
            const int ret = lcd_exact_comp_var_site_ins_merge(&cand_site, &dig_site, min_sv_len);
            if (ret < 0) {
                VariantCounts& c = variants[static_cast<size_t>(site_i)].counts;
                add_coverage(c, read.reverse, false, false, read.is_ont_palindrome);
                if (read_support_out != nullptr && chunk_region != nullptr) {
                    const auto& k = variants[static_cast<size_t>(site_i)].key;
                    read_support_out->push_back(ReadSupportRow{k.tid,
                                                               k.pos,
                                                               k.type,
                                                               k.ref_len,
                                                               k.alt,
                                                               read.qname,
                                                               0,
                                                               0,
                                                               read.reverse,
                                                               read.mapq,
                                                               chunk_region->beg,
                                                               chunk_region->end});
                }
                ++site_i;
            } else if (ret == 0) {
                const bool lq = d.low_quality || ave < min_bq;
                add_coverage(
                    variants[static_cast<size_t>(site_i)].counts, read.reverse, true, lq,
                    read.is_ont_palindrome);
                if (read_support_out != nullptr && chunk_region != nullptr) {
                    const auto& k = variants[static_cast<size_t>(site_i)].key;
                    read_support_out->push_back(
                        ReadSupportRow{k.tid,
                                       k.pos,
                                       k.type,
                                       k.ref_len,
                                       k.alt,
                                       read.qname,
                                       1,
                                       lq ? 1 : 0,
                                       read.reverse,
                                       read.mapq,
                                       chunk_region->beg,
                                       chunk_region->end});
                }
                ++site_i;
            } else {
                ++digar_i;
            }
        }

        for (; site_i < n; ++site_i) {
            if (const int stid = variants[static_cast<size_t>(site_i)].key.tid; stid != read.tid) {
                if (stid < read.tid) continue;
                break;
            }
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) break;
            add_coverage(variants[static_cast<size_t>(site_i)].counts, read.reverse, false, false,
                         read.is_ont_palindrome);
            if (read_support_out != nullptr && chunk_region != nullptr) {
                const auto& k = variants[static_cast<size_t>(site_i)].key;
                read_support_out->push_back(ReadSupportRow{k.tid,
                                                           k.pos,
                                                           k.type,
                                                           k.ref_len,
                                                           k.alt,
                                                           read.qname,
                                                           0,
                                                           0,
                                                           read.reverse,
                                                           read.mapq,
                                                           chunk_region->beg,
                                                           chunk_region->end});
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Interval utilities
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Sorts by start and merges overlapping or adjacent intervals in place.
 *
 * Merged `label` is the maximum of merged components (preserves noisy-region size metadata).
 * @param intervals In/out interval list.
 */
void merge_intervals(std::vector<Interval>& intervals) {
    if (intervals.empty()) return;
    std::sort(intervals.begin(), intervals.end(), [](const Interval& lhs, const Interval& rhs) {
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        return lhs.end < rhs.end;
    });
    size_t write_i = 0;
    for (size_t read_i = 1; read_i < intervals.size(); ++read_i) {
        if (intervals[read_i].beg <= intervals[write_i].end + 1) {
            intervals[write_i].end = std::max(intervals[write_i].end, intervals[read_i].end);
            intervals[write_i].label = std::max(intervals[write_i].label, intervals[read_i].label);
        } else {
            intervals[++write_i] = intervals[read_i];
        }
    }
    intervals.resize(write_i + 1);
}

/**
 * @brief Builds an unindexed `cgranges_t` from 1-based inclusive intervals (stored 0-based half-open).
 *
 * Returns nullptr if empty or all intervals have `beg > end`. Caller must `cr_index` before overlap queries.
 * @param intervals 1-based inclusive source intervals.
 * @return New `cgranges_t*` or nullptr.
 */
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals) {
    if (intervals.empty()) return nullptr;
    cgranges_t* cr = cr_init();
    int added = 0;
    for (const Interval& iv : intervals) {
        if (iv.beg > iv.end) continue;
        const int32_t st = static_cast<int32_t>(iv.beg - 1);
        const int32_t en = static_cast<int32_t>(iv.end);
        const int32_t label = iv.label > 0 ? iv.label : 1;
        cr_add(cr, "cr", st, en, label);
        ++added;
    }
    if (added == 0) {
        cr_destroy(cr);
        return nullptr;
    }
    // Leave unindexed: callers run cr_index once before overlap/merge.
    return cr;
}

/**
 * @brief Reads intervals from the synthetic `"cr"` contig into \a out (inverse of `intervals_to_cr`).
 * @param cr Indexed or unindexed tree using contig `"cr"`.
 * @param out Destination vector (cleared first).
 */
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out) {
    out.clear();
    if (cr == nullptr || cr->n_r == 0) return;
    const int32_t ctg = cr_get_ctg(cr, "cr");
    if (ctg < 0) return;
    const cr_ctg_t* c = &cr->ctg[ctg];
    for (int64_t i = c->off; i < c->off + c->n; ++i) {
        const hts_pos_t beg = static_cast<hts_pos_t>(cr_start(cr, i) + 1);
        const hts_pos_t end = static_cast<hts_pos_t>(cr_end(cr, i));
        out.push_back(Interval{beg, end, cr_label(cr, i)});
    }
}

// longcallD `post_process_noisy_regs`: final `cr_add(noisy_regs, "cr", noisy_reg_start[reg_i],
// noisy_reg_end[reg_i], ...)` stores `st == noisy_reg_start` (1-based genomic start), not `start-1`.
// Use these when converting that tree to/from `Interval` so `sort_noisy_regs` lens match
// `cr_end - cr_start` and step-4 bounds match `collect_noisy_vars1`.

static void intervals_from_cr_lcd_chunk_noisy_post_merge(const cgranges_t* cr, std::vector<Interval>& out) {
    out.clear();
    if (cr == nullptr || cr->n_r == 0) return;
    const int32_t ctg = cr_get_ctg(cr, "cr");
    if (ctg < 0) return;
    const cr_ctg_t* c = &cr->ctg[ctg];
    for (int64_t i = c->off; i < c->off + c->n; ++i) {
        const hts_pos_t beg = static_cast<hts_pos_t>(cr_start(cr, i));
        const hts_pos_t end = static_cast<hts_pos_t>(cr_end(cr, i));
        out.push_back(Interval{beg, end, cr_label(cr, i)});
    }
}

static cgranges_t* intervals_to_cr_lcd_chunk_noisy_post_merge(const std::vector<Interval>& intervals) {
    if (intervals.empty()) return nullptr;
    cgranges_t* cr = cr_init();
    int added = 0;
    for (const Interval& iv : intervals) {
        if (iv.beg > iv.end) continue;
        const int32_t st = static_cast<int32_t>(iv.beg);
        const int32_t en = static_cast<int32_t>(iv.end);
        const int32_t label = iv.label > 0 ? iv.label : 1;
        cr_add(cr, "cr", st, en, label);
        ++added;
    }
    if (added == 0) {
        cr_destroy(cr);
        return nullptr;
    }
    return cr;
}

/**
 * @brief Indexed `cgranges_t` for one read's `noisy_regions` (overlap tests in noisy preprocessing).
 * @param read Read whose `noisy_regions` may be empty.
 * @return Indexed tree or nullptr.
 */
cgranges_t* build_read_noisy_cr(const ReadRecord& read) {
    if (read.noisy_regions.empty()) return nullptr;
    cgranges_t* cr = intervals_to_cr(read.noisy_regions);
    if (cr != nullptr) cr_index(cr);
    return cr;
}

/**
 * @brief Reference span for overlap tests (1-based inclusive).
 *
 * Insertions yield `var_end == key.pos - 1` (zero-length span) per longcallD.
 * @param key Variant locus.
 * @param var_start Output first reference coordinate.
 * @param var_end Output last reference coordinate.
 */
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end) {
    if (key.type == VariantType::Snp) {
        var_start = var_end = key.pos;
        return;
    }
    if (key.type == VariantType::Deletion) {
        var_start = key.pos;
        var_end = key.pos + key.ref_len - 1;
        return;
    }
    // Insertions: longcallD uses [pos, pos] for overlap/noisy-ratio checks.
    var_start = key.pos;
    var_end = key.pos;
}

// ════════════════════════════════════════════════════════════════════════════
// Noisy region processing
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Expands \a start/\a end to the union with any overlapping low-complexity intervals.
 * @param low_comp_cr Indexed low-complexity tree or nullptr.
 * @param start Reference start (1-based, converted internally).
 * @param end Reference end.
 * @param new_start Output expanded start.
 * @param new_end Output expanded end.
 */
static void low_comp_cr_start_end(cgranges_t* low_comp_cr,
                                  int32_t start, int32_t end,
                                  int32_t* new_start, int32_t* new_end) {
    *new_start = start;
    *new_end = end;
    if (low_comp_cr == nullptr || low_comp_cr->n_r == 0) return;
    int64_t* low_comp_b = nullptr;
    int64_t max_low_comp_n = 0;
    const int64_t low_comp_n =
        cr_overlap(low_comp_cr, "cr", start - 1, end, &low_comp_b, &max_low_comp_n);
    for (int64_t j = 0; j < low_comp_n; ++j) {
        const int32_t s = cr_start(low_comp_cr, low_comp_b[j]) + 1;
        const int32_t e = cr_end(low_comp_cr, low_comp_b[j]);
        if (s < *new_start) *new_start = s;
        if (e > *new_end) *new_end = e;
    }
    std::free(low_comp_b);
}

/**
 * @brief Extends noisy regions into sdust intervals, then `cr_merge` (longcallD `cr_extend_noisy_regs_with_low_comp`).
 * @param chunk_noisy_regs Unindexed or indexed noisy tree (consumed/replaced when low-comp extends).
 * @param low_comp_cr Low-complexity intervals or nullptr.
 * @param merge_dis Merge distance for `cr_merge`.
 * @param min_sv_len `min_sv_len` argument to `cr_merge`.
 * @return Merged `cgranges_t*` (caller owns prior tree if replaced).
 */
static cgranges_t* cr_extend_noisy_regs_with_low_comp(cgranges_t* chunk_noisy_regs,
                                                      cgranges_t* low_comp_cr,
                                                      int merge_dis,
                                                      int min_sv_len) {
    if (chunk_noisy_regs == nullptr) return nullptr;
    cgranges_t* cur = chunk_noisy_regs;
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        cgranges_t* new_noisy_regs = cr_init();
        for (int i = 0; i < cur->n_r; ++i) {
            const int32_t start = cr_start(cur, i) + 1;
            const int32_t end = cr_end(cur, i);
            int32_t new_start = start;
            int32_t new_end = end;
            low_comp_cr_start_end(low_comp_cr, start, end, &new_start, &new_end);
            cr_add(new_noisy_regs, "cr", new_start - 1, new_end, cr_label(cur, i));
        }
        cr_index(new_noisy_regs);
        cr_destroy(cur);
        cur = new_noisy_regs;
    }
    return cr_merge(cur, -1, merge_dis, min_sv_len);
}

/**
 * @brief Categories excluded when flanking noisy regions with nearby clean variants.
 * @param c Final or intermediate category.
 */
static bool category_skipped_for_noisy_flank(VariantCategory c) {
    return c == VariantCategory::LowCoverage || c == VariantCategory::StrandBias ||
           c == VariantCategory::NonVariant;
}

/**
 * @brief Per-noisy-region flanking bounds including nearest clean candidates (longcallD `collect_noisy_reg_start_end`).
 * @param chunk_noisy_regs Indexed noisy regions.
 * @param cand Candidates (sorted as in chunk).
 * @param cand_cate Parallel category per candidate (same length as `cand`).
 * @param start_out Output left bounds per noisy region.
 * @param end_out Output right bounds per noisy region.
 */
static void collect_noisy_reg_start_end_pgphase(const cgranges_t* chunk_noisy_regs,
                                                const CandidateTable& cand,
                                                const std::vector<VariantCategory>& cand_cate,
                                                std::vector<int32_t>& start_out,
                                                std::vector<int32_t>& end_out) {
    const auto flank_span = [](const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end) {
        // Match longcallD collect_noisy_reg_start_end:
        // var_start = var->pos; var_end = var->pos + var->ref_len - 1.
        var_start = key.pos;
        var_end = key.pos + key.ref_len - 1;
    };

    const int n_noisy = chunk_noisy_regs->n_r;
    start_out.assign(n_noisy, 0);
    end_out.assign(n_noisy, 0);
    if (cand.empty()) {
        for (int reg_i = 0; reg_i < n_noisy; ++reg_i) {
            const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
            const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
            start_out[static_cast<size_t>(reg_i)] = ori_reg_start - kNoisyRegFlankLen;
            end_out[static_cast<size_t>(reg_i)] = ori_reg_end + kNoisyRegFlankLen;
        }
        return;
    }
    std::vector<int> max_left_var_i(n_noisy, -1);
    std::vector<int> min_right_var_i(n_noisy, -1);

    int reg_i = 0;
    int var_i = 0;
    while (reg_i < n_noisy && var_i < static_cast<int>(cand.size())) {
        if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(var_i)])) {
            ++var_i;
            continue;
        }
        hts_pos_t var_start = 0;
        hts_pos_t var_end = 0;
        flank_span(cand[static_cast<size_t>(var_i)].key, var_start, var_end);
        const int32_t reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t reg_end = cr_end(chunk_noisy_regs, reg_i);
        if (static_cast<int32_t>(var_start) > reg_end) {
            if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
                min_right_var_i[static_cast<size_t>(reg_i)] = var_i;
            }
            ++reg_i;
        } else if (static_cast<int32_t>(var_end) < reg_start) {
            max_left_var_i[static_cast<size_t>(reg_i)] = var_i;
            ++var_i;
        } else {
            ++var_i;
        }
    }

    for (reg_i = 0; reg_i < n_noisy; ++reg_i) {
        if (max_left_var_i[static_cast<size_t>(reg_i)] == -1) {
            max_left_var_i[static_cast<size_t>(reg_i)] =
                std::min(std::max(0, static_cast<int>(cand.size()) - 1), 0);
        }
        if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
            min_right_var_i[static_cast<size_t>(reg_i)] =
                std::max(0, static_cast<int>(cand.size()) - 1);
        }
        const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
        int32_t cur_start = ori_reg_start - kNoisyRegFlankLen;
        for (int v = max_left_var_i[static_cast<size_t>(reg_i)]; v >= 0; --v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            flank_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(ve) < cur_start - 1) break;
            if (static_cast<int32_t>(vs) - kNoisyRegFlankLen < cur_start) {
                cur_start = static_cast<int32_t>(vs) - kNoisyRegFlankLen;
            }
        }
        start_out[static_cast<size_t>(reg_i)] = cur_start;
        int32_t cur_end = ori_reg_end + kNoisyRegFlankLen;
        for (int v = min_right_var_i[static_cast<size_t>(reg_i)]; v < static_cast<int>(cand.size()); ++v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            flank_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(vs) > cur_end + 1) break;
            if (static_cast<int32_t>(ve) + kNoisyRegFlankLen > cur_end) {
                cur_end = static_cast<int32_t>(ve) + kNoisyRegFlankLen;
            }
        }
        end_out[static_cast<size_t>(reg_i)] = cur_end;
    }
}

/**
 * @brief Merges per-read noisy intervals into unified chunk-level regions and extends low-comp windows.
 *
 * Implements `longcallD`'s `pre_process_noisy_regs()` from `collect_var.c`.
 * It sweeps across individually flagged noisy segments (dense mismatches),
 * merges overlapping intervals from multiple reads via `cgranges_t` merges,
 * extends them to encompass neighboring low-complexity sequence (like homopolymers),
 * and finally drops regions that lack sufficient read support (min_alt_depth or min_af)
 * to be deemed true structural variations rather than isolated noise bursts.
 *
 * Call order matches longcallD `collect_var_main`: after deduplicated sites and allele-depth
 * collection, immediately before `classify_cand_vars`.
 *
 * @param chunk Collection of read records and current candidates in this block.
 * @param opts User-defined filtering thresholds.
 */
void pre_process_noisy_regs_pgphase(BamChunk& chunk, const Options& opts) {
    if (chunk.noisy_regions.empty()) return;

    CrangesOwner noisy_own;
    noisy_own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (noisy_own.cr == nullptr) {
        chunk.noisy_regions.clear();
        return;
    }
    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(chunk.low_complexity_regions));

    cr_index(noisy_own.cr);
    if (low_own.cr != nullptr) cr_index(low_own.cr);

    const int merge_dis = opts.noisy_reg_merge_dis;
    const int msv = opts.min_sv_len;

    // longcallD: cr_extend_noisy_regs_with_low_comp then a second cr_merge.
    cgranges_t* noisy = noisy_own.release();
    noisy = cr_extend_noisy_regs_with_low_comp(noisy, low_own.cr, merge_dis, msv);
    noisy = cr_merge(noisy, -1, merge_dis, msv);
    noisy_own.adopt(noisy);

    cgranges_t* noisy_regs = noisy_own.cr;
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    std::vector<uint8_t> skip_noisy_reg(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_total(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_noisy(static_cast<size_t>(noisy_regs->n_r), 0);

    const std::vector<int>& read_order = chunk.ordered_read_ids;
    const auto visit_read = [&](const ReadRecord& read) {
        if (read.is_skipped) return;
        const int64_t beg = read.beg;
        const int64_t end = read.end;
        const int64_t ovlp_n = cr_overlap(
            noisy_regs, "cr", static_cast<int32_t>(beg - 1), static_cast<int32_t>(end), &ovlp_b, &max_b);
        CrangesOwner read_noisy_own;
        read_noisy_own.adopt(build_read_noisy_cr(read));
        int64_t* read_ovlp_b = nullptr;
        int64_t read_max_b = 0;
        for (int64_t ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            const int r_i = static_cast<int>(ovlp_b[ovlp_i]);
            noisy_reg_to_total[static_cast<size_t>(r_i)]++;
            const int32_t noisy_reg_start = cr_start(noisy_regs, r_i) + 1;
            const int32_t noisy_reg_end = cr_end(noisy_regs, r_i);
            int64_t noisy_digar_ovlp_n = 0;
            if (read_noisy_own.cr != nullptr) {
                noisy_digar_ovlp_n =
                    cr_overlap(read_noisy_own.cr, "cr", noisy_reg_start - 1, noisy_reg_end,
                               &read_ovlp_b, &read_max_b);
            }
            if (noisy_digar_ovlp_n > 0) {
                noisy_reg_to_noisy[static_cast<size_t>(r_i)]++;
            }
        }
        std::free(read_ovlp_b);
    };
    if (read_order.empty() && !chunk.reads.empty()) {
        for (const ReadRecord& r : chunk.reads) visit_read(r);
    } else {
        for (int ord : read_order) {
            if (ord < 0 || static_cast<size_t>(ord) >= chunk.reads.size()) continue;
            visit_read(chunk.reads[static_cast<size_t>(ord)]);
        }
    }
    std::free(ovlp_b);

    const int min_noisy_reg_reads = opts.min_alt_depth;
    const float min_noisy_reg_ratio = static_cast<float>(opts.min_af);
    const int min_overlap_reads = opts.min_noisy_reg_total_depth;
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        const int n_noisy = noisy_reg_to_noisy[static_cast<size_t>(i)];
        const int n_total = noisy_reg_to_total[static_cast<size_t>(i)];
        if ((min_overlap_reads > 0 && n_total < min_overlap_reads) ||
            n_noisy < min_noisy_reg_reads ||
            (n_total > 0 &&
             static_cast<float>(n_noisy) / static_cast<float>(n_total) < min_noisy_reg_ratio)) {
            skip_noisy_reg[static_cast<size_t>(i)] = 1;
            ++n_skipped;
        }
    }

    if (n_skipped > 0) {
        cgranges_t* new_noisy_regs = cr_init();
        for (int i = 0; i < noisy_regs->n_r; ++i) {
            if (skip_noisy_reg[static_cast<size_t>(i)]) continue;
            cr_add(new_noisy_regs, "cr", cr_start(noisy_regs, i), cr_end(noisy_regs, i),
                   cr_label(noisy_regs, i));
        }
        cr_index(new_noisy_regs);
        noisy_own.adopt(new_noisy_regs);
    }

    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

/**
 * @brief Widens noisy intervals using classified candidates, then re-merges (longcallD `post_process_noisy_regs`).
 * @param chunk Chunk whose `noisy_regions` are rewritten.
 * @param cand Classified candidates (categories drive flanking).
 */
void post_process_noisy_regs_pgphase(BamChunk& chunk, const CandidateTable& cand) {
    if (chunk.noisy_regions.empty()) return;
    CrangesOwner noisy_own;
    noisy_own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (noisy_own.cr == nullptr) {
        chunk.noisy_regions.clear();
        return;
    }
    cr_index(noisy_own.cr);
    cgranges_t* chunk_noisy_regs = noisy_own.release();
    const int n_noisy_regs = chunk_noisy_regs->n_r;
    std::vector<int32_t> noisy_reg_start(static_cast<size_t>(n_noisy_regs));
    std::vector<int32_t> noisy_reg_end(static_cast<size_t>(n_noisy_regs));
    std::vector<VariantCategory> cand_cate;
    cand_cate.reserve(cand.size());
    for (const auto& cv : cand) cand_cate.push_back(cv.counts.category);
    collect_noisy_reg_start_end_pgphase(
        chunk_noisy_regs, cand, cand_cate, noisy_reg_start, noisy_reg_end);

    cgranges_t* noisy_regs = cr_init();
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        cr_add(noisy_regs, "cr",
               noisy_reg_start[static_cast<size_t>(reg_i)],
               noisy_reg_end[static_cast<size_t>(reg_i)],
               cr_label(chunk_noisy_regs, reg_i));
    }
    cr_index(noisy_regs);
    cr_destroy(chunk_noisy_regs);
    cgranges_t* merged = cr_merge(noisy_regs, 0, -1, -1);
    noisy_own.adopt(merged);
    intervals_from_cr_lcd_chunk_noisy_post_merge(noisy_own.cr, chunk.noisy_regions);
}

/**
 * @brief Traverses surviving variants and drops those entirely swallowed by a noisy block.
 *
 * Implements `longcallD`'s final masking pass equivalent (`cr_is_contained`).
 * Like classify compaction, skips sites whose bitmask intersects
 * `kLongcalldNotCandVarCate` (`LONGCALLD_NOT_CAND_VAR_CATE`).
 * Once noisy intervals are merged and widened via `post_process_noisy_regs_pgphase()`, this
 * marks surviving germline candidates fully contained in noisy spans as `NON_VAR`.
 *
 * This ensures that a single structural variance window gets reported cleanly
 * instead of hundreds of raw SNVs generated through misalignment.
 *
 * @param chunk Chunk tracking the current region's evaluated candidates.
 * @note No-op when `chunk.noisy_regions` is empty.
 */
void apply_noisy_containment_filter(BamChunk& chunk) {
    if (chunk.noisy_regions.empty()) return;
    CrangesOwner own;
    own.adopt(intervals_to_cr_lcd_chunk_noisy_post_merge(chunk.noisy_regions));
    if (own.cr == nullptr || own.cr->n_r == 0) return;
    cr_index(own.cr);
    int64_t* b = nullptr;
    int64_t m = 0;
    for (auto& cv : chunk.candidates) {
        if ((cv.lcd_var_i_to_cate & kLongcalldNotCandVarCate) != 0) continue;
        const VariantKey& k = cv.key;
        int64_t n = 0;
        if (k.type == VariantType::Insertion) {
            n = cr_is_contained(
                own.cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), &b, &m);
        } else {
            n = cr_is_contained(own.cr, "cr",
                                static_cast<int32_t>(k.pos - 1),
                                static_cast<int32_t>(k.pos + k.ref_len),
                                &b, &m);
        }
        if (n > 0) {
            cv.counts.category = VariantCategory::NonVariant;
            cv.lcd_var_i_to_cate = category_to_flag(VariantCategory::NonVariant);
        }
    }
    std::free(b);
}

static bool is_not_candidate_category(VariantCategory c) {
    return c == VariantCategory::NonVariant ||
           c == VariantCategory::LowCoverage ||
           c == VariantCategory::StrandBias;
}

static void prune_not_candidate_variants(BamChunk& chunk) {
    CandidateTable kept;
    kept.reserve(chunk.candidates.size());
    for (CandidateVariant& cv : chunk.candidates) {
        if (is_not_candidate_category(cv.counts.category)) continue;
        kept.push_back(std::move(cv));
    }
    chunk.candidates = std::move(kept);
}

/**
 * @brief Fills `chunk.low_complexity_regions` via sdust on the loaded reference slice.
 * @param chunk Must have `ref_seq` and region/ref bounds aligned with the slice.
 */
void populate_low_complexity_intervals(BamChunk& chunk) {
    chunk.low_complexity_regions.clear();
    if (chunk.ref_seq.empty()) return;

    const hts_pos_t beg = std::max(chunk.region.beg, chunk.ref_beg);
    const hts_pos_t end = std::min(chunk.region.end, chunk.ref_end);
    if (beg > end) return;

    const size_t offset = static_cast<size_t>(beg - chunk.ref_beg);
    const int len = static_cast<int>(end - beg + 1);
    int n = 0;
    uint64_t* intervals = sdust(nullptr,
                                reinterpret_cast<const uint8_t*>(chunk.ref_seq.data() + offset),
                                len, kSdustThreshold, kSdustWindow, &n);
    for (int i = 0; i < n; ++i) {
        const hts_pos_t rel_beg = static_cast<hts_pos_t>(intervals[i] >> 32);
        const hts_pos_t rel_end = static_cast<hts_pos_t>(static_cast<uint32_t>(intervals[i]));
        if (rel_end <= rel_beg) continue;
        chunk.low_complexity_regions.push_back(
            Interval{beg + rel_beg, beg + rel_end - 1, static_cast<int>(rel_end - rel_beg)});
    }
    std::free(intervals);
}

/**
 * @brief Sets `ordered_read_ids` and unions per-read `noisy_regions` into chunk-level intervals.
 * @param chunk Chunk with `reads` populated.
 */
void populate_chunk_read_indexes(BamChunk& chunk) {
    chunk.ordered_read_ids.clear();
    chunk.noisy_regions.clear();
    chunk.ordered_read_ids.reserve(chunk.reads.size());
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        chunk.ordered_read_ids.push_back(static_cast<int>(read_i));
        if (read.is_skipped) continue;
        for (const Interval& iv : read.noisy_regions) {
            if (iv.end < chunk.region.beg || iv.beg > chunk.region.end) continue;
            chunk.noisy_regions.push_back(iv);
        }
    }
    // longcallD uses qsort(comp_bam_read_sort): pos asc, end desc, NM asc, qname asc.
    // Use std::sort (non-stable) to mirror qsort tie behavior as closely as possible.
    std::sort(chunk.ordered_read_ids.begin(), chunk.ordered_read_ids.end(),
              [&](int lhs_i, int rhs_i) {
                  const ReadRecord& lhs = chunk.reads[static_cast<size_t>(lhs_i)];
                  const ReadRecord& rhs = chunk.reads[static_cast<size_t>(rhs_i)];
                  if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
                  if (lhs.end != rhs.end) return lhs.end > rhs.end;
                  if (lhs.nm != rhs.nm) return lhs.nm < rhs.nm;
                  return lhs.qname < rhs.qname;
              });
}

// ════════════════════════════════════════════════════════════════════════════
// Variant classification — math utilities
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Log of the hypergeometric PMF term (longcallD `log_hypergeometric`, `lgamma` only).
 * @param a,b,c,d Fourfold table entries (row1: a,b; row2: c,d).
 * @return Log probability mass; `-inf` if totals are non-positive.
 */
static double log_hypergeom_pg(int a, int b, int c, int d) {
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N = n1 + n2;
    if (N <= 0) return -std::numeric_limits<double>::infinity();
    if (n1 > n2) return log_hypergeom_pg(c, d, a, b);
    if (m1 > m2) return log_hypergeom_pg(b, a, d, c);
    return std::lgamma(static_cast<double>(n1 + 1)) + std::lgamma(static_cast<double>(n2 + 1)) +
           std::lgamma(static_cast<double>(m1 + 1)) + std::lgamma(static_cast<double>(m2 + 1)) -
           (std::lgamma(static_cast<double>(a + 1)) + std::lgamma(static_cast<double>(b + 1)) +
            std::lgamma(static_cast<double>(c + 1)) + std::lgamma(static_cast<double>(d + 1)) +
            std::lgamma(static_cast<double>(N + 1)));
}

/**
 * @brief Evaluates statistical probability of read variance within noisy loci.
 *
 * Implements Fisher's Exact test and log hypergeometric distribution (analogous logic
 * in `fisher_exact()` of longcallD) evaluating whether read strands deviate
 * significantly from uniform distribution.
 *
 * Calculates Strand Bias and noise correlations dynamically.
 *
 * @param a Alt-forward count.
 * @param b Alt-reverse count.
 * @param c,d Padding cells so `a+b+c+d` matches longcallD’s symmetric 2×2 layout.
 * @return Two-sided Fisher p-value in [0,1].
 */
static double fisher_exact_two_tail(int a, int b, int c, int d) {
    if (a + b + c + d <= 0) return 1.0;
    const double log_p_observed = log_hypergeom_pg(a, b, c, d);
    const double p_observed = std::exp(log_p_observed);
    double total_p = 0.0;
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    const int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    const int denom = a + b + c + d;
    const int mode_a =
        denom > 0
            ? static_cast<int>((static_cast<double>(a + b) * static_cast<double>(a + c)) /
                               static_cast<double>(denom))
            : 0;
    for (int delta = 0; delta <= max_a - min_a; ++delta) {
        int current_a = mode_a + delta;
        if (current_a <= max_a) {
            int current_b = (a + b) - current_a;
            int current_c = (a + c) - current_a;
            int current_d = (b + d) - current_b;
            if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                const double p = std::exp(log_hypergeom_pg(current_a, current_b, current_c, current_d));
                if (p <= p_observed + DBL_EPSILON) total_p += p;
            }
        }
        if (delta > 0) {
            current_a = mode_a - delta;
            if (current_a >= min_a) {
                int current_b = (a + b) - current_a;
                int current_c = (a + c) - current_a;
                int current_d = (b + d) - current_b;
                if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                    const double p =
                        std::exp(log_hypergeom_pg(current_a, current_b, current_c, current_d));
                    if (p <= p_observed + DBL_EPSILON) total_p += p;
                }
            }
        }
    }
    return total_p;
}

/**
 * @brief Maps ASCII base to 0–3 (ACGT) or 4 for other.
 */
static int nt4_from_ref_char(char ch) {
    switch (std::toupper(static_cast<unsigned char>(ch))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}

/**
 * @brief Encoded base at absolute `abs_pos` within a loaded reference substring starting at `ref_beg`.
 */
static int ref_nt4_at(const std::string& ref_seq, hts_pos_t ref_beg, hts_pos_t abs_pos) {
    if (ref_seq.empty() || abs_pos < ref_beg) return 4;
    const hts_pos_t rel = abs_pos - ref_beg;
    if (rel < 0 || static_cast<size_t>(rel) >= ref_seq.size()) return 4;
    return nt4_from_ref_char(ref_seq[static_cast<size_t>(rel)]);
}

/**
 * @brief True if a short indel sits inside or immediately beside a homopolymer or tandem repeat run.
 *
 * Only applied to indels with span ≤ `xid` (default 5 bp from `noisy_reg_max_xgaps`); larger
 * indels return false immediately. SNPs are also gated: `ref_len` must be ≤ `xid`.
 *
 * **Algorithm:** Checks up to 6 bp **downstream** of the variant for a repeat run, then up to
 * 6 bp **upstream**. For each direction, tests repeat unit lengths 1–6 bp; a unit qualifies if
 * at least **3 consecutive copies** are present in the reference. The first matching unit in
 * either direction makes the variant a homopolymer candidate and returns `true`.
 *
 * Example — homopolymer (unit = 1 bp):
 * ```
 * reference: ...C A A A A A A A G...
 *                  ^ deletion of one A
 * downstream: AAAAAA — 6 copies of "A" → homopolymer
 * ```
 *
 * Example — short tandem repeat (unit = 2 bp):
 * ```
 * reference: ...G C A C A C A C A T...
 *                    ^ deletion of "CA"
 * downstream: CACACA — 3 copies of "CA" → tandem repeat
 * ```
 *
 * @param var Variant to test (type must be Insertion or Deletion; SNPs are gated by xid).
 * @param ref_seq Reference substring for the chunk.
 * @param ref_beg Absolute 1-based coordinate of `ref_seq[0]`.
 * @param ref_end Unused (signature parity with longcallD).
 * @param xid Maximum indel span to check; variants larger than this return false.
 */
static bool var_is_homopolymer_pg(const VariantKey& var,
                                  const std::string& ref_seq,
                                  hts_pos_t ref_beg,
                                  hts_pos_t ref_end,
                                  int xid) {
    if (ref_seq.empty()) return false;
    hts_pos_t start_pos = 0;
    hts_pos_t end_pos = 0;
    if (var.type == VariantType::Snp) {
        start_pos = var.pos - 1;
        end_pos = var.pos + 1;
    } else if (var.type == VariantType::Insertion) {
        const int alt_len = static_cast<int>(var.alt.size());
        if (alt_len > xid) return false;
        start_pos = var.pos - 1;
        end_pos = var.pos;
    } else {
        if (var.ref_len > xid) return false;
        start_pos = var.pos + var.ref_len - 1;
        end_pos = var.pos;
    }
    (void)ref_end;
    int is_homopolymer = 1;
    constexpr int max_unit_len = 6;
    constexpr int n_check_copy_num = 3;
    uint8_t ref_bases[6];
    for (int i = 0; i < 6; ++i) {
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, end_pos + i));
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, end_pos + i * rep_unit_len + j) != ref_bases[j]) {
                    is_homopolymer = 0;
                    break;
                }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) break;
    }
    if (is_homopolymer) return true;
    for (int i = 0; i < 6; ++i) {
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, start_pos - i));
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, start_pos - i * rep_unit_len - j) != ref_bases[j]) {
                    is_homopolymer = 0;
                    break;
                }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) break;
    }
    return is_homopolymer != 0;
}

/**
 * @brief True if the reference flanking an indel is a tandem repeat of the deleted/inserted motif.
 *
 * Only applied to indels with span ≤ `xid` (default 5 bp); larger indels return false immediately.
 * Complementary to `var_is_homopolymer_pg`: catches micro-homology / STR contexts where the
 * alignment can be shifted left or right with equal edit cost, making the exact breakpoint ambiguous.
 *
 * **Algorithm:**
 * - **Deletion** of length L: compares the deleted sequence (fetched from `ref_seq` at `pos`)
 *   against the `3 × L` bases immediately following the deletion endpoint. If they match the
 *   deleted motif repeated 3 times, returns `true`.
 * - **Insertion** of length L: constructs a 3-copy repeat of the inserted sequence and compares
 *   it against `3 × L` bases of the reference starting at `pos`. Returns `true` on match.
 *
 * Example — deletion in a tandem repeat:
 * ```
 * delete "AT" at pos 100 (reference: ...ATATAT...)
 * check ref[102..107] = "ATATAT" == "AT"×3 → true
 * ```
 *
 * Example — insertion consistent with STR:
 * ```
 * insert "GC" at pos 200 (reference: ...GCGCGC...)
 * check ref[200..205] = "GCGCGC" == "GC"×3 → true
 * ```
 *
 * @param var Variant to test (must be Insertion or Deletion).
 * @param ref_seq Reference substring for the chunk.
 * @param ref_beg Absolute 1-based coordinate of `ref_seq[0]`.
 * @param ref_end Last valid absolute coordinate in `ref_seq`.
 * @param xid Maximum indel span to check; variants larger than this return false.
 * @return True if the flanking reference matches 3 tandem copies of the indel motif.
 */
static bool var_is_repeat_region_pg(const VariantKey& var,
                                    const std::string& ref_seq,
                                    hts_pos_t ref_beg,
                                    hts_pos_t ref_end,
                                    int xid) {
    if (ref_seq.empty()) return false;
    const hts_pos_t pos = var.pos;
    if (var.type == VariantType::Deletion) {
        const int del_len = var.ref_len;
        if (del_len > xid) return false;
        const int len = del_len * 3;
        if (pos < ref_beg || pos + del_len + len > ref_end) return false;
        const size_t off = static_cast<size_t>(pos - ref_beg);
        const size_t off2 = static_cast<size_t>(pos + del_len - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size() ||
            off2 + static_cast<size_t>(len) > ref_seq.size()) return false;
        return std::memcmp(ref_seq.data() + off, ref_seq.data() + off2, static_cast<size_t>(len)) == 0;
    }
    if (var.type == VariantType::Insertion) {
        const int ins_len = static_cast<int>(var.alt.size());
        if (ins_len > xid) return false;
        const int len = ins_len * 3;
        if (pos < ref_beg || pos + len > ref_end) return false;
        const size_t off = static_cast<size_t>(pos - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size()) return false;
        std::string ref_b = ref_seq.substr(off, static_cast<size_t>(len));
        std::string alt_b = ref_b;
        for (int j = ins_len; j < len; ++j) {
            alt_b[static_cast<size_t>(j)] = alt_b[static_cast<size_t>(j - ins_len)];
        }
        for (int j = 0; j < ins_len; ++j) {
            alt_b[static_cast<size_t>(j)] = var.alt[static_cast<size_t>(j)];
        }
        return ref_b == alt_b;
    }
    return false;
}

// ════════════════════════════════════════════════════════════════════════════
// Noisy reads ratio cache  (longcallD `build_var_noisy_reads_ratio_cache`)
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Lazily builds read coverage and merged-XID interval trees for `var_noisy_reads_ratio`.
 * @param chunk Chunk whose `var_noisy_read_*` fields are populated once per classification pass.
 */
static void build_var_noisy_read_cache(BamChunk& chunk) {
    if (chunk.var_noisy_read_cov_cr) return;

    chunk.var_noisy_read_cov_cr.reset(cr_init());
    chunk.var_noisy_read_err_cr.reset(cr_init());
    chunk.var_noisy_read_marks.assign(chunk.reads.size(), 0);
    chunk.var_noisy_read_mark_id = 0;

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& r = chunk.reads[read_i];
        if (r.is_skipped || r.beg > r.end) continue;

        cr_add(chunk.var_noisy_read_cov_cr.get(), "cr",
               static_cast<int32_t>(r.beg - 1), static_cast<int32_t>(r.end),
               static_cast<int32_t>(read_i));

        // Merge adjacent XID events per read → noisy intervals.
        hts_pos_t noisy_start = -1;
        hts_pos_t noisy_end = -1;
        bool has_noisy = false;
        for (const DigarOp& d : r.digars) {
            if (d.type != DigarType::Snp && d.type != DigarType::Insertion &&
                d.type != DigarType::Deletion) continue;
            const hts_pos_t cur_start = d.pos - 1;
            hts_pos_t cur_end = d.pos;
            if (d.type == DigarType::Snp || d.type == DigarType::Deletion) cur_end += d.len - 1;
            if (!has_noisy) {
                noisy_start = cur_start;
                noisy_end = cur_end;
                has_noisy = true;
                continue;
            }
            if (cur_start < noisy_end) {
                if (cur_end > noisy_end) noisy_end = cur_end;
                continue;
            }
            cr_add(chunk.var_noisy_read_err_cr.get(), "cr",
                   static_cast<int32_t>(noisy_start), static_cast<int32_t>(noisy_end),
                   static_cast<int32_t>(read_i));
            noisy_start = cur_start;
            noisy_end = cur_end;
        }
        if (has_noisy) {
            cr_add(chunk.var_noisy_read_err_cr.get(), "cr",
                   static_cast<int32_t>(noisy_start), static_cast<int32_t>(noisy_end),
                   static_cast<int32_t>(read_i));
        }
    }
    cr_index(chunk.var_noisy_read_cov_cr.get());
    cr_index(chunk.var_noisy_read_err_cr.get());
}

/**
 * @brief Fraction of overlapping reads that carry a noisy digar span inside `[var_start,var_end]`.
 * @param chunk Chunk with reads (cache built on first use).
 * @param var_start Variant span start (1-based).
 * @param var_end Variant span end (1-based).
 * @return Ratio in [0,1], or 0 if no overlapping reads.
 */
static double var_noisy_reads_ratio(BamChunk& chunk, hts_pos_t var_start, hts_pos_t var_end) {
    build_var_noisy_read_cache(chunk);

    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    const int64_t total_n = cr_overlap(chunk.var_noisy_read_cov_cr.get(), "cr",
                                       static_cast<int32_t>(var_start - 1),
                                       static_cast<int32_t>(var_end), &ovlp_b, &max_b);
    const int total_reads = static_cast<int>(total_n);
    int noisy_reads = 0;

    if (total_reads > 0) {
        ++chunk.var_noisy_read_mark_id;
        if (chunk.var_noisy_read_mark_id == INT_MAX) {
            std::fill(chunk.var_noisy_read_marks.begin(), chunk.var_noisy_read_marks.end(), 0);
            chunk.var_noisy_read_mark_id = 1;
        }
        int64_t* noisy_b = nullptr;
        int64_t noisy_max_b = 0;
        const int64_t noisy_n = cr_overlap(chunk.var_noisy_read_err_cr.get(), "cr",
                                           static_cast<int32_t>(var_start - 1),
                                           static_cast<int32_t>(var_end), &noisy_b, &noisy_max_b);
        for (int64_t i = 0; i < noisy_n; ++i) {
            const int ri =
                static_cast<int>(cr_label(chunk.var_noisy_read_err_cr.get(), noisy_b[i]));
            if (ri < 0 || static_cast<size_t>(ri) >= chunk.var_noisy_read_marks.size()) continue;
            if (chunk.var_noisy_read_marks[static_cast<size_t>(ri)] == chunk.var_noisy_read_mark_id)
                continue;
            chunk.var_noisy_read_marks[static_cast<size_t>(ri)] = chunk.var_noisy_read_mark_id;
            ++noisy_reads;
        }
        std::free(noisy_b);
    }
    std::free(ovlp_b);
    return total_reads > 0 ? static_cast<double>(noisy_reads) / static_cast<double>(total_reads) : 0.0;
}

/**
 * @brief Adds a variant’s reference span (expanded by low-complexity overlap) to `var_cr` when filters pass.
 * @param var_cr Accumulator for dense-variant / repeat-het intervals.
 * @param low_comp_cr Low-complexity tree or nullptr.
 * @param var Locus to add.
 * @param chunk Chunk for noisy-read ratio when \a check_noisy_reads_ratio is true.
 * @param check_noisy_reads_ratio If true, requires `var_noisy_reads_ratio` ≥ `opts.min_af`.
 * @param opts Thresholds (`min_af`).
 */
static void cr_add_var_to_noisy_cr(cgranges_t* var_cr,
                                   cgranges_t* low_comp_cr,
                                   const VariantKey& var,
                                   BamChunk& chunk,
                                   bool check_noisy_reads_ratio,
                                   const Options& opts) {
    hts_pos_t var_start = 0;
    hts_pos_t var_end = 0;
    variant_genomic_span(var, var_start, var_end);
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        int64_t* low_comp_b = nullptr;
        int64_t max_low_comp_n = 0;
        const int64_t low_comp_n = cr_overlap(low_comp_cr, "cr",
                                              static_cast<int32_t>(var_start - 1),
                                              static_cast<int32_t>(var_end),
                                              &low_comp_b, &max_low_comp_n);
        for (int64_t j = 0; j < low_comp_n; ++j) {
            const int32_t start = cr_start(low_comp_cr, low_comp_b[j]) + 1;
            const int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
            if (start < var_start) var_start = start;
            if (end > var_end) var_end = end;
        }
        std::free(low_comp_b);
    }
    if (!check_noisy_reads_ratio || var_noisy_reads_ratio(chunk, var_start, var_end) >= opts.min_af) {
        cr_add(var_cr, "cr",
               static_cast<int32_t>(var_start - 1), static_cast<int32_t>(var_end), 1);
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Variant classification
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief First-pass category from depth, ONT strand bias, AF band, and repeat context (longcallD `classify_var_cate`).
 *
 * Sets `counts.allele_fraction`. ONT runs Fisher exact on forward vs reverse alt counts.
 *
 * @param key Variant identity.
 * @param counts Allele and strand tallies (in/out).
 * @param ref_slice Loaded reference sequence for the chunk.
 * @param ref_beg Absolute coordinate of `ref_slice[0]`.
 * @param ref_end Last valid absolute coordinate for repeat/homopolymer checks.
 * @param opts Depth, AF, strand p-value, and `noisy_reg_max_xgaps`.
 * @return Initial `VariantCategory` before noisy-overlap second pass.
 */
static VariantCategory classify_variant_initial(const VariantKey& key,
                                                VariantCounts& counts,
                                                const std::string& ref_slice,
                                                hts_pos_t ref_beg,
                                                hts_pos_t ref_end,
                                                const Options& opts) {
    const int depth_with_low_quality = counts.total_cov + counts.low_qual_cov;
    counts.allele_fraction =
        counts.total_cov == 0
            ? 0.0
            : static_cast<double>(counts.alt_cov) / static_cast<double>(counts.total_cov);

    if (depth_with_low_quality < opts.min_depth) return VariantCategory::LowCoverage;
    if (counts.alt_cov < opts.min_alt_depth) return VariantCategory::LowCoverage;

    if (opts.is_ont()) {
        const int fa = counts.forward_alt;
        const int ra = counts.reverse_alt;
        const int expected = (fa + ra) / 2;
        if (expected > 0) {
            const double p = fisher_exact_two_tail(fa, ra, expected, expected);
            if (p < opts.strand_bias_pval) return VariantCategory::StrandBias;
        }
    }

    if (counts.allele_fraction < opts.min_af) return VariantCategory::LowAlleleFraction;
    if (counts.allele_fraction > opts.max_af) return VariantCategory::CleanHom;

    if ((key.type == VariantType::Insertion || key.type == VariantType::Deletion) &&
        (var_is_homopolymer_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps) ||
         var_is_repeat_region_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps))) {
        return VariantCategory::RepeatHetIndel;
    }
    return key.type == VariantType::Snp ? VariantCategory::CleanHetSnp : VariantCategory::CleanHetIndel;
}

/** longcallD-shaped `classify_cand_vars`: two passes by necessity, not two separate classifiers.
 *
 * **Pass 1** — Only *local* rules: `classify_variant_initial` fills `cats[]` and
 * `candvarcate_initial` (coverage, allele fraction, strand bias, repeat het, clean categories).
 * Every candidate that is not `LOW_COV` (and not skipped for ONT strand bias when building the
 * index) inserts its reference interval into `var_pos_cr`, which is then indexed. Density checks
 * need the full set of surviving loci in that structure before any per-site `n_ov` query.
 *
 * **Pass 2** — *Context* rules: overlap with existing `chunk_noisy` → `NON_VAR` (HiFi);
 * `RepeatHetIndel` handling adds to `noisy_var_cr`; `cr_overlap` on `var_pos_cr` with `n_ov > 1`
 * flags dense clusters for `noisy_var_cr`; then `LOW_AF` is folded to `LOW_COV` (longcallD
 * `var_i_to_cate` after this loop). Merged ranges are written back to `chunk.noisy_regions`.
 *
 * Note: longcallD `NOISY_CAND_HET` / `NOISY_CAND_HOM` (`e` / `h`) are **not** produced here; they
 * come from MSA-based recall inside noisy regions. This pass matches longcallD `classify_cand_vars`
 * (candidates remain `CLEAN_*` / `REP_HET_INDEL` until `NON_VAR` from containment, etc.).
 *
 * @param chunk Chunk with `candidates`, `ref_seq`, `low_complexity_regions`, and `noisy_regions` inputs.
 * @param opts Platform and merge options.
 */
static void classify_cand_vars_pgphase(BamChunk& chunk, const Options& opts) {
    CandidateTable& variants = chunk.candidates;
    if (variants.empty()) return;

    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(chunk.low_complexity_regions));
    cgranges_t* low_comp = low_own.cr;
    if (low_comp != nullptr) cr_index(low_comp);

    // longcallD `classify_cand_vars` does not modify cand.alt_ref_base (only init_cand_vars_based_on_sites
    // sets 4 and noisy `make_cand_vars_from_baln0` sets consensus column bytes). SNP ref_base is filled
    // here for classifier parity with existing pgphase paths.
    for (auto& cv : variants) {
        if (cv.key.type == VariantType::Snp) {
            cv.ref_base = static_cast<uint8_t>(ref_nt4_at(chunk.ref_seq, chunk.ref_beg, cv.key.pos));
        }
    }

    // Pass 1: local labels + interval index of all non–LOW_COV candidates (see function comment).
    std::vector<VariantCategory> cats;
    cats.reserve(variants.size());
    for (auto& cv : variants) {
        const VariantCategory init = classify_variant_initial(
            cv.key, cv.counts, chunk.ref_seq, chunk.ref_beg, chunk.ref_end, opts);
        cv.counts.candvarcate_initial = init;
        cats.push_back(init);
    }

    cgranges_t* var_pos_cr = cr_init();
    for (size_t i = 0; i < variants.size(); ++i) {
        const VariantCategory c = cats[i];
        if (c == VariantCategory::LowCoverage) continue;
        if (opts.is_ont() && c == VariantCategory::StrandBias) continue;
        const VariantKey& k = variants[i].key;
        if (k.type == VariantType::Insertion) {
            cr_add(var_pos_cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), 1);
        } else {
            cr_add(var_pos_cr, "cr",
                   static_cast<int32_t>(k.pos - 1),
                   static_cast<int32_t>(k.pos + k.ref_len - 1), 1);
        }
    }
    cr_index(var_pos_cr);

    cgranges_t* chunk_noisy = intervals_to_cr(chunk.noisy_regions);
    if (chunk_noisy == nullptr) {
        chunk_noisy = cr_init();
    } else {
        cr_index(chunk_noisy);
    }

    cgranges_t* noisy_var_cr = cr_init();
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    const hts_pos_t reg_beg = chunk.region.beg;
    const hts_pos_t reg_end = chunk.region.end;

    // Pass 2: noisy mask, dense overlaps via var_pos_cr, noisy_var_cr merge, LOW_AF→LOW_COV.
    for (size_t i = 0; i < variants.size(); ++i) {
        VariantCategory c = cats[i];
        // longcallD `classify_cand_vars` pass 2: skip `NON_VAR` and `STRAND_BIAS` before overlap checks.
        if (c == VariantCategory::NonVariant || c == VariantCategory::StrandBias) continue;

        if (chunk_noisy->n_r > 0) {
            int64_t n = 0;
            if (variants[i].key.type == VariantType::Insertion) {
                n = cr_overlap(chunk_noisy, "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos),
                               &ovlp_b, &max_b);
            } else {
                n = cr_overlap(chunk_noisy, "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                               &ovlp_b, &max_b);
            }
            if (n > 0) {
                cats[i] = VariantCategory::NonVariant;
                continue;
            }
        }

        if (c == VariantCategory::LowCoverage) continue;

        if (c == VariantCategory::RepeatHetIndel) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, false, opts);
            }
            continue;
        }

        int64_t n_ov = 0;
        if (variants[i].key.type == VariantType::Insertion) {
            n_ov = cr_overlap(var_pos_cr, "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos),
                              &ovlp_b, &max_b);
        } else {
            n_ov = cr_overlap(var_pos_cr, "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                              &ovlp_b, &max_b);
        }
        if (n_ov > 1) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, true, opts);
            }
        }

        // Match longcallD: LOW_AF is rewritten to LOW_COV in var_i_to_cate after this loop.
        if (c == VariantCategory::LowAlleleFraction) cats[i] = VariantCategory::LowCoverage;
    }

    std::free(ovlp_b);
    cr_destroy(var_pos_cr);

    if (noisy_var_cr->n_r > 0) {
        cr_index(noisy_var_cr);
        cgranges_t* merged =
            cr_merge2(chunk_noisy, noisy_var_cr, -1, opts.noisy_reg_merge_dis, opts.min_sv_len);
        cr_destroy(chunk_noisy);
        cr_destroy(noisy_var_cr);
        chunk_noisy = merged;
    } else {
        cr_destroy(noisy_var_cr);
    }

    intervals_from_cr(chunk_noisy, chunk.noisy_regions);
    cr_destroy(chunk_noisy);

    for (size_t i = 0; i < variants.size(); ++i) {
        variants[i].counts.category = cats[i];
        variants[i].lcd_var_i_to_cate = category_to_flag(cats[i]);
    }
}

/**
 * @brief Classifies all chunk candidates (longcallD classify path).
 * @param chunk Chunk with allele counts populated.
 * @param opts Classification options.
 * @param header Reserved (unused).
 */
void classify_chunk_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    (void)header;
    classify_cand_vars_pgphase(chunk, opts);
}


static inline int read_order_visit_count(const BamChunk& chunk) {
    return chunk.ordered_read_ids.empty() ? static_cast<int>(chunk.reads.size())
                                         : static_cast<int>(chunk.ordered_read_ids.size());
}

/** Shape of longcallD `var_site_t` for `update_read_vs_all_var_profile_from_digar` only. */
struct LcdProfileSite {
    hts_pos_t pos = 0;
    int var_type = 0;
    int ref_len = 0;
    int alt_len = 0;
    const uint8_t* alt_seq = nullptr;
};

static void lcd_make_profile_site_from_cand(const CandidateVariant& cv, LcdProfileSite* out) {
    out->pos = cv.key.pos;
    out->alt_seq = reinterpret_cast<const uint8_t*>(cv.key.alt.data());
    if (cv.key.type == VariantType::Snp) {
        out->var_type = BAM_CDIFF;
        out->ref_len = 1;
        out->alt_len = static_cast<int>(cv.key.alt.size());
    } else if (cv.key.type == VariantType::Insertion) {
        out->var_type = BAM_CINS;
        out->ref_len = 0;
        out->alt_len = static_cast<int>(cv.key.alt.size());
    } else {
        out->var_type = BAM_CDEL;
        out->ref_len = cv.key.ref_len;
        out->alt_len = 0;
    }
}

/** longcallD `make_var_site_from_digar` (`collect_var.c`), including clip rows. */
static void lcd_make_profile_site_from_digar(const DigarOp& op, LcdProfileSite* out) {
    out->pos = op.pos;
    out->var_type = BAM_CEQUAL;
    out->ref_len = 1;
    out->alt_len = op.len;
    out->alt_seq = nullptr;
    switch (op.type) {
        case DigarType::Snp:
            out->var_type = BAM_CDIFF;
            out->ref_len = 1;
            out->alt_len = static_cast<int>(op.alt.size());
            out->alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return;
        case DigarType::Insertion:
            out->var_type = BAM_CINS;
            out->ref_len = 0;
            out->alt_len = static_cast<int>(op.alt.size());
            out->alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return;
        case DigarType::Deletion:
            out->var_type = BAM_CDEL;
            out->ref_len = op.len;
            out->alt_len = 0;
            return;
        case DigarType::SoftClip:
            out->var_type = BAM_CSOFT_CLIP;
            return;
        case DigarType::HardClip:
            out->var_type = BAM_CHARD_CLIP;
            return;
        case DigarType::Equal:
        default:
            out->var_type = BAM_CEQUAL;
            return;
    }
}

static int lcd_profile_ovlp_var_site(const LcdProfileSite* var1, const LcdProfileSite* var2) {
    const int beg1 = static_cast<int>(var1->pos);
    const int end1 = static_cast<int>(var1->pos + var1->ref_len);
    const int beg2 = static_cast<int>(var2->pos);
    const int end2 = static_cast<int>(var2->pos + var2->ref_len);
    if (var1->ref_len == 0 && var2->ref_len == 0) return (beg1 == beg2) ? 1 : 0;
    if (var1->ref_len == 0) return (beg1 > beg2 && end1 < end2) ? 1 : 0;
    if (var2->ref_len == 0) return (beg2 > beg1 && end2 < end1) ? 1 : 0;
    if (beg1 >= end2 || beg2 >= end1) return 0;
    return 1;
}

static int lcd_profile_exact_comp_var_site(const LcdProfileSite* var1, const LcdProfileSite* var2) {
    hts_pos_t var1_pos = (var1->var_type == BAM_CDIFF) ? var1->pos : var1->pos - 1;
    hts_pos_t var2_pos = (var2->var_type == BAM_CDIFF) ? var2->pos : var2->pos - 1;
    if (var1_pos < var2_pos) return -1;
    if (var1_pos > var2_pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->alt_len < var2->alt_len) return -1;
    if (var1->alt_len > var2->alt_len) return 1;
    if (var1->var_type == BAM_CDIFF || var1->var_type == BAM_CINS) {
        return std::memcmp(var1->alt_seq, var2->alt_seq, static_cast<size_t>(var1->alt_len));
    }
    return 0;
}

static int lcd_profile_comp_ovlp_var_site(const LcdProfileSite* var1, const LcdProfileSite* var2, int* is_ovlp) {
    *is_ovlp = lcd_profile_ovlp_var_site(var1, var2);
    return lcd_profile_exact_comp_var_site(var1, var2);
}

/**
 * @brief Per-read candidate allele sweep (longcallD `collect_read_var_profile`).
 *
 * Visits reads in `ordered_read_ids` order when populated (same as longcallD), else `0..n-1`.
 * Populates `chunk.read_var_profile` and `chunk.read_var_cr` ([start_var_idx, end_var_idx+1) → read_i).
 */
static void collect_read_var_profile(const Options& opts, BamChunk& chunk) {
    const CandidateTable& variants = chunk.candidates;
    const int n = static_cast<int>(variants.size());
    chunk.read_var_profile.assign(chunk.reads.size(), ReadVariantProfile{});
    for (size_t i = 0; i < chunk.reads.size(); ++i)
        chunk.read_var_profile[i].read_id = static_cast<int>(i);

    cgranges_t* cr = cr_init();

    if (n == 0) {
        cr_index(cr);
        chunk.read_var_cr.reset(cr);
        return;
    }

    for (int ord = 0; ord < read_order_visit_count(chunk); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || static_cast<size_t>(read_i) >= chunk.reads.size()) continue;
        const ReadRecord& read = chunk.reads[static_cast<size_t>(read_i)];
        if (read.is_skipped || read.tid < 0) continue;

        ReadVariantProfile cur;
        cur.read_id = static_cast<int>(read_i);
        const hts_pos_t pos_start = read.beg;
        const hts_pos_t pos_end = read.end;
        int site_i = get_var_site_start(variants, pos_start, n);
        size_t digar_i = 0;
        const std::vector<DigarOp>& dig = read.digars;

        cgranges_t* noisy_cr = intervals_to_cr(read.noisy_regions);
        if (noisy_cr != nullptr) cr_index(noisy_cr);
        // longcallD `bam_utils.h` `is_in_noisy_reg(pos, noisy_regs)`: cr_overlap(..., pos, pos+1).
        const auto is_in_noisy_reg = [&](hts_pos_t pos) -> bool {
            if (noisy_cr == nullptr) return false;
            int64_t* ovlp_b = nullptr;
            int64_t max_b = 0;
            const int64_t ovlp_n =
                cr_overlap(noisy_cr, "cr", static_cast<int32_t>(pos), static_cast<int32_t>(pos + 1), &ovlp_b, &max_b);
            free(ovlp_b);
            return ovlp_n > 0;
        };

        const auto append_observation = [&](int idx, int allele, int alt_qi) {
            if (cur.start_var_idx < 0) cur.start_var_idx = idx;
            while (cur.start_var_idx + static_cast<int>(cur.alleles.size()) < idx) {
                cur.alleles.push_back(-1);
                cur.alt_qi.push_back(-1);
            }
            cur.end_var_idx = idx;
            cur.alleles.push_back(allele);
            cur.alt_qi.push_back(alt_qi);
        };

        // longcallD `update_read_vs_all_var_profile_from_digar` (`bam_utils.c`): skip NON_VAR by bitmask.
        while (site_i < n && digar_i < dig.size()) {
            if (variants[static_cast<size_t>(site_i)].lcd_var_i_to_cate == kLongcalldNonVar) {
                ++site_i;
                continue;
            }
            if (dig[digar_i].type == DigarType::Equal) {
                ++digar_i;
                continue;
            }
            LcdProfileSite var_site0{};
            lcd_make_profile_site_from_cand(variants[static_cast<size_t>(site_i)], &var_site0);
            LcdProfileSite digar_site{};
            lcd_make_profile_site_from_digar(dig[digar_i], &digar_site);

            const DigarOp& d = dig[digar_i];
            const int ave_qual = get_digar_ave_qual(d, read.qual);
            const int var_read_pos = d.qi;

            int is_ovlp = 0;
            const int ret = lcd_profile_comp_ovlp_var_site(&var_site0, &digar_site, &is_ovlp);

            if (is_ovlp == 0) {
                if (ret < 0) {
                    append_observation(site_i, 0, -1);
                    ++site_i;
                } else if (ret > 0) {
                    ++digar_i;
                } else {
                    ++site_i;
                    ++digar_i;
                }
            } else {
                if (ret == 0) {
                    int allele_i = 1;
                    if (ave_qual < opts.min_bq) allele_i = -2;
                    append_observation(site_i, allele_i, var_read_pos);
                    ++site_i;
                } else {
                    append_observation(site_i, -1, -1);
                    ++site_i;
                }
            }
        }

        // Tail loop mirrors lcd trailing `for (; var_i < n_cand_vars; ++var_i)` including `is_in_noisy_reg`.
        for (; site_i < n; ++site_i) {
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) break;
            if (is_in_noisy_reg(variants[static_cast<size_t>(site_i)].key.pos)) continue;
            append_observation(site_i, 0, -1);
        }

        if (noisy_cr != nullptr) cr_destroy(noisy_cr);

        if (cur.start_var_idx >= 0) {
            cr_add(cr, "cr", cur.start_var_idx, cur.end_var_idx + 1,
                   static_cast<int32_t>(read_i));
            chunk.read_var_profile[read_i] = std::move(cur);
        }
    }
    cr_index(cr);
    chunk.read_var_cr.reset(cr);
}

/**
 * @brief Prints longcallD-style noisy-region summary at `--verbose 2`.
 */
static void dump_all_noisy_regions(const BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    if (opts.verbose < 2) return;
    const char* chrom = "unknown";
    if (header != nullptr && chunk.region.tid >= 0 && chunk.region.tid < header->n_targets) {
        chrom = header->target_name[chunk.region.tid];
    }
    std::fprintf(stderr, "AllNoisyRegions: %s:%" PRId64 "-%" PRId64 " (%zu)\n",
                 chrom,
                 static_cast<int64_t>(chunk.region.beg),
                 static_cast<int64_t>(chunk.region.end),
                 chunk.noisy_regions.size());
    for (const Interval& region : chunk.noisy_regions) {
        std::fprintf(stderr, "NoisyRegion: %s:%" PRId64 "-%" PRId64 " (%d)\n",
                     chrom,
                     static_cast<int64_t>(region.beg - 1),
                     static_cast<int64_t>(region.end),
                     region.label);
    }
}

/**
 * @brief Per-chunk candidate collection and phasing scaffold, arranged like longcallD `collect_var_main`.
 *
 * BAM/FASTA I/O is intentionally outside this function: callers prepare `chunk.reads`,
 * reference slices, and initial noisy/low-complexity state before entering the numbered
 * steps below.
 *
 * **Phase 1 — site collection (steps 1.1–1.2)**
 * Candidate sites (SNP/INS/DEL) are extracted from read digars, then each site is visited
 * again to accumulate ref/alt/low-qual counts and optional per-read support rows.
 *
 * **Phase 2 — classification (steps 2.1–2.6)**
 * Read-level noisy intervals are merged and filtered (`pre_process_noisy_regs`), candidates
 * are labelled (CLEAN_HET_SNP, LOW_COV, STRAND_BIAS, REP_HET_INDEL, …), noisy regions are
 * widened using classified context (`post_process_noisy_regs`), and — for non-ONT reads —
 * candidates fully contained in noisy spans are demoted to NON_VAR.
 *
 * **Guard** — if both `chunk.candidates` and `chunk.noisy_regions` are empty after
 * classification the function returns early: there is nothing to phase and nothing to
 * propagate to downstream callers.
 *
 * **Phase 3 — phasing scaffold (steps 3.1–3.2)**
 * `collect_read_var_profile` builds a sparse per-read allele table: for every candidate site
 * each read covers, it records 0 (ref), 1 (alt), or -2 (low-quality alt).
 * `assign_hap_based_on_germline_het_vars_kmeans` runs **once** here with **`kCandGermlineClean`**
 * (longcallD `LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR`).
 * longcallD’s comments describe additional rounds (`+ NOISY_CAND_*`, somatic); the **second**
 * germline k-means with `LONGCALLD_CAND_GERMLINE_VAR_CATE` happens **after** noisy-region MSA (step 4),
 * which `collect_noisy_vars_step4` runs when MSA-recalled noisy candidates are merged.
 * The separate longcallD `out_somatic` path is still outside the current pgPhase model.
 *
 * **Example — two het sites phased by shared reads, one hom site**
 * @code
 *   Candidate sites:   A (pos 100, CleanHetSnp)   B (pos 200, CleanHetSnp)   C (pos 300, CleanHom)
 *
 *   Read 1: ----[A=alt]--------[B=ref]----     Agrees with (A=1, B=0) consensus
 *   Read 2: ----[A=alt]--------[B=ref]----     Reinforces (A=1, B=0)
 *   Read 3:            --------[B=alt]--[C=alt] Conflicting read if assigned to same hap; implies opposite
 *
 *   K-Means Output Phase set anchored at A:
 *     A → phase_set=100, hap_alt=1, hap_ref=2   (anchor: alt on hap1)
 *     B → phase_set=100, hap_alt=2, hap_ref=1   (polarized relative to A: alt on hap2)
 *     C → phase_set=300, hap_alt=3, hap_ref=0   (CleanHom: alt on both, no K-Means iteration involvement)
 * @endcode
 *
 * @param chunk Chunk with reads and region pre-filled by the caller.
 * @param opts Thresholds (min BQ, AF, strand p-value, SV length, verbosity, …).
 * @param header BAM header for contig names used in verbose logging.
 * @param read_support_out If non-null, per-read allele observations are appended here.
 */
/** longcallD `collect_var.c` `make_variants` lines 1481–1488 (`active_reg_beg` / `active_reg_end`). */
static bool lcd_make_variants_active_region_contains(hts_pos_t active_reg_beg,
                                                     hts_pos_t active_reg_end,
                                                     const VariantKey& key) {
    const hts_pos_t vcf_pos = (key.type == VariantType::Deletion || key.type == VariantType::Insertion)
                                  ? key.pos - 1
                                  : key.pos;
    return !(vcf_pos < active_reg_beg || vcf_pos > active_reg_end);
}

void collect_var_main(BamChunk& chunk,
                      const Options& opts,
                      const bam_hdr_t* header,
                      std::vector<ReadSupportRow>* read_support_out) {
    // 1.1. collect X/I/D candidate sites from parsed read digars.
    collect_candidate_sites_from_records(
        chunk.region, chunk.reads, chunk.candidates, opts.min_sv_len);

    // 1.2. collect reference and alternate allele support for every candidate site.
    collect_allele_counts_from_records(chunk.reads,
                                       chunk.candidates,
                                       &chunk.region,
                                       read_support_out,
                                       opts.min_bq,
                                       opts.min_sv_len);

    // 2.1. pre-process read-level noisy regions before classification.
    pre_process_noisy_regs_pgphase(chunk, opts);

    // 2.2. identify clean, repeat/noisy, low-AF, strand-biased, and low-coverage candidates.
    // 2.3. add repeat/dense candidate spans back into the noisy-region model.
    // 2.4. reserve somatic/mosaic classification for the longcallD-compatible somatic branch.
    if (opts.verbose >= 2) {
        std::fprintf(stderr,
                     "PRE_CLASSIFY noisy_regions (%zu intervals):\n",
                     chunk.noisy_regions.size());
        for (const Interval& iv : chunk.noisy_regions) {
            std::fprintf(stderr,
                         "  %" PRId64 "-%" PRId64 " (label=%d)\n",
                         static_cast<int64_t>(iv.beg),
                         static_cast<int64_t>(iv.end),
                         iv.label);
        }
    }
    classify_chunk_candidates(chunk, opts, header);

    // 2.5. post-process noisy regions using classified candidate context.
    post_process_noisy_regs_pgphase(chunk, chunk.candidates);

    // 2.6. final containment pass: candidates fully inside noisy spans become NON_VAR.
    apply_noisy_containment_filter(chunk);
    prune_not_candidate_variants(chunk);

    dump_all_noisy_regions(chunk, opts, header);
    // longcallD returns here only when no candidates and no noisy intervals; otherwise it runs
    // noisy-region MSA (step 4) — `collect_noisy_vars_step4` — even with zero candidates.
    if (chunk.candidates.empty() && chunk.noisy_regions.empty()) return;

    if (!chunk.candidates.empty()) {
        // 3.1. collect clean-region read-wise variant profiles.
        collect_read_var_profile(opts, chunk);

        // 3.2. co-phasing using clean-region SNPs + indels + clean hom (longcallD collect_var_main).
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean);
    }

    // 4. iteratively call variants in noisy regions via MSA and re-run k-means (longcallD step 4).
    collect_noisy_vars_step4(chunk, opts);

    const hts_pos_t active_reg_beg = chunk.region.beg;
    const hts_pos_t active_reg_end = chunk.region.end;
    for (CandidateVariant& cand : chunk.candidates) {
        cand.lcd_make_variants_region_pass =
            lcd_make_variants_active_region_contains(active_reg_beg, active_reg_end, cand.key);
    }
}

} // namespace pgphase_collect
