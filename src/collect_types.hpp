#ifndef PGPHASE_COLLECT_TYPES_HPP
#define PGPHASE_COLLECT_TYPES_HPP

/**
 * @file collect_types.hpp
 * @brief Shared types, defaults, and RAII I/O helpers for collect-bam-variation.
 */

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/sam.h>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Constants
// ════════════════════════════════════════════════════════════════════════════

/** @name Defaults
 *  @brief Default thresholds and window sizes (many mirror longcallD constants).
 */
///@{

constexpr int kDefaultMinMapq = 30;
constexpr int kDefaultMinBaseq = 10;
constexpr int kMinSvLen = 30;
constexpr hts_pos_t kDefaultChunkSize = 500000;
constexpr double kDefaultMaxVarRatioPerRead = 0.05;
constexpr double kDefaultMaxNoisyFracPerRead = 0.5; // longcallD LONGCALLD_MAX_NOISY_FRAC_PER_READ
constexpr int kDefaultMinDepth = 5;
constexpr int kDefaultMinAltDepth = 2;
constexpr double kDefaultMinAf = 0.20;
constexpr double kDefaultMaxAf = 0.80;
constexpr int kDefaultNoisyRegSlideWinHifi = 100; // longcallD LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN
constexpr int kDefaultNoisyRegSlideWinOnt = 25;   // longcallD LONGCALLD_NOISY_REG_ONT_SLIDE_WIN
constexpr int kDefaultNoisyRegSlideWinShortReads = 25;
constexpr int kLongClipLength = 30;
constexpr int kClipFlank = 100;
constexpr hts_pos_t kReferenceFlank = 50000;
constexpr int kSdustThreshold = 5;
constexpr int kSdustWindow = 20;
constexpr int kNoisyRegMergeDis = 500;
constexpr int kNoisyRegFlankLen = 10;
constexpr int kLongcalldMinSvLen = 30;
constexpr double kDefaultStrandBiasPvalOnt = 0.01;
constexpr int kDefaultNoisyRegMaxXgaps = 5;
/** If `>0`, a merged noisy region is kept only if at least this many (non-skipped) reads overlap it.
 * 0 = no extra gate (longcallD `pre_process_noisy_regs` has no separate total-depth check). Use 5 to match
 * the published “total read coverage” noisy-region description. */
constexpr int kDefaultMinNoisyRegTotalDepth = 0;

///@}

// ════════════════════════════════════════════════════════════════════════════
// Enumerations
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Variant class aligned with BAM CIGAR op types used in digars.
 */
enum class VariantType : uint8_t {
    Snp = 8,       // BAM_CDIFF
    Insertion = 1, // BAM_CINS
    Deletion = 2   // BAM_CDEL
};

enum class VariantCategory : uint8_t {
    /** `LONGCALLD_LOW_COV_VAR` — first-pass letter `L` (low alt depth). */
    LowCoverage,
    /**
     * `LONGCALLD_LOW_AF_VAR` — first-pass letter `l` (low allele fraction).
     * After the second classify loop, longcallD stores this as LOW_COV in `var_i_to_cate`
     * (collect_var.c ~981–983); we fold the same way before using `cats[]` as final category.
     */
    LowAlleleFraction,
    StrandBias,
    CleanHetSnp,
    CleanHetIndel,
    CleanHom,
    /**
     * longcallD `LONGCALLD_NOISY_CAND_HET_VAR` (`e`). Assigned when noisy-region MSA recalls a
     * variant seen on one consensus haplotype only (`update_cand_var_profile_from_cons_aln_str2`).
     * Digar `classify_cand_vars` does **not** set this — same as longcallD (candidates stay
     * `CLEAN_*` / `REP_HET_INDEL` until containment or later stages).
     */
    NoisyCandHet,
    /** longcallD `LONGCALLD_NOISY_CAND_HOM_VAR` (`h`); same provenance as `NoisyCandHet`. */
    NoisyCandHom,
    NoisyResolved,
    /** longcallD `LONGCALLD_REP_HET_VAR`: homopolymer / STR indels in clean AF band. */
    RepeatHetIndel,
    NonVariant
};

/**
 * @brief Per-read alignment event type in the digar stream (match, SNP, indel, clip, skip).
 */
enum class DigarType : uint8_t {
    Equal,
    Snp,
    Insertion,
    Deletion,
    SoftClip,
    HardClip,
    RefSkip
};

/**
 * @brief Sequencing preset: affects noisy window size and ONT strand-bias testing.
 */
enum class ReadTechnology : uint8_t {
    Hifi,
    Ont,
    ShortReads
};

enum class OutputAlignmentFormat : uint8_t {
    None,
    Sam,
    Bam,
    Cram
};

// ════════════════════════════════════════════════════════════════════════════
// Options & region input
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Parsed CLI state for collect-bam-variation (threads, thresholds, paths, regions).
 */
struct Options {
    int threads = 1;
    int min_mapq = kDefaultMinMapq;
    int min_bq = kDefaultMinBaseq;
    int min_depth = kDefaultMinDepth;
    int min_alt_depth = kDefaultMinAltDepth;
    hts_pos_t chunk_size = kDefaultChunkSize;
    double min_af = kDefaultMinAf;
    double max_af = kDefaultMaxAf;
    double max_var_ratio_per_read = kDefaultMaxVarRatioPerRead;
    double max_noisy_frac_per_read = kDefaultMaxNoisyFracPerRead;
    int noisy_reg_slide_win = -1; // if <0, pick default based on read_technology
    /** longcallD: `opt->noisy_reg_merge_dis` / `opt->min_sv_len` in `pre_process_noisy_regs` + `cr_merge`. */
    int noisy_reg_merge_dis = kNoisyRegMergeDis;
    int min_sv_len = kLongcalldMinSvLen;
    ReadTechnology read_technology = ReadTechnology::Hifi;
    double strand_bias_pval = kDefaultStrandBiasPvalOnt;
    int noisy_reg_max_xgaps = kDefaultNoisyRegMaxXgaps;
    int min_noisy_reg_total_depth = kDefaultMinNoisyRegTotalDepth;
    // Step 4: noisy-region MSA options (mirrors longcallD call_var_opt_t fields).
    int max_noisy_reg_len = 50000; ///< Skip regions longer than this bp (LONGCALLD_MAX_NOISY_REG_LEN).
    int max_noisy_reg_cov = 1000;  ///< Skip regions with more overlapping reads than this (LONGCALLD_MAX_NOISY_REG_COV).
    int min_hap_full_reads = 1;    ///< Min full-cover reads per hap for hap-aware MSA (LONGCALLD_MIN_HAP_FULL_READS).
    int min_hap_reads = 2;         ///< Min reads per hap total (LONGCALLD_MIN_HAP_READS).
    int min_noisy_reg_size_to_sample_reads = 10000; ///< Sample reads for re-alignment in regions >= this (LONGCALLD_MIN_NOISY_REG_SIZE_TO_SAMPLE_READS).
    int noisy_reg_flank_len = kNoisyRegFlankLen; ///< Flank used when scoring boundary deletions (LONGCALLD_NOISY_REG_FLANK_LEN).
    // Alignment scoring (mirrors longcallD call_var_opt_t / align.h defaults).
    int match      = 2;   ///< Match score (LONGCALLD_MATCH_SCORE).
    int mismatch   = 6;   ///< Mismatch penalty (LONGCALLD_MISMATCH_SCORE).
    int gap_open1  = 6;   ///< Gap-open penalty 1 (LONGCALLD_GAP_OPEN1_SCORE).
    int gap_ext1   = 2;   ///< Gap-extend penalty 1 (LONGCALLD_GAP_EXT1_SCORE).
    int gap_open2  = 24;  ///< Gap-open penalty 2 (LONGCALLD_GAP_OPEN2_SCORE).
    int gap_ext2   = 1;   ///< Gap-extend penalty 2 (LONGCALLD_GAP_EXT2_SCORE).
    int gap_aln    = 1;   ///< Gap alignment direction: 1=left-most (LONGCALLD_GAP_LEFT_ALN).
    double partial_aln_ratio = 1.1; ///< Max longer/shorter ratio for partial alignment (LONGCALLD_PARTIAL_ALN_RATIO).
    int verbose = 0;
    bool include_filtered = false;
    bool autosome = false;
    bool input_is_list = false;
    std::string region_file;
    std::vector<std::string> regions;
    std::string ref_fasta;
    std::string bam_file;
    std::vector<std::string> bam_files;
    std::string output_tsv = "output.tsv";
    std::string output_vcf;
    std::string output_phased_vcf;
    /** If non-empty, output phased SAM/BAM/CRAM with HP/PS tags (longcallD-style -S/-b/-C). */
    std::string output_aln;
    /** Output alignment format selected by -S/-b/-C (longcallD style). */
    OutputAlignmentFormat output_aln_format = OutputAlignmentFormat::None;
    /** If true, refine phased read alignments from per-read digars before writing output alignment. */
    bool refine_aln = false;
    /** If non-empty, write per-read allele observations for downstream phasing. */
    std::string read_support_tsv;
    /** If non-empty, write per-read HAP / PHASE_SET after k-means phasing (per chunk). */
    std::string phase_read_tsv;
    /** Optional pgbam sidecar file used for fallback chunk stitching when overlap reads have no signal. */
    std::string pgbam_file;
    std::string debug_site; // CHR:POS, emits per-read digar hits to stderr
    /** CLI command string used for PG:CL header field in phased SAM/BAM/CRAM output. */
    std::string command_line;

    /**
     * @brief First BAM path: `bam_files.front()` when multi-input, else legacy `bam_file`.
     */
    const std::string& primary_bam_file() const {
        return bam_files.empty() ? bam_file : bam_files.front();
    }

    /** @brief True when `--ont` mode (Fisher strand bias, shorter noisy window default). */
    bool is_ont() const {
        return read_technology == ReadTechnology::Ont;
    }

    /** @brief True when `--short-reads` mode. */
    bool is_short_reads() const {
        return read_technology == ReadTechnology::ShortReads;
    }
};

// ════════════════════════════════════════════════════════════════════════════
// Region & read data
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief User-specified inclusion interval on one contig (1-based inclusive coordinates).
 */
struct RegionFilter {
    bool enabled = false;
    std::string chrom;
    hts_pos_t beg = 1;  // 1-based inclusive
    hts_pos_t end = -1; // 1-based inclusive, -1 means contig end
};

/**
 * @brief One tiled window on a contig plus neighbor metadata for chunk-boundary overlap logic.
 */
struct RegionChunk {
    int chunk_id = -1;
    int tid = -1;
    hts_pos_t beg = 1; // 1-based inclusive
    hts_pos_t end = 0; // 1-based inclusive
    int reg_chunk_i = -1;
    int reg_i = -1;
    int prev_chunk_id = -1;
    int prev_tid = -1;
    hts_pos_t prev_beg = 0;
    hts_pos_t prev_end = 0;
    int next_chunk_id = -1;
    int next_tid = -1;
    hts_pos_t next_beg = 0;
    hts_pos_t next_end = 0;

    /** @brief True if `prev_chunk_id` points to the preceding chunk on the same contig. */
    bool has_prev_region() const { return prev_chunk_id >= 0; }
    /** @brief True if `next_chunk_id` points to the following chunk on the same contig. */
    bool has_next_region() const { return next_chunk_id >= 0; }
};

/**
 * @brief Candidate identity: contig, 1-based position, type, ref length, alternate sequence.
 */
struct VariantKey {
    int tid = -1;
    hts_pos_t pos = 0; // 1-based. Insertions are between pos - 1 and pos, matching longcallD.
    VariantType type = VariantType::Snp;
    int ref_len = 0;
    std::string alt;

    /**
     * @brief Position key for sorting (indels use `pos - 1` like longcallD).
     */
    hts_pos_t sort_pos() const {
        return type == VariantType::Snp ? pos : pos - 1;
    }
};

/**
 * @brief Closed interval on the reference with optional integer label (noisy size, merge metadata).
 */
struct Interval {
    hts_pos_t beg = 0; // 1-based inclusive
    hts_pos_t end = 0; // 1-based inclusive
    int label = 0;
};

/**
 * @brief One CIGAR/MD/CS-derived event on a read (position, type, length, bases, quality flag).
 */
struct DigarOp {
    hts_pos_t pos = 0; // 1-based reference coordinate
    DigarType type = DigarType::Equal;
    int len = 0;
    int qi = 0;
    bool low_quality = false;
    std::string alt;
};


/** @brief `std::unique_ptr` deleter calling `bam_destroy1`. */
struct AlignmentDeleter {
    void operator()(bam1_t* p) const { bam_destroy1(p); }
};

/** @brief `std::unique_ptr` deleter calling `cr_destroy`. */
struct CgrangesDeleter {
    void operator()(cgranges_t* p) const { if (p) cr_destroy(p); }
};

/**
 * @brief One line of `--read-support` TSV: read evidence at a candidate locus.
 */
struct ReadSupportRow {
    int tid = -1;
    hts_pos_t pos = 0;
    VariantType type = VariantType::Snp;
    int ref_len = 0;
    std::string alt;
    std::string qname;
    int is_alt = 0;       // 1 = read carries alt allele at site, 0 = ref
    int is_low_qual = 0;  // 1 = alt/low-qual observation (matches low_qual_cov path)
    bool reverse = false;
    int mapq = 0;
    hts_pos_t chunk_beg = 0;
    hts_pos_t chunk_end = 0;
};

/**
 * @brief One input read after parsing: coordinates, digars, qualities, noisy subregions.
 */
struct ReadRecord {
    int tid = -1;
    int input_index = 0;
    hts_pos_t beg = 0;
    hts_pos_t end = 0;
    bool reverse = false;
    int nm = 0;
    int mapq = 0;
    std::string qname;
    // Kept for future phasing/BAM-output consumers; candidate collection itself works from digars.
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment;
    // Per-base qualities are copied out because allele quality/depth accounting needs them after parsing.
    std::vector<uint8_t> qual;
    // longcallD digar equivalent: ordered per-read =/X/I/D/clip/refskip operations.
    std::vector<DigarOp> digars;
    // Dense-error intervals discovered while building digars; merged into chunk-level noisy regions.
    std::vector<Interval> noisy_regions;
    bool is_skipped = false;
    bool is_ont_palindrome = false;
    int n_clean_agree_snps = 0;    // populated during phasing (Step 2)
    int n_clean_conflict_snps = 0; // populated during phasing (Step 2)
    int total_cand_events = 0; // longcallD n_total_cand_vars (includes long-clip noisy windows)
};

/**
 * @brief Counts of reads skipped at chunk boundaries when loading from a BAM slice.
 */
struct OverlapSkipCounts {
    int upstream = 0;
    int downstream = 0;
};

/**
 * @brief Parsed pgbam sidecar map: set_id -> graph thread IDs.
 */
struct PgbamSidecarData {
    std::unordered_map<uint32_t, std::vector<uint64_t>> set_to_threads;
};

// ════════════════════════════════════════════════════════════════════════════
// Variant data
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Per-candidate depth, strand tallies, categories, and allele fraction.
 */
struct VariantCounts {
    int total_cov = 0;
    int ref_cov = 0;   // alle_covs[0] in longcallD
    int alt_cov = 0;   // alle_covs[1] in longcallD
    int low_qual_cov = 0;
    int forward_ref = 0;  // strand_to_alle_covs[0][0]
    int reverse_ref = 0;  // strand_to_alle_covs[1][0]
    int forward_alt = 0;  // strand_to_alle_covs[0][1]
    int reverse_alt = 0;  // strand_to_alle_covs[1][1]
    /** longcallD `n_uniq_alles`; width of `alle_covs` / `hap_to_alle_profile` when multi-allele. */
    int n_uniq_alles = 2;
    /**
     * longcallD `alle_covs` (optional). When empty, init consensus uses `ref_cov`/`alt_cov` as
     * alleles 0/1 only. When non-empty, `get_var_init_max_cov_allele` scans with strict `>` ties
     * (lower index wins), matching longcallD.
     */
    std::vector<int> alle_covs;
    /** After full classify_cand_vars (noisy overlap, LOW_AF→LOW_COV rewrite, etc.). */
    VariantCategory category = VariantCategory::LowCoverage;
    /** longcallD first-loop `classify_var_cate` only (matches CandVarCate before `After classify var:`). */
    VariantCategory candvarcate_initial = VariantCategory::LowCoverage;
    double allele_fraction = 0.0;
};

/**
 * @brief Key plus counts and reference-encoded bases for classification.
 */
struct CandidateVariant {
    VariantKey key;
    VariantCounts counts;
    uint8_t ref_base = 4;     // 0-3=ACGT, 4=unknown; SNPs only (longcallD cand_var_t::ref_base)
    uint8_t alt_ref_base = 4; // 0-3=ACGT, 4=unknown; INS/DEL only (longcallD cand_var_t::alt_ref_base)
    hts_pos_t phase_set = 0;
    int hap_alt = 0;
    int hap_ref = 0;
    /** True when the variant lies in a homopolymer or STR run (longcallD `is_homopolymer_indel`). */
    bool is_homopolymer_indel = false;
    /**
     * Per-haplotype allele count profiles for k-means phasing (longcallD `hap_to_alle_profile`).
     * Indexed [hap 0–2][allele i]. hap=0 unused here; haps 1–2 are diploid. Length matches
     * `counts.n_uniq_alles` / `counts.alle_covs.size()` (default 2). Reset each phasing call.
     */
    std::array<std::vector<int>, 3> hap_to_alle_profile{};
    /**
     * Consensus allele per haplotype (longcallD `hap_to_cons_alle`).
     * -1 = unknown/undecided; 0 = ref; 1 = alt. Written by k-means and read by scoring functions.
     */
    std::array<int, 3> hap_to_cons_alle{-1, -1, -1};
};

/** @brief Ordered list of candidates for one chunk or merged batch. */
using CandidateTable = std::vector<CandidateVariant>;

/**
 * @brief Per-read candidate allele profile, matching longcallD's read-wise variant profile shape.
 *
 * `alleles[i]` corresponds to candidate `start_var_idx + i`.
 * Values follow longcallD convention: 0=reference, 1=alternate, -1=other/non-informative,
 * -2=low-quality alternate observation.
 */
struct ReadVariantProfile {
    int read_id = -1;
    int start_var_idx = -1;
    int end_var_idx = -1;
    std::vector<int> alleles;
    std::vector<int> alt_qi;
};

// ════════════════════════════════════════════════════════════════════════════
// Chunk data
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Working state for one region chunk: reads, reference slice, noisy intervals, candidates.
 */
struct BamChunk {
    RegionChunk region;
    hts_pos_t ref_beg = 0;
    hts_pos_t ref_end = 0;
    std::string ref_seq;
    std::vector<Interval> low_complexity_regions;
    std::vector<int> ordered_read_ids;
    std::vector<Interval> noisy_regions;
    std::vector<ReadRecord> reads;
    std::vector<std::vector<int>> up_ovlp_read_i;
    std::vector<std::vector<int>> down_ovlp_read_i;
    std::vector<int> n_up_ovlp_skip_reads;
    std::vector<int> n_down_ovlp_skip_reads;
    std::vector<int> haps;
    std::vector<hts_pos_t> phase_sets;
    std::vector<ReadVariantProfile> read_var_profile;
    /** Interval tree [start_var_idx, end_var_idx+1) → read_i, built by `collect_read_var_profile`. */
    std::unique_ptr<cgranges_t, CgrangesDeleter> read_var_cr;
    CandidateTable candidates;
    // longcallD var_noisy_read_cov_cr / var_noisy_read_err_cr / var_noisy_read_marks:
    // interval-tree cache for var_noisy_reads_ratio(), built lazily before classification.
    std::unique_ptr<cgranges_t, CgrangesDeleter> var_noisy_read_cov_cr; // read coverage spans
    std::unique_ptr<cgranges_t, CgrangesDeleter> var_noisy_read_err_cr; // merged XID intervals per read
    std::vector<int> var_noisy_read_marks; // per-read dedup mark (size = reads.size())
    int var_noisy_read_mark_id = 0;
    int chunk_min_qual = 0;
    int chunk_first_quar_qual = 0;
    int chunk_median_qual = 0;
    int chunk_third_quar_qual = 0;
    int chunk_max_qual = 0;
    /**
     * When `erase_non_variant_candidates_fully_in_noisy_spans` removes a contained `NON_VAR` row
     * that had a clean first-pass category, we stash `candvarcate_initial` here so step-4 MSA
     * merges can restore INIT_CAT (longcallD first `CandVarCate` block) on the replacement row.
     */
    std::vector<std::pair<VariantKey, VariantCategory>> erased_clean_signal_initial;
};

// ════════════════════════════════════════════════════════════════════════════
// RAII deleters & I/O helpers
// ════════════════════════════════════════════════════════════════════════════

/** @brief `unique_ptr` deleter for `bam_hdr_t`. */
struct HeaderDeleter {
    void operator()(bam_hdr_t* p) const { sam_hdr_destroy(p); }
};

/** @brief `unique_ptr` deleter for `hts_idx_t`. */
struct IndexDeleter {
    void operator()(hts_idx_t* p) const { hts_idx_destroy(p); }
};

/** @brief `unique_ptr` deleter for `hts_itr_t`. */
struct IteratorDeleter {
    void operator()(hts_itr_t* p) const { hts_itr_destroy(p); }
};

/** @brief `unique_ptr` deleter for `faidx_t`. */
struct FaiDeleter {
    void operator()(faidx_t* p) const { fai_destroy(p); }
};

/**
 * @brief RAII `samFile*` for BAM/CRAM input with optional CRAM reference and HTS threads.
 */
class SamFile {
public:
    /**
     * @param path Alignment path (BAM or CRAM).
     * @param threads Passed to `hts_set_threads` when greater than 1.
     * @param ref_fasta Required for CRAM decode (`hts_set_fai_filename`).
     */
    SamFile(const std::string& path, int threads, const std::string& ref_fasta = "") : fp_(sam_open(path.c_str(), "r")) {
        if (fp_ == nullptr) throw std::runtime_error("failed to open BAM/CRAM: " + path);
        const htsFormat* fmt = hts_get_format(fp_);
        if (fmt == nullptr) throw std::runtime_error("failed to inspect alignment format: " + path);
        if (fmt->format != bam && fmt->format != cram) {
            throw std::runtime_error("input alignment file must be BAM or CRAM: " + path);
        }
        if (fmt->format == cram && !ref_fasta.empty() && hts_set_fai_filename(fp_, ref_fasta.c_str()) != 0) {
            throw std::runtime_error("failed to set reference for CRAM decoding: " + path);
        }
        if (threads > 1 && hts_set_threads(fp_, threads) != 0) {
            throw std::runtime_error("failed to configure htslib threads");
        }
    }

    /** @brief Closes the alignment file. */
    ~SamFile() {
        if (fp_ != nullptr) sam_close(fp_);
    }

    SamFile(const SamFile&) = delete;
    SamFile& operator=(const SamFile&) = delete;

    /** @brief Underlying HTSlib handle. */
    samFile* get() const { return fp_; }

private:
    samFile* fp_ = nullptr;
};

/**
 * @brief Lazy per-contig reference sequence cache backed by `faidx_t`.
 */
class ReferenceCache {
public:
    explicit ReferenceCache(faidx_t* fai) : fai_(fai) {}

    /**
     * @brief Uppercase ACGTN at 1-based \a one_based_pos on \a tid, or `N` if out of range.
     */
    char base(int tid, hts_pos_t one_based_pos, const bam_hdr_t* header) {
        if (tid < 0 || tid >= header->n_targets || one_based_pos < 1) return 'N';
        load_contig(tid, header);
        if (one_based_pos > static_cast<hts_pos_t>(seq_.size())) return 'N';
        return normalize_base(seq_[static_cast<size_t>(one_based_pos - 1)]);
    }

    /**
     * @brief Concatenates `len` bases starting at \a one_based_pos; returns `"."` if `len <= 0`.
     */
    std::string subseq(int tid, hts_pos_t one_based_pos, int len, const bam_hdr_t* header) {
        if (len <= 0) return ".";
        std::string out;
        out.reserve(static_cast<size_t>(len));
        for (int i = 0; i < len; ++i) {
            out.push_back(base(tid, one_based_pos + i, header));
        }
        return out;
    }

private:
    /** @brief Maps IUPAC ambiguity to `N`; uppercases ACGTN. */
    static char normalize_base(char base) {
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
        switch (base) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                return base;
            default:
                return 'N';
        }
    }

    /** @brief Fetches whole contig into `seq_` when `tid` changes. */
    void load_contig(int tid, const bam_hdr_t* header) {
        if (tid == loaded_tid_) return;
        hts_pos_t len = 0;
        char* raw = faidx_fetch_seq64(
            fai_, header->target_name[tid], 0, static_cast<hts_pos_t>(header->target_len[tid]) - 1, &len);
        if (raw == nullptr || len < 0) {
            free(raw);
            throw std::runtime_error(std::string("failed to fetch reference contig: ") + header->target_name[tid]);
        }
        seq_.assign(raw, raw + len);
        free(raw);
        loaded_tid_ = tid;
    }

    faidx_t* fai_ = nullptr;
    int loaded_tid_ = -1;
    std::string seq_;
};

/**
 * @brief Opens or builds a FASTA index (`fai_load3` with `FAI_CREATE`).
 * @param ref_fasta Path to reference `.fa`.
 * @return New `faidx_t*`; caller owns.
 * @throws std::runtime_error If the index cannot be created.
 */
inline faidx_t* load_reference_index(const std::string& ref_fasta) {
    faidx_t* fai = fai_load3(ref_fasta.c_str(), nullptr, nullptr, FAI_CREATE);
    if (fai == nullptr) throw std::runtime_error("failed to load/build reference FASTA index: " + ref_fasta);
    return fai;
}

/**
 * @brief Per-thread alignment and reference context: opens every BAM in `opts.bam_files`.
 */
struct WorkerContext {
    /**
     * @param opts Must list at least one indexed BAM and valid `ref_fasta`.
     * @throws std::runtime_error On missing index, header read failure, or empty `bam_files`.
     */
    explicit WorkerContext(const Options& opts)
        : fai(load_reference_index(opts.ref_fasta)),
          ref(fai.get()) {
        if (opts.bam_files.empty()) throw std::runtime_error("no input BAM/CRAM files provided");
        bams.reserve(opts.bam_files.size());
        headers.reserve(opts.bam_files.size());
        indexes.reserve(opts.bam_files.size());
        for (const std::string& path : opts.bam_files) {
            bams.push_back(std::make_unique<SamFile>(path, 1, opts.ref_fasta));
            headers.emplace_back(sam_hdr_read(bams.back()->get()));
            if (!headers.back()) throw std::runtime_error("failed to read BAM header: " + path);
            indexes.emplace_back(sam_index_load(bams.back()->get(), path.c_str()));
            if (!indexes.back()) throw std::runtime_error("region chunking requires an indexed BAM/CRAM: " + path);
        }
    }

    /** @brief Header of the first (primary) BAM in `bam_files`. */
    bam_hdr_t* primary_header() const { return headers.front().get(); }

    std::vector<std::unique_ptr<SamFile>> bams;
    std::vector<std::unique_ptr<bam_hdr_t, HeaderDeleter>> headers;
    std::vector<std::unique_ptr<hts_idx_t, IndexDeleter>> indexes;
    std::unique_ptr<faidx_t, FaiDeleter> fai;
    ReferenceCache ref;
};

} // namespace pgphase_collect

#endif
