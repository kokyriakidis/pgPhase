#ifndef PGPHASE_PHASING_TYPES_HPP
#define PGPHASE_PHASING_TYPES_HPP

// Pipeline-neutral types shared by both the BAM and graph pipelines:
// constants, enums, Options, region/chunk structs, candidate/variant types,
// read profiles, and PhasingChunk.

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/sam.h>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {


// ════════════════════════════════════════════════════════════════════════════
// Constants
// ════════════════════════════════════════════════════════════════════════════

// Default thresholds and window sizes for variant calling and phasing.

constexpr int kDefaultMinMapq = 30;
constexpr int kDefaultMinBaseq = 10;
constexpr int kMinSvLen = 30;
constexpr hts_pos_t kDefaultChunkSize = 500000;
constexpr double kDefaultMaxVarRatioPerRead = 0.05;
constexpr double kDefaultMaxNoisyFracPerRead = 0.5;
constexpr int kDefaultMinDepth = 5;
constexpr int kDefaultMinAltDepth = 2;
constexpr double kDefaultMinAf = 0.20;
constexpr double kDefaultMaxAf = 0.80;
constexpr int kDefaultNoisyRegSlideWinHifi = 100;
constexpr int kDefaultNoisyRegSlideWinOnt = 25;
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
// Chunk-stitch abstain margin: adjacent chunks merge only when the absolute
// flip-vote score strictly exceeds this margin.  0 reproduces the original
// behavior (merge on any non-zero vote).  Higher values abstain on weakly
// supported boundaries to avoid over-merging discordant blocks.
constexpr int kDefaultStitchMinMargin = 0;

// Graph-only het-indel anchor gating (hybrid pipeline only).  Graph het indels
// added to k-means as CleanHetIndel can mis-orient reads when the genotype is
// unreliable.  Keep an indel anchor only when its allele fraction sits within
// kDefaultGraphIndelAfMargin of 0.5 and it has at least
// kDefaultGraphIndelMinAlt supporting alt reads.  Defaults reproduce prior
// behavior (AF window 0.2-0.8, no extra alt floor beyond min_alt_depth).
constexpr double kDefaultGraphIndelAfMargin = 0.11;
constexpr int kDefaultGraphIndelMinAlt = 0;

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
    LowCoverage,        // below min_depth or min_alt_depth
    LowAlleleFraction,  // below min_af; folded to LowCoverage in final category
    StrandBias,         // ONT strand-bias filter (Fisher exact test)
    CleanHetSnp,        // het SNP passing all filters
    CleanHetIndel,      // het indel passing all filters
    CleanHom,           // homozygous (AF > max_af)
    NoisyCandHet,       // het recalled by noisy-region MSA (one haplotype consensus)
    NoisyCandHom,       // hom recalled by noisy-region MSA (both haplotype consensuses)
    NoisyResolved,      // noisy-region variant resolved after MSA refinement
    RepeatHetIndel,     // het indel in homopolymer / STR context
    NonVariant          // non-variant site (used as placeholder)
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
    // Max gap between adjacent noisy intervals before they are merged.
    int noisy_reg_merge_dis = kNoisyRegMergeDis;
    int min_sv_len = kLongcalldMinSvLen;
    ReadTechnology read_technology = ReadTechnology::Hifi;
    double strand_bias_pval = kDefaultStrandBiasPvalOnt;
    int noisy_reg_max_xgaps = kDefaultNoisyRegMaxXgaps;
    // Chunk-stitch abstain margin (see kDefaultStitchMinMargin).  Adjacent
    // chunks merge only when |flip_hap_score| > stitch_min_margin.
    int stitch_min_margin = kDefaultStitchMinMargin;
    // Graph-only het-indel anchor gates (hybrid pipeline).  See constants above.
    double graph_indel_af_margin = kDefaultGraphIndelAfMargin;
    int graph_indel_min_alt = kDefaultGraphIndelMinAlt;
    // Step 4: noisy-region MSA options.
    int max_noisy_reg_len = 50000; // skip regions longer than this bp
    int max_noisy_reg_cov = 1000;  // skip regions with more overlapping reads
    int min_hap_full_reads = 1;    // min full-cover reads per hap for hap-aware MSA
    int min_hap_reads = 2;         // min reads per hap total
    int min_noisy_reg_size_to_sample_reads = 10000; // subsample reads in regions >= this bp
    int noisy_reg_flank_len = kNoisyRegFlankLen;     // flank for scoring boundary deletions
    // WFA2 alignment scoring for noisy-region MSA.
    int match      = 2;
    int mismatch   = 6;
    int gap_open1  = 6;
    int gap_ext1   = 2;
    int gap_open2  = 24;
    int gap_ext2   = 1;
    int gap_aln    = 1;   // gap alignment direction: 1=left-most
    double partial_aln_ratio = 1.1; // max longer/shorter ratio for partial alignment
    int verbose = 0;
    // If false (default), skip VCF rows with non-ACGT REF/ALT bases.
    bool output_ambiguous_bases = false;
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
    std::string output_phased_bam;
    // If non-empty, output phased SAM/BAM/CRAM with HP/PS tags.
    std::string output_aln;
    // Output alignment format selected by -S/-b/-C.
    OutputAlignmentFormat output_aln_format = OutputAlignmentFormat::None;
    /** If true, refine phased read alignments from per-read digars before writing output alignment. */
    bool refine_aln = false;

    /** Optional pgbam sidecar file used for fallback chunk stitching when overlap reads have no signal. */
    std::string pgbam_file;
    int pgbam_primary_polarity_margin = 2;
    int pgbam_primary_min_winning_threads = 2;
    bool pgbam_cleanup_pass = true;
    int pgbam_cleanup_polarity_margin = 2;
    int pgbam_cleanup_min_winning_threads = 1;
    bool pgbam_relaxed_cleanup_pass = true;
    int pgbam_relaxed_cleanup_polarity_margin = 1;
    int pgbam_relaxed_cleanup_min_winning_threads = 1;
    /** Reference sample name for GBZ interval queries (e.g. "CHM13"; auto-derived from FASTA if empty). */
    std::string graph_sample;
    /** Optional GBZ-base database for future graph-native snarl/read queries. */
    std::string gbz_db;
    /** Optional raw GAF alignments for future graph-native read traversal support. */
    std::string gaf_file;
    /** Optional GAF-base database/cache for graph-native ReadSet queries. */
    std::string gaf_db;
    /** Optional precomputed vg deconstruct VCF used as a development/debug graph-site catalog. */
    std::string graph_sites_vcf;
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
    hts_pos_t pos = 0; // 1-based. Insertions are between pos-1 and pos.
    VariantType type = VariantType::Snp;
    int ref_len = 0;
    std::string alt;

    // Position key for sorting; indels use pos-1 so they sort before the anchor base.
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
    // Ordered per-read alignment events: =/X/I/D/clip/refskip operations.
    std::vector<DigarOp> digars;
    // Dense-error intervals discovered while building digars; merged into chunk-level noisy regions.
    std::vector<Interval> noisy_regions;
    bool is_skipped = false;
    bool is_ont_palindrome = false;
    int n_clean_agree_snps = 0;    // populated during phasing (Step 2)
    int n_clean_conflict_snps = 0; // populated during phasing (Step 2)
    int total_cand_events = 0; // total candidate variant events (includes long-clip noisy windows)
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
    int ref_cov = 0;   // reference allele depth (alle_covs[0])
    int alt_cov = 0;   // alternate allele depth (alle_covs[1])
    int low_qual_cov = 0;
    int forward_ref = 0;  // strand_to_alle_covs[0][0]
    int reverse_ref = 0;  // strand_to_alle_covs[1][0]
    int forward_alt = 0;  // strand_to_alle_covs[0][1]
    int reverse_alt = 0;  // strand_to_alle_covs[1][1]
    // Number of distinct alleles (width of alle_covs / hap_to_alle_profile).
    int n_uniq_alles = 2;
    // Per-allele coverage counts.  When empty, only ref_cov/alt_cov are used.
    // When non-empty, k-means init picks the highest-coverage allele (lower index wins ties).
    std::vector<int> alle_covs;
    // Final category after all classification passes (noisy overlap, AF→LowCov rewrite, etc.).
    VariantCategory category = VariantCategory::LowCoverage;
    // Category from the first classification pass only (before noisy-region adjustments).
    VariantCategory candvarcate_initial = VariantCategory::LowCoverage;
    double allele_fraction = 0.0;
};

/**
 * @brief Key plus counts and reference-encoded bases for classification.
 */
struct CandidateVariant {
    VariantKey key;
    VariantCounts counts;
    uint8_t ref_base = 4;     // 0-3=ACGT, 4=unknown; SNPs only
    uint8_t alt_ref_base = 4; // INS/DEL anchor base: 0-3=ACGT consensus, 4=use ref, >3=gap (skip VCF)
    hts_pos_t phase_set = 0;
    int hap_alt = 0;
    int hap_ref = 0;
    // True for indels in homopolymer context (set by MSA gap analysis, not by classification).
    bool is_homopolymer_indel = false;
    // True when the variant's VCF POS falls inside the chunk's active region.
    // Used during tiling-overlap dedup: prefer the copy that passes this gate.
    bool lcd_make_variants_region_pass = true;
    // Category bitmask for VCF output filtering.  Composed from kCand* flag
    // constants; VCF INFO CLEAN is set when this matches kCandGermlineClean.
    uint32_t lcd_var_i_to_cate = 0;

    // Per-haplotype allele count profiles for k-means phasing.
    // Indexed [hap 0–2][allele i]; hap 0 unused, haps 1–2 are diploid.
    // Length matches counts.n_uniq_alles.  Reset each phasing call.
    std::array<std::vector<int>, 3> hap_to_alle_profile{};
    // Consensus allele per haplotype after k-means.
    // -1 = unknown; 0 = ref; 1 = alt.
    std::array<int, 3> hap_to_cons_alle{-1, -1, -1};
};

/** @brief Ordered list of candidates for one chunk or merged batch. */
using CandidateTable = std::vector<CandidateVariant>;

// Per-read allele observations across a contiguous range of candidates.
// alleles[i] corresponds to candidate start_var_idx + i.
// Values: 0=reference, 1=alternate, -1=non-informative, -2=low-quality alt.
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
struct PhasingChunk {
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
    // Interval-tree cache for noisy-read ratio computation, built lazily before classification.
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
     * merges can restore the initial category on the replacement row.
     */
    std::vector<std::pair<VariantKey, VariantCategory>> erased_clean_signal_initial;
};

} // namespace pgphase_collect

#endif
