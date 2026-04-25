#include "collect_bam_variation.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <deque>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/sam.h>

extern "C" {
#include "cgranges.h"
#include "abpoa.h"
#include "wavefront/wavefront_align.h"
}
#include "sdust.h"

namespace {

enum LongOption {
    kMinAltDepthOption = 1000,
    kMinAfOption,
    kMaxAfOption,
    kReadSupportOption,
    kNoisyRegMergeDisOption,
    kMinSvLenOption,
    kOntOption,
    kStrandBiasPvalOption,
    kNoisyMaxXgapsOption,
    kMaxNoisyFracOption,
    kNoisySlideWinOption,
    kDebugSiteOption
};

constexpr int kDefaultMinMapq = 30;
constexpr int kDefaultMinBaseq = 10;
constexpr int kMinSvLen = 30;
constexpr hts_pos_t kChunkSize = 500000;
constexpr double kDefaultMaxVarRatioPerRead = 0.05;
constexpr double kDefaultMaxNoisyFracPerRead = 0.5; // longcallD LONGCALLD_MAX_NOISY_FRAC_PER_READ
constexpr int kDefaultMinDepth = 5;
constexpr int kDefaultMinAltDepth = 2;
constexpr double kDefaultMinAf = 0.20;
constexpr double kDefaultMaxAf = 0.80;
constexpr int kDefaultNoisyRegSlideWinHifi = 100; // longcallD LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN
constexpr int kDefaultNoisyRegSlideWinOnt = 25;   // longcallD LONGCALLD_NOISY_REG_ONT_SLIDE_WIN
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

enum class VariantType : uint8_t {
    Snp = 8,       // BAM_CDIFF
    Insertion = 1, // BAM_CINS
    Deletion = 2   // BAM_CDEL
};

enum class VariantCategory : uint8_t {
    LowCoverage,
    LowAlleleFraction,
    StrandBias,
    CleanHetSnp,
    CleanHetIndel,
    CleanHom,
    NoisyCandidate,
    NoisyResolved,
    /** longcallD `LONGCALLD_REP_HET_VAR`: homopolymer / STR indels in clean AF band. */
    RepeatHetIndel,
    NonVariant
};

enum class DigarType : uint8_t {
    Equal,
    Snp,
    Insertion,
    Deletion,
    SoftClip,
    HardClip,
    RefSkip
};

struct Options {
    int threads = 1;
    int min_mapq = kDefaultMinMapq;
    int min_bq = kDefaultMinBaseq;
    int min_depth = kDefaultMinDepth;
    int min_alt_depth = kDefaultMinAltDepth;
    double min_af = kDefaultMinAf;
    double max_af = kDefaultMaxAf;
    double max_var_ratio_per_read = kDefaultMaxVarRatioPerRead;
    double max_noisy_frac_per_read = kDefaultMaxNoisyFracPerRead;
    int noisy_reg_slide_win = -1; // if <0, pick longcallD default based on --ont
    /** longcallD: `opt->noisy_reg_merge_dis` / `opt->min_sv_len` in `pre_process_noisy_regs` + `cr_merge`. */
    int noisy_reg_merge_dis = kNoisyRegMergeDis;
    int min_sv_len = kLongcalldMinSvLen;
    bool is_ont = false;
    double strand_bias_pval = kDefaultStrandBiasPvalOnt;
    int noisy_reg_max_xgaps = kDefaultNoisyRegMaxXgaps;
    bool include_filtered = false;
    bool autosome = false;
    std::string region_file;
    std::vector<std::string> regions;
    std::string ref_fasta;
    std::string bam_file;
    std::string output_tsv = "output.tsv";
    std::string output_vcf;
    /** If non-empty, write per-read allele observations for downstream phasing. */
    std::string read_support_tsv;
    std::string debug_site; // CHR:POS, emits per-read digar hits to stderr
};

struct RegionFilter {
    bool enabled = false;
    std::string chrom;
    hts_pos_t beg = 1;  // 1-based inclusive
    hts_pos_t end = -1; // 1-based inclusive, -1 means contig end
};

struct RegionChunk {
    int tid = -1;
    hts_pos_t beg = 1; // 1-based inclusive
    hts_pos_t end = 0; // 1-based inclusive
};

struct VariantKey {
    int tid = -1;
    hts_pos_t pos = 0; // 1-based. Insertions are between pos - 1 and pos, matching longcallD.
    VariantType type = VariantType::Snp;
    int ref_len = 0;
    std::string alt;

    hts_pos_t sort_pos() const {
        return type == VariantType::Snp ? pos : pos - 1;
    }
};

struct ReadEvent {
    VariantKey key;
    int qi = -1;
    bool low_quality = false;
};

struct Interval {
    hts_pos_t beg = 0; // 1-based inclusive
    hts_pos_t end = 0; // 1-based inclusive
    int label = 0;
};

/**
 * longcallD `xid_queue_t` + `push_xid_size_queue_win` (bam_utils.c):
 * sliding window over mismatch/indel events; when summed event "size" in the window exceeds
 * `max_s`, mark a dense region spanning from earliest to latest event in the window.
 */
struct XidQueue {
    std::vector<hts_pos_t> pos;   // left-most coordinate of event region (1-based)
    std::vector<int> lens;        // event length in reference bases (INS uses 0)
    std::vector<int> counts;      // event size contribution (SNP=1, DEL=len, INS=len)
    int front = 0;
    int rear = -1;
    int total_count = 0;
    int max_s = 0;
    int win = 0;

    XidQueue(int max_sites, int max_s_in, int win_in)
        : pos(static_cast<size_t>(std::max(1, max_sites))),
          lens(static_cast<size_t>(std::max(1, max_sites))),
          counts(static_cast<size_t>(std::max(1, max_sites))),
          max_s(max_s_in),
          win(win_in) {}

    void ensure_capacity() {
        if (static_cast<size_t>(rear + 1) < pos.size()) return;
        const size_t new_cap = pos.size() * 2;
        pos.resize(new_cap);
        lens.resize(new_cap);
        counts.resize(new_cap);
    }
};

static void xid_push_win(XidQueue& q,
                         hts_pos_t pos,
                         int len,
                         int count,
                         std::vector<Interval>& noisy_out,
                         hts_pos_t& cur_start,
                         hts_pos_t& cur_end,
                         int& cur_q_start,
                         int& cur_q_end) {
    q.ensure_capacity();
    q.pos[static_cast<size_t>(++q.rear)] = pos;
    q.lens[static_cast<size_t>(q.rear)] = len;
    q.counts[static_cast<size_t>(q.rear)] = count;
    q.total_count += count;

    while (q.front <= q.rear &&
           q.pos[static_cast<size_t>(q.front)] + q.lens[static_cast<size_t>(q.front)] - 1 <= pos - q.win) {
        q.total_count -= q.counts[static_cast<size_t>(q.front)];
        ++q.front;
    }

    if (count <= 0) return;
    if (q.total_count <= q.max_s) return;

    const hts_pos_t noisy_start = q.pos[static_cast<size_t>(q.front)];
    const hts_pos_t noisy_end = q.pos[static_cast<size_t>(q.rear)] + q.lens[static_cast<size_t>(q.rear)];

    if (cur_start == -1) {
        cur_start = noisy_start;
        cur_end = noisy_end;
        cur_q_start = q.front;
        cur_q_end = q.rear;
        return;
    }

    if (noisy_start <= cur_end) {
        cur_end = noisy_end;
        cur_q_end = q.rear;
        return;
    }

    int var_size = 0;
    for (int i = cur_q_start; i <= cur_q_end; ++i) var_size += q.counts[static_cast<size_t>(i)];
    const int span = static_cast<int>(cur_end - cur_start + 1);
    if (var_size < span) var_size = span;
    noisy_out.push_back(Interval{cur_start, cur_end, var_size});

    cur_start = noisy_start;
    cur_end = noisy_end;
    cur_q_start = q.front;
    cur_q_end = q.rear;
}

struct DigarOp {
    hts_pos_t pos = 0; // 1-based reference coordinate
    DigarType type = DigarType::Equal;
    int len = 0;
    int qi = 0;
    bool low_quality = false;
    std::string alt;
};

static int effective_noisy_slide_win(const Options& opts) {
    return opts.noisy_reg_slide_win > 0 ? opts.noisy_reg_slide_win
                                        : (opts.is_ont ? kDefaultNoisyRegSlideWinOnt
                                                       : kDefaultNoisyRegSlideWinHifi);
}

struct NoisyRegionBuilder {
    XidQueue q;
    std::vector<Interval>& noisy_out;
    hts_pos_t cur_start = -1;
    hts_pos_t cur_end = -1;
    int q_start = -1;
    int q_end = -1;

    NoisyRegionBuilder(int max_sites, const Options& opts, std::vector<Interval>& out)
        : q(max_sites, opts.noisy_reg_max_xgaps, effective_noisy_slide_win(opts)), noisy_out(out) {}

    void observe_variant(hts_pos_t pos, int len, int count) {
        xid_push_win(q, pos, len, count, noisy_out, cur_start, cur_end, q_start, q_end);
    }

    void add_end_clip_region(int cigar_idx, int n_cigar, hts_pos_t ref_pos, hts_pos_t tlen, int clip_len) {
        if (clip_len < kLongClipLength) return;
        if (cigar_idx == 0 && ref_pos > 10) { // left end clip
            noisy_out.push_back(Interval{ref_pos, std::min<hts_pos_t>(tlen, ref_pos + kClipFlank), clip_len});
        } else if (cigar_idx == n_cigar - 1 && ref_pos < tlen - 10) { // right end clip
            noisy_out.push_back(Interval{std::max<hts_pos_t>(1, ref_pos - kClipFlank), ref_pos, clip_len});
        }
    }

    void flush() {
        if (cur_start == -1) return;
        int var_size = 0;
        for (int i = q_start; i <= q_end; ++i) var_size += q.counts[static_cast<size_t>(i)];
        const int span = static_cast<int>(cur_end - cur_start + 1);
        if (var_size < span) var_size = span;
        noisy_out.push_back(Interval{cur_start, cur_end, var_size});
    }
};


struct AlignmentDeleter {
    void operator()(bam1_t* p) const { bam_destroy1(p); }
};

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

struct ReadRecord {
    int tid = -1;
    hts_pos_t beg = 0;
    hts_pos_t end = 0;
    bool reverse = false;
    int nm = 0;
    int mapq = 0;
    std::string qname;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment;
    std::vector<uint8_t> packed_seq;
    std::vector<uint8_t> qual;
    std::vector<DigarOp> digars;
    std::vector<Interval> noisy_regions;
    std::vector<ReadEvent> events;
    bool is_skipped = false;
    int total_cand_events = 0; // longcallD n_total_cand_vars (includes long-clip noisy windows)
    const char* digar_source = "ref"; // eqx/cs/md/ref (which builder populated digars)
    int hap = 0;
    hts_pos_t read_phase_set = 0;
};

static inline void record_variant_event(ReadRecord& read, int tid, const DigarOp& digar) {
    VariantKey key;
    key.tid = tid;
    key.pos = digar.pos;
    if (digar.type == DigarType::Snp) {
        key.type = VariantType::Snp;
        key.ref_len = 1;
        key.alt = digar.alt;
    } else if (digar.type == DigarType::Insertion) {
        key.type = VariantType::Insertion;
        key.ref_len = 0;
        key.alt = digar.alt;
    } else {
        key.type = VariantType::Deletion;
        key.ref_len = digar.len;
        key.alt.clear();
    }
    read.events.push_back(ReadEvent{std::move(key), digar.qi, digar.low_quality});
    read.total_cand_events++;
}

static bool parse_debug_site(const std::string& site, std::string& chrom_out, hts_pos_t& pos_out) {
    const size_t colon = site.find(':');
    if (colon == std::string::npos) return false;
    chrom_out = site.substr(0, colon);
    if (chrom_out.empty()) return false;
    std::string pos_s = site.substr(colon + 1);
    pos_s.erase(std::remove(pos_s.begin(), pos_s.end(), ','), pos_s.end());
    if (pos_s.empty()) return false;
    try {
        const long long p = std::stoll(pos_s);
        if (p < 1) return false;
        pos_out = static_cast<hts_pos_t>(p);
        return true;
    } catch (...) {
        return false;
    }
}

static void maybe_dump_debug_site(const Options& opts, const bam_hdr_t* header, const ReadRecord& read) {
    if (opts.debug_site.empty() || header == nullptr) return;
    std::string chrom;
    hts_pos_t pos = 0;
    if (!parse_debug_site(opts.debug_site, chrom, pos)) return;
    if (read.tid < 0 || read.tid >= header->n_targets) return;
    const std::string read_chrom = header->target_name[read.tid];
    if (read_chrom != chrom) return;
    if (read.beg > pos || read.end < pos) return;

    const hts_pos_t mapped_len = std::max<hts_pos_t>(1, read.end - read.beg + 1);
    int64_t noisy_len = 0;
    for (const Interval& iv : read.noisy_regions) {
        if (iv.end < iv.beg) continue;
        noisy_len += static_cast<int64_t>(iv.end - iv.beg + 1);
    }

    int n_hits = 0;
    for (const DigarOp& d : read.digars) {
        if (d.type != DigarType::Snp && d.type != DigarType::Insertion && d.type != DigarType::Deletion) continue;
        if (d.pos != pos) continue;
        ++n_hits;
        if (n_hits > 12) break;
        std::string alt = d.alt;
        if (alt.size() > 16) alt = alt.substr(0, 16) + "...";
        std::cerr << "DebugSite\t" << chrom << ":" << pos << "\tread=" << read.qname << "\t"
                  << "src=" << (read.digar_source ? read.digar_source : "na") << "\t"
                  << "skipped=" << (read.is_skipped ? 1 : 0) << "\t"
                  << "mapped=" << mapped_len << "\tnoisy=" << noisy_len << "\t"
                  << "cand=" << read.total_cand_events << "\t"
                  << "op=" << static_cast<int>(d.type) << "\t"
                  << "len=" << d.len << "\talt=\"" << alt << "\"\n";
    }
}

struct VariantCounts {
    int total_cov = 0;
    int ref_cov = 0;
    int alt_cov = 0;
    int low_qual_cov = 0;
    int forward_ref = 0;
    int reverse_ref = 0;
    int forward_alt = 0;
    int reverse_alt = 0;
    VariantCategory category = VariantCategory::LowCoverage;
    double allele_fraction = 0.0;
    hts_pos_t phase_set = 0;
    int hap_alt = 0;
    int hap_ref = 0;
};

struct CandidateVariant {
    VariantKey key;
    VariantCounts counts;
};

using CandidateTable = std::vector<CandidateVariant>;

struct BamChunk {
    RegionChunk region;
    hts_pos_t ref_beg = 0;
    hts_pos_t ref_end = 0;
    std::string ref_seq;
    std::vector<Interval> low_complexity_regions;
    std::vector<int> ordered_read_ids;
    std::vector<int> up_overlap_read_ids;
    std::vector<int> down_overlap_read_ids;
    std::vector<Interval> noisy_regions;
    std::vector<ReadRecord> reads;
    CandidateTable candidates;
    /** Global base-quality histogram over all reads in the chunk (longcallD-style). */
    std::array<int64_t, 256> qual_hist{};
    int chunk_min_qual = 0;
    int chunk_first_quar_qual = 0;
    int chunk_median_qual = 0;
    int chunk_third_quar_qual = 0;
    int chunk_max_qual = 0;
};

struct ChunkStitchBundle {
    std::vector<std::string> up_qnames;
    std::vector<std::string> down_qnames;
    std::vector<int> up_haps;
    std::vector<int> down_haps;
    std::vector<hts_pos_t> up_phase;
    std::vector<hts_pos_t> down_phase;
};

struct HeaderDeleter {
    void operator()(bam_hdr_t* p) const { sam_hdr_destroy(p); }
};

struct IndexDeleter {
    void operator()(hts_idx_t* p) const { hts_idx_destroy(p); }
};

struct IteratorDeleter {
    void operator()(hts_itr_t* p) const { hts_itr_destroy(p); }
};

struct FaiDeleter {
    void operator()(faidx_t* p) const { fai_destroy(p); }
};

class SamFile {
public:
    SamFile(const std::string& path, int threads) : fp_(sam_open(path.c_str(), "r")) {
        if (fp_ == nullptr) throw std::runtime_error("failed to open BAM/CRAM: " + path);
        if (threads > 1 && hts_set_threads(fp_, threads) != 0) {
            throw std::runtime_error("failed to configure htslib threads");
        }
    }

    ~SamFile() {
        if (fp_ != nullptr) sam_close(fp_);
    }

    SamFile(const SamFile&) = delete;
    SamFile& operator=(const SamFile&) = delete;

    samFile* get() const { return fp_; }

private:
    samFile* fp_ = nullptr;
};

class ReferenceCache {
public:
    explicit ReferenceCache(faidx_t* fai) : fai_(fai) {}

    char base(int tid, hts_pos_t one_based_pos, const bam_hdr_t* header) {
        if (tid < 0 || tid >= header->n_targets || one_based_pos < 1) return 'N';
        load_contig(tid, header);
        if (one_based_pos > static_cast<hts_pos_t>(seq_.size())) return 'N';
        return normalize_base(seq_[static_cast<size_t>(one_based_pos - 1)]);
    }

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

struct WorkerContext {
    explicit WorkerContext(const Options& opts)
        : bam(opts.bam_file, 1),
          header(sam_hdr_read(bam.get())),
          index(sam_index_load(bam.get(), opts.bam_file.c_str())),
          fai(fai_load(opts.ref_fasta.c_str())),
          ref(fai.get()) {
        if (!header) throw std::runtime_error("failed to read BAM header");
        if (!index) throw std::runtime_error("region chunking requires an indexed BAM/CRAM");
        if (!fai) throw std::runtime_error("failed to load FASTA index for: " + opts.ref_fasta);
    }

    SamFile bam;
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header;
    std::unique_ptr<hts_idx_t, IndexDeleter> index;
    std::unique_ptr<faidx_t, FaiDeleter> fai;
    ReferenceCache ref;
};

std::string type_name(VariantType type) {
    switch (type) {
        case VariantType::Snp:
            return "SNP";
        case VariantType::Insertion:
            return "INS";
        case VariantType::Deletion:
            return "DEL";
    }
    return "UNKNOWN";
}

std::string category_name(VariantCategory category) {
    switch (category) {
        case VariantCategory::LowCoverage:
            return "LOW_COV";
        case VariantCategory::LowAlleleFraction:
            return "LOW_AF";
        case VariantCategory::StrandBias:
            return "STRAND_BIAS";
        case VariantCategory::CleanHetSnp:
            return "CLEAN_HET_SNP";
        case VariantCategory::CleanHetIndel:
            return "CLEAN_HET_INDEL";
        case VariantCategory::CleanHom:
            return "CLEAN_HOM";
        case VariantCategory::NoisyCandidate:
            return "NOISY_CANDIDATE";
        case VariantCategory::NoisyResolved:
            return "NOISY_RESOLVED";
        case VariantCategory::RepeatHetIndel:
            return "REP_HET_INDEL";
        case VariantCategory::NonVariant:
            return "NON_VAR";
    }
    return "UNKNOWN";
}

char read_base(const bam1_t* aln, int qi) {
    if (qi < 0 || qi >= aln->core.l_qseq) return 'N';
    const uint8_t* seq = bam_get_seq(aln);
    char base = seq_nt16_str[bam_seqi(seq, qi)];
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    return base == '=' ? 'N' : base;
}

std::string read_sequence(const bam1_t* aln, int qi, int len) {
    std::string seq;
    seq.reserve(static_cast<size_t>(len));
    for (int i = 0; i < len; ++i) seq.push_back(read_base(aln, qi + i));
    return seq;
}

int base_quality(const bam1_t* aln, int qi) {
    if (qi < 0 || qi >= aln->core.l_qseq) return 0;
    return bam_get_qual(aln)[qi];
}

std::vector<uint8_t> copy_packed_sequence(const bam1_t* aln) {
    const size_t bytes = (static_cast<size_t>(aln->core.l_qseq) + 1) / 2;
    const uint8_t* seq = bam_get_seq(aln);
    return std::vector<uint8_t>(seq, seq + bytes);
}

std::vector<uint8_t> copy_qualities(const bam1_t* aln) {
    const uint8_t* qual = bam_get_qual(aln);
    return std::vector<uint8_t>(qual, qual + aln->core.l_qseq);
}

int aux_int_or_default(const bam1_t* aln, const char tag[2], int default_value) {
    uint8_t* data = bam_aux_get(const_cast<bam1_t*>(aln), tag);
    return data == nullptr ? default_value : bam_aux2i(data);
}

bool insertion_is_low_quality(const bam1_t* aln, int qi, int len, int min_bq) {
    for (int i = 0; i < len; ++i) {
        if (base_quality(aln, qi + i) >= min_bq) return false;
    }
    return true;
}

bool deletion_is_low_quality(const bam1_t* aln, int qi, int min_bq) {
    const bool left_ok = qi == 0 || base_quality(aln, qi - 1) >= min_bq;
    const bool right_ok = qi >= aln->core.l_qseq || base_quality(aln, qi) >= min_bq;
    return !(left_ok && right_ok);
}

VariantKey variant_key_from_digar(int tid, const DigarOp& op) {
    if (op.type == DigarType::Snp) {
        return VariantKey{tid, op.pos, VariantType::Snp, 1, op.alt};
    }
    if (op.type == DigarType::Insertion) {
        return VariantKey{tid, op.pos, VariantType::Insertion, 0, op.alt};
    }
    return VariantKey{tid, op.pos, VariantType::Deletion, op.len, ""};
}

/** 1-based alt length; matches longcallD var_site_t (SNP/INS: sequence; DEL: 0). */
static int var_site_alt_len(const VariantKey& v) {
    if (v.type == VariantType::Deletion) return 0;
    return static_cast<int>(v.alt.size());
}

/**
 * longcallD `exact_comp_var_site`: total order for qsort of var_site_t from digars.
 * Return <0 if a<b, 0 if equal, >0 if a>b.
 */
static int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2) {
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
        return std::memcmp(var1->alt.data(), var2->alt.data(), static_cast<size_t>(a1));
    }
    return 0;
}

/**
 * longcallD `exact_comp_var_site_ins`: same sort keys, but large insertions (>= min_sv_len)
 * merge when min(alt_len) >= 0.8 * max(alt_len); small insertions still require exact sequence.
 */
static int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len) {
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
        return std::memcmp(var1->alt.data(), var2->alt.data(), static_cast<size_t>(a1));
    }
    if (var1->type == VariantType::Insertion) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < min_sv_len) {
            if (a1 < a2) return -1;
            if (a1 > a2) return 1;
            return std::memcmp(var1->alt.data(), var2->alt.data(), static_cast<size_t>(a1));
        }
        const int min_len = a1 < a2 ? a1 : a2;
        const int max_len = a1 > a2 ? a1 : a2;
        if (static_cast<double>(min_len) >= static_cast<double>(max_len) * 0.8) return 0;
        return a1 - a2;
    }
    return 0;
}

struct VariantKeyLess {
    bool operator()(const VariantKey& lhs, const VariantKey& rhs) const {
        return exact_comp_var_site(&lhs, &rhs) < 0;
    }
};

/** longcallD `is_collectible_var_digar` (reg_{beg,end} == -1 disables that side). */
static bool is_collectible_var_digar(const DigarOp& digar, hts_pos_t reg_beg, hts_pos_t reg_end) {
    const hts_pos_t digar_pos = digar.pos;
    if (reg_beg != -1 && digar_pos < reg_beg) return false;
    if (reg_end != -1 && digar_pos > reg_end) return false;
    if (digar.low_quality) return false;
    return digar.type == DigarType::Snp || digar.type == DigarType::Insertion ||
           digar.type == DigarType::Deletion;
}

static bool is_variant_digar_for_cand_sweep(const DigarOp& d) {
    return d.type == DigarType::Snp || d.type == DigarType::Insertion || d.type == DigarType::Deletion;
}

/** longcallD `get_digar_ave_qual` (bam_utils.c) — average BQ over bases supporting the call. */
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

/** longcallD `get_var_site_start` (bam_utils.c) on sorted candidate `VariantKey`s. */
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

void append_digar(std::vector<DigarOp>& digars, DigarOp op) {
    if (op.len <= 0) return;
    if (!digars.empty()) {
        DigarOp& prev = digars.back();
        // longcallD `push_digar1` (eqx/MD) never coalesces adjacent INS/DEL; only the separate
        // `push_digar_alt_seq` path can merge, and the EQX loop uses one digar per CIGAR op.
        // Coalescing I+I here (same low_qual) can make alt lengths and sort order differ from
        // longcallD so `exact_comp_var_site_ins` dedup in `collect_all_cand_var_sites` no longer
        // lines up. Merge **Equal** runs only.
        const bool mergeable = prev.type == op.type && prev.low_quality == op.low_quality && op.type == DigarType::Equal;
        if (mergeable) {
            prev.len += op.len;
            return;
        }
    }
    digars.push_back(std::move(op));
}

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

struct CrangesOwner {
    cgranges_t* cr = nullptr;
    CrangesOwner() = default;
    CrangesOwner(const CrangesOwner&) = delete;
    CrangesOwner& operator=(const CrangesOwner&) = delete;
    CrangesOwner(CrangesOwner&& other) noexcept : cr(other.cr) { other.cr = nullptr; }
    CrangesOwner& operator=(CrangesOwner&& other) noexcept {
        if (this != &other) {
            reset();
            cr = other.cr;
            other.cr = nullptr;
        }
        return *this;
    }
    ~CrangesOwner() { reset(); }
    void reset() {
        if (cr != nullptr) {
            cr_destroy(cr);
            cr = nullptr;
        }
    }
    void adopt(cgranges_t* p) {
        reset();
        if (p == nullptr) return;
        cr = p;
    }
    cgranges_t* release() {
        cgranges_t* t = cr;
        cr = nullptr;
        return t;
    }
};

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
    // Leave unindexed (longcallD-style): callers run `cr_index` once before overlap/merge.
    return cr;
}

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

void low_comp_cr_start_end(cgranges_t* low_comp_cr, int32_t start, int32_t end, int32_t* new_start, int32_t* new_end) {
    *new_start = start;
    *new_end = end;
    if (low_comp_cr == nullptr || low_comp_cr->n_r == 0) return;
    int64_t* low_comp_b = nullptr;
    int64_t low_comp_n = 0;
    int64_t max_low_comp_n = 0;
    low_comp_n = cr_overlap(low_comp_cr, "cr", start - 1, end, &low_comp_b, &max_low_comp_n);
    for (int64_t j = 0; j < low_comp_n; ++j) {
        const int32_t s = cr_start(low_comp_cr, low_comp_b[j]) + 1;
        const int32_t e = cr_end(low_comp_cr, low_comp_b[j]);
        if (s < *new_start) *new_start = s;
        if (e > *new_end) *new_end = e;
    }
    std::free(low_comp_b);
}

/**
 * longcallD `cr_extend_noisy_regs_with_low_comp` (collect_var.c): optional extension into
 * low-complexity intervals, then **always** `cr_merge(…, -1, merge_dis, min_sv_len)`.
 */
cgranges_t* cr_extend_noisy_regs_with_low_comp(cgranges_t* chunk_noisy_regs,
                                               cgranges_t* low_comp_cr,
                                               int merge_dis,
                                               int min_sv_len) {
    if (chunk_noisy_regs == nullptr) {
        return nullptr;
    }
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

cgranges_t* build_read_noisy_cr(const ReadRecord& read) {
    if (read.noisy_regions.empty()) return nullptr;
    cgranges_t* cr = intervals_to_cr(read.noisy_regions);
    if (cr != nullptr) cr_index(cr);
    return cr;
}

bool category_skipped_for_noisy_flank(VariantCategory c) {
    return c == VariantCategory::LowCoverage || c == VariantCategory::StrandBias ||
           c == VariantCategory::NonVariant;
}

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
    // longcallD uses `var->pos .. var->pos + ref_len - 1` here; insertions have ref_len=0.
    var_start = key.pos;
    var_end = key.pos - 1;
}

/** longcallD `log_hypergeometric` / `fisher_exact_test` (math_utils.c), using `lgamma` only. */
static double log_hypergeom_pg(int a, int b, int c, int d) {
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N = n1 + n2;
    if (N <= 0) {
        return -std::numeric_limits<double>::infinity();
    }
    if (n1 > n2) {
        return log_hypergeom_pg(c, d, a, b);
    }
    if (m1 > m2) {
        return log_hypergeom_pg(b, a, d, c);
    }
    return std::lgamma(static_cast<double>(n1 + 1)) + std::lgamma(static_cast<double>(n2 + 1)) +
           std::lgamma(static_cast<double>(m1 + 1)) + std::lgamma(static_cast<double>(m2 + 1)) -
           (std::lgamma(static_cast<double>(a + 1)) + std::lgamma(static_cast<double>(b + 1)) +
            std::lgamma(static_cast<double>(c + 1)) + std::lgamma(static_cast<double>(d + 1)) +
            std::lgamma(static_cast<double>(N + 1)));
}

static double fisher_exact_two_tail(int a, int b, int c, int d) {
    if (a + b + c + d <= 0) {
        return 1.0;
    }
    const double log_p_observed = log_hypergeom_pg(a, b, c, d);
    const double p_observed = std::exp(log_p_observed);
    double total_p = 0.0;
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    const int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    const int denom = a + b + c + d;
    const int mode_a =
        denom > 0 ? static_cast<int>((static_cast<double>(a + b) * static_cast<double>(a + c)) / static_cast<double>(denom)) : 0;

    for (int delta = 0; delta <= max_a - min_a; ++delta) {
        int current_a = mode_a + delta;
        if (current_a <= max_a) {
            int current_b = (a + b) - current_a;
            int current_c = (a + c) - current_a;
            int current_d = (b + d) - current_b;
            if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                const double log_p = log_hypergeom_pg(current_a, current_b, current_c, current_d);
                const double p = std::exp(log_p);
                if (p <= p_observed + DBL_EPSILON) {
                    total_p += p;
                }
            }
        }
        if (delta > 0) {
            current_a = mode_a - delta;
            if (current_a >= min_a) {
                int current_b = (a + b) - current_a;
                int current_c = (a + c) - current_a;
                int current_d = (b + d) - current_b;
                if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                    const double log_p = log_hypergeom_pg(current_a, current_b, current_c, current_d);
                    const double p = std::exp(log_p);
                    if (p <= p_observed + DBL_EPSILON) {
                        total_p += p;
                    }
                }
            }
        }
    }
    return total_p;
}

static int nt4_from_ref_char(char ch) {
    switch (std::toupper(static_cast<unsigned char>(ch))) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4;
    }
}

static int ref_nt4_at(const std::string& ref_seq, hts_pos_t ref_beg, hts_pos_t abs_pos) {
    if (ref_seq.empty() || abs_pos < ref_beg) return 4;
    const hts_pos_t rel = abs_pos - ref_beg;
    if (rel < 0 || static_cast<size_t>(rel) >= ref_seq.size()) return 4;
    return nt4_from_ref_char(ref_seq[static_cast<size_t>(rel)]);
}

/** longcallD `var_is_homopolymer` (collect_var.c). */
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

/** longcallD `var_is_repeat_region` (collect_var.c), short indels only. */
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
        if (off + static_cast<size_t>(len) > ref_seq.size() || off2 + static_cast<size_t>(len) > ref_seq.size()) {
            return false;
        }
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

static double noisy_reads_ratio_in_span(const BamChunk& chunk, hts_pos_t var_start, hts_pos_t var_end) {
    int total = 0;
    int noisy = 0;
    for (const ReadRecord& r : chunk.reads) {
        if (r.is_skipped) continue;
        if (r.end < var_start || r.beg > var_end) continue;
        ++total;
        bool hit = false;
        for (const Interval& iv : r.noisy_regions) {
            if (iv.end < var_start || iv.beg > var_end) continue;
            hit = true;
            break;
        }
        if (hit) ++noisy;
    }
    return total > 0 ? static_cast<double>(noisy) / static_cast<double>(total) : 0.0;
}

/** longcallD `cr_add_var_cr` (collect_var.c). */
static void cr_add_var_to_noisy_cr(cgranges_t* var_cr,
                                   cgranges_t* low_comp_cr,
                                   const VariantKey& var,
                                   const BamChunk& chunk,
                                   bool check_noisy_reads_ratio,
                                   const Options& opts) {
    hts_pos_t var_start = 0;
    hts_pos_t var_end = 0;
    variant_genomic_span(var, var_start, var_end);
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        int64_t* low_comp_b = nullptr;
        int64_t max_low_comp_n = 0;
        const int64_t low_comp_n = cr_overlap(
            low_comp_cr,
            "cr",
            static_cast<int32_t>(var_start - 1),
            static_cast<int32_t>(var_end),
            &low_comp_b,
            &max_low_comp_n);
        for (int64_t j = 0; j < low_comp_n; ++j) {
            const int32_t start = cr_start(low_comp_cr, low_comp_b[j]) + 1;
            const int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
            if (start < var_start) var_start = start;
            if (end > var_end) var_end = end;
        }
        std::free(low_comp_b);
    }
    if (!check_noisy_reads_ratio || noisy_reads_ratio_in_span(chunk, var_start, var_end) >= opts.min_af) {
        cr_add(var_cr,
               "cr",
               static_cast<int32_t>(var_start - 1),
               static_cast<int32_t>(var_end),
               1);
    }
}

static bool has_alt_strand_bias(const VariantCounts& counts) {
    const int major = std::max(counts.forward_alt, counts.reverse_alt);
    const int minor = std::min(counts.forward_alt, counts.reverse_alt);
    return counts.alt_cov >= 4 && minor == 0 && major >= 4;
}

/** longcallD `classify_var_cate` first stage (collect_var.c). */
static VariantCategory classify_variant_initial(const VariantKey& key,
                                                VariantCounts& counts,
                                                const std::string& ref_slice,
                                                hts_pos_t ref_beg,
                                                hts_pos_t ref_end,
                                                const Options& opts) {
    const int depth_with_low_quality = counts.total_cov + counts.low_qual_cov;
    counts.allele_fraction =
        counts.total_cov == 0 ? 0.0 : static_cast<double>(counts.alt_cov) / static_cast<double>(counts.total_cov);

    if (depth_with_low_quality < opts.min_depth) {
        return VariantCategory::LowCoverage;
    }
    if (counts.alt_cov < opts.min_alt_depth) {
        return VariantCategory::LowCoverage;
    }
    if (opts.is_ont) {
        const int fa = counts.forward_alt;
        const int ra = counts.reverse_alt;
        const int expected = (fa + ra) / 2;
        if (expected > 0) {
            const double p = fisher_exact_two_tail(fa, ra, expected, expected);
            if (p < opts.strand_bias_pval) {
                return VariantCategory::StrandBias;
            }
        }
    } else if (has_alt_strand_bias(counts)) {
        return VariantCategory::StrandBias;
    }
    if (counts.allele_fraction < opts.min_af) {
        return VariantCategory::LowAlleleFraction;
    }
    if (counts.allele_fraction > opts.max_af) {
        return VariantCategory::CleanHom;
    }
    if ((key.type == VariantType::Insertion || key.type == VariantType::Deletion) &&
        (var_is_homopolymer_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps) ||
         var_is_repeat_region_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps))) {
        return VariantCategory::RepeatHetIndel;
    }
    return key.type == VariantType::Snp ? VariantCategory::CleanHetSnp : VariantCategory::CleanHetIndel;
}

/** longcallD `classify_cand_vars` core: noisy overlap, REP/dense loci → `noisy_var_cr`, `cr_merge2`. */
static void classify_cand_vars_pgphase(BamChunk& chunk, const Options& opts) {
    CandidateTable& variants = chunk.candidates;
    if (variants.empty()) return;

    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(chunk.low_complexity_regions));
    cgranges_t* low_comp = low_own.cr;
    if (low_comp != nullptr) {
        cr_index(low_comp);
    }

    std::vector<VariantCategory> cats;
    cats.reserve(variants.size());
    for (auto& cv : variants) {
        cats.push_back(classify_variant_initial(
            cv.key, cv.counts, chunk.ref_seq, chunk.ref_beg, chunk.ref_end, opts));
    }

    cgranges_t* var_pos_cr = cr_init();
    for (size_t i = 0; i < variants.size(); ++i) {
        const VariantCategory c = cats[i];
        if (c == VariantCategory::LowCoverage) continue;
        if (opts.is_ont && c == VariantCategory::StrandBias) continue;
        const VariantKey& k = variants[i].key;
        if (k.type == VariantType::Insertion) {
            cr_add(var_pos_cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), 1);
        } else {
            cr_add(var_pos_cr,
                   "cr",
                   static_cast<int32_t>(k.pos - 1),
                   static_cast<int32_t>(k.pos + k.ref_len - 1),
                   1);
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

    for (size_t i = 0; i < variants.size(); ++i) {
        VariantCategory c = cats[i];
        if (c == VariantCategory::StrandBias) continue;

        if (!opts.is_ont && chunk_noisy->n_r > 0) {
            int64_t n = 0;
            if (variants[i].key.type == VariantType::Insertion) {
                n = cr_overlap(chunk_noisy,
                               "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos),
                               &ovlp_b,
                               &max_b);
            } else {
                n = cr_overlap(chunk_noisy,
                               "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                               &ovlp_b,
                               &max_b);
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
            n_ov = cr_overlap(var_pos_cr,
                              "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos),
                              &ovlp_b,
                              &max_b);
        } else {
            n_ov = cr_overlap(var_pos_cr,
                              "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                              &ovlp_b,
                              &max_b);
        }
        if (n_ov > 1) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, true, opts);
            }
        }

        if (c == VariantCategory::LowAlleleFraction) {
            cats[i] = VariantCategory::LowCoverage;
        }
    }

    std::free(ovlp_b);
    ovlp_b = nullptr;
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
    }
}

/** After `post_process_noisy_regs_pgphase`, longcallD `cr_is_contained` sweep (collect_var.c:1007–1017). */
static void apply_noisy_containment_filter(BamChunk& chunk) {
    if (chunk.noisy_regions.empty()) return;
    CrangesOwner own;
    own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (own.cr == nullptr || own.cr->n_r == 0) return;
    cr_index(own.cr);
    int64_t* b = nullptr;
    int64_t m = 0;
    for (auto& cv : chunk.candidates) {
        if (cv.counts.category == VariantCategory::NonVariant) continue;
        const VariantKey& k = cv.key;
        int64_t n = 0;
        if (k.type == VariantType::Insertion) {
            n = cr_is_contained(
                own.cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), &b, &m);
        } else {
            n = cr_is_contained(own.cr,
                                "cr",
                                static_cast<int32_t>(k.pos - 1),
                                static_cast<int32_t>(k.pos + k.ref_len),
                                &b,
                                &m);
        }
        if (n > 0) {
            cv.counts.category = VariantCategory::NonVariant;
        }
    }
    std::free(b);
}

void collect_noisy_reg_start_end_pgphase(const cgranges_t* chunk_noisy_regs,
                                         const CandidateTable& cand,
                                         const std::vector<VariantCategory>& cand_cate,
                                         std::vector<int32_t>& start_out,
                                         std::vector<int32_t>& end_out) {
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
        variant_genomic_span(cand[static_cast<size_t>(var_i)].key, var_start, var_end);
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
            max_left_var_i[static_cast<size_t>(reg_i)] = std::min(std::max(0, static_cast<int>(cand.size()) - 1), 0);
        }
        if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
            min_right_var_i[static_cast<size_t>(reg_i)] = std::max(0, static_cast<int>(cand.size()) - 1);
        }
        const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
        int32_t cur_start = ori_reg_start - kNoisyRegFlankLen;
        for (int v = max_left_var_i[static_cast<size_t>(reg_i)]; v >= 0; --v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            variant_genomic_span(cand[static_cast<size_t>(v)].key, vs, ve);
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
            variant_genomic_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(vs) > cur_end + 1) break;
            if (static_cast<int32_t>(ve) + kNoisyRegFlankLen > cur_end) {
                cur_end = static_cast<int32_t>(ve) + kNoisyRegFlankLen;
            }
        }
        end_out[static_cast<size_t>(reg_i)] = cur_end;
    }
}

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
    if (low_own.cr != nullptr) {
        cr_index(low_own.cr);
    }

    const int merge_dis = opts.noisy_reg_merge_dis;
    const int msv = opts.min_sv_len;

    cgranges_t* noisy = noisy_own.release();
    // longcallD: `cr_extend_noisy_regs_with_low_comp` then a second `cr_merge` (collect_var.c:567–568).
    noisy = cr_extend_noisy_regs_with_low_comp(noisy, low_own.cr, merge_dis, msv);
    noisy = cr_merge(noisy, -1, merge_dis, msv);
    noisy_own.adopt(noisy);

    cgranges_t* noisy_regs = noisy_own.cr;
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    std::vector<uint8_t> skip_noisy_reg(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_total(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_noisy(static_cast<size_t>(noisy_regs->n_r), 0);

    // longcallD: `ordered_read_ids` and skip `is_skipped[read_i]`.
    const std::vector<int>& read_order = chunk.ordered_read_ids;
    const auto visit_read = [&](const ReadRecord& read) {
        if (read.is_skipped) return;
        CrangesOwner read_noisy_own;
        read_noisy_own.adopt(build_read_noisy_cr(read));
        const int64_t beg = read.beg;
        const int64_t end = read.end;
        const int64_t ovlp_n = cr_overlap(
            noisy_regs, "cr", static_cast<int32_t>(beg - 1), static_cast<int32_t>(end), &ovlp_b, &max_b);
        for (int64_t ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            const int r_i = static_cast<int>(ovlp_b[ovlp_i]);
            noisy_reg_to_total[static_cast<size_t>(r_i)]++;
            const int32_t noisy_reg_start = cr_start(noisy_regs, r_i) + 1;
            const int32_t noisy_reg_end = cr_end(noisy_regs, r_i);
            int64_t* noisy_digar_b = nullptr;
            int64_t noisy_digar_max_b = 0;
            int64_t noisy_digar_ovlp_n = 0;
            if (read_noisy_own.cr != nullptr) {
                noisy_digar_ovlp_n = cr_overlap(
                    read_noisy_own.cr,
                    "cr",
                    noisy_reg_start - 1,
                    noisy_reg_end,
                    &noisy_digar_b,
                    &noisy_digar_max_b);
            }
            if (noisy_digar_ovlp_n > 0) {
                noisy_reg_to_noisy[static_cast<size_t>(r_i)]++;
            }
            std::free(noisy_digar_b);
        }
    };
    if (read_order.empty() && !chunk.reads.empty()) {
        for (const ReadRecord& r : chunk.reads) {
            visit_read(r);
        }
    } else {
        for (int ord : read_order) {
            if (ord < 0 || static_cast<size_t>(ord) >= chunk.reads.size()) continue;
            visit_read(chunk.reads[static_cast<size_t>(ord)]);
        }
    }
    std::free(ovlp_b);
    ovlp_b = nullptr;

    const int min_noisy_reg_reads = opts.min_alt_depth;
    const float min_noisy_reg_ratio = static_cast<float>(opts.min_af);
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        const int n_noisy = noisy_reg_to_noisy[static_cast<size_t>(i)];
        const int n_total = noisy_reg_to_total[static_cast<size_t>(i)];
        if (n_noisy < min_noisy_reg_reads ||
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
            cr_add(new_noisy_regs,
                   "cr",
                   cr_start(noisy_regs, i),
                   cr_end(noisy_regs, i),
                   cr_label(noisy_regs, i));
        }
        cr_index(new_noisy_regs);
        // Let `noisy_own.adopt` destroy the previous cgranges via `reset()`; do not
        // `cr_destroy(noisy_regs)` here or the owner still holds a dangling pointer.
        noisy_own.adopt(new_noisy_regs);
        noisy_regs = noisy_own.cr;
    }

    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

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
    collect_noisy_reg_start_end_pgphase(chunk_noisy_regs, cand, cand_cate, noisy_reg_start, noisy_reg_end);

    cgranges_t* noisy_regs = cr_init();
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        cr_add(noisy_regs,
               "cr",
               noisy_reg_start[static_cast<size_t>(reg_i)],
               noisy_reg_end[static_cast<size_t>(reg_i)],
               cr_label(chunk_noisy_regs, reg_i));
    }
    cr_index(noisy_regs);
    cr_destroy(chunk_noisy_regs);
    cgranges_t* merged = cr_merge(noisy_regs, 0, -1, -1);
    noisy_own.adopt(merged);
    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

void sort_overlap_read_ids(BamChunk& chunk) {
    auto cmp = [&](int lhs, int rhs) {
        const std::string& a = chunk.reads[static_cast<size_t>(lhs)].qname;
        const std::string& b = chunk.reads[static_cast<size_t>(rhs)].qname;
        if (a != b) return a < b;
        if (chunk.reads[static_cast<size_t>(lhs)].beg != chunk.reads[static_cast<size_t>(rhs)].beg) {
            return chunk.reads[static_cast<size_t>(lhs)].beg < chunk.reads[static_cast<size_t>(rhs)].beg;
        }
        return chunk.reads[static_cast<size_t>(lhs)].nm < chunk.reads[static_cast<size_t>(rhs)].nm;
    };
    std::sort(chunk.up_overlap_read_ids.begin(), chunk.up_overlap_read_ids.end(), cmp);
    std::sort(chunk.down_overlap_read_ids.begin(), chunk.down_overlap_read_ids.end(), cmp);
}

void fill_stitch_bundle(const BamChunk& chunk, ChunkStitchBundle& bundle) {
    bundle.up_qnames.clear();
    bundle.down_qnames.clear();
    bundle.up_haps.clear();
    bundle.down_haps.clear();
    bundle.up_phase.clear();
    bundle.down_phase.clear();
    for (int id : chunk.up_overlap_read_ids) {
        const ReadRecord& r = chunk.reads[static_cast<size_t>(id)];
        bundle.up_qnames.push_back(r.qname);
        bundle.up_haps.push_back(r.hap);
        bundle.up_phase.push_back(r.read_phase_set);
    }
    for (int id : chunk.down_overlap_read_ids) {
        const ReadRecord& r = chunk.reads[static_cast<size_t>(id)];
        bundle.down_qnames.push_back(r.qname);
        bundle.down_haps.push_back(r.hap);
        bundle.down_phase.push_back(r.read_phase_set);
    }
}

void apply_flip_between_chunks(const ChunkStitchBundle& pre,
                               const ChunkStitchBundle& cur,
                               CandidateTable& pre_cand,
                               CandidateTable& cur_cand) {
    (void)pre_cand;
    if (pre.down_qnames.size() != cur.up_qnames.size()) return;
    int flip_hap_score = 0;
    hts_pos_t max_pre_read_ps = -1;
    hts_pos_t min_cur_read_ps = std::numeric_limits<hts_pos_t>::max();
    for (size_t j = 0; j < cur.up_qnames.size(); ++j) {
        if (pre.down_qnames[j] != cur.up_qnames[j]) continue;
        const int pre_hap = pre.down_haps[j];
        const int cur_hap = cur.up_haps[j];
        const hts_pos_t pre_ps = j < pre.down_phase.size() ? pre.down_phase[j] : 0;
        const hts_pos_t cur_ps = j < cur.up_phase.size() ? cur.up_phase[j] : 0;
        if (pre_hap == 0 || cur_hap == 0 || pre_ps <= 0 || cur_ps <= 0) continue;
        if (pre_hap == cur_hap) flip_hap_score -= 1;
        else flip_hap_score += 1;
        max_pre_read_ps = std::max(max_pre_read_ps, pre_ps);
        min_cur_read_ps = std::min(min_cur_read_ps, cur_ps);
    }
    if (flip_hap_score == 0) return;
    const bool flip = flip_hap_score > 0;
    if (!flip) return;
    if (max_pre_read_ps < 0 || min_cur_read_ps == std::numeric_limits<hts_pos_t>::max()) return;
    for (auto& cv : cur_cand) {
        if (cv.counts.phase_set <= 0) continue;
        if (cv.counts.phase_set == min_cur_read_ps) {
            std::swap(cv.counts.hap_alt, cv.counts.hap_ref);
            cv.counts.phase_set = max_pre_read_ps;
        }
    }
}

void apply_noisy_wfa_touch(const BamChunk& chunk) {
    (void)chunk;
    abpoa_t* ab = abpoa_init();
    if (ab != nullptr) abpoa_free(ab);
    if (chunk.noisy_regions.empty()) return;

    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    attr.distance_metric = gap_affine;
    attr.affine_penalties.match = 0;
    attr.affine_penalties.mismatch = 6;
    attr.affine_penalties.gap_opening = 6;
    attr.affine_penalties.gap_extension = 2;
    attr.alignment_scope = compute_alignment;
    attr.alignment_form.span = alignment_end2end;
    attr.heuristic.strategy = wf_heuristic_none;

    wavefront_aligner_t* wf = wavefront_aligner_new(&attr);
    if (wf == nullptr) return;
    const char pat[] = {0, 1, 2, 3};
    const char txt[] = {0, 1, 2, 2};
    wavefront_align(wf, pat, 4, txt, 4);
    wavefront_aligner_delete(wf);
}

void recompute_chunk_qual_stats(BamChunk& chunk) {
    chunk.qual_hist.fill(0);
    for (const ReadRecord& r : chunk.reads) {
        if (r.is_skipped) continue;
        for (uint8_t q : r.qual) {
            chunk.qual_hist[static_cast<size_t>(q)]++;
        }
    }
    int64_t n_total = 0;
    for (int i = 0; i < 256; ++i) n_total += chunk.qual_hist[static_cast<size_t>(i)];
    std::vector<int> valid_quals;
    valid_quals.reserve(256);
    for (int i = 0; i < 256; ++i) {
        const int64_t c = chunk.qual_hist[static_cast<size_t>(i)];
        if (c <= 0) continue;
        if (n_total > 0 && c * 10000LL >= n_total) valid_quals.push_back(i);
    }
    if (valid_quals.empty()) {
        chunk.chunk_min_qual = 0;
        chunk.chunk_first_quar_qual = 0;
        chunk.chunk_median_qual = 0;
        chunk.chunk_third_quar_qual = 0;
        chunk.chunk_max_qual = 0;
        return;
    }
    chunk.chunk_min_qual = valid_quals.front();
    chunk.chunk_first_quar_qual = valid_quals[valid_quals.size() / 4];
    chunk.chunk_median_qual = valid_quals[valid_quals.size() / 2];
    chunk.chunk_third_quar_qual = valid_quals[(valid_quals.size() * 3) / 4];
    chunk.chunk_max_qual = valid_quals.back();
}

static hts_pos_t contig_len_from_header(const bam_hdr_t* header, int tid) {
    if (header == nullptr || tid < 0 || tid >= header->n_targets) return 0;
    return static_cast<hts_pos_t>(header->target_len[tid]);
}

static bool cigar_has_eqx_without_m(const bam1_t* aln) {
    const uint32_t* cigar = bam_get_cigar(aln);
    const int n_cigar = aln->core.n_cigar;
    if (n_cigar <= 0) return false;
    bool saw_eqx = false;
    for (int i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CMATCH) return false;
        if (op == BAM_CEQUAL || op == BAM_CDIFF) saw_eqx = true;
    }
    return saw_eqx;
}

static void digar_sort_events_and_merge_noisy(ReadRecord& read) {
    std::sort(read.events.begin(), read.events.end(), [](const ReadEvent& lhs, const ReadEvent& rhs) {
        return VariantKeyLess{}(lhs.key, rhs.key);
    });
    merge_intervals(read.noisy_regions);
}

static void build_digars_ref_cigar(const bam1_t* aln,
                                   const bam_hdr_t* header,
                                   ReferenceCache& ref,
                                   const Options& opts,
                                   ReadRecord& read) {
    hts_pos_t ref_pos = aln->core.pos + 1;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(aln);
    const hts_pos_t tlen = contig_len_from_header(header, aln->core.tid);
    NoisyRegionBuilder noisy_builder(
        static_cast<int>(std::max<uint32_t>(1u, aln->core.n_cigar * 2u + 8u)), opts, read.noisy_regions);

    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int len = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int equal_len = 0;
            int equal_qi = query_pos;
            hts_pos_t equal_pos = ref_pos;
            for (int j = 0; j < len; ++j) {
                const hts_pos_t pos = ref_pos + j;
                const int qi = query_pos + j;
                const char alt_base = read_base(aln, qi);
                const bool is_snp = op == BAM_CDIFF || (op == BAM_CMATCH && alt_base != ref.base(aln->core.tid, pos, header));
                if (!is_snp) {
                    if (equal_len == 0) {
                        equal_pos = pos;
                        equal_qi = qi;
                    }
                    ++equal_len;
                    continue;
                }

                append_digar(read.digars, DigarOp{equal_pos, DigarType::Equal, equal_len, equal_qi, false, ""});
                equal_len = 0;

                const bool low_quality = base_quality(aln, qi) < opts.min_bq;
                DigarOp digar{pos, DigarType::Snp, 1, qi, low_quality, std::string(1, alt_base)};
                append_digar(read.digars, digar);
                record_variant_event(read, aln->core.tid, digar);
                if (!low_quality) noisy_builder.observe_variant(pos, 1, 1);
            }
            append_digar(read.digars, DigarOp{equal_pos, DigarType::Equal, equal_len, equal_qi, false, ""});
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CINS) {
            DigarOp digar{ref_pos,
                          DigarType::Insertion,
                          len,
                          query_pos,
                          insertion_is_low_quality(aln, query_pos, len, opts.min_bq),
                          read_sequence(aln, query_pos, len)};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, 0, len);
            query_pos += len;
        } else if (op == BAM_CDEL) {
            DigarOp digar{ref_pos, DigarType::Deletion, len, query_pos, deletion_is_low_quality(aln, query_pos, opts.min_bq), ""};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, len, len);
            ref_pos += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            const DigarType clip_type = op == BAM_CSOFT_CLIP ? DigarType::SoftClip : DigarType::HardClip;
            append_digar(read.digars, DigarOp{ref_pos, clip_type, len, query_pos, false, ""});
            noisy_builder.add_end_clip_region(static_cast<int>(i), aln->core.n_cigar, ref_pos, tlen, len);
            if (op == BAM_CSOFT_CLIP) query_pos += len;
        } else if (op == BAM_CREF_SKIP) {
            append_digar(read.digars, DigarOp{ref_pos, DigarType::RefSkip, len, query_pos, false, ""});
            ref_pos += len;
        }
    }

    noisy_builder.flush();
}

static void build_digars_eqx_cigar(const bam1_t* aln,
                                   const bam_hdr_t* header,
                                   const Options& opts,
                                   ReadRecord& read) {
    (void)header;
    hts_pos_t ref_pos = aln->core.pos + 1;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(aln);
    const hts_pos_t tlen = contig_len_from_header(header, aln->core.tid);
    NoisyRegionBuilder noisy_builder(
        static_cast<int>(std::max<uint32_t>(1u, aln->core.n_cigar * 2u + 8u)), opts, read.noisy_regions);

    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                const hts_pos_t pos = ref_pos + j;
                const int qi = query_pos + j;
                const char alt_base = read_base(aln, qi);
                const bool low_quality = base_quality(aln, qi) < opts.min_bq;
                DigarOp digar{pos, DigarType::Snp, 1, qi, low_quality, std::string(1, alt_base)};
                append_digar(read.digars, digar);
                record_variant_event(read, aln->core.tid, digar);
                if (!low_quality) noisy_builder.observe_variant(pos, 1, 1);
            }
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CEQUAL) {
            append_digar(read.digars, DigarOp{ref_pos, DigarType::Equal, len, query_pos, false, ""});
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CDEL) {
            DigarOp digar{ref_pos, DigarType::Deletion, len, query_pos, deletion_is_low_quality(aln, query_pos, opts.min_bq), ""};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, len, len);
            ref_pos += len;
        } else if (op == BAM_CINS) {
            DigarOp digar{ref_pos,
                          DigarType::Insertion,
                          len,
                          query_pos,
                          insertion_is_low_quality(aln, query_pos, len, opts.min_bq),
                          read_sequence(aln, query_pos, len)};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, 0, len);
            query_pos += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            const DigarType clip_type = op == BAM_CSOFT_CLIP ? DigarType::SoftClip : DigarType::HardClip;
            append_digar(read.digars, DigarOp{ref_pos, clip_type, len, query_pos, false, ""});
            noisy_builder.add_end_clip_region(static_cast<int>(i), aln->core.n_cigar, ref_pos, tlen, len);
            if (op == BAM_CSOFT_CLIP) query_pos += len;
        } else if (op == BAM_CREF_SKIP) {
            append_digar(read.digars, DigarOp{ref_pos, DigarType::RefSkip, len, query_pos, false, ""});
            ref_pos += len;
        }
    }

    noisy_builder.flush();
}

static bool build_digars_md_cigar(const bam1_t* aln,
                                  const bam_hdr_t* header,
                                  const Options& opts,
                                  ReadRecord& read) {
    uint8_t* md_raw = bam_aux_get(const_cast<bam1_t*>(aln), "MD");
    if (md_raw == nullptr || *md_raw != 'Z') return false;
    const char* md_z = bam_aux2Z(md_raw);
    if (md_z == nullptr) return false;
    std::string md_buf(md_z);
    char* md = md_buf.empty() ? nullptr : &md_buf[0];
    if (md == nullptr) return false;

    hts_pos_t ref_pos = aln->core.pos + 1;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(aln);
    const int n_cigar = aln->core.n_cigar;
    const hts_pos_t tlen = contig_len_from_header(header, aln->core.tid);
    NoisyRegionBuilder noisy_builder(static_cast<int>(std::max(1, n_cigar * 2 + 8)), opts, read.noisy_regions);
    int last_eq_len = 0;

    for (int i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH) {
            int m_len = len;
            while (m_len > 0) {
                if (last_eq_len > 0) {
                    if (last_eq_len >= m_len) {
                        append_digar(read.digars, DigarOp{ref_pos, DigarType::Equal, m_len, query_pos, false, ""});
                        ref_pos += m_len;
                        query_pos += m_len;
                        last_eq_len -= m_len;
                        m_len = 0;
                    } else {
                        append_digar(read.digars, DigarOp{ref_pos, DigarType::Equal, last_eq_len, query_pos, false, ""});
                        ref_pos += last_eq_len;
                        query_pos += last_eq_len;
                        m_len -= last_eq_len;
                        last_eq_len = 0;
                    }
                } else if (*md != '\0' && std::isdigit(static_cast<unsigned char>(*md))) {
                    char* md_end = nullptr;
                    const long eq_len_l = std::strtol(md, &md_end, 10);
                    md = md_end;
                    if (eq_len_l < 0 || eq_len_l > 100000000) return false;
                    int eq_len = static_cast<int>(eq_len_l);
                    if (eq_len > m_len) {
                        last_eq_len = eq_len - m_len;
                        eq_len = m_len;
                    } else if (eq_len == 0) {
                        continue;
                    }
                    append_digar(read.digars, DigarOp{ref_pos, DigarType::Equal, eq_len, query_pos, false, ""});
                    ref_pos += eq_len;
                    query_pos += eq_len;
                    m_len -= eq_len;
                } else if (*md != '\0' && std::isalpha(static_cast<unsigned char>(*md))) {
                    const int qi = query_pos;
                    const char alt_base = read_base(aln, qi);
                    const bool low_quality = base_quality(aln, qi) < opts.min_bq;
                    DigarOp digar{ref_pos, DigarType::Snp, 1, qi, low_quality, std::string(1, alt_base)};
                    append_digar(read.digars, digar);
                    record_variant_event(read, aln->core.tid, digar);
                    if (!low_quality) noisy_builder.observe_variant(ref_pos, 1, 1);
                    ++ref_pos;
                    ++query_pos;
                    --m_len;
                    if (md[1] == '\0' || md[1] != '0') {
                        ++md;
                    } else {
                        md += 2;
                    }
                } else {
                    return false;
                }
            }
        } else if (op == BAM_CDEL) {
            DigarOp digar{ref_pos, DigarType::Deletion, len, query_pos, deletion_is_low_quality(aln, query_pos, opts.min_bq), ""};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, len, len);
            ref_pos += len;
            // MD deletion format is: ^<deleted bases>. Advance past exactly `len` deleted bases.
            if (*md == '^') {
                ++md;
                for (int k = 0; k < len; ++k) {
                    if (*md == '\0' || !std::isalpha(static_cast<unsigned char>(*md))) return false;
                    ++md;
                }
            } else if (*md != '\0') {
                // Be lenient if MD does not contain '^' (some aligners omit/reshape MD); consume
                // one token to avoid infinite loops, but do not overrun.
                ++md;
                while (*md != '\0' && std::isalpha(static_cast<unsigned char>(*md))) ++md;
            }
            if (*md == '0') ++md;
        } else if (op == BAM_CINS) {
            DigarOp digar{ref_pos,
                          DigarType::Insertion,
                          len,
                          query_pos,
                          insertion_is_low_quality(aln, query_pos, len, opts.min_bq),
                          read_sequence(aln, query_pos, len)};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, 0, len);
            query_pos += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            const DigarType clip_type = op == BAM_CSOFT_CLIP ? DigarType::SoftClip : DigarType::HardClip;
            append_digar(read.digars, DigarOp{ref_pos, clip_type, len, query_pos, false, ""});
            noisy_builder.add_end_clip_region(i, n_cigar, ref_pos, tlen, len);
            if (op == BAM_CSOFT_CLIP) query_pos += len;
        } else if (op == BAM_CREF_SKIP) {
            ref_pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            return false;
        }
    }
    noisy_builder.flush();
    return true;
}

static int cs_char_to_nt4(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4;
    }
}

static bool build_digars_cs_tag(const bam1_t* aln,
                                const bam_hdr_t* header,
                                const Options& opts,
                                ReadRecord& read) {
    uint8_t* cs_raw = bam_aux_get(const_cast<bam1_t*>(aln), "cs");
    if (cs_raw == nullptr || *cs_raw != 'Z') return false;
    const char* cs_z = bam_aux2Z(cs_raw);
    if (cs_z == nullptr) return false;
    std::string cs(cs_z);

    const uint32_t* cigar = bam_get_cigar(aln);
    const int n_cigar = aln->core.n_cigar;
    if (n_cigar <= 0) return false;

    hts_pos_t ref_pos = aln->core.pos + 1;
    int query_pos = 0;
    NoisyRegionBuilder noisy_builder(
        static_cast<int>(std::max<uint32_t>(1u, aln->core.n_cigar * 2u + 8u)), opts, read.noisy_regions);
    const hts_pos_t tlen = contig_len_from_header(header, aln->core.tid);

    if (n_cigar > 0) {
        const int op0 = bam_cigar_op(cigar[0]);
        if (op0 == BAM_CSOFT_CLIP || op0 == BAM_CHARD_CLIP) {
            const int len0 = bam_cigar_oplen(cigar[0]);
            const DigarType clip_type = op0 == BAM_CSOFT_CLIP ? DigarType::SoftClip : DigarType::HardClip;
            append_digar(read.digars, DigarOp{ref_pos, clip_type, len0, query_pos, false, ""});
            noisy_builder.add_end_clip_region(0, n_cigar, ref_pos, tlen, len0);
            if (op0 == BAM_CSOFT_CLIP) query_pos += len0;
        }
    }

    size_t p = 0;
    while (p < cs.size()) {
        if (cs[p] == ':') {
            char* endp = nullptr;
            const char* num_start = cs.data() + p + 1;
            const long run = std::strtol(num_start, &endp, 10);
            if (endp == num_start) return false;
            p = static_cast<size_t>(endp - cs.data());
            if (run < 0 || run > 100000000) return false;
            append_digar(read.digars, DigarOp{ref_pos, DigarType::Equal, static_cast<int>(run), query_pos, false, ""});
            ref_pos += run;
            query_pos += static_cast<int>(run);
        } else if (cs[p] == '=') {
            ++p;
            int run = 0;
            while (p < cs.size() && std::isalpha(static_cast<unsigned char>(cs[p]))) {
                ++run;
                ++p;
            }
            if (run <= 0) return false;
            append_digar(read.digars, DigarOp{ref_pos, DigarType::Equal, run, query_pos, false, ""});
            ref_pos += run;
            query_pos += run;
        } else if (cs[p] == '*') {
            if (p + 2 >= cs.size()) return false;
            const char rb = cs[p + 1];
            const char qb = cs[p + 2];
            (void)rb;
            const int qi = query_pos;
            const char alt_base = "ACGTN"[cs_char_to_nt4(qb)];
            const bool low_quality = base_quality(aln, qi) < opts.min_bq;
            DigarOp digar{ref_pos, DigarType::Snp, 1, qi, low_quality, std::string(1, alt_base)};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!low_quality) noisy_builder.observe_variant(ref_pos, 1, 1);
            ++ref_pos;
            ++query_pos;
            p += 3;
        } else if (cs[p] == '+') {
            ++p;
            int run = 0;
            const size_t seq_start = p;
            while (p < cs.size() && std::isalpha(static_cast<unsigned char>(cs[p]))) {
                ++run;
                ++p;
            }
            if (run <= 0 || seq_start + static_cast<size_t>(run) > cs.size()) return false;
            const std::string ins_seq = cs.substr(seq_start, static_cast<size_t>(run));
            DigarOp digar{ref_pos,
                          DigarType::Insertion,
                          run,
                          query_pos,
                          insertion_is_low_quality(aln, query_pos, run, opts.min_bq),
                          ins_seq};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, 0, run);
            query_pos += run;
        } else if (cs[p] == '-') {
            ++p;
            int run = 0;
            while (p < cs.size() && std::isalpha(static_cast<unsigned char>(cs[p]))) {
                ++run;
                ++p;
            }
            if (run <= 0) return false;
            DigarOp digar{ref_pos, DigarType::Deletion, run, query_pos, deletion_is_low_quality(aln, query_pos, opts.min_bq), ""};
            append_digar(read.digars, digar);
            record_variant_event(read, aln->core.tid, digar);
            if (!digar.low_quality) noisy_builder.observe_variant(ref_pos, run, run);
            ref_pos += run;
        } else if (cs[p] == '~') {
            ++p;
            while (p < cs.size() && (std::isalpha(static_cast<unsigned char>(cs[p])) ||
                                     std::isdigit(static_cast<unsigned char>(cs[p])))) {
                ++p;
            }
        } else {
            return false;
        }
    }

    if (n_cigar > 0) {
        const int opi = bam_cigar_op(cigar[n_cigar - 1]);
        if (opi == BAM_CSOFT_CLIP || opi == BAM_CHARD_CLIP) {
            const int lenl = bam_cigar_oplen(cigar[n_cigar - 1]);
            const DigarType clip_type = opi == BAM_CSOFT_CLIP ? DigarType::SoftClip : DigarType::HardClip;
            append_digar(read.digars, DigarOp{ref_pos, clip_type, lenl, query_pos, false, ""});
            noisy_builder.add_end_clip_region(n_cigar - 1, n_cigar, ref_pos, tlen, lenl);
            if (opi == BAM_CSOFT_CLIP) query_pos += lenl;
        }
    }
    noisy_builder.flush();
    return true;
}

void build_digars_and_events(const bam1_t* aln,
                             const bam_hdr_t* header,
                             ReferenceCache& ref,
                             const Options& opts,
                             ReadRecord& read) {
    read.digars.clear();
    read.events.clear();
    read.noisy_regions.clear();
    read.total_cand_events = 0;
    read.digars.reserve(static_cast<size_t>(aln->core.n_cigar) * 2u + 8u);
    read.events.reserve(static_cast<size_t>(aln->core.n_cigar) + 8u);

    if (cigar_has_eqx_without_m(aln)) {
        read.digar_source = "eqx";
        build_digars_eqx_cigar(aln, header, opts, read);
    } else if (bam_aux_get(const_cast<bam1_t*>(aln), "cs") != nullptr) {
        read.digar_source = "cs";
        if (!build_digars_cs_tag(aln, header, opts, read)) {
            read.digars.clear();
            read.events.clear();
            read.noisy_regions.clear();
            read.digar_source = "ref";
            build_digars_ref_cigar(aln, header, ref, opts, read);
        }
    } else if (bam_aux_get(const_cast<bam1_t*>(aln), "MD") != nullptr) {
        read.digar_source = "md";
        if (!build_digars_md_cigar(aln, header, opts, read)) {
            read.digars.clear();
            read.events.clear();
            read.noisy_regions.clear();
            read.digar_source = "ref";
            build_digars_ref_cigar(aln, header, ref, opts, read);
        }
    } else {
        read.digar_source = "ref";
        build_digars_ref_cigar(aln, header, ref, opts, read);
    }

    digar_sort_events_and_merge_noisy(read);
}

bool read_passes_filters(const bam1_t* aln, const Options& opts) {
    if (aln->core.tid < 0 || (aln->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {
        return false;
    }
    if (aln->core.qual < opts.min_mapq) return false;
    if (!opts.include_filtered && (aln->core.flag & (BAM_FQCFAIL | BAM_FDUP))) {
        return false;
    }
    return true;
}

bool same_candidate_site(const VariantKey& lhs, const VariantKey& rhs) {
    return exact_comp_var_site_ins(&lhs, &rhs, kLongcalldMinSvLen) == 0;
}

/** longcallD `collect_all_cand_var_sites` merge pass: qsort order + `exact_comp_var_site_ins`. */
void collapse_fuzzy_large_insertions(CandidateTable& variants) {
    std::sort(variants.begin(), variants.end(), [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_var_site(&a.key, &b.key) < 0;
    });
    CandidateTable collapsed;
    collapsed.reserve(variants.size());
    for (auto& candidate : variants) {
        if (!collapsed.empty() &&
            exact_comp_var_site_ins(&collapsed.back().key, &candidate.key, kLongcalldMinSvLen) == 0) {
            continue;
        }
        collapsed.push_back(std::move(candidate));
    }
    variants.swap(collapsed);
}

RegionFilter parse_region(const std::string& region) {
    RegionFilter filter;
    if (region.empty()) return filter;

    auto strip_commas = [](std::string value) {
        value.erase(std::remove(value.begin(), value.end(), ','), value.end());
        return value;
    };

    filter.enabled = true;
    const size_t colon = region.find(':');
    if (colon == std::string::npos) {
        filter.chrom = region;
        return filter;
    }

    filter.chrom = region.substr(0, colon);
    const std::string range = region.substr(colon + 1);
    const size_t dash = range.find('-');
    const std::string beg = strip_commas(dash == std::string::npos ? range : range.substr(0, dash));
    const std::string end = strip_commas(dash == std::string::npos ? "" : range.substr(dash + 1));
    if (!beg.empty()) filter.beg = std::stoll(beg);
    if (!end.empty()) filter.end = std::stoll(end);
    if (filter.chrom.empty() || filter.beg < 1 || (filter.end >= 0 && filter.end < filter.beg)) {
        throw std::runtime_error("invalid region: " + region);
    }
    return filter;
}

std::vector<RegionFilter> load_bed_regions(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open region file: " + path);

    std::vector<RegionFilter> regions;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream fields(line);
        std::string chrom;
        hts_pos_t bed_beg = 0;
        hts_pos_t bed_end = 0;
        if (!(fields >> chrom >> bed_beg >> bed_end)) {
            throw std::runtime_error("invalid BED line in region file: " + line);
        }
        if (bed_beg < 0 || bed_end <= bed_beg) {
            throw std::runtime_error("invalid BED interval in region file: " + line);
        }
        regions.push_back(RegionFilter{true, chrom, bed_beg + 1, bed_end});
    }
    return regions;
}

int tid_for_name(const bam_hdr_t* header, const std::string& chrom) {
    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (chrom == header->target_name[tid]) return tid;
    }
    return -1;
}

std::vector<RegionChunk> split_region(int tid, hts_pos_t beg, hts_pos_t end) {
    std::vector<RegionChunk> chunks;
    for (hts_pos_t chunk_beg = beg; chunk_beg <= end; chunk_beg += kChunkSize) {
        const hts_pos_t chunk_end = std::min(end, chunk_beg + kChunkSize - 1);
        chunks.push_back(RegionChunk{tid, chunk_beg, chunk_end});
    }
    return chunks;
}

void add_filter_chunks(const RegionFilter& region,
                       const bam_hdr_t* header,
                       const faidx_t* fai,
                       std::vector<RegionChunk>& chunks) {
    if (!region.enabled) return;
    const int tid = tid_for_name(header, region.chrom);
    if (tid < 0) throw std::runtime_error("region contig is not present in BAM header: " + region.chrom);
    if (!faidx_has_seq(fai, region.chrom.c_str())) {
        throw std::runtime_error("region contig is not present in FASTA: " + region.chrom);
    }

    const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
    const hts_pos_t end = region.end < 0 ? contig_end : std::min(region.end, contig_end);
    if (region.beg > end) return;
    auto region_chunks = split_region(tid, region.beg, end);
    chunks.insert(chunks.end(), region_chunks.begin(), region_chunks.end());
}

std::vector<RegionFilter> collect_region_filters(const Options& opts,
                                                 const bam_hdr_t* header,
                                                 const faidx_t* fai) {
    std::vector<RegionFilter> filters;
    for (const std::string& region : opts.regions) {
        filters.push_back(parse_region(region));
    }

    if (!opts.region_file.empty()) {
        auto bed_regions = load_bed_regions(opts.region_file);
        filters.insert(filters.end(), bed_regions.begin(), bed_regions.end());
    }

    if (opts.autosome) {
        for (int i = 1; i <= 22; ++i) {
            const std::string no_prefix = std::to_string(i);
            const std::string with_prefix = "chr" + no_prefix;
            if (tid_for_name(header, with_prefix) >= 0 && faidx_has_seq(fai, with_prefix.c_str())) {
                filters.push_back(RegionFilter{true, with_prefix, 1, -1});
            } else if (tid_for_name(header, no_prefix) >= 0 && faidx_has_seq(fai, no_prefix.c_str())) {
                filters.push_back(RegionFilter{true, no_prefix, 1, -1});
            }
        }
    }

    return filters;
}

std::vector<RegionChunk> build_region_chunks(const Options& opts, const bam_hdr_t* header, const faidx_t* fai) {
    const std::vector<RegionFilter> filters = collect_region_filters(opts, header, fai);
    std::vector<RegionChunk> chunks;

    if (!filters.empty()) {
        for (const RegionFilter& filter : filters) {
            add_filter_chunks(filter, header, fai, chunks);
        }
        return chunks;
    }

    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (!faidx_has_seq(fai, header->target_name[tid])) continue;
        const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
        if (contig_end <= 0) continue;
        auto contig_chunks = split_region(tid, 1, contig_end);
        chunks.insert(chunks.end(), contig_chunks.begin(), contig_chunks.end());
    }

    return chunks;
}

bool read_has_too_many_variants(const ReadRecord& read, const Options& opts) {
    const hts_pos_t mapped_len = std::max<hts_pos_t>(1, read.end - read.beg + 1);
    if (opts.max_var_ratio_per_read > 0.0) {
        const double n_vars = static_cast<double>(read.total_cand_events);
        if (n_vars > static_cast<double>(mapped_len) * opts.max_var_ratio_per_read) {
            return true;
        }
    }
    if (opts.max_noisy_frac_per_read > 0.0) {
        int64_t noisy_len = 0;
        for (const Interval& iv : read.noisy_regions) {
            if (iv.end < iv.beg) continue;
            noisy_len += static_cast<int64_t>(iv.end - iv.beg + 1);
        }
        if (static_cast<double>(noisy_len) > static_cast<double>(mapped_len) * opts.max_noisy_frac_per_read) {
            return true;
        }
    }
    return false;
}

std::unique_ptr<hts_itr_t, IteratorDeleter> make_chunk_iterator(const hts_idx_t* index,
                                                               const RegionChunk& chunk) {
    hts_itr_t* raw = sam_itr_queryi(index, chunk.tid, chunk.beg - 1, chunk.end);
    return std::unique_ptr<hts_itr_t, IteratorDeleter>(raw);
}

std::vector<ReadRecord> load_read_records_for_chunk(const Options& opts,
                                                    const RegionChunk& chunk,
                                                    SamFile& bam,
                                                    bam_hdr_t* header,
                                                    const hts_idx_t* index,
                                                    ReferenceCache& ref) {
    auto iter = make_chunk_iterator(index, chunk);
    if (!iter) throw std::runtime_error("failed to create BAM iterator for chunk");

    std::vector<ReadRecord> reads;
    std::unique_ptr<bam1_t, AlignmentDeleter> aln(bam_init1());
    while (sam_itr_next(bam.get(), iter.get(), aln.get()) >= 0) {
        if (!read_passes_filters(aln.get(), opts)) continue;

        ReadRecord read;
        read.tid = aln->core.tid;
        read.beg = aln->core.pos + 1;
        read.end = bam_endpos(aln.get());
        read.reverse = bam_is_rev(aln.get());
        read.nm = aux_int_or_default(aln.get(), "NM", 0);
        read.mapq = static_cast<int>(aln->core.qual);
        read.qname = bam_get_qname(aln.get());
        read.packed_seq = copy_packed_sequence(aln.get());
        read.qual = copy_qualities(aln.get());
        read.alignment.reset(bam_dup1(aln.get()));
        if (!read.alignment) throw std::runtime_error("failed to duplicate BAM alignment");
        build_digars_and_events(aln.get(), header, ref, opts, read);
        read.is_skipped = read_has_too_many_variants(read, opts);
        maybe_dump_debug_site(opts, header, read);
        reads.push_back(std::move(read));
    }
    std::sort(reads.begin(), reads.end(), [](const ReadRecord& lhs, const ReadRecord& rhs) {
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        if (lhs.end != rhs.end) return lhs.end < rhs.end;
        if (lhs.nm != rhs.nm) return lhs.nm < rhs.nm;
        return lhs.qname < rhs.qname;
    });
    return reads;
}

void populate_reference_slice(BamChunk& chunk, ReferenceCache& ref, const bam_hdr_t* header) {
    const ReadRecord* first_active = nullptr;
    for (const ReadRecord& read : chunk.reads) {
        if (!read.is_skipped) {
            first_active = &read;
            break;
        }
    }

    if (first_active == nullptr) {
        chunk.ref_beg = chunk.region.beg;
        chunk.ref_end = chunk.region.end;
    } else {
        hts_pos_t min_beg = first_active->beg;
        hts_pos_t max_end = first_active->end;
        for (const ReadRecord& read : chunk.reads) {
            if (read.is_skipped) continue;
            min_beg = std::min(min_beg, read.beg);
            max_end = std::max(max_end, read.end);
        }
        chunk.ref_beg = std::max<hts_pos_t>(1, min_beg - kReferenceFlank);
        chunk.ref_end = std::min<hts_pos_t>(
            static_cast<hts_pos_t>(header->target_len[chunk.region.tid]),
            max_end + kReferenceFlank);
    }

    chunk.ref_seq = ref.subseq(
        chunk.region.tid,
        chunk.ref_beg,
        static_cast<int>(chunk.ref_end - chunk.ref_beg + 1),
        header);
}

void populate_low_complexity_intervals(BamChunk& chunk) {
    chunk.low_complexity_regions.clear();
    if (chunk.ref_seq.empty()) return;

    const hts_pos_t beg = std::max(chunk.region.beg, chunk.ref_beg);
    const hts_pos_t end = std::min(chunk.region.end, chunk.ref_end);
    if (beg > end) return;

    const size_t offset = static_cast<size_t>(beg - chunk.ref_beg);
    const int len = static_cast<int>(end - beg + 1);
    int n = 0;
    uint64_t* intervals = sdust(
        nullptr,
        reinterpret_cast<const uint8_t*>(chunk.ref_seq.data() + offset),
        len,
        kSdustThreshold,
        kSdustWindow,
        &n);
    for (int i = 0; i < n; ++i) {
        const hts_pos_t rel_beg = static_cast<hts_pos_t>(intervals[i] >> 32);
        const hts_pos_t rel_end = static_cast<hts_pos_t>(static_cast<uint32_t>(intervals[i]));
        if (rel_end <= rel_beg) continue;
        chunk.low_complexity_regions.push_back(Interval{
            beg + rel_beg,
            beg + rel_end - 1,
            static_cast<int>(rel_end - rel_beg)});
    }
    std::free(intervals);
}

void populate_chunk_read_indexes(BamChunk& chunk) {
    chunk.ordered_read_ids.clear();
    chunk.up_overlap_read_ids.clear();
    chunk.down_overlap_read_ids.clear();
    chunk.noisy_regions.clear();

    chunk.ordered_read_ids.reserve(chunk.reads.size());
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        if (read.is_skipped) continue;
        chunk.ordered_read_ids.push_back(static_cast<int>(read_i));
        if (read.beg < chunk.region.beg) chunk.up_overlap_read_ids.push_back(static_cast<int>(read_i));
        if (read.end > chunk.region.end) chunk.down_overlap_read_ids.push_back(static_cast<int>(read_i));
        chunk.noisy_regions.insert(chunk.noisy_regions.end(), read.noisy_regions.begin(), read.noisy_regions.end());
    }
    merge_intervals(chunk.noisy_regions);
    sort_overlap_read_ids(chunk);
}

void finalize_bam_chunk(BamChunk& chunk, ReferenceCache& ref, const bam_hdr_t* header) {
    recompute_chunk_qual_stats(chunk);
    populate_reference_slice(chunk, ref, header);
    populate_low_complexity_intervals(chunk);
    populate_chunk_read_indexes(chunk);
}

void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants) {
    variants.clear();
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            variants.push_back(CandidateVariant{variant_key_from_digar(read.tid, digar), VariantCounts{}});
        }
    }

    collapse_fuzzy_large_insertions(variants);
}

void add_coverage(VariantCounts& counts, bool reverse, bool alt, bool low_quality) {
    if (low_quality) {
        counts.low_qual_cov++;
        return;
    }

    counts.total_cov++;
    if (alt) {
        counts.alt_cov++;
        reverse ? counts.reverse_alt++ : counts.forward_alt++;
    } else {
        counts.ref_cov++;
        reverse ? counts.reverse_ref++ : counts.forward_ref++;
    }
}

std::vector<ReadEvent>::const_iterator find_matching_event(const std::vector<ReadEvent>& events,
                                                           const VariantKey& key) {
    const hts_pos_t target_pos = key.sort_pos();
    auto event_it = std::lower_bound(
        events.begin(),
        events.end(),
        target_pos,
        [](const ReadEvent& event, hts_pos_t pos) {
            return event.key.sort_pos() < pos;
        });

    for (; event_it != events.end() && event_it->key.sort_pos() == target_pos; ++event_it) {
        if (same_candidate_site(key, event_it->key)) return event_it;
    }
    return events.end();
}

void infer_read_hap_phase_from_candidates(BamChunk& chunk, const CandidateTable& cand) {
    for (ReadRecord& read : chunk.reads) {
        if (read.is_skipped) continue;
        read.hap = 0;
        read.read_phase_set = 0;
        for (const ReadEvent& ev : read.events) {
            if (ev.low_quality) continue;
            for (const CandidateVariant& cv : cand) {
                if (!same_candidate_site(cv.key, ev.key)) continue;
                if (cv.counts.category != VariantCategory::CleanHetSnp &&
                    cv.counts.category != VariantCategory::CleanHetIndel) {
                    continue;
                }
                read.read_phase_set = cv.counts.phase_set;
                const auto ev_it = find_matching_event(read.events, cv.key);
                read.hap = (ev_it != read.events.end()) ? cv.counts.hap_alt : cv.counts.hap_ref;
                break;
            }
            if (read.hap != 0) break;
        }
    }
}

/**
 * longcallD `init_cand_vars_based_on_sites` + `update_cand_vars_from_digar` (bam_utils.c):
 * one pass over each read's digars, merged with the sorted `CandidateTable` using
 * `exact_comp_var_site_ins`; alt BQ = max(digar flag, mean qual < min_bq) like
 * is_low_qual || ave_qual < opt->min_bq; trailing sites with pos <= read end get
 * ref unless past pos_end.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        const RegionChunk* chunk_region,
                                        std::vector<ReadSupportRow>* read_support_out,
                                        int min_bq) {
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
            if (!is_variant_digar_for_cand_sweep(dig[digar_i])) {
                ++digar_i;
                continue;
            }
            const DigarOp& d = dig[digar_i];
            VariantKey dkey = variant_key_from_digar(read.tid, d);
            const int ave = get_digar_ave_qual(d, read.qual);
            const int ret = exact_comp_var_site_ins(
                &variants[static_cast<size_t>(site_i)].key, &dkey, kLongcalldMinSvLen);
            if (ret < 0) {
                VariantCounts& c = variants[static_cast<size_t>(site_i)].counts;
                add_coverage(c, read.reverse, false, false);
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
                    variants[static_cast<size_t>(site_i)].counts, read.reverse, true, lq);
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
                if (stid < read.tid) {
                    continue;
                }
                break;
            }
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) {
                break;
            }
            add_coverage(variants[static_cast<size_t>(site_i)].counts, read.reverse, false, false);
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

void assign_clean_phase_sets(CandidateTable& variants) {
    hts_pos_t current_phase_set = 0;
    for (CandidateVariant& candidate : variants) {
        const bool is_clean_het = candidate.counts.category == VariantCategory::CleanHetSnp ||
                                  candidate.counts.category == VariantCategory::CleanHetIndel;
        if (!is_clean_het) continue;

        if (current_phase_set == 0) current_phase_set = candidate.key.pos;
        candidate.counts.phase_set = current_phase_set;

        // This subcommand reports collection-time phase support, not final VCF genotypes.
        candidate.counts.hap_alt = candidate.counts.alt_cov >= candidate.counts.ref_cov ? 1 : 2;
        candidate.counts.hap_ref = candidate.counts.hap_alt == 1 ? 2 : 1;
    }
}

void resolve_simple_noisy_candidates(CandidateTable& variants) {
    for (CandidateVariant& candidate : variants) {
        if (candidate.counts.category != VariantCategory::NoisyCandidate &&
            candidate.counts.category != VariantCategory::RepeatHetIndel) {
            continue;
        }
        if (candidate.key.alt.size() >= static_cast<size_t>(kMinSvLen) || candidate.key.ref_len >= kMinSvLen) {
            candidate.counts.category = VariantCategory::NoisyResolved;
        }
    }
}

void assign_haplotype_refinement_light(CandidateTable& variants) {
    (void)variants;
}

void run_collect_var_pipeline(BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    (void)header;
    classify_cand_vars_pgphase(chunk, opts);
    assign_clean_phase_sets(chunk.candidates);
    resolve_simple_noisy_candidates(chunk.candidates);
    assign_haplotype_refinement_light(chunk.candidates);
}

static void load_and_prepare_chunk(BamChunk& chunk, const Options& opts, WorkerContext& context) {
    chunk.reads = load_read_records_for_chunk(
        opts, chunk.region, context.bam, context.header.get(), context.index.get(), context.ref);
    finalize_bam_chunk(chunk, context.ref, context.header.get());
}

static void collect_prephase_candidates(BamChunk& chunk,
                                        const Options& opts,
                                        std::vector<ReadSupportRow>* read_support_out) {
    pre_process_noisy_regs_pgphase(chunk, opts);
    collect_candidate_sites_from_records(chunk.region, chunk.reads, chunk.candidates);
    collect_allele_counts_from_records(
        chunk.reads, chunk.candidates, &chunk.region, read_support_out, opts.min_bq);
}

static void classify_and_filter_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    run_collect_var_pipeline(chunk, opts, header);
    post_process_noisy_regs_pgphase(chunk, chunk.candidates);
    if (!opts.is_ont) {
        apply_noisy_containment_filter(chunk);
    }
}

static void finalize_chunk_outputs(BamChunk& chunk, ChunkStitchBundle& stitch_out) {
    infer_read_hap_phase_from_candidates(chunk, chunk.candidates);
    fill_stitch_bundle(chunk, stitch_out);
    apply_noisy_wfa_touch(chunk);
}

CandidateTable collect_chunks_parallel(const Options& opts,
                                       const std::vector<RegionChunk>& chunks,
                                       std::vector<std::vector<ReadSupportRow>>* read_support_batches) {
    std::vector<CandidateTable> chunk_tables(chunks.size());
    std::vector<ChunkStitchBundle> stitch_info(chunks.size());
    if (chunks.empty()) return CandidateTable{};

    if (read_support_batches != nullptr) read_support_batches->assign(chunks.size(), {});

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), chunks.size());
    std::atomic<size_t> next_chunk{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&, worker_i]() {
            (void)worker_i;
            try {
                WorkerContext context(opts);

                while (true) {
                    const size_t chunk_i = next_chunk.fetch_add(1);
                    if (chunk_i >= chunks.size()) break;

                    BamChunk chunk;
                    chunk.region = chunks[chunk_i];
                    load_and_prepare_chunk(chunk, opts, context);
                    std::vector<ReadSupportRow>* rs_ptr = nullptr;
                    if (read_support_batches != nullptr) rs_ptr = &(*read_support_batches)[chunk_i];
                    collect_prephase_candidates(chunk, opts, rs_ptr);
                    classify_and_filter_candidates(chunk, opts, context.header.get());
                    finalize_chunk_outputs(chunk, stitch_info[chunk_i]);
                    chunk_tables[chunk_i] = std::move(chunk.candidates);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }

    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);

    for (size_t ii = 1; ii < chunks.size(); ++ii) {
        if (chunks[ii - 1].tid != chunks[ii].tid) continue;
        if (chunks[ii - 1].end + 1 != chunks[ii].beg) continue;
        apply_flip_between_chunks(
            stitch_info[ii - 1], stitch_info[ii], chunk_tables[ii - 1], chunk_tables[ii]);
    }

    CandidateTable merged;
    for (CandidateTable& table : chunk_tables) {
        merged.insert(
            merged.end(),
            std::make_move_iterator(table.begin()),
            std::make_move_iterator(table.end()));
    }
    collapse_fuzzy_large_insertions(merged);
    return merged;
}

std::vector<RegionChunk> load_region_chunks(const Options& opts) {
    SamFile bam(opts.bam_file, 1);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    std::unique_ptr<hts_idx_t, IndexDeleter> index(sam_index_load(bam.get(), opts.bam_file.c_str()));
    if (!index) throw std::runtime_error("region chunking requires an indexed BAM/CRAM");
    std::unique_ptr<faidx_t, FaiDeleter> fai(fai_load(opts.ref_fasta.c_str()));
    if (!fai) throw std::runtime_error("failed to load FASTA index for: " + opts.ref_fasta);
    return build_region_chunks(opts, header.get(), fai.get());
}

void write_read_support_tsv(const Options& opts, const std::vector<std::vector<ReadSupportRow>>& by_chunk) {
    SamFile bam(opts.bam_file, 1);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");

    std::ofstream out(opts.read_support_tsv);
    if (!out) throw std::runtime_error("failed to open read support output: " + opts.read_support_tsv);

    out << "CHROM\tPOS\tTYPE\tREF_LEN\tALT\tQNAME\tIS_ALT\tLOW_QUAL\tREVERSE\tMAPQ\tCHUNK_BEG\tCHUNK_END\n";

    for (const std::vector<ReadSupportRow>& batch : by_chunk) {
        for (const ReadSupportRow& r : batch) {
            const char* chrom = header->target_name[r.tid];
            out << chrom << '\t' << r.pos << '\t' << type_name(r.type) << '\t' << r.ref_len << '\t'
                << (r.alt.empty() ? "." : r.alt) << '\t' << r.qname << '\t' << r.is_alt << '\t' << r.is_low_qual
                << '\t' << (r.reverse ? 1 : 0) << '\t' << r.mapq << '\t' << r.chunk_beg << '\t' << r.chunk_end
                << '\n';
        }
    }
}

void write_variants(const Options& opts, faidx_t* fai, const CandidateTable& variants) {
    SamFile bam(opts.bam_file, 1);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai);

    std::ofstream out(opts.output_tsv);
    if (!out) throw std::runtime_error("failed to open output: " + opts.output_tsv);

    out << "CHROM\tPOS\tTYPE\tREF\tALT\tDP\tREF_COUNT\tALT_COUNT\tLOW_QUAL_COUNT"
        << "\tFORWARD_REF\tREVERSE_REF\tFORWARD_ALT\tREVERSE_ALT"
        << "\tAF\tCATEGORY\tPHASE_SET\tHAP_ALT\tHAP_REF\n";

    for (const CandidateVariant& candidate : variants) {
        const VariantKey& key = candidate.key;
        const VariantCounts& counts = candidate.counts;
        const std::string chrom = header->target_name[key.tid];
        std::string ref_seq = ".";
        std::string alt_seq = key.alt.empty() ? "." : key.alt;

        if (key.type == VariantType::Snp) {
            ref_seq = std::string(1, ref.base(key.tid, key.pos, header.get()));
        } else if (key.type == VariantType::Deletion) {
            ref_seq = ref.subseq(key.tid, key.pos, key.ref_len, header.get());
        }

        out << chrom << '\t' << key.pos << '\t' << type_name(key.type) << '\t' << ref_seq << '\t'
            << alt_seq << '\t' << counts.total_cov << '\t' << counts.ref_cov << '\t'
            << counts.alt_cov << '\t' << counts.low_qual_cov << '\t' << counts.forward_ref << '\t'
            << counts.reverse_ref << '\t' << counts.forward_alt << '\t' << counts.reverse_alt << '\t'
            << counts.allele_fraction << '\t' << category_name(counts.category) << '\t'
            << counts.phase_set << '\t' << counts.hap_alt << '\t' << counts.hap_ref << '\n';
    }
}

void write_variants_vcf(const Options& opts, faidx_t* fai, const CandidateTable& variants) {
    if (opts.output_vcf.empty()) return;

    SamFile bam(opts.bam_file, 1);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai);

    std::ofstream out(opts.output_vcf);
    if (!out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);

    out << "##fileformat=VCFv4.2\n";
    {
        std::time_t t = std::time(nullptr);
        std::tm* tm = std::localtime(&t);
        char date_buf[16] = {0};
        if (tm != nullptr && std::strftime(date_buf, sizeof(date_buf), "%Y%m%d", tm) > 0) {
            out << "##fileDate=" << date_buf << "\n";
        }
    }
    out << "##source=pgphase collect-bam-variation\n";
    out << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    out << "##FILTER=<ID=LowQual,Description=\"Low quality candidate\">\n";
    out << "##FILTER=<ID=RefCall,Description=\"Reference call candidate\">\n";
    out << "##FILTER=<ID=NoCall,Description=\"Site has depth=0 resulting in no call\">\n";
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
    out << "##INFO=<ID=CLEAN,Number=0,Type=Flag,Description=\"Clean-region candidate variant\">\n";
    out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    out << "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
    out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n";
    out << "##INFO=<ID=REFC,Number=1,Type=Integer,Description=\"Reference allele count\">\n";
    out << "##INFO=<ID=ALTC,Number=1,Type=Integer,Description=\"Alternate allele count\">\n";
    out << "##INFO=<ID=LQC,Number=1,Type=Integer,Description=\"Low-quality observation count\">\n";
    out << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Alternate allele fraction\">\n";
    out << "##INFO=<ID=CAT,Number=1,Type=String,Description=\"pgPhase candidate category\">\n";
    for (int32_t tid = 0; tid < header->n_targets; ++tid) {
        out << "##contig=<ID=" << header->target_name[tid] << ",length=" << header->target_len[tid] << ">\n";
    }
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    for (const CandidateVariant& candidate : variants) {
        const VariantKey& key = candidate.key;
        const VariantCounts& counts = candidate.counts;
        const std::string chrom = header->target_name[key.tid];

        hts_pos_t pos = key.pos;
        std::string ref_seq;
        std::string alt_seq;

        if (key.type == VariantType::Snp) {
            ref_seq = std::string(1, ref.base(key.tid, key.pos, header.get()));
            alt_seq = key.alt.empty() ? "." : key.alt;
        } else if (key.type == VariantType::Insertion) {
            const hts_pos_t anchor_pos = std::max<hts_pos_t>(1, key.pos - 1);
            pos = anchor_pos;
            const char anchor_base = ref.base(key.tid, anchor_pos, header.get());
            ref_seq = std::string(1, anchor_base);
            alt_seq = ref_seq + key.alt;
        } else { // Deletion
            const hts_pos_t anchor_pos = std::max<hts_pos_t>(1, key.pos - 1);
            pos = anchor_pos;
            const char anchor_base = ref.base(key.tid, anchor_pos, header.get());
            const std::string del_seq = ref.subseq(key.tid, key.pos, key.ref_len, header.get());
            ref_seq = std::string(1, anchor_base) + del_seq;
            alt_seq = std::string(1, anchor_base);
        }

        std::string filter = "PASS";
        if (counts.total_cov == 0) {
            filter = "NoCall";
        } else if (counts.category == VariantCategory::NonVariant) {
            filter = "RefCall";
        } else if (counts.category != VariantCategory::CleanHetSnp &&
                   counts.category != VariantCategory::CleanHetIndel &&
                   counts.category != VariantCategory::CleanHom) {
            filter = "LowQual";
        }

        const hts_pos_t end_pos = pos + static_cast<hts_pos_t>(ref_seq.size()) - 1;
        std::ostringstream info;
        info << "END=" << end_pos;
        if (counts.category == VariantCategory::CleanHetSnp || counts.category == VariantCategory::CleanHetIndel ||
            counts.category == VariantCategory::CleanHom) {
            info << ";CLEAN";
        }
        if (key.type == VariantType::Insertion || key.type == VariantType::Deletion) {
            const int svlen = (key.type == VariantType::Insertion) ? static_cast<int>(key.alt.size()) : -key.ref_len;
            if (std::abs(svlen) >= opts.min_sv_len) {
                info << ";SVTYPE=" << (svlen > 0 ? "INS" : "DEL");
                info << ";SVLEN=" << svlen;
            }
        }
        info << ";DP=" << counts.total_cov << ";REFC=" << counts.ref_cov << ";ALTC=" << counts.alt_cov
             << ";LQC=" << counts.low_qual_cov << ";AF=" << counts.allele_fraction
             << ";CAT=" << category_name(counts.category);

        out << chrom << '\t' << pos << "\t.\t" << ref_seq << '\t' << alt_seq << "\t.\t" << filter << '\t'
            << info.str() << '\n';
    }
}

void print_collect_help() {
    std::cout << "Usage: pgphase collect-bam-variation [options] <ref.fa> <input.bam>\n"
              << "Options:\n"
              << "  -t, --threads INT             Region worker threads [1]\n"
              << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
              << "  -B, --min-bq INT              Minimum base quality for candidate sites [10]\n"
              << "  -D, --min-depth INT           Minimum total depth for clean candidates [5]\n"
              << "      --min-alt-depth INT       Minimum alternate depth for clean candidates [2]\n"
              << "      --min-af FLOAT            Minimum allele fraction for clean het candidates [0.20]\n"
              << "      --max-af FLOAT            Maximum allele fraction for clean het candidates [0.80]\n"
              << "  -r, --region STR              Optional region; may be repeated\n"
              << "      --region-file FILE        BED file of regions to process\n"
              << "      --autosome                Process chr1-22 / 1-22 only\n"
              << "  -j, --max-var-ratio FLOAT     Skip reads above this variant/ref-span ratio [0.05]\n"
              << "      --max-noisy-frac FLOAT    Skip reads with > this fraction in noisy regions [0.5]\n"
              << "      --include-filtered        Include QC-fail and duplicate reads\n"
              << "  -o, --output FILE             Output TSV file [output.tsv]\n"
              << "  -v, --vcf-output FILE         Optional VCF output for collected candidates\n"
              << "      --read-support FILE       Per-read ref/alt observations at candidates (for phasing)\n"
              << "      --noisy-merge-dis INT     Max distance (bp) to merge noisy/SV windows [500]\n"
              << "      --min-sv-len INT          min_sv_len for noisy-region cgranges merge [30]\n"
              << "      --noisy-slide-win INT     Slide window (bp) for per-read noisy regions [HiFi 100 / ONT 25]\n"
              << "      --ont                     ONT mode: Fisher exact test for alt strand bias\n"
              << "      --strand-bias-pval FLOAT  max p-value for ONT strand filter [0.01]\n"
              << "      --noisy-max-xgaps INT     max indel len (bp) for STR/homopolymer flags [5]\n"
              << "\n"
              << "Regions can also be supplied after <input.bam>, e.g.\n"
              << "  pgphase collect-bam-variation ref.fa hifi.bam chr11:1000-2000 chr12:1-500\n";
}

} // namespace

int collect_bam_variation(int argc, char* argv[]) {
    Options opts;
    optind = 1;

    const struct option long_options[] = {
        {"threads", required_argument, nullptr, 't'},
        {"min-mapq", required_argument, nullptr, 'q'},
        {"min-bq", required_argument, nullptr, 'B'},
        {"min-depth", required_argument, nullptr, 'D'},
        {"min-alt-depth", required_argument, nullptr, kMinAltDepthOption},
        {"min-af", required_argument, nullptr, kMinAfOption},
        {"max-af", required_argument, nullptr, kMaxAfOption},
        {"region", required_argument, nullptr, 'r'},
        {"region-file", required_argument, nullptr, 'R'},
        {"autosome", no_argument, nullptr, 'a'},
        {"max-var-ratio", required_argument, nullptr, 'j'},
        {"max-noisy-frac", required_argument, nullptr, kMaxNoisyFracOption},
        {"include-filtered", no_argument, nullptr, 'f'},
        {"output", required_argument, nullptr, 'o'},
        {"vcf-output", required_argument, nullptr, 'v'},
        {"read-support", required_argument, nullptr, kReadSupportOption},
        {"noisy-merge-dis", required_argument, nullptr, kNoisyRegMergeDisOption},
        {"min-sv-len", required_argument, nullptr, kMinSvLenOption},
        {"noisy-slide-win", required_argument, nullptr, kNoisySlideWinOption},
        {"debug-site", required_argument, nullptr, kDebugSiteOption},
        {"ont", no_argument, nullptr, kOntOption},
        {"strand-bias-pval", required_argument, nullptr, kStrandBiasPvalOption},
        {"noisy-max-xgaps", required_argument, nullptr, kNoisyMaxXgapsOption},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "t:q:B:D:r:R:aj:o:v:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 't':
                opts.threads = std::stoi(optarg);
                break;
            case 'q':
                opts.min_mapq = std::stoi(optarg);
                break;
            case 'B':
                opts.min_bq = std::stoi(optarg);
                break;
            case 'D':
                opts.min_depth = std::stoi(optarg);
                break;
            case kMinAltDepthOption:
                opts.min_alt_depth = std::stoi(optarg);
                break;
            case kMinAfOption:
                opts.min_af = std::stod(optarg);
                break;
            case kMaxAfOption:
                opts.max_af = std::stod(optarg);
                break;
            case 'r':
                opts.regions.push_back(optarg);
                break;
            case 'R':
                opts.region_file = optarg;
                break;
            case 'a':
                opts.autosome = true;
                break;
            case 'j':
                opts.max_var_ratio_per_read = std::stod(optarg);
                break;
            case kMaxNoisyFracOption:
                opts.max_noisy_frac_per_read = std::stod(optarg);
                break;
            case 'f':
                opts.include_filtered = true;
                break;
            case 'o':
                opts.output_tsv = optarg;
                break;
            case 'v':
                opts.output_vcf = optarg;
                break;
            case kReadSupportOption:
                opts.read_support_tsv = optarg;
                break;
            case kNoisyRegMergeDisOption:
                opts.noisy_reg_merge_dis = std::stoi(optarg);
                break;
            case kMinSvLenOption:
                opts.min_sv_len = std::stoi(optarg);
                break;
            case kNoisySlideWinOption:
                opts.noisy_reg_slide_win = std::stoi(optarg);
                break;
            case kDebugSiteOption:
                opts.debug_site = optarg;
                break;
            case kOntOption:
                opts.is_ont = true;
                break;
            case kStrandBiasPvalOption:
                opts.strand_bias_pval = std::stod(optarg);
                break;
            case kNoisyMaxXgapsOption:
                opts.noisy_reg_max_xgaps = std::stoi(optarg);
                break;
            case 'h':
                print_collect_help();
                return 0;
            default:
                print_collect_help();
                return 1;
        }
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_bq < 0 ||
        opts.min_depth < 0 || opts.min_alt_depth < 0 || opts.noisy_reg_merge_dis < 0 || opts.min_sv_len < 0 ||
        opts.min_af < 0.0 || opts.max_af < opts.min_af || opts.strand_bias_pval < 0.0 ||
        opts.strand_bias_pval > 1.0 || opts.max_var_ratio_per_read < 0.0 || opts.max_noisy_frac_per_read < 0.0 ||
        opts.noisy_reg_max_xgaps < 0 || opts.noisy_reg_slide_win < -1) {
        std::cerr << "Error: numeric thresholds are invalid\n";
        return 1;
    }

    if (optind + 2 > argc) {
        print_collect_help();
        return 1;
    }
    opts.ref_fasta = argv[optind];
    opts.bam_file = argv[optind + 1];
    for (int arg_i = optind + 2; arg_i < argc; ++arg_i) {
        opts.regions.push_back(argv[arg_i]);
    }

    try {
        std::unique_ptr<faidx_t, FaiDeleter> fai(fai_load(opts.ref_fasta.c_str()));
        if (!fai) throw std::runtime_error("failed to load FASTA index for: " + opts.ref_fasta);

        const std::vector<RegionChunk> chunks = load_region_chunks(opts);
        std::vector<std::vector<ReadSupportRow>> read_support_batches;
        std::vector<std::vector<ReadSupportRow>>* rs_batches_ptr =
            opts.read_support_tsv.empty() ? nullptr : &read_support_batches;
        CandidateTable variants = collect_chunks_parallel(opts, chunks, rs_batches_ptr);
        write_variants(opts, fai.get(), variants);
        write_variants_vcf(opts, fai.get(), variants);
        if (!opts.read_support_tsv.empty()) {
            write_read_support_tsv(opts, read_support_batches);
        }

        std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
                  << " worker thread(s)\n";
        std::cerr << "Collected " << variants.size() << " candidate variant sites into "
                  << opts.output_tsv << "\n";
        if (!opts.output_vcf.empty()) {
            std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
        }
        if (!opts.read_support_tsv.empty()) {
            size_t n_rows = 0;
            for (const auto& b : read_support_batches) n_rows += b.size();
            std::cerr << "Wrote " << n_rows << " read x candidate observations to " << opts.read_support_tsv
                      << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
