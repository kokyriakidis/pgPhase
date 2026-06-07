#ifndef PGPHASE_BAM_IO_TYPES_HPP
#define PGPHASE_BAM_IO_TYPES_HPP

// BAM/CRAM I/O infrastructure: RAII deleters for htslib handles, SamFile
// wrapper, ReferenceCache, and per-thread WorkerContext.  Separated from
// phasing_types.hpp so that graph-only code can avoid pulling in faidx.h.

#include "phasing_types.hpp"

#include <cctype>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/sam.h>

namespace pgphase_collect {

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
