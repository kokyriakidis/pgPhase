/**
 * @file align.cpp
 * @brief WFA2-lib / abPOA alignment wrappers and noisy-region helpers.
 *
 * Mirrors longcallD align.c, limited to the functions needed for step 4
 * (iterative noisy-region MSA variant calling).  Functions that touch TE
 * detection, BAM re-writing, or somatic calling are not ported.
 */

#include "align.hpp"
#include "bam_digar.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>

// WFA2-lib (C library — needs extern "C" wrapper since headers lack __cplusplus guards).
extern "C" {
#include "wavefront/wavefront_align.h"
#include "alignment/cigar.h"
}

// longcallD uses edlib_xgaps as part of noisy-region sampling.
// Ported here for strict parity.
#include "edlib.h"

// abPOA
#include "abpoa.h"

// HTSlib
#include <htslib/sam.h>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Helpers
// ════════════════════════════════════════════════════════════════════════════

// seq_nt16_int[] maps bam_seqi nibble → 2-bit base (0=A, 1=C, 2=G, 3=T, 4=N).
// Matches longcallD's `seq_nt16_int` from seq.h.
static const uint8_t kSeqNt16Int[16] = {4,0,1,4, 2,4,4,4, 3,4,4,4, 4,4,4,4};

static inline int add_phase_set(hts_pos_t ps,
                                std::vector<hts_pos_t>& uniq,
                                int& n_uniq) {
    for (int i = 0; i < n_uniq; ++i) {
        if (uniq[static_cast<size_t>(i)] == ps) return i;
    }
    uniq[static_cast<size_t>(n_uniq)] = ps;
    return n_uniq++;
}

static bool is_homopolymer(const uint8_t* seq, int seq_len, int flank_len,
                           int* hp_start, int* hp_end, int* hp_len) {
    if (seq_len < 2 * flank_len || seq_len > 2 * flank_len + 50) return false;
    constexpr int min_hp = 5;
    *hp_start = -1; *hp_end = -1; *hp_len = 0;
    for (int i = flank_len - 1; i < seq_len - flank_len + 1; ++i) {
        if (seq[i] == seq[i - 1]) {
            if (*hp_start == -1) *hp_start = i - 2;
            (*hp_len)++;
        } else {
            if (*hp_len >= min_hp) { *hp_end = i; return true; }
            *hp_start = -1; *hp_len = 0;
        }
    }
    if (*hp_len >= min_hp) { *hp_end = *hp_start + *hp_len + 1; return true; }
    return false;
}

// Build alignment strings from WFA cigar — mirrors longcallD wfa_collect_pretty_alignment.
static int wfa_collect_pretty_alignment(cigar_t* cigar,
                                        const uint8_t* pattern, int plen,
                                        const uint8_t* text, int tlen,
                                        std::vector<uint8_t>& pat_alg,
                                        std::vector<uint8_t>& txt_alg) {
    const int max_buf = plen + tlen + 1;
    pat_alg.assign(static_cast<size_t>(max_buf), 0);
    txt_alg.assign(static_cast<size_t>(max_buf), 0);

    const char* ops   = cigar->operations;
    const int beg_off = cigar->begin_offset;
    const int end_off = cigar->end_offset;

    int alg_pos = 0, pp = 0, tp = 0;
    for (int i = beg_off; i < end_off; ++i) {
        switch (ops[i]) {
            case 'M':
            case 'X':
                pat_alg[static_cast<size_t>(alg_pos)] = pattern[pp++];
                txt_alg[static_cast<size_t>(alg_pos)] = text[tp++];
                ++alg_pos;
                break;
            case 'I':
                pat_alg[static_cast<size_t>(alg_pos)] = 5; // gap
                txt_alg[static_cast<size_t>(alg_pos)] = text[tp++];
                ++alg_pos;
                break;
            case 'D':
                pat_alg[static_cast<size_t>(alg_pos)] = pattern[pp++];
                txt_alg[static_cast<size_t>(alg_pos)] = 5; // gap
                ++alg_pos;
                break;
            default:
                break;
        }
    }
    pat_alg.resize(static_cast<size_t>(alg_pos));
    txt_alg.resize(static_cast<size_t>(alg_pos));
    return alg_pos;
}

// ════════════════════════════════════════════════════════════════════════════
// calc_read_error_rate  (mirrors longcallD seq.c)
// ════════════════════════════════════════════════════════════════════════════

double calc_read_error_rate(int len, const uint8_t* qual) {
    if (len <= 0 || qual == nullptr) return 0.0;
    double expected_errors = 0.0;
    for (int i = 0; i < len; ++i)
        expected_errors += std::pow(10.0, -static_cast<double>(qual[i]) / 10.0);
    return expected_errors / len;
}

// ════════════════════════════════════════════════════════════════════════════
// full_cover_cmp  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int full_cover_cmp(int c1, int c2) {
    if (c1 == c2) return 0;
    if (noisyIsBothCover(c1)) return 1;
    if (noisyIsBothCover(c2)) return -1;
    if (noisyIsLeftCover(c1) && noisyIsLeftCover(c2)) return 0;
    if (noisyIsRightCover(c1) && noisyIsRightCover(c2)) return 0;
    return c1 - c2;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_read_info  (mirrors longcallD align.c lines 1377–1461)
// ════════════════════════════════════════════════════════════════════════════

NoisyReadInfo collect_noisy_read_info(const Options& opts, const BamChunk& chunk,
                                      hts_pos_t reg_beg, hts_pos_t reg_end,
                                      const std::vector<int>& noisy_read_ids) {
    const int n = static_cast<int>(noisy_read_ids.size());
    NoisyReadInfo info;
    info.n_reads        = n;
    info.noisy_read_ids = noisy_read_ids;
    info.lens           .resize(static_cast<size_t>(n), 0);
    info.seqs           .resize(static_cast<size_t>(n));
    info.quals          .resize(static_cast<size_t>(n));
    info.strands        .resize(static_cast<size_t>(n), 0);
    info.fully_covers   .resize(static_cast<size_t>(n), 0);
    info.haps           .resize(static_cast<size_t>(n), 0);
    info.phase_sets     .resize(static_cast<size_t>(n), 0);
    info.read_id_to_full_covers.assign(chunk.reads.size(), 0);
    info.read_reg_beg          .assign(chunk.reads.size(), 0);
    info.read_reg_end          .assign(chunk.reads.size(), 0);

    for (int i = 0; i < n; ++i) {
        const int read_id = noisy_read_ids[static_cast<size_t>(i)];
        const ReadRecord& r = chunk.reads[static_cast<size_t>(read_id)];
        const auto& digars  = r.digars;
        const int n_digar   = static_cast<int>(digars.size());
        const bam1_t* aln   = r.alignment.get();
        const uint8_t* bseq = bam_get_seq(aln);

        info.strands[static_cast<size_t>(i)] = r.reverse ? 1 : 0;
        info.haps   [static_cast<size_t>(i)] = chunk.haps[static_cast<size_t>(read_id)];
        info.phase_sets[static_cast<size_t>(i)] = chunk.phase_sets[static_cast<size_t>(read_id)];

        // Match longcallD collect_noisy_read_info exactly:
        // reg_read_end starts from digar2qlen(read_digars)-1, then hard-clips adjust both ends.
        int reg_read_beg = 0;
        int reg_read_end = 0;
        {
            int qlen = 0;
            if (n_digar > 0) {
                const DigarOp& last = digars[static_cast<size_t>(n_digar - 1)];
                qlen = last.qi;
                if (last.type == DigarType::Equal || last.type == DigarType::Snp ||
                    last.type == DigarType::Insertion || last.type == DigarType::SoftClip ||
                    last.type == DigarType::HardClip) {
                    qlen += last.len;
                }
            }
            reg_read_end = qlen - 1;
        }
        if (n_digar > 0 && digars[0].type == DigarType::HardClip) {
            reg_read_beg = digars[0].len;
        }
        if (n_digar > 0 && digars[static_cast<size_t>(n_digar - 1)].type == DigarType::HardClip) {
            reg_read_end = digars[static_cast<size_t>(n_digar - 1)].qi - 1;
        }

        hts_pos_t reg_digar_beg = -1, reg_digar_end = -1;
        int beg_is_del = 0, end_is_del = 0;

        for (int d = 0; d < n_digar; ++d) {
            const DigarOp& dop = digars[static_cast<size_t>(d)];
            if (dop.type == DigarType::SoftClip || dop.type == DigarType::HardClip) continue;

            hts_pos_t digar_beg = dop.pos;
            hts_pos_t digar_end;
            if (dop.type == DigarType::Snp || dop.type == DigarType::Equal ||
                dop.type == DigarType::Deletion) {
                digar_end = digar_beg + static_cast<hts_pos_t>(dop.len) - 1;
            } else {
                digar_end = digar_beg; // INS / RefSkip span only one ref pos
            }

            if (digar_beg > reg_end) break;
            if (digar_end < reg_beg) continue;

            if (digar_beg <= reg_beg && digar_end >= reg_beg) {
                if (dop.type == DigarType::Deletion) {
                    reg_digar_beg = reg_beg;
                    reg_read_beg  = dop.qi; // qi is on the right side of the DEL
                    if (dop.len > opts.noisy_reg_flank_len) beg_is_del = 1;
                } else {
                    reg_digar_beg = reg_beg;
                    reg_read_beg  = dop.qi + static_cast<int>(reg_beg - digar_beg);
                }
            }
            if (digar_beg <= reg_end && digar_end >= reg_end) {
                if (dop.type == DigarType::Deletion) {
                    reg_digar_end = reg_end;
                    reg_read_end  = dop.qi - 1; // qi is left side of the DEL
                    if (dop.len > opts.noisy_reg_flank_len) end_is_del = 1;
                } else {
                    reg_digar_end = reg_end;
                    reg_read_end  = dop.qi + static_cast<int>(reg_end - digar_beg);
                }
            }
        }

        int cover = kNoisyNoCover;
        if (reg_digar_beg == reg_beg && reg_digar_end == reg_end) {
            if (!beg_is_del && !end_is_del)      cover = kNoisyLeftCover | kNoisyRightCover;
            else if (!beg_is_del && end_is_del)  cover = kNoisyLeftCover | kNoisyRightGap;
            else if (beg_is_del && !end_is_del)  cover = kNoisyLeftGap   | kNoisyRightCover;
            else                                 cover = kNoisyLeftGap   | kNoisyRightGap;
        } else if (reg_digar_beg == reg_beg) {
            cover = beg_is_del ? kNoisyLeftGap : kNoisyLeftCover;
        } else if (reg_digar_end == reg_end) {
            cover = end_is_del ? kNoisyRightGap : kNoisyRightCover;
        }

        if (reg_read_end < reg_read_beg) { reg_read_end = reg_read_beg; }

        const int rlen = reg_read_end - reg_read_beg + 1;
        info.lens[static_cast<size_t>(i)] = rlen;
        info.seqs [static_cast<size_t>(i)].resize(static_cast<size_t>(rlen));
        info.quals[static_cast<size_t>(i)].resize(static_cast<size_t>(rlen));

        for (int j = reg_read_beg; j <= reg_read_end; ++j) {
            const int ji = j - reg_read_beg;
            info.seqs [static_cast<size_t>(i)][static_cast<size_t>(ji)] =
                kSeqNt16Int[bam_seqi(bseq, j)];
            info.quals[static_cast<size_t>(i)][static_cast<size_t>(ji)] =
                r.qual[static_cast<size_t>(j)];
        }

        info.fully_covers[static_cast<size_t>(i)] = cover;
        info.read_id_to_full_covers[static_cast<size_t>(read_id)] = cover;
        info.read_reg_beg          [static_cast<size_t>(read_id)] = reg_read_beg;
        info.read_reg_end          [static_cast<size_t>(read_id)] = reg_read_end;
        if (opts.verbose >= 2) {
            std::fprintf(stderr,
                         "NoisyReadInfo\tread_id=%d\treg=%" PRId64 "-%" PRId64
                         "\tread_q=%d-%d\tlen=%d\tcover=%d\tseq=",
                         read_id,
                         static_cast<int64_t>(reg_beg),
                         static_cast<int64_t>(reg_end),
                         reg_read_beg,
                         reg_read_end,
                         rlen,
                         cover);
            const int dump_len = std::min(rlen, 24);
            for (int j = 0; j < dump_len; ++j) {
                std::fprintf(stderr, "%c",
                             "ACGTN"[info.seqs[static_cast<size_t>(i)][static_cast<size_t>(j)] < 5
                                 ? info.seqs[static_cast<size_t>(i)][static_cast<size_t>(j)]
                                 : 4]);
            }
            std::fprintf(stderr, "\n");
        }
    }
    return info;
}

static inline int digar_type_to_bam_op(DigarType t) {
    switch (t) {
        case DigarType::Equal: return BAM_CEQUAL;
        case DigarType::Snp: return BAM_CDIFF;
        case DigarType::Insertion: return BAM_CINS;
        case DigarType::Deletion: return BAM_CDEL;
        case DigarType::SoftClip: return BAM_CSOFT_CLIP;
        case DigarType::HardClip: return BAM_CHARD_CLIP;
        case DigarType::RefSkip: return BAM_CREF_SKIP;
    }
    return BAM_CEQUAL;
}

static inline DigarType bam_op_to_digar_type(int op) {
    switch (op) {
        case BAM_CEQUAL: return DigarType::Equal;
        case BAM_CDIFF: return DigarType::Snp;
        case BAM_CINS: return DigarType::Insertion;
        case BAM_CDEL: return DigarType::Deletion;
        case BAM_CSOFT_CLIP: return DigarType::SoftClip;
        case BAM_CHARD_CLIP: return DigarType::HardClip;
        case BAM_CREF_SKIP: return DigarType::RefSkip;
        default: return DigarType::Equal;
    }
}

static inline bool same_digar_for_merge(const DigarOp& a, const DigarOp& b) {
    const int aop = digar_type_to_bam_op(a.type);
    if (aop != BAM_CEQUAL && aop != BAM_CINS && aop != BAM_CDEL) return false;
    return a.type == b.type && a.low_quality == b.low_quality;
}

static inline void push_digar0(std::vector<DigarOp>& out, const DigarOp& d) {
    if (d.len <= 0) return;
    if (out.empty() || !same_digar_for_merge(out.back(), d)) {
        out.push_back(d);
        return;
    }
    if (out.back().type == DigarType::Insertion) {
        out.back().alt += d.alt;
    }
    out.back().len += d.len;
}

static inline void push_digar_alt_seq(std::vector<DigarOp>& out, const DigarOp& d) {
    if (d.len <= 0) return;
    if (out.empty() || !same_digar_for_merge(out.back(), d)) {
        out.push_back(d);
        return;
    }
    if (out.back().type == DigarType::Insertion) {
        out.back().alt += d.alt;
    }
    out.back().len += d.len;
}

/** longcallD `bam_utils.h` `double_check_digar`: returns 1 if invalid qi chain, else 0. */
static inline int double_check_digar(const std::vector<DigarOp>& digars) {
    const int n_digar = static_cast<int>(digars.size());
    if (n_digar == 0) return 0;
    for (int i = n_digar - 1; i > 0; --i) {
        const DigarOp& prev = digars[static_cast<size_t>(i - 1)];
        const int prev_op = digar_type_to_bam_op(prev.type);
        int qi = prev.qi;
        if (prev_op == BAM_CEQUAL || prev_op == BAM_CMATCH || prev_op == BAM_CDIFF ||
            prev_op == BAM_CINS || prev_op == BAM_CSOFT_CLIP || prev_op == BAM_CHARD_CLIP) {
            qi = prev.qi + prev.len;
        }
        if (qi != digars[static_cast<size_t>(i)].qi) return 1;
    }
    return 0;
}

/** longcallD `bam_utils.h` `digar2qlen` on the pre-MSA digar stream (uses last op only). */
static inline int digar2qlen_lcd(const std::vector<DigarOp>& digars) {
    if (digars.empty()) return 0;
    const DigarOp& last = digars.back();
    const int op = digar_type_to_bam_op(last.type);
    int qlen = last.qi;
    if (op == BAM_CEQUAL || op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CINS ||
        op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
        qlen += last.len;
    }
    return qlen;
}

static inline int digar_q_end(const DigarOp& d) {
    const int op = digar_type_to_bam_op(d.type);
    if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CINS) return d.qi + d.len - 1;
    return d.qi;
}

static inline hts_pos_t digar_ref_end(const DigarOp& d) {
    const int op = digar_type_to_bam_op(d.type);
    if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) return d.pos + d.len - 1;
    return d.pos;
}

static std::vector<DigarOp> collect_left_digars(const std::vector<DigarOp>& digars,
                                                 int read_noisy_beg,
                                                 hts_pos_t ref_noisy_beg) {
    std::vector<DigarOp> left;
    for (size_t i = 0; i < digars.size(); ++i) {
        const DigarOp& d = digars[i];
        const int op = digar_type_to_bam_op(d.type);
        if (i == 0 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            push_digar0(left, d);
            continue;
        }
        const int dq_end = digar_q_end(d);
        const hts_pos_t dr_end = digar_ref_end(d);
        if (d.qi >= read_noisy_beg && d.pos >= ref_noisy_beg) break;
        if (dq_end < read_noisy_beg && dr_end < ref_noisy_beg) {
            push_digar0(left, d);
        } else if (dq_end >= read_noisy_beg || dr_end >= ref_noisy_beg) {
            DigarOp c = d;
            if (op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
                c.len = read_noisy_beg - c.qi;
            } else if (op == BAM_CDEL) {
                c.len = static_cast<int>(ref_noisy_beg - c.pos);
            }
            push_digar0(left, c);
            break;
        } else {
            std::cerr << "Error: collect_left_digars: digar_ref_end: " << dr_end
                      << " ref_noisy_beg: " << ref_noisy_beg << "\n";
            std::abort();
        }
    }
    return left;
}

static std::vector<DigarOp> collect_right_digars(const std::vector<DigarOp>& digars,
                                                  int read_noisy_end,
                                                  hts_pos_t ref_noisy_end) {
    std::vector<DigarOp> right;
    for (size_t i = 0; i < digars.size(); ++i) {
        const DigarOp& d = digars[i];
        const int op = digar_type_to_bam_op(d.type);
        if (i + 1 == digars.size() && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            push_digar0(right, d);
            continue;
        }
        const int dq_end = digar_q_end(d);
        const hts_pos_t dr_end = digar_ref_end(d);
        if (dq_end <= read_noisy_end && dr_end <= ref_noisy_end) continue;
        if (d.qi > read_noisy_end && d.pos > ref_noisy_end) {
            push_digar0(right, d);
        } else if (d.qi <= read_noisy_end || d.pos <= ref_noisy_end) {
            DigarOp c = d;
            if (op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
                c.len = dq_end - read_noisy_end;
                c.qi = read_noisy_end + 1;
                if (op == BAM_CINS && !d.alt.empty()) {
                    const int shift = read_noisy_end + 1 - d.qi;
                    if (shift >= 0 && shift < static_cast<int>(d.alt.size())) {
                        c.alt = d.alt.substr(static_cast<size_t>(shift));
                    } else {
                        c.alt.clear();
                    }
                } else if (op != BAM_CINS) {
                    c.pos = ref_noisy_end + 1;
                }
            } else if (op == BAM_CDEL) {
                c.len = static_cast<int>(dr_end - ref_noisy_end);
                c.pos = ref_noisy_end + 1;
            }
            push_digar0(right, c);
        } else {
            std::cerr << "Error: collect_right_digars: digar_ref_end: " << dr_end
                      << " ref_noisy_end: " << ref_noisy_end << "\n";
            std::abort();
        }
    }
    return right;
}

// longcallD align.c `collect_full_msa_digars` (read_end unused; kept for parity with caller signature).
static std::vector<DigarOp> collect_full_msa_digars(int read_beg, int /*read_end*/, hts_pos_t ref_beg,
                                                     const AlnStr& ref_read) {
    std::vector<DigarOp> out;
    const int msa_len = ref_read.aln_len;
    if (msa_len <= 0) return out;
    int read_pos = read_beg;
    hts_pos_t ref_pos = ref_beg;
    const int left_read_start = 0;
    const int right_read_end = msa_len - 1;
    for (int i = 0; i < msa_len; ++i) {
        const uint8_t rb = ref_read.target_aln[static_cast<size_t>(i)];
        const uint8_t qb = ref_read.query_aln [static_cast<size_t>(i)];
        if (rb == 5 && qb == 5) continue;
        if (qb != 5 && rb != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                if (qb == rb) {
                    push_digar0(out, DigarOp{ref_pos, DigarType::Equal, 1, read_pos, false, ""});
                } else {
                    DigarOp d{ref_pos, DigarType::Snp, 1, read_pos, false,
                              std::string(1, "ACGTN"[qb < 5 ? qb : 4])};
                    push_digar0(out, d);
                }
            }
            ++read_pos;
            ++ref_pos;
        } else if (qb != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                push_digar0(out, DigarOp{ref_pos, DigarType::Insertion, 1, read_pos, false,
                                         std::string(1, "ACGTN"[qb < 5 ? qb : 4])});
            }
            ++read_pos;
        } else {
            if (i >= left_read_start && i <= right_read_end) {
                push_digar0(out, DigarOp{ref_pos, DigarType::Deletion, 1, read_pos, false, ""});
            }
            ++ref_pos;
        }
    }
    return out;
}

// longcallD align.c `collect_left_msa_digars` (read_end overwritten immediately in LCD; unused).
static std::vector<DigarOp> collect_left_msa_digars(int read_beg, int /*read_end*/, int qlen,
                                                     hts_pos_t ref_beg, const AlnStr& ref_read) {
    std::vector<DigarOp> out;
    const int msa_len = ref_read.aln_len;
    if (msa_len <= 0) return out;
    int read_pos = read_beg;
    int read_end_pos = read_pos - 1;
    hts_pos_t ref_pos = ref_beg;
    const int left_read_start = 0;
    int right_read_end = msa_len - 1;
    int right_skipped_read_base = 0;
    int is_covered_by_ref = 0;
    for (int i = msa_len - 1; i >= 0; --i) {
        const uint8_t rb = ref_read.target_aln[static_cast<size_t>(i)];
        const uint8_t qb = ref_read.query_aln [static_cast<size_t>(i)];
        if (rb != 5) is_covered_by_ref = 1;
        if (is_covered_by_ref && qb != 5) {
            right_read_end = i;
            break;
        }
        if (!is_covered_by_ref && qb != 5) ++right_skipped_read_base;
    }
    for (int i = 0; i < msa_len; ++i) {
        if (ref_read.query_aln[static_cast<size_t>(i)] != 5) ++read_end_pos;
    }
    for (int i = 0; i < msa_len; ++i) {
        const uint8_t rb = ref_read.target_aln[static_cast<size_t>(i)];
        const uint8_t qb = ref_read.query_aln [static_cast<size_t>(i)];
        if (rb == 5 && qb == 5) continue;
        if (qb != 5 && rb != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                if (qb == rb) {
                    push_digar0(out, DigarOp{ref_pos, DigarType::Equal, 1, read_pos, false, ""});
                } else {
                    DigarOp d{ref_pos, DigarType::Snp, 1, read_pos, false,
                              std::string(1, "ACGTN"[qb < 5 ? qb : 4])};
                    push_digar0(out, d);
                }
            }
            ++read_pos;
            ++ref_pos;
        } else if (qb != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                push_digar0(out, DigarOp{ref_pos, DigarType::Insertion, 1, read_pos, false,
                                         std::string(1, "ACGTN"[qb < 5 ? qb : 4])});
            }
            ++read_pos;
        } else {
            if (i >= left_read_start && i <= right_read_end) {
                push_digar0(out, DigarOp{ref_pos, DigarType::Deletion, 1, read_pos, false, ""});
            }
            ++ref_pos;
        }
    }
    if (read_end_pos < qlen - 1 || right_skipped_read_base > 0) {
        push_digar0(out, DigarOp{ref_pos, DigarType::SoftClip,
                                 qlen - 1 - read_end_pos + right_skipped_read_base, read_end_pos + 1, false,
                                 ""});
    }
    return out;
}

// longcallD align.c `collect_right_msa_digars` (read_beg unused in LCD).
static std::vector<DigarOp> collect_right_msa_digars(int /*read_beg*/, int read_end, hts_pos_t ref_beg,
                                                      hts_pos_t ref_end, const AlnStr& ref_read) {
    std::vector<DigarOp> out;
    const int msa_len = ref_read.aln_len;
    if (msa_len <= 0) return out;
    int read_pos = read_end + 1;
    hts_pos_t _ref_pos = ref_end + 1;
    hts_pos_t ref_pos = ref_beg;
    int left_read_start = 0;
    const int right_read_end = msa_len - 1;
    int left_skipped_read_base = 0;
    int is_covered_by_ref = 0;
    for (int i = 0; i < msa_len; ++i) {
        const uint8_t rb = ref_read.target_aln[static_cast<size_t>(i)];
        const uint8_t qb = ref_read.query_aln [static_cast<size_t>(i)];
        if (rb != 5) is_covered_by_ref = 1;
        if (is_covered_by_ref && qb != 5) {
            left_read_start = i;
            break;
        }
        if (!is_covered_by_ref && qb != 5) ++left_skipped_read_base;
    }
    for (int i = msa_len - 1; i >= 0; --i) {
        if (ref_read.target_aln[static_cast<size_t>(i)] != 5) --_ref_pos;
        if (ref_read.query_aln[static_cast<size_t>(i)] != 5) {
            --read_pos;
            ref_pos = _ref_pos;
        }
    }
    if (read_pos > 0 || left_skipped_read_base > 0) {
        push_digar0(out,
                    DigarOp{ref_pos, DigarType::SoftClip, read_pos + left_skipped_read_base, 0, false, ""});
    }
    read_pos += left_skipped_read_base;
    for (int i = left_read_start; i <= right_read_end; ++i) {
        const uint8_t rb = ref_read.target_aln[static_cast<size_t>(i)];
        const uint8_t qb = ref_read.query_aln [static_cast<size_t>(i)];
        if (rb == 5 && qb == 5) continue;
        if (qb != 5 && rb != 5) {
            if (qb == rb) {
                push_digar0(out, DigarOp{ref_pos, DigarType::Equal, 1, read_pos, false, ""});
            } else {
                DigarOp d{ref_pos, DigarType::Snp, 1, read_pos, false,
                          std::string(1, "ACGTN"[qb < 5 ? qb : 4])};
                push_digar0(out, d);
            }
            ++read_pos;
            ++ref_pos;
        } else if (qb != 5) {
            push_digar0(out, DigarOp{ref_pos, DigarType::Insertion, 1, read_pos, false,
                                     std::string(1, "ACGTN"[qb < 5 ? qb : 4])});
            ++read_pos;
        } else {
            push_digar0(out, DigarOp{ref_pos, DigarType::Deletion, 1, read_pos, false, ""});
            ++ref_pos;
        }
    }
    return out;
}

static void update_digars_from_msa1(ReadRecord& read, const AlnStr& ref_read_aln, int full_cover,
                                    hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                                    int read_beg, int read_end) {
    if (noisyIsNotCover(full_cover)) return;
    std::vector<DigarOp> left, right, msa, merged;
    if (noisyIsBothCover(full_cover) ||
        (noisyIsLeftCover(full_cover) && noisyIsRightGap(full_cover)) ||
        (noisyIsRightCover(full_cover) && noisyIsLeftGap(full_cover))) {
        left = collect_left_digars(read.digars, read_beg, noisy_reg_beg);
        right = collect_right_digars(read.digars, read_end, noisy_reg_end);
        msa = collect_full_msa_digars(read_beg, read_end, noisy_reg_beg, ref_read_aln);
        for (const DigarOp& d : left) push_digar_alt_seq(merged, d);
        for (const DigarOp& d : msa) push_digar0(merged, d);
        for (const DigarOp& d : right) push_digar_alt_seq(merged, d);
    } else if (noisyIsLeftCover(full_cover)) {
        left = collect_left_digars(read.digars, read_beg, noisy_reg_beg);
        msa = collect_left_msa_digars(read_beg, read_end, digar2qlen_lcd(read.digars), noisy_reg_beg, ref_read_aln);
        for (const DigarOp& d : left) push_digar_alt_seq(merged, d);
        for (const DigarOp& d : msa) push_digar0(merged, d);
    } else if (noisyIsRightCover(full_cover)) {
        right = collect_right_digars(read.digars, read_end, noisy_reg_end);
        msa = collect_right_msa_digars(read_beg, read_end, noisy_reg_beg, noisy_reg_end, ref_read_aln);
        for (const DigarOp& d : msa) push_digar0(merged, d);
        for (const DigarOp& d : right) push_digar_alt_seq(merged, d);
    }
    // longcallD align.c `update_digars_from_msa1`: replace digars only if `double_check_digar` returns 0.
    if (double_check_digar(merged)) {
        return;
    }
    read.digars = std::move(merged);
}

static void update_digars_from_aln_str(BamChunk& chunk,
                                       hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                                       const std::vector<int>& read_id_to_full_covers,
                                       const std::vector<int>& read_reg_beg,
                                       const std::vector<int>& read_reg_end,
                                       const std::array<std::vector<AlnStr>, 2>& aln_strs,
                                       int n_cons,
                                       const std::array<int, 2>& clu_n_seqs,
                                       const std::array<std::vector<int>, 2>& clu_read_ids) {
    for (int ci = 0; ci < n_cons; ++ci) {
        const int n_seqs = clu_n_seqs[static_cast<size_t>(ci)];
        for (int read_i = 0; read_i < n_seqs; ++read_i) {
            const int read_id = clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(read_i)];
            const int read_beg = read_reg_beg[static_cast<size_t>(read_id)];
            const int read_end = read_reg_end[static_cast<size_t>(read_id)];
            const size_t ref_read_slot = static_cast<size_t>((read_i + 1) * 2);
            if (ref_read_slot >= aln_strs[static_cast<size_t>(ci)].size()) continue;
            const AlnStr& ref_read = aln_strs[static_cast<size_t>(ci)][ref_read_slot];
            update_digars_from_msa1(chunk.reads[static_cast<size_t>(read_id)], ref_read,
                                    read_id_to_full_covers[static_cast<size_t>(read_id)],
                                    noisy_reg_beg, noisy_reg_end, read_beg, read_end);
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// sort_noisy_region_reads  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

void sort_noisy_region_reads(NoisyReadInfo& info, bool use_error_rate) {
    const int n = info.n_reads;
    std::vector<double> err(static_cast<size_t>(n), 0.0);
    if (use_error_rate) {
        for (int i = 0; i < n; ++i) {
            const auto& q = info.quals[static_cast<size_t>(i)];
            err[static_cast<size_t>(i)] = calc_read_error_rate(
                info.lens[static_cast<size_t>(i)], q.data());
        }
    }
    // bubble sort matching longcallD exactly
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int cc = full_cover_cmp(info.fully_covers[static_cast<size_t>(i)],
                                    info.fully_covers[static_cast<size_t>(j)]);
            bool do_swap = false;
            if (cc < 0) {
                do_swap = true;
            } else if (cc == 0) {
                if (use_error_rate && err[static_cast<size_t>(i)] > err[static_cast<size_t>(j)]) {
                    do_swap = true;
                } else if ((!use_error_rate ||
                            err[static_cast<size_t>(i)] == err[static_cast<size_t>(j)]) &&
                           info.lens[static_cast<size_t>(i)] < info.lens[static_cast<size_t>(j)]) {
                    do_swap = true;
                }
            }
            if (do_swap) {
                std::swap(info.noisy_read_ids[static_cast<size_t>(i)], info.noisy_read_ids[static_cast<size_t>(j)]);
                std::swap(info.lens          [static_cast<size_t>(i)], info.lens          [static_cast<size_t>(j)]);
                std::swap(info.seqs          [static_cast<size_t>(i)], info.seqs          [static_cast<size_t>(j)]);
                std::swap(info.quals         [static_cast<size_t>(i)], info.quals         [static_cast<size_t>(j)]);
                std::swap(info.strands       [static_cast<size_t>(i)], info.strands       [static_cast<size_t>(j)]);
                std::swap(info.fully_covers  [static_cast<size_t>(i)], info.fully_covers  [static_cast<size_t>(j)]);
                std::swap(info.haps          [static_cast<size_t>(i)], info.haps          [static_cast<size_t>(j)]);
                std::swap(info.phase_sets    [static_cast<size_t>(i)], info.phase_sets    [static_cast<size_t>(j)]);
                std::swap(err                [static_cast<size_t>(i)], err                [static_cast<size_t>(j)]);
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// collect_phase_set_with_both_haps  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

hts_pos_t collect_phase_set_with_both_haps(const NoisyReadInfo& info,
                                           int min_hap_full_reads,
                                           int min_hap_all_reads) {
    const int n = info.n_reads;
    std::vector<hts_pos_t> uniq(static_cast<size_t>(n), 0);
    int n_uniq = 0;

    struct PsStats { int full[2]{}; int all[2]{}; int min_full_len[2]{INT_MAX, INT_MAX}; };
    std::vector<PsStats> stats(static_cast<size_t>(n));

    for (int i = 0; i < n; ++i) {
        if (info.haps[static_cast<size_t>(i)] == 0) continue;
        const int ps_i = add_phase_set(info.phase_sets[static_cast<size_t>(i)], uniq, n_uniq);
        const int hap  = info.haps[static_cast<size_t>(i)] - 1; // 0-indexed
        const int len  = info.lens[static_cast<size_t>(i)];
        const int cov  = info.fully_covers[static_cast<size_t>(i)];

        if (noisyIsBothCover(cov)) {
            stats[static_cast<size_t>(ps_i)].full[hap]++;
            stats[static_cast<size_t>(ps_i)].all[hap]++;
            if (len < stats[static_cast<size_t>(ps_i)].min_full_len[hap])
                stats[static_cast<size_t>(ps_i)].min_full_len[hap] = len;
        } else if (noisyIsLeftCover(cov) || noisyIsRightCover(cov)) {
            if (len >= stats[static_cast<size_t>(ps_i)].min_full_len[hap])
                stats[static_cast<size_t>(ps_i)].all[hap]++;
        }
    }

    hts_pos_t max_ps = -1;
    int max_ps_i = -1;
    int max_minor = -1, max_major = -1;

    for (int i = 0; i < n_uniq; ++i) {
        const int minor = std::min(stats[static_cast<size_t>(i)].full[0], stats[static_cast<size_t>(i)].full[1]);
        const int major = std::max(stats[static_cast<size_t>(i)].full[0], stats[static_cast<size_t>(i)].full[1]);
        if (minor > max_minor ||
            (minor == max_minor && major > max_major)) {
            max_minor = minor; max_major = major;
            max_ps    = uniq[static_cast<size_t>(i)];
            max_ps_i  = i;
        }
    }

    if (max_minor < min_hap_full_reads) return -1;
    if (max_ps_i >= 0) {
        if (stats[static_cast<size_t>(max_ps_i)].all[0] < min_hap_all_reads ||
            stats[static_cast<size_t>(max_ps_i)].all[1] < min_hap_all_reads) return -1;
    }
    return max_ps;
}

// ════════════════════════════════════════════════════════════════════════════
// wfa_end2end_aln  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int wfa_end2end_aln(const uint8_t* pattern, int plen,
                    const uint8_t* text,    int tlen,
                    int gap_aln, int b, int q, int e, int q2, int e2,
                    int heuristic, int affine_gap,
                    std::vector<uint8_t>& pattern_alg, std::vector<uint8_t>& text_alg) {
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    if (affine_gap == kWfaAffine2p) {
        attr.distance_metric = gap_affine_2p;
        attr.affine2p_penalties.match       = 0;
        attr.affine2p_penalties.mismatch    = b;
        attr.affine2p_penalties.gap_opening1 = q;
        attr.affine2p_penalties.gap_extension1 = e;
        attr.affine2p_penalties.gap_opening2 = q2;
        attr.affine2p_penalties.gap_extension2 = e2;
    } else {
        attr.distance_metric = gap_affine;
        attr.affine_penalties.match     = 0;
        attr.affine_penalties.mismatch  = b;
        attr.affine_penalties.gap_opening  = q;
        attr.affine_penalties.gap_extension = e;
    }
    attr.alignment_scope = compute_alignment;
    attr.alignment_form.span = alignment_end2end;

    if (heuristic == kWfaNoHeuristic) {
        attr.heuristic.strategy = wf_heuristic_none;
    } else if (heuristic == kWfaAdaptive) {
        attr.heuristic.strategy = wf_heuristic_wfadaptive;
    } else if (heuristic == kWfaZdrop) {
        attr.heuristic.strategy = wf_heuristic_zdrop;
        attr.heuristic.zdrop = std::min(500, static_cast<int>(std::min(plen, tlen) * 0.1));
        attr.heuristic.steps_between_cutoffs = 100;
    }

    // Optionally reverse for left-gap alignment.
    std::vector<uint8_t> p_rev, t_rev;
    const uint8_t* p = pattern;
    const uint8_t* t = text;
    if (gap_aln == kGapLeftAln) {
        p_rev.resize(static_cast<size_t>(plen));
        t_rev.resize(static_cast<size_t>(tlen));
        for (int i = 0; i < plen; ++i) p_rev[static_cast<size_t>(i)] = pattern[plen - i - 1];
        for (int i = 0; i < tlen; ++i) t_rev[static_cast<size_t>(i)] = text   [tlen - i - 1];
        p = p_rev.data(); t = t_rev.data();
    }

    wavefront_aligner_t* wf = wavefront_aligner_new(&attr);
    wavefront_align(wf, reinterpret_cast<const char*>(p), plen,
                        reinterpret_cast<const char*>(t), tlen);

    int aln_len = wfa_collect_pretty_alignment(wf->cigar, p, plen, t, tlen,
                                               pattern_alg, text_alg);

    if (gap_aln == kGapLeftAln) {
        std::reverse(pattern_alg.begin(), pattern_alg.end());
        std::reverse(text_alg   .begin(), text_alg   .end());
    }

    wavefront_aligner_delete(wf);
    return aln_len;
}

// ════════════════════════════════════════════════════════════════════════════
// wfa_trim_aln_str  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

void wfa_trim_aln_str(int full_cover, AlnStr& s) {
    if (noisyIsNotCover(full_cover) || noisyIsBothCover(full_cover)) return;

    // One side cover + one side gap: keep the whole alignment as-is.
    if ((noisyIsLeftCover(full_cover)  && noisyIsRightGap(full_cover)) ||
        (noisyIsRightCover(full_cover) && noisyIsLeftGap(full_cover))) {
        s.target_beg = 0; s.target_end = s.aln_len - 1;
        s.query_beg  = 0; s.query_end  = s.aln_len - 1;
        return;
    }

    if (noisyIsLeftCover(full_cover)) {
        // Trim right: keep up to the last matching position in query.
        int target_end = -1, query_end = -1;
        for (int i = s.aln_len - 1; i >= 0; --i) {
            if (query_end == -1 && s.query_aln[static_cast<size_t>(i)] != 5 &&
                s.target_aln[static_cast<size_t>(i)] == s.query_aln[static_cast<size_t>(i)])
                query_end = i;
            if (target_end == -1 && s.target_aln[static_cast<size_t>(i)] != 5)
                target_end = i;
            if (target_end != -1 && query_end != -1) break;
        }
        if (query_end == -1) query_end = target_end;
        s.aln_len   = target_end + 1;
        s.target_beg = 0; s.target_end = target_end;
        s.query_beg  = 0; s.query_end  = query_end;
        for (int i = query_end + 1; i < s.aln_len; ++i)
            s.query_aln[static_cast<size_t>(i)] = 5;
    } else { // right-cover: trim left
        int query_start = -1, target_start = -1;
        for (int i = 0; i < s.aln_len; ++i) {
            if (query_start == -1 && s.query_aln[static_cast<size_t>(i)] != 5 &&
                s.target_aln[static_cast<size_t>(i)] == s.query_aln[static_cast<size_t>(i)])
                query_start = i;
            if (target_start == -1 && s.target_aln[static_cast<size_t>(i)] != 5)
                target_start = i;
            if (target_start != -1 && query_start != -1) break;
        }
        if (query_start == -1) query_start = target_start;
        s.aln_len = s.aln_len - target_start;
        if (target_start != 0) {
            // Shift both arrays left by target_start.
            std::vector<uint8_t> new_target(s.target_aln.begin() + target_start,
                                            s.target_aln.begin() + target_start + s.aln_len);
            std::vector<uint8_t> new_query (s.query_aln .begin() + target_start,
                                            s.query_aln .begin() + target_start + s.aln_len);
            s.target_aln = std::move(new_target);
            s.query_aln  = std::move(new_query);
        }
        s.target_beg = 0; s.target_end = s.aln_len - 1;
        s.query_beg  = query_start - target_start;
        s.query_end  = s.aln_len - 1;
        for (int i = 0; i < s.query_beg; ++i)
            s.query_aln[static_cast<size_t>(i)] = 5;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// wfa_collect_aln_str  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int wfa_collect_aln_str(const Options& opts, const uint8_t* target, int tlen,
                        const uint8_t* query,  int qlen,
                        int full_cover, int heuristic, int affine_gap, AlnStr& s) {
    if (noisyIsNotCover(full_cover)) return 0;
    s = AlnStr{};

    const int gap_aln = opts.gap_aln;
    const int b = opts.mismatch, q = opts.gap_open1, e = opts.gap_ext1;
    const int q2 = opts.gap_open2, e2 = opts.gap_ext2;

    if (noisyIsBothCover(full_cover)) {
        int aln_len = wfa_end2end_aln(target, tlen, query, qlen,
                                      gap_aln, b, q, e, q2, e2,
                                      heuristic, affine_gap,
                                      s.target_aln, s.query_aln);
        s.aln_len = aln_len;
        s.target_beg = 0; s.target_end = aln_len - 1;
        s.query_beg  = 0; s.query_end  = aln_len - 1;
    } else {
        const int ext_dir = noisyIsLeftCover(full_cover)
            ? kExtAlnLeftToRight : kExtAlnRightToLeft;
        const double ratio = opts.partial_aln_ratio;
        int _tlen = tlen, _qlen = qlen, t_start = 0, q_start = 0;

        if (ext_dir == kExtAlnLeftToRight) {
            if (tlen > static_cast<int>(qlen * ratio)) _tlen = static_cast<int>(qlen * ratio);
            else if (qlen > static_cast<int>(tlen * ratio)) _qlen = static_cast<int>(tlen * ratio);
        } else {
            if (tlen > static_cast<int>(qlen * ratio)) {
                _tlen = static_cast<int>(qlen * ratio); t_start = tlen - _tlen;
            } else if (qlen > static_cast<int>(tlen * ratio)) {
                _qlen = static_cast<int>(tlen * ratio); q_start = qlen - _qlen;
            }
        }
        // Flip gap direction for left-to-right partial alignment.
        int eff_gap_aln = gap_aln;
        if (ext_dir == kExtAlnLeftToRight)
            eff_gap_aln = (gap_aln == kGapLeftAln) ? kGapRightAln : kGapLeftAln;

        wfa_end2end_aln(target + t_start, _tlen, query + q_start, _qlen,
                        eff_gap_aln, b, q, e, q2, e2,
                        kWfaZdrop, kWfaAffine2p,
                        s.target_aln, s.query_aln);
        s.aln_len = static_cast<int>(s.target_aln.size());
        s.target_beg = 0; s.target_end = s.aln_len - 1;
        s.query_beg  = 0; s.query_end  = s.aln_len - 1;
        wfa_trim_aln_str(full_cover, s);
        s.target_beg += t_start; s.target_end += t_start;
        s.query_beg  += q_start; s.query_end  += q_start;
    }
    return 0;
}

// ════════════════════════════════════════════════════════════════════════════
// make_cons_read_aln_str  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int make_cons_read_aln_str(const uint8_t* cons_row, const uint8_t* read_row,
                           int msa_len, int full_cover, AlnStr& out) {
    out = AlnStr{};
    out.target_aln.reserve(static_cast<size_t>(msa_len));
    out.query_aln .reserve(static_cast<size_t>(msa_len));
    for (int i = 0; i < msa_len; ++i) {
        if (cons_row[i] != 5 || read_row[i] != 5) { // skip double-gap columns
            out.target_aln.push_back(cons_row[i]);
            out.query_aln .push_back(read_row[i]);
        }
    }
    out.aln_len    = static_cast<int>(out.target_aln.size());
    out.target_beg = 0; out.target_end = out.aln_len - 1;
    out.query_beg  = 0; out.query_end  = out.aln_len - 1;
    wfa_trim_aln_str(full_cover, out);
    return out.aln_len;
}

// ════════════════════════════════════════════════════════════════════════════
// make_ref_read_aln_str  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int make_ref_read_aln_str(const Options& opts, const AlnStr& ref_cons,
                          const AlnStr& cons_read, AlnStr& ref_read) {
    const int max_len = ref_cons.aln_len + cons_read.aln_len;
    ref_read = AlnStr{};
    ref_read.target_aln.reserve(static_cast<size_t>(max_len));
    ref_read.query_aln .reserve(static_cast<size_t>(max_len));

    int i = 0, j = 0;
    while (i < ref_cons.aln_len && j < cons_read.aln_len) {
        const uint8_t rc_query = ref_cons.query_aln [static_cast<size_t>(i)];
        const uint8_t cr_tgt   = cons_read.target_aln[static_cast<size_t>(j)];

        if (rc_query == 5 && cr_tgt == 5) {
            // Both gaps in the cons row: need to locally re-align ref-del vs read-ins.
            int ref_del_len = 1;
            while (i + ref_del_len < ref_cons.aln_len &&
                   ref_cons.query_aln[static_cast<size_t>(i + ref_del_len)] == 5)
                ++ref_del_len;
            int read_del_len = 1;
            while (j + read_del_len < cons_read.aln_len &&
                   cons_read.target_aln[static_cast<size_t>(j + read_del_len)] == 5)
                ++read_del_len;

            std::vector<uint8_t> ref_aln, read_aln;
            wfa_end2end_aln(ref_cons.target_aln .data() + i, ref_del_len,
                            cons_read.query_aln.data() + j, read_del_len,
                            opts.gap_aln, opts.mismatch,
                            opts.gap_open1, opts.gap_ext1,
                            opts.gap_open2, opts.gap_ext2,
                            kWfaNoHeuristic, kWfaAffine2p,
                            ref_aln, read_aln);
            for (size_t k = 0; k < ref_aln.size(); ++k) {
                ref_read.target_aln.push_back(ref_aln[k]);
                ref_read.query_aln .push_back(read_aln[k]);
            }
            i += ref_del_len; j += read_del_len;
        } else if (rc_query != 5 && cr_tgt != 5) {
            ref_read.target_aln.push_back(ref_cons .target_aln[static_cast<size_t>(i)]);
            ref_read.query_aln .push_back(cons_read.query_aln [static_cast<size_t>(j)]);
            ++i; ++j;
        } else if (rc_query == 5) {
            ref_read.target_aln.push_back(ref_cons.target_aln[static_cast<size_t>(i)]);
            ref_read.query_aln .push_back(5);
            ++i;
        } else {
            ref_read.target_aln.push_back(5);
            ref_read.query_aln .push_back(cons_read.query_aln[static_cast<size_t>(j)]);
            ++j;
        }
    }
    while (i < ref_cons.aln_len) {
        ref_read.target_aln.push_back(ref_cons.target_aln[static_cast<size_t>(i)]);
        ref_read.query_aln .push_back(5);
        ++i;
    }
    while (j < cons_read.aln_len) {
        ref_read.target_aln.push_back(5);
        ref_read.query_aln .push_back(cons_read.query_aln[static_cast<size_t>(j)]);
        ++j;
    }
    ref_read.aln_len = static_cast<int>(ref_read.target_aln.size());
    ref_read.target_beg = ref_read.target_end = ref_read.query_beg = ref_read.query_end = -1;
    return ref_read.aln_len;
}

// ════════════════════════════════════════════════════════════════════════════
// abpoa_aln_msa_cons  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

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
                       std::array<std::vector<std::vector<uint8_t>>, 2>& msa_seqs) {
    abpoa_t*      ab   = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->wb = -1;
    abpt->inc_path_score = 1;
    abpt->out_cons = 1;
    abpt->out_msa  = 1;
    abpt->cons_algrm = ABPOA_MF;
    abpt->max_n_cons = max_n_cons;
    abpt->min_freq   = static_cast<float>(opts.min_af);
    abpt->match      = opts.match;
    abpt->mismatch   = opts.mismatch;
    abpt->gap_open1  = opts.gap_open1;
    abpt->gap_ext1   = opts.gap_ext1;
    abpt->gap_open2  = opts.gap_open2;
    abpt->gap_ext2   = opts.gap_ext2;
    abpoa_post_set_para(abpt);

    // Build raw C arrays for abpoa_msa.
    std::vector<int> lens_c(read_lens.begin(), read_lens.end());
    abpoa_msa(ab, abpt, n_reads, nullptr, lens_c.data(),
              const_cast<uint8_t**>(read_seqs.data()), nullptr, nullptr);

    abpoa_cons_t* abc = ab->abc;
    int n_cons = 0;
    if (abc->n_cons > 0) {
        cons_lens.resize(static_cast<size_t>(abc->n_cons));
        cons_seqs.resize(static_cast<size_t>(abc->n_cons));
        for (int ci = 0; ci < abc->n_cons; ++ci) {
            cons_lens[static_cast<size_t>(ci)] = abc->cons_len[ci];
            cons_seqs[static_cast<size_t>(ci)].assign(
                abc->cons_base[ci], abc->cons_base[ci] + abc->cons_len[ci]);
        }

        if (abc->n_cons == 2) {
            for (int ci = 0; ci < 2; ++ci) {
                clu_n_seqs[static_cast<size_t>(ci)] = abc->clu_n_seq[ci];
                clu_read_ids[static_cast<size_t>(ci)].resize(static_cast<size_t>(abc->clu_n_seq[ci]));
                for (int k = 0; k < abc->clu_n_seq[ci]; ++k)
                    clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(k)] =
                        read_ids[static_cast<size_t>(abc->clu_read_ids[ci][k])];
            }
        } else {
            clu_n_seqs[0] = n_reads; clu_n_seqs[1] = 0;
            clu_read_ids[0].resize(static_cast<size_t>(n_reads));
            for (int k = 0; k < n_reads; ++k)
                clu_read_ids[0][static_cast<size_t>(k)] = read_ids[static_cast<size_t>(k)];
        }

        // Collect MSA rows split by cluster.
        msa_seq_lens.resize(static_cast<size_t>(abc->n_cons));
        for (int ci = 0; ci < abc->n_cons; ++ci) {
            msa_seq_lens[static_cast<size_t>(ci)] = abc->msa_len;
            const int n_read_rows = (abc->n_cons > 1) ? abc->clu_n_seq[ci] : n_reads;
            const int n_rows = n_read_rows + 1; // reads + consensus
            msa_seqs[static_cast<size_t>(ci)].resize(static_cast<size_t>(n_rows));
            if (abc->n_cons > 1) {
                for (int k = 0; k < abc->clu_n_seq[ci]; ++k) {
                    const int src = abc->clu_read_ids[ci][k];
                    msa_seqs[static_cast<size_t>(ci)][static_cast<size_t>(k)].assign(
                        abc->msa_base[src], abc->msa_base[src] + abc->msa_len);
                }
            } else {
                for (int k = 0; k < n_reads; ++k) {
                    msa_seqs[static_cast<size_t>(ci)][static_cast<size_t>(k)].assign(
                        abc->msa_base[k], abc->msa_base[k] + abc->msa_len);
                }
            }
            // consensus row is stored after all reads in abc->msa_base
            int cons_row = abc->n_seq + ci;
            msa_seqs[static_cast<size_t>(ci)][static_cast<size_t>(n_read_rows)].assign(
                abc->msa_base[cons_row], abc->msa_base[cons_row] + abc->msa_len);
        }
        n_cons = abc->n_cons;
    }
    abpoa_free_para(abpt);
    abpoa_free(ab);
    return n_cons;
}

// ════════════════════════════════════════════════════════════════════════════
// abpoa_partial_aln_msa_cons  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

// ── longcallD exact X-gap counting (edlib_xgaps) ─────────────────────────────

static int edlibAlignmentToXGAPS(const unsigned char* const alignment, const int alignmentLength) {
    // Mirrors longcallD edlibAlignmentToXGAPS:
    // mismatch counts as 1; insert/delete counts as 1 per run (gap-open events).
    int n_gaps = 0, n_mismatch = 0;
    for (int i = 0; i < alignmentLength; i++) {
        if (alignment[i] == EDLIB_EDOP_MATCH) {
            continue;
        } else if (alignment[i] == EDLIB_EDOP_MISMATCH) {
            n_mismatch++;
        } else if (alignment[i] == EDLIB_EDOP_INSERT || alignment[i] == EDLIB_EDOP_DELETE) {
            if (i == 0 || alignment[i - 1] != alignment[i]) n_gaps++;
        }
    }
    return n_mismatch + n_gaps;
}

static int edlib_xgaps(const uint8_t* target, int tlen, const uint8_t* query, int qlen) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    int n_x_gaps = edlibAlignmentToXGAPS(result.alignment, result.alignmentLength);
    edlibFreeAlignResult(result);
    return n_x_gaps;
}

static void collect_aln_beg_end_from_pretty(const int ext_direction,
                                             const int ref_len, const int read_len,
                                             const std::vector<uint8_t>& ref_aln,
                                             const std::vector<uint8_t>& read_aln,
                                             int* ref_beg, int* ref_end,
                                             int* query_beg, int* query_end) {
    // Mirrors longcallD collect_aln_beg_end, but reconstructs CIGAR ops from
    // pretty alignment columns where gaps are marked as 5.
    *ref_beg = 1; *query_beg = 1; *ref_end = ref_len, *query_end = read_len;
    const int n = static_cast<int>(ref_aln.size());

    if (ext_direction == kExtAlnLeftToRight) {
        int tmp_ref_end = 0, tmp_read_end = 0;
        for (int i = 0; i < n; ++i) {
            const uint8_t r = ref_aln[static_cast<size_t>(i)];
            const uint8_t q = read_aln[static_cast<size_t>(i)];
            const bool r_gap = (r == 5);
            const bool q_gap = (q == 5);
            if (!r_gap && !q_gap) {
                tmp_ref_end++; tmp_read_end++;
                // CEQUAL/CMATCH in longcallD corresponds to bases equal.
                if (r == q) {
                    *ref_end = tmp_ref_end;
                    *query_end = tmp_read_end;
                }
            } else if (!r_gap && q_gap) {
                // CDEL: reference consumed, query not.
                tmp_ref_end++;
            } else if (r_gap && !q_gap) {
                // CINS: query consumed, reference not.
                tmp_read_end++;
            }
        }
    } else {
        int tmp_ref_beg = ref_len + 1, tmp_read_beg = read_len + 1;
        for (int i = n - 1; i >= 0; --i) {
            const uint8_t r = ref_aln[static_cast<size_t>(i)];
            const uint8_t q = read_aln[static_cast<size_t>(i)];
            const bool r_gap = (r == 5);
            const bool q_gap = (q == 5);
            if (!r_gap && !q_gap) {
                tmp_ref_beg--; tmp_read_beg--;
                if (r == q) {
                    *ref_beg = tmp_ref_beg;
                    *query_beg = tmp_read_beg;
                }
            } else if (!r_gap && q_gap) {
                tmp_ref_beg--;
            } else if (r_gap && !q_gap) {
                tmp_read_beg--;
            }
        }
    }
}

static int cal_wfa_partial_aln_beg_end(const int ext_direction, const Options& opts,
                                       const uint8_t* _target, int _tlen,
                                       const uint8_t* _query, int _qlen,
                                       int* target_beg, int* target_end,
                                       int* query_beg, int* query_end) {
    int gap_aln = opts.gap_aln;
    int b = opts.mismatch, q = opts.gap_open1, e = opts.gap_ext1;
    int q2 = opts.gap_open2, e2 = opts.gap_ext2;
    double ratio = opts.partial_aln_ratio;
    int tlen = _tlen, qlen = _qlen;
    const uint8_t* target = _target;
    const uint8_t* query = _query;

    if (ext_direction == kExtAlnLeftToRight) {
        if (_tlen > _qlen * ratio) tlen = static_cast<int>(_qlen * ratio);
        else if (_qlen > _tlen * ratio) qlen = static_cast<int>(_tlen * ratio);
    } else if (ext_direction == kExtAlnRightToLeft) {
        if (_tlen > _qlen * ratio) {
            target = _target + _tlen - static_cast<int>(_qlen * ratio);
            tlen = static_cast<int>(_qlen * ratio);
        } else if (_qlen > _tlen * ratio) {
            query = _query + _qlen - static_cast<int>(_tlen * ratio);
            qlen = static_cast<int>(_tlen * ratio);
        }
    }

    if (ext_direction == kExtAlnLeftToRight) {
        gap_aln = (gap_aln == kGapRightAln) ? kGapLeftAln : kGapRightAln;
    }

    int min_len = std::min(tlen, qlen);
    if (ext_direction == kExtAlnLeftToRight) {
        int x_gaps = edlib_xgaps(target, min_len, query, min_len);
        if (x_gaps > min_len * 0.10) return 0;
    } else {
        int x_gaps = edlib_xgaps(target + tlen - min_len, min_len, query + qlen - min_len, min_len);
        if (x_gaps > min_len * 0.10) return 0;
    }

    std::vector<uint8_t> pat_alg, txt_alg;
    const int cigar_len = wfa_end2end_aln(target, tlen, query, qlen, gap_aln, b, q, e, q2, e2,
                                          kWfaNoHeuristic, kWfaAffine2p, pat_alg, txt_alg);
    int ret = 1;
    if (cigar_len == 0) ret = 0;
    else collect_aln_beg_end_from_pretty(ext_direction, _tlen, _qlen, pat_alg, txt_alg,
                                          target_beg, target_end, query_beg, query_end);
    return ret;
}

static int collect_partial_aln_beg_end(const Options& opts, const int sampling_reads,
                                        const uint8_t* target, int tlen, int target_full_cover,
                                        const uint8_t* query, int qlen, int query_full_cover,
                                        int* target_beg, int* target_end,
                                        int* query_beg, int* query_end) {
    *target_beg = 1; *target_end = tlen; *query_beg = 1; *query_end = qlen;
    int ret = 1;
    assert(noisyIsBothCover(target_full_cover) != 0);

    if (noisyIsBothCover(target_full_cover)) {
        if (noisyIsBothCover(query_full_cover) ||
            (noisyIsLeftCover(query_full_cover) && noisyIsRightGap(query_full_cover)) ||
            (noisyIsRightCover(query_full_cover) && noisyIsLeftGap(query_full_cover))) {
            if (sampling_reads) {
                int x_gaps = edlib_xgaps(target, tlen, query, qlen);
                if (x_gaps > std::min(tlen, qlen) * 0.10) return 0;
            }
            return 1;
        } else {
            if (noisyIsLeftCover(query_full_cover)) {
                ret = cal_wfa_partial_aln_beg_end(kExtAlnLeftToRight, opts, target, tlen, query, qlen,
                                                   target_beg, target_end, query_beg, query_end);
            } else if (noisyIsRightCover(query_full_cover)) {
                ret = cal_wfa_partial_aln_beg_end(kExtAlnRightToLeft, opts, target, tlen, query, qlen,
                                                   target_beg, target_end, query_beg, query_end);
            }
        }
    } else if (noisyIsLeftCover(target_full_cover)) {
        assert(noisyIsLeftCover(query_full_cover) != 0 && "target is left-cover but query is not left-cover");
        ret = cal_wfa_partial_aln_beg_end(kExtAlnLeftToRight, opts, target, tlen, query, qlen,
                                          target_beg, target_end, query_beg, query_end);
    } else if (noisyIsRightCover(target_full_cover)) {
        assert(noisyIsRightCover(query_full_cover) != 0 && "target is right-cover but query is not right-cover");
        ret = cal_wfa_partial_aln_beg_end(kExtAlnRightToLeft, opts, target, tlen, query, qlen,
                                          target_beg, target_end, query_beg, query_end);
    } else {
        return 0;
    }
    return ret;
}

int abpoa_partial_aln_msa_cons(const Options& opts, int sampling_reads,
                               int n_reads,
                               const std::vector<int>& read_ids,
                               const std::vector<uint8_t*>& read_seqs,
                               const std::vector<uint8_t*>& /*read_quals*/,
                               const std::vector<int>& read_lens,
                               const std::vector<int>& read_full_covers,
                               int /*max_n_cons*/,
                               int& cons_len_out,
                               std::vector<uint8_t>& cons_seq_out,
                               int& clu_n_seqs_out,
                               std::vector<int>& clu_read_ids_out,
                               int& msa_seq_len_out,
                               std::vector<std::vector<uint8_t>>& msa_seqs_out) {
    abpoa_t*      ab   = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->sub_aln = 1;
    abpt->inc_path_score = 1;
    abpt->out_cons = 1;
    abpt->out_msa  = 1;
    abpt->cons_algrm = ABPOA_MF;
    abpt->max_n_cons = 1;
    abpt->min_freq   = static_cast<float>(opts.min_af);
    abpt->match      = opts.match;
    abpt->mismatch   = opts.mismatch;
    abpt->gap_open1  = opts.gap_open1;
    abpt->gap_ext1   = opts.gap_ext1;
    abpt->gap_open2  = opts.gap_open2;
    abpt->gap_ext2   = opts.gap_ext2;
    abpoa_post_set_para(abpt);

    ab->abs->n_seq = n_reads;
    int n_added = 0;
    for (int i = 0; i < n_reads; ++i) {
        if (read_lens[static_cast<size_t>(i)] <= 0) continue;

        abpoa_res_t res{}; res.graph_cigar = nullptr; res.n_cigar = 0;
        int exc_beg = 0, exc_end = 1;
        int seq_beg_cut = 0, seq_end_cut = 0;

        if (i != 0) {
            int ref_beg = 1, ref_end = read_lens[0];
            int read_beg = 1, read_end = read_lens[static_cast<size_t>(i)];
            if (collect_partial_aln_beg_end(
                    opts, sampling_reads,
                    read_seqs[0], read_lens[0], read_full_covers[0],
                    read_seqs[static_cast<size_t>(i)], read_lens[static_cast<size_t>(i)], read_full_covers[static_cast<size_t>(i)],
                    &ref_beg, &ref_end, &read_beg, &read_end) == 0) {
                continue;
            }
            const int beg_id = ref_beg + 1;
            const int end_id = ref_end + 1;
            seq_beg_cut = read_beg - 1;
            seq_end_cut = read_lens[static_cast<size_t>(i)] - read_end;
            abpoa_subgraph_nodes(ab, abpt, beg_id, end_id, &exc_beg, &exc_end);
        }
        uint8_t* seq_ptr = read_seqs[static_cast<size_t>(i)] + seq_beg_cut;
        const int seq_use_len  = read_lens[static_cast<size_t>(i)] - seq_beg_cut - seq_end_cut;
        if (seq_use_len <= 0) continue;

        abpoa_align_sequence_to_subgraph(ab, abpt, exc_beg, exc_end,
                                         seq_ptr, seq_use_len, &res);
        abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end,
                                     seq_ptr, nullptr, seq_use_len, nullptr,
                                     res, i, n_reads, 0);
        if (res.n_cigar) free(res.graph_cigar);
        ++n_added;
    }

    abpoa_output(ab, abpt, nullptr);
    abpoa_cons_t* abc = ab->abc;
    int n_cons = 0;

    if (abc->n_cons > 0) {
        cons_len_out = abc->cons_len[0];
        cons_seq_out.assign(abc->cons_base[0], abc->cons_base[0] + abc->cons_len[0]);
        if (abc->n_cons > 1) {
            clu_n_seqs_out = abc->clu_n_seq[0];
            clu_read_ids_out.resize(static_cast<size_t>(abc->clu_n_seq[0]));
            for (int k = 0; k < abc->clu_n_seq[0]; ++k) {
                clu_read_ids_out[static_cast<size_t>(k)] =
                    read_ids[static_cast<size_t>(abc->clu_read_ids[0][k])];
            }
        } else {
            clu_n_seqs_out = n_reads;
            clu_read_ids_out.resize(static_cast<size_t>(n_reads));
            for (int k = 0; k < n_reads; ++k) {
                clu_read_ids_out[static_cast<size_t>(k)] = read_ids[static_cast<size_t>(k)];
            }
        }

        msa_seq_len_out = abc->msa_len;
        const int n_read_rows = (abc->n_cons > 1) ? abc->clu_n_seq[0] : n_reads;
        const int n_rows = n_read_rows + 1;
        msa_seqs_out.resize(static_cast<size_t>(n_rows));
        if (abc->n_cons > 1) {
            for (int k = 0; k < abc->clu_n_seq[0]; ++k) {
                const int src = abc->clu_read_ids[0][k];
                msa_seqs_out[static_cast<size_t>(k)].assign(
                    abc->msa_base[src], abc->msa_base[src] + abc->msa_len);
            }
        } else {
            for (int k = 0; k < n_reads; ++k) {
                msa_seqs_out[static_cast<size_t>(k)].assign(
                    abc->msa_base[k], abc->msa_base[k] + abc->msa_len);
            }
        }
        int cons_row = abc->n_seq;
        msa_seqs_out[static_cast<size_t>(n_read_rows)].assign(
            abc->msa_base[cons_row], abc->msa_base[cons_row] + abc->msa_len);
        n_cons = 1;
    }

    abpoa_free_para(abpt);
    abpoa_free(ab);
    return n_cons;
}

// ════════════════════════════════════════════════════════════════════════════
// wfa_collect_noisy_aln_str_no_ps_hap  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int wfa_collect_noisy_aln_str_no_ps_hap(const Options& opts, NoisyReadInfo& info,
                                         const uint8_t* ref_seq, int ref_seq_len,
                                         bool collect_ref_read_aln_str,
                                         std::array<int, 2>& clu_n_seqs,
                                         std::array<std::vector<int>, 2>& clu_read_ids,
                                         std::array<std::vector<AlnStr>, 2>& aln_strs) {
    // Collect full-cover reads.
    std::vector<int>       full_ids;
    std::vector<uint8_t*>  full_seqs;
    std::vector<int>       full_lens;
    std::vector<int>       full_covers;

    const int n = info.n_reads;
    for (int i = 0; i < n; ++i) {
        if (info.lens[static_cast<size_t>(i)] <= 0 ||
            !noisyIsBothCover(info.fully_covers[static_cast<size_t>(i)])) continue;
        full_ids   .push_back(i);
        full_seqs  .push_back(info.seqs[static_cast<size_t>(i)].data());
        full_lens  .push_back(info.lens[static_cast<size_t>(i)]);
        full_covers.push_back(info.fully_covers[static_cast<size_t>(i)]);
    }
    if (full_ids.empty()) return 0;
    if (full_lens[0] >= opts.max_noisy_reg_len) return 0;

    // Run de-novo abPOA.
    std::vector<int>                           cons_lens;
    std::vector<std::vector<uint8_t>>          cons_seqs;
    std::vector<int>                           msa_seq_lens;
    std::array<std::vector<std::vector<uint8_t>>, 2> msa_seqs;

    const int n_full = static_cast<int>(full_ids.size());
    int n_cons = abpoa_aln_msa_cons(opts, n_full, full_ids, full_seqs, full_lens,
                                    2, cons_lens, cons_seqs, clu_n_seqs, clu_read_ids,
                                    msa_seq_lens, msa_seqs);
    if (n_cons == 0) return 0;
    if (opts.verbose >= 2) std::fprintf(stderr, "n_cons: %d\n", n_cons);

    for (int ci = 0; ci < n_cons; ++ci) {
        // aln_strs[ci][0] = ref-vs-cons.
        aln_strs[static_cast<size_t>(ci)].resize(
            static_cast<size_t>(1 + clu_n_seqs[static_cast<size_t>(ci)] * 2));

        wfa_collect_aln_str(opts, ref_seq, ref_seq_len,
                            cons_seqs[static_cast<size_t>(ci)].data(),
                            cons_lens[static_cast<size_t>(ci)],
                            kNoisyBothCover, kWfaNoHeuristic, kWfaAffine2p,
                            aln_strs[static_cast<size_t>(ci)][0]);

        for (int k = 0; k < clu_n_seqs[static_cast<size_t>(ci)]; ++k) {
            const int read_i = clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(k)];
            const int read_id = info.noisy_read_ids[static_cast<size_t>(read_i)];
            clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(k)] = read_id;
            const int msa_len  = msa_seq_lens[static_cast<size_t>(ci)];
            const uint8_t* cons_row = msa_seqs[static_cast<size_t>(ci)]
                                               [static_cast<size_t>(clu_n_seqs[static_cast<size_t>(ci)])].data();
            const uint8_t* read_row = msa_seqs[static_cast<size_t>(ci)]
                                               [static_cast<size_t>(k)].data();
            // cons-vs-read aln: slot (k+1)*2 - 1.
            const size_t cons_read_slot = static_cast<size_t>((k + 1) * 2 - 1);
            make_cons_read_aln_str(cons_row, read_row, msa_len,
                                  info.fully_covers[static_cast<size_t>(read_i)],
                                   aln_strs[static_cast<size_t>(ci)][cons_read_slot]);
            const size_t ref_read_slot = static_cast<size_t>((k + 1) * 2);
            if (collect_ref_read_aln_str) {
                make_ref_read_aln_str(opts,
                                      aln_strs[static_cast<size_t>(ci)][0],
                                      aln_strs[static_cast<size_t>(ci)][cons_read_slot],
                                      aln_strs[static_cast<size_t>(ci)][ref_read_slot]);
            }
            if (opts.verbose >= 2) {
                std::fprintf(stderr, "%d:\n", read_id);
                std::fprintf(stderr, "cons:\t");
                for (int kk = 0; kk < aln_strs[static_cast<size_t>(ci)][cons_read_slot].aln_len; ++kk) {
                    std::fprintf(stderr, "%c",
                                 "ACGTN-"[aln_strs[static_cast<size_t>(ci)][cons_read_slot].target_aln[static_cast<size_t>(kk)]]);
                }
                std::fprintf(stderr, "\nread:\t");
                for (int kk = 0; kk < aln_strs[static_cast<size_t>(ci)][cons_read_slot].aln_len; ++kk) {
                    std::fprintf(stderr, "%c",
                                 "ACGTN-"[aln_strs[static_cast<size_t>(ci)][cons_read_slot].query_aln[static_cast<size_t>(kk)]]);
                }
                std::fprintf(stderr, "\n");
            }
        }
    }
    return n_cons;
}

// ════════════════════════════════════════════════════════════════════════════
// wfa_collect_noisy_aln_str_with_ps_hap  (mirrors longcallD align.c)
// ════════════════════════════════════════════════════════════════════════════

int wfa_collect_noisy_aln_str_with_ps_hap(const Options& opts, bool sampling_reads,
                                           NoisyReadInfo& info,
                                           hts_pos_t ps,
                                           int /*min_hap_full_reads*/, int /*min_hap_all_reads*/,
                                           const uint8_t* ref_seq, int ref_seq_len,
                                           bool collect_ref_read_aln_str,
                                           std::array<int, 2>& clu_n_seqs,
                                           std::array<std::vector<int>, 2>& clu_read_ids,
                                           std::array<std::vector<AlnStr>, 2>& aln_strs) {
    const int n = info.n_reads;
    // Check if region is a homopolymer — if so, require full-cover reads only.
    int hp_start, hp_end, hp_len;
    const bool use_non_full =
        !is_homopolymer(ref_seq, ref_seq_len, opts.noisy_reg_flank_len,
                        &hp_start, &hp_end, &hp_len);

    int n_cons = 0;
    std::array<std::vector<uint8_t>, 2>          cons_seqs{};
    std::array<int, 2>                           cons_lens{0, 0};
    std::array<std::vector<std::vector<uint8_t>>, 2> msa_seqs{};
    std::array<int, 2>                           msa_seq_lens{0, 0};

    for (int hap = 1; hap <= 2; ++hap) {
        const int ci = hap - 1;
        // Collect reads for this haplotype.
        std::vector<int>      ps_ids;
        std::vector<uint8_t*> ps_seqs;
        std::vector<uint8_t*> ps_quals;
        std::vector<int>      ps_lens;
        std::vector<int>      ps_covers;

        for (int i = 0; i < n; ++i) {
            if (info.lens[static_cast<size_t>(i)] <= 0) continue;
            if (info.phase_sets[static_cast<size_t>(i)] != ps) continue;
            if (info.haps[static_cast<size_t>(i)] != hap) continue;
            if (!use_non_full && !noisyIsBothCover(info.fully_covers[static_cast<size_t>(i)])) continue;

            ps_ids   .push_back(info.noisy_read_ids[static_cast<size_t>(i)]);
            ps_seqs  .push_back(info.seqs [static_cast<size_t>(i)].data());
            ps_quals .push_back(info.quals[static_cast<size_t>(i)].data());
            ps_lens  .push_back(info.lens [static_cast<size_t>(i)]);
            ps_covers.push_back(info.fully_covers[static_cast<size_t>(i)]);
        }
        if (ps_ids.empty()) continue;
        if (ps_lens[0] >= opts.max_noisy_reg_len) break;

        int out_clu_n  = 0;
        std::vector<int> out_clu_ids;
        int out_msa_len = 0;
        std::vector<std::vector<uint8_t>> out_msa;

        const int n_ps = static_cast<int>(ps_ids.size());
        // longcallD passes global read ids into abPOA, not 0..n_ps-1.
        // This matters because abPOA cluster read-ids are propagated downstream.
        int ret = abpoa_partial_aln_msa_cons(opts, sampling_reads, n_ps,
                                             ps_ids, ps_seqs, ps_quals, ps_lens, ps_covers,
                                             1,
                                             cons_lens[static_cast<size_t>(ci)], cons_seqs[static_cast<size_t>(ci)],
                                             out_clu_n, out_clu_ids,
                                             out_msa_len, out_msa);
        if (ret) {
            clu_n_seqs  [static_cast<size_t>(ci)] = out_clu_n;
            clu_read_ids[static_cast<size_t>(ci)] = out_clu_ids;
            msa_seq_lens[static_cast<size_t>(ci)] = out_msa_len;
            msa_seqs    [static_cast<size_t>(ci)] = std::move(out_msa);
            ++n_cons;
        }
    }

    if (n_cons != 2) return 0;

    // Build alignment strings for each haplotype.
    for (int hap = 1; hap <= 2; ++hap) {
        const int ci = hap - 1;
        const int n_hap = clu_n_seqs[static_cast<size_t>(ci)];
        aln_strs[static_cast<size_t>(ci)].resize(static_cast<size_t>(1 + n_hap * 2));

        // [0] = ref-vs-cons.
        wfa_collect_aln_str(opts, ref_seq, ref_seq_len,
                            cons_seqs[static_cast<size_t>(ci)].data(),
                            cons_lens[static_cast<size_t>(ci)],
                            kNoisyBothCover, kWfaNoHeuristic, kWfaAffine2p,
                            aln_strs[static_cast<size_t>(ci)][0]);

        int k_out = 0;
        for (int i = 0; i < n; ++i) {
            if (info.lens[static_cast<size_t>(i)] <= 0) continue;
            if (info.phase_sets[static_cast<size_t>(i)] != ps) continue;
            if (info.haps[static_cast<size_t>(i)] != hap) continue;
            if (!use_non_full && !noisyIsBothCover(info.fully_covers[static_cast<size_t>(i)])) continue;

            const int msa_len = msa_seq_lens[static_cast<size_t>(ci)];
            const uint8_t* cons_row = msa_seqs[static_cast<size_t>(ci)]
                                               [static_cast<size_t>(clu_n_seqs[static_cast<size_t>(ci)])].data();
            const uint8_t* read_row = msa_seqs[static_cast<size_t>(ci)]
                                               [static_cast<size_t>(k_out)].data();
            const size_t cons_read_slot = static_cast<size_t>((k_out + 1) * 2 - 1);
            make_cons_read_aln_str(cons_row, read_row, msa_len,
                                   info.fully_covers[static_cast<size_t>(i)],
                                   aln_strs[static_cast<size_t>(ci)][cons_read_slot]);
            const size_t ref_read_slot = static_cast<size_t>((k_out + 1) * 2);
            if (collect_ref_read_aln_str) {
                make_ref_read_aln_str(opts,
                                      aln_strs[static_cast<size_t>(ci)][0],
                                      aln_strs[static_cast<size_t>(ci)][cons_read_slot],
                                      aln_strs[static_cast<size_t>(ci)][ref_read_slot]);
            }
            if (opts.verbose >= 2) {
                std::fprintf(stderr, "%d\tps=%" PRId64 "\thap=%d\tfull_cover=%d\n",
                             info.noisy_read_ids[static_cast<size_t>(i)],
                             static_cast<int64_t>(ps),
                             hap,
                             info.fully_covers[static_cast<size_t>(i)]);
                std::fprintf(stderr, "cons:\t");
                for (int kk = 0; kk < aln_strs[static_cast<size_t>(ci)][cons_read_slot].aln_len; ++kk) {
                    std::fprintf(stderr, "%c",
                                 "ACGTN-"[aln_strs[static_cast<size_t>(ci)][cons_read_slot].target_aln[static_cast<size_t>(kk)]]);
                }
                std::fprintf(stderr, "\nread:\t");
                for (int kk = 0; kk < aln_strs[static_cast<size_t>(ci)][cons_read_slot].aln_len; ++kk) {
                    std::fprintf(stderr, "%c",
                                 "ACGTN-"[aln_strs[static_cast<size_t>(ci)][cons_read_slot].query_aln[static_cast<size_t>(kk)]]);
                }
                std::fprintf(stderr, "\n");
            }
            ++k_out;
        }
    }
    return 2;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_reg_aln_strs  (mirrors longcallD align.c / align.h)
// ════════════════════════════════════════════════════════════════════════════

int collect_noisy_reg_aln_strs(const Options& opts, BamChunk& chunk,
                                hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                                const std::vector<int>& noisy_read_ids,
                                const std::vector<uint8_t>& ref_seq_vec,
                                std::array<int, 2>& clu_n_seqs,
                                std::array<std::vector<int>, 2>& clu_read_ids,
                                std::array<std::vector<AlnStr>, 2>& aln_strs) {
    if (noisy_read_ids.empty()) return 0;
    const int ref_seq_len = static_cast<int>(ref_seq_vec.size());

    // 1. Extract per-read sequences / cover info.
    NoisyReadInfo info = collect_noisy_read_info(opts, chunk,
                                                  noisy_reg_beg, noisy_reg_end,
                                                  noisy_read_ids);

    // 2. Sort: BOTH_COVER first, then error rate (sampling mode), then length desc.
    const bool sampling = (noisy_reg_end - noisy_reg_beg + 1) >=
                          static_cast<hts_pos_t>(opts.min_noisy_reg_size_to_sample_reads);
    sort_noisy_region_reads(info, sampling);

    // 3. Find phase set with reads from both haplotypes.
    const hts_pos_t ps = collect_phase_set_with_both_haps(info,
                                                           opts.min_hap_full_reads,
                                                           opts.min_hap_reads);

    int n_full_reads = 0;
    for (int i = 0; i < info.n_reads; ++i)
        if (noisyIsBothCover(info.fully_covers[static_cast<size_t>(i)])) ++n_full_reads;

    int n_cons = 0;
    if (opts.verbose >= 1) {
        std::fprintf(stderr, "BranchSelect ps=%" PRId64 " n_full_reads=%d n_reads=%d\n",
                     static_cast<int64_t>(ps), n_full_reads, info.n_reads);
    }
    // longcallD align.c: collect_ref_read_aln_str = ((refine_bam && out_aln_fp) || out_somatic); pgPhase has no somatic port.
    const bool collect_ref_read_aln_str =
        (opts.refine_aln && !opts.output_aln.empty());

    if (ps > 0) {
        n_cons = wfa_collect_noisy_aln_str_with_ps_hap(
            opts, sampling, info, ps,
            opts.min_hap_full_reads, opts.min_hap_reads,
            ref_seq_vec.data(), ref_seq_len, collect_ref_read_aln_str,
            clu_n_seqs, clu_read_ids, aln_strs);
    } else if (n_full_reads >= opts.min_depth) {
        n_cons = wfa_collect_noisy_aln_str_no_ps_hap(
            opts, info, ref_seq_vec.data(), ref_seq_len, collect_ref_read_aln_str,
            clu_n_seqs, clu_read_ids, aln_strs);
    }

    // longcallD: if (n_cons > 0 && ((refine_bam && out_aln_fp) || out_somatic)) update_digars_from_aln_str(...)
    if (n_cons > 0 && opts.refine_aln && !opts.output_aln.empty()) {
        update_digars_from_aln_str(chunk, noisy_reg_beg, noisy_reg_end,
                                   info.read_id_to_full_covers,
                                   info.read_reg_beg, info.read_reg_end,
                                   aln_strs, n_cons, clu_n_seqs, clu_read_ids);
    }
    return n_cons;
}

} // namespace pgphase_collect
