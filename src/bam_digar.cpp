#include "bam_digar.hpp"

#include "collect_var.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Noisy sliding window
// ════════════════════════════════════════════════════════════════════════════

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

/** Return the per-technology default slide window, or opts.noisy_reg_slide_win if overridden. */
static int effective_noisy_slide_win(const Options& opts) {
    if (opts.noisy_reg_slide_win > 0) return opts.noisy_reg_slide_win;
    if (opts.is_ont()) return kDefaultNoisyRegSlideWinOnt;
    if (opts.is_short_reads()) return kDefaultNoisyRegSlideWinShortReads;
    return kDefaultNoisyRegSlideWinHifi;
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

    bool suppress_left_clip = false;
    bool suppress_right_clip = false;

    void add_end_clip_region(int cigar_idx, int n_cigar, hts_pos_t ref_pos, hts_pos_t tlen, int clip_len) {
        if (clip_len < kLongClipLength) return;
        if (cigar_idx == 0 && ref_pos > 10) { // left end clip
            if (!suppress_left_clip)
                noisy_out.push_back(Interval{ref_pos, std::min<hts_pos_t>(tlen, ref_pos + kClipFlank), clip_len});
        } else if (cigar_idx == n_cigar - 1 && ref_pos < tlen - 10) { // right end clip
            if (!suppress_right_clip)
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
// ════════════════════════════════════════════════════════════════════════════
// Read sequence & quality utilities
// ════════════════════════════════════════════════════════════════════════════

/** Increment read.total_cand_events for one SNP/INS/DEL.
 *  total_cand_events counts all candidate-quality events including long-clip noisy windows,
 *  used by read_has_too_many_variants to decide whether to skip the read. */
static inline void record_variant_event(ReadRecord& read, int tid, const DigarOp& digar) {
    (void)tid;
    (void)digar;
    read.total_cand_events++;
}

/** Parse a "CHR:POS" debug-site string; returns false if malformed. */
inline bool parse_debug_site(const std::string& site, std::string& chrom_out, hts_pos_t& pos_out) {
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

/** If opts.debug_site matches this read's span, dump per-digar hits to stderr for manual inspection. */
inline void maybe_dump_debug_site(const Options& opts, const bam_hdr_t* header, const ReadRecord& read) {
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
                  << "skipped=" << (read.is_skipped ? 1 : 0) << "\t"
                  << "mapped=" << mapped_len << "\tnoisy=" << noisy_len << "\t"
                  << "cand=" << read.total_cand_events << "\t"
                  << "op=" << static_cast<int>(d.type) << "\t"
                  << "len=" << d.len << "\talt=\"" << alt << "\"\n";
    }
}


/** Return the base at query position qi as an uppercase char, 'N' if ambiguous or out of range. */
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

std::vector<uint8_t> copy_qualities(const bam1_t* aln) {
    const uint8_t* qual = bam_get_qual(aln);
    return std::vector<uint8_t>(qual, qual + aln->core.l_qseq);
}

int aux_int_or_default(const bam1_t* aln, const char tag[2], int default_value) {
    uint8_t* data = bam_aux_get(const_cast<bam1_t*>(aln), tag);
    return data == nullptr ? default_value : bam_aux2i(data);
}

/** An insertion is low-quality only if ALL of its bases are below min_bq (longcallD behaviour). */
bool insertion_is_low_quality(const bam1_t* aln, int qi, int len, int min_bq) {
    for (int i = 0; i < len; ++i) {
        if (base_quality(aln, qi + i) >= min_bq) return false;
    }
    return true;
}

/** A deletion is low-quality if BOTH flanking bases (qi-1 and qi) are below min_bq. */
bool deletion_is_low_quality(const bam1_t* aln, int qi, int min_bq) {
    const bool left_ok = qi == 0 || base_quality(aln, qi - 1) >= min_bq;
    const bool right_ok = qi >= aln->core.l_qseq || base_quality(aln, qi) >= min_bq;
    return !(left_ok && right_ok);
}
// ════════════════════════════════════════════════════════════════════════════
// Digar building
// ════════════════════════════════════════════════════════════════════════════

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
void recompute_chunk_qual_stats(BamChunk& chunk) {
    int64_t qual_hist[256] = {};
    for (const ReadRecord& r : chunk.reads) {
        if (r.is_skipped) continue;
        for (uint8_t q : r.qual) {
            qual_hist[static_cast<size_t>(q)]++;
        }
    }
    int64_t n_total = 0;
    for (int i = 0; i < 256; ++i) n_total += qual_hist[i];
    std::vector<int> valid_quals;
    valid_quals.reserve(256);
    for (int i = 0; i < 256; ++i) {
        const int64_t c = qual_hist[i];
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

/** Contig length from header; 0 if header/tid is invalid. Used to bound end-clip noisy windows. */
static hts_pos_t contig_len_from_header(const bam_hdr_t* header, int tid) {
    if (header == nullptr || tid < 0 || tid >= header->n_targets) return 0;
    return static_cast<hts_pos_t>(header->target_len[tid]);
}

/** True if CIGAR uses =/X ops (no M): selects the EQX parser path which needs no reference. */
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

/** Merge adjacent noisy intervals after digar construction. */
static void digar_merge_noisy(ReadRecord& read) {
    merge_intervals(read.noisy_regions);
}

/** Build digars by comparing each M base against the reference (slowest path; no MD/cs tag). */
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
    noisy_builder.suppress_left_clip  = read.is_ont_palindrome && bam_is_rev(aln);
    noisy_builder.suppress_right_clip = read.is_ont_palindrome && !bam_is_rev(aln);

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

/** Build digars from =/X CIGAR ops; no reference lookup needed since mismatches are explicit. */
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
    noisy_builder.suppress_left_clip  = read.is_ont_palindrome && bam_is_rev(aln);
    noisy_builder.suppress_right_clip = read.is_ont_palindrome && !bam_is_rev(aln);

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

/** Build digars from M-CIGAR + MD tag; faster than ref lookup but requires a valid MD string.
 *  Returns false (caller falls back to ref path) if MD is absent, malformed, or inconsistent. */
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
    noisy_builder.suppress_left_clip  = read.is_ont_palindrome && bam_is_rev(aln);
    noisy_builder.suppress_right_clip = read.is_ont_palindrome && !bam_is_rev(aln);
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

/** Build digars from the minimap2 cs tag (handles :, =, *, -, +, ~ cs operations).
 *  Returns false (caller falls back to ref path) if cs is absent or malformed. */
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
    noisy_builder.suppress_left_clip  = read.is_ont_palindrome && bam_is_rev(aln);
    noisy_builder.suppress_right_clip = read.is_ont_palindrome && !bam_is_rev(aln);
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
    read.noisy_regions.clear();
    read.total_cand_events = 0;
    read.digars.reserve(static_cast<size_t>(aln->core.n_cigar) * 2u + 8u);

    if (cigar_has_eqx_without_m(aln)) {
        build_digars_eqx_cigar(aln, header, opts, read);
    } else if (bam_aux_get(const_cast<bam1_t*>(aln), "cs") != nullptr) {
        if (!build_digars_cs_tag(aln, header, opts, read)) {
            read.digars.clear();
            read.noisy_regions.clear();
            build_digars_ref_cigar(aln, header, ref, opts, read);
        }
    } else if (bam_aux_get(const_cast<bam1_t*>(aln), "MD") != nullptr) {
        if (!build_digars_md_cigar(aln, header, opts, read)) {
            read.digars.clear();
            read.noisy_regions.clear();
            build_digars_ref_cigar(aln, header, ref, opts, read);
        }
    } else {
        build_digars_ref_cigar(aln, header, ref, opts, read);
    }

    digar_merge_noisy(read);
}

// ════════════════════════════════════════════════════════════════════════════
// Read filtering
// ════════════════════════════════════════════════════════════════════════════

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

static bool overlaps_region(int read_tid,
                            hts_pos_t read_beg,
                            hts_pos_t read_end,
                            int region_tid,
                            hts_pos_t region_beg,
                            hts_pos_t region_end) {
    if (region_tid < 0 || read_tid != region_tid) return false;
    return !(read_end < region_beg || read_beg > region_end);
}

bool read_overlaps_prev_region(const RegionChunk& chunk, int read_tid, hts_pos_t read_beg, hts_pos_t read_end) {
    return overlaps_region(read_tid, read_beg, read_end, chunk.prev_tid, chunk.prev_beg, chunk.prev_end);
}

bool read_overlaps_next_region(const RegionChunk& chunk, int read_tid, hts_pos_t read_beg, hts_pos_t read_end) {
    return overlaps_region(read_tid, read_beg, read_end, chunk.next_tid, chunk.next_beg, chunk.next_end);
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

// ════════════════════════════════════════════════════════════════════════════
// ONT palindrome detection
// ════════════════════════════════════════════════════════════════════════════

/** longcallD `is_ont_palindrome_clip` (bam_utils.c): detect palindromic ONT reads via SA tag.
 *  A read is palindromic when a supplementary alignment on the same contig overlaps the primary
 *  span, indicating the read folded back on itself and mapped twice to the same region. */
static bool detect_ont_palindrome(const bam1_t* aln, const bam_hdr_t* header, const Options& opts) {
    if (!opts.is_ont()) return false;
    const uint8_t* sa_raw = bam_aux_get(aln, "SA");
    if (!sa_raw) return false;
    const char* sa = bam_aux2Z(sa_raw);
    if (!sa || !*sa) return false;

    const hts_pos_t primary_pos = aln->core.pos + 1;
    const hts_pos_t primary_end = bam_endpos(aln);
    const int primary_tid = aln->core.tid;
    const char* primary_chrom = (header && primary_tid >= 0 && primary_tid < header->n_targets)
                                    ? header->target_name[primary_tid] : nullptr;

    char rname[512], cigar_buf[512];
    char strand = '+';
    int sa_pos = 0, mapq = 0, nm = 0;
    const char* p = sa;
    while (*p) {
        const char* entry_end = std::strchr(p, ';');
        if (!entry_end) entry_end = p + std::strlen(p);
        rname[0] = cigar_buf[0] = '\0';
        if (std::sscanf(p, "%511[^,],%d,%c,%511[^,],%d,%d",
                        rname, &sa_pos, &strand, cigar_buf, &mapq, &nm) == 6) {
            if (!primary_chrom || std::strcmp(rname, primary_chrom) == 0) {
                hts_pos_t sa_end = static_cast<hts_pos_t>(sa_pos);
                for (const char* cp = cigar_buf; *cp; ) {
                    const int len = static_cast<int>(std::strtol(cp, const_cast<char**>(&cp), 10));
                    if (*cp == 'M' || *cp == 'D' || *cp == '=' || *cp == 'X') sa_end += len;
                    if (*cp) ++cp;
                }
                if (sa_end >= primary_pos && static_cast<hts_pos_t>(sa_pos) <= primary_end) {
                    return true;
                }
            }
        }
        if (!*entry_end) break;
        p = entry_end + 1;
    }
    return false;
}

// ════════════════════════════════════════════════════════════════════════════
// BAM loading
// ════════════════════════════════════════════════════════════════════════════

std::vector<ReadRecord> load_read_records_for_chunk(const Options& opts,
                                                    const RegionChunk& chunk,
                                                    int input_index,
                                                    SamFile& bam,
                                                    bam_hdr_t* header,
                                                    const hts_idx_t* index,
                                                    ReferenceCache& ref,
                                                    OverlapSkipCounts* overlap_skip_counts) {
    auto iter = make_chunk_iterator(index, chunk);
    if (!iter) throw std::runtime_error("failed to create BAM iterator for chunk");

    std::vector<ReadRecord> reads;
    reads.reserve(4096);
    std::unique_ptr<bam1_t, AlignmentDeleter> aln(bam_init1());
    while (sam_itr_next(bam.get(), iter.get(), aln.get()) >= 0) {
        const int read_tid = aln->core.tid;
        const hts_pos_t read_beg = aln->core.pos + 1;
        const hts_pos_t read_end = bam_endpos(aln.get());
        const bool overlaps_prev = read_overlaps_prev_region(chunk, read_tid, read_beg, read_end);
        const bool overlaps_next = read_overlaps_next_region(chunk, read_tid, read_beg, read_end);

        if (!read_passes_filters(aln.get(), opts)) {
            if (overlap_skip_counts != nullptr) {
                if (overlaps_prev) ++overlap_skip_counts->upstream;
                if (overlaps_next) ++overlap_skip_counts->downstream;
            }
            continue;
        }

        ReadRecord read;
        read.tid = read_tid;
        read.input_index = input_index;
        read.beg = read_beg;
        read.end = read_end;
        read.reverse = bam_is_rev(aln.get());
        read.nm = aux_int_or_default(aln.get(), "NM", 0);
        read.mapq = static_cast<int>(aln->core.qual);
        read.qname = bam_get_qname(aln.get());
        read.qual = copy_qualities(aln.get());
        read.alignment.reset(bam_dup1(aln.get()));
        if (!read.alignment) throw std::runtime_error("failed to duplicate BAM alignment");
        read.is_ont_palindrome = detect_ont_palindrome(aln.get(), header, opts);
        build_digars_and_events(aln.get(), header, ref, opts, read);
        read.is_skipped = read_has_too_many_variants(read, opts);
        maybe_dump_debug_site(opts, header, read);
        reads.push_back(std::move(read));
    }
    return reads;
}

// ════════════════════════════════════════════════════════════════════════════
// Chunk finalization
// ════════════════════════════════════════════════════════════════════════════

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

void finalize_bam_chunk(BamChunk& chunk, ReferenceCache& ref, const bam_hdr_t* header) {
    recompute_chunk_qual_stats(chunk);
    populate_reference_slice(chunk, ref, header);
    populate_low_complexity_intervals(chunk);
    populate_chunk_read_indexes(chunk);
}

} // namespace pgphase_collect
