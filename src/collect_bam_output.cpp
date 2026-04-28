#include "collect_bam_output.hpp"

#include <algorithm>
#include <cinttypes>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>

#include <htslib/kstring.h>

namespace pgphase_collect {

namespace {

[[noreturn]] inline void fatal_bam_output(const std::string& msg) {
    std::cerr << msg << "\n";
    std::abort();
}

inline int sam_realloc_bam_data_compat(bam1_t* b, int desired_size) {
    if (desired_size <= static_cast<int>(b->m_data)) return 0;
    const int new_cap = desired_size;
    uint8_t* p = static_cast<uint8_t*>(realloc(b->data, static_cast<size_t>(new_cap)));
    if (p == nullptr) return -1;
    b->data = p;
    b->m_data = static_cast<uint32_t>(new_cap);
    return 0;
}

inline int digar_query_len(const std::vector<DigarOp>& digars) {
    int qlen = 0;
    for (const DigarOp& d : digars) {
        switch (d.type) {
            case DigarType::Equal:
            case DigarType::Snp:
            case DigarType::Insertion:
            case DigarType::SoftClip:
                qlen += d.len;
                break;
            default:
                break;
        }
    }
    return qlen;
}

inline uint32_t digar_to_cigar_op(DigarType type) {
    switch (type) {
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

inline int cigar_is_identical(const uint32_t* cigar, int n_cigar, const uint32_t* new_cigar, int new_n_cigar) {
    if (n_cigar != new_n_cigar) return 0;
    for (int i = 0; i < n_cigar; ++i) {
        if (cigar[i] != new_cigar[i]) return 0;
    }
    return 1;
}

inline int get_nm_from_digar(const std::vector<DigarOp>& digars) {
    int nm = 0;
    for (const DigarOp& d : digars) {
        if (d.type == DigarType::Snp || d.type == DigarType::Insertion || d.type == DigarType::Deletion) {
            nm += d.len;
        }
    }
    return nm;
}

inline char nt_char_upper(char b) {
    switch (std::toupper(static_cast<unsigned char>(b))) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return std::toupper(static_cast<unsigned char>(b));
        default:
            return 'N';
    }
}

inline char nt_char_lower(char b) {
    char u = nt_char_upper(b);
    if (u == 'N') return 'n';
    return static_cast<char>(std::tolower(static_cast<unsigned char>(u)));
}

inline kstring_t* get_md_from_digar(const std::vector<DigarOp>& digars,
                                    const std::string& ref_seq,
                                    hts_pos_t ref_beg,
                                    hts_pos_t ref_end) {
    kstring_t* md = static_cast<kstring_t*>(calloc(1, sizeof(kstring_t)));
    hts_pos_t pos = 0;
    int eq_len = 0;
    for (const DigarOp& d : digars) {
        const int len = d.len;
        if (d.type == DigarType::Equal) {
            eq_len += len;
        } else if (d.type == DigarType::Snp) {
            kputw(eq_len, md);
            for (int j = 0; j < len; ++j) {
                pos = d.pos + j;
                char ref_base = 'N';
                if (pos >= ref_beg && pos < ref_end) ref_base = nt_char_upper(ref_seq[static_cast<size_t>(pos - ref_beg)]);
                kputc(ref_base, md);
            }
            eq_len = 0;
        } else if (d.type == DigarType::Deletion) {
            kputw(eq_len, md);
            kputc('^', md);
            for (int j = 0; j < len; ++j) {
                pos = d.pos + j;
                char ref_base = 'N';
                if (pos >= ref_beg && pos < ref_end) ref_base = nt_char_upper(ref_seq[static_cast<size_t>(pos - ref_beg)]);
                kputc(ref_base, md);
            }
            eq_len = 0;
        }
    }
    if (eq_len > 0) kputw(eq_len, md);
    return md;
}

inline kstring_t* get_cs_from_digar(const std::vector<DigarOp>& digars,
                                    const std::string& ref_seq,
                                    hts_pos_t ref_beg,
                                    hts_pos_t ref_end) {
    kstring_t* cs = static_cast<kstring_t*>(calloc(1, sizeof(kstring_t)));
    hts_pos_t pos = 0;
    for (const DigarOp& d : digars) {
        const int len = d.len;
        if (d.type == DigarType::Equal) {
            kputc(':', cs);
            kputw(len, cs);
        } else if (d.type == DigarType::Snp) {
            for (int j = 0; j < len; ++j) {
                kputc('*', cs);
                char ref_base = 'n';
                if (d.pos + j >= ref_beg && d.pos + j < ref_end) {
                    ref_base = nt_char_lower(ref_seq[static_cast<size_t>(d.pos - ref_beg + j)]);
                }
                kputc(ref_base, cs);
                char alt_base = 'n';
                if (!d.alt.empty() && j < static_cast<int>(d.alt.size())) alt_base = nt_char_lower(d.alt[static_cast<size_t>(j)]);
                kputc(alt_base, cs);
            }
        } else if (d.type == DigarType::Insertion) {
            kputc('+', cs);
            for (int j = 0; j < len; ++j) {
                char ins_base = 'n';
                if (!d.alt.empty() && j < static_cast<int>(d.alt.size())) ins_base = nt_char_lower(d.alt[static_cast<size_t>(j)]);
                kputc(ins_base, cs);
            }
        } else if (d.type == DigarType::Deletion) {
            kputc('-', cs);
            for (int j = 0; j < len; ++j) {
                pos = d.pos + j;
                char ref_base = 'n';
                if (pos >= ref_beg && pos < ref_end) ref_base = nt_char_lower(ref_seq[static_cast<size_t>(pos - ref_beg)]);
                kputc(ref_base, cs);
            }
        }
    }
    return cs;
}

inline int update_bam1_tags(bam1_t* b,
                            const std::vector<DigarOp>& digars,
                            const std::string& ref_seq,
                            hts_pos_t ref_beg,
                            hts_pos_t ref_end) {
    uint8_t* nm_tag = bam_aux_get(b, "NM");
    if (nm_tag != nullptr) {
        const int nm = get_nm_from_digar(digars);
        if (nm != bam_aux2i(nm_tag)) {
            bam_aux_del(b, nm_tag);
            bam_aux_append(b, "NM", 'i', 4, reinterpret_cast<const uint8_t*>(&nm));
        }
    }
    uint8_t* md_tag = bam_aux_get(b, "MD");
    if (md_tag != nullptr) {
        kstring_t* md = get_md_from_digar(digars, ref_seq, ref_beg, ref_end);
        if (strcmp(md->s, reinterpret_cast<char*>(bam_aux2Z(md_tag))) != 0) {
            bam_aux_del(b, md_tag);
            bam_aux_append(b, "MD", 'Z', md->l + 1, reinterpret_cast<uint8_t*>(md->s));
        }
        free(md->s);
        free(md);
    }
    uint8_t* cs_tag = bam_aux_get(b, "cs");
    if (cs_tag != nullptr) {
        kstring_t* cs = get_cs_from_digar(digars, ref_seq, ref_beg, ref_end);
        if (strcmp(cs->s, reinterpret_cast<char*>(bam_aux2Z(cs_tag))) != 0) {
            bam_aux_del(b, cs_tag);
            bam_aux_append(b, "cs", 'Z', cs->l + 1, reinterpret_cast<uint8_t*>(cs->s));
        }
        free(cs->s);
        free(cs);
    }
    return 0;
}

inline int refine_bam1(const ReadRecord& read,
                       const std::string& ref_seq,
                       hts_pos_t ref_beg,
                       hts_pos_t ref_end,
                       bam1_t* b) {
    const hts_pos_t new_pos = read.digars.empty() ? 0 : read.digars.front().pos;
    if (new_pos < 1) {
        std::fprintf(stderr, "Error: Invalid position %" PRId64 " for read %s\n", new_pos, bam_get_qname(b));
        return -1;
    }
    if (new_pos != b->core.pos + 1) b->core.pos = new_pos - 1;

    uint32_t* new_cigar = static_cast<uint32_t*>(malloc(read.digars.size() * sizeof(uint32_t)));
    int new_n_cigar = 0;
    int new_m_cigar = static_cast<int>(read.digars.size());
    for (const DigarOp& d : read.digars) {
        const uint32_t op = digar_to_cigar_op(d.type);
        if (new_n_cigar > 0 && bam_cigar_op(new_cigar[new_n_cigar - 1]) == op) {
            const uint32_t len = bam_cigar_oplen(new_cigar[new_n_cigar - 1]) + static_cast<uint32_t>(d.len);
            new_cigar[new_n_cigar - 1] = bam_cigar_gen(len, op);
        } else {
            if (new_n_cigar == new_m_cigar) {
                new_m_cigar = std::max(2 * new_m_cigar, 1);
                new_cigar = static_cast<uint32_t*>(realloc(new_cigar, static_cast<size_t>(new_m_cigar) * sizeof(uint32_t)));
            }
            new_cigar[new_n_cigar++] = bam_cigar_gen(static_cast<uint32_t>(d.len), op);
        }
    }
    if (!(b->core.flag & BAM_FSUPPLEMENTARY) && new_n_cigar > 0) {
        if ((new_cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
            new_cigar[0] = (new_cigar[0] & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
        }
        if ((new_cigar[new_n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
            new_cigar[new_n_cigar - 1] = (new_cigar[new_n_cigar - 1] & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
        }
    }

    uint32_t* old_cigar = bam_get_cigar(b);
    int old_n_cigar = b->core.n_cigar;
    if (cigar_is_identical(old_cigar, old_n_cigar, new_cigar, new_n_cigar)) {
        free(new_cigar);
        return 0;
    }
    {
        const int qname_len = b->core.l_qname;
        const int old_data_len = b->l_data;
        const int new_data_len = old_data_len + (new_n_cigar - old_n_cigar) * static_cast<int>(sizeof(uint32_t));
        if (new_data_len > static_cast<int>(b->m_data)) {
            sam_realloc_bam_data_compat(b, new_data_len);
        }
        uint8_t* new_data = static_cast<uint8_t*>(malloc(static_cast<size_t>(b->m_data)));
        memcpy(new_data, b->data, static_cast<size_t>(qname_len));
        memcpy(new_data + qname_len, new_cigar, static_cast<size_t>(new_n_cigar) * sizeof(uint32_t));
        const int rest_len = old_data_len - (qname_len + old_n_cigar * static_cast<int>(sizeof(uint32_t)));
        memcpy(new_data + qname_len + new_n_cigar * static_cast<int>(sizeof(uint32_t)),
               b->data + qname_len + old_n_cigar * static_cast<int>(sizeof(uint32_t)),
               static_cast<size_t>(rest_len));
        free(b->data);
        free(new_cigar);
        b->data = new_data;
        b->core.n_cigar = static_cast<uint32_t>(new_n_cigar);
        b->l_data = new_data_len;
    }
    update_bam1_tags(b, read.digars, ref_seq, ref_beg, ref_end);
    return 0;
}

inline void update_hp_ps_tags(bam1_t* rec, int hap, hts_pos_t ps) {
    uint8_t* hp_tag = bam_aux_get(rec, "HP");
    if (hap != 0) {
        if (hp_tag != nullptr) {
            if (hap != bam_aux2i(hp_tag)) {
                bam_aux_del(rec, hp_tag);
                bam_aux_append(rec, "HP", 'i', 4, reinterpret_cast<const uint8_t*>(&hap));
            }
        } else {
            bam_aux_append(rec, "HP", 'i', 4, reinterpret_cast<const uint8_t*>(&hap));
        }
    } else {
        if (hp_tag != nullptr) bam_aux_del(rec, hp_tag);
    }

    uint8_t* ps_tag = bam_aux_get(rec, "PS");
    if (ps > 0) {
        if (ps_tag != nullptr) {
            if (ps != bam_aux2i(ps_tag)) {
                bam_aux_del(rec, ps_tag);
                bam_aux_append(rec, "PS", 'i', 4, reinterpret_cast<const uint8_t*>(&ps));
            }
        } else {
            bam_aux_append(rec, "PS", 'i', 4, reinterpret_cast<const uint8_t*>(&ps));
        }
    } else {
        if (ps_tag != nullptr) bam_aux_del(rec, ps_tag);
    }
}

inline void strip_hp_ps_tags(bam1_t* rec) {
    uint8_t* hp = bam_aux_get(rec, "HP");
    if (hp != nullptr) bam_aux_del(rec, hp);
    uint8_t* ps = bam_aux_get(rec, "PS");
    if (ps != nullptr) bam_aux_del(rec, ps);
}

} // namespace

PhasedAlignmentWriter::PhasedAlignmentWriter(const Options& opts, const bam_hdr_t* header) : opts_(opts) {
    if (opts_.output_aln.empty()) return;
    if (opts_.bam_files.empty()) fatal_bam_output("No input BAM/CRAM files provided.");
    const char* mode = "wb";
    if (opts_.output_aln_format == OutputAlignmentFormat::Sam) mode = "w";
    if (opts_.output_aln_format == OutputAlignmentFormat::Cram) mode = "wc";
    out_ = hts_open(opts_.output_aln.c_str(), mode);
    if (out_ == nullptr) fatal_bam_output("Failed to open alignment output file '" + opts_.output_aln + "'");
    if (opts_.output_aln_format == OutputAlignmentFormat::Cram && hts_set_fai_filename(out_, opts_.ref_fasta.c_str()) != 0) {
        hts_close(out_);
        out_ = nullptr;
        fatal_bam_output("Failed to set reference file for output CRAM encoding: " + opts_.ref_fasta);
    }
    hts_set_threads(out_, opts_.threads);

    header_.reset(sam_hdr_dup(header));
    if (!header_) {
        hts_close(out_);
        out_ = nullptr;
        fatal_bam_output("Failed to duplicate BAM header.");
    }
    const char* cl = opts_.command_line.empty() ? "pgphase collect-bam-variation" : opts_.command_line.c_str();
    if (sam_hdr_add_pg(header_.get(), "pgphase", "VN", "0.0.0", "CL", cl, nullptr) < 0) {
        hts_close(out_);
        out_ = nullptr;
        fatal_bam_output("Fail to add PG line to bam header.");
    }
    if (sam_hdr_write(out_, header_.get()) < 0) {
        hts_close(out_);
        out_ = nullptr;
        fatal_bam_output("Failed to write BAM header.");
    }

    in_bams_.reserve(opts_.bam_files.size());
    in_headers_.reserve(opts_.bam_files.size());
    in_indexes_.reserve(opts_.bam_files.size());
    for (const std::string& path : opts_.bam_files) {
        in_bams_.push_back(std::make_unique<SamFile>(path, 1, opts_.ref_fasta));
        in_headers_.emplace_back(sam_hdr_read(in_bams_.back()->get()));
        if (!in_headers_.back()) fatal_bam_output("Failed to read alignment file header '" + path + "'");
        in_indexes_.emplace_back(sam_index_load(in_bams_.back()->get(), path.c_str()));
        if (!in_indexes_.back()) fatal_bam_output("Failed to load index for '" + path + "'");
    }
}

PhasedAlignmentWriter::~PhasedAlignmentWriter() {
    if (out_ != nullptr) hts_close(out_);
}

int PhasedAlignmentWriter::write_chunks(const std::vector<BamChunk>& chunks) {
    if (out_ == nullptr) return 0;
    int n_out_reads = 0;
    for (const BamChunk& chunk : chunks) {
        const hts_pos_t reg_beg = chunk.region.beg;
        const hts_pos_t reg_end = chunk.region.end;
        const int tid = chunk.region.tid;
        int global_read_i = 0;
        for (size_t bi = 0; bi < in_bams_.size(); ++bi) {
            bam_hdr_t* in_header = in_headers_[bi].get();
            hts_idx_t* in_idx = in_indexes_[bi].get();
            hts_itr_t* iter = sam_itr_queryi(in_idx, tid, reg_beg - 1, reg_end);
            if (iter == nullptr) {
                const char* tname = (tid >= 0 && tid < in_header->n_targets) ? in_header->target_name[tid] : ".";
                std::cerr << "Failed to create iterator for region: "
                          << tname << ":" << reg_beg << "-" << reg_end << ". Skipping.\n";
                continue;
            }
            bam1_t* read = bam_init1();
            int read_i = 0;
            int skip_read_i = 0;
            while (sam_itr_next(in_bams_[bi]->get(), iter, read) >= 0) {
                if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) ||
                    read->core.qual < opts_.min_mapq) {
                    if (skip_read_i >= chunk.n_up_ovlp_skip_reads[bi]) {
                        strip_hp_ps_tags(read);
                        if (sam_write1(out_, in_header, read) < 0) {
                            const char* tname = (tid >= 0 && tid < in_header->n_targets) ? in_header->target_name[tid] : ".";
                            std::cerr << "Failed to write BAM record. "
                                      << tname << ":" << (read->core.pos + 1) << " " << bam_get_qname(read) << "\n";
                            std::abort();
                        }
                        ++n_out_reads;
                    }
                    ++skip_read_i;
                    continue;
                }
                if (read_i >= static_cast<int>(chunk.up_ovlp_read_i[bi].size())) {
                    const ReadRecord& rr = chunk.reads[static_cast<size_t>(global_read_i)];
                    const int hap = static_cast<size_t>(global_read_i) < chunk.haps.size() ? chunk.haps[static_cast<size_t>(global_read_i)] : 0;
                    const hts_pos_t ps = static_cast<size_t>(global_read_i) < chunk.phase_sets.size() ? chunk.phase_sets[static_cast<size_t>(global_read_i)] : static_cast<hts_pos_t>(-1);
                    if (opts_.refine_aln) {
                        if (digar_query_len(rr.digars) != read->core.l_qseq) {
                            const char* tname = (tid >= 0 && tid < in_header->n_targets) ? in_header->target_name[tid] : ".";
                            std::cerr << "Read length mismatch when writing BAM record: "
                                      << tname << ":" << (read->core.pos + 1) << " " << bam_get_qname(read)
                                      << ", digar_qlen: " << digar_query_len(rr.digars)
                                      << ", bam_qlen: " << read->core.l_qseq << "\n";
                            std::abort();
                        }
                        if (refine_bam1(rr, chunk.ref_seq, chunk.ref_beg, chunk.ref_end, read) != 0) {
                            const char* tname = (tid >= 0 && tid < in_header->n_targets) ? in_header->target_name[tid] : ".";
                            std::cerr << "Failed to refine BAM record: "
                                      << tname << ":" << (read->core.pos + 1) << " " << bam_get_qname(read) << "\n";
                            std::abort();
                        }
                    }
                    update_hp_ps_tags(read, hap, ps);
                    if (sam_write1(out_, in_header, read) < 0) {
                        const char* tname = (tid >= 0 && tid < in_header->n_targets) ? in_header->target_name[tid] : ".";
                        std::cerr << "Failed to write BAM record. "
                                  << tname << ":" << (read->core.pos + 1) << " " << bam_get_qname(read) << "\n";
                        std::abort();
                    }
                    ++n_out_reads;
                }
                ++read_i;
                ++global_read_i;
            }
            bam_destroy1(read);
            hts_itr_destroy(iter);
        }
    }
    return n_out_reads;
}

} // namespace pgphase_collect
