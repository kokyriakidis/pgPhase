/**
 * @file collect_output.cpp
 * @brief TSV, optional VCF, and read-support writers for collect-bam-variation candidates.
 *
 * @details Outputs describe pre-phasing candidates (not diploid genotypes). Opening the primary BAM
 * for SQ names may throw `std::runtime_error` on I/O failure.
 */

#include "collect_output.hpp"

#include "collect_phase.hpp"

#include <ctime>
#include <fstream>
#include <sstream>

namespace pgphase_collect {

/**
 * @brief Returns the string representation for a VariantType.
 *
 * Used for SNP / INS / DEL labels in TSV, VCF INFO, and read-support TYPE columns.
 *
 * @param type The variant type enum value.
 * @return A string literal ("SNP", "INS", "DEL", or "UNKNOWN").
 */
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

/**
 * @brief Returns the string representation for a VariantCategory.
 *
 * Provides LongcallD-style category tokens (e.g., CLEAN_HET_SNP, LOW_COV)
 * for TSV and VCF INFO.CAT columns.
 *
 * @param category The variant category enum value.
 * @return A string corresponding to the category label.
 */
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
        case VariantCategory::NoisyCandHet:
            return "NOISY_CAND_HET";
        case VariantCategory::NoisyCandHom:
            return "NOISY_CAND_HOM";
        case VariantCategory::NoisyResolved:
            return "NOISY_RESOLVED";
        case VariantCategory::RepeatHetIndel:
            return "REP_HET_INDEL";
        case VariantCategory::NonVariant:
            return "NON_VAR";
    }
    return "UNKNOWN";
}

/**
 * @brief Writes the TSV header for read support output.
 *
 * Defines the columns for the `--read-support` output format, which consists
 * of one line per read x candidate observation.
 *
 * @param out Output stream to write the header to.
 */
void write_read_support_header(std::ostream& out) {
    out << "CHROM\tPOS\tTYPE\tREF_LEN\tALT\tQNAME\tIS_ALT\tLOW_QUAL\tREVERSE\tMAPQ\tCHUNK_BEG\tCHUNK_END\n";
}

/**
 * @brief Writes body lines for read support observations.
 *
 * Outputs one chunk's read x site observations including the sequence name
 * (CHROM) resolved from the BAM header.
 *
 * @param out Output stream to write the rows to.
 * @param header BAM header used to resolve reference sequence IDs to names.
 * @param rows A vector of `ReadSupportRow` structures containing the observation data.
 */
void write_phase_read_tsv_header(std::ostream& out) {
    out << "CHUNK_ID\tREG_CHUNK_I\tCHROM\tCHUNK_BEG\tCHUNK_END\t"
        << "QNAME\tREAD_CHROM\tINPUT_IDX\tREAD_BEG\tREAD_END\tMAPQ\tREVERSE\tSKIPPED\t"
        << "HAP\tPHASE_SET\n";
}

void write_phase_read_tsv_rows(std::ostream& out, const bam_hdr_t* header, const BamChunk& chunk) {
    const RegionChunk& reg = chunk.region;
    const char* chunk_chrom =
        (reg.tid >= 0 && reg.tid < header->n_targets) ? header->target_name[reg.tid] : ".";
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        const ReadRecord& r = chunk.reads[i];
        const char* read_chrom =
            (r.tid >= 0 && r.tid < header->n_targets) ? header->target_name[r.tid] : ".";
        const int hap = i < chunk.haps.size() ? chunk.haps[i] : 0;
        const hts_pos_t ps = i < chunk.phase_sets.size() ? chunk.phase_sets[i] : static_cast<hts_pos_t>(-1);
        out << reg.chunk_id << '\t' << reg.reg_chunk_i << '\t' << chunk_chrom << '\t' << reg.beg << '\t'
            << reg.end << '\t' << r.qname << '\t' << read_chrom << '\t' << r.input_index << '\t' << r.beg
            << '\t' << r.end << '\t' << r.mapq << '\t' << (r.reverse ? 1 : 0) << '\t' << (r.is_skipped ? 1 : 0)
            << '\t' << hap << '\t' << ps << '\n';
    }
}

void write_read_support_rows(std::ostream& out,
                             const bam_hdr_t* header,
                             const std::vector<ReadSupportRow>& rows) {
    for (const ReadSupportRow& r : rows) {
        const char* chrom = header->target_name[r.tid];
        out << chrom << '\t' << r.pos << '\t' << type_name(r.type) << '\t' << r.ref_len << '\t'
            << (r.alt.empty() ? "." : r.alt) << '\t' << r.qname << '\t' << r.is_alt << '\t' << r.is_low_qual
            << '\t' << (r.reverse ? 1 : 0) << '\t' << r.mapq << '\t' << r.chunk_beg << '\t' << r.chunk_end
            << '\n';
    }
}

/**
 * @brief Writes the full read support TSV file from all processed chunks.
 *
 * Concatenates per-chunk batches of read observations into the final TSV file specified in `opts`.
 * Order matches the streaming chunk processing.
 *
 * @param opts Configuration options containing the primary BAM file and read support output path.
 * @param by_chunk A vector of vectors, where each inner vector represents a chunk's read observations.
 * @throws std::runtime_error If the BAM header cannot be read or the output file cannot be opened.
 */
void write_read_support_tsv(const Options& opts, const std::vector<std::vector<ReadSupportRow>>& by_chunk) {
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");

    std::ofstream out(opts.read_support_tsv);
    if (!out) throw std::runtime_error("failed to open read support output: " + opts.read_support_tsv);

    write_read_support_header(out);
    for (const std::vector<ReadSupportRow>& batch : by_chunk) {
        write_read_support_rows(out, header.get(), batch);
    }
}

/**
 * @brief Writes the TSV header for the main candidate variant table.
 *
 * Defines columns for ref/alt sequence, depth, strand counts, allele fraction (AF),
 * and category. Trailing columns carry the k-means phasing scaffold (`PHASE_SET`, `HAP_ALT`, `HAP_REF`)
 * when `collect_var_main` has run phasing; otherwise they are zero.
 *
 * @param out Output stream to write the header line to.
 */
void write_variants_tsv_header(std::ostream& out) {
    out << "CHROM\tPOS\tTYPE\tREF\tALT\tDP\tREF_COUNT\tALT_COUNT\tLOW_QUAL_COUNT"
        << "\tFORWARD_REF\tREVERSE_REF\tFORWARD_ALT\tREVERSE_ALT"
        << "\tAF\tCATEGORY\tINIT_CAT\tPHASE_SET\tHAP_ALT\tHAP_REF\n";
}

/**
 * @brief Serializes the mapped candidate variant evaluations to TSV.
 *
 * Writes one row per CandidateVariant, obtaining REF/ALT sequences from the
 * `ReferenceCache` for SNPs and deletions. Category and count fields mirror internal
 * classifications for cross-checks against longcallD-style debugging output.
 *
 * @param out Open file stream targeting a `.tsv`.
 * @param header BAM header for contig names.
 * @param ref Reference cache used to extract sequence for SNPs and deletions.
 * @param variants Validated list of categorized variants.
 */
void write_variants_tsv_records(std::ostream& out,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants) {
    for (const CandidateVariant& candidate : variants) {
        const VariantKey& key = candidate.key;
        const VariantCounts& counts = candidate.counts;
        const std::string chrom = header->target_name[key.tid];
        std::string ref_seq = ".";
        std::string alt_seq = key.alt.empty() ? "." : key.alt;

        if (key.type == VariantType::Snp) {
            ref_seq = std::string(1, ref.base(key.tid, key.pos, header));
        } else if (key.type == VariantType::Deletion) {
            ref_seq = ref.subseq(key.tid, key.pos, key.ref_len, header);
        }

        out << chrom << '\t' << key.pos << '\t' << type_name(key.type) << '\t' << ref_seq << '\t'
            << alt_seq << '\t' << counts.total_cov << '\t' << counts.ref_cov << '\t'
            << counts.alt_cov << '\t' << counts.low_qual_cov << '\t' << counts.forward_ref << '\t'
            << counts.reverse_ref << '\t' << counts.forward_alt << '\t' << counts.reverse_alt << '\t'
            << counts.allele_fraction << '\t' << category_name(counts.category) << '\t'
            << category_name(counts.candvarcate_initial) << '\t' << candidate.phase_set << '\t'
            << candidate.hap_alt << '\t' << candidate.hap_ref << '\n';
    }
}

/**
 * @brief Writes the full set of categorized variants to a TSV file.
 *
 * Opens the primary BAM file to obtain reference sequence names via the header,
 * configures the FASTA cache, and streams the full merged candidate set to `opts.output_tsv`.
 *
 * @param opts Program options containing I/O paths.
 * @param fai FASTA index for the reference genome.
 * @param variants The complete table of evaluated candidate variants.
 * @throws std::runtime_error If the BAM header cannot be read or `opts.output_tsv` cannot be opened.
 */
void write_variants(const Options& opts, faidx_t* fai, const CandidateTable& variants) {
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai);

    std::ofstream out(opts.output_tsv);
    if (!out) throw std::runtime_error("failed to open output: " + opts.output_tsv);

    write_variants_tsv_header(out);
    write_variants_tsv_records(out, header.get(), ref, variants);
}

/**
 * @brief Writes a valid VCF v4.2 header for candidate variants.
 *
 * Emits standard VCF pragmas, contig definitions from the BAM header,
 * and FILTER/INFO lines for candidate (not final genotype) semantics.
 *
 * @param out Output stream to write the VCF header to.
 * @param opts Program options (reserved; currently unused in the header text).
 * @param header BAM header containing reference names and lengths.
 */
void write_variants_vcf_header(std::ostream& out, const Options& opts, const bam_hdr_t* header) {
    (void)opts;
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
    out << "##FILTER=<ID=LowQual,Description=\"Low quality variant\">\n";
    out << "##FILTER=<ID=RefCall,Description=\"Reference call candidate\">\n";
    out << "##FILTER=<ID=NoCall,Description=\"Site has depth=0 resulting in no call\">\n";
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
    out << "##INFO=<ID=CLEAN,Number=0,Type=Flag,Description=\"Clean-region variant (SNP or simple indel in non-repetitive region)\">\n";
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
}

void write_phased_variants_vcf_header(std::ostream& out, const Options& opts, const bam_hdr_t* header) {
    (void)opts;
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
    out << "##FILTER=<ID=LowQual,Description=\"Low quality variant\">\n";
    out << "##FILTER=<ID=RefCall,Description=\"Reference call candidate\">\n";
    out << "##FILTER=<ID=NoCall,Description=\"Site has depth=0 resulting in no call\">\n";
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
    out << "##INFO=<ID=CLEAN,Number=0,Type=Flag,Description=\"Clean-region variant (SNP or simple indel in non-repetitive region)\">\n";
    out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    out << "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
    out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n";
    out << "##INFO=<ID=REFC,Number=1,Type=Integer,Description=\"Reference allele count\">\n";
    out << "##INFO=<ID=ALTC,Number=1,Type=Integer,Description=\"Alternate allele count\">\n";
    out << "##INFO=<ID=LQC,Number=1,Type=Integer,Description=\"Low-quality observation count\">\n";
    out << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Alternate allele fraction\">\n";
    out << "##INFO=<ID=CAT,Number=1,Type=String,Description=\"pgPhase candidate category\">\n";
    out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    out << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set anchor coordinate\">\n";
    for (int32_t tid = 0; tid < header->n_targets; ++tid) {
        out << "##contig=<ID=" << header->target_name[tid] << ",length=" << header->target_len[tid] << ">\n";
    }
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}

namespace {
struct VcfRecordCore {
    hts_pos_t pos = 0;
    std::string ref_seq;
    std::string alt_seq;
    std::string filter;
    std::string info;
};

/** longcallD `vcf_utils.c` `write_var_to_vcf`: skip when `opt->out_amb_base == 0` and any nt ≥ 4 (non-ACGT). */
static bool lcd_vcf_seq_has_non_acgt(const std::string& s) {
    for (unsigned char uc : s) {
        const char c = static_cast<char>(uc);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return true;
    }
    return false;
}

/**
 * longcallD collect_var.c `make_variants`: for INS/DEL, when cand.alt_ref_base != 4 the first
 * emitted ALT byte is the raw consensus anchor (`alt_ref_base`); when it equals 4, the anchor is
 * taken from reference (`nst_nt4_table`). Values > 3 (e.g. abPOA gap 5 or '-' mapped to 5 in LCD's
 * nst_nt4_table) are copied into alt_bases and rejected by vcf_utils.c write_var_to_vcf
 * (`alt_bases[j][k] >= 4` with verbose \"Invalid alt base\") — those variants never appear in LCD VCF.
 */
static bool passes_lcd_write_var_alt_ref_base_gate(const CandidateVariant& candidate) {
    if (candidate.key.type != VariantType::Insertion && candidate.key.type != VariantType::Deletion)
        return true;
    if (candidate.alt_ref_base != 4 && candidate.alt_ref_base > 3) return false;
    return true;
}

static bool passes_longcalld_vcf_amb_base_gate(const Options& opts, const VcfRecordCore& core) {
    if (opts.output_ambiguous_bases) return true;
    if (lcd_vcf_seq_has_non_acgt(core.ref_seq)) return false;
    if (lcd_vcf_seq_has_non_acgt(core.alt_seq)) return false;
    return true;
}

static bool is_longcalld_germline_output_category(VariantCategory category) {
    switch (category) {
        case VariantCategory::CleanHetSnp:
        case VariantCategory::CleanHetIndel:
        case VariantCategory::CleanHom:
        case VariantCategory::NoisyCandHet:
        case VariantCategory::NoisyCandHom:
            return true;
        default:
            return false;
    }
}

/**
 * longcallD `vcf_utils.c` `write_var_to_vcf`: skip when `var.DP < opt->min_dp` or `var.AD[1] < opt->min_alt_dp`
 * (germline branch). `make_variants` sets `var.DP = cand_vars[cand_i].total_cov` only — not including
 * low_qual_cov — so VCF emission must match that gate even though `classify_var_cate` uses
 * `total_cov + low_qual_cov` for LOW_COV classification.
 */
static bool passes_longcalld_vcf_depth_gates(const CandidateVariant& candidate, const Options& opts) {
    return candidate.counts.total_cov >= opts.min_depth &&
           candidate.counts.alt_cov >= opts.min_alt_depth;
}

/** longcallD `make_variants`: `is_clean = (var_i_to_cate & LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE) != 0`. */
static bool lcd_make_variants_is_clean(const CandidateVariant& candidate) {
    return (candidate.lcd_var_i_to_cate & kCandGermlineClean) != 0;
}

static std::pair<int, int> derive_hap_alt_ref_from_consensus(const CandidateVariant& candidate) {
    int c1 = candidate.hap_to_cons_alle[1];
    int c2 = candidate.hap_to_cons_alle[2];
    if (c1 == -1 && c2 == -1) {
        c1 = c2 = candidate.hap_to_cons_alle[0];
    }
    if (c1 == -1) c1 = 0;
    if (c2 == -1) c2 = 0;
    const bool h1_alt = (c1 != 0);
    const bool h2_alt = (c2 != 0);
    if (h1_alt && h2_alt) return {3, 0};
    if (h1_alt && !h2_alt) return {1, 2};
    if (!h1_alt && h2_alt) return {2, 1};
    return {0, 0};
}

static bool is_alt_genotype(const CandidateVariant& candidate) {
    const auto [hap_alt, hap_ref] = derive_hap_alt_ref_from_consensus(candidate);
    (void)hap_ref;
    return hap_alt == 1 || hap_alt == 2 || hap_alt == 3;
}

static std::vector<const CandidateVariant*> project_longcalld_like_vcf_candidates(
    const CandidateTable& variants,
    const Options& opts,
    bool require_alt_genotype) {
    std::vector<const CandidateVariant*> projected;
    projected.reserve(variants.size());
    for (const CandidateVariant& candidate : variants) {
        if (!is_longcalld_germline_output_category(candidate.counts.category)) continue;
        if (!passes_longcalld_vcf_depth_gates(candidate, opts)) continue;
        if (require_alt_genotype && !is_alt_genotype(candidate)) continue;
        if (!candidate.lcd_make_variants_region_pass) continue;
        projected.push_back(&candidate);
    }
    return projected;
}

static VcfRecordCore build_vcf_record_core(const CandidateVariant& candidate,
                                           const Options& opts,
                                           const bam_hdr_t* header,
                                           ReferenceCache& ref) {
    const VariantKey& key = candidate.key;
    const VariantCounts& counts = candidate.counts;
    (void)header;

    VcfRecordCore core;
    core.pos = key.pos;
    if (key.type == VariantType::Snp) {
        core.ref_seq = std::string(1, ref.base(key.tid, key.pos, header));
        core.alt_seq = key.alt.empty() ? "." : key.alt;
    } else if (key.type == VariantType::Insertion) {
        const hts_pos_t anchor_pos = std::max<hts_pos_t>(1, key.pos - 1);
        core.pos = anchor_pos;
        const char anchor_base = ref.base(key.tid, anchor_pos, header);
        // Mirrors longcallD make_variants INS anchor branch (raw nt index when != 4).
        const char alt_anchor_base = (candidate.alt_ref_base != 4)
                                         ? static_cast<char>(
                                               "ACGTN"[static_cast<size_t>(candidate.alt_ref_base)])
                                         : anchor_base;
        core.ref_seq = std::string(1, anchor_base);
        core.alt_seq = std::string(1, alt_anchor_base) + key.alt;
    } else { // Deletion
        const hts_pos_t anchor_pos = std::max<hts_pos_t>(1, key.pos - 1);
        core.pos = anchor_pos;
        const char anchor_base = ref.base(key.tid, anchor_pos, header);
        const char alt_anchor_base = (candidate.alt_ref_base != 4)
                                         ? static_cast<char>(
                                               "ACGTN"[static_cast<size_t>(candidate.alt_ref_base)])
                                         : anchor_base;
        const std::string del_seq = ref.subseq(key.tid, key.pos, key.ref_len, header);
        core.ref_seq = std::string(1, anchor_base) + del_seq;
        core.alt_seq = std::string(1, alt_anchor_base);
    }

    core.filter = "PASS";
    if (counts.total_cov == 0) {
        core.filter = "NoCall";
    } else if (counts.category == VariantCategory::NonVariant) {
        core.filter = "RefCall";
    } else if (counts.category == VariantCategory::LowCoverage ||
               counts.category == VariantCategory::LowAlleleFraction ||
               counts.category == VariantCategory::StrandBias) {
        // Keep obvious pre-call failures in LowQual; clean/noisy called candidates remain PASS.
        core.filter = "LowQual";
    }

    const hts_pos_t end_pos = core.pos + static_cast<hts_pos_t>(core.ref_seq.size()) - 1;
    std::ostringstream info;
    info << "END=" << end_pos;
    if (lcd_make_variants_is_clean(candidate)) {
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
    core.info = info.str();
    return core;
}
} // namespace

/**
 * @brief Writes VCF body records for candidate variants.
 *
 * Produces left-normalized records: SNP at `key.pos`; insertions and deletions use the
 * anchor base at POS-1 per VCF convention. FILTER is PASS for called candidate categories,
 * RefCall for non-variants, NoCall when depth is zero, and LowQual for obvious pre-call
 * failures (LOW_COV / LOW_AF / STRAND_BIAS). INFO includes END, optional CLEAN (clean-region
 * categories only), depth and allele fields, AF, CAT, and SVTYPE/SVLEN when |SVLEN| >=
 * `opts.min_sv_len`.
 *
 * @param out Output stream for the `.vcf` body.
 * @param opts Options (e.g. `min_sv_len` for tagging large indels).
 * @param header BAM header for contig names.
 * @param ref Reference sequence cache for REF/ALT bases.
 * @param variants Candidate variants to emit.
 */
void write_variants_vcf_records(std::ostream& out,
                                const Options& opts,
                                const bam_hdr_t* header,
                                ReferenceCache& ref,
                                const CandidateTable& variants) {
    const std::vector<const CandidateVariant*> projected =
        project_longcalld_like_vcf_candidates(variants, opts, true);
    for (const CandidateVariant* candidate_ptr : projected) {
        const CandidateVariant& candidate = *candidate_ptr;
        if (!passes_lcd_write_var_alt_ref_base_gate(candidate)) continue;
        const VariantKey& key = candidate.key;
        const std::string chrom = header->target_name[key.tid];
        const VcfRecordCore core = build_vcf_record_core(candidate, opts, header, ref);
        if (!passes_longcalld_vcf_amb_base_gate(opts, core)) continue;
        out << chrom << '\t' << core.pos << "\t.\t" << core.ref_seq << '\t' << core.alt_seq
            << "\t.\t" << core.filter << '\t' << core.info << '\n';
    }
}

void write_phased_variants_vcf_records(std::ostream& out,
                                       const Options& opts,
                                       const bam_hdr_t* header,
                                       ReferenceCache& ref,
                                       const CandidateTable& variants) {
    const std::vector<const CandidateVariant*> projected =
        project_longcalld_like_vcf_candidates(variants, opts, true);
    for (const CandidateVariant* candidate_ptr : projected) {
        const CandidateVariant& candidate = *candidate_ptr;
        if (!passes_lcd_write_var_alt_ref_base_gate(candidate)) continue;
        const VariantKey& key = candidate.key;
        const std::string chrom = header->target_name[key.tid];
        const VcfRecordCore core = build_vcf_record_core(candidate, opts, header, ref);
        if (!passes_longcalld_vcf_amb_base_gate(opts, core)) continue;

        const auto [hap_alt, hap_ref] = derive_hap_alt_ref_from_consensus(candidate);
        std::string gt = "./.";
        if (hap_alt == 1 && hap_ref == 2) gt = "1|0";
        else if (hap_alt == 2 && hap_ref == 1) gt = "0|1";
        else if (hap_alt == 3) gt = "1|1";
        else if (candidate.counts.category == VariantCategory::NonVariant) gt = "0/0";

        std::string ps = ".";
        if ((gt == "1|0" || gt == "0|1") && candidate.phase_set > 0) {
            ps = std::to_string(candidate.phase_set);
        }
        out << chrom << '\t' << core.pos << "\t.\t" << core.ref_seq << '\t' << core.alt_seq
            << "\t.\t" << core.filter << '\t' << core.info << "\tGT:PS\t" << gt << ':' << ps << '\n';
    }
}

/**
 * @brief Writes optional candidate-variant VCF output.
 *
 * If `opts.output_vcf` is non-empty, opens the path, writes the VCF header and all records;
 * otherwise does nothing.
 *
 * @param opts Program options (`output_vcf`, BAM path, reference FASTA).
 * @param fai FASTA index for reference bases.
 * @param variants Categorized candidates to serialize.
 * @throws std::runtime_error If the BAM header cannot be read or the VCF path cannot be opened.
 */
void write_variants_vcf(const Options& opts, faidx_t* fai, const CandidateTable& variants) {
    if (opts.output_vcf.empty()) return;

    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai);

    std::ofstream out(opts.output_vcf);
    if (!out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);

    write_variants_vcf_header(out, opts, header.get());
    write_variants_vcf_records(out, opts, header.get(), ref, variants);
}

} // namespace pgphase_collect
