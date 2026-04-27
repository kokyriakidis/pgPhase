/**
 * @file collect_output.cpp
 * @brief TSV, optional VCF, and read-support writers for collect-bam-variation candidates.
 *
 * @details Outputs describe pre-phasing candidates (not diploid genotypes). Opening the primary BAM
 * for SQ names may throw `std::runtime_error` on I/O failure.
 */

#include "collect_output.hpp"

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
}

/**
 * @brief Writes VCF body records for candidate variants.
 *
 * Produces left-normalized records: SNP at `key.pos`; insertions and deletions use the
 * anchor base at POS-1 per VCF convention. FILTER is PASS for clean categories, RefCall
 * for non-variants, NoCall when depth is zero, and LowQual otherwise. INFO includes END,
 * optional CLEAN, depth and allele fields, AF, CAT, and SVTYPE/SVLEN when |SVLEN| >=
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
    for (const CandidateVariant& candidate : variants) {
        const VariantKey& key = candidate.key;
        const VariantCounts& counts = candidate.counts;
        const std::string chrom = header->target_name[key.tid];

        hts_pos_t pos = key.pos;
        std::string ref_seq;
        std::string alt_seq;

        if (key.type == VariantType::Snp) {
            ref_seq = std::string(1, ref.base(key.tid, key.pos, header));
            alt_seq = key.alt.empty() ? "." : key.alt;
        } else if (key.type == VariantType::Insertion) {
            const hts_pos_t anchor_pos = std::max<hts_pos_t>(1, key.pos - 1);
            pos = anchor_pos;
            const char anchor_base = ref.base(key.tid, anchor_pos, header);
            ref_seq = std::string(1, anchor_base);
            alt_seq = ref_seq + key.alt;
        } else { // Deletion
            const hts_pos_t anchor_pos = std::max<hts_pos_t>(1, key.pos - 1);
            pos = anchor_pos;
            const char anchor_base = ref.base(key.tid, anchor_pos, header);
            const std::string del_seq = ref.subseq(key.tid, key.pos, key.ref_len, header);
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
