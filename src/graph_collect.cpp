#include "graph_collect.hpp"

#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_pipeline.hpp"
#include "collect_types.hpp"
#include "collect_var.hpp"
#include "gbz_ffi.h"
#include "graph_bam_adapter.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#include <htslib/faidx.h>
#include <htslib/sam.h>

namespace pgphase_collect {

namespace {

// Mirrors graph_bam_adapter.cpp site_key (static there; duplicated here to avoid coupling).
static std::string graph_site_key_str(const GraphSite& site, size_t /*index*/) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
}

static bool graph_site_starts_in_interval(const GraphSite& site,
                                          const std::string& contig,
                                          hts_pos_t beg,
                                          hts_pos_t end) {
    const std::string site_contig = site.ref_contig.empty() ? site.chrom : site.ref_contig;
    if (!contig.empty() && site_contig != contig) return false;
    const hts_pos_t site_beg0 = (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1;
    return site_beg0 >= beg && site_beg0 < end;
}

static GraphSiteCatalog filter_graph_catalog_for_interval(const GraphSiteCatalog& catalog,
                                                          const std::string& contig,
                                                          hts_pos_t beg,
                                                          hts_pos_t end) {
    GraphSiteCatalog out;
    for (const GraphSite& site : catalog.sites) {
        if (graph_site_starts_in_interval(site, contig, beg, end)) out.sites.push_back(site);
    }
    return out;
}

static bam_hdr_t* build_synthetic_header(faidx_t* fai) {
    bam_hdr_t* hdr = sam_hdr_init();
    if (!hdr) throw std::runtime_error("failed to allocate synthetic BAM header");
    const int n_seq = faidx_nseq(fai);
    for (int i = 0; i < n_seq; ++i) {
        const char* name = faidx_iseq(fai, i);
        const hts_pos_t len = faidx_seq_len(fai, name);
        const std::string len_str = std::to_string(len);
        if (sam_hdr_add_line(hdr, "SQ", "SN", name, "LN", len_str.c_str(), nullptr) < 0) {
            sam_hdr_destroy(hdr);
            throw std::runtime_error(std::string("failed to add SQ line for ") + name);
        }
    }
    return hdr;
}

// Converts phased multi-allelic graph candidates to biallelic CandidateVariants compatible
// with the existing TSV/VCF writers. One output row per passing alt allele per site.
static CandidateTable graph_chunks_to_candidate_table(
    const std::vector<GraphBamChunkBuildResult>& graph_chunks,
    const std::unordered_map<std::string, const GraphSite*>& site_id_to_site,
    const std::unordered_map<std::string, int>& contig_to_tid,
    const Options& opts)
{
    CandidateTable result;

    for (const GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        const BamChunk& chunk = graph_chunk.chunk;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& mcand = chunk.candidates[ci];

            auto site_it = site_id_to_site.find(graph_chunk.site_ids[ci]);
            if (site_it == site_id_to_site.end()) {
                // Biallelic pair from multiallelic site has ID "base_id:orig_alt"; strip suffix.
                const auto colon_pos = graph_chunk.site_ids[ci].rfind(':');
                if (colon_pos == std::string::npos) continue;
                site_it = site_id_to_site.find(graph_chunk.site_ids[ci].substr(0, colon_pos));
                if (site_it == site_id_to_site.end()) continue;
            }
            const GraphSite& site = *site_it->second;

            const std::string& chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
            auto tid_it = contig_to_tid.find(chrom);
            if (tid_it == contig_to_tid.end()) continue;
            const int fai_tid = tid_it->second;

            const std::vector<int>* orig_idx =
                (ci < graph_chunk.site_allele_orig_idx.size())
                    ? &graph_chunk.site_allele_orig_idx[ci]
                    : nullptr;

            const std::vector<int>& alle_covs = mcand.counts.alle_covs;
            const int n_new_alleles = static_cast<int>(alle_covs.size());

            for (int new_a = 1; new_a < n_new_alleles; ++new_a) {
                // Map surviving allele index back to original walk index in GraphSite.
                int orig_walk_idx = new_a;
                if (orig_idx != nullptr && new_a < static_cast<int>(orig_idx->size())) {
                    orig_walk_idx = (*orig_idx)[new_a];
                }

                // orig_walk_idx: 0 = ref walk, 1 = first alt walk → site.alts[0], etc.
                const int alt_idx = orig_walk_idx - 1;
                if (alt_idx < 0 || alt_idx >= static_cast<int>(site.alts.size())) continue;
                const std::string& alt_seq = site.alts[static_cast<size_t>(alt_idx)];
                if (alt_seq.empty() || alt_seq == "*") continue;

                const std::string& ref_seq = site.ref;
                if (ref_seq.empty()) continue;

                CandidateVariant cand;
                cand.key.tid = fai_tid;

                // Derive VariantKey from VCF-anchored alleles matching BAM-path convention:
                //   SNP : key.pos = site.pos, key.alt = alt base, ref_len = 1
                //   INS : key.pos = site.pos+1 (after anchor), key.alt = inserted bases
                //   DEL : key.pos = site.pos+1 (first deleted base), key.alt = deleted bases
                if (ref_seq.size() == 1 && alt_seq.size() == 1) {
                    // SNP
                    cand.key.type = VariantType::Snp;
                    cand.key.pos = site.pos;
                    cand.key.alt = alt_seq;
                    cand.key.ref_len = 1;
                } else if (alt_seq.size() > ref_seq.size() && ref_seq[0] == alt_seq[0]) {
                    // Left-anchored insertion: strip shared prefix; pos after anchor.
                    cand.key.type = VariantType::Insertion;
                    cand.key.pos = site.pos + 1;
                    cand.key.alt = alt_seq.substr(ref_seq.size());
                    cand.key.ref_len = 0;
                } else if (ref_seq.size() > alt_seq.size() && ref_seq[0] == alt_seq[0]) {
                    // Left-anchored deletion: empty alt matches BAM-path convention.
                    cand.key.type = VariantType::Deletion;
                    cand.key.pos = site.pos + static_cast<hts_pos_t>(alt_seq.size());
                    cand.key.alt = "";
                    cand.key.ref_len = static_cast<int>(ref_seq.size() - alt_seq.size());
                } else {
                    // Complex / MNP: no shared anchor base — classify by net length change.
                    cand.key.pos = site.pos;
                    cand.key.alt = alt_seq;
                    cand.key.ref_len = static_cast<int>(ref_seq.size());
                    if (alt_seq.size() > ref_seq.size()) {
                        cand.key.type = VariantType::Insertion;
                    } else if (ref_seq.size() > alt_seq.size()) {
                        cand.key.type = VariantType::Deletion;
                    } else {
                        cand.key.type = VariantType::Snp;  // MNP (equal length)
                    }
                }

                const int ref_cov = alle_covs.empty() ? 0 : alle_covs[0];
                const int alt_cov = alle_covs[static_cast<size_t>(new_a)];
                const int total_cov = ref_cov + alt_cov;
                cand.counts.ref_cov = ref_cov;
                cand.counts.alt_cov = alt_cov;
                cand.counts.total_cov = total_cov;
                cand.counts.forward_ref = ref_cov;
                cand.counts.forward_alt = alt_cov;
                cand.counts.allele_fraction =
                    total_cov > 0 ? static_cast<double>(alt_cov) / total_cov : 0.0;
                cand.counts.n_uniq_alles = 2;
                cand.counts.alle_covs = {ref_cov, alt_cov};

                const bool is_hom_alt = (ref_cov == 0 && alt_cov >= opts.min_alt_depth);
                if (alt_cov < opts.min_alt_depth || total_cov < opts.min_depth) {
                    cand.counts.category = VariantCategory::LowCoverage;
                    cand.counts.candvarcate_initial = VariantCategory::LowCoverage;
                    cand.lcd_var_i_to_cate = kLongcalldLowCovVar;
                } else if (is_hom_alt) {
                    cand.counts.category = VariantCategory::CleanHom;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHom;
                    cand.lcd_var_i_to_cate = kCandCleanHom;
                } else if (cand.counts.allele_fraction < opts.min_af ||
                           cand.counts.allele_fraction > opts.max_af) {
                    cand.counts.category = VariantCategory::LowAlleleFraction;
                    cand.counts.candvarcate_initial = VariantCategory::LowAlleleFraction;
                    cand.lcd_var_i_to_cate = kLongcalldLowAfVar;
                } else if (cand.key.type == VariantType::Snp) {
                    cand.counts.category = VariantCategory::CleanHetSnp;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHetSnp;
                    cand.lcd_var_i_to_cate = kCandCleanHetSnp;
                } else {
                    cand.counts.category = VariantCategory::CleanHetIndel;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
                    cand.lcd_var_i_to_cate = kCandCleanHetIndel;
                }

                // Translate multi-allelic hap_to_cons_alle to biallelic for this alt.
                // Homozygous alt: both haplotypes carry the alt allele, no phase set.
                const int hap1 = mcand.hap_to_cons_alle[1];
                const int hap2 = mcand.hap_to_cons_alle[2];
                cand.hap_to_cons_alle[0] = -1;
                if (is_hom_alt) {
                    cand.hap_to_cons_alle[1] = 1;
                    cand.hap_to_cons_alle[2] = 1;
                    cand.hap_alt = 1;
                    cand.hap_ref = 1;
                    cand.phase_set = -1;
                } else {
                    cand.hap_to_cons_alle[1] = (hap1 == new_a) ? 1 : 0;
                    cand.hap_to_cons_alle[2] = (hap2 == new_a) ? 1 : 0;
                    cand.hap_alt = cand.hap_to_cons_alle[1];
                    cand.hap_ref = cand.hap_to_cons_alle[2];
                    cand.phase_set = mcand.phase_set;
                }
                cand.alt_ref_base = 4;  // use FASTA anchor (BAM-path default)
                cand.lcd_make_variants_region_pass = true;

                result.push_back(std::move(cand));
            }
        }
    }

    std::stable_sort(result.begin(), result.end(),
                     [](const CandidateVariant& a, const CandidateVariant& b) {
                         return exact_comp_cand_var(&a, &b) < 0;
                     });

    return result;
}

// Processes one batch of graph chunks in parallel (one thread pool per reg_chunk_i batch,
// mirroring collect_chunk_batch_parallel in collect_pipeline.cpp).
// Each worker queries overlapping reads via the gbz-base FFI (one SQLite connection
// per thread), then build_graph_bam_chunk + assign_hap k-means.
// Peak memory = threads × reads_per_chunk (like BAM).
// After all workers join: populate_graph_chunk_overlaps + stitch_chunk_haps.
static std::vector<GraphBamChunkBuildResult> process_graph_chunk_batch(
    const GraphSiteCatalog& catalog,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const GraphQueryConfig& qconfig,
    const std::string& ref_sample,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const Options& opts)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphBamChunkBuildResult> graph_chunks(batch_size);

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            // Each thread opens its own GBZ + GAF database connection
            // (SQLite requires one connection per thread).
            char* err = nullptr;
            void* gbz_h = pgphase_gbz_open(qconfig.gbz_db.c_str(), &err);
            if (!gbz_h) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("failed to open GBZ: " + msg));
                return;
            }
            void* gaf_h = pgphase_gaf_open(qconfig.gaf_db.c_str(), &err);
            if (!gaf_h) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                pgphase_gbz_close(gbz_h);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("failed to open GAF-base: " + msg));
                return;
            }
            // Validate GBZ/GAF compatibility (same check the query binary does).
            if (pgphase_gbz_gaf_validate(gbz_h, gaf_h, &err) != 0) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                pgphase_gaf_close(gaf_h);
                pgphase_gbz_close(gbz_h);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("GBZ/GAF-base incompatible: " + msg));
                return;
            }
            // RAII cleanup for FFI handles.
            struct HandleCleanup {
                void* gbz; void* gaf;
                ~HandleCleanup() {
                    if (gaf) pgphase_gaf_close(gaf);
                    if (gbz) pgphase_gbz_close(gbz);
                }
            } cleanup{gbz_h, gaf_h};

            try {
                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    const RegionChunk& region = chunks[batch_begin + offset];
                    const std::string contig_name = header->target_name[region.tid];

                    // Strip pangenome sample prefix:
                    // "CHM13#0#chr20" → "chr20"; plain "chr20" stays unchanged.
                    const std::string query_contig = [&]() -> std::string {
                        auto it = fai_full_to_suffix.find(contig_name);
                        return (it != fai_full_to_suffix.end()) ? it->second : contig_name;
                    }();
                    GraphSiteCatalog chunk_catalog = filter_graph_catalog_for_interval(
                        catalog, contig_name, region.beg - 1, region.end);

                    std::vector<GraphReadAllele> chunk_rows;
                    if (!chunk_catalog.sites.empty()) {
                        chunk_rows = query_gbz_interval_gaf_ffi(
                            gbz_h, gaf_h, ref_sample, query_contig,
                            region.beg - 1, region.end,
                            chunk_catalog, qconfig.min_mapq);
                    }

                    graph_chunks[offset] = build_graph_bam_chunk(
                        chunk_catalog,
                        chunk_rows,
                        contig_name,
                        region.beg - 1,
                        region.end,
                        region.chunk_id,
                        opts);

                    assign_hap_based_on_germline_het_vars_kmeans(
                        graph_chunks[offset].chunk, opts, kCandGermlineClean);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);

    // populate_graph_chunk_overlaps finds reads shared between adjacent chunks by name,
    // setting down_ovlp_read_i / up_ovlp_read_i for stitch_chunk_haps.
    populate_graph_chunk_overlaps(graph_chunks);

    // stitch_chunk_haps needs a flat BamChunk vector (mirrors phase_graph_bam_chunks).
    std::vector<BamChunk> bam_chunks;
    bam_chunks.reserve(batch_size);
    for (GraphBamChunkBuildResult& gc : graph_chunks)
        bam_chunks.push_back(std::move(gc.chunk));
    stitch_chunk_haps(bam_chunks, &opts, nullptr);
    for (size_t i = 0; i < batch_size; ++i)
        graph_chunks[i].chunk = std::move(bam_chunks[i]);

    return graph_chunks;
}

// Per-chunk tabix queries on an indexed GAF file.  Each worker seeks directly
// to the overlapping region — only the relevant reads are decompressed and
// parsed, making this efficient even for very large (100+ GB) GAF files.
static std::vector<GraphBamChunkBuildResult> process_graph_chunk_batch_indexed_gaf(
    const GraphSiteCatalog& catalog,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const std::string& gaf_file,
    int min_mapq,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const Options& opts)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphBamChunkBuildResult> graph_chunks(batch_size);

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            try {
                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    const RegionChunk& region = chunks[batch_begin + offset];
                    const std::string contig_name = header->target_name[region.tid];

                    // The tabix index uses reference-suffix names ("chr20"),
                    // while the header may use pangenome names ("CHM13#0#chr20").
                    // scan_indexed_gaf_chunk handles this via
                    // tbx_seq_tid_with_pangenome_fallback, but we also try
                    // the suffix directly for a faster first lookup.
                    const std::string query_contig = [&]() -> std::string {
                        auto it = fai_full_to_suffix.find(contig_name);
                        return (it != fai_full_to_suffix.end()) ? it->second : contig_name;
                    }();

                    GraphSiteCatalog chunk_catalog = filter_graph_catalog_for_interval(
                        catalog, contig_name, region.beg - 1, region.end);

                    std::vector<GraphReadAllele> chunk_rows;
                    if (!chunk_catalog.sites.empty()) {
                        chunk_rows = scan_indexed_gaf_chunk(
                            gaf_file, query_contig,
                            region.beg - 1, region.end,
                            chunk_catalog, min_mapq);
                    }

                    graph_chunks[offset] = build_graph_bam_chunk(
                        chunk_catalog,
                        chunk_rows,
                        contig_name,
                        region.beg - 1,
                        region.end,
                        region.chunk_id,
                        opts);

                    assign_hap_based_on_germline_het_vars_kmeans(
                        graph_chunks[offset].chunk, opts, kCandGermlineClean);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);

    populate_graph_chunk_overlaps(graph_chunks);

    std::vector<BamChunk> bam_chunks;
    bam_chunks.reserve(batch_size);
    for (GraphBamChunkBuildResult& gc : graph_chunks)
        bam_chunks.push_back(std::move(gc.chunk));
    stitch_chunk_haps(bam_chunks, &opts, nullptr);
    for (size_t i = 0; i < batch_size; ++i)
        graph_chunks[i].chunk = std::move(bam_chunks[i]);

    return graph_chunks;
}

void run_collect_graph_variation(const Options& opts) {
    const bool use_indexed_gaf = !opts.gaf_file.empty();
    if (!use_indexed_gaf) {
        if (opts.gbz_db.empty())
            throw std::runtime_error("--gbz-db is required for collect-graph-variation without --gaf");
        if (opts.gaf_db.empty())
            throw std::runtime_error("--gaf-db is required for collect-graph-variation without --gaf");
    } else {
        // The --gaf path requires a tabix-indexed, bgzip-compressed GAF with
        // annotated coordinate columns so per-chunk region queries are efficient.
        require_indexed_gaf(opts.gaf_file);
        if (!opts.gbz_db.empty() || !opts.gaf_db.empty()) {
            std::cerr << "Warning: --gaf provided; ignoring --gbz-db/--gaf-db\n";
        }
    }

    // 1. Reference FASTA index first — needed to resolve autosome contig names for
    //    region filters, which are then passed to the VCF loader so only sites in the
    //    requested regions are parsed (tabix-assisted for bgzipped + indexed VCFs).
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(build_synthetic_header(fai.get()));

    // 2. Build a chrom alias map between the FASTA and VCF naming conventions.
    //    Pangenome FASTAs use "SAMPLE#HAP#CHROM" (e.g. "CHM13#0#chr20") while graph VCFs
    //    use the plain reference name ("chr20"), or vice versa.  We inspect the FAI names
    //    and build a bidirectional suffix map so both directions resolve automatically.
    //
    //    fai_suffix_to_full : "chr20" → "CHM13#0#chr20"  (used when VCF is short, FAI is full)
    //    fai_full_to_suffix : "CHM13#0#chr20" → "chr20"  (used when VCF is full, FAI is short)
    std::unordered_map<std::string, std::string> fai_suffix_to_full;
    std::unordered_map<std::string, std::string> fai_full_to_suffix;
    {
        const int nseq = faidx_nseq(fai.get());
        fai_suffix_to_full.reserve(static_cast<size_t>(nseq));
        for (int i = 0; i < nseq; ++i) {
            const std::string full(faidx_iseq(fai.get(), i));
            const size_t h = full.rfind('#');
            if (h != std::string::npos) {
                fai_suffix_to_full.emplace(full.substr(h + 1), full);
                fai_full_to_suffix.emplace(full, full.substr(h + 1));
            }
        }
    }

    // Resolves a contig name to its canonical FAI form (the name that exists in the header).
    // Handles both "chr20" → "CHM13#0#chr20" and "CHM13#0#chr20" → "chr20".
    auto resolve_contig = [&](const std::string& name) -> std::string {
        if (faidx_has_seq(fai.get(), name.c_str())) return name;
        // Short name → full pangenome name ("chr20" → "CHM13#0#chr20")
        auto it = fai_suffix_to_full.find(name);
        if (it != fai_suffix_to_full.end()) return it->second;
        // Full pangenome name → short name ("CHM13#0#chr20" → "chr20")
        auto it2 = fai_full_to_suffix.find(name);
        if (it2 != fai_full_to_suffix.end() && faidx_has_seq(fai.get(), it2->second.c_str()))
            return it2->second;
        return name; // unchanged; add_filter_chunks will throw a clear error
    };

    std::vector<RegionFilter> region_filters;
    for (const std::string& r : opts.regions) {
        RegionFilter f = parse_region(r);
        f.chrom = resolve_contig(f.chrom);
        region_filters.push_back(std::move(f));
    }
    if (!opts.region_file.empty()) {
        auto bed = load_bed_regions(opts.region_file);
        for (RegionFilter& f : bed) f.chrom = resolve_contig(f.chrom);
        region_filters.insert(region_filters.end(), bed.begin(), bed.end());
    }
    if (opts.autosome) {
        // With pangenome FASTAs, "chr1"–"chr22" won't be found directly; resolve_contig
        // maps them to the full name (e.g. "CHM13#0#chr1").
        for (int i = 1; i <= 22; ++i) {
            for (const std::string& candidate :
                     {"chr" + std::to_string(i), std::to_string(i)}) {
                const std::string resolved = resolve_contig(candidate);
                if (faidx_has_seq(fai.get(), resolved.c_str())) {
                    region_filters.push_back(RegionFilter{true, resolved, 1, -1});
                    break;
                }
            }
        }
    }

    // When no region filters were specified, build whole-chromosome filters for every
    // FASTA contig so the VCF loader can use tabix instead of streaming the entire file.
    if (region_filters.empty()) {
        const int nseq = faidx_nseq(fai.get());
        region_filters.reserve(static_cast<size_t>(nseq));
        for (int i = 0; i < nseq; ++i)
            region_filters.push_back(
                RegionFilter{true, std::string(faidx_iseq(fai.get(), i)), 1, -1});
    }

    // 3. Load graph site catalog — region-filtered if any filters are active.
    //    For bgzipped + tabix-indexed VCFs this seeks directly to each region.
    // Both paths (indexed GAF and GBZ-DB FFI) need allele traversal strings
    // for compact site index construction / structured callback matching.
    GraphSiteCatalog catalog = load_graph_site_catalog_from_vcf(opts.graph_sites_vcf,
                                                                 region_filters,
                                                                 true);
    if (opts.verbose >= 1) {
        std::cerr << "Loaded " << catalog.sites.size() << " graph sites from "
                  << opts.graph_sites_vcf << "\n";
    }

    // 4a. Normalize site chrom/ref_contig to match the FAI naming convention.
    //     After this, all downstream comparisons (site_starts_in_interval, contig_to_tid,
    //     build_site_to_chunk_map) use exact string matching with no special cases.
    //     The remap covers both directions detected above.
    {
        // Build the VCF-name → FAI-name map from the two suffix tables.
        std::unordered_map<std::string, std::string> chrom_remap;
        // VCF uses short name, FAI uses full ("chr20" → "CHM13#0#chr20")
        for (const auto& kv : fai_suffix_to_full) chrom_remap.emplace(kv.first, kv.second);
        // VCF uses full name, FAI uses short ("CHM13#0#chr20" → "chr20")
        for (const auto& kv : fai_full_to_suffix) {
            if (faidx_has_seq(fai.get(), kv.second.c_str()))
                chrom_remap.emplace(kv.first, kv.second);
        }
        for (GraphSite& site : catalog.sites) {
            auto it = chrom_remap.find(site.chrom);
            if (it != chrom_remap.end()) site.chrom = it->second;
            if (!site.ref_contig.empty()) {
                auto it2 = chrom_remap.find(site.ref_contig);
                if (it2 != chrom_remap.end()) site.ref_contig = it2->second;
            }
        }
    }

    // 5. contig name → FAI-order tid (used when remapping candidate tids for output).
    std::unordered_map<std::string, int> contig_to_tid;
    contig_to_tid.reserve(static_cast<size_t>(header->n_targets));
    for (int32_t tid = 0; tid < header->n_targets; ++tid) {
        contig_to_tid[header->target_name[tid]] = tid;
    }

    // 6. Tile genome into RegionChunks using the resolved region_filters (which may have
    //    contig names like "CHM13#0#chr20" resolved from a user-supplied "chr20").
    const std::vector<RegionChunk> chunks =
        build_region_chunks(opts, header.get(), fai.get(), region_filters);
    if (chunks.empty()) {
        std::cerr << "No region chunks to process\n";
        return;
    }
    if (opts.verbose >= 1) {
        std::cerr << "Tiled genome into " << chunks.size() << " chunks ("
                  << opts.chunk_size << " bp, " << opts.threads << " thread(s))\n";
    }

    // 7. Build per-chunk query config for the GBZ/GAF-base FFI path.
    GraphQueryConfig qconfig;
    qconfig.gbz_db   = opts.gbz_db;
    qconfig.gaf_db   = opts.gaf_db;
    qconfig.min_mapq  = opts.min_mapq;

    // Derive the reference sample name for GBZ interval queries.
    // For pangenome FASTAs ("CHM13#0#chr20") extract "CHM13".
    // For plain FASTAs ("chr20") leave empty so query uses its GENERIC_SAMPLE default.
    std::string ref_sample = opts.graph_sample;
    if (ref_sample.empty()) {
        const int nseq = faidx_nseq(fai.get());
        for (int i = 0; i < nseq && ref_sample.empty(); ++i) {
            const std::string name(faidx_iseq(fai.get(), i));
            const size_t h = name.find('#');
            if (h != std::string::npos && h > 0)
                ref_sample = name.substr(0, h);
        }
    }
    if (!use_indexed_gaf && opts.verbose >= 1 && !ref_sample.empty())
        std::cerr << "Using reference sample \"" << ref_sample << "\" for GBZ interval queries\n";

    // 8. site_id → GraphSite* (for biallelic conversion after phasing).
    std::unordered_map<std::string, const GraphSite*> site_id_to_site;
    site_id_to_site.reserve(catalog.sites.size());
    for (size_t i = 0; i < catalog.sites.size(); ++i) {
        site_id_to_site[graph_site_key_str(catalog.sites[i], i)] = &catalog.sites[i];
    }

    // 9. Open output streams.
    std::ofstream variant_out(opts.output_tsv);
    if (!variant_out) throw std::runtime_error("failed to open output: " + opts.output_tsv);
    write_variants_tsv_header(variant_out);

    std::ofstream vcf_out;
    if (!opts.output_vcf.empty()) {
        vcf_out.open(opts.output_vcf);
        if (!vcf_out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);
        write_variants_vcf_header(vcf_out, opts, header.get());
    }
    std::ofstream phased_vcf_out;
    if (!opts.output_phased_vcf.empty()) {
        phased_vcf_out.open(opts.output_phased_vcf);
        if (!phased_vcf_out)
            throw std::runtime_error("failed to open phased VCF output: " + opts.output_phased_vcf);
        write_phased_variants_vcf_header(phased_vcf_out, opts, header.get());
    }

    ReferenceCache ref(fai.get());

    // Phased BAM: open output file and write header.
    struct SamFileCloser { void operator()(samFile* fp) const { if (fp) hts_close(fp); } };
    std::unique_ptr<samFile, SamFileCloser> phased_bam_fp;
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> phased_bam_hdr(nullptr, &sam_hdr_destroy);
    std::unordered_map<std::string, PhaseReadOutputRow> phased_bam_rows;
    std::unordered_set<std::string> phased_bam_emitted;
    const bool emit_phased_bam = !opts.output_phased_bam.empty();
    if (emit_phased_bam) {
        samFile* raw_fp = hts_open(opts.output_phased_bam.c_str(), "wb");
        if (!raw_fp)
            throw std::runtime_error("failed to open phased BAM: " + opts.output_phased_bam);
        phased_bam_fp.reset(raw_fp);
        sam_hdr_t* hdr = sam_hdr_init();
        if (!hdr) throw std::runtime_error("failed to allocate phased BAM header");
        phased_bam_hdr.reset(hdr);
        if (sam_hdr_write(phased_bam_fp.get(), phased_bam_hdr.get()) < 0)
            throw std::runtime_error("failed to write phased BAM header");
    }

    // 10. Process chunks in reg_chunk_i batches (one contig per batch), streaming output.
    //     Mirrors run_collect_bam_variation's batch loop exactly.
    size_t n_variants = 0;
    size_t batch_begin = 0;
    while (batch_begin < chunks.size()) {
        size_t batch_end = batch_begin + 1;
        while (batch_end < chunks.size() &&
               chunks[batch_end].reg_chunk_i == chunks[batch_begin].reg_chunk_i) {
            ++batch_end;
        }

        std::vector<GraphBamChunkBuildResult> graph_chunks =
            use_indexed_gaf
                ? process_graph_chunk_batch_indexed_gaf(
                      catalog, chunks, batch_begin, batch_end, header.get(),
                      opts.gaf_file, opts.min_mapq, fai_full_to_suffix, opts)
                : process_graph_chunk_batch(
                      catalog, chunks, batch_begin, batch_end, header.get(),
                      qconfig, ref_sample, fai_full_to_suffix, opts);

        CandidateTable variants =
            graph_chunks_to_candidate_table(graph_chunks, site_id_to_site, contig_to_tid, opts);

        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty())
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        if (!opts.output_phased_vcf.empty())
            write_phased_variants_vcf_records(phased_vcf_out, opts, header.get(), ref, variants);

        // Phased BAM: accumulate post-stitch read assignments, flush reads
        // that won't appear in the next batch.
        if (emit_phased_bam) {
            for (const GraphBamChunkBuildResult& gc : graph_chunks)
                merge_graph_chunk_into_read_rows(phased_bam_rows, gc);
            // Collect read names in the next batch to avoid flushing overlap reads early.
            std::unordered_set<std::string> next_batch_qnames;
            if (batch_end < chunks.size()) {
                // Peek: next batch shares the same reg_chunk_i? If not, no overlap.
                // Cross-contig batches never share reads, so flush all.
            }
            flush_graph_phase_bam_after_merge(
                phased_bam_fp.get(), phased_bam_hdr.get(),
                phased_bam_rows, nullptr, phased_bam_emitted);
        }

        batch_begin = batch_end;
    }

    std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
              << " worker thread(s)\n";
    std::cerr << "Collected " << n_variants << " candidate variant sites into "
              << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty())
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    if (!opts.output_phased_vcf.empty())
        std::cerr << "Wrote phased candidate VCF to " << opts.output_phased_vcf << "\n";
    if (emit_phased_bam)
        std::cerr << "Wrote phased BAM to " << opts.output_phased_bam << "\n";
}

static void print_graph_collect_help() {
    std::cout
        << "Usage: pgphase collect-graph-variation [options]\n"
        << "\n"
        << "Required:\n"
        << "      --ref FILE                Reference FASTA (indexed)\n"
        << "      --sites FILE              Sites VCF from build-snarl-catalog (bgzipped + tabix-indexed)\n"
        << "\n"
        << "Read input:\n"
        << "      --gaf FILE                Raw GAF alignments; builds a fast pgphase observation index\n"
        << "      --gbz-db FILE             GBZ graph database (legacy --gaf-db path)\n"
        << "      --gaf-db FILE             GAF-base read alignment database (legacy path)\n"
        << "\n"
        << "Options:\n"
        << "  -o, --output FILE             Output TSV [output.tsv]\n"
        << "  -v, --vcf-output FILE         Candidate VCF output\n"
        << "      --phased-vcf-out FILE     Phased VCF with GT:DP:AD:VAF:GQ:PS\n"
        << "      --phased-bam-out FILE     Unaligned BAM with HP/PS tags per read\n"
        << "  -t, --threads INT             Worker threads [1]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
        << "  -D, --min-depth INT           Minimum total depth [5]\n"
        << "      --min-alt-depth INT       Minimum alt depth [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction [0.80]\n"
        << "      --min-sv-len INT          Min SV length for SVTYPE/SVLEN tags [30]\n"
        << "      --chunk-size INT          Region chunk size in bp [500000]\n"
        << "  -r, --region STR              Restrict to region (may be repeated)\n"
        << "      --region-file FILE        BED file of regions\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "      --sample NAME             Reference sample name for GBZ interval queries\n"
        << "                                (auto-derived from FASTA if not provided)\n"
        << "      --gbz-query-bin FILE      Path to GBZ-base query binary [query]\n"
        << "      --ont                     ONT read mode (enables strand-bias filter)\n"
        << "  -V, --verbose INT             Verbosity level [0]\n"
        << "  -h, --help                    Print this help\n"
        << "\n"
        << "Examples:\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites sites.vcf.gz \\\n"
        << "      --gaf reads.gaf \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      --phased-bam-out phased.bam \\\n"
        << "      -t 8\n"
        << "\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites sites.vcf.gz \\\n"
        << "      --gbz-db reads.gaf.db \\\n"
        << "      --ont \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      -t 16\n";
}

enum GraphCollectOption {
    kGcMinAltDepth = 1000,
    kGcMinAf,
    kGcMaxAf,
    kGcMinSvLen,
    kGcChunkSize,
    kGcPhasedVcf,
    kGcGbzDb,
    kGcGafFile,
    kGcGafDb,
    kGcRegionFile,
    kGcAutosome,
    kGcSample,
    kGcGbzQueryBin,
    kGcOnt,
    kGcPhasedBam,
    kGcRef,
    kGcSites,
};

} // namespace

} // namespace pgphase_collect

int collect_graph_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;

    {
        std::ostringstream cmd;
        cmd << "pgphase collect-graph-variation";
        for (int i = 1; i < argc; ++i) cmd << ' ' << argv[i];
        opts.command_line = cmd.str();
    }

    optind = 1;
    const struct option long_options[] = {
        {"output",            required_argument, nullptr, 'o'},
        {"vcf-output",        required_argument, nullptr, 'v'},
        {"phased-vcf-out",    required_argument, nullptr, kGcPhasedVcf},
        {"phased-bam-out",   required_argument, nullptr, kGcPhasedBam},
        {"threads",           required_argument, nullptr, 't'},
        {"min-mapq",          required_argument, nullptr, 'q'},
        {"min-depth",         required_argument, nullptr, 'D'},
        {"min-alt-depth",     required_argument, nullptr, kGcMinAltDepth},
        {"min-af",            required_argument, nullptr, kGcMinAf},
        {"max-af",            required_argument, nullptr, kGcMaxAf},
        {"min-sv-len",        required_argument, nullptr, kGcMinSvLen},
        {"chunk-size",        required_argument, nullptr, kGcChunkSize},
        {"region",            required_argument, nullptr, 'r'},
        {"region-file",       required_argument, nullptr, kGcRegionFile},
        {"autosome",          no_argument,       nullptr, kGcAutosome},
        {"gaf",               required_argument, nullptr, kGcGafFile},
        {"gaf-file",          required_argument, nullptr, kGcGafFile},
        {"gbz-db",            required_argument, nullptr, kGcGbzDb},
        {"gaf-db",            required_argument, nullptr, kGcGafDb},
        {"sample",            required_argument, nullptr, kGcSample},
        {"gbz-query-bin",     required_argument, nullptr, kGcGbzQueryBin},
        {"ont",               no_argument,       nullptr, kGcOnt},
        {"ref",               required_argument, nullptr, kGcRef},
        {"sites",             required_argument, nullptr, kGcSites},
        {"verbose",           required_argument, nullptr, 'V'},
        {"help",              no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "o:v:t:q:D:r:V:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'o': opts.output_tsv = optarg; break;
            case 'v': opts.output_vcf = optarg; break;
            case kGcPhasedVcf:    opts.output_phased_vcf = optarg; break;
            case kGcPhasedBam:    opts.output_phased_bam = optarg; break;
            case 't': opts.threads = std::stoi(optarg); break;
            case 'q': opts.min_mapq = std::stoi(optarg); break;
            case 'D': opts.min_depth = std::stoi(optarg); break;
            case kGcMinAltDepth:  opts.min_alt_depth = std::stoi(optarg); break;
            case kGcMinAf:        opts.min_af = std::stod(optarg); break;
            case kGcMaxAf:        opts.max_af = std::stod(optarg); break;
            case kGcMinSvLen:     opts.min_sv_len = std::stoi(optarg); break;
            case kGcChunkSize:    opts.chunk_size = std::stoll(optarg); break;
            case 'r': opts.regions.push_back(optarg); break;
            case kGcRegionFile:   opts.region_file = optarg; break;
            case kGcAutosome:     opts.autosome = true; break;
            case kGcGafFile:      opts.gaf_file = optarg; break;
            case kGcGbzDb:        opts.gbz_db = optarg; break;
            case kGcGafDb:        opts.gaf_db = optarg; break;
            case kGcSample:       opts.graph_sample = optarg; break;
            case kGcGbzQueryBin:  opts.gbz_query_bin = optarg; break;
            case kGcOnt:          opts.read_technology = ReadTechnology::Ont; break;
            case kGcRef:          opts.ref_fasta = optarg; break;
            case kGcSites:        opts.graph_sites_vcf = optarg; break;
            case 'V': opts.verbose = std::stoi(optarg); break;
            case 'h': print_graph_collect_help(); return 0;
            default:  print_graph_collect_help(); return 1;
        }
    }

    if (opts.ref_fasta.empty() || opts.graph_sites_vcf.empty()) {
        std::cerr << "Error: --ref and --sites are required\n";
        print_graph_collect_help();
        return 1;
    }

    require_graph_site_vcf_tabix_index(opts.graph_sites_vcf);

    if (opts.gaf_file.empty() && (opts.gbz_db.empty() || opts.gaf_db.empty())) {
        std::cerr << "Error: provide --gaf, or provide --gbz-db and --gaf-db\n";
        print_graph_collect_help();
        return 1;
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_depth < 0 ||
        opts.min_alt_depth < 0 || opts.min_af < 0.0 || opts.max_af < opts.min_af ||
        opts.min_sv_len < 0 || opts.chunk_size < 1 || opts.verbose < 0) {
        std::cerr << "Error: invalid numeric threshold\n";
        return 1;
    }

    try {
        run_collect_graph_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
