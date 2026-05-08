#include "graph_pipeline.hpp"

#include "graph_bam_adapter.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <algorithm>
#include <atomic>
#include <exception>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <optional>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kstring.h>

namespace {

struct GraphOptions {
    std::string gaf_file;
    std::string graph_sites_vcf;
    std::string graph_sites_tsv;
    std::string graph_read_support_tsv;
    std::string graph_site_counts_tsv;
    std::string graph_read_profile_tsv;
    std::string graph_phase_sites_tsv = "graph_phase_sites.tsv";
    std::string graph_phase_reads_tsv = "graph_phase_reads.tsv";
    std::string graph_phase_bam;
    std::string graph_filtered_sites_tsv;
    std::string graph_contig;
    hts_pos_t graph_interval_beg = -1;
    hts_pos_t graph_interval_end = -1;
    hts_pos_t graph_chunk_size = 500000;
    int threads = 1;
    int verbose = 0;
};

enum LongOption {
    kGafFileOption = 1000,
    kGraphSitesVcfOption,
    kGraphSitesTsvOption,
    kGraphReadSupportTsvOption,
    kGraphSiteCountsTsvOption,
    kGraphReadProfileTsvOption,
    kGraphPhaseSitesTsvOption,
    kGraphPhaseReadsTsvOption,
    kGraphFilteredSitesTsvOption,
    kGraphPhaseBamOption,
    kGraphContigOption,
    kGraphIntervalOption,
    kGraphChunkSizeOption,
    kThreadsOption
};

void print_help() {
    std::cout
        << "Usage: pgphase phase-graph [options]\n"
        << "\n"
        << "Required:\n"
        << "      --gaf FILE                pggaf indexed annotated GAF (.gaf.gz + .tbi; rc/rb/re required)\n"
        << "      --graph-sites-vcf FILE    Decomposed graph-site VCF (from build-snarl-catalog)\n"
        << "\n"
        << "Options:\n"
        << "      --contig NAME             Restrict to this reference contig\n"
        << "      --interval BEG..END       0-based half-open coordinate bounds (default: full catalog)\n"
        << "      --chunk-size INT          Phasing chunk size in bp [500000]\n"
        << "      --graph-phase-sites FILE  Output: graph-site hap consensus TSV [graph_phase_sites.tsv]\n"
        << "      --graph-phase-reads FILE  Output: read HAP / PHASE_SET TSV [graph_phase_reads.tsv]\n"
        << "      --graph-phase-bam FILE    Output: unaligned BAM with HP/PS tags per read\n"
        << "      --graph-sites-tsv FILE    Diagnostic: parsed graph-site catalog\n"
        << "      --graph-read-support FILE Diagnostic: read->site allele support TSV\n"
        << "      --graph-site-counts FILE  Diagnostic: site allele depth/count TSV\n"
        << "      --graph-read-profile FILE Diagnostic: sparse read x site allele profile TSV\n"
        << "      --graph-filtered-sites FILE Diagnostic: sites removed by depth/AF filters\n"
        << "  -t, --threads INT             Worker threads [1]\n"
        << "  -V, --verbose INT             Verbosity level [0]\n"
        << "  -h, --help                    Print this help\n";
}

bool parse_half_open_interval(const std::string& text, hts_pos_t& beg, hts_pos_t& end) {
    const size_t sep = text.find("..");
    if (sep == std::string::npos) return false;
    const std::string left = text.substr(0, sep);
    const std::string right = text.substr(sep + 2);
    if (left.empty() || right.empty()) return false;
    beg = static_cast<hts_pos_t>(std::stoll(left));
    end = static_cast<hts_pos_t>(std::stoll(right));
    return beg >= 0 && end > beg;
}

std::vector<std::pair<hts_pos_t, hts_pos_t>>
split_graph_interval(hts_pos_t beg, hts_pos_t end, hts_pos_t chunk_size) {
    std::vector<std::pair<hts_pos_t, hts_pos_t>> chunks;
    for (hts_pos_t chunk_beg = beg; chunk_beg < end; chunk_beg += chunk_size) {
        chunks.emplace_back(chunk_beg, std::min(end, chunk_beg + chunk_size));
    }
    return chunks;
}

void append_catalog(pgphase_collect::GraphSiteCatalog& dst,
                    const pgphase_collect::GraphSiteCatalog& src) {
    dst.sites.insert(dst.sites.end(), src.sites.begin(), src.sites.end());
}


struct GraphChunkWorkResult {
    pgphase_collect::GraphSiteCatalog catalog;
    std::vector<pgphase_collect::GraphReadAllele> rows;
    pgphase_collect::GraphBamChunkBuildResult build_result;
};

// VCF+indexed-GAF fast path: catalogs are pre-sliced; rows are retrieved per
// chunk with tabix from pggaf's rc/rb/re coordinate index.
// Processes the half-open index range [batch_beg, batch_end) in parallel, returning
// results[0..batch_size-1] corresponding to chunk indices [batch_beg..batch_end-1].
// chunk_id_offset is added to each local chunk index so IDs are globally monotonic
// across contigs.
// Peak RAM = n_workers × one_chunk_rows + batch_size × BamChunk.
static std::vector<GraphChunkWorkResult> run_vcf_gaf_chunks_parallel(
    const std::vector<pgphase_collect::GraphSiteCatalog>& chunk_catalogs,
    const std::string& indexed_gaf_file,
    const std::vector<std::pair<hts_pos_t, hts_pos_t>>& intervals,
    size_t batch_beg,
    size_t batch_end,
    int n_threads,
    const std::string& contig,
    const pgphase_collect::Options& filter_opts,
    int chunk_id_offset = 0) {
    using namespace pgphase_collect;
    const size_t batch_size = batch_end - batch_beg;
    std::vector<GraphChunkWorkResult> results(batch_size);
    const size_t n_workers = std::min<size_t>(static_cast<size_t>(n_threads), batch_size);
    std::atomic<size_t> next{0};
    std::exception_ptr first_error;
    std::mutex err_mutex;
    std::vector<std::thread> workers;
    workers.reserve(n_workers);
    for (size_t wi = 0; wi < n_workers; ++wi) {
        workers.emplace_back([&]() {
            try {
                while (true) {
                    const size_t offset = next.fetch_add(1);
                    if (offset >= batch_size) break;
                    const size_t i = batch_beg + offset;
                    std::vector<GraphReadAllele> chunk_rows = scan_indexed_gaf_chunk(
                        indexed_gaf_file, contig,
                        intervals[i].first, intervals[i].second,
                        chunk_catalogs[i], filter_opts.min_mapq);
                    results[offset].build_result = build_graph_bam_chunk(
                        chunk_catalogs[i], chunk_rows, contig,
                        intervals[i].first, intervals[i].second,
                        chunk_id_offset + static_cast<int>(i), filter_opts);
                    results[offset].catalog = chunk_catalogs[i];
                    results[offset].rows = std::move(chunk_rows);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(err_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);
    return results;
}

// GBZ path: catalog-fetch (via gbz-base query) + per-chunk GAF DB lookup + build.
// GetCatalog(chunk_i, beg, end) → GraphSiteCatalog.
template<typename GetCatalog>
static std::vector<GraphChunkWorkResult> run_gbz_chunks_parallel(
    const std::vector<std::pair<hts_pos_t, hts_pos_t>>& intervals,
    int n_threads,
    const pgphase_collect::GraphQueryConfig& query_config,
    const std::string& contig,
    const pgphase_collect::Options& filter_opts,
    GetCatalog get_catalog) {
    using namespace pgphase_collect;
    const size_t n = intervals.size();
    std::vector<GraphChunkWorkResult> results(n);
    const size_t n_workers = std::min<size_t>(static_cast<size_t>(n_threads), n);
    std::atomic<size_t> next{0};
    std::exception_ptr first_error;
    std::mutex err_mutex;
    std::vector<std::thread> workers;
    workers.reserve(n_workers);
    for (size_t wi = 0; wi < n_workers; ++wi) {
        workers.emplace_back([&]() {
            try {
                while (true) {
                    const size_t i = next.fetch_add(1);
                    if (i >= n) break;
                    GraphSiteCatalog chunk_catalog =
                        get_catalog(i, intervals[i].first, intervals[i].second);
                    std::vector<GraphReadAllele> chunk_rows =
                        collect_graph_read_alleles_for_catalog(query_config, chunk_catalog);
                    results[i].build_result = build_graph_bam_chunk(
                        chunk_catalog, chunk_rows, contig,
                        intervals[i].first, intervals[i].second,
                        static_cast<int>(i), filter_opts);
                    results[i].catalog = std::move(chunk_catalog);
                    results[i].rows    = std::move(chunk_rows);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(err_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);
    return results;
}

static bool parse_uint_field(std::string_view text, int& out) {
    if (text.empty()) return false;
    int value = 0;
    for (char c : text) {
        if (c < '0' || c > '9') return false;
        value = value * 10 + (c - '0');
    }
    out = value;
    return true;
}

static bool parse_gaf_name_mapq_from_columns(const std::string& line,
                                             int first_gaf_column,
                                             std::string& read_name,
                                             int& mapq) {
    const int qname_col = first_gaf_column;
    const int mapq_col = first_gaf_column + 11;
    int field = 0;
    size_t pos = 0;
    while (pos <= line.size()) {
        const size_t tab = line.find('\t', pos);
        const size_t end = tab == std::string::npos ? line.size() : tab;
        if (field == qname_col) read_name.assign(line, pos, end - pos);
        if (field == mapq_col) {
            return parse_uint_field(std::string_view(line.data() + pos, end - pos), mapq);
        }
        if (tab == std::string::npos) break;
        pos = tab + 1;
        ++field;
    }
    return false;
}

// Collect unique, MAPQ-filtered read names from raw or pggaf indexed GAF.
static std::vector<std::string> collect_gaf_read_names(const std::string& gaf_file, int min_mapq) {
    std::vector<std::string> names;
    htsFile* raw_fp = hts_open(gaf_file.c_str(), "r");
    if (raw_fp == nullptr) throw std::runtime_error("failed to open GAF: " + gaf_file);
    struct HtsFileGuard {
        htsFile* fp = nullptr;
        ~HtsFileGuard() { if (fp != nullptr) hts_close(fp); }
    } fp_guard{raw_fp};

    kstring_t line = KS_INITIALIZE;
    struct KStringGuard {
        kstring_t* line = nullptr;
        ~KStringGuard() { if (line != nullptr) ks_free(line); }
    } line_guard{&line};

    while (hts_getline(raw_fp, '\n', &line) >= 0) {
        if (line.l == 0 || line.s[0] == '#') continue;
        const std::string text(line.s, line.l);
        int mapq = 255;
        std::string read_name;
        if (!parse_gaf_name_mapq_from_columns(text, 3, read_name, mapq) &&
            !parse_gaf_name_mapq_from_columns(text, 0, read_name, mapq)) {
            continue;
        }
        if (mapq >= min_mapq && !read_name.empty()) names.push_back(std::move(read_name));
    }
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

hts_pos_t catalog_end_bound(const pgphase_collect::GraphSiteCatalog& catalog) {
    hts_pos_t end = 1;
    for (const pgphase_collect::GraphSite& site : catalog.sites) {
        end = std::max(end, site.ref_end > 0 ? site.ref_end : site.pos);
    }
    return end;
}

} // namespace

int phase_graph(int argc, char* argv[]) {
    using namespace pgphase_collect;
    GraphOptions opts;

    const struct option long_options[] = {
        {"gaf",                 required_argument, nullptr, kGafFileOption},
        {"gaf-file",            required_argument, nullptr, kGafFileOption},
        {"graph-sites-vcf",     required_argument, nullptr, kGraphSitesVcfOption},
        {"graph-sites-tsv",     required_argument, nullptr, kGraphSitesTsvOption},
        {"graph-read-support",  required_argument, nullptr, kGraphReadSupportTsvOption},
        {"graph-site-counts",   required_argument, nullptr, kGraphSiteCountsTsvOption},
        {"graph-read-profile",  required_argument, nullptr, kGraphReadProfileTsvOption},
        {"graph-phase-sites",   required_argument, nullptr, kGraphPhaseSitesTsvOption},
        {"graph-phase-reads",   required_argument, nullptr, kGraphPhaseReadsTsvOption},
        {"graph-filtered-sites",required_argument, nullptr, kGraphFilteredSitesTsvOption},
        {"graph-phase-bam",     required_argument, nullptr, kGraphPhaseBamOption},
        {"contig",              required_argument, nullptr, kGraphContigOption},
        {"interval",            required_argument, nullptr, kGraphIntervalOption},
        {"chunk-size",          required_argument, nullptr, kGraphChunkSizeOption},
        {"threads",             required_argument, nullptr, 't'},
        {"verbose",             required_argument, nullptr, 'V'},
        {"help",                no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    optind = 1;
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "t:V:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case kGafFileOption:            opts.gaf_file = optarg; break;
            case kGraphSitesVcfOption:      opts.graph_sites_vcf = optarg; break;
            case kGraphSitesTsvOption:      opts.graph_sites_tsv = optarg; break;
            case kGraphReadSupportTsvOption:opts.graph_read_support_tsv = optarg; break;
            case kGraphSiteCountsTsvOption: opts.graph_site_counts_tsv = optarg; break;
            case kGraphReadProfileTsvOption:opts.graph_read_profile_tsv = optarg; break;
            case kGraphPhaseSitesTsvOption: opts.graph_phase_sites_tsv = optarg; break;
            case kGraphPhaseReadsTsvOption:      opts.graph_phase_reads_tsv = optarg; break;
            case kGraphFilteredSitesTsvOption:   opts.graph_filtered_sites_tsv = optarg; break;
            case kGraphPhaseBamOption:           opts.graph_phase_bam = optarg; break;
            case kGraphContigOption:        opts.graph_contig = optarg; break;
            case kGraphIntervalOption:
                if (!parse_half_open_interval(optarg, opts.graph_interval_beg, opts.graph_interval_end)) {
                    std::cerr << "Error: --interval must be BEG..END with 0-based half-open coordinates\n";
                    return 1;
                }
                break;
            case kGraphChunkSizeOption: opts.graph_chunk_size = std::stoll(optarg); break;
            case 't': opts.threads = std::stoi(optarg); break;
            case 'V': opts.verbose = std::stoi(optarg); break;
            case 'h': print_help(); return 0;
            default: print_help(); return 1;
        }
    }

    if (opts.gaf_file.empty()) {
        std::cerr << "Error: --gaf is required\n";
        print_help();
        return 1;
    }
    if (opts.graph_sites_vcf.empty()) {
        std::cerr << "Error: --graph-sites-vcf is required\n";
        print_help();
        return 1;
    }
    if (opts.graph_chunk_size < 1 || opts.threads < 1 || opts.verbose < 0) {
        std::cerr << "Error: numeric thresholds are invalid\n";
        return 1;
    }

    try {
        pgphase_collect::Options filter_opts;
        filter_opts.verbose = opts.verbose;
        filter_opts.touch_read_phase = true;

        Options phase_opts;
        phase_opts.read_technology = ReadTechnology::Hifi;
        phase_opts.output_aln = "graph";
        phase_opts.verbose = opts.verbose;

        // ── Determine which contigs to process ────────────────────────────────
        // Single contig: use --contig directly (--interval applies only here).
        // Whole genome: discover all contigs from VCF (streaming, one at a time).
        std::vector<std::string> contigs;
        if (!opts.graph_contig.empty()) {
            contigs.push_back(opts.graph_contig);
        } else {
            contigs = load_graph_site_contig_names(opts.graph_sites_vcf);
            if (contigs.empty())
                throw std::runtime_error("no contigs found in VCF: " + opts.graph_sites_vcf);
            if (opts.verbose >= 1)
                std::cerr << "Discovered " << contigs.size() << " contig(s) in VCF\n";
        }

        // ── Open all output files once; write headers before the contig loop ──
        auto open_out = [](const std::string& path) -> std::ofstream {
            if (path.empty()) return {};
            std::ofstream f(path);
            if (!f) throw std::runtime_error("failed to open output: " + path);
            return f;
        };
        std::ofstream variants_out      = open_out(opts.graph_phase_sites_tsv);
        std::ofstream site_counts_out   = open_out(opts.graph_site_counts_tsv);
        std::ofstream profiles_out      = open_out(opts.graph_read_profile_tsv);
        std::ofstream filtered_out      = open_out(opts.graph_filtered_sites_tsv);
        std::ofstream read_support_out  = open_out(opts.graph_read_support_tsv);

        if (variants_out)     write_graph_bam_variants_tsv_header(variants_out);
        if (site_counts_out)  write_graph_bam_site_counts_tsv_header(site_counts_out);
        if (profiles_out)     write_graph_bam_read_profiles_tsv_header(profiles_out);
        if (filtered_out)     write_graph_bam_filtered_sites_tsv_header(filtered_out);
        if (read_support_out) write_graph_read_alleles_tsv_header(read_support_out);

        // ── Global accumulators across all contigs ────────────────────────────
        GraphSiteCatalog diagnostic_catalog;  // for graph_sites_tsv
        std::unordered_map<std::string, PhaseReadOutputRow> rows_by_read;
        size_t total_rows   = 0;
        size_t total_chunks = 0;
        int    chunk_id_offset = 0;

        const size_t nw = static_cast<size_t>(std::max(1, opts.threads));
        require_indexed_gaf(opts.gaf_file);

        // ── Per-contig streaming loop ─────────────────────────────────────────
        // For each contig:
        //   1. Load only that contig's sites from the VCF (tabix if indexed).
        //   2. Build per-chunk catalogs in the contig's own [0, end) position space.
        //   3. Query the indexed annotated GAF per chunk with tabix.
        //   4. Run the streaming build → assign → stitch → emit loop.
        // Peak catalog RAM = one contig's sites (not the whole genome at once).
        // Peak BamChunk RAM = O(n_threads × one_chunk).
        for (const std::string& contig : contigs) {
            // Build a RegionFilter that limits loading to this contig.
            // For the single-contig case, also respect --interval bounds.
            RegionFilter rf;
            rf.enabled = true;
            rf.chrom   = contig;
            if (contigs.size() == 1 && opts.graph_interval_beg >= 0) {
                // Convert 0-based half-open to the filter's 0-based-min-pos convention.
                rf.beg = opts.graph_interval_beg;
                rf.end = (opts.graph_interval_end > opts.graph_interval_beg)
                         ? opts.graph_interval_end : -1;
            } else {
                rf.beg = 0;
                rf.end = -1;
            }

            GraphSiteCatalog full_catalog =
                load_graph_site_catalog_from_vcf(opts.graph_sites_vcf, {rf});

            if (full_catalog.sites.empty()) {
                if (opts.verbose >= 2)
                    std::cerr << "Contig " << contig << ": no sites, skipping\n";
                continue;
            }

            if (opts.verbose >= 1)
                std::cerr << "Contig " << contig << ": querying indexed GAF for "
                          << full_catalog.sites.size() << " graph sites...\n";

            // Coordinate bounds for this contig (always 0-based half-open).
            const hts_pos_t interval_beg = rf.beg;
            const hts_pos_t interval_end = (rf.end > rf.beg)
                                           ? rf.end
                                           : catalog_end_bound(full_catalog);
            const auto intervals = split_graph_interval(
                interval_beg, interval_end, opts.graph_chunk_size);
            const size_t n = intervals.size();

            // Build per-chunk catalogs in O(n_sites) — one arithmetic pass.
            std::vector<GraphSiteCatalog> chunk_catalogs(n);
            for (const GraphSite& site : full_catalog.sites) {
                const hts_pos_t site_beg0 =
                    (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1;
                if (site_beg0 < interval_beg || site_beg0 >= interval_end) continue;
                const size_t ci = static_cast<size_t>(
                    (site_beg0 - interval_beg) / opts.graph_chunk_size);
                if (ci >= n) continue;
                chunk_catalogs[ci].sites.push_back(site);
            }

            // Build per-contig site_id → GraphSite* lookup (valid until end of
            // this loop iteration; freed with full_catalog).
            std::unordered_map<std::string, const GraphSite*> site_map;
            site_map.reserve(full_catalog.sites.size());
            for (const GraphSite& site : full_catalog.sites) {
                if (!site.id.empty() && site.id != ".")
                    site_map.emplace(site.id, &site);
            }

            // Streaming build → assign → stitch → emit loop (same as before,
            // but scoped to this contig).
            const size_t n_per_batch = nw;
            std::optional<GraphBamChunkBuildResult> prev_chunk;
            size_t contig_rows = 0;

            for (size_t batch_beg = 0; batch_beg < n; batch_beg += n_per_batch) {
                const size_t batch_end = std::min(batch_beg + n_per_batch, n);
                auto results = run_vcf_gaf_chunks_parallel(
                    chunk_catalogs, opts.gaf_file, intervals,
                    batch_beg, batch_end,
                    opts.threads, contig, filter_opts,
                    chunk_id_offset);

                for (size_t bi = 0; bi < results.size(); ++bi) {
                    append_catalog(diagnostic_catalog, results[bi].catalog);
                    contig_rows += results[bi].rows.size();
                    if (read_support_out)
                        write_graph_read_alleles_tsv_rows(read_support_out, results[bi].rows);
                    GraphBamChunkBuildResult cur = std::move(results[bi].build_result);

                    assign_graph_chunk_hap(cur, phase_opts);

                    // Build-time outputs: written before stitch.
                    if (site_counts_out)
                        write_graph_bam_site_counts_tsv_rows(site_counts_out, cur);
                    if (profiles_out)
                        write_graph_bam_read_profiles_tsv_rows(profiles_out, cur);
                    if (filtered_out)
                        write_graph_bam_filtered_sites_tsv_rows(filtered_out, cur);

                    // Stitch cur with prev; prev becomes fully finalized.
                    if (prev_chunk.has_value()) {
                        stitch_graph_chunk_pair(*prev_chunk, cur, phase_opts);
                        if (variants_out)
                            write_graph_bam_variants_tsv_rows(
                                variants_out, *prev_chunk, site_map);
                        merge_graph_chunk_into_read_rows(rows_by_read, *prev_chunk);
                    }

                    prev_chunk = std::move(cur);
                }
            }

            // Finalize the last chunk of this contig.
            if (prev_chunk.has_value()) {
                if (variants_out)
                    write_graph_bam_variants_tsv_rows(
                        variants_out, *prev_chunk, site_map);
                merge_graph_chunk_into_read_rows(rows_by_read, *prev_chunk);
                prev_chunk.reset();
            }

            total_rows += contig_rows;
            if (opts.verbose >= 1)
                std::cerr << "Contig " << contig << ": " << contig_rows
                          << " read-allele rows, " << n << " chunk(s)\n";
            chunk_id_offset += static_cast<int>(n);
            total_chunks    += n;
        } // end per-contig loop

        // ── Post-loop outputs ─────────────────────────────────────────────────

        // graph_sites_tsv: write once from the accumulated diagnostic catalog.
        if (!opts.graph_sites_tsv.empty()) {
            std::ofstream out(opts.graph_sites_tsv);
            if (!out) throw std::runtime_error(
                "failed to open graph-site output: " + opts.graph_sites_tsv);
            write_graph_site_catalog_tsv(out, diagnostic_catalog);
        }

        // Phase BAM / reads TSV: written from the accumulated per-read map.
        // A second GAF pass collects all MAPQ-passing read names so reads with
        // no site observations appear as unphased rows in the output.
        if (!opts.graph_phase_reads_tsv.empty() || !opts.graph_phase_bam.empty()) {
            const auto all_read_names =
                collect_gaf_read_names(opts.gaf_file, filter_opts.min_mapq);

            if (!opts.graph_phase_reads_tsv.empty()) {
                std::ofstream out(opts.graph_phase_reads_tsv);
                if (!out) throw std::runtime_error(
                    "failed to open graph phase-read output: " + opts.graph_phase_reads_tsv);
                write_graph_bam_phase_reads_tsv_from_rows(out, rows_by_read, all_read_names);
            }
            if (!opts.graph_phase_bam.empty()) {
                write_graph_bam_phase_bam_from_rows(
                    opts.graph_phase_bam, rows_by_read, all_read_names);
            }
        }

        if (opts.verbose >= 1) {
            std::cerr << "Graph phasing used " << diagnostic_catalog.sites.size()
                      << " graph sites across " << contigs.size() << " contig(s), "
                      << total_rows << " read-allele rows, and "
                      << total_chunks << " chunk(s)\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
