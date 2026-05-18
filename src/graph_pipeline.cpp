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
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <optional>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

namespace {

struct SamFileCloser {
    void operator()(samFile* fp) const {
        if (fp != nullptr) hts_close(fp);
    }
};

struct GraphOptions {
    std::string gaf_file;
    std::string graph_sites_vcf;
    std::string graph_sites_tsv = "graph_sites.tsv";
    std::string graph_phase_bam = "graph_phased.bam";
    std::string graph_filtered_sites_tsv;
    std::string graph_contig;
    hts_pos_t graph_interval_beg = -1;
    hts_pos_t graph_interval_end = -1;
    hts_pos_t graph_chunk_size = 500000;
    int threads = 1;
    int verbose = 0;
    bool ont = false;
};

enum LongOption {
    kGafFileOption = 1000,
    kGraphSitesVcfOption,
    kGraphSitesTsvOption,
    kGraphPhaseBamOption,
    kGraphFilteredSitesTsvOption,
    kGraphContigOption,
    kGraphIntervalOption,
    kGraphChunkSizeOption,
    kThreadsOption,
    kOntOption
};

void print_help() {
    std::cout
        << "Usage: pgphase phase-graph [options]\n"
        << "\n"
        << "Required:\n"
        << "      --gaf FILE                pggaf indexed annotated GAF (.gaf.gz + .tbi; rc/rb/re required)\n"
        << "      --graph-sites-vcf FILE    Decomposed graph-site VCF (bgzip + tabix .tbi/.csi required)\n"
        << "\n"
        << "Outputs (defaults shown):\n"
        << "      --graph-sites-tsv FILE      Site scaffold TSV (PHASE_SET, HAP1_ALLELE, HAP2_ALLELE, ...) [graph_sites.tsv]\n"
        << "      --graph-phase-bam FILE        Unaligned BAM with HP/PS tags per read [graph_phased.bam]\n"
        << "      --graph-filtered-sites FILE   Optional audit TSV: drops with explicit filter_reason "
           "[none]\n"
        << "\n"
        << "Options:\n"
        << "      --contig NAME             Restrict to this reference contig\n"
        << "      --interval BEG..END     0-based half-open [BEG,END); default = whole contig\n"
        << "                              (from ##contig length). Length <= --chunk-size => 1 chunk.\n"
        << "      --chunk-size INT          Target chunk width in bp [500000]\n"
        << "  -t, --threads INT             Worker threads [1]\n"
        << "      --ont                     ONT read mode (enables strand-bias filter)\n"
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
    // Half-open tiling [beg, end). Same recipe whether bounds come from --interval or
    // full contig: span wider than chunk_size → multiple chunks; span ≤ chunk_size → one chunk.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> chunks;
    for (hts_pos_t chunk_beg = beg; chunk_beg < end; chunk_beg += chunk_size) {
        chunks.emplace_back(chunk_beg, std::min(end, chunk_beg + chunk_size));
    }
    return chunks;
}

struct GraphChunkWorkResult {
    size_t read_allele_rows = 0;
    pgphase_collect::GraphBamChunkBuildResult build_result;
};

// Indexed graph-site VCF + indexed GAF: tabix site load per chunk + tabix GAF scan.
// One GraphSitesTabixReader per worker. Processes chunk indices [batch_beg, batch_end).
static std::vector<GraphChunkWorkResult> run_vcf_gaf_chunks_parallel(
    const std::string& sites_vcf_path,
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
                GraphSitesTabixReader reader(sites_vcf_path);
                while (true) {
                    const size_t offset = next.fetch_add(1);
                    if (offset >= batch_size) break;
                    const size_t i = batch_beg + offset;
                    GraphSiteCatalog chunk_catalog =
                        reader.load_half_open_interval(
                            contig,
                            intervals[i].first,
                            intervals[i].second,
                            true);
                    std::vector<GraphReadAllele> chunk_rows = scan_indexed_gaf_chunk(
                        indexed_gaf_file, contig,
                        intervals[i].first, intervals[i].second,
                        chunk_catalog, filter_opts.min_mapq);
                    results[offset].read_allele_rows = chunk_rows.size();
                    results[offset].build_result = build_graph_bam_chunk(
                        chunk_catalog, chunk_rows, contig,
                        intervals[i].first, intervals[i].second,
                        chunk_id_offset + static_cast<int>(i), filter_opts);
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

} // namespace

int phase_graph(int argc, char* argv[]) {
    using namespace pgphase_collect;
    GraphOptions opts;

    const struct option long_options[] = {
        {"gaf",                required_argument, nullptr, kGafFileOption},
        {"gaf-file",           required_argument, nullptr, kGafFileOption},
        {"graph-sites-vcf",    required_argument, nullptr, kGraphSitesVcfOption},
        {"graph-sites-tsv",    required_argument, nullptr, kGraphSitesTsvOption},
        {"graph-phase-bam",    required_argument, nullptr, kGraphPhaseBamOption},
        {"graph-filtered-sites", required_argument, nullptr, kGraphFilteredSitesTsvOption},
        {"contig",             required_argument, nullptr, kGraphContigOption},
        {"interval",           required_argument, nullptr, kGraphIntervalOption},
        {"chunk-size",         required_argument, nullptr, kGraphChunkSizeOption},
        {"threads",            required_argument, nullptr, 't'},
        {"ont",                no_argument,       nullptr, kOntOption},
        {"verbose",            required_argument, nullptr, 'V'},
        {"help",               no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    optind = 1;
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "t:V:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case kGafFileOption:           opts.gaf_file = optarg; break;
            case kGraphSitesVcfOption:     opts.graph_sites_vcf = optarg; break;
            case kGraphSitesTsvOption:     opts.graph_sites_tsv = optarg; break;
            case kGraphPhaseBamOption:     opts.graph_phase_bam = optarg; break;
            case kGraphFilteredSitesTsvOption:
                opts.graph_filtered_sites_tsv = optarg;
                break;
            case kGraphContigOption:       opts.graph_contig = optarg; break;
            case kGraphIntervalOption:
                if (!parse_half_open_interval(optarg, opts.graph_interval_beg, opts.graph_interval_end)) {
                    std::cerr << "Error: --interval must be BEG..END with 0-based half-open coordinates\n";
                    return 1;
                }
                break;
            case kGraphChunkSizeOption: opts.graph_chunk_size = std::stoll(optarg); break;
            case kOntOption: opts.ont = true; break;
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
    if (opts.graph_sites_tsv.empty() || opts.graph_phase_bam.empty()) {
        std::cerr << "Error: output paths must be non-empty\n";
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
        phase_opts.read_technology = opts.ont ? ReadTechnology::Ont : ReadTechnology::Hifi;
        phase_opts.output_aln = "graph";
        phase_opts.verbose = opts.verbose;

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

        std::ofstream sites_out(opts.graph_sites_tsv);
        if (!sites_out)
            throw std::runtime_error("failed to open graph sites output: " + opts.graph_sites_tsv);
        write_graph_bam_phase_sites_tsv_header(sites_out);

        std::unordered_set<std::string> phase_reads_emitted;
        std::unique_ptr<samFile, SamFileCloser> phase_bam_fp;
        std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> phase_bam_hdr(nullptr,
                                                                             sam_hdr_destroy);

        samFile* raw_fp = hts_open(opts.graph_phase_bam.c_str(), "wb");
        if (raw_fp == nullptr)
            throw std::runtime_error("failed to open output BAM: " + opts.graph_phase_bam);
        phase_bam_fp.reset(raw_fp);
        sam_hdr_t* hdr = sam_hdr_init();
        if (hdr == nullptr)
            throw std::runtime_error("failed to allocate phase-graph BAM header");
        phase_bam_hdr.reset(hdr);
        if (sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "unknown", nullptr) < 0)
            throw std::runtime_error("failed to build phase-graph BAM header");
        if (sam_hdr_write(phase_bam_fp.get(), phase_bam_hdr.get()) < 0)
            throw std::runtime_error("failed to write phase-graph BAM header");

        std::ofstream filtered_out;
        const bool emit_filtered_audit = !opts.graph_filtered_sites_tsv.empty();
        if (emit_filtered_audit) {
            filtered_out.open(opts.graph_filtered_sites_tsv);
            if (!filtered_out)
                throw std::runtime_error("failed to open filtered-site audit: " +
                                         opts.graph_filtered_sites_tsv);
            write_graph_bam_filtered_sites_tsv_header(filtered_out);
        }
        std::unordered_map<std::string, size_t> filter_reason_hist;

        std::unordered_map<std::string, PhaseReadOutputRow> rows_by_read;
        size_t total_rows         = 0;
        size_t total_chunks       = 0;
        size_t total_distinct_ids = 0;
        int    chunk_id_offset    = 0;

        const size_t nw = static_cast<size_t>(std::max(1, opts.threads));
        require_indexed_gaf(opts.gaf_file);
        require_graph_site_vcf_tabix_index(opts.graph_sites_vcf);

        // Per-contig: stitch chunks like collect-bam-variation; stream scaffold site rows + phase-BAM.
        for (const std::string& contig : contigs) {
            RegionFilter rf;
            rf.enabled = true;
            rf.chrom   = contig;
            if (contigs.size() == 1 && opts.graph_interval_beg >= 0) {
                rf.beg = opts.graph_interval_beg;
                rf.end = (opts.graph_interval_end > opts.graph_interval_beg)
                         ? opts.graph_interval_end : -1;
            } else {
                rf.beg = 0;
                rf.end = -1;
            }

            const hts_pos_t interval_beg = rf.beg;

            const hts_pos_t interval_end = (rf.end > rf.beg)
                                               ? rf.end
                                               : graph_site_contig_end_bp_from_vcf_header(
                                                     opts.graph_sites_vcf, contig);

            const auto intervals =
                split_graph_interval(interval_beg, interval_end, opts.graph_chunk_size);
            const size_t n = intervals.size();
            if (n == 0) continue;

            if (opts.verbose >= 1) {
                std::cerr << "Contig " << contig << ": streaming indexed sites, " << n
                          << " chunk(s)\n";
            }

            const size_t n_per_batch = nw;
            std::optional<GraphBamChunkBuildResult> prev_chunk;
            size_t contig_rows = 0;
            size_t contig_distinct_sites = 0;

            for (size_t batch_beg = 0; batch_beg < n; batch_beg += n_per_batch) {
                const size_t batch_end = std::min(batch_beg + n_per_batch, n);
                std::vector<GraphChunkWorkResult> results = run_vcf_gaf_chunks_parallel(
                    opts.graph_sites_vcf, opts.gaf_file, intervals,
                    batch_beg, batch_end, opts.threads, contig, filter_opts,
                    chunk_id_offset);

                std::stable_sort(results.begin(), results.end(),
                                 [](const GraphChunkWorkResult& a, const GraphChunkWorkResult& b) {
                                     return a.build_result.chunk.region.beg <
                                            b.build_result.chunk.region.beg;
                                 });

                for (size_t bi = 0; bi < results.size(); ++bi) {
                    contig_distinct_sites += results[bi].build_result.site_ids.size();
                    contig_rows += results[bi].read_allele_rows;
                    GraphBamChunkBuildResult cur = std::move(results[bi].build_result);

                    assign_graph_chunk_hap(cur, phase_opts);

                    for (const auto& fs : cur.filtered_sites)
                        ++filter_reason_hist[fs.filter_reason];
                    if (emit_filtered_audit)
                        write_graph_bam_filtered_sites_tsv_rows(filtered_out, cur);

                    if (prev_chunk.has_value()) {
                        stitch_graph_chunk_pair(*prev_chunk, cur, phase_opts);
                        write_graph_bam_phase_sites_tsv_rows(sites_out, *prev_chunk);
                        merge_graph_chunk_into_read_rows(rows_by_read, *prev_chunk);
                        std::unordered_set<std::string> next_qnames;
                        next_qnames.reserve(cur.chunk.reads.size() * 2 + 8);
                        for (const ReadRecord& rr : cur.chunk.reads)
                            next_qnames.insert(rr.qname);
                        flush_graph_phase_bam_after_merge(
                            phase_bam_fp.get(),
                            phase_bam_hdr.get(),
                            rows_by_read,
                            &next_qnames,
                            phase_reads_emitted);
                    }

                    prev_chunk = std::move(cur);
                }
            }

            if (prev_chunk.has_value()) {
                write_graph_bam_phase_sites_tsv_rows(sites_out, *prev_chunk);
                merge_graph_chunk_into_read_rows(rows_by_read, *prev_chunk);
                prev_chunk.reset();
            }
            flush_graph_phase_bam_after_merge(
                phase_bam_fp.get(),
                phase_bam_hdr.get(),
                rows_by_read,
                nullptr,
                phase_reads_emitted);

            total_rows += contig_rows;
            total_distinct_ids += contig_distinct_sites;
            if (opts.verbose >= 1)
                std::cerr << "Contig " << contig << ": " << contig_rows
                          << " read-allele rows, " << contig_distinct_sites
                          << " phased site row(s) (chunk-local)\n";
            chunk_id_offset += static_cast<int>(n);
            total_chunks += n;
        }

        const std::vector<std::string> all_read_names =
            collect_gaf_read_names(opts.gaf_file, filter_opts.min_mapq);
        std::vector<std::string> sorted_all(all_read_names.begin(), all_read_names.end());
        std::sort(sorted_all.begin(), sorted_all.end());
        write_graph_phase_bam_for_unobserved(
            phase_bam_fp.get(),
            phase_bam_hdr.get(),
            phase_reads_emitted,
            sorted_all);

        if (opts.verbose >= 1 && !filter_reason_hist.empty()) {
            size_t n_filtered_records = 0;
            for (const auto& kv : filter_reason_hist) n_filtered_records += kv.second;
            std::cerr << "Filtered-site audit (explicit drop rows): " << n_filtered_records
                      << " across all chunks\n";
            std::vector<std::pair<std::string, size_t>> pairs(filter_reason_hist.begin(),
                                                              filter_reason_hist.end());
            std::sort(pairs.begin(), pairs.end(),
                      [](const std::pair<std::string, size_t>& a,
                         const std::pair<std::string, size_t>& b) {
                          return a.second != b.second ? a.second > b.second : a.first < b.first;
                      });
            for (const auto& kv : pairs)
                std::cerr << "  " << kv.first << '\t' << kv.second << '\n';
        }

        if (opts.verbose >= 1) {
            std::cerr << "Graph phasing used " << total_distinct_ids
                      << " phased site row(s) across " << contigs.size()
                      << " contig(s), " << total_rows << " read-allele rows, and "
                      << total_chunks << " chunk(s)\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
