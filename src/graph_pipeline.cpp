#include "graph_pipeline.hpp"

#include "graph_bam_adapter.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <exception>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <cstring>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

// ---------- Per-chunk disk shard infrastructure (mirrors graph_collect.cpp) ----------

// Temp-file store: per_worker × per_chunk TSV files, written during GAF scan,
// read per-chunk during parallel build. Peak RAM = n_workers × one_chunk_rows.
struct ChunkRowShard {
    std::string dir;
    // paths[worker][chunk] = shard file path
    std::vector<std::vector<std::string>> paths;

    ChunkRowShard() = default;
    ChunkRowShard(const ChunkRowShard&) = delete;
    ChunkRowShard& operator=(const ChunkRowShard&) = delete;
    ~ChunkRowShard() {
        for (const auto& wp : paths)
            for (const auto& p : wp)
                std::remove(p.c_str());
        if (!dir.empty()) rmdir(dir.c_str());
    }

    // Load and return all rows for chunk_index from every worker's shard file.
    std::vector<pgphase_collect::GraphReadAllele> load_chunk(size_t chunk_index) const {
        using namespace pgphase_collect;
        std::vector<GraphReadAllele> rows;
        for (const auto& wp : paths) {
            if (chunk_index >= wp.size()) continue;
            std::ifstream in(wp[chunk_index]);
            if (!in) continue;
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty()) continue;
                // site_id \t chrom \t pos \t read_name \t mapq \t allele
                const char* s = line.c_str();
                const char* e = s + line.size();
                auto next = [&](const char* from) -> std::pair<std::string_view, const char*> {
                    const char* t = static_cast<const char*>(
                        std::memchr(from, '\t', static_cast<size_t>(e - from)));
                    if (!t) t = e;
                    return {{from, static_cast<size_t>(t - from)}, t < e ? t + 1 : e};
                };
                auto [site_id, p1] = next(s);
                auto [chrom,   p2] = next(p1);
                auto [pos_sv,  p3] = next(p2);
                auto [rname,   p4] = next(p3);
                auto [mapq_sv, p5] = next(p4);
                auto [ale_sv,  p6] = next(p5);
                rows.push_back(GraphReadAllele{
                    std::string(site_id),
                    std::string(chrom),
                    static_cast<hts_pos_t>(std::stoll(std::string(pos_sv))),
                    std::string(rname),
                    std::stoi(std::string(ale_sv)),
                    std::string(),
                    std::stoi(std::string(mapq_sv))
                });
            }
        }
        return rows;
    }
};

// LRU-cached append writer: keeps up to kMaxOpen file descriptors open at once.
class ChunkRowWriter {
public:
    explicit ChunkRowWriter(const std::vector<std::string>& paths) : paths_(paths) {}

    void write(size_t chunk_idx, const pgphase_collect::GraphReadAllele& row) {
        std::ofstream& out = stream_for(chunk_idx);
        out << row.site_id << '\t' << row.chrom << '\t' << row.pos << '\t'
            << row.read_name << '\t' << row.mapq << '\t' << row.allele << '\n';
    }

private:
    static constexpr size_t kMaxOpen = 64;
    struct Entry { std::ofstream stream; std::list<size_t>::iterator lru_it; };

    std::ofstream& stream_for(size_t idx) {
        auto it = open_.find(idx);
        if (it != open_.end()) {
            lru_.erase(it->second.lru_it);
            lru_.push_front(idx);
            it->second.lru_it = lru_.begin();
            return it->second.stream;
        }
        if (open_.size() >= kMaxOpen) {
            open_.erase(lru_.back());
            lru_.pop_back();
        }
        lru_.push_front(idx);
        Entry e{std::ofstream(paths_[idx], std::ios::app), lru_.begin()};
        if (!e.stream)
            throw std::runtime_error("failed to open chunk shard: " + paths_[idx]);
        auto [ins_it, ok] = open_.emplace(idx, std::move(e));
        (void)ok;
        return ins_it->second.stream;
    }

    const std::vector<std::string>& paths_;
    std::unordered_map<size_t, Entry> open_;
    std::list<size_t> lru_;
};

// ---------------------------------------------------------------------------------

struct GraphOptions {
    std::string gaf_file;
    std::string graph_sites_vcf;
    std::string graph_sites_tsv;
    std::string graph_read_support_tsv;
    std::string graph_site_counts_tsv;
    std::string graph_read_profile_tsv;
    std::string graph_phase_sites_tsv = "graph_phase_sites.tsv";
    std::string graph_phase_reads_tsv = "graph_phase_reads.tsv";
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
        << "      --gaf FILE                Raw GAF graph alignments\n"
        << "      --graph-sites-vcf FILE    Decomposed graph-site VCF (from build-snarl-catalog)\n"
        << "\n"
        << "Options:\n"
        << "      --contig NAME             Restrict to this reference contig\n"
        << "      --interval BEG..END       0-based half-open coordinate bounds (default: full catalog)\n"
        << "      --chunk-size INT          Phasing chunk size in bp [500000]\n"
        << "      --graph-phase-sites FILE  Output: graph-site hap consensus TSV [graph_phase_sites.tsv]\n"
        << "      --graph-phase-reads FILE  Output: read HAP / PHASE_SET TSV [graph_phase_reads.tsv]\n"
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

// VCF+GAF fast path: catalogs are pre-sliced; rows are loaded per-chunk from disk shards.
// Each worker loads only the rows for the chunk it is currently building, keeping
// peak RAM proportional to (n_workers × one_chunk_rows) instead of all rows at once.
static std::vector<GraphChunkWorkResult> run_vcf_gaf_chunks_parallel(
    const std::vector<pgphase_collect::GraphSiteCatalog>& chunk_catalogs,
    const ChunkRowShard& shard,
    const std::vector<std::pair<hts_pos_t, hts_pos_t>>& intervals,
    int n_threads,
    const std::string& contig,
    const pgphase_collect::Options& filter_opts) {
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
                    const std::vector<GraphReadAllele> chunk_rows = shard.load_chunk(i);
                    results[i].build_result = build_graph_bam_chunk(
                        chunk_catalogs[i], chunk_rows, contig,
                        intervals[i].first, intervals[i].second,
                        static_cast<int>(i), filter_opts);
                    results[i].catalog = chunk_catalogs[i];
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

void write_required_outputs(const GraphOptions& opts,
                            const pgphase_collect::GraphSiteCatalog& catalog,
                            const ChunkRowShard& shard,
                            const std::vector<pgphase_collect::GraphBamChunkBuildResult>& graph_chunks,
                            int min_mapq) {
    using namespace pgphase_collect;

    if (!opts.graph_sites_tsv.empty()) {
        std::ofstream out(opts.graph_sites_tsv);
        if (!out) throw std::runtime_error("failed to open graph-site output: " + opts.graph_sites_tsv);
        write_graph_site_catalog_tsv(out, catalog);
    }
    if (!opts.graph_read_support_tsv.empty()) {
        std::ofstream out(opts.graph_read_support_tsv);
        if (!out) throw std::runtime_error("failed to open graph read-support output: " + opts.graph_read_support_tsv);
        std::vector<GraphReadAllele> all_rows;
        const size_t n_chunks = shard.paths.empty() ? 0 : shard.paths[0].size();
        for (size_t ci = 0; ci < n_chunks; ++ci) {
            const std::vector<GraphReadAllele> chunk_rows = shard.load_chunk(ci);
            all_rows.insert(all_rows.end(), chunk_rows.begin(), chunk_rows.end());
        }
        write_graph_read_alleles_tsv(out, all_rows);
    }
    if (!opts.graph_site_counts_tsv.empty()) {
        std::ofstream out(opts.graph_site_counts_tsv);
        if (!out) throw std::runtime_error("failed to open graph site-count output: " + opts.graph_site_counts_tsv);
        write_graph_bam_site_counts_tsv(out, graph_chunks);
    }
    if (!opts.graph_read_profile_tsv.empty()) {
        std::ofstream out(opts.graph_read_profile_tsv);
        if (!out) throw std::runtime_error("failed to open graph read-profile output: " + opts.graph_read_profile_tsv);
        write_graph_bam_read_profiles_tsv(out, graph_chunks);
    }
    if (!opts.graph_filtered_sites_tsv.empty()) {
        std::ofstream out(opts.graph_filtered_sites_tsv);
        if (!out) throw std::runtime_error("failed to open graph filtered-sites output: " + opts.graph_filtered_sites_tsv);
        write_graph_bam_filtered_sites_tsv(out, graph_chunks);
    }

    {
        std::ofstream out(opts.graph_phase_sites_tsv);
        if (!out) throw std::runtime_error("failed to open graph phase-site output: " + opts.graph_phase_sites_tsv);
        write_graph_bam_variants_tsv(out, graph_chunks, catalog);
    }
    {
        // Fast second pass: collect every read name that passed the MAPQ filter,
        // so reads with no site observations can be emitted as hap=0.
        std::vector<std::string> all_read_names;
        {
            std::ifstream gaf(opts.gaf_file);
            std::string line;
            while (std::getline(gaf, line)) {
                if (line.empty() || line[0] == '#') continue;
                // GAF col 0 = read name, col 11 (0-indexed) = MAPQ.
                size_t pos = 0, field = 0;
                int mapq = 255;
                std::string read_name;
                while (true) {
                    const size_t tab = line.find('\t', pos);
                    const size_t len = (tab == std::string::npos) ? line.size() - pos : tab - pos;
                    if (field == 0) read_name.assign(line, pos, len);
                    else if (field == 11) { mapq = std::atoi(line.c_str() + pos); break; }
                    if (tab == std::string::npos) break;
                    pos = tab + 1;
                    ++field;
                }
                if (mapq >= min_mapq && !read_name.empty())
                    all_read_names.push_back(std::move(read_name));
            }
            std::sort(all_read_names.begin(), all_read_names.end());
            all_read_names.erase(
                std::unique(all_read_names.begin(), all_read_names.end()),
                all_read_names.end());
        }
        std::ofstream out(opts.graph_phase_reads_tsv);
        if (!out) throw std::runtime_error("failed to open graph phase-read output: " + opts.graph_phase_reads_tsv);
        write_graph_bam_phase_reads_tsv(out, graph_chunks, all_read_names);
    }
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

        GraphSiteCatalog catalog;
        std::vector<GraphBamChunkBuildResult> graph_chunks;

        // Load catalog; split into chunk-size intervals; pre-route GAF rows per chunk.
        GraphSiteCatalog full_catalog = load_graph_site_catalog_from_vcf(opts.graph_sites_vcf);
        if (opts.verbose >= 1)
            std::cerr << "Scanning GAF for " << full_catalog.sites.size() << " graph sites...\n";

        const hts_pos_t interval_beg = opts.graph_interval_beg >= 0 ? opts.graph_interval_beg : 0;
        const hts_pos_t interval_end = (opts.graph_interval_end > opts.graph_interval_beg)
                                       ? opts.graph_interval_end
                                       : catalog_end_bound(full_catalog);
        const auto intervals = split_graph_interval(interval_beg, interval_end, opts.graph_chunk_size);
        const size_t n = intervals.size();
        const size_t nw = static_cast<size_t>(std::max(1, opts.threads));

        // Build per-chunk catalogs in O(n_sites) — one pass that assigns each site by
        // arithmetic to its chunk interval, avoiding n_chunks × n_sites iterations.
        std::vector<GraphSiteCatalog> chunk_catalogs(n);
        std::unordered_map<std::string, size_t> site_to_chunk;
        site_to_chunk.reserve(full_catalog.sites.size());
        for (size_t si = 0; si < full_catalog.sites.size(); ++si) {
            const GraphSite& site = full_catalog.sites[si];
            if (!opts.graph_contig.empty()) {
                const std::string& sc = site.ref_contig.empty() ? site.chrom : site.ref_contig;
                if (sc != opts.graph_contig) continue;
            }
            const hts_pos_t site_beg0 = (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1;
            if (site_beg0 < interval_beg || site_beg0 >= interval_end) continue;
            const size_t ci = static_cast<size_t>(
                (site_beg0 - interval_beg) / opts.graph_chunk_size);
            if (ci >= n) continue;
            chunk_catalogs[ci].sites.push_back(site);
            {
                const std::string key = (!site.id.empty() && site.id != ".")
                    ? site.id
                    : site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
                site_to_chunk.emplace(key, ci);
            }
        }

        // Release walk string storage from chunk_catalogs — walk sequences were needed
        // only to build the compact index for the GAF scan. Walk vector sizes (allele
        // counts) are preserved so n_alleles checks in build_graph_bam_chunk still work.
        for (auto& cc : chunk_catalogs) {
            for (GraphSite& site : cc.sites) {
                for (GraphWalk& w : site.allele_walks) w.clear();
                site.allele_traversals.clear();
            }
        }

        // Create per-(worker, chunk) shard files in a temp directory.
        // During the GAF scan each worker writes only to its own files — no locking.
        // During the parallel build each worker loads its chunk's files from disk,
        // so peak RAM = n_workers × one_chunk_rows instead of all rows at once.
        static std::atomic<int> shard_counter{0};
        const std::string shard_dir =
            "/tmp/pgphase_graph_shard_" + std::to_string(getpid()) +
            "_" + std::to_string(++shard_counter);
        if (mkdir(shard_dir.c_str(), 0700) != 0)
            throw std::runtime_error("failed to create shard dir: " + shard_dir);

        ChunkRowShard shard;
        shard.dir = shard_dir;
        shard.paths.resize(nw);
        for (size_t wi = 0; wi < nw; ++wi) {
            shard.paths[wi].resize(n);
            for (size_t ci = 0; ci < n; ++ci)
                shard.paths[wi][ci] = shard_dir + "/w" + std::to_string(wi) +
                                      "_c" + std::to_string(ci) + ".tsv";
        }

        // Scan GAF once; each worker writes matched rows to its per-chunk shard file.
        // Using the releasing variant frees walk string storage from full_catalog during
        // the scan, reducing peak RAM (chunk_catalogs retain their own walk copies).
        std::vector<std::unique_ptr<ChunkRowWriter>> writers(nw);
        for (size_t wi = 0; wi < nw; ++wi)
            writers[wi] = std::make_unique<ChunkRowWriter>(shard.paths[wi]);

        std::atomic<size_t> row_counter{0};
        scan_gaf_for_catalog_emit_parallel_releasing_walks(
            opts.gaf_file, full_catalog, filter_opts.min_mapq, nw,
            [&](size_t worker_id, GraphReadAllele&& row) {
                auto it = site_to_chunk.find(row.site_id);
                if (it != site_to_chunk.end()) {
                    writers[worker_id]->write(it->second, row);
                    row_counter.fetch_add(1, std::memory_order_relaxed);
                }
            });
        writers.clear();  // flush and close all shard files

        const size_t total_rows = row_counter.load();
        if (opts.verbose >= 1)
            std::cerr << "GAF scan produced " << total_rows << " read-allele rows\n";

        auto results = run_vcf_gaf_chunks_parallel(chunk_catalogs, shard, intervals,
                                                   opts.threads, opts.graph_contig, filter_opts);
        graph_chunks.reserve(results.size());
        for (auto& r : results) {
            append_catalog(catalog, r.catalog);
            graph_chunks.push_back(std::move(r.build_result));
        }

        Options phase_opts;
        phase_opts.read_technology = ReadTechnology::Hifi;
        phase_opts.output_aln = "graph";
        phase_opts.verbose = opts.verbose;
        phase_graph_bam_chunks(graph_chunks, phase_opts);
        write_required_outputs(opts, catalog, shard, graph_chunks, filter_opts.min_mapq);

        if (opts.verbose >= 1) {
            std::cerr << "Graph phasing used " << catalog.sites.size()
                      << " graph sites, " << total_rows
                      << " exact read-allele rows, and "
                      << graph_chunks.size() << " graph chunk(s)\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
