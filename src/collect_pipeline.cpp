#include "collect_pipeline.hpp"

#include "bam_digar.hpp"
#include "collect_output.hpp"
#include "collect_var.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace pgphase_collect {

/**
 * Region construction, parallel chunk execution, streaming writer, and CLI for
 * `collect-bam-variation`. Coordinates are 1-based inclusive (BED is converted in
 * load_bed_regions). Chunks are batched by `reg_chunk_i` in run_collect_bam_variation
 * so TSV/VCF/read-support can be emitted without storing all candidates at once.
 */

// ════════════════════════════════════════════════════════════════════════════
// Region chunking
// ════════════════════════════════════════════════════════════════════════════

/** Parse `chr`, `chr:start`, or `chr:start-end` (commas allowed). Empty string → disabled filter. */
RegionFilter parse_region(const std::string& region) {
    RegionFilter filter;
    if (region.empty()) return filter;

    auto strip_commas = [](std::string value) {
        value.erase(std::remove(value.begin(), value.end(), ','), value.end());
        return value;
    };

    filter.enabled = true;
    const size_t colon = region.find(':');
    if (colon == std::string::npos) {
        filter.chrom = region;
        return filter;
    }

    filter.chrom = region.substr(0, colon);
    const std::string range = region.substr(colon + 1);
    const size_t dash = range.find('-');
    const std::string beg = strip_commas(dash == std::string::npos ? range : range.substr(0, dash));
    const std::string end = strip_commas(dash == std::string::npos ? "" : range.substr(dash + 1));
    if (!beg.empty()) filter.beg = std::stoll(beg);
    if (!end.empty()) filter.end = std::stoll(end);
    if (filter.chrom.empty() || filter.beg < 1 || (filter.end >= 0 && filter.end < filter.beg)) {
        throw std::runtime_error("invalid region: " + region);
    }
    return filter;
}

/** Load BED chrom/start/end; convert 0-based half-open to 1-based inclusive [bed_beg+1, bed_end]. */
std::vector<RegionFilter> load_bed_regions(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open region file: " + path);
    std::vector<RegionFilter> regions;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string chrom;
        hts_pos_t bed_beg = 0;
        hts_pos_t bed_end = 0;
        if (!(fields >> chrom >> bed_beg >> bed_end)) {
            throw std::runtime_error("invalid BED line in region file: " + line);
        }
        if (bed_beg < 0 || bed_end <= bed_beg) {
            throw std::runtime_error("invalid BED interval in region file: " + line);
        }
        regions.push_back(RegionFilter{true, chrom, bed_beg + 1, bed_end});
    }
    return regions;
}

/** BAM SQ name lookup; -1 if absent. */
static int tid_for_name(const bam_hdr_t* header, const std::string& chrom) {
    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (chrom == header->target_name[tid]) return tid;
    }
    return -1;
}

/** Tile [beg,end] into windows of up to chunk_size bp (last may be shorter). */
static std::vector<RegionChunk> split_region(int tid,
                                             hts_pos_t beg,
                                             hts_pos_t end,
                                             hts_pos_t chunk_size) {
    std::vector<RegionChunk> chunks;
    for (hts_pos_t chunk_beg = beg; chunk_beg <= end; chunk_beg += chunk_size) {
        const hts_pos_t chunk_end = std::min(end, chunk_beg + chunk_size - 1);
        RegionChunk chunk;
        chunk.tid = tid;
        chunk.beg = chunk_beg;
        chunk.end = chunk_end;
        chunks.push_back(chunk);
    }
    return chunks;
}

/** Sort by contig/start; set chunk_id, reg_chunk_i, reg_i, prev/next neighbour geometry. */
static void annotate_chunk_neighbors(std::vector<RegionChunk>& chunks) {
    std::sort(chunks.begin(), chunks.end(), [](const RegionChunk& lhs, const RegionChunk& rhs) {
        if (lhs.tid != rhs.tid) return lhs.tid < rhs.tid;
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        return lhs.end < rhs.end;
    });

    int reg_chunk_i = -1;
    int reg_i = 0;
    int last_tid = -1;
    for (size_t i = 0; i < chunks.size(); ++i) {
        RegionChunk& chunk = chunks[i];
        chunk.chunk_id = static_cast<int>(i);
        if (i == 0 || chunk.tid != last_tid) {
            ++reg_chunk_i;
            reg_i = 0;
        }
        chunk.reg_chunk_i = reg_chunk_i;
        chunk.reg_i = reg_i++;
        last_tid = chunk.tid;
    }
    for (size_t i = 0; i < chunks.size(); ++i) {
        RegionChunk& chunk = chunks[i];
        if (i > 0 && chunks[i - 1].tid == chunk.tid) {
            chunk.prev_chunk_id = chunks[i - 1].chunk_id;
            chunk.prev_tid = chunks[i - 1].tid;
            chunk.prev_beg = chunks[i - 1].beg;
            chunk.prev_end = chunks[i - 1].end;
        }
        if (i + 1 < chunks.size() && chunks[i + 1].tid == chunk.tid) {
            chunk.next_chunk_id = chunks[i + 1].chunk_id;
            chunk.next_tid = chunks[i + 1].tid;
            chunk.next_beg = chunks[i + 1].beg;
            chunk.next_end = chunks[i + 1].end;
        }
    }
}

/** Append chunks for [filter.beg, filter.end] clipped to contig (end=-1 → full contig). */
static void add_filter_chunks(const RegionFilter& region,
                              const bam_hdr_t* header,
                              const faidx_t* fai,
                              hts_pos_t chunk_size,
                              std::vector<RegionChunk>& chunks) {
    if (!region.enabled) return;
    const int tid = tid_for_name(header, region.chrom);
    if (tid < 0)
        throw std::runtime_error("region contig is not present in BAM header: " + region.chrom);
    if (!faidx_has_seq(fai, region.chrom.c_str()))
        throw std::runtime_error("region contig is not present in FASTA: " + region.chrom);

    const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
    const hts_pos_t end = region.end < 0 ? contig_end : std::min(region.end, contig_end);
    if (region.beg > end) return;
    auto region_chunks = split_region(tid, region.beg, end, chunk_size);
    chunks.insert(chunks.end(), region_chunks.begin(), region_chunks.end());
}

/** Combine -r regions, optional --region-file BED, and optional --autosome filters. */
static std::vector<RegionFilter> collect_region_filters(const Options& opts,
                                                        const bam_hdr_t* header,
                                                        const faidx_t* fai) {
    std::vector<RegionFilter> filters;
    for (const std::string& region : opts.regions) {
        filters.push_back(parse_region(region));
    }
    if (!opts.region_file.empty()) {
        auto bed_regions = load_bed_regions(opts.region_file);
        filters.insert(filters.end(), bed_regions.begin(), bed_regions.end());
    }
    if (opts.autosome) {
        for (int i = 1; i <= 22; ++i) {
            const std::string no_prefix = std::to_string(i);
            const std::string with_prefix = "chr" + no_prefix;
            if (tid_for_name(header, with_prefix) >= 0 && faidx_has_seq(fai, with_prefix.c_str())) {
                filters.push_back(RegionFilter{true, with_prefix, 1, -1});
            } else if (tid_for_name(header, no_prefix) >= 0 &&
                       faidx_has_seq(fai, no_prefix.c_str())) {
                filters.push_back(RegionFilter{true, no_prefix, 1, -1});
            }
        }
    }
    return filters;
}

/**
 * With filters: only those intervals. Without: every indexed FASTA contig. Always annotate_chunk_neighbors.
 */
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai) {
    const std::vector<RegionFilter> filters = collect_region_filters(opts, header, fai);
    std::vector<RegionChunk> chunks;

    if (!filters.empty()) {
        for (const RegionFilter& filter : filters) {
            add_filter_chunks(filter, header, fai, opts.chunk_size, chunks);
        }
        annotate_chunk_neighbors(chunks);
        return chunks;
    }

    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (!faidx_has_seq(fai, header->target_name[tid])) continue;
        const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
        if (contig_end <= 0) continue;
        auto contig_chunks = split_region(tid, 1, contig_end, opts.chunk_size);
        chunks.insert(chunks.end(), contig_chunks.begin(), contig_chunks.end());
    }

    annotate_chunk_neighbors(chunks);
    return chunks;
}

/** Primary BAM + index + reference FAI, then build_region_chunks. */
std::vector<RegionChunk> load_region_chunks(const Options& opts) {
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    std::unique_ptr<hts_idx_t, IndexDeleter> index(
        sam_index_load(bam.get(), opts.primary_bam_file().c_str()));
    if (!index) throw std::runtime_error("region chunking requires an indexed BAM/CRAM");
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));
    return build_region_chunks(opts, header.get(), fai.get());
}

// ════════════════════════════════════════════════════════════════════════════
// Pipeline
// ════════════════════════════════════════════════════════════════════════════

/**
 * Per-chunk arrays for downstream phasing hooks: overlap read indices per input BAM, skip counts,
 * placeholder haplotype / phase_set state (longcallD-shaped; not used by candidate-only output).
 */
static void initialize_longcalld_chunk_state(BamChunk& chunk, size_t n_bams) {
    chunk.up_ovlp_read_i.assign(n_bams, {});
    chunk.down_ovlp_read_i.assign(n_bams, {});
    chunk.n_up_ovlp_skip_reads.assign(n_bams, 0);
    chunk.n_down_ovlp_skip_reads.assign(n_bams, 0);
    chunk.haps.assign(chunk.reads.size(), 0);
    chunk.phase_scores.assign(chunk.reads.size(), 0);
    chunk.phase_sets.assign(chunk.reads.size(), -1);

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        if (read.input_index < 0 || static_cast<size_t>(read.input_index) >= n_bams) continue;
        if (read_overlaps_prev_region(chunk.region, read.tid, read.beg, read.end)) {
            chunk.up_ovlp_read_i[read.input_index].push_back(static_cast<int>(read_i));
        }
        if (read_overlaps_next_region(chunk.region, read.tid, read.beg, read.end)) {
            chunk.down_ovlp_read_i[read.input_index].push_back(static_cast<int>(read_i));
        }
    }
}

/**
 * @brief Loads all input BAMs for one region chunk and finishes digar/reference setup.
 *
 * Queries each BAM for overlapping reads, merges and sorts `chunk.reads`, initializes
 * per-read overlap indices and placeholder phasing fields, then calls `finalize_bam_chunk`.
 *
 * @param chunk Chunk state; `chunk.region` must already be set.
 * @param opts Read filters, reference path, and technology flags.
 * @param context Open BAM handles, indexes, headers, and shared reference cache.
 */
static void load_and_prepare_chunk(BamChunk& chunk, const Options& opts, WorkerContext& context) {
    chunk.reads.clear();
    std::vector<OverlapSkipCounts> overlap_skip_counts(context.bams.size());
    for (size_t input_i = 0; input_i < context.bams.size(); ++input_i) {
        std::vector<ReadRecord> reads = load_read_records_for_chunk(
            opts,
            chunk.region,
            static_cast<int>(input_i),
            *context.bams[input_i],
            context.headers[input_i].get(),
            context.indexes[input_i].get(),
            context.ref,
            &overlap_skip_counts[input_i]);
        chunk.reads.reserve(chunk.reads.size() + reads.size());
        chunk.reads.insert(chunk.reads.end(),
                           std::make_move_iterator(reads.begin()),
                           std::make_move_iterator(reads.end()));
    }
    std::sort(chunk.reads.begin(), chunk.reads.end(), [](const ReadRecord& lhs, const ReadRecord& rhs) {
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        if (lhs.end != rhs.end) return lhs.end > rhs.end;
        if (lhs.nm != rhs.nm) return lhs.nm < rhs.nm;
        return lhs.qname < rhs.qname;
    });
    initialize_longcalld_chunk_state(chunk, context.bams.size());
    for (size_t input_i = 0; input_i < overlap_skip_counts.size(); ++input_i) {
        chunk.n_up_ovlp_skip_reads[input_i] = overlap_skip_counts[input_i].upstream;
        chunk.n_down_ovlp_skip_reads[input_i] = overlap_skip_counts[input_i].downstream;
    }
    finalize_bam_chunk(chunk, context.ref, context.primary_header());
}

/** Noisy-region prep → unique candidate sites → allele/read-support counts. */
static void collect_prephase_candidates(BamChunk& chunk,
                                        const Options& opts,
                                        std::vector<ReadSupportRow>* read_support_out) {
    pre_process_noisy_regs_pgphase(chunk, opts);
    collect_candidate_sites_from_records(chunk.region, chunk.reads, chunk.candidates);
    collect_allele_counts_from_records(
        chunk.reads, chunk.candidates, &chunk.region, read_support_out, opts.min_bq);
}

/** classify_chunk_candidates → post_process noisy → optional noisy containment (non-ONT). */
static void classify_and_filter_candidates(BamChunk& chunk,
                                           const Options& opts,
                                           const bam_hdr_t* header) {
    classify_chunk_candidates(chunk, opts, header);
    post_process_noisy_regs_pgphase(chunk, chunk.candidates);
    if (!opts.is_ont()) apply_noisy_containment_filter(chunk);
}

/**
 * @brief Sequentially collects variants from candidate regions.
 *
 * Implements the core multi-stage worker thread map-reduce pattern:
 * It pulls bam portions by region coordinates, identifies individual variants via 
 * parsed CIGAR operations/CS tags, and assigns them categories (HET, STRAND BIAS, 
 * LOW COVERAGE, etc) prior to merging the sub-chunks together.
 *
 * This implements the overarching control flow found in `longcallD`'s 
 * top-level `collect_var_main()`.
 *
 * @param region Sub-region coordinate struct mapping a piece of genome.
 * @param opts Options describing behavior and qualities.
 * @param context Thread-local variables passing caching pointers safely.
 * @return CandidateTable A self-contained list of evaluated sites.
 */
static BamChunk process_chunk(const RegionChunk& region,
                              const Options& opts,
                              WorkerContext& context,
                              std::vector<ReadSupportRow>* read_support_out) {
    BamChunk chunk;
    chunk.region = region;
    load_and_prepare_chunk(chunk, opts, context);
    collect_prephase_candidates(chunk, opts, read_support_out);
    classify_and_filter_candidates(chunk, opts, context.primary_header());
    return chunk;
}

/**
 * @brief Simple memory-efficient concatenation of threaded variant results.
 *
 * Emulates the array `memcpy` structure when chunks synchronize. Rather than 
 * maintaining overlapping complex locks over a giant data structure, each thread 
 * independently allocates `BamChunk` variants and then moves (`std::make_move_iterator`) 
 * them deterministically and sequentially so genomic order is perfectly maintained.
 *
 * Does not re-evaluate `fuzzy_collapse` boundary conflicts (which mimics 
 * exactly `longcallD`'s static logic of collapsing purely inside chunks).
 *
 * @param chunks The completed chunk outputs evaluated by thread workers.
 * @return A unified table mapping the whole runtime block.
 */
static CandidateTable merge_chunk_candidates(std::vector<BamChunk>& chunks) {
    CandidateTable merged;
    for (BamChunk& chunk : chunks) {
        CandidateTable& table = chunk.candidates;
        merged.insert(merged.end(),
                      std::make_move_iterator(table.begin()),
                      std::make_move_iterator(table.end()));
    }
    return merged;
}

/**
 * @brief Outputs of parallel chunk processing for one `reg_chunk_i` batch.
 */
struct ChunkBatchResult {
    std::vector<BamChunk> chunks;
    std::vector<std::vector<ReadSupportRow>> read_support_batches;
};

/**
 * @brief Runs `process_chunk` on a contiguous slice of region chunks with a thread pool.
 *
 * Each worker constructs its own `WorkerContext` (per-thread BAM + FAI handles). On failure,
 * the first exception is stored and rethrown after all workers join.
 *
 * @param opts Thread count and I/O options.
 * @param chunks Full chunk list; only indices `[batch_begin, batch_end)` are processed.
 * @param batch_begin First index in `chunks` (inclusive).
 * @param batch_end One past the last index (exclusive).
 * @param collect_read_support If true, each slot in `read_support_batches` receives that chunk's rows.
 * @return Per-chunk `BamChunk` results in offset order, plus optional read-support batches.
 */
static ChunkBatchResult collect_chunk_batch_parallel(const Options& opts,
                                                     const std::vector<RegionChunk>& chunks,
                                                     size_t batch_begin,
                                                     size_t batch_end,
                                                     bool collect_read_support) {
    if (batch_begin > batch_end || batch_end > chunks.size()) {
        throw std::runtime_error("invalid chunk batch range");
    }
    const size_t batch_size = batch_end - batch_begin;
    ChunkBatchResult result;
    result.chunks.resize(batch_size);
    if (collect_read_support) result.read_support_batches.resize(batch_size);
    if (batch_size == 0) return result;

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&, worker_i]() {
            (void)worker_i;
            try {
                WorkerContext context(opts);
                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    std::vector<ReadSupportRow>* rs_ptr = nullptr;
                    if (collect_read_support) rs_ptr = &result.read_support_batches[offset];
                    result.chunks[offset] =
                        process_chunk(chunks[batch_begin + offset], opts, context, rs_ptr);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }

    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);

    return result;
}

/**
 * @brief Dispatcher driving the threadpool executing BAM collection logic.
 *
 * Implements a modern C++ thread dispatch architecture managing identical loops
 * defined by `collect_ref_seq_bam_main()` nested looping over regions.
 *
 * Divides linearly scheduled `RegionChunk` boundaries across an active set of `std::thread`,
 * spawning atomic threads parsing discrete BAM boundaries with unique indices without locking.
 * Collects returned processed vectors of variants, reassembles them serially 
 * removing arbitrary split overlaps in `merge_chunk_candidates()`.
 *
 * Handles exception safety seamlessly, ensuring memory cleans up on thread death.
 *
 * @param opts Run flags config struct.
 * @param chunks A pre-split list of genomic ranges to distribute.
 * @param read_support_batches Output buffer to stash parsed variant-to-read links.
 */
CandidateTable collect_chunks_parallel(
    const Options& opts,
    const std::vector<RegionChunk>& chunks,
    std::vector<std::vector<ReadSupportRow>>* read_support_batches) {
    if (chunks.empty()) return CandidateTable{};
    ChunkBatchResult result =
        collect_chunk_batch_parallel(opts, chunks, 0, chunks.size(), read_support_batches != nullptr);
    if (read_support_batches != nullptr) *read_support_batches = std::move(result.read_support_batches);
    return merge_chunk_candidates(result.chunks);
}

/**
 * @brief End-to-end collect-bam-variation driver with streaming output.
 *
 * Groups chunks by `reg_chunk_i`, processes each batch in parallel, merges candidates in memory
 * only within the batch, then appends TSV rows and optional VCF / read-support lines. Does not
 * hold the full genome candidate set in RAM.
 *
 * @param opts Output paths, reference, BAM list, and threading configuration.
 */
void run_collect_bam_variation(const Options& opts) {
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));

    const std::vector<RegionChunk> chunks = load_region_chunks(opts);
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai.get());

    std::ofstream variant_out(opts.output_tsv);
    if (!variant_out) throw std::runtime_error("failed to open output: " + opts.output_tsv);
    write_variants_tsv_header(variant_out);

    std::ofstream vcf_out;
    if (!opts.output_vcf.empty()) {
        vcf_out.open(opts.output_vcf);
        if (!vcf_out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);
        write_variants_vcf_header(vcf_out, opts, header.get());
    }

    std::ofstream read_support_out;
    if (!opts.read_support_tsv.empty()) {
        read_support_out.open(opts.read_support_tsv);
        if (!read_support_out) {
            throw std::runtime_error("failed to open read support output: " + opts.read_support_tsv);
        }
        write_read_support_header(read_support_out);
    }

    size_t n_variants = 0;
    size_t n_read_support_rows = 0;
    size_t batch_begin = 0;
    while (batch_begin < chunks.size()) {
        size_t batch_end = batch_begin + 1;
        while (batch_end < chunks.size() &&
               chunks[batch_end].reg_chunk_i == chunks[batch_begin].reg_chunk_i) {
            ++batch_end;
        }

        ChunkBatchResult batch = collect_chunk_batch_parallel(
            opts, chunks, batch_begin, batch_end, !opts.read_support_tsv.empty());
        CandidateTable variants = merge_chunk_candidates(batch.chunks);
        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty()) {
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        }
        if (!opts.read_support_tsv.empty()) {
            for (const std::vector<ReadSupportRow>& rows : batch.read_support_batches) {
                n_read_support_rows += rows.size();
                write_read_support_rows(read_support_out, header.get(), rows);
            }
        }

        batch_begin = batch_end;
    }

    std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
              << " worker thread(s)\n";
    std::cerr << "Collected " << n_variants << " candidate variant sites into "
              << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty()) {
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    }
    if (!opts.read_support_tsv.empty()) {
        std::cerr << "Wrote " << n_read_support_rows << " read x candidate observations to "
                  << opts.read_support_tsv << "\n";
    }
}

} // namespace pgphase_collect

// ════════════════════════════════════════════════════════════════════════════
// CLI entry point
// ════════════════════════════════════════════════════════════════════════════

namespace pgphase_collect {

enum LongOption {
    kMinAltDepthOption = 1000,
    kMinAfOption,
    kMaxAfOption,
    kReadSupportOption,
    kNoisyRegMergeDisOption,
    kMinSvLenOption,
    kChunkSizeOption,
    kHifiOption,
    kOntOption,
    kShortReadsOption,
    kStrandBiasPvalOption,
    kNoisyMaxXgapsOption,
    kMinNoisyRegTotalDepthOption,
    kMaxNoisyFracOption,
    kNoisySlideWinOption,
    kDebugSiteOption,
    kInputIsListOption
};

/**
 * @brief Reads a newline-separated list of BAM/CRAM paths.
 *
 * Strips surrounding whitespace, skips blank lines and `#` comments. Requires at least one path.
 *
 * @param path Text file passed with `--input-is-list` / `-L`.
 * @return Non-empty list of alignment file paths.
 */
static std::vector<std::string> load_bam_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open BAM/CRAM list: " + path);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(in, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;
        const size_t last = line.find_last_not_of(" \t\r\n");
        files.push_back(line.substr(first, last - first + 1));
    }
    if (files.empty()) throw std::runtime_error("BAM/CRAM list is empty: " + path);
    return files;
}

/**
 * @brief Prints usage and option summary for `collect-bam-variation` to stdout.
 */
static void print_collect_help() {
    std::cout
        << "Usage: pgphase collect-bam-variation [options] <ref.fa> <input.bam|bam.list> [region ...]\n"
        << "Options:\n"
        << "  -L, --input-is-list          Treat input path as a list of BAM/CRAM files\n"
        << "  -X, --extra-bam FILE         Extra input BAM/CRAM file; may be repeated\n"
        << "  -t, --threads INT             Region worker threads [1]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
        << "  -B, --min-bq INT              Minimum base quality for candidate sites [10]\n"
        << "  -D, --min-depth INT           Minimum total depth for clean candidates [5]\n"
        << "      --min-alt-depth INT       Minimum alternate depth for clean candidates [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction for clean het candidates [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction for clean het candidates [0.80]\n"
        << "  -r, --region STR              Optional region; may be repeated\n"
        << "      --region-file FILE        BED file of regions to process\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "  -j, --max-var-ratio FLOAT     Skip reads above this variant/ref-span ratio [0.05]\n"
        << "      --max-noisy-frac FLOAT    Skip reads with > this fraction in noisy regions [0.5]\n"
        << "      --include-filtered        Include QC-fail and duplicate reads\n"
        << "  -o, --output FILE             Output TSV file [output.tsv]\n"
        << "  -v, --vcf-output FILE         Optional VCF output for collected candidates\n"
        << "      --read-support FILE       Per-read ref/alt observations at candidates (for phasing)\n"
        << "      --chunk-size INT          Region chunk size in bp [500000]\n"
        << "      --noisy-merge-dis INT     Max distance (bp) to merge noisy/SV windows [500]\n"
        << "      --min-sv-len INT          min_sv_len for noisy-region cgranges merge [30]\n"
        << "      --noisy-slide-win INT     Slide window (bp) for per-read noisy regions [HiFi 100 / ONT/short reads 25]\n"
        << "      --hifi                    HiFi mode: 100 bp noisy window [default; no ONT Fisher strand test]\n"
        << "      --ont                     ONT mode: 25 bp window + Fisher exact test for alt strand bias\n"
        << "      --short-reads             Short-read mode: 25 bp noisy window (no ONT Fisher strand test)\n"
        << "      --strand-bias-pval FLOAT  max p-value for ONT strand filter [0.01]\n"
        << "      --noisy-max-xgaps INT     max indel len (bp) for STR/homopolymer flags [5]\n"
        << "      --min-noisy-reg-total-depth INT   Min reads overlapping a noisy region to keep it [0 off; 5 = paper]\n"
        << "\n"
        << "Regions can also be supplied after <input.bam|bam.list>, e.g.\n"
        << "  pgphase collect-bam-variation ref.fa hifi.bam chr11:1000-2000 chr12:1-500\n";
}

} // namespace pgphase_collect

/**
 * @brief CLI entry for the `collect-bam-variation` subcommand.
 *
 * Parses GNU long options into `Options`, validates numeric thresholds, resolves the BAM list,
 * then calls `run_collect_bam_variation`. Expects `argv` with the `collect-bam-variation` token
 * already removed by the caller.
 *
 * @param argc Argument count.
 * @param argv Argument vector (reference FASTA, input BAM or list, optional region strings).
 * @return 0 on success, 1 on usage or validation error, 1 if `run_collect_bam_variation` throws.
 */
int collect_bam_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;
    std::vector<std::string> extra_bam_files;
    optind = 1;

    const struct option long_options[] = {
        {"threads",                   required_argument, nullptr, 't'},
        {"min-mapq",                  required_argument, nullptr, 'q'},
        {"min-bq",                    required_argument, nullptr, 'B'},
        {"min-depth",                 required_argument, nullptr, 'D'},
        {"min-alt-depth",             required_argument, nullptr, kMinAltDepthOption},
        {"min-af",                    required_argument, nullptr, kMinAfOption},
        {"max-af",                    required_argument, nullptr, kMaxAfOption},
        {"region",                    required_argument, nullptr, 'r'},
        {"region-file",               required_argument, nullptr, 'R'},
        {"autosome",                  no_argument,       nullptr, 'a'},
        {"max-var-ratio",             required_argument, nullptr, 'j'},
        {"max-noisy-frac",            required_argument, nullptr, kMaxNoisyFracOption},
        {"include-filtered",          no_argument,       nullptr, 'f'},
        {"output",                    required_argument, nullptr, 'o'},
        {"vcf-output",                required_argument, nullptr, 'v'},
        {"read-support",              required_argument, nullptr, kReadSupportOption},
        {"chunk-size",                required_argument, nullptr, kChunkSizeOption},
        {"noisy-merge-dis",           required_argument, nullptr, kNoisyRegMergeDisOption},
        {"min-sv-len",                required_argument, nullptr, kMinSvLenOption},
        {"noisy-slide-win",           required_argument, nullptr, kNoisySlideWinOption},
        {"debug-site",                required_argument, nullptr, kDebugSiteOption},
        {"extra-bam",                 required_argument, nullptr, 'X'},
        {"input-is-list",             no_argument,       nullptr, 'L'},
        {"hifi",                      no_argument,       nullptr, kHifiOption},
        {"ont",                       no_argument,       nullptr, kOntOption},
        {"short-reads",               no_argument,       nullptr, kShortReadsOption},
        {"strand-bias-pval",          required_argument, nullptr, kStrandBiasPvalOption},
        {"noisy-max-xgaps",           required_argument, nullptr, kNoisyMaxXgapsOption},
        {"min-noisy-reg-total-depth", required_argument, nullptr, kMinNoisyRegTotalDepthOption},
        {"help",                      no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    bool read_technology_was_set = false;
    bool read_technology_conflict = false;
    const auto set_read_technology = [&](ReadTechnology tech) {
        if (read_technology_was_set && opts.read_technology != tech) {
            read_technology_conflict = true;
            return;
        }
        opts.read_technology = tech;
        read_technology_was_set = true;
    };

    while ((opt = getopt_long(argc, argv, "t:q:B:D:r:R:aj:o:v:hX:L", long_options, &long_index)) != -1) {
        switch (opt) {
            case 't': opts.threads = std::stoi(optarg); break;
            case 'q': opts.min_mapq = std::stoi(optarg); break;
            case 'B': opts.min_bq = std::stoi(optarg); break;
            case 'D': opts.min_depth = std::stoi(optarg); break;
            case kMinAltDepthOption:    opts.min_alt_depth = std::stoi(optarg); break;
            case kMinAfOption:          opts.min_af = std::stod(optarg); break;
            case kMaxAfOption:          opts.max_af = std::stod(optarg); break;
            case 'r': opts.regions.push_back(optarg); break;
            case 'R': opts.region_file = optarg; break;
            case 'a': opts.autosome = true; break;
            case 'j': opts.max_var_ratio_per_read = std::stod(optarg); break;
            case kMaxNoisyFracOption:   opts.max_noisy_frac_per_read = std::stod(optarg); break;
            case 'f': opts.include_filtered = true; break;
            case 'o': opts.output_tsv = optarg; break;
            case 'v': opts.output_vcf = optarg; break;
            case kReadSupportOption:    opts.read_support_tsv = optarg; break;
            case kChunkSizeOption:      opts.chunk_size = std::stoll(optarg); break;
            case kNoisyRegMergeDisOption: opts.noisy_reg_merge_dis = std::stoi(optarg); break;
            case kMinSvLenOption:       opts.min_sv_len = std::stoi(optarg); break;
            case kNoisySlideWinOption:  opts.noisy_reg_slide_win = std::stoi(optarg); break;
            case kDebugSiteOption:      opts.debug_site = optarg; break;
            case 'X': extra_bam_files.push_back(optarg); break;
            case 'L': opts.input_is_list = true; break;
            case kHifiOption:           set_read_technology(ReadTechnology::Hifi); break;
            case kOntOption:            set_read_technology(ReadTechnology::Ont); break;
            case kShortReadsOption:     set_read_technology(ReadTechnology::ShortReads); break;
            case kStrandBiasPvalOption: opts.strand_bias_pval = std::stod(optarg); break;
            case kNoisyMaxXgapsOption:  opts.noisy_reg_max_xgaps = std::stoi(optarg); break;
            case kMinNoisyRegTotalDepthOption:
                opts.min_noisy_reg_total_depth = std::stoi(optarg);
                break;
            case 'h': print_collect_help(); return 0;
            default:  print_collect_help(); return 1;
        }
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_bq < 0 || opts.chunk_size < 1 ||
        opts.min_depth < 0 || opts.min_alt_depth < 0 || opts.noisy_reg_merge_dis < 0 ||
        opts.min_sv_len < 0 || opts.min_af < 0.0 || opts.max_af < opts.min_af ||
        opts.strand_bias_pval < 0.0 || opts.strand_bias_pval > 1.0 ||
        opts.max_var_ratio_per_read < 0.0 || opts.max_noisy_frac_per_read < 0.0 ||
        opts.noisy_reg_max_xgaps < 0 || opts.min_noisy_reg_total_depth < 0 ||
        opts.noisy_reg_slide_win < -1) {
        std::cerr << "Error: numeric thresholds are invalid\n";
        return 1;
    }
    if (read_technology_conflict) {
        std::cerr << "Error: choose only one of --hifi, --ont, or --short-reads\n";
        return 1;
    }
    if (optind + 2 > argc) {
        print_collect_help();
        return 1;
    }

    opts.ref_fasta = argv[optind];
    const std::string input_path = argv[optind + 1];
    if (opts.input_is_list) {
        opts.bam_files = load_bam_list(input_path);
    } else {
        opts.bam_files.push_back(input_path);
    }
    opts.bam_files.insert(opts.bam_files.end(), extra_bam_files.begin(), extra_bam_files.end());
    opts.bam_file = opts.bam_files.front();
    for (int arg_i = optind + 2; arg_i < argc; ++arg_i) {
        opts.regions.push_back(argv[arg_i]);
    }

    try {
        run_collect_bam_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
