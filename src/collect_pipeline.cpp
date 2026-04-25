#include "collect_pipeline.hpp"

#include "bam_digar.hpp"
#include "candidate_collection.hpp"
#include "collect_output.hpp"
#include "collect_regions.hpp"
#include "noisy_regions.hpp"
#include "variant_classification.hpp"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace pgphase_collect {

static void load_and_prepare_chunk(BamChunk& chunk, const Options& opts, WorkerContext& context) {
    chunk.reads.clear();
    for (size_t input_i = 0; input_i < context.bams.size(); ++input_i) {
        std::vector<ReadRecord> reads = load_read_records_for_chunk(
            opts,
            chunk.region,
            *context.bams[input_i],
            context.headers[input_i].get(),
            context.indexes[input_i].get(),
            context.ref);
        chunk.reads.insert(
            chunk.reads.end(),
            std::make_move_iterator(reads.begin()),
            std::make_move_iterator(reads.end()));
    }
    std::sort(chunk.reads.begin(), chunk.reads.end(), [](const ReadRecord& lhs, const ReadRecord& rhs) {
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        if (lhs.end != rhs.end) return lhs.end < rhs.end;
        if (lhs.nm != rhs.nm) return lhs.nm < rhs.nm;
        return lhs.qname < rhs.qname;
    });
    finalize_bam_chunk(chunk, context.ref, context.primary_header());
}

static void collect_prephase_candidates(BamChunk& chunk,
                                        const Options& opts,
                                        std::vector<ReadSupportRow>* read_support_out) {
    pre_process_noisy_regs_pgphase(chunk, opts);
    collect_candidate_sites_from_records(chunk.region, chunk.reads, chunk.candidates);
    collect_allele_counts_from_records(
        chunk.reads, chunk.candidates, &chunk.region, read_support_out, opts.min_bq);
}

static void classify_and_filter_candidates(BamChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    classify_chunk_candidates(chunk, opts, header);
    post_process_noisy_regs_pgphase(chunk, chunk.candidates);
    if (!opts.is_ont()) {
        apply_noisy_containment_filter(chunk);
    }
}

static CandidateTable process_chunk(const RegionChunk& region,
                                    const Options& opts,
                                    WorkerContext& context,
                                    std::vector<ReadSupportRow>* read_support_out) {
    BamChunk chunk;
    chunk.region = region;
    load_and_prepare_chunk(chunk, opts, context);
    collect_prephase_candidates(chunk, opts, read_support_out);
    classify_and_filter_candidates(chunk, opts, context.primary_header());
    return std::move(chunk.candidates);
}

static CandidateTable merge_chunk_candidates(std::vector<CandidateTable>& chunk_tables) {
    CandidateTable merged;
    for (CandidateTable& table : chunk_tables) {
        merged.insert(
            merged.end(),
            std::make_move_iterator(table.begin()),
            std::make_move_iterator(table.end()));
    }
    collapse_fuzzy_large_insertions(merged);
    return merged;
}

CandidateTable collect_chunks_parallel(const Options& opts,
                                       const std::vector<RegionChunk>& chunks,
                                       std::vector<std::vector<ReadSupportRow>>* read_support_batches) {
    std::vector<CandidateTable> chunk_tables(chunks.size());
    if (chunks.empty()) return CandidateTable{};

    if (read_support_batches != nullptr) read_support_batches->assign(chunks.size(), {});

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), chunks.size());
    std::atomic<size_t> next_chunk{0};
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
                    const size_t chunk_i = next_chunk.fetch_add(1);
                    if (chunk_i >= chunks.size()) break;

                    std::vector<ReadSupportRow>* rs_ptr = nullptr;
                    if (read_support_batches != nullptr) rs_ptr = &(*read_support_batches)[chunk_i];
                    chunk_tables[chunk_i] = process_chunk(chunks[chunk_i], opts, context, rs_ptr);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }

    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);

    return merge_chunk_candidates(chunk_tables);
}

void run_collect_bam_variation(const Options& opts) {
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));

    const std::vector<RegionChunk> chunks = load_region_chunks(opts);
    std::vector<std::vector<ReadSupportRow>> read_support_batches;
    std::vector<std::vector<ReadSupportRow>>* rs_batches_ptr =
        opts.read_support_tsv.empty() ? nullptr : &read_support_batches;
    CandidateTable variants = collect_chunks_parallel(opts, chunks, rs_batches_ptr);
    write_variants(opts, fai.get(), variants);
    write_variants_vcf(opts, fai.get(), variants);
    if (!opts.read_support_tsv.empty()) {
        write_read_support_tsv(opts, read_support_batches);
    }

    std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
              << " worker thread(s)\n";
    std::cerr << "Collected " << variants.size() << " candidate variant sites into "
              << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty()) {
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    }
    if (!opts.read_support_tsv.empty()) {
        size_t n_rows = 0;
        for (const auto& b : read_support_batches) n_rows += b.size();
        std::cerr << "Wrote " << n_rows << " read x candidate observations to " << opts.read_support_tsv
                  << "\n";
    }
}

} // namespace pgphase_collect
