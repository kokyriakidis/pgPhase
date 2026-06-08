#ifndef PGPHASE_GRAPH_QUERY_HPP
#define PGPHASE_GRAPH_QUERY_HPP

#include "collect_types.hpp"
#include "graph_sites.hpp"

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <string>
#include <vector>

namespace pgphase_collect {

struct GraphReadAllele {
    std::string site_id;
    std::string chrom;
    hts_pos_t pos = 0;
    std::string read_name;
    int allele = kGraphAlleleMissing;
    int mapq = 255;
    bool reverse = false; // true when the read traverses the snarl in reverse-complement
};

struct GraphQueryConfig {
    std::string gbz_db;
    std::string gaf_db;
    int min_mapq = kDefaultMinMapq;
};

// RAII handle for a pggaf coordinate-indexed GAF file (bgzipped + tabix).
// Opened once per thread, reused across chunks to avoid repeated file
// open / tabix index load overhead.
struct IndexedGafHandle {
    htsFile* fp  = nullptr;
    tbx_t*   tbx = nullptr;

    IndexedGafHandle() = default;
    explicit IndexedGafHandle(const std::string& path);
    ~IndexedGafHandle();
    IndexedGafHandle(const IndexedGafHandle&) = delete;
    IndexedGafHandle& operator=(const IndexedGafHandle&) = delete;
    IndexedGafHandle(IndexedGafHandle&& o) noexcept;
    IndexedGafHandle& operator=(IndexedGafHandle&& o) noexcept;
    explicit operator bool() const { return fp && tbx; }
};

// Per-chunk query using a pre-opened IndexedGafHandle.
// The GAF must have pggaf's three leading coordinate columns (rc rb re).
std::vector<GraphReadAllele>
scan_indexed_gaf_chunk(IndexedGafHandle& handle,
                       const std::string& contig,
                       hts_pos_t beg,
                       hts_pos_t end,
                       const GraphSiteCatalogView& catalog,
                       int min_mapq);

// Validates that the GAF file is bgzip-compressed and has a .tbi index.
void require_indexed_gaf(const std::string& indexed_gaf_file);

// FFI-based interval query: uses pre-opened GBZ/GAF handles instead of
// spawning a subprocess per chunk.
std::vector<GraphReadAllele>
query_gbz_interval_gaf_ffi(void* gbz_handle,
                            void* gaf_handle,
                            const std::string& sample,
                            const std::string& contig,
                            hts_pos_t beg,
                            hts_pos_t end,
                            const GraphSiteCatalogView& catalog,
                            int min_mapq);

// Genome-coordinate range covered by a reference path in the GBZ.
// For subgraphed GBZ files the path may start at a non-zero offset.
// `end` is a conservative lower bound (the true path may extend up to
// ~1000 bp beyond).
struct PathRange {
    hts_pos_t start = 0;
    hts_pos_t end   = 0;
    bool valid      = false;
};

// Query the genome-coordinate range of a reference path.
PathRange query_gbz_path_range(void* gbz_handle,
                               const std::string& sample,
                               const std::string& contig);

} // namespace pgphase_collect

#endif
