# `phase-graph` Implementation Description

This document describes the implementation of the `pgphase phase-graph` command from command-line
input to final phased output. The command takes a decomposed graph-site VCF (produced by
`build-snarl-catalog` / `vg deconstruct`) and a raw GAF alignment file and produces per-read
haplotype assignments and phased graph-site calls without requiring a reference FASTA or coordinate-
sorted BAM. It is the graph-native complement to `collect-bam-variation`: instead of calling
candidates from BAM digars and phasing them, it reads pre-called snarl genotype candidates from the
VCF, matches raw graph-walk alignments to those candidates, and runs the same k-means + chunk-
stitching engine used by the BAM path.

The implementation spans the following source files:

```text
graph_pipeline.cpp    CLI parser, contig discovery, per-contig streaming orchestration,
                      shard infrastructure, parallel GAF scan, parallel chunk build loop,
                      streaming output coordination, final output writers.

graph_bam_adapter.cpp Build GraphReadAllele rows into a BamChunk (allele filtering, depth/AF
                      gating, read-profile construction); per-chunk hap assignment; overlap
                      population; chunk-pair stitching; per-chunk and end-of-run output writers.

graph_bam_adapter.hpp Public API for the adapter: GraphBamChunkBuildResult, PhaseReadOutputRow,
                      FilteredGraphSite; streaming build/assign/stitch/merge/write functions;
                      header writers and row writers for each output file.

graph_sites.hpp/.cpp  GraphSite / GraphSiteCatalog data types; VCF catalog loader
                      (load_graph_site_catalog_from_vcf, with tabix fast-path and streaming
                      fallback); contig name discovery (load_graph_site_contig_names);
                      walk parsing and string round-trip helpers.

graph_query.hpp/.cpp  GraphReadAllele type; compact site index (CompactGraphSite,
                      CompactGraphSiteIndex); node-handle encoding; GAF streaming and walk-
                      match logic; scan_gaf_for_catalog_emit_parallel_releasing_walks; per-
                      chunk and diagnostic row writers.

collect_phase.cpp     Shared k-means hap assignment (assign_hap_based_on_germline_het_vars_kmeans)
                      and chunk-boundary stitching (stitch_chunk_haps / flip_chunk_hap) used
                      identically by both the BAM and graph paths.

collect_types.hpp     Shared Options struct (including touch_read_phase), BamChunk, and other
                      types consumed by collect_phase.cpp.
```

### Reading guide

Sections §1–§14 describe the pipeline in execution order. §15 covers the memory model and the
current I/O tradeoff between RAM and GAF scan count. §16 summarises the full pipeline.

Coordinate conventions:

- Graph VCF `POS` is **1-based** (standard VCF).
- All interval arithmetic inside the pipeline uses **0-based half-open** positions derived from
  `site.pos - 1` (or `site.ref_beg - 1` when `ref_beg` is set).
- `chunk_size` is specified in base-pairs of reference coordinate space.

---

## 1. Command-Line Configuration

```text
pgphase phase-graph [options]
  --gaf FILE                Raw GAF graph alignments (required)
  --graph-sites-vcf FILE    Decomposed graph-site VCF (required)
  --contig NAME             Restrict to one reference contig
  --interval BEG..END       0-based half-open bounds (with --contig only)
  --chunk-size INT          Phasing chunk size in bp [500000]
  --graph-phase-sites FILE  Phased graph-site calls TSV [graph_phase_sites.tsv]
  --graph-phase-reads FILE  Per-read HAP/PHASE_SET TSV [graph_phase_reads.tsv]
  --graph-phase-bam FILE    Unaligned BAM with HP/PS aux tags per read
  --graph-sites-tsv FILE    Diagnostic: parsed graph-site catalog
  --graph-read-support FILE Diagnostic: raw read→site allele observation TSV
  --graph-site-counts FILE  Diagnostic: per-site allele depth counts TSV
  --graph-read-profile FILE Diagnostic: sparse read×site allele profile TSV
  --graph-filtered-sites FILE Diagnostic: sites removed by depth/AF filters
  -t, --threads INT         Worker threads [1]
  -V, --verbose INT         Verbosity level [0]
```

Important defaults and semantics:

- `--chunk-size` follows the same tiling logic as `collect-bam-variation`: the genome is divided
  into non-overlapping windows of this width and each window is processed as one `BamChunk`.
  Adjacent chunks share overlap reads for boundary stitching (§12).
- `--contig` restricts the entire run to a single named contig. Without it, all contigs present
  in the VCF are processed in sequence (§3).
- `--interval` is only valid when `--contig` is also supplied. It is silently ignored for multi-
  contig runs.
- `--graph-phase-sites` and `--graph-phase-reads` both default to non-empty paths, so phased
  output is always produced unless explicitly redirected to `/dev/null`.

---

## 2. Key Data Structures

### 2.1 `GraphSite` (`graph_sites.hpp`)

One decomposed snarl from the VCF:

```cpp
struct GraphSite {
    std::string chrom;         // raw VCF CHROM (may be "CHM13#0#chr20")
    std::string ref_contig;    // normalised contig ("chr20"), derived from CHROM or RC INFO
    hts_pos_t   pos;           // 1-based VCF POS
    std::string id;            // VCF ID (site key used in all shard routing)
    std::string ref;           // REF allele sequence
    std::vector<std::string> alts;
    std::vector<std::string> allele_traversals; // raw AT strings (released after index build)
    std::vector<GraphWalk>   allele_walks;       // parsed walk vectors (released after scan)
    hts_pos_t ref_beg, ref_end;                  // from END= INFO field
    int level;                                   // LV= nesting level
    std::string parent, root;                    // PS=/RS= parent/root snarl IDs
    std::vector<int> conditional_parent_alleles; // PA=
    bool has_spanning_deletion;
    bool eligible;
    std::string skip_reason;
};
```

`GraphWalk` is `std::vector<GraphWalkStep>` where each step holds a node ID string and an
orientation flag. Walk index 0 is always the reference traversal; subsequent indices are alt
traversals.

### 2.2 `GraphReadAllele` (`graph_query.hpp`)

One observation of a read traversing a site:

```cpp
struct GraphReadAllele {
    std::string site_id;    // matches GraphSite::id (or constructed chrom:pos:ref key)
    std::string chrom;
    hts_pos_t   pos;
    std::string read_name;
    int         allele;     // 0 = ref, 1+ = alt index; kGraphAlleleMissing / kGraphAlleleAmbiguous
    std::string walk;       // diagnostic only; empty in shard files
    int         mapq;
};
```

### 2.3 `GraphBamChunkBuildResult` (`graph_bam_adapter.hpp`)

Output of one `build_graph_bam_chunk` call — wraps a fully-constructed `BamChunk` together with
the site metadata needed for downstream writers:

```cpp
struct GraphBamChunkBuildResult {
    BamChunk                    chunk;                // k-means input/output state
    std::vector<std::string>    site_ids;             // sites present in this chunk
    std::vector<std::vector<int>> site_allele_orig_idx; // surviving allele → original allele index
    std::vector<FilteredGraphSite> filtered_sites;    // sites removed by depth/AF gates
};
```

`BamChunk` (from `collect_types.hpp`) is the same struct used by the BAM path: it holds read
profiles, candidate variant rows, hap assignments, phase sets, and overlap bookkeeping. It is
non-copyable (owns a `cgranges_t*` via `unique_ptr`); all moves use `std::move`.

### 2.4 `PhaseReadOutputRow` (`graph_bam_adapter.hpp`)

Compact per-read accumulator maintained across chunks for the whole run:

```cpp
struct PhaseReadOutputRow {
    std::string read_name;
    int         chunk_id = -1;
    int         hap = 0;
    hts_pos_t   phase_set = -1;
    bool        has_phased_assignment = false;
    int         copies = 0;
    int         best_assignment_obs = -1;
    int         assignment_chunk_id = std::numeric_limits<int>::max();
    std::unordered_map<std::string, int> allele_by_site;
};
```

`assignment_chunk_id` is used for conflict resolution when a read appears in multiple chunks:
the assignment from the chunk with the smallest global `chunk_id` (i.e., earlier in the genome)
wins. `best_assignment_obs` tracks the k-means observation count so that a better-supported
assignment in a later chunk can still override a weak earlier one.

---

## 3. Contig Discovery

When `--contig` is supplied the pipeline processes exactly that contig. When it is absent,
`load_graph_site_contig_names` enumerates all distinct reference contigs in the VCF:

```cpp
std::vector<std::string> load_graph_site_contig_names(const std::string& vcf_path);
```

**Tabix fast-path** (bgzipped VCF with `.tbi` or `.csi` index): reads the sequence-name table
from the index with `tbx_seqnames` — no data lines are scanned.

**Streaming fallback** (plain-text or bgzipped VCF without index): streams all data lines,
extracts the first tab-delimited field, and collects distinct values. Result is sorted
lexicographically.

In both paths, pangenome-prefixed contig names (`CHM13#0#chr20`) are normalised to the suffix
after the last `#` (`chr20`). Deduplication preserves VCF-header order for the tabix path and
sorts for the fallback path.

The returned list drives the outer contig loop in `phase_graph()`. Each contig is processed
independently with no shared mutable state between iterations; `rows_by_read` and
`diagnostic_catalog` are the only accumulators that grow across iterations.

---

## 4. Catalog Loading (Per-Contig)

For each contig, the pipeline loads only that contig's sites:

```cpp
GraphSiteCatalog load_graph_site_catalog_from_vcf(
    const std::string& path,
    const std::vector<RegionFilter>& filters = {},
    bool keep_allele_traversal_strings = true);
```

A `RegionFilter` specifying the current contig name is passed. When a tabix index is present the
loader queries only the relevant sequence from the index; without an index it streams all data
lines and skips non-matching `CHROM` values. Both paths converge to a stable-sorted catalog
(`chrom`, then `pos`, then `id`).

After loading, site metadata is re-validated:

- `graph_site_validation_skip_reason` checks that the reference walk (allele 0) has at least two
  nodes and that no structural problems (missing ID, missing walks, mismatched allele counts) are
  present. Sites that fail receive `eligible = false` and a `skip_reason` string.
- Spanning-deletion sites (`ALT = *`) set `has_spanning_deletion = true` and are excluded from
  phasing (not placed in any chunk catalog).

The catalog retains allele walk strings (`allele_traversals`, `allele_walks`) at this point
because they are needed to build the compact index for the GAF scan (§7). They are released
during or after the scan to reclaim memory.

---

## 5. Chunk Interval Splitting

```cpp
std::vector<std::pair<hts_pos_t, hts_pos_t>>
split_graph_interval(hts_pos_t beg, hts_pos_t end, hts_pos_t chunk_size);
```

Returns a vector of contiguous half-open intervals `[beg, beg+chunk_size), [beg+chunk_size,
beg+2×chunk_size), …` covering `[beg, end)`. The last interval is clipped to `end`.

`end` is computed as `catalog_end_bound(full_catalog)` — the maximum `ref_end` (or `pos`) across
all sites in the current contig's catalog — unless overridden by `--interval`.

Each interval becomes one chunk index. The relationship between a site and its chunk is:

```text
site_beg0 = (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1   // 0-based
ci        = (site_beg0 - interval_beg) / chunk_size
```

This is an O(1) arithmetic assignment — no scan over intervals is needed.

---

## 6. Per-Chunk Catalog Construction

After interval splitting, a single O(n_sites) pass routes each site to its chunk:

```cpp
std::vector<GraphSiteCatalog> chunk_catalogs(n);  // one per interval
std::unordered_map<std::string, size_t> site_to_chunk;

for (const GraphSite& site : full_catalog.sites) {
    const hts_pos_t site_beg0 = ...;
    const size_t ci = (site_beg0 - interval_beg) / chunk_size;
    chunk_catalogs[ci].sites.push_back(site);
    site_to_chunk.emplace(site_key(site), ci);
}
```

`site_key` returns `site.id` when the VCF ID is not `.`, otherwise constructs
`"chrom:pos:ref"`. This key is used identically during the GAF scan to route matched
`GraphReadAllele` rows to the correct chunk's shard file.

After the site-to-chunk assignment, **walk strings are released from `chunk_catalogs`**:

```cpp
for (auto& cc : chunk_catalogs)
    for (GraphSite& s : cc.sites) {
        for (GraphWalk& w : s.allele_walks) w.clear();
        s.allele_traversals.clear();
    }
```

The chunk catalogs now retain only site metadata (id, chrom, pos, ref, alts, allele counts).
This reduces RAM before the parallel build phase. The walk strings are still present in
`full_catalog` at this point and remain there until the GAF scan releases them (§7).

---

## 7. Shard Infrastructure

Each contig run creates a fresh per-(worker × chunk) temporary directory:

```text
/tmp/pgphase_graph_shard_<pid>_<counter>/
    w0_c0.tsv   w0_c1.tsv   …   w0_cN.tsv
    w1_c0.tsv   w1_c1.tsv   …   w1_cN.tsv
    …
    wW_c0.tsv   …           …   wW_cN.tsv
```

`ChunkRowShard` is an RAII wrapper: its destructor removes all shard files and calls `rmdir` on
the directory.

`ChunkRowWriter` is an LRU-caching append writer that keeps up to 64 file descriptors open
simultaneously. When a worker writes to chunk index `ci`, the writer appends one TSV row
(site_id, chrom, pos, read_name, mapq, allele) to `paths[worker_id][ci]`. Workers never share
a writer, so no locking is needed at the file level.

`ChunkRowShard::load_chunk(ci)` reconstitutes a `std::vector<GraphReadAllele>` from all
workers' shard files for chunk `ci`. This is called once per chunk during the parallel build
phase (§9).

The shard mechanism decouples the GAF scan (streaming write, all threads) from the chunk build
(random-access read, all threads) without holding all rows in RAM simultaneously. Peak RAM from
shard rows = `n_workers × one_chunk_rows` rather than all rows across all chunks.

---

## 8. GAF Scanning and Walk Matching

### 8.1 Compact Site Index

Before scanning the GAF, `scan_gaf_for_catalog_emit_parallel_releasing_walks` builds a
`CompactGraphSiteIndex` from the contig's `full_catalog`:

```cpp
struct CompactGraphSite {
    size_t   original_index;   // back-reference into full_catalog.sites
    std::string site_id;
    std::string chrom;
    hts_pos_t   pos;
    CompactHandle left, right;           // ref-walk boundary node handles (forward)
    CompactHandle left_rev, right_rev;  // boundary node handles (reverse complement)
    std::vector<std::vector<CompactHandle>> allele_walks;  // all alleles, all nodes
};

struct CompactGraphSiteIndex {
    std::vector<CompactGraphSite> sites;
    std::unordered_map<CompactHandle, std::vector<size_t>> boundary_to_sites;
};
```

`CompactHandle` is a `uint32_t` encoding: `(node_id << 1) | orientation_bit`. Encoding with
`graph_walk_to_compact_handles` converts all walk steps. The `boundary_to_sites` map indexes
every snarl's left and right boundary node (both orientations) to the set of compact sites that
share that boundary. This is the primary lookup used during GAF scanning.

Only sites for which `graph_site_is_queryable` returns true and whose reference walk has at
least two nodes are included in the compact index. Sites with malformed or empty walks are
silently skipped.

After building the compact index the function optionally releases allele walk strings from
`full_catalog` (the `_releasing_walks` variant used here). This frees the last remaining walk
string storage from the per-contig catalog, since the compact index has its own `uint32`
representation.

### 8.2 Parallel GAF Streaming

The GAF file is scanned with `n_workers` threads. Each thread reads a disjoint range of lines:
the file is divided into roughly equal byte ranges and each thread seeks to the start of the
first complete line within its range.

For each GAF line, `parse_gaf_core_fields` extracts three fields without heap allocation using
`string_view`:

```text
field  0 → read_name
field  5 → path (node walk, e.g. ">114818865>114818866>114818867")
field 11 → mapq
```

Lines with `mapq < min_mapq` are skipped. Lines whose walk cannot be parsed as compact handles
fall through to a string-walk path.

### 8.3 Walk Matching

For each parsed read walk, the scanner searches for sites whose boundary nodes appear in the
read's path:

1. Iterate every handle in the read's compact walk.
2. Look up the handle in `boundary_to_sites`. For each candidate site:
   a. Find the position of the site's `left` boundary in the read walk.
   b. Find the position of the site's `right` boundary at or after `left`.
   c. Extract the sub-walk `[left_pos, right_pos]` from the read.
   d. Compare the sub-walk against each allele walk with `compact_span_matches_allele`.
   e. Emit a `GraphReadAllele` with the matching allele index (0 = ref, 1+ = alt).
   f. If no allele matches, emit `kGraphAlleleAmbiguous`; if the span is not found, skip.

A read can match multiple sites (it traverses multiple snarls). Each match is emitted as a
separate `GraphReadAllele` observation.

### 8.4 Routing to Shard Files

The emitter callback in `phase_graph()` receives `(worker_id, GraphReadAllele&&)` and routes
the row using `site_to_chunk`:

```cpp
auto it = site_to_chunk.find(row.site_id);
if (it != site_to_chunk.end())
    writers[worker_id]->write(it->second, row);
```

`site_to_chunk` was built in §6 from the current contig's sites. Observations for unknown site
IDs (sites outside the coordinate window or already-released sites) are dropped silently.

After all threads complete, `writers.clear()` flushes and closes all shard files.

---

## 9. Parallel Chunk Build

`run_vcf_gaf_chunks_parallel` processes a batch of `[batch_beg, batch_end)` chunk indices in
parallel. Each worker:

1. Calls `shard.load_chunk(i)` — reads all workers' TSV rows for chunk `i` from disk.
2. Calls `build_graph_bam_chunk(chunk_catalogs[i], chunk_rows, contig, beg, end, chunk_id, opts)`.

### 9.1 `build_graph_bam_chunk`

This function constructs a `BamChunk` from the raw `GraphReadAllele` observations for one chunk.

**Phase 1 — Alt-only allele filtering**: sites where all observations are reference allele
(`allele == 0`) are dropped. For sites that survive, only allele indices that were actually
observed in at least one read are retained. The `site_allele_orig_idx` mapping records the
original VCF allele index for each surviving allele position.

**Phase 2 — Depth and AF filtering**: for each remaining site, counts of reference and
non-reference observations are tallied. Sites that fail any of:
- total depth < `min_depth`
- allele fraction (non-ref / total) < `min_af`
- allele fraction > `max_af`

are moved to `filtered_sites` (a `FilteredGraphSite` record with the failure reason) and
excluded from phasing.

**Read profile construction**: surviving observations are used to populate
`BamChunk::read_var_profiles` — the same `ReadVariantProfile` / `hap_to_alle_profile` struct
used by the BAM path. Each read gets one entry per site it was observed at, with the allele
index stored as the profile value. The chunk's candidate variant rows (`BamChunk::cand_vars`)
are populated from the surviving site list.

The resulting `BamChunk` is indistinguishable to the phasing engine from one produced by the
BAM path.

---

## 10. Per-Chunk Hap Assignment

```cpp
void assign_graph_chunk_hap(GraphBamChunkBuildResult& gc, const Options& opts);
```

This is a thin wrapper around `phase_graph_bam_chunks` restricted to a single chunk. It calls
`assign_hap_based_on_germline_het_vars_kmeans` (from `collect_phase.cpp`) on the chunk's
`BamChunk`, runs the same k-means read-clustering algorithm used by the BAM path (§18 of the
BAM path doc), and stores hap and phase-set assignments into the chunk's read profiles and
candidate variants.

No cross-chunk communication happens here. Stitching is a separate step (§12).

---

## 11. Build-Time Streaming Outputs

Three diagnostic outputs are written immediately after `assign_graph_chunk_hap`, before any
stitching. They reflect pre-stitch per-chunk state.

### 11.1 `--graph-site-counts` TSV

Columns: `CHUNK_ID SITE_ID CHROM POS REF_COV ALT_COV TOTAL_COV AF`

One row per surviving site per chunk. Depths come from the `BamChunk` candidate allele counts
populated in §9.1.

### 11.2 `--graph-read-profile` TSV

Columns: `CHUNK_ID READ SITE_ID ALLELE`

One row per (read, site) observation that survived both filtering phases in §9.1. This is the
sparse read×site matrix fed to k-means.

### 11.3 `--graph-filtered-sites` TSV

Columns: `CHUNK_ID SITE_ID CHROM POS REF_COV ALT_COV TOTAL_COV AF FILTER_REASON`

One row per site removed in Phase 2 depth/AF filtering (`ref_only`, `low_depth`, `high_af`,
`low_af`).

These three files have their headers written once before the contig loop begins.

---

## 12. Chunk-Pair Overlap and Stitching

### 12.1 Overlap Population

```cpp
void populate_graph_chunk_pair_overlaps(GraphBamChunkBuildResult& pre,
                                        GraphBamChunkBuildResult& cur);
```

Builds the `BamChunk::down_ovlp_read_i` / `up_ovlp_read_i` overlap lists expected by
`stitch_chunk_haps`. For the graph path, two adjacent chunks share reads that were assigned
observations in both windows (a read that spans a chunk boundary will appear in both chunks'
`GraphReadAllele` rows, which become entries in both `BamChunk::read_var_profiles`).

The function builds a name-to-index map for each chunk, then intersects them:

```cpp
const auto pre_reads = read_index_by_name(pre.chunk);
const auto cur_reads = read_index_by_name(cur.chunk);
for (const auto& [name, cur_idx] : cur_reads) {
    auto it = pre_reads.find(name);
    if (it == pre_reads.end()) continue;
    pre.chunk.down_ovlp_read_i[0].push_back(it->second);
    cur.chunk.up_ovlp_read_i[0].push_back(cur_idx);
}
```

The overlap lists use a single slot (index 0) because each chunk-pair has exactly one boundary.
This mirrors the BAM path's `build_bam_chunk_overlaps` but without coordinate-based positional
logic — shard membership determines overlap, not genomic window containment.

### 12.2 `stitch_graph_chunk_pair`

```cpp
void stitch_graph_chunk_pair(GraphBamChunkBuildResult& pre,
                             GraphBamChunkBuildResult& cur,
                             const Options& opts);
```

1. Calls `populate_graph_chunk_pair_overlaps(pre, cur)`.
2. Moves `pre.chunk` and `cur.chunk` into a two-element `std::vector<BamChunk>` (using
   `push_back(std::move(...))` — not an initializer list, since BamChunk is non-copyable).
3. Calls `stitch_chunk_haps(pair, &opts, nullptr)` — the same stitching function used by the
   BAM path.
4. Moves the chunks back.

`stitch_chunk_haps` applies `flip_chunk_hap`: if the common-read vote across the boundary
disagrees with the current hap assignment of `cur`, it flips `cur`'s read hap assignments,
candidate hap assignments, and phase-set polarities.

`opts.touch_read_phase` must be `true` for flip to propagate to per-read hap values (not just
candidate variants). The pipeline sets this unconditionally:

```cpp
filter_opts.touch_read_phase = true;
```

This ensures that `merge_graph_chunk_into_read_rows` (§13) sees correctly stitched per-read
hap assignments after every stitch.

### 12.3 Stitch-Time Output: `--graph-phase-sites` TSV

Columns: `CHROM POS TYPE REF ALT DP REF_COV ALT_COV AF PHASE_SET HAP_ALT HAP_REF`

Written for `prev_chunk` immediately after `stitch_graph_chunk_pair` finalises it. REF and ALT
sequences come from the `site_map` (site_id → `const GraphSite*` built from `full_catalog`
before the streaming loop). The header is written once before the contig loop.

---

## 13. Per-Read Accumulator

```cpp
void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphBamChunkBuildResult& gc);
```

After a chunk is fully stitched and finalized it is merged into the global `rows_by_read` map.
For each read in the chunk's `BamChunk::read_var_profiles`:

- If the read has not been seen before: insert a new `PhaseReadOutputRow`.
- If it has been seen: update if the new chunk provides a better assignment. "Better" is defined
  as: lower `assignment_chunk_id` (earlier in the genome), or equal `assignment_chunk_id` and
  higher `best_assignment_obs` (more k-means evidence).

The map grows monotonically across all chunks and all contigs. It is the single source of truth
for end-of-run phase BAM and reads TSV output.

---

## 14. End-of-Run Outputs

### 14.1 `--graph-sites-tsv` Diagnostic

The `diagnostic_catalog` accumulates one `GraphSiteCatalog` entry per chunk throughout the run
(via `append_catalog`). After all contigs complete, the catalog is serialised by
`write_graph_site_catalog_tsv`. This is a tab-delimited dump of every site that was processed,
regardless of whether it survived depth/AF filtering, useful for inspecting the raw snarl
catalog in tabular form.

### 14.2 `--graph-read-support` Diagnostic

Written per-contig by re-scanning that contig's shard files while they are still alive (before
the shard RAII destructor removes them):

```cpp
if (read_support_out) {
    for (size_t ci = 0; ci < n_shards; ++ci)
        write_graph_read_alleles_tsv_rows(read_support_out, shard.load_chunk(ci));
}
```

The header is written once before the contig loop. This is the raw (pre-filtering)
read→site→allele observation table produced by the GAF scan.

### 14.3 `--graph-phase-reads` TSV and `--graph-phase-bam`

Both require a complete view of all reads across all contigs, so they are written after the
contig loop from the accumulated `rows_by_read` map.

A second GAF pass (`collect_gaf_read_names`) collects all MAPQ-passing read names so that reads
with no site observations appear in the output as unphased rows (HAP=0, PHASE_SET=-1). This
second pass reads only fields 0 and 11 (read name and MAPQ) from each GAF line.

**Phase reads TSV** columns:
`READ HAP PHASE_SET COPIES BEST_OBS CHUNK_ID ASSIGNMENT_CHUNK_ID ALLELE_BY_SITE`

**Phase BAM** (`write_graph_bam_phase_bam_from_rows`): produces an unaligned name-sorted BAM.
For each read in `rows_by_read` (or in `all_read_names`):
- If `has_phased_assignment`: set `HP` aux tag (1 or 2) and `PS` aux tag (phase-set value).
- If unphased: record is written with no HP/PS tags.

---

## 15. Memory Model and I/O Tradeoffs

### 15.1 Peak RAM Components

```text
Per-contig catalog (walks present):   O(n_contig_sites × avg_walk_length)
                                      Released during GAF scan.
Per-contig catalog (walks released):  O(n_contig_sites × metadata_size)
                                      Released at end of contig loop iteration.
Compact index (during scan):          O(n_contig_sites × avg_allele_count × avg_walk_nodes × 4B)
                                      Released after scan.
site_to_chunk routing map:            O(n_contig_sites × site_id_string)
                                      Released at end of contig loop iteration.
chunk_catalogs (walks released):      O(n_contig_sites × metadata_size)
                                      Released at end of contig loop iteration.
Shard in-flight rows (per worker):    O(1 shard file per chunk, disk)
BamChunk batch:                       O(n_threads × one_chunk_reads × profile_width)
                                      Released after each batch.
rows_by_read accumulator:             O(n_observed_reads × ~60 bytes + allele_map)
                                      Grows across all contigs, released at end of run.
diagnostic_catalog:                   O(n_total_sites × metadata_size)
                                      Grows across all contigs; skip --graph-sites-tsv if not needed.
```

The dominant cost on the graph path is the per-contig catalog with walk strings, which is held
from catalog load until the GAF scan's compact-index build releases it. For typical human
chromosomes (~50 K–500 K sites), this is manageable; for a whole-genome VCF without `--contig`
the contig loop ensures only one chromosome's walks are live at a time.

### 15.2 GAF I/O: One Scan Per Contig

The current implementation performs one complete GAF scan per contig. For a whole-genome run
with 25 autosomes + sex chromosomes, this means 25–26 full passes over the GAF file.

**Why one scan per contig rather than one scan total**: the GAF scanner matches reads by
graph-node handles (§8.3). Building the compact index requires the allele walk data for all
sites in the catalog. Loading all chromosomes' walk data simultaneously would require holding
the entire genome-wide catalog in RAM — the RAM cost this architecture is designed to avoid.

**Future optimisation (node-routing single scan)**: because every graph node belongs to exactly
one chromosome, the first node in a read's walk handle determines the contig unambiguously.
A lightweight first pass over the VCF can build a `boundary_node_handle → ref_contig` routing
table (one entry per snarl, two handles, no walk strings). A single GAF pass using this table
would route raw GAF lines to per-contig intermediate files; the per-contig allele matching would
then run against those pre-filtered lines. This reduces GAF I/O from O(n_contigs × GAF_size) to
O(GAF_size + sum_contig_reads) and enables parallel contig phasing, at the cost of intermediate
disk space proportional to GAF size.

---

## 16. Summary of the Complete Pipeline

```text
phase_graph()
│
├── 3. Discover contigs (tabix seqnames / stream VCF)
│
├── Open all output files; write TSV headers
│
└── for each contig:
    │
    ├── §4  Load contig catalog from VCF (tabix or stream, RegionFilter)
    ├── §5  Split [0, catalog_end_bound) into chunk intervals
    ├── §6  Assign each site to its chunk (O(n_sites), arithmetic)
    │       Release walk strings from chunk_catalogs
    ├── §7  Create per-(worker × chunk) shard dir (RAII ChunkRowShard)
    │
    ├── §8  Parallel GAF scan
    │       Build CompactGraphSiteIndex from full_catalog (releases walk strings)
    │       n_workers stream GAF → match node handles → emit GraphReadAllele
    │       Callback routes each row to per-chunk shard file via site_to_chunk
    │       Writers flushed and closed
    │
    ├── §9–10 Streaming build → assign loop (batch of n_threads chunks at a time)
    │   │
    │   └── for each batch:
    │       │
    │       ├── run_vcf_gaf_chunks_parallel
    │       │   Each worker: load_chunk(i) → build_graph_bam_chunk → BamChunk
    │       │
    │       └── for each chunk result cur:
    │           ├── assign_graph_chunk_hap(cur)        [k-means]
    │           ├── §11 write site_counts / profiles / filtered_sites rows
    │           ├── if prev exists:
    │           │   ├── stitch_graph_chunk_pair(prev, cur)
    │           │   ├── §12.3 write variants_tsv rows for prev
    │           │   └── §13  merge_graph_chunk_into_read_rows(rows_by_read, prev)
    │           └── prev = cur
    │
    ├── Finalize last chunk (write variants + merge into rows_by_read)
    ├── §14.2 Write graph_read_support rows from shards
    └── Shard RAII destructor: delete temp files and dir
        Release full_catalog, chunk_catalogs, site_to_chunk, site_map
│
├── §14.1 Write graph_sites_tsv from diagnostic_catalog
│
└── §14.3 Second GAF pass (collect_gaf_read_names)
    Write graph_phase_reads_tsv from rows_by_read
    Write graph_phase_bam from rows_by_read
```
