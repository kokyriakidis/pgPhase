# `collect-bam-variation` Implementation Description

This document describes the implementation of the `pgphase collect-bam-variation` command from command-line input to final candidate output. The goal of this command is to collect and classify germline, non-mosaic **pre-phasing** candidate variation from one or more indexed BAM/CRAM files. The output is **not** a final phased callset: it is a structured table (TSV and optionally VCF) of candidate SNP, insertion, and deletion sites with read-level support counts, strand tallies, and quality categories aligned with longcallD’s candidate stage. That table is intended for downstream phasing or for parity checks against longcallD’s internal classifications. After classification, each chunk runs longcallD-shaped **read profiling** and **k-means haplotype clustering** (§18); the main TSV carries **`PHASE_SET`**, **`HAP_ALT`**, and **`HAP_REF`** from that scaffold. Adjacent chunks are stitched (§19) to merge phase blocks across boundaries; full noisy-region MSA (longcallD step 4) is not implemented.

The implementation is organized as a **staged pipeline** so that memory stays bounded on whole-genome runs and so each stage has a clear responsibility:

- **`collect_pipeline.cpp`** — End-to-end orchestration: parse regions and BED, build `RegionChunk` tiles with neighbour metadata (`reg_chunk_i`, prev/next overlap hints), stream batches to disk, and expose the `collect_bam_variation` CLI. This is the “control plane” of the tool.
- **`collect_var.cpp`** — The “biology core” ported from longcallD: gather candidate keys from digars, sweep reads for allele depths, then merge and filter chunk-level noisy regions (`cgranges`, `pre_process_noisy_regs` timing), run Fisher strand bias in ONT mode, run the two-pass `classify_cand_vars`-shaped classification, build read-level allele profiles (`collect_read_var_profile`), and co-phase clean germline het markers iteratively (`assign_hap_based_on_germline_het_vars_kmeans`).
- **`bam_digar.cpp`** — The “alignment plane”: stream BAM/CRAM into `ReadRecord` rows, apply read filters, build digars from EQX / MD / `cs` / reference comparison, and detect **per-read** noisy intervals using the same sliding-window rule as longcallD (see §8).
- **`collect_output.cpp`** — Serializers only: TSV header/body, optional candidate VCF (left-normalized, FILTER/INFO semantics for candidates), and optional read×site read-support TSV. No classification logic lives here.
- **`collect_phase.cpp` / `collect_phase.hpp`** — K-means read–haplotype clustering (`assign_hap_based_on_germline_het_vars_kmeans`), ported from longcallD `assign_hap.c`; category bitmask flags (`kCandCleanHetSnp`, …) for selecting which `VariantCategory` values participate in phasing.
- **`collect_types.hpp`** — Shared structs (`Options`, `BamChunk`, `ReadRecord`, `CandidateVariant`, `ReadVariantProfile`, per-read `haps` / `phase_sets`, …) and small RAII helpers (`SamFile`, `ReferenceCache`, HTS deleters).

Together, these files implement “longcallD-shaped” candidate collection in C++17 with explicit streaming and multi-BAM pooling.

### Source documentation (Doxygen)

The `collect_*` sources and headers are written for **Doxygen-style** extraction:

- Each major `.cpp` begins with an `@file` / `@brief` (and sometimes `@details`) block describing scope and conventions.
- Public and internal helpers use `@brief`, `@param`, `@return`, `@throws`, and `@note` where that clarifies contracts (for example, which functions reopen the BAM for SQ names, or when `std::runtime_error` is thrown on I/O failure).
- Header files (`collect_pipeline.hpp`, `collect_var.hpp`, `collect_output.hpp`, `collect_types.hpp`) duplicate or summarize API intent so a reader browsing headers alone sees the same contracts as the implementation.

Generating HTML from the tree is optional; the comments are also meant to be read in-place in the IDE.

### Vendored `cgranges` (longcallD fork, not vanilla lh3)

The project ships **`src/cgranges.c`** and **`src/cgranges.h`** as **vendored** sources (they live in-repo, not as a git submodule). For the current longcallD `main` branch on GitHub, these files are **byte-identical** to longcallD’s `src/cgranges.*`.

They are **not** the same as the minimal upstream library [lh3/cgranges](https://github.com/lh3/cgranges) alone: longcallD’s copy adds APIs this pipeline relies on, including:

- **`cr_merge` / `cr_merge2`** — merge nearby intervals with fixed and dynamic distance parameters (noisy-region geometry after extension).
- **`cr_is_contained`** — test whether a query span is fully covered by some stored interval (non-ONT containment filter for candidates inside noisy blocks).

**Why vendored instead of a submodule?** A plain clone plus `make` works everywhere; the exact C snapshot stays pinned to what was validated against longcallD; and release tarballs do not depend on recursive submodule checkout. If longcallD updates `cgranges` upstream, this tree should be refreshed deliberately (diff and test) rather than floating on a submodule pointer.

The main source files are:

```text
collect_pipeline.cpp   region parsing, BED/autosome filters, chunk neighbours, streaming run,
                       getopt CLI (collect_bam_variation), WorkerContext batch workers
collect_var.cpp        candidate sites, allele counts, noisy-region prep/post-process,
                       Fisher/strand (ONT), classify_cand_vars, collect_read_var_profile,
                       assign_hap_based_on_germline_het_vars_kmeans
bam_digar.cpp          BAM/CRAM iteration, filters, EQX/MD/cs/ref digars, per-read noisy windows
collect_output.cpp     TSV, optional VCF, optional read-support TSV
collect_phase.cpp      assign_hap_based_on_germline_het_vars_kmeans (longcallD assign_hap.c),
                       stitch_chunk_haps / flip_chunk_hap (longcallD stitch_var_main / flip_variant_hap)
collect_phase.hpp      kCand* flags (incl. kCandGermlineClean, kCandGermlineVarCate), assign_hap decl.,
                       stitch_chunk_haps decl.
collect_types.hpp      Options, BamChunk, ReadRecord, candidates, ReadVariantProfile,
                       hap_to_alle_profile, alle_covs / n_uniq_alles, ReferenceCache, …
main.cpp               dispatches `pgphase collect-bam-variation` to collect_bam_variation()
```

## 1. Command-Line Configuration

The command is invoked as:

```text
pgphase collect-bam-variation [options] <ref.fa> <input.bam|bam.list> [region ...]
```

The command-line parser builds an `Options` object. Important defaults are:

```text
threads = 1
min_mapq = 30
min_bq = 10
min_depth = 5
min_alt_depth = 2
min_af = 0.20
max_af = 0.80
chunk_size = 500000
max_var_ratio_per_read = 0.05
max_noisy_frac_per_read = 0.5
noisy_reg_merge_dis = 500
min_sv_len = 30
noisy_reg_max_xgaps = 5
min_noisy_reg_total_depth = 0 (0 = no extra minimum on overlapping reads; use 5 for the published “total coverage” rule)
verbose = 0 (use 2 for longcallD style AllNoisyRegions output)
HiFi noisy slide window = 100 bp
ONT noisy slide window = 25 bp
short-read noisy slide window = 25 bp
ONT strand-bias p-value threshold = 0.01
```

The input can be a single BAM/CRAM, a list of BAM/CRAM files, or a primary BAM/CRAM plus extra BAM/CRAM files:

```text
pgphase collect-bam-variation ref.fa sample.bam
pgphase collect-bam-variation --input-is-list ref.fa bam_list.txt
pgphase collect-bam-variation ref.fa part1.bam -X part2.bam -X part3.bam
```

Multiple input BAM/CRAM files are treated as files from the same sample. For each genomic chunk, reads from all input files are loaded and pooled before candidate collection. This means the candidate table represents the combined evidence of the sample, not independent per-file calls.

Two option semantics are worth making explicit:

```text
--include-filtered
  controls whether BAM reads marked QC-fail or duplicate are loaded.
  It does not mean "include filtered candidate categories in the output";
  the TSV writer already emits all final candidate categories.

--hifi
  selects HiFi-oriented behavior: a 100 bp default noisy-read window. Like
  longcallD `classify_var_cate`, the Fisher / strand-bias test is not applied
  in non-ONT modes. This is the default if no read technology flag is supplied.

--ont
  selects ONT-oriented behavior: a 25 bp default noisy-read window and
  Fisher exact testing for alternate-strand imbalance (longcallD
  `var_is_strand_bias` / `classify_var_cate`).

--short-reads
  selects short-read behavior: a 25 bp default noisy-read window. Strand-bias
  classification is ONT-only, same as longcallD: HiFi/short modes do not run
  the Fisher strand test.

Only one of --hifi, --ont, and --short-reads may be supplied.
```

The command supports explicit region restriction:

```text
pgphase collect-bam-variation ref.fa sample.bam chr11:1000-2000
pgphase collect-bam-variation ref.fa sample.bam -r chr11:1000-2000 -r chr12:1-500
pgphase collect-bam-variation ref.fa sample.bam --region-file regions.bed
pgphase collect-bam-variation ref.fa sample.bam --autosome
```

Region arguments are parsed as 1-based inclusive intervals. BED intervals are converted from 0-based half-open BED coordinates to 1-based inclusive coordinates by adding one to the BED start and keeping the BED end as the inclusive end.

Internally, most genomic positions stored in `RegionChunk`, `VariantKey`, `DigarOp`, `ReadRecord`, and `Interval` are 1-based inclusive coordinates. HTSlib iterator calls and the `cgranges` interval library use 0-based half-open intervals, so the code converts at the boundary:

```text
1-based inclusive interval: chr11:1001-1100
HTSlib/cgranges interval:   start = 1000, end = 1100
```

Insertions use the longcallD convention: an insertion at position `P` is between `P - 1` and `P`. This is why insertion interval operations often use `[P - 1, P)` in 0-based half-open form.

## 2. Reference and Alignment Resource Setup

Reference FASTA loading uses HTSlib's `fai_load3` with index creation enabled. This supports ordinary local FASTA files and remote FASTA paths supported by HTSlib. The `ReferenceCache` object fetches and caches one contig sequence at a time, normalizing bases to `A`, `C`, `G`, `T`, or `N`.

Each BAM/CRAM is opened through `SamFile`, which validates that the file format is BAM or CRAM. For CRAM input, the reference FASTA path is provided to HTSlib through `hts_set_fai_filename`, allowing CRAM decoding to retrieve reference bases.

For each worker thread, `WorkerContext` opens all input alignment files, reads their headers, loads their indexes, and owns a shared reference cache for that worker. This avoids sharing HTSlib file handles across worker threads.

Example:

```text
Input files:
  part1.bam
  part2.bam

Worker 1 opens:
  part1.bam + index
  part2.bam + index
  reference FASTA index

Worker 2 opens:
  part1.bam + index
  part2.bam + index
  reference FASTA index
```

This design makes chunk processing parallel while keeping each worker's HTSlib state independent.

## 3. Region and Chunk Construction

The pipeline first determines which genomic regions to process.

If the user supplied regions, a BED file, or `--autosome`, those filters are used. If no region filter is supplied, every contig present in both the BAM header and the FASTA index is processed.

The `--autosome` option checks chromosomes `1` through `22`. For each chromosome number, it first looks for a `chr`-prefixed contig such as `chr11`; if that is absent, it looks for the unprefixed form such as `11`. A chromosome is used only if it is present in both the BAM header and the FASTA index.

Each region is split into fixed-size chunks. The default chunk size is `500000` bp, configurable with `--chunk-size`.

Example:

```text
Requested region:
  chr11:1-1,200,000

With default chunk_size = 500000:
  chunk 1: chr11:1-500000
  chunk 2: chr11:500001-1000000
  chunk 3: chr11:1000001-1200000
```

Chunking provides two advantages. First, it limits the amount of read and candidate state held in memory at once. Second, it creates independent work units for parallel processing.

Reads can span chunk boundaries. Because each chunk queries reads overlapping that chunk, a long read may be loaded for two adjacent chunks. Candidate collection is still restricted to Digar events inside the current chunk boundaries. This means a boundary-spanning read can contribute evidence to candidates on both sides, while each chunk only creates candidate sites for its own interval. Fuzzy large-insertion deduplication runs **only inside each chunk** (longcallD `collect_all_cand_var_sites`); batch output **concatenates** per-chunk tables without a second contig-wide collapse, so near-boundary fuzzy duplicates across chunks are not merged at the candidate stage (same as longcallD).

## 4. Parallel Chunk Processing and Streaming Output

The central orchestrator is `run_collect_bam_variation` (`collect_pipeline.cpp`). It **streams** results to disk so the full merged candidate table for an entire genome need not sit in memory at once. The streaming loop groups chunks by `reg_chunk_i`, runs a thread pool over each group, merges only within the batch, and appends to open output streams—so peak RAM scales with **one contig’s batch** plus worker scratch space, not with every variant on the genome. For readers navigating the code, `collect_pipeline.cpp` / `.hpp` document each step (`load_region_chunks`, `collect_chunk_batch_parallel`, CLI parsing) with Doxygen-style comments.

High-level flow:

```text
load reference index (FAI)
build sorted RegionChunk list (`annotate_chunk_neighbors`: chunk_id, reg_chunk_i, prev/next overlaps)
open TSV (and optional VCF / read-support) streams; write headers
for each batch of consecutive chunks sharing the same reg_chunk_i:
    run collect_chunk_batch_parallel (atomic chunk claiming within the batch)
    merge_chunk_candidates(batch) → concat chunk tables only (fuzzy collapse already per chunk)
    append write_variants_tsv_records (and optional VCF / read-support rows)
close streams; print summary counts
```

**`reg_chunk_i` batches:** After `annotate_chunk_neighbors`, every chunk on the same reference sequence (contig) shares one `reg_chunk_i`; the batch loop processes **all chunks of one contig** together, concatenates their candidate tables, and appends. If the user supplies multiple disjoint region filters on different chromosomes, chunk order is sorted by contig and start, so each contig’s chunks still form their own batch. Batching is for streaming I/O and ordering; it does **not** re-run fuzzy collapse across chunks (longcallD parity).

**Workers:** Inside a batch, each worker owns a `WorkerContext`, claims the next chunk index with `fetch_add`, runs `process_chunk` (load reads → noisy prep → candidates → classify), and stores a `BamChunk` in a **fixed offset** so completion order does not scramble indices.

**`collect_chunks_parallel`:** A small helper still exists for “all chunks in one shot” (e.g. tests); the default CLI path uses the streaming loop above.

Example (same contig, five chunks, two threads inside one batch):

```text
Chunks on chr11 (same reg_chunk_i): C0, C1, C2, C3, C4

Workers steal indices until the batch is done:
  worker may process C0, C2, C4 / C1, C3 (order of claiming varies)

Merge for this batch:
  candidates(C0) ∥ … ∥ candidates(C4) → append to TSV (each Ci already fuzzy-collapsed internally)
```

## 5. Read Loading for a Chunk

For each chunk, the pipeline loads reads from every input BAM/CRAM using the BAM/CRAM index. The iterator query uses the chunk's contig id and half-open HTSlib coordinates:

```text
chunk: chr11:100001-200000
HTSlib query: tid(chr11), start=100000, end=200000
```

Each alignment record is filtered before conversion:

```text
discard if unmapped
discard if secondary
discard if supplementary
discard if MAPQ < min_mapq
discard QC-fail or duplicate unless --include-filtered is set
```

The exclusion of supplementary alignments is important. Many DNA structural-variant aligners represent split-read evidence using supplementary records and `SA` tags. This candidate collector intentionally does not assemble split alignments at this stage; it collects SNPs and indels from each primary alignment record's local alignment representation. Reference skips (`N` CIGAR operations) are also not treated as deletions.

For each accepted alignment, the code stores a `ReadRecord` with (among others):

```text
reference contig id
input BAM index (multi-input pooling)
1-based alignment start and end
strand
NM tag value if present
mapping quality
query name
base qualities (copy)
duplicate of the BAM alignment (bam1_t) for downstream use
Digar operations
read-level noisy regions
ONT palindrome flag, overlap-skip bookkeeping, phasing placeholders (haps, phase_sets, …)
read skip status
```

After loading reads from all input files, reads are sorted by:

```text
alignment start (ascending)
alignment end (descending — longer span first when start matches)
NM tag (ascending)
query name (ascending)
```

This deterministic sorting makes output stable across thread counts and input iteration details.

When more than one BAM/CRAM is provided, this sorting happens after reads from all inputs have been pooled for the chunk. The candidate collector therefore sees a single read collection per chunk:

```text
chunk reads = reads_from_part1 + reads_from_part2 + ... + reads_from_partN
sort pooled reads deterministically
process pooled reads once
```

## 6. Digar Construction

A `DigarOp` is the internal representation of an alignment event. It records:

```text
reference position
operation type
operation length
query/read position
low-quality flag
alternate sequence for SNPs and insertions
```

The supported Digar types are:

```text
Equal      aligned match block
Snp        single-base mismatch
Insertion  inserted read sequence relative to reference
Deletion   deleted reference sequence relative to read
SoftClip   soft-clipped read sequence
HardClip   hard-clipped read sequence
RefSkip    CIGAR N reference skip, e.g. spliced RNA alignment
```

The code can build Digars from several alignment encodings:

```text
1. EQX CIGAR if the CIGAR uses = and X without M
2. cs tag if present
3. MD tag if present
4. reference comparison fallback for M CIGAR operations
```

This priority matters because different aligners encode mismatches differently. If a `cs` or `MD` path fails to parse safely, the code falls back to reference comparison.

The parser priority also preserves information when the aligner has already emitted explicit mismatch operations. For example, an EQX CIGAR has separate `=` and `X` operations, so a mismatch can be emitted without comparing every aligned base to the reference. If only ordinary `M` operations are available, the implementation needs either a `cs` tag, an `MD` tag, or the reference sequence to determine which aligned bases are matches versus SNPs.

Adjacent `Equal` Digars are merged to keep match runs compact. Adjacent insertions and deletions are deliberately not coalesced in the generic append path, because longcallD keeps those events separate in its CIGAR/MD-derived representation. Preserving this behavior keeps candidate sorting and large-insertion collapsing compatible with longcallD.

Examples of tag-derived parsing:

```text
MD tag example:
  CIGAR: 10M
  MD:    4A5

Interpretation:
  4 matching reference bases
  1 mismatch where the reference base is A
  5 matching reference bases

The read base at the mismatch position becomes the SNP alternate allele.
```

```text
cs tag example:
  cs: :4*ag:5

Interpretation:
  :4    four matches
  *ag   reference a, query g, therefore SNP alt G
  :5    five matches
```

For insertions and deletions, the `cs` tag uses:

```text
+ACG  insertion of ACG in the read
-TTA  deletion of TTA from the reference
```

The code still checks the CIGAR and read coordinates while parsing these tags, because the tag must be projected onto the read and reference positions.

### 6.1 Example: SNP from Reference Comparison

Suppose the reference and read are:

```text
Reference: A C G T A
Read:      A C T T A
                 ^
```

At the third base, the reference has `G` and the read has `T`. The Digar emitted is:

```text
type = Snp
pos = reference position of G
len = 1
alt = "T"
low_quality = base_quality(T) < min_bq
```

If the base quality is at least `min_bq`, this SNP also contributes one unit to the read-level noisy-event sliding window.

For SNPs, low quality is determined from the base quality of the read base carrying the alternate allele.

### 6.2 Example: Insertion

CIGAR:

```text
10M3I20M
```

If the insertion occurs after reference position `1009` and the inserted read bases are `TGA`, the internal representation is:

```text
type = Insertion
pos = 1010
ref_len = 0
alt = "TGA"
```

Insertions use the position convention inherited from longcallD: an insertion at `pos` is between `pos - 1` and `pos`. The inserted sequence is stored as the alternate allele.

An insertion is marked low-quality only if all inserted bases have quality below `min_bq`. If at least one inserted base meets the base-quality threshold, the insertion Digar is not marked low-quality at construction time.

### 6.3 Example: Deletion

CIGAR:

```text
10M4D20M
```

If the deletion starts at reference position `1010`, the Digar is:

```text
type = Deletion
pos = 1010
ref_len = 4
alt = ""
```

The deleted reference sequence is not stored in the `VariantKey`; it is fetched later from the reference when TSV or VCF output is written.

A deletion has no read bases inside the deleted reference interval. Therefore its low-quality check uses the qualities of the read bases flanking the deletion. The deletion is treated as high-quality only when the left and right anchors are available and pass the base-quality threshold, with special handling at the start or end of the read.

### 6.4 Example: Reference Skip

CIGAR:

```text
50M1000N50M
```

The `N` operation means the alignment skips `1000` reference bases without consuming read bases. This is common in spliced RNA alignments and is not treated as a deletion candidate. It is stored as `RefSkip`, which prevents it from being counted as a genomic deletion.

The difference between deletion and reference skip is:

```text
D: consumes reference bases and represents missing sequence in the read
N: consumes reference bases but represents an alignment jump, usually a splice gap
```

Only deletions become deletion candidates. Reference skips are ignored by candidate collection.

## 7. Per-Read Quality and Skip Decisions

After Digars are built, the read is checked for excessive variation or excessive noisy-region coverage.

The first skip rule is based on the density of candidate events:

```text
total_cand_events > mapped_length * max_var_ratio_per_read
```

With the default `max_var_ratio_per_read = 0.05`, a 10,000 bp read is skipped if it has more than 500 candidate events.

Example:

```text
read length = 10000
candidate events = 700
threshold = 10000 * 0.05 = 500
700 > 500, so the read is skipped
```

The second skip rule is based on the fraction of the read covered by read-level noisy regions:

```text
noisy_region_length_on_read > mapped_length * max_noisy_frac_per_read
```

With the default `max_noisy_frac_per_read = 0.5`, a 10,000 bp read is skipped if more than 5,000 bp are covered by noisy regions.

Skipped reads remain in the chunk's read list, but downstream collection and noisy-region support counting ignore them. This preserves debug visibility while preventing poor alignments from dominating candidate discovery.

This distinction creates two levels of filtering:

```text
alignment-level filtering:
  read is not loaded at all if it is unmapped, secondary, supplementary, low MAPQ, etc.

read-level skipping after Digar construction:
  read was loaded and parsed, but is ignored for candidate/noisy support because it is too variant-dense or too noisy
```

The second level requires Digar construction first, because the program cannot know whether a read has too many events or too much noisy coverage until the alignment has been converted into Digar operations.

## 8. Read-Level Noisy-Region Detection

Noisy-region detection begins **while** Digars are being built for each read. The implementation walks the alignment once; whenever a non–low-quality SNP, insertion, or deletion is recorded, it feeds a **sliding-window accumulator** that sums event “weight” inside a fixed genomic width on the reference. If that sum exceeds a threshold, the read is considered locally unreliable and a **read-level noisy interval** is opened, extended, or merged—**before** chunk-level union in §9; that list is refined after allele counting in §12.

This design mirrors longcallD’s intent: dense mismatch and indel signal in a short window is a proxy for local misalignment or excessive micro-errors, so simple pileup-style calling there is untrustworthy.

### Parity with longcallD (`xid_queue_t` / `push_xid_size_queue_win`)

In longcallD, the same algorithm lives in `bam_utils.c` as `xid_queue_t` and `push_xid_size_queue_win`. In pgPhase, the equivalent logic is **`XidQueue`** and **`xid_push_win`** in `bam_digar.cpp`, wrapped by **`NoisyRegionBuilder`** (constructor picks window width from technology, `observe_variant` pushes events, `flush` emits a trailing open interval).

**The sliding-window math is the same:** push `(pos, len, count)`, subtract evicted events from the front while `pos[front] + len[front] - 1 <= pos - win`, and when `count > 0` and the running sum **`> max_s`** (`noisy_reg_max_xgaps`), compute `noisy_start = pos[front]`, `noisy_end = pos[rear] + lens[rear]`, then merge or flush intervals using the same adjacency rule and the same **`var_size = max(sum of counts, genomic span)`** label.

**Intentional structural differences** (output-equivalent for interval geometry):

| Aspect | longcallD | pgPhase |
|--------|-----------|---------|
| **Where intervals go** | Immediately `cr_add` into a per-read `cgranges_t` (0-based half-open storage). | Append **`Interval{start,end,label}`** to `ReadRecord::noisy_regions` (**1-based inclusive**). Chunk code later converts via `intervals_to_cr` / merge when it needs trees. |
| **`is_dense[]` flag** | Allocated and set to 1 for queue indices when noisy; **never read** elsewhere in that file. | **Omitted** (dead state in longcallD). |
| **Eviction loop** | `while` without `front <= rear` guard. | **`while (front <= rear && …)`** so the loop never indexes past an empty queue (defensive port). |

So the **detection rule and merge behaviour match longcallD**; pgPhase differs in **C++ data layout** (vectors on `ReadRecord`, defer building `cgranges` until chunk merge) and in **dropping unused bookkeeping**.

The effective window is:

```text
--noisy-slide-win if provided
25 bp in ONT mode
25 bp in short-read mode
100 bp in HiFi mode
```

The maximum allowed event load in the window is `noisy_reg_max_xgaps`, default `5`. Event contributions are:

```text
SNP: 1
Insertion: inserted length
Deletion: deleted length
```

Example:

```text
HiFi window = 100 bp
noisy_reg_max_xgaps = 5

Events on a read:
  chr11:1000 SNP, contribution 1
  chr11:1010 INS length 2, contribution 2
  chr11:1030 DEL length 3, contribution 3

Total contribution in 100 bp = 1 + 2 + 3 = 6
6 > 5, so a read-level noisy interval is started.
```

The noisy interval spans from the earliest event currently in the window to the latest event. Adjacent or overlapping noisy intervals on the same read are merged.

The interval also carries a numeric label. When a dense event window is flushed, the label is the larger of:

```text
sum of event contributions in the window
span length of the interval
```

This label is passed through the interval data structure and is used by the merge routines as interval metadata. Conceptually, it records how large or severe the noisy interval is.

Long terminal clips also create noisy regions. A soft or hard clip of at least `30 bp` at a read end creates a noisy interval extending `100 bp` into the reference from the clip edge.

Example:

```text
CIGAR: 40S160M
alignment starts at chr11:5000

40S is a long left-end clip.
Add noisy interval approximately chr11:5000-5100.
```

This captures places where the read alignment begins or ends abruptly, which may indicate local misalignment, structural variation, or unresolved sequence.

Only terminal clips are used for this long-clip noisy-region rule. Internal clips are not treated the same way here. A left-end clip creates a region starting at the current reference position and extending to the right; a right-end clip creates a region extending leftward into the reference.

## 9. Chunk Finalization

After reads have been loaded and sorted for a chunk, the chunk is finalized:

```text
compute base-quality histogram
load reference slice with flanking sequence
detect low-complexity reference intervals
collect and merge read-level noisy intervals into chunk-level noisy intervals
```

The reference slice spans active reads plus a `50000 bp` flank on each side, clipped to the contig boundary. If all reads are skipped, the slice falls back to the chunk interval itself.

Low-complexity sequence is detected with `sdust` over the reference slice. These intervals are later used to extend noisy regions and to classify repeat-associated indels.

The current `sdust` parameters are:

```text
threshold = 5
window = 20
```

The low-complexity intervals are stored in genomic coordinates relative to the chunk reference slice. They are not candidate variants by themselves; they are context annotations used by later noisy-region and repeat-indel logic.

Example:

```text
chunk: chr11:100000-200000
active reads span chr11:95000-210000
reference flank = 50000

reference slice:
  chr11:45000-260000
```

The chunk-level noisy-region list is initially formed by concatenating noisy intervals from non-skipped reads and merging overlapping or adjacent intervals.

## 10. Candidate Site Collection

After chunk finalization (§9), the collector builds a raw candidate site list from non-skipped reads. This matches longcallD `collect_all_cand_var_sites` (including fuzzy large-insertion collapse within the chunk).

For each read and each Digar:

```text
ignore skipped reads
ignore non-variant Digars such as Equal, clips, and RefSkip
ignore low-quality variant Digars
ignore variant Digars outside the chunk boundaries
convert SNP/INS/DEL Digar to VariantKey
append CandidateVariant with empty counts
```

Low-quality variant Digars are not allowed to create new candidate sites. However, if another high-quality read already created the candidate site, a low-quality matching observation can still be counted later as `low_qual_cov` during allele counting. This avoids discovering sites from weak evidence while still preserving information about low-quality observations at discovered sites.

Example raw events:

```text
read1: chr11:1000 SNP A>G
read2: chr11:1000 SNP A>G
read3: chr11:1005 DEL length 2
read4: chr11:1010 INS T
```

The raw candidate table initially contains all four observations. It is then sorted and collapsed into unique candidate sites:

```text
chr11:1000 SNP A>G
chr11:1005 DEL length 2
chr11:1010 INS T
```

The duplicated SNP appears once as a site. Allele counts are computed in §11.

Large insertions use fuzzy collapsing. Insertions of at least `30 bp` at the same position are considered equivalent if their lengths are similar:

```text
min(length1, length2) >= 0.8 * max(length1, length2)
```

Example:

```text
read1: chr11:2000 INS length 52
read2: chr11:2000 INS length 55

52 >= 0.8 * 55 = 44
These are collapsed into one candidate site.
```

Small insertions require exact sequence equality:

```text
chr11:2000 INS A
chr11:2000 INS T

These remain separate candidate sites.
```

**Deletions (no fuzzy length merge):** longcallD `exact_comp_var_site_ins` applies the 0.8 length rule **only to large insertions**. Deletions are the same site only when they share the same breakpoint and **`ref_len`** (deletion length on the reference). Two deletions at the same position with **different** lengths (e.g. 40 bp vs 50 bp) remain **separate** candidate rows—there is no longcallD-style fuzzy merge for large DELs, and pgPhase matches that.

The sort order is deterministic and longcallD-compatible. It compares:

```text
contig id
sort position
variant type
reference length
alternate length
alternate sequence for SNPs and small insertions
```

For SNPs, the sort position is the SNP position. For insertions and deletions, the sort position is `pos - 1`, matching the insertion-between-bases convention.

## 11. Allele Count Collection

Once unique candidate sites are known (§10), the pipeline performs a second pass over non-skipped reads to count support for each site (`collect_cand_vars` in longcallD; `collect_allele_counts_from_records` here). **Next**, chunk-level noisy intervals are refined (§12), and only then does variant classification (§13) run.

For each read, the algorithm starts near the first candidate overlapping the read. It then walks the sorted candidate table and the read's variant Digars together.

For each candidate:

```text
if the read has a matching Digar:
    count an alternate observation
else if the read spans the candidate position:
    count a reference observation
else:
    no observation
```

The matching test uses the same longcallD-compatible comparison used for candidate collapsing. Therefore, a read carrying a large insertion can support a candidate large insertion even if the inserted sequence is not exactly identical, as long as the large-insertion fuzzy comparison considers them equivalent. This is important for long reads, where large inserted alleles may be represented with slightly different lengths or sequences across alignments.

Reference observations are counted only when the read reaches the candidate position. A read that ends before the site, starts after the site, maps to another contig, or is skipped contributes no observation for that candidate.

Alternate observations are marked low-quality if:

```text
the Digar was already low-quality
or the average base quality supporting the Digar is less than min_bq
```

Low-quality alternate observations increment `low_qual_cov` and do not increment `total_cov` or `alt_cov`. This means the classifier can use `total_cov` for high-quality allele fraction while still using `total_cov + low_qual_cov` when checking whether the site had enough overall read evidence.

Example:

```text
Candidate: chr11:1000 SNP A>G

10 high-quality reads support ref
1 high-quality read supports alt
3 low-quality reads support alt

total_cov = 11
ref_cov = 10
alt_cov = 1
low_qual_cov = 3
depth_with_low_quality = 14
```

This candidate still fails `min_alt_depth = 2`, because only high-quality alternate observations contribute to `alt_cov`.

For deletions, the average quality used in this second pass is computed over the same anchor bases that represent the deleted interval in the read. For insertions and SNPs, it is computed over the read bases carrying the alternate sequence.

Counts tracked per candidate are:

```text
total_cov
ref_cov
alt_cov
low_qual_cov
forward_ref
reverse_ref
forward_alt
reverse_alt
allele_fraction
```

Example:

```text
Candidate: chr11:1000 SNP A>G

read1 spans 1000 and has G: alt, forward
read2 spans 1000 and has A: ref, reverse
read3 spans 1000 and has G: alt, reverse
read4 does not cover 1000: ignored

total_cov = 3
ref_cov = 1
alt_cov = 2
forward_alt = 1
reverse_alt = 1
allele_fraction = 2 / 3 = 0.667
```

If `--read-support` is enabled, the same pass writes per-read observations. This output is useful for downstream phasing because it records whether each read supports the reference or alternate allele at each candidate site.

## 12. Pre-Processing Noisy Regions

After candidate sites and allele depths are in place (§10–§11), chunk-level noisy regions are refined (`pre_process_noisy_regs`). The provisional list is the read-union from §9; this step extends through low-complexity sequence, merges nearby intervals, and drops intervals with insufficient read support. Classification (§13) uses the refined mask, matching longcallD `collect_var_main`.

The preprocessing order is:

```text
raw read-supported noisy intervals
extend through overlapping low-complexity intervals
merge nearby noisy intervals
count read support for each merged interval
discard noisy intervals with insufficient support
```

### 12.1 Extension Through Low-Complexity Regions

If a noisy interval overlaps a low-complexity interval, it is extended to include that low-complexity sequence.

Example:

```text
noisy interval:        chr11:1000-1030
low-complexity repeat: chr11:1020-1100

extended noisy region: chr11:1000-1100
```

This is useful because indel alignment around homopolymers and short tandem repeats is often unstable. A small noisy signal near the edge of a repeat may represent the entire repeat region rather than only the original small interval.

### 12.2 Merge Nearby Noisy Regions

After extension, nearby noisy intervals are merged using `noisy_reg_merge_dis`, default `500 bp`, and `min_sv_len`, default `30 bp`.

Example:

```text
noisy A: chr11:1000-1100
noisy B: chr11:1300-1400
distance = 199 bp

Because 199 <= 500, merge:
  chr11:1000-1400
```

This prevents a complex locus from being split into many small intervals.

### 12.3 Read-Support Filter for Noisy Regions

The code then counts support for each merged noisy interval.

For each non-skipped read:

```text
if the read overlaps any part of the noisy interval:
    increment n_total
    if the read has a read-level noisy interval overlapping it:
        increment n_noisy
```

The noisy region is kept only if the following hold:

```text
n_noisy >= min_alt_depth
n_noisy / n_total >= min_af
```

If `min_noisy_reg_total_depth` is greater than zero, the implementation also requires:

```text
n_total >= min_noisy_reg_total_depth
```

The default is `0`, which does **not** add an extra check (this matches the longcallD `pre_process_noisy_regs` logic, which only uses the two conditions above with `n_noisy` and the ratio of `n_noisy` to `n_total`). The published LongcallD methods description also lists a **minimum total overlapping read count**; set `--min-noisy-reg-total-depth 5` to enforce that third gate. When `0`, there is no separate “total read coverage” minimum.

With defaults (when `min_noisy_reg_total_depth` is 0 or unset in code):

```text
min_alt_depth = 2
min_af = 0.20
```

Example dropped region:

```text
merged noisy region: chr11:1000-1100
20 reads overlap any part of it
1 read has noisy evidence there

n_noisy = 1 < 2
The interval is removed from chunk.noisy_regions.
```

Example kept region:

```text
merged noisy region: chr11:1000-1100
20 reads overlap any part of it
6 reads have noisy evidence there

n_noisy = 6
n_noisy / n_total = 6 / 20 = 0.30
0.30 >= 0.20
The interval remains a trusted noisy region.
```

This step prevents one poor read from causing an entire locus to be treated as noisy.

Because this support check happens after low-complexity extension and noisy-region merging, the denominator and numerator are measured on the merged interval, not on each original read-level interval. In other words, the code asks whether the final candidate noisy locus is supported by enough reads.

## 13. Initial Variant Classification

Each candidate is classified using coverage, allele fraction, and local sequence context. A Fisher
strand-imbalance test is applied only in ONT mode, matching longcallD (§13.2).

The first computed value is:

```text
depth_with_low_quality = total_cov + low_qual_cov
allele_fraction = alt_cov / total_cov if total_cov > 0 else 0
```

The classification order is important.

### 13.1 Low Coverage

A candidate is `LOW_COV` if:

```text
depth_with_low_quality < min_depth
or alt_cov < min_alt_depth
```

Example:

```text
min_depth = 5
min_alt_depth = 2

candidate depth = 4
alt_cov = 3

depth 4 < 5, so category = LOW_COV
```

Another example:

```text
depth = 20
alt_cov = 1

alt_cov 1 < 2, so category = LOW_COV
```

### 13.2 Strand Bias (ONT only)

This matches longcallD `classify_var_cate`: the Fisher / strand-bias check runs
**only** when the read-technology mode is ONT (`opt->is_ont` in longcallD).

In ONT mode, the code applies a Fisher exact test to alternate forward and
reverse support. It compares observed alternate strand counts to an expected
balanced split. If the two-tailed p-value is below `strand_bias_pval`, the
candidate is `STRAND_BIAS`.

Example:

```text
forward_alt = 10
reverse_alt = 0
```

This is suspicious because all alternate observations are on one strand; with
sufficient support the two-tailed p-value can fall below the threshold, giving
`STRAND_BIAS`.

In **HiFi** and **short-read** modes, there is no separate strand-bias
category from this test: the implementation follows longcallD and does not run
Fisher (or a substitute heuristic) for non-ONT.

### 13.3 Low Allele Fraction

If allele fraction is below `min_af`, the initial category is `LOW_AF`.

Example:

```text
alt_cov = 3
total_cov = 30
AF = 0.10
min_af = 0.20

category = LOW_AF
```

Later, this is converted to `LOW_COV` in the longcallD-compatible classification pass.

### 13.4 Clean Homozygous Candidate

If allele fraction is above `max_af`, the category is `CLEAN_HOM`.

Example:

```text
alt_cov = 19
total_cov = 20
AF = 0.95
max_af = 0.80

category = CLEAN_HOM
```

This indicates that nearly all reads support the alternate allele.

### 13.5 Repeat-Associated Heterozygous Indel

For insertions and deletions, the classifier checks whether the indel is in a homopolymer or short tandem repeat context. This only applies to short indels controlled by `noisy_reg_max_xgaps`, default `5`.

The homopolymer/repeat decision is made from the reference sequence around the candidate, not from the read names or read sequences alone. For short indels, the classifier checks two related patterns.

First, it tests for homopolymer or short-repeat context by scanning both sides of the indel for repeated units up to length `6`. The code checks whether at least three copies of the same unit are visible next to the variant.

Example homopolymer with unit length 1:

```text
Reference:  C A A A A A A A G
Unit:         A
Copies:       A A A ...
Candidate:     delete one A
```

Because the nearby reference consists of repeated `A` bases, the deletion is considered homopolymer-associated.

Example short tandem repeat with unit length 2:

```text
Reference:  G C A C A C A C A T
Unit:         C A
Copies:       CA CA CA ...
Candidate:       insert or delete CA
```

Second, the repeat-region check compares the reference sequence before and after the indel. For a deletion of length `L`, it compares approximately `3 * L` bases starting at the deletion with `3 * L` bases after the deletion. For an insertion of length `L`, it asks whether inserting the alternate sequence would preserve a local repeated pattern. These tests identify short tandem repeat contexts where equivalent alignments can be shifted left or right.

Example homopolymer:

```text
Reference: AAAAAAAA
Candidate: delete one A
```

Such a candidate is classified as:

```text
REP_HET_INDEL
```

The biological reason is that short indels in homopolymers and tandem repeats are frequently represented inconsistently by read alignments.

### 13.6 Clean Heterozygous Candidate

If no previous filter applies, the candidate becomes:

```text
CLEAN_HET_SNP for SNPs
CLEAN_HET_INDEL for insertions and deletions
```

Example:

```text
candidate: chr11:1000 A>G
depth = 30
alt_cov = 15
AF = 0.50
forward_alt = 8
reverse_alt = 7
not in noisy region

category = CLEAN_HET_SNP
```

## 14. Noisy-Region Feedback During Classification

Classification also feeds information back into the noisy-region model.

First, in non-ONT mode, including HiFi and short-read mode, if a candidate overlaps a trusted chunk noisy region, it is classified as `NON_VAR`.

Example:

```text
trusted noisy region: chr11:1000-1100
candidate SNP: chr11:1050 A>G

category = NON_VAR
```

This means the candidate is not trusted as a clean pre-phasing site.

Second, repeat-associated indels are added back to a noisy-region interval set. If the repeat indel overlaps low-complexity sequence, the noisy interval is extended to include the low-complexity span.

Example:

```text
candidate: chr11:2000 deletion of A
reference context: AAAAAAAAAA

Add candidate span, extended through the homopolymer, to noisy regions.
```

Third, dense overlapping candidate sites are treated as noisy. The code builds an interval index of candidate positions. If a candidate overlaps more than one candidate position, the locus is considered dense.

Example:

```text
chr11:3000 SNP A>G
chr11:3000 INS T
chr11:3001 DEL length 1
```

These overlapping signals are unlikely to be a simple isolated variant. The region can be added to the noisy interval set, provided the local noisy-read ratio is high enough.

The newly generated noisy intervals are merged with the previous chunk noisy intervals using `cr_merge2`.

The local noisy-read ratio check for dense loci is the same conceptual ratio used elsewhere:

```text
reads with read-level noisy evidence in the candidate span
---------------------------------------------------------
all non-skipped reads overlapping the candidate span
```

The dense locus is added to the noisy interval set only if this ratio is at least `min_af`. Repeat-associated indels are added without this extra ratio check, because the repeat classification itself is already sequence-context evidence that the locus is alignment-ambiguous.

## 15. Post-Processing Noisy Regions

After classification, noisy intervals are expanded slightly and merged again.

Each noisy region starts with a default flank of `10 bp` on each side. The algorithm then **walks outward** from both ends through the sorted candidate list, extending the boundary as long as it finds eligible candidates adjacent to the current edge. The walk stops only when a gap of more than 1 bp separates the next candidate from the current boundary.

Categories skipped for this extension are:

```text
LOW_COV
STRAND_BIAS
NON_VAR
```

All other categories — `CLEAN_HET_SNP`, `CLEAN_HET_INDEL`, `CLEAN_HOM`, `REP_HET_INDEL`, `LOW_AF`, `NOISY_CAND_HET`, `NOISY_CAND_HOM`, `NOISY_RESOLVED` — are eligible to pull the boundary further.

**Chain expansion detail:** This is not a single-step "nearest neighbor" lookup. The boundary keeps moving as long as the next eligible candidate is within 1 bp of it, so a dense cluster of variants adjacent to a noisy region is absorbed entirely. Only a genuine sequence gap (> 1 bp between the current boundary and the next candidate) breaks the expansion. For each candidate that qualifies, the boundary is pulled to `variant_start - 10` (left walk) or `variant_end + 10` (right walk).

Example — single nearby candidate (does not touch):

```text
noisy region: chr11:1000-1100
initial extension: chr11:990-1110

right walk — next eligible candidate at chr11:1115-1117
  vs=1115 > cur_end+1=1111 → gap detected, walk stops
  cur_end stays 1110
```

Example — chain of adjacent candidates:

```text
noisy region: chr11:1000-1100
initial extension: chr11:990-1110

right walk:
  candidate at 1108-1110: vs=1108 ≤ 1111 → cur_end = 1110+10 = 1120
  candidate at 1118-1120: vs=1118 ≤ 1121 → cur_end = 1120+10 = 1130
  candidate at 1130-1135: vs=1130 ≤ 1131 → cur_end = 1135+10 = 1145
  next candidate at 1300:  vs=1300 > 1146 → gap, walk stops
  → final noisy region extends to 1145
```

After extension, overlapping intervals are merged.

Post-processing is deliberately performed after candidate classification, because the final noisy-region boundaries can depend on nearby candidate categories. A clean nearby variant may extend a noisy locus, while `LOW_COV`, `STRAND_BIAS`, and `NON_VAR` sites do not extend it.

## 16. Final Noisy-Containment Filter

In non-ONT mode, including HiFi and short-read mode, the pipeline performs a final containment sweep. Any candidate contained in a finalized noisy region is marked `NON_VAR`.

Example:

```text
final noisy region: chr11:5000-5100
candidate: chr11:5050 SNP C>T

category becomes NON_VAR
```

For insertions, containment is tested over the insertion's anchor interval:

```text
insertion at pos P is treated as interval [P-1, P)
```

This final filter is not applied in ONT mode in the current implementation, matching the ported longcallD behavior for this part of the pipeline.

## 17. Resolving Large Noisy Candidates

After classification, candidates with category `NOISY_CAND_HET`, `NOISY_CAND_HOM`, or `REP_HET_INDEL` can be converted to `NOISY_RESOLVED` if they are large enough:

```text
insertion length >= 30
or deletion length >= 30
```

Example:

```text
candidate: chr11:10000 deletion length 45
category before resolution: REP_HET_INDEL

45 >= 30
category after resolution: NOISY_RESOLVED
```

This distinguishes large events that originated from noisy or repeat-associated logic but are large enough to be retained as resolved structural-like candidates.

In the current **digar collect** flow, only `REP_HET_INDEL` is typically promoted here: longcallD `classify_cand_vars` does not assign `NOISY_CAND_HET` / `NOISY_CAND_HOM` (`e` / `h`) to BAM-sweep candidates (those labels are for MSA-recalled variants inside noisy regions). The initial classifier most commonly produces `LOW_COV`, `STRAND_BIAS`, `LOW_AF`, `CLEAN_HOM`, `REP_HET_INDEL`, `CLEAN_HET_SNP`, `CLEAN_HET_INDEL`, and `NON_VAR`.


## 18. Phasing Scaffold (Steps 3.1–3.2)

Per-chunk biology is driven by **`collect_var_main`** (`collect_var.cpp`), which mirrors longcallD’s numbered `collect_var_main`: steps **1.x** (sites and allele counts), **2.x** (noisy prep, `classify_chunk_candidates`, noisy post-process, non-ONT containment), then **3.1–3.2** (read profiles and k-means phasing). The worker entry point is **`process_chunk`** → **`collect_var_main`** (`collect_pipeline.cpp`).

Phasing runs **only when** `chunk.candidates` is non-empty after classification (same guard pattern as longcallD: no candidates ⇒ no profile or k-means work).

### Step 3.1: `collect_read_var_profile` (static, `collect_var.cpp`)

For each non-skipped read, this walks digars and the sorted candidate table in lockstep (same merge logic as the allele-depth sweep: `exact_comp_var_site_ins`, `kLongcalldMinSvLen`). It fills:

- **`ReadVariantProfile`** on `chunk.read_var_profile[read_i]`: sparse `alleles` / `alt_qi` for variant indices from `start_var_idx` through `end_var_idx`.
- Allele codes match longcallD profiling: **`0`** = ref, **`1`** = alt (if base quality passes `min_bq` and the digar is not marked low-quality), **`-2`** = low-quality alt observation (skipped in hap scoring), and trailing sites inside the read span with no alt digar match get **ref (`0`)** via the tail loop.

It also builds **`chunk.read_var_cr`**: a `cgranges` interval index with half-open **`[start_var_idx, end_var_idx + 1)`** per read, labeled by `read_i`. That index supports **`cr_overlap`** queries “reads covering variant *i*” used in the pivot sweep and in adjacent-het linkage counts—same geometry as longcallD `collect_read_var_profile`.

Reads are visited in **`ordered_read_ids`** order when that vector is populated (same as longcallD’s loop over `chunk->ordered_read_ids[i]`); otherwise indices **`0 … n_reads-1`** are used (equivalent after **`load_and_prepare_chunk`** sorts reads by start, end (desc), NM, qname).

### Step 3.2: `assign_hap_based_on_germline_het_vars_kmeans` (`collect_phase.cpp`)

Public entry declared in **`collect_phase.hpp`**. It implements longcallD **`assign_hap_based_on_germline_het_vars_kmeans`**.

**Single call in phase 3 (`collect_var_main` when `n_cand_vars > 0`).** longcallD runs **`assign_hap_based_on_germline_het_vars_kmeans`** once with
**`LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR`**
(= **`LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE`** / **`kCandGermlineClean`**). Comments in longcallD describe further rounds (`+ NOISY_CAND_*`, somatic); the k-means that includes **`LONGCALLD_CAND_GERMLINE_VAR_CATE`** is issued **later**, after noisy-region MSA (step 4) when new variants exist — not in the same `if (n_cand_vars > 0)` block. pgPhase now mirrors that germline/noisy step-4 path: MSA-recalled `NOISY_CAND_*` variants are merged into the chunk profile and trigger a second k-means pass with **`kCandGermlineVarCate`**. The separate longcallD `out_somatic` branch remains outside the current pgPhase model.

**Other masks.** **`kCandHetVarCate`** matches **`LONGCALLD_CAND_HET_VAR_CATE`**. **`kCandGermlineVarCate`** matches **`LONGCALLD_CAND_GERMLINE_VAR_CATE`** for use if a future MSA step mirrors longcallD’s follow-up k-means.

**Category mask (general).** Callers pass `kCand*` flags; see **`collect_phase.hpp`**.

**Per-variant state (on `CandidateVariant`).**

- **`hap_to_alle_profile[1]` / `[2]`**: per-hap read counts per allele index. Width is **`variant_allele_slots`**: default **2** (ref/alt from `VariantCounts`). If **`VariantCounts::alle_covs`** is populated, length follows **`n_uniq_alles`** / vector size so multi-allele sites can match longcallD’s `n_uniq_alles` (read observations remain 0/1 today unless extended).
- **`hap_to_cons_alle[0..2]`**: consensus allele indices; het sites start at `-1` for haps 1–2 until profiles accumulate.
- **`get_var_init_max_cov_allele`**: scans **`alle_covs`** with strict **`>`** (ties keep lower index), else **`ref_cov` / `alt_cov`**, matching longcallD’s `get_var_init_max_cov_allele` tie behavior.

**Algorithm phases (same structure as `assign_hap.c`).**

1. **Phase 1 — pivot sweep:** **`select_init_var`** picks the deepest `CLEAN_HET_SNP`, else `CLEAN_HET_INDEL`, else noisy het SNP/indel (homopolymer indels excluded from noisy-indel pivot). Variants are visited in an outward order from that pivot index within the **valid** (flag-filtered) list. For each variant, overlapping reads with **`haps[read_i] == 0`** get **`init_assign_read_hap`**, then **`update_var_hap_profile_cons_alle`**. Reads that still have no informative score are seeded as hap **1** (longcallD behavior). Hom categories are skipped in this round.
2. **Phase 2 — up to 10 iterations:** **`iter_update_var_hap_cons_phase_set`** counts spanning-read agreement vs conflict between **adjacent phased het** variants (via `read_var_cr` overlap on global variant indices). Weak linkage (`n_agree < 2` and `n_conflict < 2`) starts a new **phase_set** anchor at **`VariantKey::sort_pos()`** (SNP → `pos`, indel → `pos - 1`, matching longcallD). When conflicts dominate, **`flip`** toggles; the consensus swap uses longcallD’s loop **`for (hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap)`** swapping `[hap]` with `[3-hap]`—at diploid ploidy **2** this is **two swaps** and leaves **`hap_to_cons_alle[1]`/`[2]`** unchanged (binary parity with released `assign_hap.c`). **`iter_update_var_hap_to_cons_alle`** clears profiles, re-assigns every read in **`ordered_read_ids`** order, rebuilds profiles, and updates consensus; stops early if haps 1–2 consensus is unchanged.
3. **Phase 3:** **`update_read_phase_set`** sets each read’s **`phase_sets[read_i]`** from the first phased het variant in profile order.
4. **Phase 4:** Fill **`hap_alt` / `hap_ref`** on every candidate from finalized **`hap_to_cons_alle`**, including **`-1` → ref** and hom-index fallback when both haps are unknown.

**ONT:** Homopolymer indels use the **67%** read-support guard when updating consensus alleles, as in longcallD **`update_var_hap_to_cons_alle`**.

**Outputs.** Per-candidate **`phase_set`**, **`hap_alt`**, **`hap_ref`** are written to the main TSV (§21). Per-read **`haps`** and **`phase_sets`** live on **`BamChunk`** for downstream use; they are not serialized to the candidate TSV today.

### 18.1 Mid-Free: Releasing Intermediates After K-means (`mid_free_chunk`)

After `collect_var_main` completes for a chunk, `process_chunk` calls **`mid_free_chunk`** to release heavy per-chunk fields before stitching holds all chunks of the contig in memory simultaneously.

This mirrors longcallD **`bam_chunk_mid_free`** (`bam_utils.c`), called from within each parallel worker before `stitch_var_main`. Both tools free the same subset of fields; correspondences are:

| pgPhase field | longcallD field | freed by |
|---|---|---|
| `r.digars` | `bam_read->digars` | `bam_chunk_free_digar` |
| `r.qual` | `bam_read->qual` | `bam_chunk_free_digar` |
| `r.noisy_regions` | `bam_read->noisy_regs` (cgranges) | `bam_chunk_free_digar` |
| `chunk.low_complexity_regions` | `low_comp_cr` | `cr_destroy` |
| `chunk.var_noisy_read_cov_cr` | `var_noisy_read_cov_cr` | `bam_chunk_clear_var_noisy_read_cache` |
| `chunk.var_noisy_read_err_cr` | `var_noisy_read_err_cr` | `bam_chunk_clear_var_noisy_read_cache` |
| `chunk.var_noisy_read_marks` | `var_noisy_read_marks` | `bam_chunk_clear_var_noisy_read_cache` |
| `chunk.var_noisy_read_mark_id` → 0 | `var_noisy_read_mark_id` → 0 | `bam_chunk_clear_var_noisy_read_cache` |
| `chunk.noisy_regions` | `chunk_noisy_regs` + `noisy_reg_to_reads/n_reads` | `cr_destroy` + free |
| `chunk.read_var_cr` | `read_var_cr` | `cr_destroy` |

Fields deliberately **not** freed (matching longcallD `bam_chunk_mid_free`):

| Kept field | Reason |
|---|---|
| `r.alignment` (BAM record) | Future phased-BAM output; longcallD keeps `bam_read->bam_record` |
| `chunk.ref_seq` | Freed later in post-free equivalent; longcallD defers to `bam_chunk_post_free` |
| `chunk.ordered_read_ids` | longcallD iterates it in `update_chunk_read_hap_phase_set1` |
| `chunk.read_var_profile` | longcallD never frees in mid or post_free (only in full `bam_chunk_free`) |

**Memory model.** All chunks for one contig stay in memory concurrently during stitching (§19). `mid_free_chunk` reduces per-chunk footprint before that window so peak RAM scales with **`n_chunks × (candidates + reads[haps/phase_sets/alignment])`** rather than the much larger **`n_chunks × (candidates + reads[all fields])`**.

## 19. Chunk-Boundary Stitching (`stitch_chunk_haps`)

After all chunks in a contig batch are independently phased by k-means (§18), their haplotype assignments are **local**: hap 1 in one chunk may label the same physical haplotype as hap 2 in the next. `stitch_chunk_haps` (`collect_phase.cpp`) corrects this by inspecting reads that span chunk boundaries and flipping assignments where the majority disagree. Phase-set anchors are also merged so that consecutive, consistently oriented chunks form one continuous phase block.

This mirrors longcallD **`stitch_var_main`** / **`flip_variant_hap`** (`collect_var.c`). pgPhase matches the same algorithm with one intentional difference: **qname lookup** instead of j-th-index correspondence (see below).

### 19.1 Inner Function: `flip_chunk_hap(pre, cur)`

For each adjacent pair of chunks on the same contig, `flip_chunk_hap` runs these steps:

**Step 1 — build overlap map from the previous chunk.**
Iterate `pre.down_ovlp_read_i` (reads that extend into the downstream boundary), skip skipped reads, and record `qname → (hap, phase_set)` in a hash map. If the map is empty, stop (no boundary reads on the pre side).

**Step 2 — early exit when either chunk has no candidates.**
If `pre.candidates` or `cur.candidates` is empty, no k-means ran and no flip is needed — stop.

**Step 3 — score boundary reads in the current chunk.**
Iterate `cur.up_ovlp_read_i` (reads that reach into the upstream boundary of `cur`). For each non-skipped read with `hap ≠ 0` that also appears in the pre-map and has `pre_hap ≠ 0`:

```text
if pre_hap == cur_hap:  flip_score -= 1   (consistent — counts against a flip)
if pre_hap != cur_hap:  flip_score += 1   (inconsistent — counts toward a flip)
track max_pre_ps = max over matched reads of pre.phase_sets[ri]
track min_cur_ps = min over matched reads of cur.phase_sets[ri]
```

If `flip_score == 0` (tied or no informative boundary reads), stop.

**Step 4 — apply flip and phase-set merge.**
`do_flip = (flip_score > 0)` (majority of boundary reads were inconsistent).

For every `CandidateVariant v` in `cur` whose `phase_set == min_cur_ps`:
- If `do_flip`: swap `hap_to_cons_alle[1] ↔ [2]`; remap `hap_alt` and `hap_ref` (1↔2).
- If `max_pre_ps >= 0`: rewrite `v.phase_set = max_pre_ps`.

For every read `i` in `cur` whose `phase_sets[i] == min_cur_ps`:
- If `do_flip && haps[i] != 0`: `haps[i] = 3 - haps[i]` (1→2, 2→1).
- If `max_pre_ps >= 0`: `phase_sets[i] = max_pre_ps`.

**`stitch_chunk_haps`** iterates pairs `(chunks[0], chunks[1])`, `(chunks[1], chunks[2])`, … left to right, calling `flip_chunk_hap` for each.

### 19.2 Phase-Set Anchor Semantics

`min_cur_ps` is the earliest phase-set anchor among boundary-spanning reads in the current chunk; `max_pre_ps` is the latest anchor among the matching reads in the previous chunk. After stitching, all variants and reads that carried `min_cur_ps` now carry `max_pre_ps`, effectively **extending the previous chunk's phase block** into the current one and merging them into a single continuous block.

### 19.3 Worked Example — No Flip Needed

```text
Chunks: C0 (chr11:1–500000)  →  C1 (chr11:500001–1000000)

Boundary reads (span 490000–510000):
  readA: C0 hap=1 ps=450000  →  C1 hap=1 ps=500000
  readB: C0 hap=2 ps=450000  →  C1 hap=2 ps=500000

flip_score = -1 - 1 = -2   (both consistent)
do_flip = false
max_pre_ps = 450000,  min_cur_ps = 500000

C1 variants/reads with ps=500000: phase_set rewritten to 450000.
C1 candidates: hap_to_cons_alle unchanged; hap_alt/hap_ref unchanged.
Result: C0 and C1 share one phase block anchored at 450000.
```

### 19.4 Worked Example — Flip Needed

```text
Chunks: C0 (chr11:1–500000)  →  C1 (chr11:500001–1000000)

Boundary reads (span 490000–510000):
  readA: C0 hap=1 ps=450000  →  C1 hap=2 ps=500000   (inconsistent)
  readB: C0 hap=2 ps=450000  →  C1 hap=1 ps=500000   (inconsistent)
  readC: C0 hap=1 ps=450000  →  C1 hap=2 ps=500000   (inconsistent)

flip_score = +1 +1 +1 = +3
do_flip = true
max_pre_ps = 450000,  min_cur_ps = 500000

C1 variants with ps=500000:
  hap_to_cons_alle[1] ↔ [2] swapped
  hap_alt: 1→2, 2→1;  hap_ref: 1→2, 2→1
  phase_set: 500000 → 450000

C1 reads with ps=500000:
  haps: 1→2, 2→1
  phase_sets: 500000 → 450000

Result: C1 labels are coherent with C0; both in one phase block anchored at 450000.
```

### 19.5 Difference from longcallD: qname vs j-th Index

longcallD's `flip_variant_hap` uses a **j-th-index loop**: it assumes the j-th read in `down_ovlp_read_i` corresponds to the j-th read in `up_ovlp_read_i` by BAM-loading order. pgPhase uses **qname lookup**: it hashes `pre.down_ovlp_read_i` reads by query name and looks them up in `cur.up_ovlp_read_i`. This is more robust when reads are pooled from multiple input BAM files or when worker completion order differs.

### 19.6 Integration in the Batch Loop

`stitch_chunk_haps` is called once per contig batch, after all workers have finished `mid_free_chunk` (§18.1), before candidate-table concatenation and TSV/VCF output:

```text
collect_chunk_batch_parallel(batch, ...)     # parallel k-means + mid_free_chunk per chunk
stitch_chunk_haps(batch.chunks)              # sequential boundary flip, left to right
merge_chunk_candidates(batch)               # concat candidate tables
append write_variants_tsv_records(...)       # TSV + optional VCF / read-support
```

This matches longcallD's `kt_pipeline` pipeline-2 stage where `stitch_var_main` runs sequentially after the parallel phase-3 workers. Both tools hold all chunks of one contig in memory during this sequential pass; `mid_free_chunk` makes that feasible by releasing large per-chunk intermediates beforehand.

## 20. Merging Candidate Tables Across Chunks

Each chunk produces a candidate table with **`collapse_fuzzy_large_insertions` already applied inside that chunk** (same as longcallD `collect_all_cand_var_sites`). For each **`reg_chunk_i` batch** (see §4), the pipeline **concatenates** those tables in chunk order and appends rows to the TSV (and optional VCF / read-support) streams. There is **no** second `collapse_fuzzy_large_insertions` on the concatenated list, matching longcallD (candidate-site deduplication is per BAM/region chunk only; chunk-boundary stitching (§19) does not merge candidate sites).

Example:

```text
chunk 1: chr11:1-500000
chunk 2: chr11:500001-1000000

If two fuzzy-equivalent large insertion rows appeared in different chunks (unusual but possible
near boundaries), they remain separate rows in the output, as they would in longcallD’s
per-chunk candidate lists.
```

## 21. TSV Output

The primary output is a TSV file. Each row is a candidate variant with counts and category:

```text
CHROM
POS
TYPE
REF
ALT
DP
REF_COUNT
ALT_COUNT
LOW_QUAL_COUNT
FORWARD_REF
REVERSE_REF
FORWARD_ALT
REVERSE_ALT
AF
CATEGORY
INIT_CAT
PHASE_SET
HAP_ALT
HAP_REF
```

- **`CATEGORY`**: final label after the full longcallD-shaped classify pass (noisy overlap, `LOW_AF`→`LOW_COV` rewrite, post-process, containment filter where applicable). Use this for phasing filters and for VCF projection eligibility (`--vcf-output` / `--phased-vcf-output`) plus `INFO.CAT`.
- **`INIT_CAT`**: first-pass `classify_var_cate` only (matches longcallD’s first `CandVarCate-` block before `After classify var:`). Use for parity checks against `longcallD -V 2` stderr (`scripts/compare_candidates.py --category-stage initial`).

Example row:

```text
chr11  1000  SNP  A  G  30  15  15  0  8  7  9  6  0.5  CLEAN_HET_SNP  CLEAN_HET_SNP  1000  1  2
```

- **`PHASE_SET`**: `CandidateVariant::phase_set` from k-means — the `sort_pos()` of the first variant in the block (`pos` for SNPs, `pos - 1` for indels); **0** if phasing did not assign a set. The value is genomic position only, with no chromosome embedded. The canonical phase block identifier is the `(CHROM, PS)` pair; downstream tools (WhatsHap, GATK) interpret it this way, so the same numeric value on two different chromosomes is not a collision.
- **`HAP_ALT` / `HAP_REF`**: polarized haplotype ids (**1** / **2** for het scaffolding, **3** for hom alt-on-both, **0** when unset), from `assign_hap_based_on_germline_het_vars_kmeans` Phase 4.

The row is still a **candidate** table, not a full diploid VCF genotype: these columns are the **read-clustering scaffold** from §18 with boundary stitching from §19, not full MSA or global assembly phasing. **`collect_output.hpp`** / **`collect_output.cpp`** document I/O contracts (Doxygen): VCF `FILTER` and `INFO.CAT` follow **final** `CATEGORY`; writers may **throw** `std::runtime_error` if the primary BAM header cannot be read or an output path cannot be opened.

For deletions, the reference allele is fetched from the FASTA. For insertions, the TSV uses `REF = "."` and stores the inserted sequence in `ALT`.

The TSV writer emits every candidate present in the final candidate table, including filtered categories such as `LOW_COV`, `STRAND_BIAS`, and `NON_VAR`. The **`CATEGORY`** column is therefore essential for downstream phasing: choose which categories are usable rather than assuming every TSV row is clean.

Typical phase-informative categories are:

```text
CLEAN_HET_SNP
CLEAN_HET_INDEL
```

Potentially useful but not heterozygous-marker categories include:

```text
CLEAN_HOM
NOISY_RESOLVED
```

Categories such as `LOW_COV`, `STRAND_BIAS`, and `NON_VAR` should generally be excluded from ordinary germline phasing markers.

## 22. Optional VCF Output

If `--vcf-output` is provided, the command writes a **final-call projected VCF** (longcallD-style output
surface), not a dump of every candidate TSV row.

SNP representation:

```text
candidate: chr11:1000 A>G
VCF:
  POS = 1000
  REF = A
  ALT = G
```

Insertion representation uses the previous reference base as an anchor:

```text
insertion between 999 and 1000, inserted TGA
anchor base at 999 = C

VCF:
  POS = 999
  REF = C
  ALT = CTGA
```

Deletion representation also uses the previous reference base as an anchor:

```text
delete bases 1000-1002, anchor base at 999 = C, deleted sequence = ATG

VCF:
  POS = 999
  REF = CATG
  ALT = C
```

Only projected germline/noisy-call categories are emitted to VCF. Rows that remain candidate-only
(`LOW_COV`, `LOW_AF`, `STRAND_BIAS`, `NON_VAR`) stay visible in TSV but are omitted from VCF.

For emitted rows, FILTER follows longcallD-style call output:

```text
PASS    emitted variant rows
RefCall reserved for non-variant calls (typically omitted by projection)
LowQual reserved for low-quality filtered calls (typically omitted by projection)
NoCall  reserved for depth=0 calls (typically omitted by projection)
```

`--vcf-output` remains site-level (no FORMAT/sample columns), but emitted sites now follow final-call
projection gates (category + depth/alt-depth) for longcallD parity.

The `INFO` field includes:

```text
END
CLEAN for clean categories
SVTYPE and SVLEN for large indels (|SVLEN| >= min_sv_len, default 30)
DP
REFC
ALTC
LQC
AF
CAT
```

**SVLEN sign convention (VCF spec):** `SVLEN` is the difference in length between ALT and REF (`ALT_len - REF_len`). For deletions REF is longer, so `SVLEN` is negative (e.g. a 82 bp deletion has `SVLEN=-82`). For insertions ALT is longer, so `SVLEN` is positive. This is the standard VCF 4.2 convention; tools like bcftools and VEP expect this sign.

## 22.1 Optional Phased VCF Output (`--phased-vcf-output`)

If `--phased-vcf-output FILE` is given, the pipeline writes a projected phased VCF (same projected site set as
`--vcf-output`) with `GT:PS` FORMAT fields derived from the k-means scaffold (§18). This is the output visible
in files like `pgphase_phased.vcf`.

**GT** (Genotype) encodes the diploid allele assignment at each site:

```text
1|0   hap1 = ALT, hap2 = REF  (phased, hap_alt=1, hap_ref=2)
0|1   hap1 = REF, hap2 = ALT  (phased, hap_alt=2, hap_ref=1)
1|1   homozygous ALT on both   (hap_alt=3)
0/0   homozygous REF           (not emitted by projected phased VCF)
./.   unresolved by k-means    (not emitted by projected phased VCF)
```

The pipe `|` separator signals a phased genotype; the slash `/` separator signals an unphased one. This follows the VCF 4.2 specification.

**PS** (Phase Set) is the anchor coordinate of the phase block. It is set to the `PHASE_SET` value from the k-means output — the `sort_pos()` of the first variant in the block (see §21). It is written as `.` for unphased calls (`0/0` and `./.`). Two records sharing the same `(CHROM, PS)` value were phased together in one k-means run; a change in PS on the same chromosome marks a phase-set break caused by insufficient spanning reads between adjacent het variants.

Example phased VCF record:

```text
chr11  1000  .  A  G  .  PASS  END=1000;CLEAN;DP=30;...  GT:PS  1|0:1000
```

The `GT:PS` fields come from the same `hap_alt`/`hap_ref`/`phase_set` fields that populate the TSV
`HAP_ALT`/`HAP_REF`/`PHASE_SET` columns. `--vcf-output` does **not** include FORMAT/sample columns; only
`--phased-vcf-output` does.

## 23. Optional Read-Support Output

If `--read-support` is provided, the command writes one row per read-candidate observation:

```text
CHROM
POS
TYPE
REF_LEN
ALT
QNAME
IS_ALT
LOW_QUAL
REVERSE
MAPQ
CHUNK_BEG
CHUNK_END
```

Example:

```text
chr11  1000  SNP  1  G  read_42  1  0  1  60  1  500000
```

This means `read_42` supports the alternate allele `G` at `chr11:1000`, the observation is not low-quality, the read is on the reverse strand, and the observation came from the chunk `1-500000`.

For phasing, this file is valuable because it converts the candidate set into read-by-site allele observations.

The read-support output uses the same candidate table as the TSV. Therefore, it can include observations for candidates that are later classified as filtered categories. A phasing stage should combine read-support rows with the candidate category table and select the categories appropriate for the intended phasing model.

## 24. Determinism and Reproducibility

Several implementation choices are designed to make output deterministic:

```text
reads are sorted by start, then end (desc), NM, and query name
candidate sites are sorted by a stable longcallD-compatible comparator
within a batch, each worker writes results into a fixed per-chunk slot (chunk index order)
streaming append order follows the global RegionChunk vector order, batch by batch
```

This means that, for the same input files and options, changing `--threads` should not change the final TSV candidate set or row order. The validation script checks this property by comparing single-threaded and multi-threaded HiFi output (hashed file identity).

## 25. Non-Goals and Current Boundaries

This command intentionally stops at **candidate** collection, classification, the **k-means read-clustering scaffold** (§18), and **chunk-boundary stitching** (§19). It does **not** perform:

```text
somatic noisy-region calling from longcallD's `out_somatic` path
final genotype likelihood (GL) models or production diploid genotyping
supplementary-alignment stitching
split-read structural-variant reconstruction
sample-level FORMAT genotype emission in VCF
```

These boundaries matter for interpretation. A `NON_VAR` or `LOW_COV` candidate in the TSV is not a final biological assertion that no event exists at the locus. It means the collector did not consider that candidate reliable enough under its evidence and classification rules. **`PHASE_SET` / `HAP_*`** summarize the scaffold from §18 and stitching from §19; they do not replace downstream MSA or global phasing when longcallD runs its full caller.

## 26. Summary of the Complete Pipeline

The full implementation can be summarized as:

```text
parse command line (collect_bam_variation in collect_pipeline.cpp)
select read technology mode: HiFi, ONT, or short reads
load reference FAI; build region chunks (parse_region, BED, autosome, split, neighbours)
open TSV / optional VCF / optional read-support streams; write headers
for each reg_chunk_i batch (typically one contig’s chunks):
    parallel for each chunk in the batch:
        for each input BAM/CRAM:
            load overlapping reads
            filter alignments
            build Digars from EQX, cs, MD, or reference comparison
            detect read-level noisy regions
            skip reads with excessive event/noisy density
        pool and sort reads across inputs
        load reference slice
        detect low-complexity reference intervals (sdust)
        merge read-level noisy intervals into chunk noisy intervals
        collect unique SNP/INS/DEL candidate sites (fuzzy large-insertion collapse within chunk)
        count ref/alt/low-quality/strand support
        extend and filter noisy regions by read support (pre_process; longcallD order)
        classify candidates (collect_var.cpp)
        add repeat/dense candidates back into noisy regions
        post-process noisy regions
        mark contained non-ONT candidates as NON_VAR when applicable
        collect_read_var_profile (3.1); assign_hap k-means (3.2, kCandGermlineClean) when candidates non-empty
        mid_free_chunk: release digars, noisy intervals, interval trees, read_var_cr
    stitch_chunk_haps: flip hap labels and merge phase-set anchors at chunk boundaries (sequential)
    concatenate candidate tables for this batch (no second fuzzy collapse; longcallD parity)
    append TSV rows; append optional VCF / read-support rows
close streams; log chunk count and observation totals
```

The central design principle is separation between site discovery and evidence counting. The first pass identifies possible candidate sites from alignment events. The second pass counts how reads support those sites. The classification stage then decides whether each site is clean, low-support, strand-biased, repeat-associated, noisy, resolved, homozygous-like, or non-variant.

This makes the output appropriate for downstream phasing: clean heterozygous SNPs and indels are readily usable as phase-informative markers, while ambiguous loci remain visible through their categories and counts instead of being silently discarded.

**Source documentation:** The `collect_*` and `collect_phase` translation units (`collect_pipeline`, `collect_var`, `collect_phase`, `collect_output`, `collect_types`) carry Doxygen-style `@file` / `@brief` / `@param` / `@return` comments so behavior matches this document at the symbol level. The introduction’s **Vendored `cgranges`** subsection explains the longcallD fork and why it is vendored; **§8** documents parity for the per-read noisy sliding window (`bam_digar` vs longcallD `xid_queue_t`); **§18.1** documents `mid_free_chunk` field-by-field parity with longcallD `bam_chunk_mid_free`; **§19** documents `stitch_chunk_haps` / `flip_chunk_hap` parity with longcallD `stitch_var_main` / `flip_variant_hap`. Re-vendor `src/cgranges.{c,h}` only after an intentional diff against longcallD if upstream changes.

## 27. 2026 Strict Parity Updates (Line-by-Line longcallD Alignment)

This section records the strict parity work that was applied after the initial C++ port so `pgphase`
matches longcallD phased VCF behavior on both ONT and HiFi fixtures.

### 27.1 Candidate/Profile Matching and Multi-Allelic Handling

- `collect_read_var_profile` now follows longcallD matching order for germline read profiling:
  overlap by `ovlp_var_site(...)`, then strict site equality by `exact_comp_var_site(...)`.
- A C++ `ovlp_var_site` helper was ported from longcallD and used in read-profile assignment.
- For this read-profile path, insertion-specific fuzzy matching (`exact_comp_var_site_ins`) is not
  used; strict compare is used for all variant types, matching longcallD behavior at this stage.
- This fixed incorrect per-read allele assignment at same-position multi-allelic insertion sites
  (for example `C->CA` / `C->CAA` separation and downstream hap support).

### 27.2 Haplotype Projection Parity (`hap_to_cons_alle` -> `HAP_ALT/HAP_REF` -> `GT`)

- `collect_phase.cpp` and `collect_output.cpp` now project ALT/REF haplotypes using longcallD logic:
  any non-zero consensus allele index is treated as ALT in hap polarization (`c != 0`).
- `variant_allele_slots(...)` / `get_var_init_max_cov_allele(...)` were aligned to longcallD slot
  sizing semantics: prefer `n_uniq_alles` (when > 0), otherwise use `alle_covs.size()`.
- These changes removed incorrect `0/0` or unresolved outcomes at sites where longcallD emits phased
  ALT genotypes.

### 27.3 Noisy-Region Depth/Count Semantics

- In noisy MSA profile updates (`update_cand_var_profile_from_cons_aln_str*`), `total_cov` now
  increments per full-cover read exactly as longcallD does.
- Depth finalization no longer overwrites a previously collected non-zero noisy `total_cov`.
- This aligned `DP/REFC/ALTC/LQC/AF` fields used in final category gating and VCF INFO.

### 27.4 Homopolymer Indel Veto Parity

- `var_is_homopolymer_indel` in `collect_phase_noisy.cpp` was replaced with a literal longcallD port.
- Insertion homopolymer evaluation uses the VCF anchor convention (`ref_pos-1`) as in longcallD
  downstream behavior.
- This corrected false homopolymer tagging that previously suppressed phasing for valid noisy indels.

### 27.5 Noisy MSA Alignment/Boundary Parity (`align.cpp`)

- The previous approximation path for partial alignment boundaries was removed.
- longcallD functions were ported directly: `edlib_xgaps`, `cal_wfa_partial_aln_beg_end`,
  `collect_partial_aln_beg_end`.
- In both PS and no-PS abPOA paths, cluster/read-id plumbing was aligned to longcallD:
  global read ids are passed through to abPOA helpers, and incorrect local->global remapping was removed.
- One-consensus handling now mirrors longcallD cluster population behavior.

### 27.6 Noisy Region Control-Flow/Coordinate Parity (Final 1 bp Fix)

- `collect_noisy_vars1` call order now matches longcallD:
  1) clip region via `collect_reg_ref_bseq`, then 2) collect noisy-region reads.
- Final step-4 coordinate parity fix:
  `collect_noisy_vars1` now enters noisy calling with longcallD start semantics (`cr_start`), which in
  this C++ interval representation corresponds to `noisy_reg_beg = reg.beg - 1`.
- This removed the remaining 1 bp ONT drift (`1435486` vs `1435485`) in no-PS noisy calling.

### 27.7 VCF Projection and Emission Contract Parity

- Optional VCF outputs are now longcallD-style projected callsets, not candidate dumps.
- Candidate-only rows (`NON_VAR`, `LOW_COV`, `LOW_AF`, `STRAND_BIAS`) are retained in TSV but excluded
  from projected VCF call emission.
- FILTER/INFO semantics were aligned to longcallD call output behavior, including `CAT` and `CLEAN`.
- Left-normalization/rotation was not retained in writer output because longcallD does not do this in
  its VCF writer path.

### 27.8 Compare Script Parity/Validation Updates

- `scripts/compare_phased_vcf.sh` now passes the selected technology mode to `pgphase`
  (`--ont` / `--hifi`) so comparisons run against the intended model.
- Multi-allelic same-position records are keyed with ALT (for non-SVLEN records), preventing false
  collapsed comparisons at shared positions.

### 27.9 Current Parity Status (Fixtures)

Using the repository parity script and test fixtures:

```text
bash scripts/compare_phased_vcf.sh chr11 --ont
  MATCH=400 DIFF=0 ONLY_PG=0 ONLY_LCD=0

bash scripts/compare_phased_vcf.sh chr11 --hifi
  MATCH=385 DIFF=0 ONLY_PG=0 ONLY_LCD=0
```

This is the current strict parity baseline for phased-het `GT:PS` output between `pgphase` and
longcallD on the provided `chr11` ONT/HiFi test datasets.
