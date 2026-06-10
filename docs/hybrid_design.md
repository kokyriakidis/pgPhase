# Hybrid BAM+Graph Phasing: Design Document

## Goal

Combine BAM pipeline's de novo variant calling (83k het sites, 97% accuracy)
with graph pipeline's read coverage (98% reads aligned) to phase more reads
than either pipeline alone.

## Data flow comparison

### BAM pipeline (per chunk)

```
BAM reads → load_read_records_for_chunk
         → collect_candidate_sites_from_records     (de novo site discovery)
         → collect_allele_counts_from_records        (per-read allele counts)
         → classify_chunk_candidates                 (filter/categorize)
         → collect_read_var_profile                  (sparse read×site matrix)
         → assign_hap_based_on_germline_het_vars_kmeans (k-means phasing)
         → collect_noisy_vars_step4                  (MSA recall + 2nd k-means)
```

Key structures:
- `PhasingChunk.reads[]` — ReadRecord with CIGAR-derived digars, positions
- `PhasingChunk.candidates[]` — CandidateVariant with VariantKey, counts
- `PhasingChunk.read_var_profile[]` — ReadVariantProfile (read_id, start/end var idx, alleles[])
- `PhasingChunk.read_var_cr` — cgranges interval tree for read↔candidate overlap

### Graph pipeline (per chunk)

```
Sites VCF → load_sites_for_region                   (load snarl sites)
GAF reads → query_gbz_interval_gaf_ffi / scan_indexed_gaf_chunk
         → build_graph_chunk                         (allele walk matching)
           → add_graph_candidate                     (site → CandidateVariant)
           → dedup read observations
           → depth/AF/strand filtering
           → biallelic decomposition
           → add_read_profile                        (read → ReadVariantProfile)
           → rebuild_read_var_cr
         → classify_graph_candidates                 (categorize)
         → apply_graph_noise_filter                  (homopolymer/repeat/SDUST)
         → assign_hap_based_on_germline_het_vars_kmeans
```

Key difference: graph pipeline doesn't do de novo site discovery — sites come
from the VCF.  Read allele observations come from graph walk matching, not
CIGAR parsing.

### Shared phasing core

Both pipelines produce the same `PhasingChunk` structure.  The k-means phasing
code (`assign_hap_based_on_germline_het_vars_kmeans`) operates identically on
both — it only sees `candidates[]`, `read_var_profile[]`, and `read_var_cr`.

## Hybrid approach

### Per-chunk flow

```
1. BAM variant calling (existing)
   BAM reads → candidate discovery → classification → CandidateTable

2. Graph site augmentation (new)
   Sites VCF → load snarl sites for chunk region
   For each snarl site NOT already in BAM CandidateTable:
     If site has ≥2 alleles and passes basic filters:
       Add as new CandidateVariant at end of table

3. Graph read injection (new)
   GAF reads → query allele observations for chunk region
   For each graph read:
     For each allele observation at a snarl site:
       Find matching CandidateVariant (by position + allele match)
       If found AND read not already in chunk.reads:
         Add ReadRecord (synthetic, with graph-derived positions)
         Build ReadVariantProfile spanning observed candidates

4. Rebuild read_var_cr with augmented reads

5. Unified k-means (existing)
   assign_hap_based_on_germline_het_vars_kmeans on augmented chunk

6. Noisy-region MSA (existing, BAM reads only)
   collect_noisy_vars_step4 — only BAM reads have CIGAR for MSA
```

### Candidate matching

A snarl site matches a BAM candidate when:
- Same GRCh38 position (after normalizing pangenome CHROM → contig suffix)
- Same REF and ALT alleles (exact match)

The 56k exact-match sites from the overlap analysis are the bridgeable sites.
Position-only matches (different alleles) are skipped — allele mismatch would
inject noise.

### Read identity

A graph read may already exist in `chunk.reads` (if it mapped to GRCh38 and
was loaded from the BAM).  To avoid double-counting:
- Build a set of read names from BAM-loaded reads
- Only add graph reads whose names are NOT in the set
- For reads present in both: they already have BAM-derived allele observations;
  graph observations at additional (graph-only) sites extend their profile

### Graph-only candidates

Sites in the graph VCF with no BAM counterpart (~4k het sites) are added as
new candidates.  Their allele counts come entirely from graph read observations.
They participate in k-means like any other candidate.

These sites are valuable because they're in regions where GRCh38 has mapping
issues — exactly where the pangenome helps most.

### Read profiles for graph reads

Graph reads don't have CIGAR-derived digars or base qualities.  Their
ReadRecord is synthetic:
- `tid` = chunk tid
- `beg/end` = min/max candidate positions they observe
- `mapq` = GAF mapping quality
- `digars` = empty (no CIGAR)
- `noisy_regions` = empty (no XID intervals)

Their ReadVariantProfile is built the same way as in `add_read_profile`:
alleles[] spans from first to last observed candidate index, with -1 for
unobserved positions.

### Phase set and haplotype output

After k-means, every read (BAM and graph) has a haplotype assignment.
- BAM reads: written to phased BAM with HP/PS tags (existing logic)
- Graph-only reads: written to a supplementary unaligned BAM with HP/PS tags,
  or appended to the main BAM as unmapped records with HP/PS tags

## Implementation plan

### New files
- `src/hybrid_collect.cpp` — orchestrates the hybrid per-chunk flow
- `src/hybrid_collect.hpp` — declarations

### Modified files
- `src/main.cpp` — add `collect-hybrid-variation` command
- `Makefile` — add new source files

### Reused components
- `collect_pipeline.cpp` — chunking, BAM I/O, chunk stitching, output
- `collect_var.cpp` — `collect_var_main` (BAM variant calling)
- `graph_sites.cpp` — `load_sites_for_region` (snarl site loading)
- `graph_query.cpp` — GAF querying
- `graph_bam_adapter.cpp` — `add_read_profile`, `rebuild_read_var_cr`
- `collect_phase.cpp` — k-means phasing
- `noise_filter.cpp` — reference-based noise detection

### Key functions to implement

1. `find_matching_candidate(CandidateTable, pos, ref, alt) → index or -1`
   Binary search by position, then exact REF/ALT match.

2. `augment_candidates_with_graph_sites(chunk, graph_sites, ref_seq)`
   Add graph-only sites as new CandidateVariant entries.

3. `inject_graph_reads(chunk, graph_rows, site_to_candidate_map)`
   Add graph-only reads and extend profiles for existing reads.

4. `hybrid_process_chunk(region, opts, bam_context, graph_context)`
   Orchestrate: BAM calling → graph augmentation → unified k-means.
