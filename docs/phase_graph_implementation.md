# `phase-graph` implementation description

This document is an **in-depth** description of `pgphase phase-graph`: how tiled **tabix-indexed** graph-site VCF windows and **tabix-indexed pggaf GAF** ranges are joined into `BamChunk` scaffolds, how phasing and stitching reuse `collect_phase.cpp`, and how **streaming** outputs preserve bounded memory on whole-genome runs.

For BAM alignment semantics shared with this path (k-means internals, `stitch_chunk_haps` voting/propagation), see **`docs/collect_bam_variation_implementation.md`**.

---

## 1. Goals and non-goals

**Goals**

- Tile each reference contig into half-open genomic intervals, load only overlapping VCF sites per tile, and scan only overlapping GAF records per tile.
- Represent graph-supported alleles as **`CandidateVariant`** + **`ReadVariantProfile`** rows so **`assign_hap_based_on_germline_het_vars_kmeans`** and **`stitch_chunk_haps`** run **unchanged** from the alignment pipeline (modulo category bitmask).
- Emit a **phased site scaffold TSV** (`PHASE_SET`, hap consensus alleles) and an **unaligned BAM** with optional **`HP` / `PS`** aux tags.
- Record every site dropped before or during het filtering with an explicit **`filter_reason`** (optional TSV + stderr histogram).

**Non-goals**

- `phase-graph` does **not** mirror all graph diagnostics on `collect-bam-variation` (catalog dump, read-profile TSV, variant-style phase dump as separate default outputs, etc.). Those live under **`collect_pipeline.cpp`** and related writers.
- Indexed-GAF mode does **not** shell out to `query` / GBZ here; that subprocess path exists for other commands (`query_gbz_interval_gaf` in `graph_query.cpp`).

---

## 2. Entrypoint and global control flow

**File:** `src/graph_pipeline.cpp`, function **`phase_graph`**.

High-level loop:

1. Parse CLI; validate paths and numeric bounds.
2. Build two **`Options`** bundles:
   - **`filter_opts`** — passed into **`build_graph_bam_chunk`** and GAF scan thresholding: copies **`verbose`**, sets **`touch_read_phase = true`**. All depth / AF / alt-depth / MAPQ thresholds use **`collect_types.hpp`** defaults (`min_mapq`, `min_depth`, `min_af`, `max_af`, `min_alt_depth`, …) because **`phase-graph` does not expose those flags on the CLI**.
   - **`phase_opts`** — passed into **`assign_graph_chunk_hap`** / **`stitch_graph_chunk_pair`**: **`ReadTechnology::Hifi`**, **`output_aln = "graph"`**, **`verbose`**.
3. Open **`--graph-sites-tsv`** and write **`write_graph_bam_phase_sites_tsv_header`** once.
4. Open **`--graph-phase-bam`**, allocate a minimal SAM header (`HD VN=1.6 SO=unknown`), write header.
5. Optionally open **`--graph-filtered-sites`** and write **`write_graph_bam_filtered_sites_tsv_header`**; allocate **`filter_reason_hist`**.
6. **`require_indexed_gaf`** + **`require_graph_site_vcf_tabix_index`** — both inputs must be BGZF-backed with **`.tbi`** (see §5).
7. For each logical contig (either **`--contig`** alone or every name from **`load_graph_site_contig_names`**):
   - Resolve **`[interval_beg, interval_end)`**: from **`--interval`** when exactly one contig is selected; otherwise **`[0, graph_site_contig_end_bp_from_vcf_header(...))`**.
   - **`split_graph_interval`** → vector of half-open chunk spans **`[chunk_beg, chunk_end)`**.
   - Process chunks in **batches of width `threads`** (see §7).
8. After **all** contigs: **`collect_gaf_read_names`** reads the **entire** GAF file (not tabix-sliced) to build the sorted universe of MAPQ-passing QNAMEs; **`write_graph_phase_bam_for_unobserved`** appends BAM records for reads that never entered **`rows_by_read`** / **`emitted_read_names`** (§13).
9. If **`verbose >= 1`**, print filtered histogram and totals.

---

## 3. Coordinate conventions

| Layer | Convention | Notes |
|-------|------------|--------|
| **`--interval`**, tiling **`split_graph_interval`** | **0-based half-open** `[beg, end)` | Parser rejects `beg < 0` or `end <= beg`. |
| **`BamChunk.region.beg` / `region.end`** | **1-based inclusive start**, end consistent with collect pipeline | Set as **`beg + 1`** and **`end`** from the tile’s half-open `(beg, end)` when building the chunk in **`build_graph_bam_chunk`**. |
| **`site_starts_in_interval`** | Compares **0-based** site start **`site_beg0 = (ref_beg > 0 ? ref_beg : site.pos) - 1`** against tile **`[beg, end)`** | A site belongs to **exactly one** chunk: the tile whose half-open interval contains its **start**. Snarls that span chunk boundaries do **not** duplicate the site across chunks. |
| **Tabix catalog filter** (`load_half_open_tabix_catalog`) | **`RegionFilter`** uses **`f.beg = beg0 + 1`**, **`f.end = end0`** | Bridges half-open tiling coordinates into the same filtering helper used elsewhere for VCF extraction. |
| **Indexed GAF `tbx_itr_queryi`** | **`[beg, end)`** half-open | Matches BAM-style interval query on **`rc`/`rb`/`re`** coordinates. |

---

## 4. Module map

| Module | Role |
|--------|------|
| **`graph_pipeline.cpp`** | CLI; per-contig tiling; parallel batch dispatch; **`stable_sort`** batch results by **`region.beg`**; ordered **`stitch_graph_chunk_pair`**; streaming site TSV + BAM flush; histogram; final unobserved pass. |
| **`graph_bam_adapter.hpp` / `.cpp`** | **`GraphBamChunkBuildResult`**, **`build_graph_bam_chunk`**, overlap **`populate_graph_chunk_pair_overlaps`**, **`assign_graph_chunk_hap`**, **`stitch_graph_chunk_pair`**, **`merge_graph_chunk_into_read_rows`**, **`merge_phase_read_assignment`**, **`flush_graph_phase_bam_after_merge`**, TSV/BAM writers. |
| **`graph_query.hpp` / `.cpp`** | **`GraphReadAllele`**, **`scan_indexed_gaf_chunk`** (tabix interval + compact walk matching), **`require_indexed_gaf`**. |
| **`graph_sites.hpp` / `.cpp`** | **`GraphSitesTabixReader::load_half_open_interval`** → **`load_half_open_tabix_catalog`** + **`finalize_graph_site_catalog_inplace`**. |
| **`collect_phase.cpp`** | **`assign_hap_based_on_germline_het_vars_kmeans`** (via adapter), **`stitch_chunk_haps`**. |
| **`collect_types.hpp`** | **`BamChunk`**, **`CandidateVariant`**, **`ReadRecord`**, **`ReadVariantProfile`**, **`Options`**. |

---

## 5. Input contracts

### 5.1 Graph-site VCF

- Must be **bgzip-compressed** with **tabix index** (`.tbi` or `.csi`) — enforced by **`require_graph_site_vcf_tabix_index`**.
- **`GraphSitesTabixReader`** holds an **`htsFile*`** + **`tbx_t*`** per **worker thread** (each worker constructs its own reader inside **`run_vcf_gaf_chunks_parallel`**).
- **`load_half_open_interval(contig, beg0, end0, keep_strings=true)`** for phase-graph passes **`keep_allele_traversal_strings = true`** so walks remain available for indexed-GAF matching **unless** released earlier by other paths (see queryable checks in **`build_graph_bam_chunk`**).

### 5.2 Indexed GAF (pggaf)

- Expected layout: **three leading columns** **`rc` `rb` `re`** (reference contig + 0-based half-open interval), followed by the **standard GAF columns** starting at column index **`3`** — see **`scan_indexed_gaf_chunk`** calling **`scan_gaf_line_compact(..., first_gaf_column = 3, ...)`**.
- Must be **bgzip + tabix** — **`require_indexed_gaf`** opens the path and loads **`tbx_index_load`**.
- **`tbx_seq_tid_with_pangenome_fallback`** resolves **`logical_chrom`** against index SQ lines (supports pangenome-style naming split when needed).

---

## 6. Contig selection and empty-output pitfalls

- With **`--contig NAME`**, only that string drives tiling and **`site_starts_in_interval`** equality against **`site.ref_contig`** (or **`site.chrom`** when ref contig empty).
- With no **`--contig`**, contig list comes from **`load_graph_site_contig_names`** (VCF header / index sequence dictionary).
- **Mismatch:** VCF sites tagged with **`chr20`** but **`--contig CHM13#0#chr20`** yields **no** candidates — sites never pass **`site_starts_in_interval`**. Tabix may still load lines by physical chrom name internally, but the adapter drops them when **`site_contig != contig`**.

---

## 7. Parallelism, batching, and ordering

**Worker pool:** **`run_vcf_gaf_chunks_parallel`** spins **`min(threads, batch_size)`** threads sharing an atomic **`next`** index into the current batch.

Each worker iteration:

1. **`GraphSitesTabixReader reader(sites_vcf)`** — **thread-local** tabix reader (HTS handles not shared across threads).
2. **`load_half_open_interval`** → **`GraphSiteCatalog`** for **`intervals[i]`**.
3. **`scan_indexed_gaf_chunk(gaf, contig, beg, end, catalog, filter_opts.min_mapq)`** → **`std::vector<GraphReadAllele>`**.
4. **`build_graph_bam_chunk(..., chunk_id_offset + i, filter_opts)`** → **`GraphBamChunkBuildResult`**.

**Batch sizing:** Outer loop advances **`batch_beg`** by **`n_per_batch = nw`** where **`nw = max(1, threads)`**. So default **`threads = 1`** processes **one chunk at a time** (still parallel-ready).

**Ordering:** Results return in **completion order**. Before stitching, **`phase_graph` stable-sorts** each batch by **`build_result.chunk.region.beg`**. This restores **genomic order** so **`prev_chunk` ↔ `cur`** adjacency matches physical tiles even when chunk *k* finishes before *k−1*.

**Chunk IDs:** **`chunk_id_offset`** increments by **chunk count per contig** so IDs stay unique across the whole run when multiple contigs are processed.

---

## 8. Site catalog → candidate scaffold (`build_graph_bam_chunk`)

**Output shell:** **`GraphBamChunkBuildResult`** carries:

- **`chunk`** — fully populated **`BamChunk`** (region, candidates, reads, profiles, overlap vectors, `haps`, `phase_sets`, **`read_var_cr`**).
- **`site_ids`** — parallel to **`chunk.candidates`** (after Phase 2 possibly **multiple IDs per original multiallelic site**).
- **`site_allele_orig_idx`** — maps surviving biallelic allele indices back to original walk indices (for ALT strings / diagnostics).
- **`filtered_sites`** — audit rows.
- **`variant_emit_rows`** — CHROM/POS/REF/ALT snapshots for variant-style writers (not the default `phase-graph` outputs).

### 8.1 Site inclusion (pre-candidate)

Iterate **`catalog.sites`**. For each site whose **start** lies in **`[beg,end)`** **and** matches **`contig`**:

- **`precandidate_ineligible`** — **`!site.eligible`** after finalize.
- **`precandidate_not_queryable`** — walk strings not released **and** **`!graph_site_is_queryable(site)`** (indexed GAF requires numeric node walks in compact index).
- **`precandidate_monoallelic`** — fewer than two **`allele_walks`**.

Eligible sites become **multi-allele Phase-1 candidates** via **`add_graph_candidate`** (synthetic **`tid`** from a local `contig → tid` map; single-contig runs effectively use **`tid = 0`**).

**Site identity string:** **`site_key`** prefers **`site.id`** when non-empty and not `"."`; otherwise **`chrom:pos:ref`**.

### 8.2 Parent / conditional-child gating

If **`site.parent`** is set and both parent and child are present as Phase-1 candidates, **`parent_candidate[child] = parent_idx`**.

When aggregating per-read observations, **conditional_parent_alleles** restrict counting: the child observation is counted only if the same read shows the parent allele in the allowed parent-allele set (`std::find` on **`conditional_parent_alleles`**).

### 8.3 From `GraphReadAllele` rows to allele counts

- Index **`site_to_candidate`** maps **`site_id` → Phase-1 candidate index**.
- **`allele_by_read_site`**: for each read and site, keep one allele; **conflicting** duplicate observations for the same read×site → **`allele = -1`** (discarded downstream).
- **`read_obs`**: after conditional-parent filtering, increment **`allele_counts[site][allele]`** and append **`GraphProfileObservation{site_i, allele}`**.

### 8.4 Phase 1 — per-site alt noise trim

Mirror BAM het classification comment in source: drop alt walks with **`count < opts.min_alt_depth`**. Build **`allele_remap[site][old_allele]`** (‑1 = dropped). **Ref walk index 0** always survives. Recompute **`CandidateVariant` counts** from condensed **`alle_covs`**.

Remap **conditional_parent_alleles** through the **parent’s** **`allele_remap`** so child gates stay consistent after parent alleles disappear.

**`allele_orig_idx`** records, for each surviving **new** alt column, which **original** walk index it came from (used for multiallelic **`SITE_ID:walkIdx`** suffixes in Phase 2).

### 8.5 Phase 2 — biallelic decomposition + gates

Each Phase-1 site with **`≥ 2`** surviving alleles is expanded into **one `CandidateVariant` per ref-vs-alt pair**:

- Pair **`SITE_ID`**: if only one alt survives, reuse base **`site_ids[i]`**; if multiple alts, append **`":" + orig_walk_index`** for disambiguation.
- **`FilteredGraphSite`** reasons on rejected pairs:
  - **`ref_only`** — fewer than two alleles after Phase 1 (`ac.size() < 2`).
  - **`low_depth`** — **`ref_c + alt_c < opts.min_depth`**.
  - **`high_af`** — AF **`> opts.max_af`**.
  - **`low_af`** — AF **`< opts.min_af`**.
- Surviving pairs become **`new_cands`** with **`n_uniq_alles = 2`**, updated **`allele_fraction`**, etc.

Parent pointers **`parent_candidate`** / **`conditional_parent_alleles`** are **remapped** into the **new pair index space** by scanning **`old_site_to_pairs`** (first matching parent pair wins).

### 8.6 Phase 3 — remap read profiles to pair indices

For each read’s observations:

- **Biallelic fast path:** map **ref** → allele 0 of the single pair; map matching **alt** → allele 1.
- **Multiallelic:** **ref** observation fans out to **allele 0 on every** surviving pair from that site; **alt *j*** maps to **allele 1** only on the pair whose **`old_alt_phase1`** matches.

Sort observations per read by **`site_index`**, deduplicate consecutive same-site entries with conflicting alleles → **`allele = -1`**, drop negatives, then **`add_read_profile`** builds **`ReadRecord`** + **`ReadVariantProfile`** (interval spans first/last variant **`sort_pos`**).

Finally:

- **`chunk.haps` / `chunk.phase_sets`** sized to **`reads.size()`**, initialized to **`0` / -1**.
- **`rebuild_read_var_cr`** builds **`cgranges`** index over profiles for overlap queries inside phasing.

---

## 9. Per-chunk haplotype assignment

**`assign_graph_chunk_hap`** calls **`assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean)`**.

Only candidates whose category bitmask intersects **`kCandGermlineClean`** participate in clustering; graph-built rows start as **`CleanHetSnp`** / matching **`lcd_var_i_to_cate`**.

No cross-chunk information is used here — stitching happens **after** each chunk has standalone hap labels and **`phase_set`** on **`CandidateVariant`** / reads.

---

## 10. Adjacent-chunk stitching

**Streaming path:** **`stitch_graph_chunk_pair(pre, cur, phase_opts)`**:

1. **`populate_graph_chunk_pair_overlaps`** — match **`ReadRecord.qname`** between **`pre.chunk`** and **`cur.chunk`**; fill **`down_ovlp_read_i`** / **`up_ovlp_read_i`** (batch-of-one overlap slot **`[0]`**).
2. Move both **`BamChunk`**s into **`std::vector<BamChunk> pair`**.
3. **`stitch_chunk_haps(pair, &phase_opts, nullptr)`** — **`nullptr`** disables `.pgbam` sidecar augmentation.
4. Move results back into **`pre.chunk`** and **`cur.chunk`**.

**Batch alternative:** **`phase_graph_bam_chunks`** (used elsewhere) assigns hap for every chunk, **`populate_graph_chunk_overlaps`** on the full vector, then **one** **`stitch_chunk_haps`** over all chunks — equivalent stitching semantics, different scheduling.

---

## 11. Streaming site TSV (`--graph-sites-tsv`)

**Header** (`write_graph_bam_phase_sites_tsv_header`):

`CHUNK_ID`, `SITE_INDEX`, `SITE_ID`, `POS`, `N_ALLELES`, `DEPTH`, `ALLELE_COUNTS`, `PHASE_SET`, `HAP1_ALLELE`, `HAP2_ALLELE`.

**When rows are written**

- After **`stitch_graph_chunk_pair(prev, cur)`**, emit **`write_graph_bam_phase_sites_tsv_rows(sites_out, *prev_chunk)`** — the **previous** tile is now **final** through its downstream boundary.
- At **contig end**, emit **`write_graph_bam_phase_sites_tsv_rows`** for the **last** held chunk.

Thus interior chunks are emitted **once**, **after** stitching with their successor; the terminal chunk is emitted **without** a following stitch.

**Content:** Rows reflect **`candidate.phase_set`** and **`candidate.hap_to_cons_alle[1/2]`** after **`stitch_chunk_haps`** updated **`BamChunk`** state.

---

## 12. Per-read accumulator (`PhaseReadOutputRow`)

Fields (**`graph_bam_adapter.hpp`**):

- **`read_name`**, **`chunk_id`** / **`copies`** — diagnostics for multi-chunk appearances; **`chunk_id` becomes `-1`** if the read span touches multiple chunk IDs during merges.
- **`hap`**, **`phase_set`**, **`has_phased_assignment`** — latest stitched assignment (§13).
- **`allele_by_site`** — merged observations **`merge_phase_read_observation`**: first allele wins; conflicting later allele → **`-1`** for that site key.

**`merge_graph_chunk_into_read_rows`** walks **`read_var_profile`**:

- For each site observation with **`allele >= 0`**, merge **site_id + allele**.
- **`merge_phase_read_assignment`** — **always overwrites** **`hap` / `phase_set`** from the **current** chunk’s **`ReadRecord` index**, and sets **`has_phased_assignment`** iff **`hap ∈ {1,2}`** and **`phase_set ≥ 0`**.

This tracks **`PhasedAlignmentWriter::write_chunks`** semantics: HP/PS reflect the **owning** chunk at emission, which for overlapping tiles is the **later** (genomic downstream) chunk after stitching.

---

## 13. Streaming phased BAM (`--graph-phase-bam`)

### 13.1 Record shape

**`write_bam_record`** emits:

- **`FLAG = FUNMAP`** (unaligned),
- **No RNAME/POS** (`tid = -1`, `pos = -1`),
- Optional **`HP:i`** and **`PS:i`** aux tags when phased.

### 13.2 Incremental flush

After merging **`prev_chunk`** into **`rows_by_read`**, **`flush_graph_phase_bam_after_merge`**:

- If **`next_chunk_qnames`** is non-null, **retain** rows whose **`read_name` appears in the next chunk’s read set** (they may receive a fresher **`merge_phase_read_assignment`** later).
- **Drain** all others: sort drained keys lexicographically for deterministic output order, write BAM, insert into **`emitted_read_names`**.

**End of contig:** call with **`next_chunk_qnames = nullptr`** to flush **all** remaining reads for that contig.

### 13.3 Unobserved reads

**`collect_gaf_read_names`** scans the **whole** GAF (supports both column layouts offset **`0`** or **`3`** for pggaf-prefixed lines), collects MAPQ-filtered unique names, sorts.

**`write_graph_phase_bam_for_unobserved`** appends **unmapped** records **without HP/PS** for names never emitted — reads that never produced graph observations merged into **`rows_by_read`**.

---

## 14. Filtered-site audit and verbosity

**Optional `--graph-filtered-sites`:** for each chunk, **`write_graph_bam_filtered_sites_tsv_rows`** dumps **`filtered_sites`**:

`site_id`, `ref_cov`, `alt_cov`, `total_cov`, `allele_fraction`, `filter_reason`.

**`-V 1` stderr:**

- Per-contig chunk count and **sum of `read_allele_rows`** (GAF-derived row count).
- Line **`chunk-local phased site row(s)`** — note this sums **`site_ids.size()`** per chunk across chunks; it is **not** the same as globally distinct scaffold rows (sites can only appear in one chunk by construction, but the metric is literal per-chunk candidate counts summed).
- Global histogram: **`filter_reason` → count**, sorted by descending count. Total filtered rows must equal lines in the audit file (excluding header).

---

## 15. Error propagation

**`run_vcf_gaf_chunks_parallel`** catches worker exceptions under mutex and **`std::rethrow_exception`** after **`join`** — first failure wins.

---

## 16. Tests and parity tooling

| Artifact | Coverage |
|----------|----------|
| **`src/test_graph_bam_adapter.cpp`** | Chunk build, AF thresholds, phase-site / phase-read writers (batch-style). |
| **`src/test_graph_sites.cpp`**, **`src/test_graph_phase.cpp`** | Catalog / matrix helpers. |
| **`test_phase_block_stitch`** | **`flip_chunk_hap`** / propagate semantics depended on by **`stitch_chunk_haps`**. |

---

## 17. Quick reference — pipeline ordering inside `phase_graph` inner loop

For each **`cur`** chunk result (already sorted within batch):

1. **`assign_graph_chunk_hap(cur, phase_opts)`**
2. Accumulate **`filter_reason_hist`**; optionally append filtered TSV rows
3. If **`prev_chunk`** exists:
   - **`stitch_graph_chunk_pair(*prev_chunk, cur, phase_opts)`**
   - **`write_graph_bam_phase_sites_tsv_rows`** for **`prev_chunk`**
   - **`merge_graph_chunk_into_read_rows(rows_by_read, *prev_chunk)`**
   - Build **`next_qnames`** from **`cur.chunk.reads`**
   - **`flush_graph_phase_bam_after_merge(..., &next_qnames, emitted)`**
4. **`prev_chunk = move(cur)`**

After all batches on contig: write **last** chunk sites, merge last chunk rows, **`flush(..., nullptr)`**.

---

## 18. Related documentation

- **`docs/collect_bam_variation_implementation.md`** — full collect pipeline; stitching § references **`stitch_chunk_haps`** behaviour in detail.
