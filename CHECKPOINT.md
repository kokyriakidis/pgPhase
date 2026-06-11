# pgphase Evaluation Checkpoint

## Test Setup

**Sample:** HG002 chr20 HiFi reads  
**Reference:** CHM13 T2T v2.0 (contig `CHM13#0#chr20`, 66.2 Mbp)  
**Truth:** HG002 v1.1 diploid assembly (chr20 slices only — `HG002_1.1_MATERNAL.chr20.fasta` / `HG002_1.1_PATERNAL.chr20.fasta`)  
**Evaluator:** diplinator (haplotype-aware read assignment) + pgphase accuracy script

### Key files

| File | Description |
|------|-------------|
| `test_data/HG002.chr20.fq.gz` | HiFi reads (272,016 total) |
| `test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam` | Reads aligned to CHM13 chr20 |
| `test_data/chm13v2.0.chr20.renamed.fa` | CHM13 chr20 reference (contig `CHM13#0#chr20`) |
| `test_data/chr20.sites.striped.vcf.gz` | 977,275 graph sites (AT field with numeric node IDs) |
| `test_data/HG002.chr20.annotated.coord.gaf.gz` | bgzipped + tabix-indexed GAF for graph pipeline |
| `test_data/chr20.phase_graph.bam` | Phased BAM from graph pipeline |
| `test_data/bam_eval_chr20/phased.bam` | Phased BAM from BAM pipeline |

---

## Results — Graph Pipeline (indexed GAF)

**Command:**
```
pgphase collect-graph-variation \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --sites test_data/chr20.sites.striped.vcf.gz \
    --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
    --phased-bam-out test_data/chr20.phase_graph.bam \
    -o test_data/phase_graph_smoke/reads.tsv \
    --phased-vcf-out test_data/phase_graph_smoke/phased.vcf \
    -t 16
```

> Note: `--sites-vcf` → `--sites`, `--out-bam` → `--phased-bam-out`, `--ref` now required (commit `e751d83`).

**Eval dir:** `test_data/diplinator_eval_chr20/`

### v1 — before noise filter

| Metric | Value |
|--------|-------|
| Total reads | 231,382 |
| Phased reads | 226,831 (98.0%) |
| Phase sets | 421 |
| Phase block N50 | 912 Kbp |
| Largest block | 2.1 Mbp |
| Overall accuracy | 89.73% |
| Hamming error rate | 10.27% |
| Switch error rate | 6.14% |
| Perfect phase sets | 78 / 411 (19.0%) |
| Phaseable PS (>60%) | 308 PS, 207,571 reads, 93.00% |
| Unphaseable PS (≤60%) | 103 PS, 19,232 reads |

### v2 — with reference-based noise filter (commit `e751d83`)

| Metric | Value |
|--------|-------|
| Total reads | 231,382 |
| Phased reads | 216,130 (93.4%) |
| Phase sets | 425 |
| Phase block N50 | 898 Kbp |
| Largest block | 2.1 Mbp |
| Overall accuracy | 92.77% |
| Hamming error rate | 7.23% |
| Switch error rate | 1.04% |
| Switchflip error rate | 2.51% |
| Switch errors | 2,251 |
| Flip errors | 3,167 |
| Switchflip errors | 5,418 |
| Switch opportunities | 215,693 |
| Perfect phase sets | 96 / 419 (22.9%) |
| Phaseable PS (>60%) | 327 PS, 202,396 reads, 95.36% |
| Unphaseable PS (≤60%) | 92 PS, 13,716 reads |

---

## Results — BAM Pipeline (no pgbam sidecar)

**Command:**
```
pgphase collect-bam-variation \
    --hifi \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    --out-bam test_data/bam_eval_chr20/phased.bam \
    --phased-vcf-out test_data/bam_eval_chr20/phased.vcf \
    -o test_data/bam_eval_chr20/reads.tsv \
    -t 16
```

> Note: `--ref` and `--bam` named options added in commit `1614581`.

**Eval dir:** `test_data/bam_eval_chr20/`

| Metric | Value |
|--------|-------|
| Total reads | 272,016 |
| Phased reads | 220,377 (81.0%) |
| Phase sets | 431 |
| Phase block N50 | 937 Kbp |
| Largest block | 53.2 Mbp |
| Overall accuracy | 97.02% |
| Hamming error rate | 2.98% |
| Switch error rate | 0.62% |
| Switchflip error rate | 1.45% |
| Switch errors | 1,356 |
| Flip errors | 1,826 |
| Switchflip errors | 3,182 |
| Switch opportunities | 219,946 |
| Perfect phase sets | 133 / 430 (30.9%) |
| Phaseable PS (>60%) | 366 PS, 214,808 reads, 98.14% |
| Unphaseable PS (≤60%) | 64 PS, 5,568 reads |

---

## Head-to-Head Comparison

| Metric | Graph v1 (no filter) | Graph v2 (noise filter) | BAM (no sidecar) |
|--------|---------------------|------------------------|-----------------|
| Phased reads | 98.0% | 93.4% | 81.0% |
| Phase sets | 421 | 425 | 431 |
| N50 | 912 Kbp | 898 Kbp | **937 Kbp** |
| auN | — | 981 Kbp | **9,828 Kbp** |
| Accuracy | 89.73% | 92.77% | **97.02%** |
| Hamming error rate | 10.27% | 7.23% | **2.98%** |
| Switch error rate | 6.14% | 1.04% | **0.62%** |
| Switchflip error rate | — | 2.51% | **1.45%** |
| Switch errors | — | 2,251 | **1,356** |
| Flip errors | — | 3,167 | **1,826** |
| Switchflip errors | — | 5,418 | **3,182** |
| Switch opportunities | — | 215,693 | 219,946 |
| Perfect PS | 19.0% | 22.9% | **30.9%** |

**Key observations:**
- Noise filter improves graph accuracy +3 pp (89.7% → 92.8%) and cuts switch error 6× (6.14% → 1.04%)
- Noise filter costs ~4.6% phased reads (noisy indels excluded from k-means)
- BAM pipeline leads on all accuracy metrics: Hamming (2.98%), switch (0.62%), switchflip (1.45%)
- BAM auN (9.8 Mbp) far exceeds graph (981 Kbp), driven by the large 53 Mbp block
- Graph pipeline phases more reads (93% vs 81%), compensating for lower per-block accuracy

---

## Site Overlap Analysis: Graph Sites vs BAM Het Variants

Measured with `scripts/analyze_site_overlap.py` on chr20 data.

| Metric | Value |
|--------|-------|
| Graph sites (total positions) | ~977K |
| BAM het variants | 83,013 |
| Shared positions | 62,600 (75.5% of BAM het) |
| Exact allele match | 56,590 (90.4% of shared) |
| Position-only match (diff alleles) | 6,010 (9.6% of shared) |
| BAM-only (no graph site) | 20,413 (24.5% of BAM het) |

### Interpretation

75.5% of BAM het sites have a graph snarl site at the same GRCh38 position, and 90.4% of
those are exact allele matches. This means graph reads traversing those snarls carry the
same allele information that the BAM pipeline uses for phasing — they just come from
pangenome-aligned reads that may not map to GRCh38.

The 24.5% BAM-only sites are likely private variants or rare alleles not represented in the
pangenome reference panel. These sites can only be phased using BAM-mapped reads.

### Hybrid bridging approach

**Goal:** Combine BAM pipeline's variant calling (83k het sites, 97% accuracy) with graph
pipeline's read mapping (98% reads aligned, including reads unmappable on GRCh38).

**Approach:** Run BAM pipeline as primary caller. For the ~19% unphased reads, use their
graph allele observations at the ~56k bridgeable sites to assign haplotypes and extend
phase blocks.

**Implementation phases:**

1. **Read haplotype assignment** — For each GAF read not in the BAM phased output, find its
   graph allele observations at bridgeable snarl sites. Look up the BAM haplotype consensus
   at those sites. Assign the read to the most consistent haplotype.

2. **Enhanced chunk stitching** — Use newly-phased reads as connectors between BAM phase
   blocks. Two blocks that couldn't be stitched (no shared GRCh38-mapped reads) might now
   be connected through pangenome-located reads that span both regions. The pangenome gives
   read locality information that GRCh38 can't — two reads unmapped on GRCh38 might be
   clearly co-located in the pangenome (traversing the same snarls in the same region).

**Expected gains:**
- Phase the missing ~19% of reads (currently unphased by BAM pipeline)
- Improve N50 by connecting fragmented BAM phase blocks through pangenome-bridged reads
- Maintain BAM pipeline's 97% accuracy (graph reads are assigned to existing haplotypes,
  not used to re-phase)

**Risks:**
- Allele matching imprecision at the 9.6% position-only sites could inject noise
- Reads unmapped on GRCh38 may have ambiguous pangenome positions (multi-mapping)
- Phase set merging logic needs to handle cross-pipeline block boundaries

---

## Merge Feasibility: BAM + Graph Phased BAMs

Measured with `scripts/merge_feasibility.py` on chr20 data.

### Read category breakdown

| Category | Reads | % of total |
|----------|-------|-----------|
| A) Phased by both | 208,092 | 76.5% |
| B) BAM-only phased | 12,285 | 4.5% |
| C) Graph-only (rescue candidates) | 8,038 | 3.0% |
| D) Neither phased | 43,601 | 16.0% |

### Haplotype agreement (category A — doubly-phased reads)

| Metric | Value |
|--------|-------|
| Doubly-phased reads | 208,092 |
| Same HP label | 128,104 (61.6%) |
| Different HP label | 79,988 (38.4%) |
| Overlapping block pairs | 639 |
| Concordant pairs (≥90% same, ≥2 reads) | 249 — 118,177 reads |
| Flipped pairs (≥90% opposite, ≥2 reads) | 169 — 67,159 reads |
| Discordant pairs (<90% majority) | 197 — 22,732 reads |

### Rescue potential

| Metric | Value |
|--------|-------|
| BAM phased reads | 220,377 (81.0%) |
| After graph rescue (potential) | 228,415 (84.0%) |
| Reads rescued | +8,038 (+3.0 pp) |
| Remaining unphased | 43,601 (16.0%) |
| Graph-only in BAM (unphased) | 8,038 |
| Graph-only not in BAM at all | 0 |

### Assessment

**RISKY** — 197 / 615 block pairs (32%) are discordant (neither clearly same nor flipped orientation). A simple polarity-based merge would inject errors at those block boundaries. The 43,601 D-category reads (16%) are absent from both pipelines and represent the hard ceiling for any hybrid approach.

The concordant+flipped pairs do cover the majority of doubly-phased reads (185,336 reads, 89%), suggesting that orientation *is* resolvable for most large blocks — the discordant pairs are likely small blocks with insufficient read overlap. A read-count-weighted merge may still be feasible.

---

## VCF-Level Metrics (chr20, no whatshap/hiphase needed)

Phase block stats extracted directly from phased VCFs; NC50 computed via `bedtools subtract`
of `switch_positions.bed` from phase block spans.

| Metric | Graph v2 (noise filter) | BAM (no sidecar) |
|--------|------------------------|-----------------|
| Het variants in VCF | 60,304 | 83,013 |
| Phased het variants | 60,304 (100%) | 83,013 (100%) |
| Phase blocks (VCF) | 427 | 431 |
| VCF N50 | 426 Kbp | 416 Kbp |
| NC50 (switch-corrected N50) | 167 Kbp | **214 Kbp** |
| Corrected phase blocks | 5,660 | **3,688** |
| Total corrected span | 50.6 Mbp | **55.2 Mbp** |

**Notes:**
- NC50 = N50 of phase blocks after subtracting switch-error positions (analogous to HiPhase paper NGC50 but at chr20 scale)
- Graph calls fewer het variants (60K vs 83K) because the graph pipeline uses graph snarls as sites, not de-novo variant calling from alignments
- BAM's lower corrected block count (3,688 vs 5,660) confirms fewer switch errors chopping up its blocks

---

## Noise Filtering: BAM vs Graph Pipeline

### What we ported

The BAM pipeline's `classify_variant_initial` applies three reference-context checks to
indel candidates:

1. **Homopolymer context** — indel sits in a run of ≥3 copies of a 1–6 bp repeat unit
2. **Tandem repeat context** — deleted/inserted motif matches ≥3 tandem copies in flanking reference
3. **SDUST low-complexity** — variant position falls inside a low-complexity interval

Indels flagged by (1) or (2) become `RepeatHetIndel` and are excluded from k-means phasing
(`lcd_var_i_to_cate = kLongcalldRepHetVar`, which doesn't match `kCandGermlineClean`).

These three checks were extracted into a shared module (`src/noise_filter.hpp/cpp`) and wired
into the graph pipeline via `apply_graph_noise_filter` in `graph_bam_adapter.cpp`. Each worker
thread opens its own `faidx_t`, fetches the chunk's reference slice, and reclassifies
`CleanHetIndel` → `RepeatHetIndel` for noisy sites. The flag behavior matches the BAM pipeline
exactly: `RepeatHetIndel` candidates are excluded from k-means.

### What we did NOT port (and why)

The BAM pipeline has additional filtering stages that depend on **per-read CIGAR error
intervals** — information that doesn't exist in the graph pipeline's GAF input:

| BAM feature | Why it doesn't apply to graph |
|---|---|
| Read-level noisy regions (XID clustering from CIGAR) | GAF allele observations are node-level, not base-aligned. No per-read error intervals. |
| Dense overlap check (`var_pos_cr` n_ov > 1) | Snarl parent-child gating already handles nested/overlapping sites. |
| Noisy-reads-ratio gate | Requires per-read CIGAR error intervals. |
| `post_process_noisy_regs` (widen noisy spans) | Depends on read-level noisy regions. |
| `apply_noisy_containment_filter` (demote contained candidates) | Depends on widened noisy regions. |
| Noisy-region MSA recall (`collect_noisy_vars_step4`) | Requires BAM reads for MSA realignment. |
| Second k-means with `kCandGermlineVarCate` | Only useful after MSA recall produces `NoisyCandHet`. |

### SNPs in low-complexity regions

`is_noisy_site` checks `pos_in_low_complexity` for all variant types, including SNPs. However,
`apply_graph_noise_filter` only reclassifies `CleanHetIndel` candidates — SNPs are not touched.
This matches the BAM pipeline, where SDUST low-complexity intervals are **never used to directly
reclassify any variant**. Instead, they serve a structural role:

1. **Extending noisy region boundaries** — adjacent low-complexity intervals are absorbed into
   read-level noisy spans during `pre_process_noisy_regs_pgphase`.
2. **Extending repeat-indel spans** — when a `RepeatHetIndel` is added to the noisy region tree,
   its genomic span is widened to include overlapping low-complexity intervals.

A SNP can end up filtered in the BAM pipeline, but only **indirectly**: if the low-complexity
region was already part of or adjacent to a noisy region built from read-level error clustering,
and the widened noisy span fully contains the SNP, then `apply_noisy_containment_filter` demotes
it to `NonVariant`. A SNP in a low-complexity region with no nearby read-level noise stays
`CleanHetSnp`.

If we reclassified graph SNPs in low-complexity regions, we'd be **more aggressive** than the
BAM pipeline — potentially removing legitimate het SNPs that are reliable phasing anchors.

### Future work

If the graph pipeline gains access to base-level alignment information (e.g. via surjected BAM
or GAF CIGAR strings), the read-level noisy region pipeline could be ported. Until then, the
reference-context checks are the only applicable noise filter.

---

## Hybrid Model — Phase 0 Baseline (measurement harness)

Phase 0 of the hybrid-model plan: build a measurement harness and establish a
baseline for the existing `collect-hybrid-variation` command before changing it.

### Deliverables

- `scripts/bench_hybrid.sh` — runs hybrid mode and, with `--compare`, also the
  BAM-only and graph-only pipelines on the same input into sibling dirs for A/B.
  Per-pipeline BAM flag differs: graph uses `--phased-bam-out`, bam/hybrid use
  `--out-bam`.
- `scripts/phase_block_stats.py` — VCF-level metrics (het/hom, phased fraction,
  block count, N50) from a phased VCF alone. No truth set needed. For
  switch/accuracy use `evaluate_phase_accuracy.sh` (needs diplinator).

### Environment note

The chr20 HG002 dataset and diplinator from the original eval are **not present**
in this environment, so the full switch/accuracy gate from the plan could not be
reproduced here. Phase 0 instead used the checked-in `test_data/graph/` fixture
(MICB/KIR3DL1, chr6+chr19) by reconstructing reads from the GAF `cs` tags and
aligning them with minimap2 to produce a matched BAM+sites+GAF triple.

> **Caveat:** the fixture is ultra-low-depth (~0.04× mean, 151 bp short reads),
> not the HiFi/ONT long-read target. Numbers below validate the *harness and code
> behavior*, not phasing quality. The real gate still requires the chr20 data.

### A/B result (fixture)

| Pipeline | Candidates | Phased het | Blocks | N50 |
|---|---|---|---|---|
| BAM-only | 1 | 1 | 1 | 0 |
| Graph-only | 292 | 279 | 52 | 1062 |
| Hybrid | 918 | 68 | 17 | 146 |

### Findings (confirm the static analysis)

- **`bridged=0` on every chunk.** No graph site matched a BAM candidate, so the
  intended bridging never fired. Partly because BAM found ~1 candidate at this
  depth, but it also exercises the `find_matching_candidate` sort-assumption gap.
- **Hybrid phases *fewer* het variants (68) than graph-only (279)** despite adding
  917 graph-only candidates. The unconditionally-`CleanHet` graph candidates
  (no AF/depth/het gate) flood the unified k-means and degrade it rather than
  augment it — exactly the Phase 1 hazard.
- **Verbose stats only print when `added>0` AND a worker logs them**; counts are
  per-chunk and the same `added=` value repeats across overlapping chunks,
  suggesting overlap-region sites are added in both chunks (double-add).

### Decision gate

The intended "BAM-primary, graph-rescue" behavior is **not** what the code does
(it co-phases and can regress below graph-only). This confirms the plan's
direction: proceed to **Phase 1** (graph candidate quality gates) and the
**Phase 2** BAM-frozen rescue rather than tuning the current co-phasing path.
The chr20 switch-error gate must still be run once that data is available before
any hybrid path is declared shippable.

---

## Hybrid Model — Phase 1 (graph candidate quality gates)

Phase 1 hardened graph-only candidate handling so graph sites can no longer
enter k-means as unconditional clean-het anchors.

### Changes

- **Deferred classification** (`add_graph_only_candidate`): graph sites are now
  added unclassified (category `LowCoverage`, flag 0) instead of stamped
  `CleanHet`. They stay out of k-means until gated.
- **Quality gate** (`classify_graph_only_candidates`): after counts are final
  (post `inject_graph_reads`), each graph-only candidate is run through the BAM
  pipeline's `classify_variant_initial` and `category_to_flag`. Only sites
  passing the het band (`min_depth`, `min_alt_depth`, `min_af ≤ AF ≤ max_af`)
  become `CleanHet` and enter het k-means. `classify_variant_initial` was
  exposed in `collect_var.hpp` for reuse (was file-static).
- **Ordering fix** (`process_chunk_hybrid`): the indel noise filter
  (`apply_hybrid_noise_filter`, which only acts on `CleanHetIndel`) now runs
  *after* the gate, not before. Verbose log adds `promoted=N`.
- **Sort-precondition fix** (`find_matching_candidate`): the bridge binary
  search now takes an explicit `n` = original BAM candidate count and only
  searches `candidates[0, n)`. `inject_graph_sites` asserts that range is
  sorted, so a future ordering change fails loudly instead of silently missing
  bridges.
- **Unit test** `src/test_hybrid_inject.cpp` (wired into `make unit-tests`):
  clean het promoted; homozygous, low-depth, low-alt-depth, and low-AF sites
  excluded from the het mask. All unit tests green.

### A/B result (same fixture as Phase 0)

| Pipeline | Candidates | Phased het | Promoted |
|---|---|---|---|
| BAM-only | 1 | 1 | — |
| Graph-only | 292 | 279 | — |
| Hybrid (Phase 0) | 918 | 68 | (ungated) |
| **Hybrid (Phase 1)** | 918 | **2** | **2 / 917** |

### Verification

- BAM-only path unchanged from Phase 0 (386 het / 182 hom / N50 194368) — no
  regression from exposing `classify_variant_initial`.
- `make -j` builds with zero warnings; `make unit-tests` all green.
- Asserts are active (no `-DNDEBUG`), so the sort precondition is enforced.

### Correction: the "AF=0 bug" was a synthetic-data artifact

The low-depth fixture suggested hybrid allele counting was broken (graph het
sites showing REFc=N, ALTc=0, AF=0). **This did not reproduce on real chr20
data** (see next section) — there, the same sites show correct AF≈0.4–0.5. The
AF=0 came from the Phase 0 read-reconstruction harness (`gaf_to_fastq.py` +
minimap2 short-read alignment did not carry variant alleles faithfully), not
from the hybrid pipeline. Lesson: validate hybrid behavior only on real aligned
reads, not reconstructed ones.

---

## Hybrid Model — chr20 25M validation (REAL data)

Real HG002 HiFi test data was added at `test_data/graph_chr20/`: a 500 kb
chr20 slice (25,000,001–25,500,000) with a matched BAM, coordinate GAF (2,048
reads), sites VCF (6,963 snarl sites), and N-padded ref. This is the first
hybrid run on genuine aligned reads.

> Note: `samtools index` the BAM first (`.bam.bai` is gitignored, not shipped).

### Expected-output assertions (both pass with Phase 1 build)

| Pipeline | Candidates | Expected |
|---|---|---|
| BAM-only | 1108 | 1108 ✓ |
| Graph-only | 509 (13,490 filtered) | 509 ✓ |

Phase 1 changes preserve both verified pipeline outputs exactly.

### Hybrid A/B (Phase 1)

| Pipeline | Phased het | Blocks | N50 | Largest block |
|---|---|---|---|---|
| BAM-only | 527 | 2 | 296,271 | 296,271 |
| Graph-only | 446 | 2 | 393,964 | 393,964 |
| **Hybrid** | **522** | **1** | **498,281** | 498,281 |

Hybrid bridging works on real data: **791 sites bridged**, 13,135 graph-only
sites added, **72 promoted** to clean-het by the gate, 1,692 BAM reads extended
with graph observations (`reads_injected=0` — the GAF and BAM are the same HG002
reads, so there are no graph-only reads to inject here).

Result: hybrid phases nearly as many het as BAM (522 vs 527) but collapses the
two BAM blocks into **one block with N50 498 kb** — larger than either
standalone pipeline. This is the pangenome-bridging stitching benefit the plan
predicted, appearing for the first time on real data.

### Gate behavior on real data (sanity)

Graph-only candidate categories after gating: 5,875 LOW_COV, 598 CLEAN_HOM,
397 CLEAN_HET_SNP, 88 NOISY_HET, 76 REP_HET_INDEL, 38 NOISY_HOM,
37 CLEAN_HET_INDEL, 32 LOW_AF. Only the clean-het sites enter het k-means; hom,
low-cov, low-AF, and repeat sites are correctly excluded. Spot-check vs graph
pipeline at chr20:25001841 — graph AF 0.52, hybrid AF 0.36 (REFc=7/ALTc=4),
both CLEAN_HET_SNP: alleles are counted correctly.

### Open items

- The **switch-error / accuracy** gate against the diploid truth still needs
  `diplinator` (not in this environment). The N50 gain is promising but must be
  confirmed not to come with a switch-error regression before shipping.
- Phase 2 (BAM-frozen rescue) is most valuable where the GAF contains reads
  absent from the BAM; this fixture has none, so that path needs a dataset with
  graph-only reads to exercise.

---

## Hybrid Model — Phase 3 (stitch abstain margin)

Phase 3 adds an abstain rule to chunk stitching so weakly supported block
boundaries are not merged.  This is the mitigation for the CHECKPOINT "32%
discordant pairs" over-merge risk: a wrong merge injects a switch error, so
when the evidence is thin it is safer to leave two blocks.

### Change

- `flip_chunk_hap` now merges adjacent chunks only when
  `|flip_hap_score| > opts.stitch_min_margin`.  The previous rule
  (`flip_hap_score == 0 → no merge`) is exactly `margin = 0`, the default.
- `Options::stitch_min_margin` (default `kDefaultStitchMinMargin = 0`).
- CLI flag `--stitch-min-margin INT` added to `collect-hybrid-variation`
  **only**.  The BAM and graph commands do not expose it and keep margin 0, so
  their behavior is unchanged.

### Safety: shared code, default-preserving

`flip_chunk_hap`/`stitch_chunk_haps` is shared by all pipelines.  Because the
default margin is 0 and the abstain test reduces to the original `== 0` check,
BAM and graph outputs are byte-identical (re-verified: 1108 / 509 candidates,
527/446 het, same blocks/N50 as Phase 1).

### chr20 margin sweep (real data)

Default chunk size (500 kb) puts the whole 500 kb data slice in one chunk, so
stitching is not exercised.  Forcing `--chunk-size 100000` (≈5 data chunks):

| Margin | Blocks | N50 |
|---|---|---|
| 0 | 2 | 406,166 |
| 1–2 | 2 | 406,166 |
| 4 | 3 | 305,381 |
| 8 | 4 | 297,133 |
| 16 | 5 | 98,172 |

Monotonic and correct: higher margin → more abstentions → more, smaller blocks.
The boundaries here need margin ≥ 4 before any abstains, i.e. they are
**well-supported merges**, not discordant over-merges.  At the default 500 kb
chunk size the hybrid single 498 kb block survives margins up to 128+,
confirming that stitch is genuine.

### Tests

`src/test_phase_block_stitch.cpp` gains two cases: a single agreeing overlap
read (score +1) merges at margin 0 and abstains at margin 1.  All unit tests
green; build has zero warnings from project code (one pre-existing abPOA
third-party `-Wunused-function` warning is unrelated).

### Open

- Still pending the diplinator switch-error gate to prove abstain margins trade
  N50 for accuracy as intended.  The margin gives a knob to *tune* that
  trade-off once truth-based evaluation is available.

---

## Hybrid Model — graph SNP low-complexity "gap" (REVERTED — was a bug)

An earlier change made `apply_hybrid_noise_filter` demote graph-only het SNPs in
SDUST low-complexity regions to `NonVariant`, on the premise that "this is what
the BAM pipeline does to density-noisy SNPs." **That premise was wrong, and the
change has been reverted.**

### Why it was wrong

Direct measurement on chr20 25M at the 10 demoted positions showed both source
pipelines *keep* these as real het calls:

- **BAM pipeline**: recalls them as `NOISY_CAND_HET` via noisy-region MSA (step 4)
  — AF 0.40–0.60, depth 46–78. The BAM pipeline does **not** discard them.
- **Standalone graph pipeline**: emits them as `CLEAN_HET_SNP` and phases them.

So the demotion dropped 10 genuine het SNPs that both pipelines call. The
`NOISY_CAND_HET` recall happens during phasing (after the bridge step), so these
sites are not in the BAM candidate table at graph-injection time and therefore
appear as graph-only — but they are still real variants, not noise.

### Resolution

`apply_hybrid_noise_filter` again leaves SNPs untouched (matching the standalone
graph pipeline's `apply_graph_noise_filter`); only indels are screened for
homopolymer / repeat / low-complexity noise.

---

## Hybrid Model — candidate-leak bug (graph-only non-calls not pruned)

While verifying the SNP revert, hybrid emitted **7141** candidates on chr20 25M
versus 1108 (BAM) and 509 (graph).

### Root cause

The BAM pipeline drops non-call candidates (`LowCoverage` / `NonVariant` /
`StrandBias`) via `prune_not_candidate_variants` at the end of
`collect_var_classify`, and folds `LowAlleleFraction` → `LowCoverage` in
classification pass 2 so those are pruned too. The hybrid pipeline appends
graph-only candidates *after* `collect_var_classify` has already run, and
`classify_graph_only_candidates` only runs pass-1 (`classify_variant_initial`).
So graph-only sites that failed the gate kept `LOW_COV` / `LOW_AF` / `NON_VAR`
categories and flowed straight to output — `merge_chunk_candidates` does not
filter by category. The standalone graph pipeline avoids this because its output
builder (`graph_chunks_to_candidate_table`) only emits clean/noisy categories.

The 7141 breakdown: 1234 real calls + 5875 `LOW_COV` + 32 `LOW_AF` + 11
`NON_VAR`.

### Fix

- `classify_graph_only_candidates` folds `LowAlleleFraction` → `LowCoverage`,
  matching the BAM pipeline's pass-2 rewrite.
- `process_chunk_hybrid` calls `prune_not_candidate_variants(chunk)` after
  k-means phasing. Pruning after phasing is safe: per-candidate read profiles
  (`read_var_profile`, `read_var_cr`) are consumed entirely during phasing and
  not used afterward; stitching, merge, and TSV/VCF/phased-BAM output use
  read-level state, not candidate indices. No post-phasing step (k-means or
  step-4 MSA) ever assigns a prunable category, so pruning after rather than
  before phasing removes exactly the same candidates the BAM ordering would.
- Exposed `prune_not_candidate_variants` (was file-static).

### chr20 25M result

- Hybrid: **7141 → 1234**. Categories are now only real calls
  (CLEAN_HOM / CLEAN_HET_SNP / NOISY_CAND_HET / REP_HET_INDEL / NOISY_CAND_HOM /
  CLEAN_HET_INDEL); no LOW_COV / LOW_AF / NON_VAR leak.
- BAM-only (1108) and graph-only (509) outputs byte-identical (MD5 verified).
- BAM↔hybrid site parity: of 1108 BAM sites, only 2 differed — both DEL
  encoding differences (`AA>.` vs `AA>A`). These were later traced to a real
  deletion key-collision bug and fixed (see "deletion key collision" section
  below); parity is now 1108/1108.
- Phased VCF / candidate VCF / phased BAM (HP+PS tags) all emit correctly.
- `make unit-tests` all green; added a `prune_not_candidate_variants` test and
  updated the noise-filter test to assert SNPs are kept.

---

## Hybrid Model — deletion key collision overwrote BAM calls

The 2 BAM↔hybrid site differences noted above (`AA>.` vs `AA>A`) were not a
cosmetic encoding choice — they were a real bug: the hybrid path silently
**overwrote a verified BAM deletion** with a graph deletion of a different
length.

### Root cause

`vcf_to_variant_key` (hybrid_inject.cpp) built deletion keys by stripping only a
**single** anchor base from the VCF REF/ALT, instead of the full shared prefix.
For a homopolymer site like `25412767 TAA→TA` (a 1 bp deletion) it produced
`pos=25412768, ref_len=2, alt="A"` — an internally inconsistent key claiming a
2-base ref span while keeping a residual `alt="A"`.

Two downstream consequences:

1. **No bridging.** `find_matching_candidate` requires `k.alt == target.alt`.
   BAM deletions always have `alt=""` (`variant_key_from_digar`), so a graph
   deletion key with `alt="A"` could never match — graph deletions were always
   added as separate candidates instead of bridging.
2. **Silent overwrite.** `exact_comp_var_site` **ignores `alt` for deletions**
   (it keys deletions on `pos` and `ref_len` only). The inflated `ref_len=2`
   made the graph 1 bp deletion compare *equal* to the genuine BAM **2 bp**
   deletion at the same locus. In `merge_chunk_candidates` the graph candidate
   (`lcd_make_variants_region_pass=true`) overwrote the BAM call, so only the
   graph allele survived.

This is why the VCF site `TAA→TA,T` (genuinely biallelic: a 1 bp and a 2 bp
deletion) collapsed to a single hybrid row that did not match BAM.

### Fix

Strip the **full** shared prefix in the deletion branch of
`vcf_to_variant_key`, matching both the BAM path (`variant_key_from_digar`) and
the standalone graph path (`graph_collect.cpp`): `pos = first deleted base`,
`ref_len = deleted span`, `alt = ""`. For a single-base anchor this reduces to
the previous behavior, so only homopolymer / repeat deletions change.

Hybrid-only change; the BAM and graph pipelines do not call this function.

### chr20 25M result

- BAM↔hybrid parity is now **complete**: all **1108/1108** BAM calls (and all
  **108/108** BAM deletions) are present in the hybrid output (was 1106 / 106).
- The biallelic homopolymer sites now emit **both** deletion alleles as distinct
  rows (e.g. `25412768 AA>.` 2 bp + `25412769 A>.` 1 bp), consistent across
  TSV / candidate VCF / phased VCF.
- BAM-only (1108) and graph-only (509) outputs remain byte-identical (MD5).
- Added a `vcf_to_variant_key` unit test asserting SNP / INS / DEL normalization
  and that 1 bp vs 2 bp homopolymer deletions encode distinct keys.

---

## Lessons / Pitfalls

- **pgbam sidecar causes over-stitching:** Running BAM pipeline with `--pgbam-file` merged all 49 initial phase sets into 2 chr20-spanning blocks with near-random accuracy (55%). Do not use the sidecar without understanding its signal quality.
- **Wrong reference = 0 reads:** `hg002v1.1.fasta` (maternal/paternal contigs) does not match the `CHM13#0#chr20` contig name in the aligned BAM. Must use `chm13v2.0.chr20.renamed.fa` as ref for `collect-bam-variation`.
- **Use chr20-only truth refs:** Full-genome maternal/paternal assemblies cause a small fraction of reads to leak to chr1/6/16 during truth alignment. Use `HG002_1.1_MATERNAL.chr20.fasta` / `HG002_1.1_PATERNAL.chr20.fasta` for clean results.
- **Per-chunk site loading:** Graph pipeline uses tabix-indexed per-chunk site queries (commit `5878aad`) — loading all 977K sites upfront was the prior bottleneck.
- **Non-numeric node IDs:** `build_compact_graph_site_index` now skips sites with non-numeric node IDs instead of aborting (commit `5c01954`).

---

## Overlap-Read Phase Propagation at Chunk Boundaries

### The problem

The genome is tiled into overlapping chunks for parallel phasing. A long read that spans a
chunk boundary exists as independent copies in both the upstream and downstream chunks. Each
chunk phases its copy using its own candidates. The BAM writer emits each read from its
**owning** (upstream) chunk only.

Example: chunks A (0–5 Mb) and B (4.5–10 Mb) overlap at 4.5–5 Mb. Read R starts at 4.8 Mb
and ends at 5.2 Mb. R appears in both chunks:

```
Chunk A (owns R):  candidates at 4.0, 4.3, 4.6 Mb
Chunk B:           candidates at 4.9, 5.1, 5.5, 5.8 Mb
                         ▲         ▲
                         R covers these in chunk B
```

In chunk A, R's variant profile covers the candidate at 4.6 Mb but no het sites are nearby —
R gets hap=0 (unphased). In chunk B, R covers candidates at 4.9 and 5.1 Mb, both het SNPs —
R gets hap=1, PS=4900000.

After stitching, `flip_chunk_hap` uses overlap reads (including R's phased copy in B) to
vote on whether to flip/merge chunks A and B. If the vote succeeds (flip_hap_score != 0),
chunk B's PS is rewritten to match chunk A's PS, and hap labels are flipped if needed. Now
R's copy in chunk B has HP/PS consistent with chunk A's phasing.

But the BAM writer emits R from chunk A, where it is still hap=0, PS=-1. R is written
unphased despite being successfully phased in chunk B. LongcallD has this same gap.

### The original (broken) fix

Commit `f094cea` added `propagate_overlap_read_phase_to_output_owner`, which copies HP/PS
from the downstream chunk to unphased upstream overlap reads. This recovered ~1,300 reads
but propagated **unconditionally** — including for unmerged pairs where `flip_hap_score == 0`.
For unmerged pairs the relative phase between chunks is unknown, so:

- The downstream hap=1 may correspond to the upstream hap=2 (wrong haplotype).
- The downstream PS doesn't exist in the upstream chunk (orphan phase set).

This caused pgphase to diverge from longcallD: +1,300 phased reads, +2 phase sets, with
some reads potentially assigned to the wrong haplotype.

### The fix

`stitch_chunk_haps` now tracks which adjacent pairs were successfully merged by
`flip_chunk_hap` via a `pair_stitched` vector. Propagation only fires for merged pairs,
where the downstream PS has been rewritten to match the upstream PS and hap labels have
been flipped if needed. Unmerged pairs are skipped.

### Evaluation (HG002 HiFi chr20)

| Metric | longcallD | pgphase (unconditional) | pgphase (merged-only) |
|---|---|---|---|
| Phased reads | 219,090 (80.5%) | 220,377 (81.0%) | 219,925 (80.9%) |
| Phase sets | 429 | 431 | 429 |
| Phase block N50 | 937,648 bp | 937,121 bp | 937,648 bp |
| Overall accuracy | 97.03% | 97.02% | 97.03% |
| Switch error rate | 0.61% | 0.62% | 0.61% |
| Phaseable accuracy | 98.15% | 98.14% | 98.15% |

The merged-only propagation recovers 835 reads (+0.4%) over longcallD with identical phase
set count, block structure, accuracy, and switch error rate.
