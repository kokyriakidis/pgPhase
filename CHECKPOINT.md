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

### Data provenance (full chr20 evaluation inputs)

The full-chr20 evaluation inputs are **not committed** (5.1 GB; gitignored under
`/test_data/chr20_eval/`). Re-download from the shared Google Drive folder:

```
https://drive.google.com/drive/folders/1DxDcuJ7uIugkQX_ZIbEMI5NdDg_lYT-5
```

Fetch with `gdown` (in a venv to satisfy PEP 668):

```bash
python3 -m venv /tmp/gdrive-venv
/tmp/gdrive-venv/bin/pip install -q gdown
/tmp/gdrive-venv/bin/gdown --folder \
  "https://drive.google.com/drive/folders/1DxDcuJ7uIugkQX_ZIbEMI5NdDg_lYT-5" \
  -O test_data/chr20_eval
```

| File (in `test_data/chr20_eval/`) | Role | Used by |
|------|------|---------|
| `chm13v2.0.chr20.renamed.fa` (+`.fai`) | CHM13 chr20 reference (`CHM13#0#chr20`) | BAM, graph, hybrid |
| `HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam` (+`.bai`) | HiFi reads aligned to CHM13 chr20 (272,016 reads) — graph `hs/hb/he` tags, **no** diplinator truth tags | BAM, hybrid |
| `chr20.sites.striped.vcf.gz` (+`.tbi`) | snarl sites | graph, hybrid |
| `HG002.chr20.annotated.coord.gaf.gz` (+`.tbi`) | graph alignments (GAF) | graph, hybrid |
| `HG002.chr20.fq.gz` | HiFi reads (FASTQ) | diplinator truth-BAM generation only |
| `HG002_1.1_MATERNAL.chr20.fasta` (+`.fai`) | truth maternal assembly | diplinator evaluation only |
| `HG002_1.1_PATERNAL.chr20.fasta` (+`.fai`) | truth paternal assembly | diplinator evaluation only |

Notes:
- The aligned BAM is "annotated" with graph `hs/hb/he` tags, **not** diplinator
  `HO:Z:`/`hq:i:` truth tags. The truth BAM that `evaluate_phase_accuracy.py`
  scores against must still be generated (align FASTQ to MAT/PAT assemblies →
  diplinator).
- The Drive folder ships **inputs only** — no pre-phased pipeline outputs. Both
  the BAM and hybrid phased BAMs must be produced by running the pipelines here.

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
- `Options::stitch_min_margin` (struct default `kDefaultStitchMinMargin = 0`).
- The hybrid command overrides the struct default to
  `kHybridDefaultStitchMinMargin = 10` in `collect_hybrid_variation`;
  `--stitch-min-margin INT` still overrides it per run.
- CLI flag `--stitch-min-margin INT` added to `collect-hybrid-variation`
  **only**.  The BAM and graph commands do not expose it and keep margin 0, so
  their behavior is unchanged.

### Hybrid default raised to 10 (matched-eval result)

Matched-eval on HG002 chr20 (identical evaluator, truth, and binary for BAM
and hybrid) showed that switch and flip rates are flat (~0.64% / ~0.93%) across
all configurations including pure BAM — they are intrinsic to the shared
k-means + stitcher core, not a hybrid regression.  The metrics that actually
separate the pipelines are Hamming and auN.  On those, raising the hybrid
stitch margin from 0 to 10 is a strict improvement:

| Config (af-indel 0.11) | Hamming | auN | switch | flip | perfect | covered |
|---|---|---|---|---|---|---|
| BAM baseline      | 3.09 | 9.83M  | 0.64 | 0.933 | 29.1% | 313.5M |
| hybrid margin 0   | 3.05 | 12.08M | 0.66 | 0.968 | 22.2% | 312.0M |
| hybrid margin 10  | 3.04 | 11.90M | 0.65 | 0.966 | 23.1% | 317.0M |

Margin 10 lowers Hamming, raises genome covered (+5 Mb) and perfect phase sets
(+0.9 pt) versus margin 0 at no switch/flip cost, so it is the hybrid default.

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

## Hybrid Model — graph het indels over-demoted (hybrid ≈ BAM)

### Symptom

Full-chromosome benchmarks showed the hybrid barely beating BAM (accuracy
96.86% vs 97.03%, phased reads 81.2% vs 80.9%) while the standalone graph
pipeline phased **93.4%** of reads. The hybrid was reproducing BAM instead of
harvesting the graph's phasing power, and phase-block auN regressed (11.8M →
8.9M).

### Diagnosis (chr20 25M slice)

Phasing in the hybrid runs only `collect_var_run_phasing`, whose first (and
only) k-means uses `kCandGermlineClean` — **clean het SNP + clean het indel +
clean hom**. So the count of *clean* het anchors drives phasing. Measured:

| anchors into k-means | BAM | GRAPH | HYBRID (old) |
|---|---|---|---|
| clean het SNP + clean het indel | 384 | **509** | 438 |

The graph's advantage is its **113 clean het indels**; the old hybrid kept only
**41**, demoting **29 → RepeatHetIndel** and reclassifying the rest. Two causes,
both from running BAM-conservative logic on graph-only candidates:

1. **classify_graph_only_candidates** called `classify_variant_initial` (the BAM
   classifier), which demotes homopolymer/repeat het indels to RepeatHetIndel
   **inline**, keyed on the normalized `VariantKey`. The standalone graph
   pipeline's `classify_graph_candidates` has no such inline demotion.
2. **apply_hybrid_noise_filter** reconstructed vcf_ref/vcf_alt from the
   normalized key, so `is_noisy_site` saw a different string/position than the
   graph pipeline's `apply_graph_noise_filter` (which screens on the original
   catalog ref/alt) — a different homopolymer verdict for the same site.

### Fix (hybrid-only)

Mirror the graph pipeline's two-step exactly:
`classify_graph_candidates` (depth/AF/het/hom, **no** inline indel demotion)
followed by `apply_graph_noise_filter` (screens on catalog strings).

- `classify_graph_only_candidates` now classifies graph-only candidates with the
  graph pipeline's depth/AF/het/hom logic (replicated inline; the shared
  `classify_variant_initial` is left untouched for the frozen BAM pipeline) and
  does **no** homopolymer demotion. LowAlleleFraction is still folded to
  LowCoverage to match BAM pass-2 pruning.
- `inject_graph_sites` now records each graph-only candidate's original VCF
  (pos, ref, alt) in a `GraphOnlyVcfAlleles` side-map (remapped through the
  post-injection sort permutation).
- `apply_hybrid_noise_filter` uses those original strings when available, so its
  homopolymer/repeat verdict matches the standalone graph pipeline; it falls
  back to key-reconstruction only when unavailable.

### Result (chr20 25M slice)

- Clean het anchors **438 → 448**; clean het indels **41 → 51**;
  RepeatHetIndel **54 → 44**.
- Phased het variants (phased VCF) **527 → 535**.
- The remaining gap to the graph's 509 clean anchors is **BAM-origin** sites
  whose linear pileup is genuinely noisy (`NOISY_CAND_HET`): of the graph
  clean-het sites the hybrid still marks noisy, **all 20 are bridged BAM-origin**
  and **0 are graph-only** — i.e. real BAM/graph disagreements, not a bug.
  Touching them would alter the frozen BAM classification, so they are out of
  scope.
- Invariants held: BAM (1108) and graph (509) outputs byte-identical (MD5);
  BAM call parity 1108/1108; `make unit-tests` green. Added a noise-filter test
  asserting original VCF strings override key-based reconstruction.

This slice is a single phase block (saturated N50), so the block-contiguity /
auN effect must be confirmed on a full chromosome with truth-based evaluation
(diplinator/whatshap unavailable in this environment).

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

---

## Three-Pipeline Comparison — HG002 HiFi chr20 (current)

Full chr20 evaluation using HG002 HiFi reads aligned to CHM13 chr20, evaluated
against HG002 T2T diploid assembly (diplinator). All numbers reflect the current
binary after all fixes above, including the indel noise-filter fixes (content-base
anchoring, nt4 insertion compare, non-minimal graph-allele trimming).

"Hybrid (old)" = BAM-base + graph-augment with the step-4 noisy-candidate k-means
re-orientation still on. "Hybrid (new)" = current default after disabling that
re-orientation (`skip_noisy_kmeans`, on by default for hybrid; `--keep-noisy-kmeans`
restores the old behaviour). See "Hybrid step-4 re-orientation" below.

| Metric | BAM | Graph | Hybrid (old) | **Hybrid (new)** |
|---|---|---|---|---|
| **Completeness** | | | | |
| Input reads | 272,016 | 231,382 | 272,016 | 272,016 |
| Phased reads | 219,925 (80.9%) | 203,729 (88.0%) | 220,725 (81.1%) | 212,484 (78.1%) |
| Phase sets | 429 | 368 | 415 | 348 |
| Perfect phase sets | 125 (29.1%) | 142 (39.1%) | 96 (23.5%) | **114 (33.4%)** |
| **Contiguity** | | | | |
| Phase block N50 | 937,648 bp | 917,581 bp | 945,071 bp | **1,041,015 bp** |
| Phase block auN | 9,834,145 bp | 989,988 bp | 11,876,987 bp | **13,677,722 bp** |
| Largest block | 53,165,843 bp | 2,080,875 bp | 59,316,535 bp | 58,616,492 bp |
| Median block span | 867,222 bp | 871,853 bp | 870,769 bp | 867,209 bp |
| Genome covered | 313,517,513 bp | 221,616,803 bp | 317,945,255 bp | 266,888,983 bp |
| **Accuracy** | | | | |
| Concordant reads | 213,127 | 202,076 | 213,609 | 210,724 |
| Discordant reads | 6,798 | 1,639 | 7,099 | 1,741 |
| Overall accuracy | 96.91% | 99.20% | 96.78% | **99.18%** |
| Hamming error rate | 3.09% | 0.80% | 3.22% | **0.82%** |
| Switch errors | 1,393 | 367 | 1,386 | **319** |
| Flip errors | 2,052 | 858 | 2,135 | 953 |
| Switchflip errors | 3,445 | 1,225 | 3,521 | **1,272** |
| Switch opportunities | 219,496 | 203,352 | 220,299 | 212,124 |
| Switch error rate | 0.63% | 0.18% | 0.63% | **0.15%** |
| Switchflip error rate | 1.57% | 0.60% | 1.60% | **0.60%** |

### Latest verified run (commit `c36026e`, full metrics + timings)

Re-run of all three pipelines after the snarl-duplicate dedup (`d3ac898`), the
CLI arg-parse crash fix (`49e4a49`), and the `--exclude-bed` revert (`c36026e`).
HG002 HiFi chr20 → CHM13, diplinator truth, `--min-reads 5`, 8 threads. Runtime
is wall-clock for the collect+phase run (excludes eval and BAM index/sort).

**Accuracy & errors**

| Metric | BAM | Graph | **Hybrid** |
|---|---|---|---|
| **Runtime** | 98 s | **53 s** | 326 s |
| Accuracy | 96.909% | 99.210% | **99.232%** |
| Hamming | 3.091% | 0.790% | **0.768%** |
| Switch errors | 1,393 | 364 | **291** |
| Flip errors | 2,052 | 856 | 906 |
| Switchflip | 3,445 | 1,220 | **1,197** |
| Switch rate | 0.635% | 0.179% | **0.137%** |
| Switchflip rate | 1.570% | 0.600% | **0.565%** |
| Switch opportunities | 219,496 | 203,357 | 211,900 |
| Concordant reads | 213,127 | 202,110 | 210,610 |
| Discordant reads | 6,798 | 1,610 | 1,630 |

**Yield & phase blocks**

| Metric | BAM | Graph | **Hybrid** |
|---|---|---|---|
| Total input reads | 272,016 | 231,382 | 272,016 |
| Phased reads | 219,925 | 203,734 | 212,259 |
| Fraction phased | 80.85% | **88.05%** | 78.03% |
| Reads evaluated | 219,925 | 203,720 | 212,240 |
| Phase sets (eval) | 429 | 363 | 340 |
| Perfect PS | 125 | **142** | 115 |
| Perfect PS % | 29.14% | **39.12%** | 33.82% |
| Candidates | 120,124 | 73,627 | 134,646 |
| Block N50 (bp) | 937,648 | 917,581 | **993,831** |
| Block auN (bp) | 9,834,145 | 989,984 | 7,153,008 |
| Block median span (bp) | 867,222 | 871,853 | 866,966 |
| Block max span (bp) | 53,165,843 | 2,080,875 | 39,388,202 |
| Genome covered (bp) | 313,517,513 | 221,615,836 | 246,734,240 |

**Phaseable vs unphaseable split** (confidence threshold 0.6)

| Metric | BAM | Graph | **Hybrid** |
|---|---|---|---|
| Phaseable PS / reads | 366 / 214,444 | 343 / 202,881 | 324 / 211,704 |
| Phaseable accuracy | 98.02% | **99.40%** | 99.34% |
| Unphaseable PS / reads | 63 / 5,481 | 20 / 839 | 16 / 536 |
| Unphaseable accuracy | 53.62% | 54.11% | 55.04% |

- **No regression.** Hybrid holds at 99.232% / 0.768% / switch 291 (matches the
  pre-/post-dedup verification). Graph holds the Fix-B no-pool gain: switch
  367 → 364, accuracy 99.195% → 99.210%.
- **Hybrid wins every error metric** (best accuracy, Hamming, switch, switch+flip
  rate) and the longest block N50 (993,831 bp), but is the slowest (326 s).
- **Graph is the accuracy-per-second winner** (99.210% in 53 s), highest fraction
  phased (88.05%), most perfect phase sets (142 / 39.12%), and best phaseable
  accuracy (99.40%). Its auN is far lower because it produces many uniform ~0.9 Mb
  blocks rather than a few huge spanning blocks.
- **BAM** has the largest auN/max-span (one 53 Mb block) and most genome covered,
  but the worst accuracy (3.091% Hamming) — the long blocks accumulate errors. It
  is the fallback when no graph/GAF inputs are available.
- **auN caveat:** BAM/hybrid's large auN comes from a single chromosome-spanning
  block; for switch-error quality the per-read accuracy metrics above are the
  meaningful comparison, not auN.

### Observations

- **The indel noise-filter fix transforms the graph pipeline.** Graph overall
  accuracy jumped from 92.70% to **99.20%** (Hamming 7.30% → **0.80%**) and switch
  error rate fell from 1.07% to **0.18%**. The fix demotes 17,509 homopolymer/STR
  het indels to `REP_HET_INDEL` so they no longer enter k-means as phasing anchors.
  Phased reads drop (216k → 204k) because those noisy indels previously phased
  many reads incorrectly — the pipeline now trades a small completeness loss for a
  large accuracy gain. Perfect phase sets rose from 22.0% to **39.1%**.
- **Hybrid (new) reaches graph-pipeline accuracy at far higher contiguity.** After
  disabling the step-4 noisy-candidate k-means re-orientation, hybrid jumps from
  96.78% to **99.18%** (Hamming 3.22% → **0.82%**, switch rate 0.63% → **0.15%**).
  It now matches standalone graph accuracy (99.20% / 0.18%) while keeping BAM-class
  contiguity (auN 13.7M vs graph's 0.99M — 14× longer blocks) and phasing 8.8k more
  reads than graph (212k vs 204k). It even has *fewer* switch errors than graph
  (319 vs 367). This is the single biggest hybrid improvement to date — see
  "Hybrid step-4 re-orientation" below for the mechanism.

### Hybrid step-4 re-orientation: the root cause of hybrid's accuracy gap

The post-fix marginal analysis (below) localised 96% of hybrid errors to the
**BAM core**, not the graph layer. The cause was found by A/B toggle: the BAM
pipeline's **step-4 noisy-region recall** (`collect_noisy_vars_step4`) runs an MSA
to recall variants inside noisy regions, then **re-runs k-means over
`kCandGermlineVarCate`** (which includes the recalled `NOISY_CAND_HET` candidates)
to re-assign read haplotypes. The standalone graph pipeline never does this — it
runs a single k-means over `kCandGermlineClean` only.

On hybrid, this second k-means was **actively harmful**:

- It phased **8,243 extra reads at 65% error** (5,358 wrong, 2,885 right).
- It poisoned the shared core: BAM-shared Hamming went **0.71% → 3.11%** when on.
- It also *fragmented* blocks (auN 13.7M with it off vs 11.9M on; perfect PS 114 vs 96).

**The fix (`skip_noisy_kmeans`, default on for hybrid):** keep the MSA variant
recall (noisy variants still appear in the output VCF, `NOISY_CAND_HET=22,240`) but
skip *only* the noisy-candidate k-means re-run. Confirmed by two equivalent toggles
giving identical results: skipping all of step 4 (`exp_skip_step4`) and skipping
just the re-run both yield 99.181% / 0.819% — proving the **re-orientation, not the
recall, was the entire harm**. `--keep-noisy-kmeans` restores the old behaviour.
The BAM pipeline is unchanged (the re-run still runs there; only hybrid defaults to
skipping it, because hybrid's graph-anchored first k-means already orients reads
well and the noisy re-run only adds error).

**Knob sweeps on the clean base were saturated** (no further gain):
`--graph-indel-af-margin` {0.11→0.50} *worsened* accuracy (admitting off-center
graph indels adds noise); `--stitch-min-margin` {4,10,15,20} and `--stitch-rule`
{1,3} moved nothing on accuracy (the both-strands rule dominates merges, not the
margin). `--stitch-rule 3 --stitch-min-margin 15` adds +6 perfect PS (120 vs 114)
and +5M genome covered at identical accuracy — a valid optional contiguity knob,
not promoted to default.

### Result is deterministic, not seed-dependent

The k-means has no RNG (greedy seeded init); repeated runs of the same config give
bit-identical metrics (319 switches every time). The "convergence variance"
referenced in earlier notes was about a separate experimental optimizer, not the
production k-means.

> **Remaining caveat (graph hybrid-allele trimming):** the non-minimal
> multiallelic-allele trimming added to `apply_graph_noise_filter` is still NOT in
> `apply_hybrid_noise_filter`. Low priority now — hybrid (new) error is dominated by
> the graph-influenced hard tail (1,660 reads at 14.6%), and the BAM-shared core is
> already at 0.71%. Revisit only if graph-only repeat indels are shown to drive the
> remaining hybrid errors.

### Post-fix error decomposition: is graph's accuracy transferable to hybrid? (NO)

`split_marginal.py` on the post-fix hybrid per-read results, splitting discordant
reads by origin (graph-influenced = phased by hybrid but NOT by standalone BAM;
bam-shared = phased by both):

| Read subset | Reads | Discordant | Hamming |
|---|---:|---:|---:|
| graph-influenced (marginal) | 1,704 | 284 | **16.67%** |
| bam-shared | 219,004 | 6,815 | 3.11% |
| ALL hybrid-phased | 220,708 | 7,099 | 3.22% |

> **SUPERSEDED conclusion — read with the step-4 finding above.** This section
> concluded "96% of hybrid errors are in the BAM core, which is unfixable, so you
> cannot get graph's low Hamming at BAM coverage." The *diagnosis* (96% BAM-core)
> was correct and led directly to the fix; the *prognosis* (unfixable) was WRONG.
> The BAM-core errors were caused by the step-4 noisy k-means re-orientation, and
> disabling it took hybrid to 99.18% / 0.82% at 212k reads — graph accuracy at far
> higher contiguity. The "denominator effect" framing below is therefore only
> partly true: graph's accuracy is real (see shared-read analysis) AND hybrid can
> now match it. Kept for the diagnostic trail.

Two conclusions for the "lower Hamming/switch" goal (the second is now superseded):

- **96% of hybrid errors are in the BAM core** (6,815 / 7,099). ✅ Correct, and the
  key clue: it pointed the search at a BAM-core operation (step-4 re-orientation),
  not the graph layer. The BAM-shared Hamming (3.11%) matched standalone BAM (3.09%)
  *with the re-orientation on*; with it off, BAM-shared drops to **0.71%**.
- **~~Graph's 0.80% is a pure denominator effect; you cannot get it at BAM
  coverage.~~** SUPERSEDED — hybrid (new) gets 0.82% at 212k reads. Graph's
  marginal reads do error at 16.67% in hybrid (the hard tail is real), but the bulk
  of the gap was the fixable step-4 re-orientation, not an intrinsic ceiling.

All rows in the table above were regenerated from the same current HEAD binary
(BAM/graph/hybrid all run with shipped defaults; hybrid = af0.11 + stitch margin
10 + both-strands rule), scored against the same diplinator truth BAM.

### Shared-read accuracy: graph IS better on the same reads (corrects earlier claim)

Comparing per-read concordance on the **201,202 reads both graph and BAM phase**
(`eval_graph/per_read.tsv.gz` ∩ `eval_bam/per_read.tsv.gz`, same truth BAM):

| On 201,202 shared reads | Graph | BAM |
|---|---:|---:|
| Discordant | 1,247 (**0.62%**) | 2,896 (**1.44%**) |

Per-read disagreement on those shared reads:

| Outcome | Reads |
|---|---:|
| graph RIGHT, bam WRONG | **2,082** |
| bam RIGHT, graph WRONG | 433 |
| both WRONG | 814 |

**Graph wins disagreements 4.8 : 1.** This means graph's accuracy is *not* purely a
denominator effect — on reads both attempt, graph is ~2.3× more concordant. The
earlier "pure denominator effect" framing (marginal-split section) was too strong.

**Confound (unresolved):** graph blocks are far shorter (auN 0.99M vs 9.83M; 368
vs 429 PS). A read in a shorter block spans fewer adjacent-het transitions, so it
is *mechanically* less likely to be scored discordant. Graph's 0.62% is therefore
**part genuinely-better-phasing, part shorter-blocks**; the two are not separated
by this measurement. Any strategy that lengthens graph blocks (e.g. by adding BAM
sites) would partly surrender the short-block component.

### Incremental (non-shared) reads — the hard tail

| Reads only one pipeline phases | Count | Discordant |
|---|---:|---:|
| graph-only | 2,513 | 15.6% |
| **bam-only** | **18,723** | **20.8%** |

The 18,723 reads BAM phases but graph does not error at 20.8% — graph lacks sites
there precisely in the **hard regions** (segdups), the same windows (37–38 Mb)
where prior analysis found "pulling hard reads into blocks injects switches."

### Proposed strategy: invert hybrid (graph-base + BAM-fill)

Idea (user, this session): since graph now performs so well, flip the hybrid base
from BAM to graph, and bring **validated BAM sites into the graph where the graph
has no sites**, instead of the current BAM-base + graph-augment design.

**Why it has merit:** graph is genuinely more accurate on shared reads (above), so
keeping graph as the core and only filling gaps could retain graph's clean phasing
while recovering coverage.

**Optimistic ceiling:** if the graph core keeps ~0.80% and the ~18,723 bam-only
reads are added at BAM accuracy: ≈5,541 errors / 222,438 reads ≈ **2.5% Hamming**
— would beat BAM (3.09%) and current hybrid (3.22%) at *higher* coverage.

**Risks that could collapse the ceiling:**
1. Adding BAM sites lengthens graph blocks → reawakens the short-block confound
   (some of graph's 0.62% edge was mechanical) and the segdup-switch penalty
   established in the rejected-lever ledger.
2. The bam-only fill reads error at 20.8% — the gain is bounded by how many are
   *recoverable* and whether merging them keeps graph's core clean.

**Feasibility probes to run BEFORE building (offline, no C++ — matches this
project's measure-first track record):**
1. **GAF membership of the 18,723 bam-only reads.** If present in the GAF but
   unphased for lack of sites → recoverable. If absent from the GAF → unreachable
   regardless. (Established globally: GAF 271,930 ⊂ BAM 272,016; confirm for *this*
   subset.)
2. **Block-length-fair accuracy.** Re-score graph vs BAM on shared reads within
   comparable block spans, to isolate real phasing quality from the short-block
   artifact and get a realistic (not optimistic) ceiling.

**Status:** not yet started — awaiting decision on probe-first (recommended) vs
build-first. This is the live thread to resume.

### Probe results: graph-core + BAM-gap-fill (the narrower strategy)

Refined goal (user): keep graph's phased core untouched; for reads graph leaves
**unphased**, phase them using BAM-validated sites the graph lacks. Purely additive
gap-fill — no merging/re-orienting of graph's clean blocks (avoids the short-block
confound). Both probes run offline on the post-fix per-read results:

**Probe 1 — reachability (PASS).** Gap-fill candidate set = reads BAM phases but
graph does not = **18,723 reads**. Of these, **18,668 (99.7%) are present in the
GAF** — graph sees the reads, it just has no phasing sites there. Only 55 are
unreachable. So the strategy is mechanically feasible.

**Probe 2 — accuracy ceiling (the catch).** The fill reads error at **20.84%**
under BAM, because graph lacks sites there *precisely because they are the hard
regions* (segdups / low-complexity; concentrated in 45–46, 21–22, 44–45, 13–14,
36–37 Mb windows). Projected combined:

| Method | Reads | Hamming |
|---|---:|---:|
| Graph (core only) | 203,715 | 0.80% |
| **Graph + BAM-fill (optimistic)** | **222,438** | **2.49%** |
| BAM | 219,925 | 3.09% |
| Hybrid (current) | 220,708 | 3.22% |

**Conclusion:** the strategy *works* and gives the best accuracy/coverage combo of
any method — **~2.49% Hamming at 222k reads** (beats BAM 3.09% and current hybrid
3.22%, at higher coverage). BUT it does **not** preserve graph's 0.80%: the gap
reads are intrinsically the hardest, so adding them at ~21% error drags the
combined rate up ~3×. Graph's advantage is *diluted, not inherited*. And 2.49% is
the **optimistic ceiling** — it assumes fill reads phase at full BAM accuracy and
that seam-stitching graph blocks via BAM sites injects no new switches (the
dominant real-world failure in the rejected-lever ledger). Realistic outcome sits
between 2.49% and current hybrid's 3.22%.

**Recommendation:** worth building IF the target metric rewards
accuracy-at-high-coverage (2.49% < 3.09% BAM is a real gain). NOT worth it if the
goal is graph's headline 0.80% — that is unreachable at full coverage by
construction, because the missing reads are the hard ones. Decision pending.

### Built and measured: hybrid-core + gap-fill (additive, no re-stitch)

The gap-fill strategy above was prototyped and evaluated on real chr20 data
(binary `9dcc509`, HG002 HiFi → CHM13, diplinator truth). Implementation
(`run_all3_fixed/gapfill.py`): take the hybrid phased BAM as the untouched core;
for the 9,169 reads BAM phases but hybrid leaves unphased, copy BAM's `HP` onto
the hybrid record with `PS += 1e9` (a disjoint phase-set namespace). Purely
additive — no hybrid block is renumbered, merged, or re-oriented.

**Provenance (important):** the bam / graph / hybrid pipeline outputs are all
phasings of the **same vg-giraffe alignment**, not different aligners. The BAM
is the linear surjection of the graph alignment (`@PG vg` → `pgbam annotate`);
the GAF is its graph-space form. So a "fill read" is not from another aligner —
it is a read that a different pgphase *pipeline* managed to phase. Gap-fill is a
cross-pipeline HP/PS label transfer by read name, not a re-alignment.

The hybrid phased BAM is the correct substrate because it carries all 272,016
reads as records. The **fill** BAM is read only for its HP/PS labels, so it may
be unaligned: the graph pipeline's phased BAM has no `@SQ` header, but
`samtools view` still yields its tags, and the labels land on the core's
already-aligned records. (Earlier notes claimed graph fill "needs realignment"
— that was wrong; it conflated using the graph BAM as the *core* with reading
its labels.)

**Which reads hybrid drops vs BAM (measured on this run):** 9,169 reads. They
error at **36.35%** even under the pipeline that phases them (BAM), are all
mapq 60 (not a mapping artifact), and concentrate in segdup/low-complexity
windows (44–46 Mb ≈ 3,805 reads at 38–53% error; 13–14, 35–36, 43 Mb). Graph
phases only 625 of them, at 64% accuracy. These are the deliberately-skipped
step-4 noisy reads — the `skip_noisy_kmeans` default drops them on purpose.

**Result:**

| Metric | BAM | Graph | Hybrid | **Hybrid+gapfill** |
|---|---:|---:|---:|---:|
| Phased reads | 219,925 | 203,734 | 212,259 | **221,428** |
| Overall accuracy | 96.91% | 99.21% | **99.23%** | 97.85% |
| Hamming | 3.091% | 0.790% | **0.768%** | 2.149% |
| Switch errors | 1,395 | 364 | **291** | 1,171 |
| Switch rate | 0.636% | 0.179% | **0.137%** | 0.530% |
| Flip errors | 2,050 | 856 | **906** | 1,958 |
| Perfect PS % | 29.14% | **39.12%** | 33.82% | 28.79% |
| Phase sets | 429 | 368 | 347 | 499 |
| N50 | 937,648 | 917,581 | **993,831** | 927,359 |
| auN | 9.83M | 0.99M | 7.15M | **11.25M** |
| Genome covered | 313.5M | 221.6M | 246.7M | **340.2M** |

**Findings:**
- **Hybrid+gapfill Pareto-beats pure BAM on every axis:** more reads phased
  (+1,503), lower Hamming (2.15% vs 3.09%, ~30% relative), fewer switch errors
  (1,171 vs 1,395), higher auN (11.25M vs 9.83M), more genome covered (340M vs
  314M). If a workflow ships/compares against BAM for high coverage, gap-fill is
  a strict upgrade.
- **The disjoint-PS design works:** switch errors rose only to 1,171 (below BAM's
  1,395), confirming the fill reads create their own small phase sets without
  contaminating the hybrid core. Measured Hamming (2.149%) **beat** the offline
  projection (2.25%) for this reason.
- **It does NOT preserve hybrid's 0.77% Hamming.** Adding the hard segdup reads
  (36% error) drags Hamming 0.77% → 2.15%. Unavoidable: the oracle best-of-3
  per-read ceiling is 2.16% Hamming at 222.8k reads (4,804 reads are wrong in all
  three pipelines), so 2.15% is essentially at the achievable floor for that
  coverage.

**Two-source gap-fill (hybrid core + BAM fill + graph fill) — measured but
deliberately NOT shipped:** chaining a second `gapfill.py` pass with
`ps_offset=2e9` adds the **1,374 reads only the graph pipeline phased** (neither
hybrid nor BAM). All 1,374 are already present as records in the hybrid BAM — no
realignment — and phase at 81% in isolation. **Decision (final): gap-fill will
only ever recover reads present in the projected BAM** (the linear surjection of
the graph alignment that the hybrid pipeline already loads). The graph-only reads
are out of scope by design, so this row is kept for the record only — it is not a
candidate for promotion to the binary.

| Metric | BAM | Hybrid | +gapfill (1-src) | **+gapfill2 (2-src)** |
|---|---:|---:|---:|---:|
| Phased reads | 219,925 | 212,259 | 221,428 | **222,802** |
| Overall accuracy | 96.91% | **99.23%** | 97.85% | 97.77% |
| Hamming | 3.091% | **0.768%** | 2.149% | 2.228% |
| Switch errors | 1,395 | **291** | 1,171 | 1,244 |
| auN | 9.83M | 7.15M | **11.25M** | 10.62M |

222,802 is the full union of all three pipelines (+2,877 over BAM, the maximum
phaseable at this coverage). Adding the 1,374 graph-pipeline reads costs +73
switch errors (their ~19% error on 1,374 reads ≈ +258 discordant), so 2-source
maximizes coverage while 1-source is marginally cleaner. Both Pareto-beat BAM.

**Decision matrix (all measured):**

| Goal | Method | Hamming | Reads |
|---|---|---:|---:|
| Lowest error | Hybrid (default) | 0.77% | 212.3k |
| Best accuracy *and* coverage vs BAM | **Hybrid+gapfill** | 2.15% | 221.4k |
| Max coverage (full union, offline only — not shipped) | Hybrid+gapfill2 | 2.23% | 222.8k |
| Raw BAM-class behavior | BAM / `--keep-noisy-kmeans` | 3.09% | 219.9k |

**Verification of `--keep-noisy-kmeans` (the in-pipeline alternative):** ran for
real — recovers reads to 220,633 but regresses to 96.79% / 3.21% Hamming /
1,371 switch (essentially BAM), confirming step-4 phases its extra reads at ~65%
accuracy. Additive gap-fill is strictly better than `--keep-noisy-kmeans` (2.15%
vs 3.21% Hamming at comparable coverage) because it does not let the noisy reads
re-orient the clean core.

**Status:** shipped natively as `collect-hybrid-variation --gap-fill`. The
1-source variant (hybrid core + BAM-pipeline fill → 221,428 reads / 2.15%
Hamming) is built into the binary: when `--gap-fill` is set alongside the
hybrid `skip_noisy_kmeans` default, `gap_fill_unphased_reads`
(`src/collect_var.cpp`) re-runs the `kCandGermlineVarCate` k-means into a
scratch buffer and adopts its haplotype only for reads the clean core left at
`hap==0`, writing them to `chunk.gap_haps`/`gap_phase_sets` with
`PS += kGapFillPsOffset` (1e9, disjoint namespace). The core read-phase and
per-candidate consensus state are snapshotted and byte-for-byte restored, so the
clean core is never renumbered, merged, or re-oriented — the in-binary form of
the validated `scripts/gapfill.py` 1-source pass. `collect_bam_output.cpp` emits
the gap-fill HP/PS for any read the core left unphased. Off by default.

The 2-source variant (adding the 1,374 graph-only reads → 222,802 / 2.23%) is a
**deliberate non-goal** and will not be promoted to the binary: gap-fill recovers
only reads present in the projected BAM, and the graph-only reads are out of scope
by design. It exists solely as a post-process for offline analysis via
`scripts/bench_hybrid.sh --gapfill 2`. The prototype (`scripts/gapfill.py`,
`scripts/bench_hybrid.sh --gapfill [N]`) is retained for reproducing both BAMs.

---

## Graph Het-Indel AF Gate (resolves hybrid accuracy/contiguity regression)

### Problem

Commit `af6b089` ("Stop over-demoting graph het indels") promoted graph-only het
indels to `CleanHetIndel` so they enter k-means as phasing anchors. On a full
chr20 truth evaluation this **raised** hybrid Hamming error to 4.52% (vs BAM
3.09%) while raising auN to 13.6M. The commit message itself flagged that the
accuracy recovery was never confirmed on a full chromosome — it was not realized.

### Investigation

Truth-based read-level analysis (HG002 diplinator tags) localized the damage:

- **STABLE phase-block interiors are untouched** (+27 discordant on 140,861
  shared reads). The joint k-means does *not* corrupt confident BAM regions.
- **All regression lives at block-merge seams.** ~2,701 reads are *wholesale
  orientation flips*: a BAM-clean sub-block inverted en masse after a graph
  anchor mis-oriented the merge.
- A `--stitch-min-margin` sweep (0→8) barely moved Hamming (4.52→4.37%): the
  wrong merges are **confident, not thin-margin**. The defect is upstream anchor
  selection, not the stitch decision.
- Reverting `af6b089` dropped Hamming to 3.27% but collapsed auN to 8.93M (below
  BAM). The graph het indels are a **mixture**: near-0.5-AF indels are reliable
  bridges (the auN gain); off-center-AF indels (mis-genotyped / repeat) are the
  wrong-orientation anchors. All-or-nothing loses on one axis.

### Fix

`classify_graph_only_candidates` (hybrid path only) now admits a graph het indel
as a `CleanHetIndel` anchor only when its allele fraction sits within
`graph_indel_af_margin` of 0.5; off-center indels are demoted to `LowCoverage`
(kept out of k-means, pruned from output). Two tunable gates added, defaults
reproduce nothing-changed for SNPs:

- `--graph-indel-af-margin` (default `kDefaultGraphIndelAfMargin = 0.11`)
- `--graph-indel-min-alt` (default 0 — swept, no additional effect over the
  existing `min_alt_depth` floor)

The BAM pipeline (`classify_variant_initial`) and graph-only pipeline
(`classify_graph_candidates`) are not touched; BAM candidate output verified
byte-for-byte identical.

### Result (chr20, full-chromosome truth eval)

AF-margin sweep, best at **0.11** (keep graph het indels with AF ∈ [0.39, 0.61]):

| Config | Hamming | auN | Largest block | Phased reads |
|---|---|---|---|---|
| BAM baseline | 3.091% | 9.83M | 53.2M | 219,925 |
| Hybrid af6b089 (margin 0.30) | 4.520% | 13.57M | 62.7M | 223,065 |
| Hybrid revert af6b089 | 3.265% | 8.93M | 50.2M | 220,741 |
| **Hybrid + af-gate 0.11** | **3.053%** | **12.09M** | **59.3M** | **220,952** |

The 0.11 gate **Pareto-beats BAM on every axis**: lower Hamming (−0.038pp,
~830 reads), +23% auN, +11% largest block, +1,027 phased reads. Sweep was
non-monotone with a clear optimum at 0.11 (0.12 already crosses back above BAM
accuracy), confirming a genuine edge rather than a tuning endpoint.

---

## Chunk-Stitch Rule: both-strands-bridged (hybrid default)

### Problem

The genome is phased in 500 kb chunks, then adjacent chunks are stitched by an
overlap-read vote in `flip_chunk_hap` (`collect_phase.cpp`). The original rule
merges when the net flip-vote magnitude exceeds `stitch_min_margin`. The net
margin conflates two very different seams: "1 clean read, 0 conflicts" and
"11 reads splitting 6-vs-5" both yield margin 1. Raising the margin to abstain
on the coin-flip case also abstains on the clean low-coverage case, leaving
blocks unmerged that could have been joined safely.

### Rules tested

`--stitch-rule` selects the decision (named constants in `phasing_types.hpp`):

- **0 net-margin** — original: merge when `|flip − noflip| > margin`.
- **1 both-strands-bridged (Reading A)** — merge only when the winning
  orientation has ≥1 read on *each* of its two haplotype links (no-flip needs
  `pre1→cur1` AND `pre2→cur2`; flip needs `pre1→cur2` AND `pre2→cur1`). Guards
  against merging on evidence from a single haplotype.
- **2 literal (Reading B)** — merge when both orientations have ≥1 read
  (stitches contested seams; diagnostic only).
- **3 both-strands + margin (Reading C)** — Reading A AND net vote > margin.

### Result (chr20, af0.11, full-chromosome truth eval)

| rule | Hamming hq0 | Hamming hq10 | auN | phased | blocks/chrom |
|---|---|---|---|---|---|
| 0 net-margin (control) | 3.0401 | 1.8550 | 11.91M | 220,734 | 411 |
| **1 both-strands (default)** | 3.0415 | 1.8541 | **12.00M** | **220,897** | **404** |
| 2 literal | 2.6530 | 1.4490 | 10.66M | 220,340 | 488 |
| 3 both-strands + m2 | 3.0385 | 1.8539 | 11.97M | 220,853 | 407 |

**Rule 1 improves contiguity at no accuracy cost**: auN +0.8%, +163 reads
phased, blocks/chrom 411 → 404, Hamming flat (hq0 +0.001, hq10 −0.001),
switch/flip unchanged. The control (rule 0) reproduced the baseline exactly,
confirming the refactor is behavior-preserving.

**Rule 2's apparent Hamming win is an artifact**, not an improvement: it only
merges contested seams and refuses clean unanimous ones, so blocks/chrom jumps
to 488 (genome fragmented into more, smaller blocks). Hamming is measured
within-block, so shorter blocks mechanically score lower while auN drops to
10.66M and 394 reads are lost. Diagnostic only.

### Decision

Rule 1 is the **hybrid default** (`kHybridDefaultStitchRule`), set in
`hybrid_collect.cpp` alongside the stitch-margin default. BAM/graph pipelines
keep rule 0. `--stitch-rule` overrides. This is the only lever found that adds
contiguity for free — every other lever (below) traded accuracy or auN away.

### Rule 4 confidence-weighted stitch (rejected)

The seam-merge regression traces to a minority of *confident-wrong* merges —
boundaries where the overlap vote is lopsided the wrong way (the margin sweep
0→8 barely moved Hamming, so the bad votes are not thin-margin). Hypothesis: a
cluster of reads all phased off one weak anchor casts many votes that are
really one piece of evidence; weighting each overlap read's vote by its k-means
assignment confidence `max(0,(agree−conflict)/(agree+conflict))` (from the
persisted clean-het SNP tallies, using `min(conf_pre,conf_cur)` per read) would
shrink such a seam's effective strand weight below threshold and abstain. Added
as `--stitch-rule 4` + `--stitch-min-strand-weight` (both-strands structure on
weighted sums).

Full chr20 sweep (diplinator truth eval; rule-1 control byte-identical to HEAD):

| config | Hamming% | auN | N50 | phased | perfect PS |
|---|---|---|---|---|---|
| rule 1 (shipped) | 3.0415 | 12.00M | 951,564 | 220,897 | 91 |
| rule 4 w≤1.0 | 3.0415 | 12.00M | 951,564 | 220,897 | 91 |
| rule 4 w=1.5–2.0 | 3.0416 | 12.00M | 951,564 | 220,861 | 92 |
| rule 4 w=3.0 | 3.0396 | 11.93M | 947,261 | 220,776 | 94 |
| rule 4 w=5.0 | 3.0405 | 11.87M | 947,261 | 220,678 | 96 |

**No-op at low weight, net loss at high weight.** At w≤1.0 the output is
byte-identical to rule 1: the bridging reads are *confidently* phased
(conf≈1.0), so weighting by confidence equals counting — confirming the
confident-wrong merges are not driven by weakly-assigned reads. Only at w≥3.0
does the gate bite, and then it shows the same accuracy/contiguity trap as every
other lever: a tiny Hamming gain (−0.0019) for real contiguity loss (auN −63 to
−127 kb, N50 951k→947k, −121 to −219 reads). Perfect-PS rises (91→96) but
overall it is not a win. Reverted from code; kept here so it is not re-tried.
The standing conclusion holds: the residual seam errors are confident and
upstream of the stitch decision — no reweighting or thresholding of the
overlap vote separates them from correct merges.

### BAM-CIGAR-only AF gate on graph indel anchors (rejected)

Since stitch-side fixes all failed, this attacked the genotype at its source.
The standard graph het-indel AF window uses **mixed BAM+GAF counts**: GAF
graph-read alt support can make an indel look het that BAM CIGAR evidence alone
shows as homozygous or thin-alt. The gate (`--graph-indel-bam-af-gate`)
snapshots BAM-CIGAR-only allele counts *before* `inject_graph_reads` adds GAF
observations, and admits a het-indel anchor only if its BAM-only AF is also
within the margin of 0.5 (for candidates with ≥ min-depth BAM reads). It reuses
the trusted BAM-CIGAR genotyping the BAM pipeline already does, upstream of
k-means.

Full chr20 sweep (diplinator truth eval; gate-off control byte-identical to HEAD):

| config | Hamming% | auN | N50 | phased | switchflip% | perfPS | anchors |
|---|---|---|---|---|---|---|---|
| off (shipped) | **3.0415** | 12.00M | 951,564 | 220,897 | 1.6206 | 91 | 134,616 |
| gate depth=3 | 3.2110 | 11.93M | 947,575 | 220,765 | 1.5790 | 99 | 128,830 |
| gate depth=5 | 3.2110 | 11.93M | 947,575 | 220,765 | 1.5790 | 99 | 128,910 |
| gate depth=8 | 3.2110 | 11.93M | 947,575 | 220,765 | 1.5790 | 99 | 128,987 |

**Net loss.** The gate is *not* timid — it removes ~5,700 graph indel anchors
(many graph indels genuinely have BAM-CIGAR AF inconsistent with their mixed
het call) and the depth knob barely changes the eval. Switchflip drops (1.62 →
1.58) and perfect-PS rises sharply (91 → 99), so it *does* remove
mis-genotyped, seam-mis-orienting anchors. But read-level **Hamming rises 3.04
→ 3.21%** and auN drops 62 kb: the removed anchors were, on net, contributing
more correct read assignments than wrong orientations. Even the trusted
BAM-CIGAR caller cannot separate the good from the bad — the indels it
disagrees with are still net-positive for phasing.

**Fourth confirmation of the same result.** Margin sweeps, the linkage gate,
confidence weighting, and now source-level BAM re-genotyping all fail
identically: any operation that *removes or down-weights* the suspect anchors
trades a real read-assignment/contiguity loss for a switchflip/perfect-PS gain
that does not net out. The mis-oriented seams and the contiguity-providing
bridges are the **same population** of off-center/GAF-supported indels; no
1-D signal (AF, LD, confidence, BAM-AF) separates them. Reverted from code;
kept here so it is not re-tried. If there is a win left it requires *correcting*
an anchor's orientation rather than *dropping* it — a fundamentally harder change
than any gate.

---

## Rejected graph-native levers (kept for the record)

Two graph-specific anchor-quality gates were prototyped and evaluated against
the shipped config (af0.11 + margin10), then **dropped** — neither helped on
HiFi without a worse cost elsewhere. Documented here so they are not re-tried.

### Graph read-divergence gate (`--graph-max-dv`)

The GAF carries a per-alignment divergence tag `dv:f` (the graph analog of the
BAM pipeline's XID/noisy-read filter, which the graph path lacks). Idea: drop a
read's allele votes when `dv` is high, before allele-fraction aggregation, so
noisy reads can't skew a site's AF off-center and trip the indel gate.

Result: at dv=0.005 it gave a marginal hq0 Hamming improvement (3.040 → 3.029)
but flat hq10 (1.855 → 1.854) and slightly lower auN/reads. The effect is tiny
because **HiFi reads are too accurate** — only 1.4% have dv>0.01, 0.26% >0.05
(median 0.0002, p99 0.0167). The filter is correct (extreme dv collapsed
candidate counts; dv=1.0 reproduced the baseline byte-for-byte) but there is
almost no noise to remove on HiFi. **Expected to matter on ONT**, where the
divergence tail is large. Not shipped.

### Graph SNP AF-window gate (`--graph-snp-af-margin`)

Mirrors the indel AF gate for graph het SNPs (default 0.5 = admit all). A 2×2
matched eval (af0.11, margin10):

| config | Hamming hq0 | Hamming hq10 | auN |
|---|---|---|---|
| neither (shipped) | 3.040 | 1.855 | 11.91M |
| dv-only (0.005) | 3.029 | 1.854 | 11.86M |
| sam-only (0.10) | 2.986 | 1.800 | 9.59M |
| both | 2.975 | 1.799 | 9.57M |

sam=0.10 drives the Hamming win but **collapses auN to 9.59M, below the BAM
baseline (9.83M)** — the same accuracy/contiguity trap as a too-wide indel
gate. Filtering off-center SNP anchors removes the long-range bridges that give
the graph its contiguity edge. Not shipped; default-off (admit all) remains
Pareto-best for SNPs.

### Graph anchor BAM-linkage gate (`--graph-anchor-min-concordance`)

The block-merge regression was localized to graph-only het-indel anchors that
mis-orient a seam (see "Graph Het-Indel AF Gate"). The AF window (0.11) is a
1-D proxy for "is this anchor's genotype reliable." This lever tried to replace
the proxy with a **direct** signal: a real graph het indel co-segregates with
neighboring confident BAM het SNPs on the reads spanning both (same physical
molecules), so per-anchor linkage concordance against BAM SNPs should separate
reliable bridges from mis-genotyped anchors better than AF alone. The gate runs
before k-means (so it uses `read_var_profile` linkage, not post-k-means
haplotype labels), computing each indel's best phase-agnostic 2×2 concordance
(`max(n00+n11, n01+n10)/total`) over BAM het SNP partners sharing ≥4 informative
reads, and demoting anchors below the threshold to LowCoverage.

Full chr20 sweep (diplinator truth eval; control byte-identical to HEAD and
reproduces the shipped hybrid row exactly):

| threshold | Hamming% | auN | N50 | phased | switchflip% | perfect PS |
|---|---|---|---|---|---|---|
| off (shipped) | **3.0415** | 12.00M | 951,564 | 220,897 | 1.6206 | 91 |
| 0.6 | 3.1869 | 12.00M | 951,564 | 220,828 | 1.6034 | 93 |
| 0.7 | 3.1657 | 12.03M | 954,036 | 220,794 | 1.5905 | 94 |
| 0.8 | 3.1410 | 12.00M | 951,564 | 220,714 | 1.5725 | 95 |
| 0.9 | 3.1485 | 11.97M | 951,564 | 220,692 | 1.5722 | 95 |

**Net loss at every threshold**: Hamming rises 3.04 → 3.14–3.19% while auN stays
flat. The signal *does* work in the intended direction — switchflip rate drops
monotonically (1.62 → 1.57) and perfect-PS count rises (91 → 95), so the gate
removes some genuinely mis-oriented seam anchors. But overall read-level Hamming
*rises*: the demoted anchors were, on balance, contributing more correct read
assignments than wrong orientations. Linkage concordance flags real
mis-genotypes **and** good-but-low-LD anchors indiscriminately — the same
all-or-nothing trap the AF window was tuned to 0.11 to avoid. It is not a
cleaner separator than AF; it is a noisier one. This is consistent with the
earlier finding that the wrong merges are **confident, not thin-margin**: a
filter cannot distinguish a confidently-wrong anchor from a confidently-right
one when both show strong LD with their neighbors. Reverted from code (default
0.0 reproduced the baseline byte-for-byte); kept here so it is not re-tried.
Possibly worth revisiting on ONT, where genotyping noise is larger and the
good/bad anchor populations may separate more cleanly.

# SYNTHESIS — Why the hybrid cannot beat BAM on switches/perfect-PS (read this first)

This section indexes a multi-session investigation into closing the hybrid's two
remaining deficits vs BAM (switch errors, perfect phase sets). **Conclusion: no
method tested beats the current default hybrid; the deficits are explained, not
fixable with the available data.** Every detailed section below is preserved
chronologically; this is the map.

### Current verified 3-way (chr20, vs truth)

| metric | BAM | graph | hybrid (default) |
|---|---|---|---|
| Hamming | 3.09% | 7.30% | **3.04%** ✅ |
| auN | 9.83 Mb | 0.98 Mb | **12.00 Mb** ✅ |
| N50 | 938 kb | 898 kb | **952 kb** ✅ |
| phased reads | 219,925 | 216,130 | **220,897** ✅ |
| switches | **1,393** | — | 1,440 |
| perfect PS | **125** | 92 | 91 |

The hybrid wins the metrics that matter most (accuracy, contiguity,
completeness) and loses only switches/perfect-PS.

### The one root cause behind every dead end

**GAF reads are a strict subset of BAM reads: 271,930 ⊂ 272,016, with 0 unique
reads.** The graph is a *second alignment of the same physical HiFi bases*. It
therefore carries **no evidence independent of BAM**. Its only genuine
contribution is the **3,541 extra easy SITES** BAM missed (already used by the
hybrid). Any method that re-weights or re-interprets the graph's view of *shared
reads* is information-free and capped at single-digit switches.

### Why the switch/perfect-PS deficit is structural, not a bug

Switches are counted *within* a phase set; perfect-PS is all-or-nothing per
block. The hybrid merges blocks BAM leaves split (404 vs 429 PS; auN 12.0 vs 9.8
Mb). Merging **exposes more within-block transitions** and **destroys perfect-PS
credit** when a clean block joins an imperfect neighbor — even though it adds no
discordant reads (the hybrid has *fewer*: 6,718 vs 6,798). The +47 switch gap is
localized (windows 37–38 Mb alone = +36) and traces to reads pulled into large
blocks in segdups, not to bad merge orientation (fixing all 21 wrong merges =
−3 switches, confirmed by real run).

### Ledger of everything tested (all REJECTED)

| # | Lever | Why it failed | Ceiling |
|---|---|---|---|
| 1 | k-means flip margin | no-op; switches are confident not bare-majority | 0 |
| 2 | BAM-vs-graph allele disagreement | not independent; measures alignment noise | ~0 |
| 3 | per-read LD changepoint gain | real 9× signal, but crossovers ≠ consensus switches | −29 sw, kills completeness |
| 4 | graph-as-voter / arch flip | graph 3× less accurate; no selective rule >50% | net −8,452 reads |
| 5 | stricter stitch gate (rule3+margin) | real-run −3 sw, +6 perfect-PS, −0.19 Mb auN | trade-off, not win |
| 6 | (D) graph alleles as coverage | GAF ⊂ BAM; double-counting | 0 |
| 7 | (A) graph-confirmed k-means weights | only 7 transitions fixable; loses real hets | ~3 sw |
| 8 | (C) graph-only sites as stitch bridges | only 5/40 seams reachable | <1 sw |
| 9 | (B) graph path co-occurrence linkage | recovers truth 60%/window vs BAM 97% | coin-flip |

### What would actually beat BAM

Only **data BAM does not contain**: trio/parental reads (HG002 has them — true
independent inheritance signal) or a second orthogonal sequencing technology
(e.g. ONT — independent error modes). Not another way to process the graph
alignment of the same reads.

> **PARTIALLY SUPERSEDED (post noise-filter-fix session).** This ledger was built
> *before* the indel noise-filter fixes that took graph accuracy 92.70% → 99.20%.
> Two of its premises now need re-examination:
> - "Graph is 3× less accurate" (lever #4) is **no longer true** — post-fix graph
>   is *more* accurate than BAM on shared reads (0.62% vs 1.44% discordant; wins
>   disagreements 4.8:1). See "Shared-read accuracy" in the three-pipeline section.
> - This changes the calculus for a **graph-base + BAM-fill** hybrid inversion
>   (proposed this session, see "Proposed strategy: invert hybrid"). The ledger's
>   levers all assumed BAM-base + graph-augment; an inverted architecture was never
>   tested. The "GAF ⊂ BAM / no independent read evidence" point still holds and
>   still caps gains on *shared* reads — but graph-base changes *which* pipeline
>   phases the easy core, which the ledger did not evaluate.
> The trio/ONT conclusion below still stands for beating BAM on the *shared hard
> reads*; it does not rule out the graph-base inversion, which is a different lever.

### A real, optional trade-off knob (not promoted to default)

`--stitch-rule 3 --stitch-min-margin 15` yields **+6 perfect PS** (91→97) and −3
switches at a cost of −0.19 Mb auN / −323 phased reads. Use only if perfect-PS is
the priority metric; the default maximizes the completeness win.

---

## Cross-pipeline read superset / seam-bridging (REJECTED — graph adds sites, not reads)

Hypothesis: the pangenome is more linear/contiguous and reads map better there,
so we could augment graph variation sites with confirmed BAM ones into a
cross-chunk superset, and/or let graph reads bridge the fixed-stride chunk seams
where the BAM stitch is weak. Two probes were scoped (A: widen chunk overlap so
more reads are phased on both sides of a seam; B: let graph-only bridging reads
vote in the seam stitch). Both were dropped after a data-level measurement,
**before** writing the merge/overlap code, because the enabling premise — that
the graph provides reads the BAM lacks — does not hold on HiFi.

### The graph supplies no new bridging reads

The GAF and the BAM are the **same HiFi molecules aligned two ways** (linearly
to CHM13 → BAM; to the pangenome → GAF). Genome-wide qname comparison on chr20:

| source | distinct reads |
|---|---|
| GAF | 271,930 |
| BAM (primary) | 272,016 |
| **GAF-only** | **0** (GAF ⊂ BAM) |

So the graph's contribution to the hybrid pipeline is extra **sites**, not extra
**reads**. Probe B has no new voters to add to a seam, and Probe A's premise
(seams are read-starved) is false: per 500 kb boundary there are a median of 67
bridging reads, only 2/132 boundaries below 20 crossers, and 98% have a
phaseable site within 20 kb on both sides. Adding more BAM reads or a wider
overlap recruits the same molecules that already vote.

This is config-independent (a property of the input read sets, not the phasing
parameters). Kept here so the cross-pipeline read-superset and
graph-reads-bridge-seam ideas are not re-attempted on same-molecule HiFi inputs.
May differ when the GAF and BAM come from genuinely different read sets (e.g.
graph-aligned reads that fail linear mapping entirely), which did not occur here.

> Note: an earlier write-up of this probe also reported a seam-flip "regression"
> (≈1,811 wrong-oriented reads, hybrid Hamming 4.52%, auN 13.57M). That analysis
> used a **stale `run/hybrid.phased.bam` produced by an experimental binary**
> (its auN exceeded every documented config — a tell that it was not the shipped
> default). Re-running the shipped hybrid (HEAD `034528b`, af0.11 + margin10 +
> both-strands rule 1) reproduces the documented row exactly (Hamming 3.0415%,
> auN 11,996,328, N50 951,564, 220,897 phased) and **beats BAM** (Hamming 3.04 vs
> 3.09%, 6,718 vs 6,798 discordant, auN +22%, +972 reads phased). The shipped
> hybrid is not regressed; the seam-flip figures above do not describe it. Always
> regenerate the phased BAM from the current binary before diagnosing.

## Switch-error gate on newly-phased reads (REJECTED — positive oracle ceiling, unreachable signal)

The hybrid loses to BAM on exactly two secondary metrics: switch errors
(1,410 vs 1,355, +55) and perfect phase sets (91 vs 125). Goal: close the switch
gap without giving back the hybrid's Hamming/auN advantage.

### Where the hybrid–BAM accuracy delta comes from (decomposition)

Per-read accounting that closes exactly to the −80 discordant-read delta:

| component | discordant reads |
|---|---|
| shared reads hybrid **fixed** (disc→conc) | −1,167 |
| shared reads hybrid **broke** (conc→disc) | +1,029 |
| reads only hybrid phased (the +972 net), of which 380 wrong | +380 |
| reads only BAM phased, of which 322 wrong | −322 |
| **net** | **−80 (hybrid better)** |

So the hybrid's per-read win is the sum of two large offsetting effects. Two
populations are net-harmful: the **380 discordant newly-phased reads** (real BAM
reads that BAM left unphased and hybrid phased wrong) and the **~1,012
hybrid-broke reads**.

### Oracle ceiling is real

A truth-based oracle that perfectly removes these reads clears BAM on switches:

| config | switches | vs BAM | discordant | phased |
|---|---|---|---|---|
| BAM | 1,355 | — | 6,660 | 219,610 |
| hybrid base | 1,410 | +55 (lose) | 6,629 | 220,681 |
| drop 390 disc newly-phased | 1,319 | **−36 (win)** | 6,239 | 220,291 |
| revert 1,012 hybrid-broke | 1,223 | −132 (win) | 5,617 | 220,681 |
| both | 1,113 | −242 (win) | 5,227 | 220,291 |

Dropping just the 390 discordant newly-phased reads is a win on *every* axis
(switches −91, discordant −390, only −390 phased). The headroom exists.

### But no implementable signal reaches it

The only per-read confidence available at output time is the k-means
`n_clean_agree_snps` / `n_clean_conflict_snps` (verified to survive
`mid_free_chunk`). Instrumented dump of all 220,897 phased reads, correlated
against the truth-defined harmful set:

- Harmful reads do skew low-confidence (median agree 3 vs 18 for good reads) but
  the means barely differ (18.8 vs 22.1, conflict 1.00 vs 0.08).
- **Base-rate kills it:** good reads outnumber harmful 156:1 (219,279 vs 1,402).
  The best gate precision is ~10% (`conflict≥5`: catches 48 harmful at the cost
  of 424 good); best recall ~55% destroys 37,885 good reads. No threshold or
  combination (conflict count, agree−conflict margin, low-agree) separates them.

This is the same **confident-wrong** signature as the rejected AF/linkage/stitch
levers: the harmful reads look as confidently phased as the good ones. The
oracle proves the win is *possible in principle*; the confidence signal proves it
is *not reachable* with current per-read information. Reverted (instrumentation
removed; src clean vs HEAD). Kept here so the "gate low-confidence newly-phased
reads to cut switches" idea is not re-tried without a *new* discriminative
signal (e.g. graph-vs-BAM allele agreement per read, or per-read LD with
neighbors — neither currently computed).

## Intra-chunk k-means flip margin (REJECTED — no-op; switches are confident, not bare-majority)

The hybrid's switch errors (1,440 vs BAM 1,355) are 93% intra-chunk (not at
stitch seams — switches are spread uniformly across the chunk grid, 7.1%
near-seam ≈ the 8% a uniform distribution predicts) and 95% clustered in
low-accuracy regions. Of the ~117 *isolated* (one-off) switches, 57 sit in
clean/good blocks (acc ≥0.90) and looked threshold-fixable.

The per-chunk k-means **already iterates** (10 Lloyd rounds to convergence,
`iter_update_var_hap_to_cons_alle`) and **already detects/corrects switches**
(`iter_update_var_hap_cons_phase_set` walks adjacent het-var pairs and flips the
downstream orientation when spanning-read `conflict > agree`, breaking the block
when `agree<2 && conflict<2`). So "add iterative refinement" was already done;
the isolated switches survive an iterating, switch-correcting system.

The one remaining knob was the flip *decision*: it fires on a **bare majority**
(`conflict > agree`). Hypothesis: an isolated noisy adjacent het-var pair induces
a spurious flip in an otherwise-clean block, fixable by requiring a margin
(`conflict > agree + k`). Added `--kmeans-flip-margin` (default 0); the control
(margin=0) reproduced HEAD candidates byte-for-byte (md5 4d261704). Sweep:

| margin | switch err | flip err | Hamming | auN |
|---|---|---|---|---|
| 0 (HEAD) | 1,440 | 2,133 | 3.042% | 11,996,328 |
| 1 | 1,440 | 2,132 | 3.042% | 12,026,597 |
| 2 | 1,440 | 2,132 | 3.042% | 12,026,597 |
| 3 | 1,440 | 2,132 | 3.042% | 12,026,597 |

Switches **unchanged** (1,440) at every margin; flips moved by 1. The hypothesis
is wrong: at the switch points `conflict` exceeds `agree` by far more than 3, so
the flip fires regardless of the margin. The flip rule is **confidently** making
the wrong call — same confident-wrong signature as the AF/linkage/stitch/newly-
phased-read levers. A bare-majority threshold tweak cannot help. Reverted (no-op
removed; src clean vs HEAD). Kept here so the "raise the k-means flip margin to
cut switches" idea is not re-tried — the switches are not bare-majority slips.

## BAM-vs-graph per-read allele disagreement (REJECTED — not independent; measures alignment noise, not haplotype ambiguity)

After six threshold/margin levers all hit the *confident-wrong* wall, the next
hypothesis sought a **genuinely new discriminative signal**: for each read mapped
by both pipelines, compare the BAM-derived allele to the graph-derived allele at
shared sites. The intuition was that reads whose two alignments *disagree* on
allele calls are alignment-ambiguous, and alignment ambiguity might predict the
harmful (switch-inducing) reads that k-means confidence cannot flag.

Instrumented `inject_graph_reads` (env-gated `PGPHASE_DUMP_XCHECK`) to dump, per
doubly-mapped read, the BAM allele vs graph allele at every shared site, then
joined the per-read disagreement rate against the eval `per_read.tsv.gz`
discordant/concordant labels:

| read class | mean disagree | median | % with any disagreement |
|---|---|---|---|
| concordant (good) | 0.48% | 0.00% | 4.6% |
| discordant (harmful) | 1.32% | 0.00% | 6.1% |

Both medians are 0.00%; the distributions overlap almost completely. Best gate
(disagree-rate ≥ 0.67) reaches **11.8% precision** — catches 13 of 1,069 harmful
reads while sacrificing 97 good ones. No threshold separates the classes, same
as the k-means confidence signal it was meant to replace.

**Why it fails (root cause, not just a bad threshold):** BAM and graph consume
the **same bases from the same HiFi read**, just aligned two ways. They are not
independent observations, so their disagreement measures *alignment noise*
(where the two aligners place indels/SNVs differently), which is orthogonal to
the cause of switches — *haplotype ambiguity* in segmental duplications, where
the read's bases genuinely match the other haplotype. The harmful reads are
**biologically confidently wrong**, not technically noisy, so an
alignment-noise signal cannot find them. Reverted (instrumentation removed; src
clean vs HEAD). Kept here so the "cross-check BAM vs graph alleles per read"
idea is not re-tried — the two views are not independent.

## Per-read LD changepoint gain (REJECTED — real signal, but oracle ceiling still loses on switches)

Unlike the prior six levers (all aggregate agreement counts that flag nothing),
this used a **positional** statistic. For each read, walk its phased het SNV
sites in genomic order to get a haplotype-support sequence (e.g. `1 1 1 2 2 2`),
then compute the **single-changepoint gain**: how many fewer errors a two-block
`A…A B…B` model makes vs the best one-block model. A clean within-read switch
(segdup crossover) scores high gain; scattered sequencing noise scores ~0. The
aggregate minority count used by every earlier lever cannot distinguish these.

**This is the first signal in seven attempts that actually separates the
classes:**

| signal | discordant | concordant | ratio |
|---|---|---|---|
| changepoint gain > 0 | 10.5% of reads | 1.2% of reads | ~9× |
| mean gain | 1.333 | 0.047 | ~28× |

(vs the prior levers where discordant and concordant medians/rates were
identical.) But the population is ~20:1 concordant:discordant, so even a 9×
enrichment yields only 7–26% weighted precision across the gain sweep.

**Decisive oracle (drop the flagged reads, recompute switches):**

| drop set | reads dropped | switches | flips | vs BAM (1,355) |
|---|---|---|---|---|
| baseline | 0 | 1,442 | 2,132 | +87 |
| gain ≥ 2 | 2,835 | 1,421 (−21) | 2,097 | +66 |
| gain ≥ 1 | 4,363 | 1,413 (−29) | 2,090 | +58 |

Even the most aggressive perfect-oracle drop (4,363 reads) closes only 29 of the
87-switch gap to BAM — still +58 — while sacrificing ~4,400 phased reads, which
would erase the hybrid's completeness/auN advantage (the entire reason it wins).

**Why a real signal still fails on switches:** a coincidence test showed the
within-read changepoints do **not** land at consensus switch positions — at a
500 bp window only 4–6% of high-gain discordant reads sit near a consensus
switch (≈ the concordant rate). The consensus majority vote **already absorbs**
these crossovers: the read's internal switch is outvoted before it becomes a
consensus switch, so removing the read removes its vote but rarely removes a
switch. The hybrid's remaining switch gap is a property of the **consensus over
many reads in segdups**, not of identifiable individual bad reads. No per-read
filter — however discriminative — can close it without destroying the
completeness that makes the hybrid win. Validated entirely offline (no C++
written); temp scripts removed. This closes the per-read filtering avenue: to
beat BAM on switches one must change the **consensus/phasing model in segdups**,
not filter reads.

## Graph-as-voter / graph-base-BAM-confirm architecture (REJECTED — graph is 3x less accurate, no selective override exists)

The current hybrid is **BAM-base, graph-augment**: it runs full BAM variant
calling, then adds graph snarl sites only where BAM was silent. At any **shared**
site (BAM also called it) the graph allele is **discarded**
(`extend_bam_profile_with_graph_obs` skips non-graph-only candidates). Scope on
chr20:

| site class | count | role |
|---|---|---|
| shared (BAM ∩ graph) | 63,299 | graph allele discarded |
| graph-only (extra) | 3,541 | graph's real contribution |
| BAM-only | 1,583 | graph doesn't confirm |

The proposal was to make the graph a **first-class voter** at shared sites (or
flip to graph-base/BAM-confirm), per the original intent "help phasing in the
graph with BAM-confirmed variants + extra from graph." Validated offline by
joining per-read truth labels from the standalone BAM-only and graph-only
pipelines on the 207,687 reads phased by **both**:

| | count |
|---|---|
| BAM ✓ & graph ✓ | 191,916 |
| BAM ✓ & graph ✗ (graph BREAKS) | 11,347 |
| BAM ✗ & graph ✓ (graph FIXES) | 2,895 |
| BAM ✗ & graph ✗ | 1,529 |

**Trusting graph over BAM blindly = +2,895 − 11,347 = net −8,452 reads.** The
graph is ~3× less accurate per read (6.2% vs 2.1% discordant), and the 8,443
reads only the graph phases are just 65% concordant. A graph-base flip would
wreck accuracy — which is exactly why the existing code is BAM-base.

The only way the idea survives is a **selective** override: a rule identifying
*which* shared sites/reads to trust the graph on. Exhaustive search over every
available axis — per-read mapq bands, graph hapq bands, and the eval's
BAM-bad-region BED, plus all combinations — found **no rule with >50% fix-rate at
meaningful volume**:

| stratum | fix-rate |
|---|---|
| all disagreements | 20% |
| BAM-bad region (any hapq) | 37% |
| BAM-bad region & graph hapq=60 | 18% |
| graph hapq=60 (graph max-confident) | 16% |

Even in BAM's *worst* regions the graph fixes only 37% of disagreements (net
−595). Graph's own hapq is *inversely* useful (hapq=60 → 16% fix-rate): the
graph is confidently wrong in the same hard regions BAM is.

**Root cause:** GAF reads are a strict subset of BAM reads (only 3.9% of
graph-phased reads are absent from the BAM set). The graph and BAM consume the
**same physical HiFi bases** — the graph is just a second alignment of them. So
the graph carries **no independent evidence** to overrule BAM; where BAM is wrong
(segdup, read maps to wrong copy with mapq 60) the graph inherits the same wrong
read and is wrong too. The graph's genuine value is purely the 3,541 extra
*sites* it adds (already exploited) — not a corrective vote. Validated entirely
offline (no C++ written). This closes the graph-as-voter / architecture-flip
avenue: a better signal must come from data BAM does **not** already have (e.g.
trio/parental reads, a second orthogonal sequencing technology), not from
re-weighting the graph view of the same reads.

## Wrong-merge investigation (the switch gap is NOT from bad stitch orientation)

Re-opened the "switches are structural" claim by auditing the merge decisions
against truth. The actual switch gap on chr20 is **BAM 1,393 vs hybrid 1,440 =
+47** (the earlier 1,355 figure was a stale comparison; per-phase-set tables give
1,393/1,440).

**Step 1 — are cross-block merges oriented wrong?** Identified 73 hybrid phase
sets that merge ≥2 BAM phase sets (≥5 reads each). Scored each sub-block's
orientation against truth: **21 of 73 merges are oriented WRONG.** The
discriminating signal is real and runtime-available in part:

| merge class | smallest sub-block size (median) | truth purity (median) |
|---|---|---|
| WRONG (21) | 10 reads | 0.83 |
| RIGHT (52) | 24.5 reads | 0.99 |

So wrong merges decide orientation on a small, impure set of overlap reads —
matching the vote patterns (`...(1,560,31) | (6,10)`: a 1,560-read block joined
to a 6-vs-10 coin-flip block).

**Step 2 — do the wrong merges cause the gap? NO.** Two oracles, both with truth:

- *Re-orient* every wrong sub-block to its correct truth polarity (zero
  contiguity loss): switches 1,444 → **1,441 (−3)**, flips +7.
- *Split out* small/impure weak sub-blocks: switches barely move (1,444 →
  1,437 best case), at the cost of 50–80 broken blocks.

Even a **perfect** fix of all 21 wrong merges removes only ~3 switches. Wrong
merges mostly create a small discordant *island* (a flip, consumed as 2
transitions), not a persistent switch.

**Step 3 — where is the +47 gap really?** Binned switches into 1 Mb windows.
The gap is **localized, not uniform**: windows 37–38 Mb alone account for **+36
of the +47**. Inspecting those phase sets, the mechanism is concrete: the hybrid
**adds reads/sites into already-large blocks in difficult regions**, e.g.

- PS 38393470: BAM 3,030 reads / 6 switches → hybrid 3,275 reads / **20
  switches** (the +245 reads in a segdup tripled the switches).
- PS 37005134 (hybrid): merges BAM's perfect PS 36895001 (1,676 r, 0 sw) with
  messy neighbors → 1,569 reads / **11 switches**.

Attribution of the extra switches in same-id phase sets: **325 occur in blocks
the hybrid GREW (>5% more reads); only 61 in blocks of ~same size.** The hybrid
added 14,391 reads into existing BAM blocks.

**Conclusion (corrected):** the switch gap is **not** a fixable stitch-rule bug —
fixing every wrong merge buys ~3 switches. It comes from the hybrid pulling
**more reads into large blocks within difficult/segdup regions**, where those
extra reads are genuinely harder to phase (the same biologically-confident-wrong
reads from the per-read-filter dead ends). This is the *price of the
completeness/auN win*, concentrated in a few hard windows (37–38 Mb). The only
lever that helps switches here is the same one already ruled out: stop pulling
hard reads into hard blocks — which sacrifices the completeness that makes the
hybrid win. A small, safe stitch improvement (require the weaker sub-block to
have ≥~15 reads / purity before merging) is *legitimate* but worth only a handful
of flips, not the switch gap. Validated entirely offline; no C++ written.

### Real-run confirmation: stricter stitch gate (rule 3 + margin) — oracle was exact

Two further facts narrowed the gate's reach before any run:
- Only **13 of 21** wrong-merge seams sit at 500 kb **chunk boundaries** (where
  `flip_chunk_hap` acts). The other **27 are intra-chunk** — created by k-means
  inside a chunk, which the stitch step never sees. So a stitch-gate change can
  touch at most ~13 of the wrong merges.
- The existing binary already exposes the gate via `--stitch-rule 3`
  (both-strands + margin) and `--stitch-min-margin`, so no code was needed to
  test it.

Ran the full hybrid pipeline + eval (chr20, truth BAM) at two gate settings vs
the HEAD default (rule 1, margin 10):

| metric | BAM | HEAD r1m10 | r3 m8 | r3 m15 |
|---|---|---|---|---|
| switches | 1,393 | 1,440 | 1,438 | 1,437 |
| flips | 2,052 | 2,133 | 2,131 | 2,131 |
| perfect PS | 125 | 91 | 94 | 97 |
| auN | 9.83 Mb | 12.00 Mb | 11.93 Mb | 11.81 Mb |
| Hamming | 3.09% | 3.04% | 3.04% | 3.04% |
| phased reads | 219,925 | 220,897 | 220,767 | 220,574 |

The real run matches the oracle to the read: the strictest gate moves switches by
only **−3** (1,440 → 1,437) and flips by −2. It does buy **+6 perfect phase
sets** (91 → 97, partway to BAM's 125) but at a monotonic cost in auN (−0.19 Mb)
and phased reads (−323). This is a contiguity-for-perfect-PS trade, **not** a
switch fix.

**Final verdict on the merge rule:** it is not the lever for the switch gap. The
gap is dominated by reads pulled into large blocks in a few hard windows (37–38
Mb), which no boundary stitch rule can reach. The HEAD default (rule 1, margin
10) remains the best all-round operating point — it maximizes auN/completeness,
which is the hybrid's reason to exist. If perfect-PS were the priority metric,
`--stitch-rule 3 --stitch-min-margin 15` is a valid alternative (+6 perfect PS,
−3 switches) at a small contiguity cost, but it is not promoted to default
because it sacrifices the completeness win. No source change made; sweep
artifacts removed.

## Full sweep of graph→BAM integration surfaces (all four families, REJECTED)

To answer "have we tested every way to inject graph signal into the BAM
pipeline?", enumerated the four code surfaces where graph data can enter and
validated each offline. The live hybrid already exploits surfaces 1–2 (add
graph-only sites; extend profiles at those sites). The four *untested* methods:

**D — graph alleles at shared sites as extra coverage/AF support. DEAD (no
oracle needed).** GAF reads are a strict subset of BAM: 271,930 ⊂ 272,016, with
**0** GAF-only reads. Every graph observation at a shared site is a read BAM
already counted, so "extra coverage" is double-counting — it inflates DP/AF with
zero new evidence and cannot rescue a borderline het.

**A — graph-confirmed sites as k-means anchor weights.** Real but tiny signal:
graph-unconfirmed (bam-only) het SNP sites are **1.3–1.6× more likely** to sit at
a switch than graph-confirmed shared sites — so confirmation *is* mildly
protective. But bam-only sites are only 1,583 of 64,882 (2.4%), and of 5,706
switch transitions only **7** sit where an unconfirmed site is the sole nearby
anchor (42% are at confirmed sites where weighting can't help; most are near
neither). Oracle ceiling ≈ **3 switches**, and down-weighting would discard real
het information. Not worth building.

**C — graph-only sites spanning chunk boundaries as stitch bridges. DEAD.**
Targeted the 27 intra-chunk wrong merges the stitch gate can't reach. Only
**5 of 40** wrong-merge seams have a graph-only site within 10 kb to bridge, and
fixing *all* 40 seams is worth only −3 switches (proven earlier). Ceiling < 1
switch. Root cause: bridging reads are the same physical reads (GAF ⊂ BAM), so
they replay the evidence that caused the wrong merge.

**B — graph path co-occurrence as long-range linkage (the only candidate that
could be independent). REJECTED.** The GAF node path (col 9) is the graph's
topology, not re-read bases, so it *could* carry information BAM's linear
alignment lacks. Tested by partitioning reads on graph node-bubbles and scoring
against truth across 50 windows: graph-path partition recovers the true
haplotype at only **60% per window** (median 59%, min 50% = coin flip), vs BAM's
**97%**. Concordance with BAM's own partition is 68% — but the 32% "disagreement"
is *noise*, not independent signal (it doesn't match truth). This is consistent
with the standalone graph pipeline's known accuracy (Hamming 7.3% / switch 2.58%,
~2× worse than BAM): the graph path *is* that pipeline's signal. In hard regions
the node path bifurcates on alignment ambiguity, not clean allele difference.

**Unifying conclusion.** All four fail for the **same root cause that governs the
entire investigation: GAF reads are a strict subset of BAM reads (0 unique).**
The graph is a *second alignment of the same physical HiFi bases*, so it carries
no evidence independent of BAM — not as coverage (D), not as a corrective vote
(earlier), not as linkage (B). Its only genuine contribution is the **3,541 extra
easy SITES** BAM missed, which the hybrid already uses. Surfaces that re-weight or
re-interpret the graph's *view of shared reads* (A, C, D, B, plus the earlier
graph-as-voter) are all capped at single-digit switches because there is no new
information to extract. **To beat BAM further requires data BAM does not contain
(trio/parental reads, or a second orthogonal technology) — not another way to
process the graph alignment of the same reads.** All validated offline; no C++
written; temp files removed.

## BAM-native optimizer alternatives (REJECTED — production k-means already beats MEC and transition-penalty on real matrices)

After exhausting graph→BAM integration, the question pivoted to the BAM
pipeline's **own** phasing algorithm: *is there a better optimizer than the
greedy k-means?* Tested directly against the real per-read × per-variant matrix.
**Answer: no.** Production k-means beats both MEC local search and a
transition-penalized optimizer. Details below, including a retracted earlier
claim.

### RETRACTED: the "confidence smoother" win was built on a truth-side signal

An earlier draft of this section claimed a confidence-aware smoother gated on
`hapq` cut switches −33% and Hamming +0.5pp. **That result is withdrawn.** The
`hapq` field is the `hq:i:` tag from the **diplinator truth BAM** (eval script
`evaluate_phase_accuracy.py:196-199`), not a quantity pgphase computes. It cannot
gate anything at phasing time, and "`hapq=0` ⇒ uncertain" partly describes
truth-side mapping ambiguity, not k-means tie-breaks. The lever was not
realizable; the apparent win was a measurement artifact. Lesson re-learned (same
as the LD/graph-voter dead ends): **validate the signal exists at decision time,
not just in the eval dump.**

### Honest test: dump the real matrix, run candidate optimizers offline

Added `--dump-phase-matrix PREFIX` (debug-only) to emit the exact read×variant
allele matrix k-means consumes, per chunk per flag-set, before any state
mutation (`dump_phase_matrix` in collect_phase.cpp). Joined to truth by qname.
Scoring is **PS-aware**: production assigns multiple phase sets per chunk (up to
9), each independently orientable — scoring a chunk under one global orientation
falsely inflates disc to ~80%. Within each production PS, pick the orientation
minimizing switches, then count.

Baseline reproduced correctly: **production = 1,656 switches / 3.376% disc** on
chr20 (matches the real eval-harness numbers).

### MEC local search from production init — WORSE

Starting from production's own assignment and iterating MEC consensus
re-assignment (no transition term) moved 411 reads and **increased** errors:
**+37 switches, +679 disc**. Production is already past the naive-MEC optimum —
its consensus-allele model + noisy second pass encode more than the bare matrix
energy.

### Transition-penalized optimizer — WORSE at every strength

The actual hypothesis (an HMM/WhatsHap-style switch cost). Added a contiguity
prior rewarding agreement with position-adjacent reads within each PS, swept
λ ∈ {1,2,3,5,8} from production init:

| λ | Δswitches | Δdisc | reads moved |
|---|---|---|---|
| 1 | **+123** | +1,182 | 1,835 |
| 2 | +130 | +1,166 | 1,802 |
| 3 | +537 | +3,880 | 6,313 |
| 5 | +657 | +6,629 | 9,134 |
| 8 | +715 | +7,087 | 10,796 |

Every λ is strictly worse. A contiguity prior pulls reads across **legitimate**
haplotype-block boundaries (the segdup long runs with median hapq 60 = confident
and genuinely on the other haplotype), manufacturing errors faster than it fixes
short runs. This is the same conclusion the segdup analysis predicted: the
switches that remain are **data-driven**, not optimizer tie-breaks.

### Soft EM (probabilistic assignment) — WORSE

The third candidate: replace hard assignment with soft responsibilities (EM over
a per-variant per-hap Bernoulli allele model). Initialized from production
responsibilities (not random — the earlier random-init EM was degenerate and is
not a fair test), refined to convergence, scored PS-aware:

| EM error rate ε | Δswitches | Δdisc | reads moved |
|---|---|---|---|
| 0.02 | +498 | +4,341 | 4,494 |
| 0.05 | +510 | +4,572 | 4,672 |
| 0.10 | +512 | +4,555 | 4,585 |
| 0.20 | +466 | +2,797 | 4,600 |

Worse at every error rate. The tell is a near-frozen sanity run (init
confidence 0.999, **1 round**): it already moves 3,681 reads and adds +497
switches — a single EM E-step disagrees with production on ~3,700 reads and is
net wrong. So it is not harness drift; EM genuinely disagrees and production is
right. Reason: EM's generative model (independent Bernoulli per variant per hap,
one global error rate) **discards what production encodes** — the clean-het
category weights (CleanHetSnp/Indel count double), HOM handling, and the noisy
second pass. The matrix carries the alleles but not these priors, so a model that
sees only the matrix underperforms.

### Conclusion

The greedy pivot-sweep k-means is **not** the bottleneck. On the real allele
matrix it already dominates MEC, transition-penalized, and EM/soft alternatives. The
`--dump-phase-matrix` flag is kept for future offline optimizer experiments. To
reduce BAM switches further requires *new data* (trio/parental, second
technology), consistent with the whole-investigation conclusion — not a better
optimizer over the same reads. Validated on real chr20 matrices; the only C++
added is the debug dump.

---

### (Historical, RETRACTED) original confidence-smoother write-up

The text below is preserved for the record but its conclusion is **wrong** (see
retraction above): it used the truth-side `hapq` tag as if it were a pipeline
signal.

### Root cause: BAM k-means has no transition penalty

`assign_hap_based_on_germline_het_vars_kmeans` (collect_phase.cpp:509) assigns
each read to the hap with the higher **summed** allele-match score, then runs
≤10 Lloyd rounds. Assignment is **per-read and independent** — there is no cost
for creating two transitions in the read ordering. State-of-the-art phasers
(HMM, MEC/WhatsHap) add a switch cost so a read only flips orientation when its
*own* allele evidence outweighs the cost of breaking local contiguity. BAM's
greedy hard assignment has no such term, so a read with **zero evidence margin**
(a tie) is assigned arbitrarily and can land opposite its neighbors.

### Switch breakpoints split cleanly by confidence

Decomposing BAM discordant runs by length and the pipeline's own per-read
`hapq` (haplotype quality, already computed) on chr20:

| bucket | n | median hapq | %hapq<10 | switches contributed |
|---|---|---|---|---|
| concordant | 213,127 | 60 | 2% | — |
| flip (len 1) | 1,777 | 18 | 44% | 0 (flips, not switches) |
| **short run (2–4)** | **2,061** | **0** | **61%** | **~1,551** |
| long run (5+) | 2,960 | 60 | 26% | 448 |

The dividing line is sharp:

- **Short runs (2–4 reads): median hapq = 0.** The algorithm had *no allele
  evidence* — it broke a tie and happened to flip away from neighbors. These are
  **algorithm-fixable**.
- **Long runs (5+ reads): median hapq = 60**, same as concordant reads. These
  are **confident but wrong** = true data ambiguity (segdups, where the read
  genuinely matches the other haplotype's local sequence). No optimizer fixes
  these; forcing them would *add* Hamming error. This matches the segdup
  localization (37–38 Mb) from the wrong-merge investigation.

### Truth-free smoother: −33% switches AND +0.53pp Hamming

Simulated a transition-penalizing optimizer **without using truth**: for each
read with `hapq < HQ`, reassign its orientation to the majority orientation of
its `W` nearest *confident* (`hapq ≥ HQ`) neighbors. High-confidence reads (98%
of all reads, at hapq 60) are never touched. Scored the result against truth:

| HQ< | W | switches | ΔSW | Hamming | ΔHam |
|---|---|---|---|---|---|
| baseline | | 1,999 | — | 3.091% | — |
| 20 | 6 | **1,340** | **−659 (−33%)** | **2.556%** | **+0.535pp** |
| 15 | 6 | 1,447 | −552 | 2.65% | +0.44pp |
| 10 | 6 | 1,468 | −531 | 2.745% | +0.346pp |

Safety audit (HQ<20, W6): **1,481 fixes vs 304 breaks = 4.9:1**. Every
configuration tested holds a **~4–5:1 fix:break ratio** — for each already-correct
low-confidence read wrongly flipped, 4–5 wrong reads are corrected. Because only
sub-threshold reads are ever reassigned, the confident backbone is untouched.

### Why this is different from every prior lever

Every earlier idea (per-read LD, graph-as-voter, stitch gate, the four
integration families) was a **trade-off**: fewer switches cost completeness or
Hamming. This one improves both at once because it targets reads the pipeline
*already knows* it is uncertain about (hapq≈0) — it is not introducing new
information, it is **stopping the greedy optimizer from making arbitrary
tie-breaks against local contiguity**. The signal (`hapq`) is already computed
and the operation needs only read order, so it is directly realizable in C++
inside the k-means post-pass.

**Status: validated offline, not yet implemented.** Next step (if pursued) is a
neighbor-consensus / transition-penalty post-pass in
`assign_hap_based_on_germline_het_vars_kmeans` gated on a `hapq` threshold, with
the smoother confined to sub-threshold reads. Recommended starting params from
the sweep: `HQ=20, W=6`. Must re-validate on the *real* pgphase run (synthetic
switch counts here differ from eval-harness counts; the **ratios** are what
transfer) and confirm perfect-PS and auN/N50 do not regress.

## Where hybrid's extra switches come from (within-chunk k-means, not merging)

Question: hybrid has more switches than BAM (1,393 → 1,440 in the default eval).
Is that purely the different merge rule, or do the extra graph sites change the
phasing itself?

### Test 1 — match the merge rule

BAM uses `--stitch-rule 0 --stitch-min-margin 0`; hybrid defaults to rule 1 +
margin 10. Re-ran hybrid with **BAM's exact stitch rule**:

| run | switches | flips | disc | Hamming | perfect-PS | N50 |
|---|---|---|---|---|---|---|
| BAM | 1,393 | 2,052 | 6,798 | 3.09% | 125 | 938k |
| hybrid (default rule) | 1,781 | 2,530 | 10,082 | 4.52% | 70 | 975k |
| hybrid (BAM rule) | 1,448 | 2,137 | 6,744 | 3.05% | 89 | 954k |

Matching the merge rule closes most of the gap (1,781 → 1,448) and Hamming
actually beats BAM (3.05% vs 3.09%). **Most of the default-hybrid switch excess
was the aggressive stitch rule, not the graph.** But a residual **+55 switches**
remains with identical merging.

### Test 2 — locate the residual switches

Classified the switch breakpoints present in hybrid-same but **not** in BAM
(persistent switches only, flips excluded), by distance to the nearest 500 kb
chunk boundary:

- 392 newly-introduced switch breakpoints.
- **390 (99%) are chunk-interior**, median 122 kb from any boundary.
- Only 2 are within 2 kb of a boundary.

So the residual is **not** a merge artifact — it is created **inside** chunks,
during k-means, exactly the hypothesis that adding sites changes the within-chunk
partition.

### Test 3 — mechanism: partition instability, not graph-site errors

Two further checks rule out the naive "graph sites place errors" story:

1. **No colocation with graph sites.** New switches within W bp of a graph-only
   site: 4% (150 bp), 7% (500 bp), 23% (2 kb) — all **at or below** the random
   baseline (4.7% / 15.5% / 62%). The switches do not sit at the added sites.
2. **Churn is two-way and large.** Hybrid (BAM rule) vs BAM:
   - ADDS 392 switch breakpoints
   - REMOVES 335 switch breakpoints
   - net **+57** (≈ the +55 summary delta); total churn **727**.

The net (+57) is small against the churn (727). Adding ~8k–10k graph sites does
not inject errors at specific loci; it **perturbs the global k-means consensus**,
which reshuffles ~700 breakpoints genome-wide and happens to land slightly
net-worse on switches (while improving Hamming and contiguity). This is
**optimizer sensitivity to the site set**, not a graph-quality problem.

### Answer to "would the same merge rule make it the same?"

Almost, but not exactly. Same merge rule removes ~85% of the switch gap and makes
hybrid match/beat BAM on Hamming. The remaining **+55 switches are within-chunk**:
the extra sites change the k-means partition (two-way churn of ~700, net +57).
They are not at the graph sites and not at chunk boundaries — they are a
second-order effect of re-running the same greedy k-means over a denser, slightly
different variant set. Consistent with the earlier finding that the greedy
k-means is itself the limiting factor: it is **not stable** under changes to the
candidate set, so even strictly-more sites can move switches in both directions.

## Inject graph sites only between unmergeable BAM blocks (REJECTED for switches — contiguity-only, structural)

Idea: instead of full hybrid (which churns the within-chunk partition, +57
switches), inject graph-only sites **only in the gaps between BAM phase blocks
that stayed separate** — the surgical version that should keep the contiguity
benefit without perturbing phasing. Distinct from family C (which targeted *wrong*
merge seams); this targets *clean* unmergeable gaps.

### Ceiling (offline, oracle orientation)

- BAM produces **429 phase blocks** with **128 inter-block gaps** (median 4.4 kb).
- **27 of 128 gaps contain a graph-only site** (median 2 sites/gap) — the only
  gaps a graph-site injection could bridge.
- Merging all 27 with **truth-optimal orientation** (best possible case):

| metric | baseline | oracle-merge 27 gaps | Δ |
|---|---|---|---|
| phase blocks | 429 | 402 | −27 (better contiguity) |
| within-PS switches | 1,586 | 1,737 | **+151** |

### Why merging can only hurt switches

Switch error is measured **within** a phase set. Two blocks stayed separate
because a hard region (low coverage / repeat) sits between them with no signal to
orient across. Bridging them does not remove any error — it **re-exposes that
hard region inside one phase set**, converting what were free block-boundary
discordances into counted within-PS switches. Even with perfect orientation the
operation is +switches / −blocks: a pure **contiguity-for-switches trade**.

### Conclusion

Confirms the mechanism from the hybrid-switch analysis: the graph's extra sites
buy contiguity (N50, fewer blocks) and **cannot** buy switch reduction, even when
applied surgically only between unmergeable blocks with oracle orientation. The
+151 here is the same structural effect as the full hybrid's +57, just isolated.
If contiguity is the goal, gap-injection is a cleaner lever than full hybrid
(touches only 27 sites, no global partition churn); if switches are the goal, it
cannot help. Validated offline; no production change made.

## Does injecting gap sites + re-running k-means actually merge the blocks? (NO — empirical)

The oracle test above forced merges. This tests what the **real k-means** does
when graph sites are present in the gaps — i.e. exactly what the hybrid pipeline
already runs. Compared BAM blocks against hybrid's per-read PS assignment for the
27 BAM gap-pairs that have graph-only sites in the gap.

### Result: k-means does not bridge

| outcome | count |
|---|---|
| BAM gap-pairs with graph sites in gap | 27 |
| hybrid **MERGED** the pair | **1** |
| hybrid kept **SEPARATE** | 26 |
| the 1 merge was **correctly oriented** | 0 |
| the 1 merge was **mis-oriented** (made a switch) | 1 |

The single merge was the smallest gap (1,258 bp) and k-means oriented it
**wrong**, creating a switch. The other 26 — including 8 of 9 gaps narrow enough
(<15 kb) for a HiFi read to span — were left unmerged.

### Why re-running k-means can't help

The gaps are **low-coverage / low-signal holes, not missing-site holes**. Of the
9 spannable gaps, most have only **2–3 reads inside the gap**; one has 36 reads
but a single site and still didn't merge. The block boundary exists *because* the
reads crossing it carry no usable heterozygous signal (low coverage, repeat, or a
homozygous stretch). Adding graph site **positions** does not add phasing
**signal**: those sites are genotyped by the same GAF⊂BAM reads, and if those
reads had clean spanning het alleles BAM would already have phased through. So
k-means correctly abstains on no signal — and the one time it committed, it
guessed wrong (50/50).

### Conclusion

"Inject in the gap + one more k-means round → good merge" **does not occur**:
1/27 merged, and that one was mis-oriented. This is the empirical confirmation of
the oracle finding — and of the whole investigation's root cause. Bridging
unmergeable blocks needs **independent linking evidence** (spanning reads with
het signal, i.e. new data), which graph sites over the same reads cannot provide.
No production change made; validated offline.

## Are there regions the graph phases BETTER than BAM? (YES — but not truth-free identifiable)

User intuition: the "GAF ⊂ BAM, 0 unique reads" finding is about read *identity*;
it does not preclude *regions* where the graph phases better (a read can be in
both but soft-clipped/misplaced in BAM yet clean in the graph). Tested directly.

### Methodology fix (important)

First attempts had a scoring bug: best-flipping a 100 kb *window* as a whole vs
best-flipping each *phase set* gave contradictory answers (a window with two
opposite-oriented PS scores 50% under window-flip but 0% under per-PS flip). The
"159 windows graph wins" and "BAM 50%/graph 0%" intermediate numbers were
artifacts and are **retracted**. Correct metric: best-flip **each PS
independently**, label every read correct/incorrect (flip-invariant), then
compare pipelines on shared reads.

### Result: graph-favorable regions are real

Reads phased by **both** pipelines (207,687):

- GRAPH right & BAM wrong: **25,745 reads**
- BAM right & GRAPH wrong: **37,556 reads**

Net favors BAM (matches the aggregate Hamming 3.1% vs 7.3%), but the two-way
disagreement is large and there is a **genuine graph-favorable subset**. Per
100 kb window (≥20 shared reads): **142 windows graph wins, 205 BAM wins, 293
tie**. In the graph-win windows BAM is internally noisy (a mid-block switch);
the graph holds a clean single block. So the intuition is **correct**: the graph
does phase some regions better.

Separately, **8,443 reads are phased by graph but not BAM at all**; of those,
~1,154 (39 PS) are well-phased (<10% disc), the rest poorly. So the graph also
adds a small amount of correctly-phased coverage BAM misses.

### But the regions are not truth-free identifiable — so not exploitable

The blocker for a selective hybrid: can we tell *which* regions to trust the
graph on **without** truth? Tested BAM's own internal consistency (fraction of
window reads agreeing with the window-majority HP) as a router signal:

| window class | n | median BAM internal consistency |
|---|---|---|
| graph wins | 142 | 0.64 |
| BAM wins | 205 | 0.58 |
| tie | 293 | 0.55 |

The signal **does not separate** graph-win from BAM-win windows (0.64 vs 0.58,
both near the 0.5 split floor). BAM looks equally unsure in regions where it is
actually right and where it is actually wrong. So a selective hybrid has no
reliable truth-free rule to defer to the graph — blindly merging pulls in the
graph's 205 losses with its 142 wins, which is exactly why every whole-pipeline
merge tested lands net-worse on switches.

### Conclusion

Both can be true at once, and both are: (1) the graph phases a real subset of
regions better than BAM (~142 windows / ~25k reads), validating the intuition;
(2) those regions cannot be identified from BAM-side signal without truth, so the
gain is not capturable by a router. The remaining path to exploit them is an
**independent confidence signal** — e.g. the graph pipeline emitting a calibrated
per-PS phasing-quality score that is trustworthy where it disagrees with BAM —
which would require validating that the graph *knows* when it is right (its own
confidence vs truth). That is the concrete next experiment if this thread is
pursued. Validated offline; scoring bug noted and corrected; no production change.

## Can a confidence signal route reads to the graph where it wins? (NO — gain is 3 lucky regions, not a generalizable signal)

Follow-up to the previous section: the graph phases ~142 windows / ~25k reads
better than BAM, but those regions were not separable by a single BAM-side
signal. This tests the full question — can **any** truth-free signal (graph-side,
BAM-side, or a learned combination) identify where to trust the graph, enough to
build a confidence-routed hybrid?

### Setup

Trusted labels: best-flip each PS, per-read correct/incorrect (flip-invariant).
Focus on the **63,301 disagreement reads** (graph and BAM give different
correctness) — the only reads a router decision affects. Base rate: 40.7%
graph-right. Target: predict graph-right without truth, route those reads to the
graph.

### Single signals — all useless

Truth-free signals tested (per graph PS and per read): PS size, het-site density
(5k/20k windows, graph and BAM), allele-frequency-near-0.5, mean depth, read
mapq, read coverage, BAM/graph orientation agreement, position offset.

| signal | AUC |
|---|---|
| read mapq | 0.50 |
| graph PS size | 0.44 |
| BAM PS size | 0.52 |
| het-density diff | 0.53 |
| all others | 0.48–0.53 |

Every single signal is at chance. Best linear combination (logistic regression,
held-out): **AUC 0.60** — below the 0.65 usability bar.

### Gradient boosting — 0.894 was leakage

A gradient-boosted model scored **AUC 0.894** with read-level 4-fold CV — but
this is **train/test leakage**: reads from one PS share a PS-size feature and
(mostly) one label, so reads of the same block leak across folds. Re-run with
**GroupKFold over phase sets** (no PS in both train and test): **AUC 0.59 ±
0.13**. The boosted model was memorizing which specific block is right, not
learning a transferable signal. (Sanity: the top features `gps_n`/`bps_n` have
individual AUC 0.44/0.52 — useless alone, confirming the 0.894 came from
group structure, not signal.)

### The gain is 3 regions, not a rule

Translating the grouped out-of-fold predictions into routing outcomes:

| threshold | routed | fixed | broke | net | precision |
|---|---|---|---|---|---|
| 0.5 | 20,280 | 9,978 | 10,302 | −324 | 0.49 |
| 0.7 | 7,038 | 4,074 | 2,964 | +1,110 | 0.58 |
| 0.8 | 3,536 | 2,646 | 890 | +1,756 | 0.75 |
| 0.9 | 1,180 | 1,177 | 3 | +1,174 | 1.00 |

The high-precision tail looks great — until you check **where it comes from**.
At thr 0.9, **1,151 of 1,318 reads are a single phase set** (PS 2717436, graph
disc 0.00). Decomposing net gain per PS at thr 0.8:

- **Top 3 PS: +1,963 net = 158% of the total** (1,244).
- **All other 24 contributing PS: −719 net** (net negative).
- **AUC excluding the top 3 PS: 0.530** (chance).

So the entire apparent gain is **3 lucky large blocks** that happen to be
graph-correct and present in the data. Everywhere else the router breaks more
reads than it fixes. This is memorization of specific regions, not a
generalizable confidence signal.

### Conclusion

**The graph does not know when it is right.** No truth-free signal — single,
linear, or gradient-boosted with PS-grouped validation — generalizes beyond a
handful of specific regions (AUC 0.53 once those are removed). A confidence-routed
hybrid is therefore not realizable from the signals available in the current
pipeline outputs: it would either need (a) a genuinely calibrated per-PS quality
score computed inside the graph phaser and validated to transfer across regions,
or (b) independent linking data. This closes the last open thread: the graph's
region-level wins are **real but not capturable** without truth. Validated
offline (sklearn GroupKFold); no production change. Methodological note for future
work: **always group-split by phase set** — read-level CV leaks and inflated this
AUC from 0.59 to 0.89.

## Why does the graph phase some regions better? (NOT more sites — k-means convergence coin-flip)

Having established the graph wins ~141 windows but those wins aren't routable,
the remaining question is *why* it wins them. Tested every data-side hypothesis on
the three window classes (graph-win=141, bam-win=204, neutral=293, 100 kb,
PS-aware, ≥30 shared reads).

### It is NOT more sites, coverage, or quality

| class | graph het (med) | bam het (med) | g−b het | bam cov | graph cov | bam mapq | graph mapq |
|---|---|---|---|---|---|---|---|
| GRAPH-WIN | 103 | 98 | **−1** | 332 | 328 | 60 | 60 |
| BAM-WIN | 108 | 108 | −3 | 350 | 352 | 60 | 60 |
| NEUTRAL | 120 | 123 | −2 | 372 | 372 | 60 | 60 |

In graph-win windows the graph has **equal or slightly fewer** het sites, the same
coverage, and the same mapq as BAM. The graph never has more sites anywhere — it
has marginally fewer everywhere. **Site count, coverage, and read quality are
ruled out.** The two pipelines see the same reads at the same sites with the same
quality, yet phase differently.

### The mechanism: each pipeline's k-means coin-flips convergence

Counting within-window switches per pipeline:

| class | BAM switches/win | graph switches/win |
|---|---|---|
| GRAPH-WIN | 9 (clean) | 60 (noisy) → wait, graph WINS here |
| BAM-WIN | 68 (noisy) | 10 (clean) |
| NEUTRAL | 60 | 60 |

(The disc-rate winner is what defines the class; the switch counts show the
*loser* in each class has a mid-block switch.) The pattern is **symmetric**: in
graph-win windows BAM carries a switch the graph avoided; in bam-win windows the
graph carries a switch BAM avoided; in neutral windows both are equally noisy.

Decomposing the **loser's** error shape:

- GRAPH-WIN: BAM's error is a **clean switch** (one contiguous wrong block ≈ a
  half-window orientation flip) in **83/141** windows; scattered noise in 58.
- BAM-WIN: graph's error is a clean switch in **120/204**; scattered in 84.

A "clean switch" means the k-means oriented half the window the wrong way — a
**convergence failure**, not ambiguous data (the other pipeline got the same
evidence right). 

### Conclusion

The graph's region-level wins are **not a data advantage** — same sites, same
reads, same coverage, same quality. They are the **same fragile greedy k-means
landing a different coin-flip**: in each hard window, each pipeline independently
either converges to the clean haplotype partition or to a half-flipped one, and
"graph-win" is simply where the graph's flip landed right and BAM's landed wrong
(and symmetrically for bam-win). This is the same k-means instability behind the
+57 hybrid churn, now shown to drive the region-level differences too. It explains
both earlier findings at once: the wins are **real** (genuine convergence
differences) but **not routable** (a coin-flip has no truth-free tell). The only
fixes are a **more stable optimizer** (shown earlier that MEC/EM/transition-penalty
don't beat the current k-means on the real matrix) or **independent data**.
Validated offline; no production change.

---

## Does a more stable optimizer beat production? (multi-restart + MEC selection)

**Motivation.** The coin-flip mechanism above suggests the win is restart-luck.
If so, running many random-pivot restarts and **selecting by a truth-free
objective** (MEC energy) should land the good basin more often and cut switches.
Tested on the real dumped matrix (133 chunks, flags140, PS-aware best-flip
scoring per phase set). Production baseline = **1,656 switches / 3.376% disc**.

### 1. Convergence variance is real and large
8 random pivot inits per chunk on the first 20 chunks: per-chunk switch count
ranges **26–92** by init alone (e.g. chunk 112: 57→149); total-switch range 941
across 20 chunks. MEC energy also varies per run (spread 944–3934), so a
truth-free objective *can* in principle distinguish runs.

### 2. But restarting the reimplemented k-means is far worse
12 restarts, pick min-MEC, scored from scratch:

| optimizer | switches | disc |
|---|---|---|
| production | 1,656 | 3.376% |
| multi-restart ×12, min-MEC | 5,473 | 43.3% |
| oracle (min-**switch** run) | 4,432 | — |

The reimplemented `pivot_init`/`kmeans` is a **degenerate baseline** — even its
oracle (4,432) can't reach production (1,656). Production's category weights,
HOM handling, and locked-pivot Phase-1 init are doing real work that a plain
Lloyd restart discards. Restarting a worse optimizer can't beat a better one.

### 3. MEC is a weak proxy for switch count
Within-chunk Spearman ρ(MEC, switches) over 40 chunks: **mean 0.565**, median
0.594, positive in 97%. But argmin(MEC) == argmin(switches) in only **13/40**
chunks. MEC correlates with truth switches but selects the truly-best run only
~⅓ of the time — too weak to be a reliable selector.

### 4. Deployment-realistic test: production never gets overridden
Include production's own assignment as one candidate alongside 12 restarts,
select min-MEC:

| | switches |
|---|---|
| production | 1,656 |
| prod + 12 restarts, min-MEC | 1,658 (**+2**) |

MEC swapped away from production in only **2/133** chunks, and one of those two
made switches worse. Net effect ≈ zero. **A more stable optimizer provides no
truth-scored gain** — production already sits at the good basin MEC would pick.

### 5. Ensemble instability does not flag errors
Per-read vote fraction across 12 aligned restarts; "unstable" = vote fraction in
(0.2, 0.8). Unstable reads are discordant **48.6%** vs stable **46.3%** —
essentially no separation (flag precision 0.486 ≈ base rate). Restart-instability
carries **no truth-free signal** about which reads are mis-phased, so a
consensus/confidence scheme built on it cannot route or correct errors.

### Conclusion
A more stable optimizer does **not** beat production k-means. (1) Variance is
real, but (2) the only optimizer that beats restarts *is* production, (3) MEC is
too weak a selector to improve on it, (4) when offered production as a candidate
MEC keeps it (net +2), and (5) restart-instability is not a usable error signal.
Combined with the earlier MEC/EM/transition-penalty results, every truth-free
optimizer variant tried fails to beat the current k-means. The coin-flip wins
are real but not capturable without **independent data**. Validated offline; no
production change. The `--dump-phase-matrix` debug flag remains the only `src/`
change and alters no phasing behavior.

---

## Can graph read-mappings link BAM phaseblocks across gaps?

**Idea (independent-data angle).** Phase with BAM, find the gaps between adjacent
BAM phaseblocks that couldn't be merged, and ask whether reads' **graph (GAF)
mappings** span those gaps even where their **BAM alignments** break. A read whose
linear BAM alignment clips at a repeat/SV but whose graph path threads through
could supply linkage BAM lacks — turning two blocks into one. This is distinct
from earlier gap-injection (which injected graph *sites* into k-means); here we
test graph *read connectivity* directly.

### Setup
- BAM phaseblocks from `bam.phased.vcf`: **431 PS blocks, 430 inter-block gaps,
  ~11 Mb total gap.** Largest gaps cluster in the peri-centromere (27–32 Mb).
- GAF: 271,930 alignments with linear ref projection (contig, ref_start, ref_end).
- A read "spans" a gap a→b if some alignment has ref_start≤a and ref_end≥b.
- "graph-only" spanner = spans in GAF but **no** BAM alignment of that read spans.

### Coordinate spanning exists, but is a centromeric artifact
287/430 gaps have ≥1 GAF spanner. Classifying by region (centromere = 26.0–29.5 Mb):

| region | n_gaps | gaps with graph-only linkage | total graph-only reads | total GAF spanners |
|---|---|---|---|---|
| **arm** | 400 | **5** | **65** | 6,288 |
| **cent** | 30 | 17 | **1,534** | 1,885 |

**96% of all graph-only linkage (1,534/1,599 reads) is inside the centromere.**
There, "spanning" is a linear-projection artifact: the graph path threads a
tandem-repeat array that the BAM aligner correctly refused to align across, and
truth alignment itself fails (the diplinator can't place these reads either).
Spanners also scatter across many production PS labels rather than anchoring the
two specific flanking blocks (e.g. 393 spanners of the 26.567→26.591 Mb gap, only
6 with any truth hap, only 3 assigned to a flanking PS).

### On the arms, there is no usable linkage
Only **4 arm gaps** have ≥3 graph-only spanners. Examining each for truth-consistent
orientation (a merge needs spanners that agree on the relative phase of the two
blocks):

| arm gap | graph-only | truth-hap split | no-truth / unaligned |
|---|---|---|---|
| 49,661,443→49,661,903 (460 bp) | 23 | **9 MAT / 9 PAT** | 5 |
| 31,765,624→31,766,466 (842 bp) | 20 | 3 MAT | **17** |
| 31,780,986→31,909,439 (128 kb) | 13 | 3 MAT | **10** |
| 29,506,920→29,510,999 (4 kb) | 7 | 1 PAT | **6** |

- The one gap with both haplotypes present (49.66 Mb) splits **9 MAT / 9 PAT** —
  perfectly mixed, which gives **zero** orientation signal (linkage requires the
  spanners to consistently tie the blocks in one relative phase, not 50/50).
- The other three are dominated by **no-truth/unaligned** reads: the graph spans
  linearly but truth can't place the reads, so the "linkage" crosses a repeat
  where the projection is unreliable.

### Conclusion
Mapping BAM phaseblocks to the graph and inspecting inter-block read mappings does
**not** yield a usable merge signal. The graph-only spanning is real but ~entirely
centromeric, where it reflects a linear-projection artifact over repeats rather
than true connectivity — exactly the regions where BAM's refusal to span is
*correct* and truth itself is undefined. On the chromosome arms, where phasing is
meaningful, graph-only spanners are too few (4 gaps) and carry no consistent
orientation (mixed haplotypes or unaligned reads). This is consistent with the
established GAF ⊂ BAM finding: the graph is a second alignment of the same HiFi
bases, so outside repeat regions it cannot connect what BAM cannot. The remaining
graph wins stay attributable to k-means coin-flips, not extra linkage. Validated
offline; no production change.

---

## Do reads map *better* in the graph (CHM13 ≠ sample)? Alignment-quality rescue

**The objection.** The reads are HG002 HiFi mapped to CHM13 — a reference that is
**not** the sample. Where HG002 diverges from CHM13 (SVs, divergent/segdup
haplotypes), the *linear* BAM alignment may be present but **wrong** (clipped,
mis-placed), while the graph threads the correct alternate path. "GAF ⊂ BAM" is
about read *identity*; it doesn't rule out the graph *aligning the same read
better*. This is a sharper test than the earlier coordinate-spanning one, which
only checked whether BAM physically reached across a gap, not whether its
alignment there was correct.

### Methodology note (a real bug I made and caught)
First pass read the wrong GAF columns: this file is a **custom format with 3
extra leading columns** (contig, ref_start, ref_end) before the standard GAF
block, so matches/blocklen/mapq are at cols **13/14/15**, not 10/11/12. The wrong
read produced a fake "18,063 reads where BAM mapq<20 but GAF mapq≥60." After
fixing the column offset, the result inverted (see below). Lesson: always sanity-
check that a "mapq" column actually ranges 0–60.

### GAF and BAM mapq are byte-for-byte identical
With the **correct** columns:

> **GAF mapq == BAM mapq for 271,930 / 271,930 reads (100.0%).**

The graph alignment is not an independent remap — it carries the *same minimap2
mapq* as the BAM, because it is the same alignment re-projected onto the graph.
So there is **no** read where the graph is more *confidently* placed than BAM
(0 reads with GAF mapq ≥ BAM mapq + 20, and 0 the other way). This is the
strongest form of GAF ⊂ BAM yet: same reads, same alignments, same confidence.

### One real residual: BAM-clipped, graph-clean reads
The graph *can* still differ in **alignment extent**: 792 reads are soft-clipped
>15% in BAM yet align ≥98% identity in GAF (637 on the chromosome arms). These
are real CHM13-vs-HG002 divergence loci where the graph's alt path lets the read
align full-length. Tracing them through the pipelines and truth:

| | count |
|---|---|
| BAM-clipped>15% & GAF identity>98%, on arm | 637 |
| …phased **concordantly** by the graph pipeline | 184 (100% concordant) |
| …phased by graph but **not** by BAM at all | **39** (all concordant) |
| …among reads BAM did phase, BAM **discordant** / graph correct | 4 |
| …of the 39 graph-only rescues, captured by current **hybrid** | **1** |

So the objection is **partially right**: there is a genuine, truth-validated set
of ~39 arm reads that the graph phases correctly because its alt path aligns them
where BAM clips — and BAM phases none of them. They sit at **31–32 Mb**, exactly
among the large BAM phaseblock gaps (128 kb, 151 kb, 170 kb …) in the segdup-rich
peri-centromeric edge, the region where CHM13≠HG002 divergence is highest.

### But the magnitude is tiny and the hybrid already misses it
- **39 reads** total, all in one ~1 Mb segdup band — against 1,393 BAM switches
  and ~223k phased reads. Even perfectly exploited, this cannot move the chr20
  switch/Hamming numbers measurably.
- The current hybrid pipeline captures **1 of 39**, so the mechanism is real but
  not wired in. Exploiting it would mean preferring the graph's alignment extent
  (not its mapq, which is identical) for clipped reads in divergent loci.

### Conclusion
The "reads map better in the graph" claim is **correct in kind but small in
degree** on this CHM13/HG002 chr20 data. Mapping *confidence* is identical (same
minimap2 mapq, 100%), so the graph never rescues via better mapq. The only real
edge is **alignment extent**: ~39 arm reads that BAM soft-clips but the graph
aligns full-length and phases correctly, clustered in the 31–32 Mb segdup band
near gaps BAM cannot close. This is a genuine independent-data signal — the first
found in this whole investigation — but at 39/223,000 reads it is far too small to
shift aggregate accuracy, and the current hybrid captures only 1. It would matter
more on a sample/region with heavier structural divergence from the reference than
chr20 offers here. Validated offline; no production change.

---

## Can the graph add correct phaseblocks where BAM phases nothing?

**Idea.** Don't fix BAM regions — *add* phaseblocks in regions BAM leaves
unphased. Even a few extra correct blocks is a net win (improves N50/auN, phases
reads BAM drops). Test: find graph PS blocks lying entirely outside BAM's phased
coverage, then check them against truth.

### Graph-only blocks exist
Merging BAM's 431 phased blocks into coverage intervals, **45 graph PS blocks
fall entirely outside** BAM coverage; **17 carry ≥2 phased het vars** (the rest
are singletons). Several are substantial (e.g. 26.87 Mb / 36 vars, 32.23 Mb / 23
vars, 44.03 Mb / 15 kb span).

### But raw graph-only blocks are only ~⅓ correct
Scoring each against truth (best-flip disc over its phased reads):

> **5 / 17 are perfect (100% concordant); the rest range down to 50–60%
> accuracy — i.e. coin-flip garbage.**

The bad blocks are the same fragile-k-means failure seen elsewhere: high-read
blocks (146, 66, 63 reads) at ~55–60% accuracy are half-flipped partitions, not
real haplotypes. Admitting them blindly would *add wrong phase*, which is worse
than leaving the region unphased.

### A truth-free filter recovers a clean subset
Per-block signals vs truth accuracy revealed the discriminators:

| signal | finding |
|---|---|
| **indel fraction** | every coin-flip block (acc 0.53–0.60) is **100% indels** with high DP (33–77) — repeat/homopolymer false-het piles |
| **degenerate span** | several bad blocks have span 0–1 bp with many reads (var/kb 999–2000) — multiple vars stacked at one position |
| **segdup band** | the remaining SNP-rich bad blocks all cluster in **32.2–32.5 Mb**, a known CHM13 peri-centromeric segmental-duplication band where paralogs create false hets |
| ~~strand bias, AF~~ | **RETRACTED — see bug below.** AF doesn't separate paralog-het from real het, but the "strand bias useless" claim was wrong: it was an artifact of a strand-accounting bug, not biology |

Applying **"≥2 SNPs AND not in the 32.2–32.5 Mb segdup band"**:

> **4 extra correct phaseblocks, 1 bad admitted, 15 reads newly phased.**

(Recall of the perfect blocks is 4/5; the SNP filter alone gives recall 5/5 but
admits 4 bad, so the segdup exclusion is what makes it usable.)

### Conclusion
The idea **works in principle and yields a real, if small, win**: ~4 extra
truth-correct phaseblocks on chr20 that BAM does not phase at all, recoverable
with a simple truth-free filter (require ≥2 SNPs, exclude the pure-indel and
segdup-band blocks). This is the same independent-data edge as the clipped-read
rescue — concentrated in the 26–32 Mb peri-centromeric region where CHM13≠HG002
divergence is highest — and again **tiny in magnitude** (4 blocks / 15 reads
against 431 BAM blocks). Unlike the optimizer and linkage avenues (which produced
*nothing*), this one is a genuine, defensible net gain and the most promising
hybrid lever found: graph-only blocks are additive (they can't introduce switches
into existing BAM blocks since they occupy disjoint regions), so the only risk is
admitting a wrong block, which the filter controls. A production implementation
would: (1) emit graph PS blocks, (2) keep only those disjoint from BAM coverage
with ≥2 SNPs and outside known segdup/centromere bands, (3) append them to the BAM
phasing. Worth doing for N50/auN even though aggregate switch/Hamming barely move.
Validated offline; no production change yet.

### BUG FOUND (RESOLVED in 8f93774): graph candidate strand counts were always forward (REVERSE=0)

> **Status: FIXED.** Resolved by commit 8f93774 in `graph_collect.cpp:148-153`
> (copy `mcand.counts.{forward,reverse}_{ref,alt}` instead of deriving
> forward-only from `alle_covs`). Verified on current HEAD: all 74,119 rows of
> `graph.candidates.tsv` satisfy `FWD+REV == COUNT` for both ref and alt, including
> multiallelic sites. The strand-bias filter remains correctly `is_ont()`-gated, so
> it does not run on the HiFi path. The historical investigation below is retained
> for the diagnostic trail.


While testing a strand-bias filter for the graph-only blocks, every graph
candidate showed `REVERSE_REF=0` and `REVERSE_ALT=0`. Investigation showed this
is **not biology** (HiFi reads map to both strands — the GAF path orientation is a
real ~49/51 split, 139,013 `>` forward / 132,917 `<` reverse) but a **bug in the
output-table rebuild**:

- Strand is parsed correctly from **path orientation** (`graph_query.cpp:56`,
  `reverse = orient == '<'`), *not* the GAF strand column (which is always `+`).
- Chunk-level accumulation is correct: `graph_bam_adapter.cpp:444-466` splits
  forward/reverse counts properly, and the **strand-bias filter does run** on them
  (`graph_bam_adapter.cpp:179-190`, Fisher-exact, `is_ont()`-gated — mirrors the
  BAM filter at `collect_var.cpp:1475-1482`).
- **The defect:** `graph_chunks_to_candidate_table` (`graph_collect.cpp:146-147`)
  rebuilds candidates from strand-**merged** `alle_covs` and sets
  `forward_ref = ref_cov; forward_alt = alt_cov`, **never setting reverse**. So the
  emitted `graph.candidates.tsv` zeroes all reverse counts.

This is the **same lossy-rebuild pattern** as the `RepeatHetIndel` re-promotion bug
(also in `graph_chunks_to_candidate_table`): the function discards category and
strand information computed upstream.

**Consequence for the analysis above:** my offline "strand bias is useless"
conclusion was drawn on the corrupted TSV (every variant looked single-strand), so
it is **invalid**. The strand-bias filter that runs *during* phasing uses correct
counts, but I could not evaluate strand bias as a graph-only-block discriminator
because the dumped data was wrong. The indel-fraction and segdup-band signals
stand; the strand-bias verdict is retracted pending a re-test on fixed output.

**Fix (one line):** in `graph_collect.cpp:146-147`, copy
`mcand.counts.forward_ref/reverse_ref/forward_alt/reverse_alt` (already strand-
correct after biallelic pairing) instead of deriving forward-only from
`alle_covs`. This is a debug-output/VCF-annotation correctness fix; it does not
change k-means (which uses the chunk counts, not the rebuilt table) but it does
affect the emitted strand-bias-derived VCF fields and any downstream strand
filtering on the TSV.

### APPLIED + RE-TESTED: strand fix and the corrected strand-bias verdict

Applied the one-line fix and regenerated `graph.candidates.tsv` + the graph eval
dump. Verification:

- Reverse counts now populated; only **0.43%** of het variants are zero-reverse
  (genuine single-strand sites), down from **100%**.
- Graph strand counts now **match the BAM pipeline** at shared positions (e.g.
  pos 59818: both Fref=8 Rref=5 Falt=7 Ralt=13).
- Site count **unchanged** (74,119) — confirms the fix is output-only and does
  not alter phasing.

**Strand-bias verdict (now valid):** Re-running per-variant Fisher-exact strand
bias on the **corrected** counts across all 17 graph-only blocks flags **zero**
variants (every block min-p > 0.09). So strand bias is **not** a discriminator for
the graph-only coin-flip blocks. The earlier "strand bias useless" conclusion was
right in outcome but had been drawn on corrupted data; it is now confirmed on
correct data. (Strand bias remains a valid filter in general — it just doesn't
separate *these* blocks, which fail for a different reason.)

### The real discriminator: the RepeatHetIndel re-promotion bug

Inspecting the worst coin-flip blocks (acc 0.53–0.60 at 40–44 Mb) directly in the
VCF shows they are **pure homopolymer / STR indels** emitted as `CLEAN`:

```
40,685,510  C>CA, C>CAA              (poly-A insertion)
42,042,659  A>ATA, A>ATATATATA       (TA-repeat expansion)
43,797,938  C>CT, C>CTT              (poly-T)
44,300,880  C>CAA, C>CAAAA, C>CAAAAAA, C>CAAAAAAA  (poly-A run, same pos)
44,032,221  AA>A, GG>G               (homopolymer deletions)
```

These are textbook `is_homopolymer_indel` / `is_repeat_indel` hits. The noise
filter **does** demote them to `RepeatHetIndel` during phasing (so they're
excluded from k-means and never properly phased — hence the ~50% coin-flip
accuracy), but `graph_chunks_to_candidate_table` rebuilds the output category
from scratch (`graph_collect.cpp:170-180`, classifying purely by type/AF/depth)
and **re-promotes them to `CleanHetIndel`**, so they pass the VCF germline gate
and surface as phased het indels with a (default) PS. This is the same lossy-
rebuild defect as the strand bug, in the same function.

**This is the actual filter the indel-fraction signal was detecting.** The
"100%-indel coin-flip block" pattern is precisely the re-promoted repeat indels.
The principled fix is not a strand or indel-fraction heuristic but **preserving
the `RepeatHetIndel` demotion** through the rebuild — then these 7 bad blocks
never enter the graph VCF at all.

### Updated recommendation for graph-only-block rescue
Breakdown of the 17 multi-var graph-only blocks by failure mode:

| failure mode | blocks | fix |
|---|---|---|
| pure repeat/homopolymer indel (acc ~0.5–0.6) | ~7 | preserve `RepeatHetIndel` demotion in rebuild |
| SNP-rich, in 32.2–32.5 Mb segdup band (acc 0.56–0.78) | 3 | segdup/centromere band exclusion |
| clean SNP-bearing, correct (acc 1.0) | ~5 | **keep — these are the win** |

So the production path is cleaner than the earlier "≥2 SNPs + segdup" heuristic:
(1) **fix the RepeatHetIndel re-promotion** (removes the pure-indel coin-flips at
the source, a correctness fix that also benefits the main graph VCF), (2) **exclude
known segdup/centromere bands**, (3) **append remaining graph-only blocks** to the
BAM phasing. Strand bias is *not* needed for this. The strand fix is still a valid
independent correctness fix and is applied. Validated offline.

### APPLIED: RepeatHetIndel re-promotion fix (+ two underlying noise-filter bugs)

Implemented the fix. It required correcting **three** distinct defects, all of
which let homopolymer/STR indels reach the phased graph VCF:

**Bug A — category dropped on rebuild (`graph_collect.cpp`).**
`graph_chunks_to_candidate_table` re-classified every indel het as
`CleanHetIndel` from scratch, discarding the noise filter's `RepeatHetIndel`
demotion. Fixed by preserving the demotion: if the chunk candidate is
`RepeatHetIndel`, carry it (and `kLongcalldRepHetVar`) through the rebuild instead
of re-promoting.

**Bug B — noise filter never fired (`graph_bam_adapter.cpp` + `noise_filter.cpp`).**
Even before the rebuild, the demotion was rarely set, for two reasons:
  1. *Non-minimal graph alleles.* Graph emits indels with the full repeat run on
     both sides (e.g. `CAAAAAAAAAAA > CAAAAAAAA` for a 3 bp deletion).
     `is_noisy_site` derives the indel length from allele sizes, so the length
     exceeded `max_xgaps` and the homopolymer scan was skipped. **Fix:** trim each
     allele pair to minimal VCF form (shared suffix then prefix) before the check.
  2. *Off-by-one in the shared scan.* `is_homopolymer_indel` / `is_repeat_indel`
     read the reference context starting at the VCF **anchor** base (`end_pos =
     pos`) instead of the inserted/deleted content one base later. For `C>CAAA`
     over reference `...c|aaaa…`, the scan saw the anchor `c` and missed the
     poly-A run. Verified empirically (`end_pos=pos` → not detected; `end_pos=pos+1`
     → detected) for all failing cases. **Fix:** start the downstream scan at
     `pos+1` (insertion) / `pos+del_len+1` (deletion), upstream at `pos`; also made
     the `is_repeat_indel` insertion comparison nt4-encoded so soft-masked
     (lowercase) reference matches. *This shared function is used only by the graph
     and hybrid-inject paths; the BAM pipeline has its own local
     `var_is_homopolymer_pg`, so this change does not affect BAM.*

**Bug C — only the first alt checked (`graph_bam_adapter.cpp`).**
Multiallelic sites (e.g. pos 33,887,895 with 9 poly-A alts) were tested only on
`meta.alts[0]`, which after decomposition need not be the phased allele — and the
longest alt could exceed `max_xgaps`. **Fix:** check **every** alt and demote the
site if any is noisy.

**Result on the graph pipeline (chr20):**

| metric | before (buggy) | after (fixed) |
|---|---|---|
| Hamming error | 1.03% | **0.81%** |
| Switch error | 0.24% | **0.18%** |
| switches | 492 | **368** |
| flips | 1,025 | **859** |
| perfect phase sets | 37.8% | **39.4%** |
| variants demoted to `REP_HET_INDEL` | (re-promoted to clean) | **16,841** |

Graph-only multi-var blocks dropped from 17 → 13; the four egregious pure
poly-A/poly-T/TA-repeat coin-flip blocks at 40–44 Mb (acc ~0.5–0.6) are gone, and
all five truth-correct clean SNP blocks remain. **7 of 8** previously-leaking
repeat-indel blocks are now filtered.

**One residual — RE-DIAGNOSED as a normalization bug, now FIXED:** pos 49,031,428
was previously described here as a "~40 bp AGGG STR expansion." That was wrong.
Trimming the catalog alleles to minimal VCF form shows it is a **single het SNP**,
`chr20:49031440 A>G` — exactly what the BAM pipeline calls. The catalog emitted it
three times (from three overlapping snarls) wrapped in 18–49 bp of equal-length
AGGG-repeat context. Because the `VariantKey` derivation in `graph_collect.cpp`
only stripped a single-base prefix for ins/del and never trimmed equal-length
alleles, padded SNPs fell into the MNP fallback with a misleading multi-bp ref/alt
and the wrong POS, and the emitted VCF showed bogus `A→AAGG…AGGC` insertions. See
"Graph allele normalization (padded-SNP bug)" below.

The 3 remaining SNP-rich bad blocks are all in the 32.2–32.5 Mb segdup band and
still require the segdup/centromere-band exclusion described above. Unit tests pass;
site count unchanged (74,119). Validated offline.

## Hybrid Model — minimal-VCF trimming in apply_hybrid_noise_filter (WIN, now default)

**Goal context:** the standing objective is lower Hamming + switch error on the
hybrid pipeline. After the step-4 re-orientation fix (hybrid → 99.18% / 0.82%), the
next lever was a normalization gap between the two noise filters.

**The gap.** `apply_graph_noise_filter` (src/graph_bam_adapter.cpp) trims each
catalog allele to **minimal VCF form** — strip the shared suffix, then the shared
prefix — *before* calling `is_noisy_site`. `apply_hybrid_noise_filter`
(src/hybrid_inject.cpp) did **not**: it ran `is_noisy_site` on the raw catalog
`(ref, alt)`. Catalog (snarl) alleles carry the **full repeat run on both flanks**,
so the derived indel length is inflated. With `max_xgaps` (default 5) bounding the
repeat scan, a genuine 1–2 bp het indel sitting inside a long homopolymer/STR looks
*longer than the window* and escapes detection in one direction but trips it in the
other — the untrimmed hybrid filter **over-demoted** real het indels to
`REP_HET_INDEL` (excluded from k-means), losing them as phasing anchors. The
standalone graph pipeline kept the same sites as `CLEAN_HET_INDEL` because it
trimmed first.

**Fix.** Port the suffix-then-prefix trim into `apply_hybrid_noise_filter`,
adjusting `noisy_pos` by the prefix shift, before the `is_noisy_site` call. Behind a
toggle (`Options::exp_hybrid_trim`), defaulted **on** for hybrid
(collect_hybrid_variation); `--no-hybrid-trim` restores the untrimmed behaviour for
diagnostics.

**Result (chr20, diplinator truth, both toggles from the same HEAD binary):**

| metric | `--no-hybrid-trim` (old) | trim (new default) | Δ |
|---|---|---|---|
| reads evaluated | 212,465 | 212,240 | −225 |
| accuracy | 99.181% | **99.232%** | +0.051 pp |
| Hamming error | 0.819% | **0.768%** | −6.2% rel |
| switch errors | 319 | **291** | −28 (−8.8%) |
| flip errors | 953 | **906** | −47 (−4.9%) |
| switch+flip | 1,272 | **1,197** | −75 (−5.9%) |
| switch rate | 0.150% | **0.137%** | |
| perfect phase sets | 114 | **115** | +1 |
| `CLEAN_HET_INDEL` | 4,652 | **4,994** | +342 |
| `REP_HET_INDEL` | 6,136 | **5,794** | −342 |

**Mechanism.** Exactly 342 graph-only het indels that the untrimmed filter
mislabelled `REP_HET_INDEL` are now correctly `CLEAN_HET_INDEL` — i.e. recovered as
k-means anchors, matching the standalone graph verdict. More correct anchors →
fewer switch/flip events. The slight read-count drop (−225) is the expected effect
of a few sites changing category; net accuracy and all error metrics improve.

This is a *correctness* fix (the two filters now agree on the same site), not a
threshold tweak, so it generalizes rather than overfitting chr20. SNPs are still
exempt (the existing SNP carve-out is unchanged). Unit tests pass
(`test_hybrid_inject`, `test_noise_filter`, all suites). Validated offline.

## Graph allele normalization (padded-SNP bug, FIXED)

**Symptom.** The "large AGGG STR" residual recorded earlier (pos 49,031,428) was
a misdiagnosis. The graph pipeline emitted three records at 49031428 / 49031431 /
49031439 with bizarre equal-length 18–49 bp ref/alt strings, all of which are the
**same single het SNP** `chr20:49031440 A>G` that the BAM pipeline calls cleanly.
The output VCF showed fake insertions like `A → AAGGAAGG…AGGC`.

**Root cause.** The `VariantKey` derivation in `graph_collect.cpp` did not
normalize catalog alleles to minimal VCF form. The ins/del branches only stripped
a single shared **prefix** base (requiring `ref[0]==alt[0]`) and never trimmed a
shared **suffix**; equal-length alleles fell through to the MNP fallback, which
kept the full untrimmed ref/alt and labelled them `SNP` when lengths matched. A
real SNP padded with equal-length repeat context (the AGGG array) therefore got
the wrong POS, a multi-bp ref/alt, and — because the catalog wraps the same site
in several overlapping snarls — was emitted multiple times.

**Scope (chr20).** 5,320 → 1,123 non-minimal alleles after the fix (the remaining
1,123 are genuine complex/MNP variants with differing bases on both ends that have
no single-anchor minimal form). 3,957 of the trimmed cases were single-base SNPs
buried in repeat padding; ~339 catalog sites were duplicate emissions of the same
post-minimization variant.

**Fix.** Normalize each `(pos, ref, alt)` to minimal VCF form (trim shared suffix,
then shared prefix, advance pos) **before** key derivation, via a new shared helper
`trim_to_minimal_vcf` in `noise_filter.{hpp,cpp}`. The same helper now backs the
two pre-existing inline trims in `apply_graph_noise_filter` and
`apply_hybrid_noise_filter` (DRY; behaviour unchanged). `site_meta` (the catalog
form used elsewhere) is left untouched.

**Impact: representation-only, zero phasing change.**

| pipeline | before | after |
|---|---|---|
| graph accuracy / Hamming | 99.195% / 0.805% | **99.195% / 0.805%** (identical) |
| graph switch / flip / perfect-PS | 367 / 858 / 142 | **367 / 858 / 142** (identical) |
| hybrid accuracy / Hamming | 99.232% / 0.768% | **99.232% / 0.768%** (identical) |

Accuracy is unchanged because k-means already classified these equal-length sites
as `CleanHetSnp` (the MNP→SNP branch) and phased them correctly; the bug only
corrupted the **emitted allele strings/POS** in the candidate TSV and VCF. The fix
makes graph/hybrid VCF output match the BAM pipeline's canonical variant
representation. Covered by `test_trim_to_minimal_vcf` in `test_noise_filter.cpp`.
Validated offline.

---

## Duplicate Variant Records & k-means Double-Counting (graph/hybrid)

**Symptom.** After minimal-VCF normalization (see previous section), the graph
pipeline still emitted some variants more than once, and the same physical variant
could enter k-means as several independent anchors.

**Root cause.** The snarl catalog wraps a single physical locus in multiple
overlapping/nested snarls. In `build_graph_chunk` candidates are keyed by snarl
`order_pos()`, so overlapping snarls that decompose to the *same* minimal-VCF
variant look distinct. Two consequences:

1. **Output:** the candidate TSV/VCF carried the same variant several times.
   475 duplicate loci on chr20, 471 sharing a phase set.
2. **k-means:** each duplicate cast an independent germline-het vote. 317
   k-means-eligible duplicates produced 329 redundant anchor votes (0.58% of
   eligible anchors), slightly distorting the clustering.

**Fix A — output dedup (`graph_collect.cpp`).** After the existing
`stable_sort`, a single pass collapses adjacent records with
`exact_comp_cand_var()==0`, keeping the highest-coverage copy and reporting any
conflicting haplotype calls under `--verbose`. Graph VCF duplicates 319→0 within a
chunk (a lone cross-chunk dup remains and is the reason this backstop stays even
with Fix B). Candidate rows 74,119→73,627.

**Fix B — k-means anchor dedup (`graph_bam_adapter.cpp`).** In Phase 2 each
biallelic pair is reduced to its minimal-VCF identity `(pos, ref, alt)` and
deduplicated *across snarl sites* before becoming a k-means anchor:

- Identity is computed with zero allocation: `minimal_vcf_id` trims to offset
  ranges over the existing `meta.ref`/`meta.alts` strings (views, not copies).
- A 64-bit FNV-1a fingerprint (`minimal_vcf_fingerprint`) buckets candidates in an
  `unordered_map<uint64_t, vector<CanonEntry>>`, reserved to the pair upper bound
  to avoid rehashing. Exact `MinimalVcfId` comparison inside the bucket resolves
  fingerprint collisions, so distinct variants are never merged.
- Same-site guard: `CanonEntry` records the source site index; the two alts of a
  single multiallelic snarl are never collapsed into each other (they are distinct
  variants by construction even if their minimal-VCF strings coincide degenerately
  when node sequences are unavailable). Only cross-site matches dedup.
- **No-pool semantics (critical).** A duplicate routes its read observations to the
  canonical anchor but does **not** pool counts. Per-read observation dedup in
  Phase 3 still collapses a read seen at both snarls into one vote.

**Pooled-vs-no-pool A/B (graph pipeline).** Pooling coverage across overlapping
snarls *regressed* phasing: it over-weights the surviving anchor and net-degrades
Hamming. Removing the extra votes *without* distorting the surviving anchor's
weight is what helps.

| variant | accuracy | Hamming | switch | flip |
|---|---|---|---|---|
| NORM-only (baseline) | 99.1954% | 0.8046% | 367 | 858 |
| Fix B pooled | 99.0443% ❌ | 0.9557% ❌ | 355 | 848 |
| **Fix B no-pool (chosen)** | **99.2097%** ✅ | **0.7903%** ✅ | **364** ✅ | **856** ✅ |

Fix B eliminated all 329 within-chunk redundant anchor votes
(eligible-distinct-loci 56,281→56,278, redundant 329→0).

**Hybrid impact: none (verified).** `build_graph_chunk` is shared by graph and
hybrid, so hybrid was re-evaluated before/after Fix A+B on the full chr20 harness:

| pipeline | before | after |
|---|---|---|
| hybrid accuracy / Hamming | 99.232% / 0.768% | **99.232% / 0.768%** (identical) |
| hybrid switch / flip / perfect-PS | 291 / 906 / 115 | **291 / 906 / 115** (identical) |
| hybrid candidate categories | — | identical |

Hybrid is unaffected because its BAM base already anchors phasing; the graph
augment only contributes at sites the BAM pass did not resolve, where the
redundant snarl votes did not change the cluster assignment. The graph-only
pipeline gains (99.195→99.210%) while hybrid holds steady.

Covered by the existing `test_graph_bam_adapter` multiallelic-decomposition tests
(which guard the same-site case) and `make unit-tests`.

---

## DeepVariant calling on HG003 GRCh38: surjection cost and the refined-pbmm2 + graph-HP win

A separate evaluation track on **HG003 chr20** (GIAB v4.2.1, GRCh38) measures how
pgphase's graph-derived phasing affects **downstream DeepVariant variant calling**,
rather than internal phasing accuracy. This is a different sample/reference than the
HG002/CHM13 work above; inputs live under the eval harness (`hg003_vg/`), not in
`test_data/`.

**Setup (all arms identical except where noted):** HG003 Revio SPRQ HiFi, 32×, chr20.
DeepVariant 1.10.0 PACBIO model, `--disable_small_model` in every arm so they share
one calling path (main CNN); the only variable is the alignment substrate and the
phasing source. hap.py vs GIAB HG003 v4.2.1 `noinconsistent` BED (chr20), vcfeval,
PASS-only. Two read alignments of the same reads:

- **pbmm2 (native):** the PacBio case-study minimap2/pbmm2 linear GRCh38 alignment.
- **graph-surject:** vg-giraffe alignment to the HPRC pangenome, surjected to GRCh38.

HP-aware arms feed pgphase HP tags to DV with
`--make_examples_extra_args "phase_reads=false,sort_by_haplotypes=true"`.

### Vanilla control reproduces the published case study

DV on the native pbmm2 BAM (DV-internal phasing) matches Google's published number
to within 0.0001 INDEL F1 (our 0.99359 vs published 0.99368; identical TP 10,561),
validating the harness. The graph-surject pipeline scores **lower** on raw DV calling:

| Config (best arm) | Alignment | INDEL F1 | SNP F1 |
|---|---|---|---|
| Vanilla pbmm2 (DV-internal) | pbmm2 | 0.99359 | 0.99910 |
| graph-surject hybrid_hp (mq20 af0.12) | surject | 0.98724 | 0.99877 |

The ~0.0064 INDEL-F1 gap is the **surjection penalty**, not the phasing. SNPs are
nearly immune (~0.0003), because the cost is almost entirely indel re-representation:
projecting a graph path onto linear GRCh38 makes CIGAR/left-shift decisions that
differ from minimap2's, and the PACBIO model was trained on minimap2 alignments.

### Threshold sweeps are saturated

`--min-af` (0.20→0.12) and `--min-mapq` (20→10→1) were swept on the surject pipeline.
AF 0.20→0.12 gained +0.00145 INDEL F1 (the only real mover) by recruiting low-VAF
candidates; MAPQ moved INDEL F1 by ≤0.0004 across the whole 20→1 range (plateaued).
Lowering either threshold recruits **mostly SNPs, not indels** (e.g. mq20→mq1 added
+5,887 SNP vs +827 indel candidates). You cannot tune past the surjection penalty
with candidate thresholds. Best surject config: **mq20 af0.12 / dv_hybrid_hp**.

### The fix: native pbmm2 substrate + graph phasing

Putting graph phasing on a **native** alignment closes the entire gap. Three
constructions, vs the vanilla baseline (TP/FP are PASS counts):

| Config | INDEL F1 | INDEL FP | SNP F1 | SNP FP |
|---|---|---|---|---|
| Vanilla pbmm2 (DV-internal) | 0.99359 | 72 | 0.99910 | 68 |
| pbmm2-hybrid, dv_hybrid_hp | **0.99420** | 67 | 0.99789 | **241** ❌ |
| **refined-pbmm2 + graph-surject HP** | 0.99396 | 69 | **0.99917** | **59** ✅ |
| graph-surject hybrid_hp | 0.98724 | — | 0.99877 | — |

- **pbmm2-hybrid** = run `collect-hybrid-variation` with the native pbmm2 BAM as
  `--bam` plus the giraffe GAF/sites (`--refine-aln --min-mapq 20 --min-af 0.12`).
  This is coordinate-safe: the GAF is queried by GRCh38 interval (its own surjection
  projection) and graph observations join to BAM reads by **read name only**, so a
  different BAM aligner still merges. It gives the **best INDEL F1 (0.99420, +8 TP)**
  but regresses SNP precision (FP 68→241).

- **refined-pbmm2 + graph-surject HP** = take the BAM-pipeline refined pbmm2 reads
  (`run_pgphase_pbmm2hybrid_af12/bam/phased.bam` — fixed indel CIGARs, native
  coords), strip their HP/PS, and stamp the **graph-surject** hybrid HP tags
  (`run_pgphase_refine_mq20_af12/hybrid/phased.bam`) by read name (`gapfill.py` with
  core=refined-pbmm2, fill=graph-surject). This is a **strict Pareto win over
  vanilla** — every metric improves, both types:

| | INDEL | SNP |
|---|---|---|
| TP | 10,561 → 10,566 (+5) | 70,107 → 70,108 (+1) |
| FN | 67 → 62 (−5) | 59 → 58 (−1) |
| FP | 72 → 69 (−3) | 68 → 59 (−9) |
| Recall | 0.99370 → 0.99417 | 0.99916 → 0.99917 |
| Precision | 0.99348 → 0.99376 | 0.99903 → 0.99916 |
| F1 | 0.99359 → 0.99396 | 0.99910 → 0.99917 |

### What caused the SNP regression (and what did not)

The +173 SNP FP in the pbmm2-hybrid arm is **not** from `--refine-aln`. Two arms on
the *same refined substrate* isolate it: pbmm2-hybrid HP (graph fed in pbmm2 coords)
→ SNP FP 241; graph-surject HP (graph phased in graph coords, then surjected) → SNP
FP 59. The damage comes entirely from the **HP source**, not from refine. The
in-pbmm2-coordinate hybrid phasing mislabels reads at some SNP sites (the untested
no-conflict profile-extension path, `hybrid_inject.cpp:295`, where pbmm2 and giraffe
disagree on placement); graph-native phasing does not. `dv_default` on the refined
pbmm2 BAM keeps SNP FP at 65, confirming refine itself is clean.

### Recommendations

- **Graph phasing is not the problem; the surjected alignment is.** Graph-derived HP
  labels are accurate enough to match (and slightly beat) DV's own phasing — when
  applied to a native alignment.
- **For a no-downside config:** refined-pbmm2 substrate + graph-surject hybrid HP.
  Beats vanilla DV on both indels and SNPs with no regression.
- **For max indel recall (accepting SNP cost):** pbmm2-hybrid dv_hybrid_hp.
- **Use the hybrid (clean-core) HP source, not gapfill** — gapfill's hard reads cost
  precision on the native substrate, same as on surject.

Eval scripts: `hg003_vg/{26..43}_*.sbatch` (pgphase / DeepVariant / hap.py per
config) and `bench_hybrid_refine_mq{1,10,20}_af12.sh`. Results under
`hg003_vg/dv_results_*/`. These artifacts live on the eval cluster, not in-repo.
