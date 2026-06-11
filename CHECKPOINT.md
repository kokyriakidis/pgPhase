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

## Lessons / Pitfalls

- **pgbam sidecar causes over-stitching:** Running BAM pipeline with `--pgbam-file` merged all 49 initial phase sets into 2 chr20-spanning blocks with near-random accuracy (55%). Do not use the sidecar without understanding its signal quality.
- **Wrong reference = 0 reads:** `hg002v1.1.fasta` (maternal/paternal contigs) does not match the `CHM13#0#chr20` contig name in the aligned BAM. Must use `chm13v2.0.chr20.renamed.fa` as ref for `collect-bam-variation`.
- **Use chr20-only truth refs:** Full-genome maternal/paternal assemblies cause a small fraction of reads to leak to chr1/6/16 during truth alignment. Use `HG002_1.1_MATERNAL.chr20.fasta` / `HG002_1.1_PATERNAL.chr20.fasta` for clean results.
- **Per-chunk site loading:** Graph pipeline uses tabix-indexed per-chunk site queries (commit `5878aad`) — loading all 977K sites upfront was the prior bottleneck.
- **Non-numeric node IDs:** `build_compact_graph_site_index` now skips sites with non-numeric node IDs instead of aborting (commit `5c01954`).

---

## BAM Pipeline Parity: `propagate_overlap_read_phase_to_output_owner`

### Observation

Side-by-side evaluation of pgphase BAM pipeline vs longcallD on the same input showed pgphase
phasing ~1,300 more reads (220,377 vs 219,090) and producing 2 extra phase sets (431 vs 429).
All other metrics (N50, auN, largest block, accuracy, switch/hamming error) were near-identical.

### Root cause

`propagate_overlap_read_phase_to_output_owner` (commit `f094cea`) is a pgphase-only addition
in `stitch_chunk_haps`. After chunk-boundary stitching, it copies HP/PS from the downstream
chunk to overlap reads that are unphased (hap=0) in the upstream (owning) chunk. LongcallD
does not have this function.

### Why longcallD omits it

Overlap reads exist in both adjacent chunks as independent copies. Each chunk phases its copy
using its own candidates and phasing context. The BAM writer emits each read from its owning
(upstream) chunk. `apply_chunk_flip_and_merge` only rewrites the **downstream** chunk's PS
values — it never touches the upstream chunk's reads.

When `flip_chunk_hap` merges two chunks (flip_hap_score != 0), the downstream copy's PS is
rewritten to `max_pre_ps` (the upstream PS) and hap labels are flipped if needed. In this case,
propagation would be correct: the downstream copy's HP/PS is now consistent with the upstream
chunk's phasing. The upstream copy just happens to be unphased because its variant profile in
the owning chunk had no informative candidates.

However, when `flip_hap_score == 0` (equal agree/disagree votes), the chunks are **not merged**
and their relative phase is unknown. `propagate_overlap_read_phase_to_output_owner` does not
distinguish between merged and unmerged pairs — it propagates unconditionally. This is
incorrect for unmerged pairs:

1. **Wrong haplotype.** The downstream hap=1 may correspond to the upstream hap=2. Propagating
   without adjusting for the unknown relative phase can assign the read to the wrong haplotype.

2. **Orphan phase sets.** The propagated PS value comes from the downstream chunk's unmerged
   phase block. This PS doesn't correspond to any candidate in the upstream chunk, creating
   phase sets that inflate the count without improving phasing.

### Decision

Conditional propagation: `stitch_chunk_haps` now tracks which adjacent pairs were successfully
merged by `flip_chunk_hap` (via a `pair_stitched` vector) and only calls
`propagate_overlap_read_phase_to_output_owner` for those pairs. For merged pairs the downstream
PS has been rewritten to `max_pre_ps` and hap labels flipped if needed, so propagation is safe.
Unmerged pairs (flip_hap_score == 0) are skipped — no propagation, no wrong-haplotype risk.

This recovers the ~1,300 overlap reads that longcallD leaves unphased at chunk boundaries,
without introducing orphan phase sets or haplotype errors. Expected effect vs longcallD:
slightly more phased reads, same phase set count, same or better accuracy.
