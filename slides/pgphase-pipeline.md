---
marp: true
theme: roche
paginate: true
---

<!-- _class: title -->

# pgphase

## Graph-Based Variant Calling & Haplotype Phasing from Pangenome Long Reads

collect-graph-variation pipeline · Technical Overview

---

<!-- _class: content -->

# Pipeline Inputs

Three required inputs for the graph-based pipeline

| | **GBZ Graph** | **GAF Alignments** | **Sites VCF** |
|---|---|---|---|
| **What** | Pangenome variation graph with GBWT haplotype index | Long reads aligned to the pangenome graph | Snarl decomposition of graph variant sites |
| **Structure** | Nodes = sequence segments, Edges = adjacencies, Threads = haplotype paths | Read paths through graph nodes, HiFi / ONT long reads | Each snarl = variant site, multi-allelic genotypes |
| **Format** | GBZ (vg toolkit) | Graph Alignment Format | VCF with reference-relative coords |

> **`pgphase collect-graph-variation`** — Queries graph + alignments per genomic chunk

---

<!-- _class: content -->

# Step 1: Chunked Graph Traversal

Process the genome in overlapping chunks for parallelism

| Chunk 1 | Chunk 2 | Chunk 3 | Chunk 4 | Chunk 5 | Chunk 6 |
|:---:|:---:|:---:|:---:|:---:|:---:|
| chr20:0-5M | chr20:4-9M | chr20:8-13M | chr20:12-17M | chr20:16-21M | chr20:20-25M |

*Overlapping regions enable phase block stitching across chunk boundaries*

**① Query** graph nodes + read alignments in chunk — *gbz-base SQLite queries via Rust FFI*

**② Collect** variant candidates from snarl traversals — *multi-allelic → biallelic candidates*

**③ Classify** regions: clean SNPs vs noisy indels — *allele fraction filter (AF 0.20–0.80)*

**④ Phase** per-chunk phasing and variant calling — *k-means on clean het SNPs, MSA on noisy*

---

<!-- _class: two-column -->

# Step 2: Variant Discovery from Snarls

**Snarl = Bubble in the Graph**

A snarl is a subgraph between two boundary nodes (S → E) with multiple traversal paths:

- **REF path**: ACGT (reference allele)
- **ALT1 path**: ACGGT (insertion)
- **ALT2 path**: AT (deletion)

Multi-allelic sites are converted to biallelic candidates. Read traversals are counted per allele path. Allele fraction filter: 0.20 < AF < 0.80.

**Region Classification**

**Clean Regions**
- Germline het SNPs (AF 0.20–0.80)
- High-confidence single-base variants
- Used as input for k-means phasing

**Noisy Regions**
- Complex indels, MNVs, structural variants
- Multiple overlapping variants
- Resolved via haplotype-aware MSA (abPOA)
- WFA2 / edlib realignment

---

<!-- _class: content -->

# Step 3: Phasing & Stitching

K-means clustering on het SNPs, then stitch chunks into chromosome-wide phase blocks

| | **K-means Phasing (k=2)** | **Phase Block Stitching** |
|---|---|---|
| **Input** | Reads × germline het SNP genotypes | Overlapping reads between adjacent chunks |
| **Method** | Each read = binary vector (REF=0, ALT=1); 2-means clustering assigns reads to HP1 or HP2 | Vote on haplotype consistency across boundary; flip HP1/HP2 if majority disagrees |
| **Scope** | Per-chunk: independent phasing | Cross-chunk: chromosome-wide phase blocks (PS tag) |
| **Notes** | Clean SNPs only (noisy sites excluded) | Optional: pgbam sidecar thread concordance |

### Noisy Region Consensus (Haplotype-Aware)

1. Partition reads by HP assignment from k-means phasing
2. Per-haplotype multiple sequence alignment using **abPOA**
3. Generate consensus sequence per haplotype
4. Realign reads to consensus with **WFA2** (or edlib fallback)
5. Call variants from haplotype-specific realignment

---

<!-- _class: content -->

# Step 4: Output

Final outputs from the graph-based pipeline

| Output | Description | Details |
|---|---|---|
| **Phased uBAM** | HP:i: and PS:i: tags per read | Haplotype assignment for every aligned read |
| **VCF** | Germline variant calls | SNPs and indels with GT, AF, DP fields |
| **Phase Blocks** | Chromosome-wide phase sets | Stitched from per-chunk k-means results |

### Pipeline Summary

> GBZ + GAF + VCF → **chunk** → **snarl variants** → **k-means phase** → **stitch** → **consensus** → VCF + uBAM

---

<!-- _class: content -->

# Full Pipeline Flow

| Stage | Step | Tool / Method | Output |
|---|---|---|---|
| **Input** | Load GBZ + GAF + Sites VCF | vg / gbz-base | Graph + alignments + sites |
| **Chunking** | Divide genome into overlapping windows | collect-graph-variation | Per-chunk data |
| **Graph Query** | Query nodes, edges, threads per chunk | gbz-base FFI (Rust → SQLite) | Local subgraph |
| **Variants** | Traverse snarls for variant candidates | Snarl decomposition | Biallelic candidates |
| **Classify** | Clean SNPs vs noisy regions | AF filter (0.20–0.80) | Region labels |
| **Phase** | K-means clustering (k=2) | 2-means on het SNP vectors | HP1 / HP2 per read |
| **Stitch** | Merge phase blocks across chunks | Overlap voting + pgbam | Chromosome-wide PS |
| **Consensus** | Haplotype-aware MSA for noisy regions | abPOA + WFA2 / edlib | Refined variants |
| **Output** | Write phased BAM + VCF | collect-graph-variation | Final deliverables |

---

<!-- _class: content -->

# Evaluation: Overall Phase Accuracy

HG002 T2T diploid assembly · chr20 · HiFi 30× · diplinator truth set

| Metric | pgphase (graph) | longcallD (linear) | WhatsHap |
|---|:---:|:---:|:---:|
| **Overall accuracy** | **XX.XX%** | XX.XX% | XX.XX% |
| **Hamming error rate** | **X.XX%** | X.XX% | X.XX% |
| **Switch error rate** | **X.XX%** | X.XX% | X.XX% |
| **Phase block N50** | **X.X Mb** | X.X Mb | X.X Mb |
| **Perfect phase sets** | XX / XX (XX%) | XX / XX (XX%) | XX / XX (XX%) |
| **Reads evaluated** | XX,XXX | XX,XXX | XX,XXX |

> *Replace XX values with actual benchmark results from `summary.json`*

---

<!-- _class: content -->

# Evaluation: Per-Category Accuracy

Accuracy breakdown by genomic region (censat, segdup, unique)

| Category | Reads | Accuracy | Hamming Error | Switch Error |
|---|:---:|:---:|:---:|:---:|
| **unique** | XX,XXX | XX.XX% | X.XX% | X.XX% |
| **segdup** | X,XXX | XX.XX% | X.XX% | X.XX% |
| **censat** | X,XXX | XX.XX% | X.XX% | X.XX% |
| **censat+segdup** | XXX | XX.XX% | X.XX% | X.XX% |

### Phaseable vs Unphaseable

| | Phase Sets | Reads | Accuracy |
|---|:---:|:---:|:---:|
| **Phaseable (>60%)** | XX | XX,XXX | XX.XX% |
| **Unphaseable (≤60%)** | XX | X,XXX | XX.XX% |

> *Phaseable threshold: phase sets with >60% accuracy. Unphaseable blocks often lack het sites in the graph.*

---

<!-- _class: content -->

# Evaluation: Graph vs Linear Pipeline

Key advantages of the graph-based approach

| Aspect | Graph (pgphase) | Linear (longcallD) |
|---|---|---|
| **Reference bias** | Eliminated — reads traverse haplotype paths | Present — mismatches at variant sites |
| **Minimizer seeding** | Not affected by known variants | Windows shift at every mismatch |
| **Variant discovery** | Snarl topology defines sites | Pileup-based candidate detection |
| **Complex regions** | Graph encodes known structural variation | Relies on realignment heuristics |
| **Phasing** | Same k-means approach | Same k-means approach |
| **Noisy regions** | abPOA + WFA2 consensus | abPOA + WFA2 consensus |

### Expected Improvements

- Higher accuracy in **segdup** and **censat** regions where reference bias is strongest
- Better variant discovery at **multi-allelic sites** encoded as graph snarls
- Reduced false negatives from **minimizer window shifts** at variant-dense regions

---

<!-- _class: title -->

# Thank You

## Questions?

pgphase · collect-graph-variation pipeline
