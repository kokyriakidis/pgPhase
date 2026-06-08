# pgphase

Variant calling and haplotype phasing for long reads, with two pipelines:

- **`collect-graph-variation`** — call and phase variants from a pangenome graph (GBZ + GAF).
- **`collect-bam-variation`** — call and phase SNPs/indels from BAM/CRAM alignments (HiFi, ONT, short reads).
- **`build-snarl-catalog`** — preprocess a GBZ pangenome graph into a phasing site catalog VCF for `collect-graph-variation`.

## Prerequisites

- C++17 compiler (GCC ≥ 8 or Clang ≥ 7)
- [htslib](https://github.com/samtools/htslib) (development headers and library)
- zlib, pthreads
- Rust toolchain (for the `gbz-base` graph query tools; install via [rustup](https://rustup.rs/))
- [vg](https://github.com/vgteam/vg) (only needed for `build-snarl-catalog`)
- [samtools](https://github.com/samtools/samtools) (optional; needed for `--refine-aln` coordinate sorting)

## Building

```bash
git clone --recursive https://github.com/kokyriakidis/pgphase.git
cd pgphase

# Build third-party libraries (WFA2, abPOA)
make third-party-libs

# Build the gbz-base Rust tools (query, gaf2db, gbz2db)
make gbz-base

# Build pgphase
make -j$(nproc)
```

To build a portable self-contained bundle (binary + shared libraries + launcher):

```bash
make portable-bundle
# Output: dist/pgphase-<os>-<arch>/
```

## Quick start

### Graph pipeline

1. Build a site catalog from a pangenome GBZ:

```bash
# Precompute snarls once
vg snarls -t 16 full.gbz > full.snarls.pb

# Build per-chromosome catalog VCFs
for chr in chr{1..22} chrX chrY; do
    vg chunk --gbz --contig $chr -x full.gbz -o /tmp/chunk_${chr}
    pgphase build-snarl-catalog \
        --ref-sample CHM13 \
        --contig $chr \
        --snarls full.snarls.pb \
        --threads 8 \
        -o ${chr}.sites.vcf.gz \
        /tmp/chunk_${chr}_graph_0_${chr}.gbz &
done
wait
```

2. Run graph-based phasing:

```bash
pgphase collect-graph-variation \
    --ref ref.fa \
    --sites chr20.sites.vcf.gz \
    --gaf reads.gaf \
    --phased-vcf-out phased.vcf \
    --phased-bam-out phased.bam \
    -t 8
```

The sites VCF must be bgzip-compressed with a tabix index (`.tbi` or `.csi`).

### BAM pipeline (HiFi)

```bash
pgphase collect-bam-variation \
    --ref ref.fa \
    --bam hifi.bam \
    --hifi \
    --phased-vcf-out phased.vcf \
    -o candidates.tsv \
    -t 8
```

### BAM pipeline (ONT)

```bash
pgphase collect-bam-variation \
    --ref ref.fa \
    --bam ont.bam \
    --ont \
    --phased-vcf-out phased.vcf \
    -o candidates.tsv \
    -t 16
```

## Subcommands

### `collect-graph-variation`

Collects and phases variants using a pangenome graph site catalog and GAF read alignments.

| Option | Description | Default |
|---|---|---|
| `--ref FILE` | Reference FASTA (indexed) | required |
| `--sites FILE` | Sites VCF (bgzipped + tabix-indexed) | required |
| `--gaf FILE` | Coordinate-indexed GAF from [pggaf](https://github.com/kokyriakidis/pggaf) (bgzipped + tabix) | — |
| `--gbz-db FILE` | GBZ graph database (from `gbz2db`) | — |
| `--gaf-db FILE` | GAF-base read alignment database (from `gaf2db`) | — |
| `--sample NAME` | Reference sample name for GBZ queries | auto-detected |
| `--hifi` | HiFi read mode | default |
| `--ont` | ONT read mode (enables strand-bias filter) | off |
| `--strand-bias-pval FLOAT` | Max p-value for ONT strand-bias filter | 0.01 |
| `-q INT` | Minimum mapping quality | 30 |
| `--phased-vcf-out FILE` | Phased VCF | — |
| `--phased-bam-out FILE` | Unaligned BAM with HP/PS tags | — |
| `-r STR` | Restrict to region (repeatable) | whole genome |
| `-t INT` | Worker threads | 1 |

Provide either `--gaf` (coordinate-indexed GAF from [pggaf](https://github.com/kokyriakidis/pggaf), supports `--pgbam-file` stitching) or `--gbz-db` + `--gaf-db` (pre-built databases from gbz-base).

### `collect-bam-variation`

Collects SNP/indel candidates from BAM/CRAM alignments and phases them using k-means clustering on read haplotype support.

Key options:

| Option | Description | Default |
|---|---|---|
| `--ref FILE` | Reference FASTA (indexed) | required |
| `--bam FILE` | Input BAM/CRAM | required |
| `--hifi` / `--ont` / `--short-reads` | Read technology mode | `--hifi` |
| `-t INT` | Worker threads | 1 |
| `-o FILE` | Candidate TSV output | `output.tsv` |
| `--phased-vcf-out FILE` | Phased VCF (GT:DP:AD:VAF:GQ:PS) | — |
| `--out-bam FILE` | Phased BAM with HP/PS tags | — |
| `-r STR` | Restrict to region (repeatable) | whole genome |
| `--autosome` | Process chr1–22 only | off |
| `-q INT` | Minimum mapping quality | 30 |
| `-D INT` | Minimum depth | 5 |
| `--min-af FLOAT` | Minimum allele fraction | 0.20 |
| `--pgbam-file FILE` | `.pgbam` sidecar for chunk stitching | — |

Run `pgphase collect-bam-variation --help` for the full option list.

### `build-snarl-catalog`

Runs `vg deconstruct -a` on a GBZ graph and produces a bgzipped, tabix-indexed sites-only VCF for use with `collect-graph-variation`.

| Option | Description | Default |
|---|---|---|
| `--ref-sample STR` | Reference sample name (e.g. `CHM13`) | required |
| `--contig STR` | Process one contig (recommended) | all |
| `-o FILE` | Output bgzipped VCF | required |
| `--snarls FILE` | Pre-built snarls file | recomputed |
| `--vg-bin PATH` | Path to `vg` binary | `vg` |
| `-t INT` | Threads for `vg deconstruct` | 4 |

## Testing

```bash
# Run validation gates (requires test_data/)
make check

# Run unit tests
make unit-tests
```

## Release

```bash
# Build + portable bundle + tarball
make release

# Strict mode: fails if output differs from golden fixtures
make release-strict
```

Output: `dist/pgphase-<os>-<arch>.tar.gz`

## License

See repository for license details.
