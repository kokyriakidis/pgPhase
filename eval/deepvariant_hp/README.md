# DeepVariant HP-tag evaluation

Does pgphase read phasing improve DeepVariant small-variant calling?

DeepVariant's PacBio model phases reads internally before calling. This harness
substitutes that internal phasing with pgphase's HP tags and measures the effect
on variant-calling accuracy (hap.py vs GIAB truth).

## The two arms

Both arms call variants from the **same reads/alignment** with DeepVariant
`PACBIO` (r1.6.1). Only the source of read phasing differs.

| Arm | Input BAM | DeepVariant phasing | Isolates |
|---|---|---|---|
| `dv_default` | raw BAM (no HP tags) | internal (`phase_reads=true`, default) | stock DeepVariant baseline |
| `dv_our_hp` | pgphase phased BAM (HP/PS tags) | **off** (`phase_reads=false`), uses our HP | pgphase phasing vs DeepVariant's own |

`dv_our_hp` keeps `sort_by_haplotypes=true` and the PacBio-default
`add_hp_channel`/`parse_sam_aux_fields`, so DeepVariant reads the HP tags already
present in the pgphase BAM and feeds them into the pileup HP channel — instead of
computing its own phasing. If `dv_our_hp` F1 > `dv_default` F1, pgphase phasing
helps.

The relevant PacBio `make_examples` defaults are taken from DeepVariant r1.6
`run_deepvariant.py`: `add_hp_channel=true`, `sort_by_haplotypes=true`,
`phase_reads=true`, `parse_sam_aux_fields=true`, `realign_reads=false`.

## Why this cannot run in the Ona dev container

DeepVariant and hap.py ship only as containers, and the dev container has no
container runtime (no docker daemon, singularity, or apptainer) and no GPU.
Run this on the cluster. Everything needed is downloaded by `fetch_resources.sh`.

## Prerequisites (cluster)

- A container engine: `singularity`/`apptainer` (recommended on HPC) or `docker`
- `curl`, `samtools`, `python3` (stdlib only)
- ~30 GB disk for reference + images + intermediates
- The two input BAMs (see below)

## Inputs you provide

Both must cover the same reads aligned to **GRCh38** (the truth set is GRCh38):

- `--raw-bam` — aligned, indexed BAM **without** pgphase HP tags (arm A).
- `--our-bam` — pgphase phased BAM **with** HP/PS tags (arm B), e.g. produced by:

  ```bash
  pgphase collect-hybrid-variation \
    --ref GRCh38.fa --bam raw.bam \
    --graph-sites sites.vcf.gz --gaf reads.coord.gaf.gz \
    --hifi -o cand.tsv -b our_phased.bam      # add --gap-fill to test gap-fill
  samtools index our_phased.bam
  ```

> Contig names must match GRCh38 (`chr20`, …). The bundled `test_data/graph_chr20`
> slice is N-padded synthetic coordinates with no truth set and **cannot** be used
> here — use a real GRCh38-aligned HG002 BAM.

## Run

```bash
cd eval/deepvariant_hp

# 1. One-time: download reference, GIAB HG002 v4.2.1 truth, and pull images.
./fetch_resources.sh                 # add --engine singularity|docker to force

# 2. Run both DeepVariant arms + hap.py + the comparison.
./run_eval.sh \
    --raw-bam /path/to/raw.bam \
    --our-bam /path/to/our_phased.bam \
    --outdir  results \
    --regions chr20 \
    --threads "$(nproc)"
```

Re-running skips arms whose VCF / hap.py summary already exist, so an
interrupted run resumes cleanly.

## Output

```
results/
  dv_default/deepvariant.vcf.gz   dv_default/happy.*    # baseline
  dv_our_hp/deepvariant.vcf.gz    dv_our_hp/happy.*     # pgphase HP
  comparison.txt                                        # side-by-side
```

`comparison.txt` (also printed) tabulates SNP and INDEL recall / precision / F1
for both arms and the delta, e.g.:

```
Type   Metric     DV default (self-phased)   DV + pgphase HP        Delta
SNP    F1                         0.999658          0.999758    +0.000100
INDEL  F1                         0.993259          0.994101    +0.000842

SNP: pgphase HP improves DeepVariant F1 (+0.000100)
INDEL: pgphase HP improves DeepVariant F1 (+0.000842)
```

(numbers above are illustrative)

## Files

| File | Purpose |
|---|---|
| `fetch_resources.sh` | Download reference, GIAB truth, DV + hap.py images |
| `run_eval.sh` | Run both DV arms, benchmark with hap.py, print comparison |
| `summarize_happy.py` | Parse two `happy.summary.csv` and tabulate the A/B delta |

## Notes

- **Versions are pinned** (DeepVariant 1.6.1, hap.py v0.3.12, GIAB v4.2.1) for
  reproducibility. Change them at the top of the scripts if needed.
- **hap.py uses the `vcfeval` engine** and `--pass-only`, matching the
  DeepVariant PacBio case study, so numbers are comparable to published results.
- For **ONT** data, change `MODEL_TYPE` to `ONT_R104` in `run_eval.sh`; the
  arm-B override (`phase_reads=false,sort_by_haplotypes=true`) is the same.
- `--regions` defaults to `chr20`, the standard fast benchmark. Pass a different
  contig (or a comma-separated list) to evaluate elsewhere; the GIAB v4.2.1 truth
  covers chr1–22.
