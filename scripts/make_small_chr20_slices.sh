#!/usr/bin/env bash
# Build tiny BAMs from whole-chr20 alignments for fast pgphase / longcallD loops.
#
# Default window: chr20 15.0–15.5 Mb (500 kb). Override with START / END (1-based inclusive).
#
# Outputs under test_data/chr20_quick/ (override OUTDIR):
#   HG003_GRCh38_chr20_<START>_<END>.bam(.bai)  — from deepVariantTest HG003 chr20 BAM
#   HG002_CHM13_chr20_<START>_<END>.bam(.bai)   — from CHM13 HG002 chr20 BAM if present
#   GRCh38_chr20_only.fa(.fai)                  — single-contig chr20 from GRCh38 (pairs with HG003 slice)
#
# Pair for runs (coordinates stay on full chr20):
#   FA="$ROOT/test_data/chr20_quick/GRCh38_chr20_only.fa"
#   BAM="$ROOT/test_data/chr20_quick/HG003_GRCh38_chr20_<START>_<END>.bam"
#   REGION="chr20:<START>-<END>"
#
# CHM13 HG002 (reference stays external unless you copy chr20 FA yourself):
#   BAM="$ROOT/test_data/chr20_quick/HG002_CHM13_chr20_<START>_<END>.bam"
#   FA="<path>/chm13v2.0.chr20.renamed.fa"
#   REGION="CHM13#0#chr20:<START>-<END>"
#
# Env:
#   SAMTOOLS, OUTDIR, START, END
#   CHM13_CHR20_BAM — default /home/kokyriakidis/Downloads/pgbam-experiments/HG002.chr20.bam
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${OUTDIR:-$ROOT/test_data/chr20_quick}"
SAMTOOLS="${SAMTOOLS:-samtools}"
START="${START:-15000001}"
END="${END:-15500000}"
CHM13_CHR20_BAM="${CHM13_CHR20_BAM:-/home/kokyriakidis/Downloads/pgbam-experiments/HG002.chr20.bam}"

fail() { echo "error: $*" >&2; exit 2; }
command -v "$SAMTOOLS" >/dev/null || fail "samtools not found"

mkdir -p "$OUTDIR"

slice_bam() {
  local src="$1" region="$2" dst="$3"
  [[ -f "$src" ]] || fail "missing BAM: $src"
  "$SAMTOOLS" view -b "$src" "$region" >"$dst"
  "$SAMTOOLS" index "$dst"
}

HG003_BAM="$ROOT/test_data/deepVariantTest/HG003.novaseq.pcr-free.35x.vg-1.55.0.chr20.bam"
GRCH38_FA="$ROOT/test_data/deepVariantTest/GRCh38_no_alt_analysis_set.fasta"

suffix="${START}_${END}"
GRCH38_REGION="chr20:${START}-${END}"
CHM13_REGION="CHM13#0#chr20:${START}-${END}"

if [[ -f "$HG003_BAM" ]]; then
  echo "slicing HG003 GRCh38 BAM -> $OUTDIR/HG003_GRCh38_chr20_${suffix}.bam ($GRCH38_REGION)"
  slice_bam "$HG003_BAM" "$GRCH38_REGION" "$OUTDIR/HG003_GRCh38_chr20_${suffix}.bam"
else
  echo "skip: HG003 BAM not found at $HG003_BAM"
fi

if [[ -f "$CHM13_CHR20_BAM" ]]; then
  echo "slicing HG002 CHM13 BAM -> $OUTDIR/HG002_CHM13_chr20_${suffix}.bam ($CHM13_REGION)"
  slice_bam "$CHM13_CHR20_BAM" "$CHM13_REGION" "$OUTDIR/HG002_CHM13_chr20_${suffix}.bam"
else
  echo "skip: CHM13 chr20 BAM not found at $CHM13_CHR20_BAM (set CHM13_CHR20_BAM)"
fi

if [[ -f "$GRCH38_FA" ]]; then
  echo "writing chr20-only GRCh38 FASTA -> $OUTDIR/GRCh38_chr20_only.fa"
  "$SAMTOOLS" faidx "$GRCH38_FA" "chr20" >"$OUTDIR/GRCh38_chr20_only.fa"
  "$SAMTOOLS" faidx "$OUTDIR/GRCh38_chr20_only.fa"
else
  echo "skip: GRCh38 FASTA not found at $GRCH38_FA"
fi

echo "done. Listing $OUTDIR:"
ls -lh "$OUTDIR"
