#!/usr/bin/env bash
#
# bench_hiphase.sh — Variant-call with DeepVariant, phase with HiPhase,
#                    and produce a haplotagged BAM.
#
# Pipeline:
#   1. DeepVariant — call variants from aligned BAM
#   2. hiphase     — phase + haplotag in a single pass
#
# Environment setup (run once):
#   micromamba create -f envs/hiphase.yaml -y
#   micromamba activate bench-hiphase
#
# Usage:
#   ./scripts/bench_hiphase.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --outdir results/hiphase \
#       --threads 16

set -euo pipefail

REF=""
BAM=""
VCF=""
SV_VCF=""
OUTDIR=""
THREADS=16
DV_MODEL="PACBIO"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --ref FILE       Reference FASTA (indexed)
  --bam FILE       Aligned BAM/CRAM (indexed)
  --outdir DIR     Output directory

Optional:
  --vcf FILE       Pre-called VCF (skip DeepVariant step)
  --sv-vcf FILE    SV VCF for joint phasing (e.g. from pbsv)
  --threads INT    Threads [16]
  --dv-model STR   DeepVariant model type [PACBIO]
  -h, --help       Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)        REF="$2";        shift 2 ;;
        --bam)        BAM="$2";        shift 2 ;;
        --vcf)        VCF="$2";        shift 2 ;;
        --sv-vcf)     SV_VCF="$2";     shift 2 ;;
        --outdir)     OUTDIR="$2";     shift 2 ;;
        --threads)    THREADS="$2";    shift 2 ;;
        --dv-model)   DV_MODEL="$2";   shift 2 ;;
        -h|--help)    usage ;;
        *)            echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "$REF" ]]    && { echo "Error: --ref required" >&2; usage; }
[[ -z "$BAM" ]]    && { echo "Error: --bam required" >&2; usage; }
[[ -z "$OUTDIR" ]] && { echo "Error: --outdir required" >&2; usage; }

for f in "$REF" "$BAM"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done
for cmd in hiphase samtools; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

DV_VCF="$OUTDIR/deepvariant.vcf.gz"
PHASED_VCF="$OUTDIR/phased.vcf.gz"
HAPLOTAGGED_BAM="$OUTDIR/haplotagged.bam"
STATS_TSV="$OUTDIR/stats.tsv"
BLOCKS_TSV="$OUTDIR/blocks.tsv"
SUMMARY_TSV="$OUTDIR/summary.tsv"
TIMING_LOG="$OUTDIR/timing.log"

# ── Step 1: Call variants with DeepVariant ────────────────────────────
if [[ -n "$VCF" ]]; then
    echo "[1/2] Using provided VCF: $VCF"
    DV_VCF="$VCF"
elif [[ ! -f "$DV_VCF" ]]; then
    echo "[1/2] Calling variants with DeepVariant (model: ${DV_MODEL}) ..."
    command -v run_deepvariant >/dev/null 2>&1 || \
        { echo "Error: run_deepvariant not found. Activate bench-hiphase env." >&2; exit 1; }
    { time run_deepvariant \
        --model_type="${DV_MODEL}" \
        --ref="$REF" \
        --reads="$BAM" \
        --output_vcf="$DV_VCF" \
        --num_shards="$THREADS" ; } 2> >(tee "$TIMING_LOG" >&2)
else
    echo "[1/2] DeepVariant VCF exists, skipping."
fi

# ── Step 2: Phase + haplotag with HiPhase ─────────────────────────────
SV_ARGS=""
SV_OUT_ARGS=""
if [[ -n "$SV_VCF" ]]; then
    [[ -f "$SV_VCF" ]] || { echo "Error: SV VCF not found: $SV_VCF" >&2; exit 1; }
    SV_ARGS="--vcf $SV_VCF"
    SV_OUT_ARGS="--output-vcf $OUTDIR/phased_sv.vcf.gz"
fi

if [[ ! -f "$HAPLOTAGGED_BAM" ]]; then
    echo "[2/2] Running HiPhase (phase + haplotag) ..."
    { time hiphase \
        --reference "$REF" \
        --bam "$BAM" \
        --vcf "$DV_VCF" \
        $SV_ARGS \
        --output-bam "$HAPLOTAGGED_BAM" \
        --output-vcf "$PHASED_VCF" \
        $SV_OUT_ARGS \
        --stats-file "$STATS_TSV" \
        --blocks-file "$BLOCKS_TSV" \
        --summary-file "$SUMMARY_TSV" \
        --threads "$THREADS" ; } 2>> "$TIMING_LOG"

    samtools index -@ "$THREADS" "$HAPLOTAGGED_BAM"
else
    echo "[2/2] Haplotagged BAM exists, skipping."
fi

echo ""
echo "Done. Outputs in $OUTDIR/"
echo "  DeepVariant VCF: $DV_VCF"
echo "  Phased VCF:      $PHASED_VCF"
echo "  Haplotagged BAM: $HAPLOTAGGED_BAM"
echo "  Stats:           $STATS_TSV"
echo "  Blocks:          $BLOCKS_TSV"
echo "  Summary:         $SUMMARY_TSV"
echo "  Timing log:      $TIMING_LOG"
