#!/usr/bin/env bash
#
# bench_longphase.sh — Variant-call with DeepVariant, phase with LongPhase,
#                      and produce a haplotagged BAM.
#
# Pipeline:
#   1. DeepVariant (via Docker) — call variants from aligned BAM
#   2. longphase phase          — phase variants using reads
#   3. longphase haplotag       — tag each read with HP/PS
#
# Prerequisites:
#   - docker (for DeepVariant)
#   - longphase (https://github.com/twolinin/longphase/releases)
#   - samtools, bgzip, tabix
#
# Usage:
#   ./scripts/bench_longphase.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --outdir results/longphase \
#       --threads 16

set -euo pipefail

REF=""
BAM=""
VCF=""
SV_VCF=""
OUTDIR=""
THREADS=16
PLATFORM="pb"
INDELS=true
DV_VERSION="1.10.0"
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
  --sv-vcf FILE    SV VCF for co-phasing (from Sniffles or CuteSV)
  --ont            Use ONT mode instead of PacBio HiFi [default: --pb]
  --threads INT    Threads [16]
  --no-indels      Skip indel co-phasing
  --dv-version STR DeepVariant Docker tag [1.10.0]
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
        --ont)        PLATFORM="ont";  shift ;;
        --no-indels)  INDELS=false;    shift ;;
        --dv-version) DV_VERSION="$2"; shift 2 ;;
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
for cmd in longphase samtools bgzip tabix; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

REF_ABS="$(realpath "$REF")"
BAM_ABS="$(realpath "$BAM")"
REF_DIR="$(dirname "$REF_ABS")"
BAM_DIR="$(dirname "$BAM_ABS")"

DV_VCF="$OUTDIR/deepvariant.vcf.gz"
PREFIX="$OUTDIR/phased"
PHASED_VCF="${PREFIX}.vcf"
PHASED_VCF_GZ="${PREFIX}.vcf.gz"
HAPLOTAGGED_BAM="$OUTDIR/haplotagged.bam"
TIMING_LOG="$OUTDIR/timing.log"

INDEL_FLAG=""
[[ "$INDELS" == true ]] && INDEL_FLAG="--indels"

# ── Step 1: Call variants with DeepVariant ────────────────────────────
if [[ -n "$VCF" ]]; then
    echo "[1/3] Using provided VCF: $VCF"
    DV_VCF="$VCF"
elif [[ ! -f "$DV_VCF" ]]; then
    echo "[1/3] Calling variants with DeepVariant ${DV_VERSION} (model: ${DV_MODEL}) ..."
    command -v docker >/dev/null 2>&1 || { echo "Error: docker not found" >&2; exit 1; }
    OUTDIR_ABS="$(realpath "$OUTDIR")"
    { time docker run \
        -v "${REF_DIR}":"/ref_dir" \
        -v "${BAM_DIR}":"/bam_dir" \
        -v "${OUTDIR_ABS}":"/output" \
        google/deepvariant:"${DV_VERSION}" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="${DV_MODEL}" \
        --ref="/ref_dir/$(basename "$REF_ABS")" \
        --reads="/bam_dir/$(basename "$BAM_ABS")" \
        --output_vcf="/output/deepvariant.vcf.gz" \
        --num_shards="$THREADS" ; } 2> >(tee "$TIMING_LOG" >&2)
else
    echo "[1/3] DeepVariant VCF exists, skipping."
fi

# ── Step 2: Phase variants with LongPhase ─────────────────────────────
SV_PHASE_ARGS=""
SV_TAG_ARGS=""
if [[ -n "$SV_VCF" ]]; then
    [[ -f "$SV_VCF" ]] || { echo "Error: SV VCF not found: $SV_VCF" >&2; exit 1; }
    SV_PHASE_ARGS="--sv-file $SV_VCF"
    SV_TAG_ARGS="--sv-file ${PREFIX}_SV.vcf"
fi

if [[ ! -f "$PHASED_VCF" ]] && [[ ! -f "$PHASED_VCF_GZ" ]]; then
    echo "[2/3] Phasing variants with LongPhase ..."
    { time longphase phase \
        -s "$DV_VCF" \
        -b "$BAM" \
        -r "$REF" \
        -t "$THREADS" \
        -o "$PREFIX" \
        $INDEL_FLAG \
        $SV_PHASE_ARGS \
        --$PLATFORM ; } 2>> "$TIMING_LOG"
else
    echo "[2/3] Phased VCF exists, skipping."
fi

if [[ -f "$PHASED_VCF" ]] && [[ ! -f "$PHASED_VCF_GZ" ]]; then
    bgzip -c "$PHASED_VCF" > "$PHASED_VCF_GZ"
    tabix -p vcf "$PHASED_VCF_GZ"
fi

# ── Step 3: Haplotag reads ────────────────────────────────────────────
if [[ ! -f "$HAPLOTAGGED_BAM" ]]; then
    echo "[3/3] Haplotagging reads ..."
    { time longphase haplotag \
        -s "$PHASED_VCF" \
        -b "$BAM" \
        -r "$REF" \
        -t "$THREADS" \
        -o "$OUTDIR/haplotagged" \
        $SV_TAG_ARGS ; } 2>> "$TIMING_LOG"

    if [[ -f "$OUTDIR/haplotagged.bam" ]]; then
        samtools index -@ "$THREADS" "$OUTDIR/haplotagged.bam"
    fi
else
    echo "[3/3] Haplotagged BAM exists, skipping."
fi

echo ""
echo "Done. Outputs in $OUTDIR/"
echo "  DeepVariant VCF: $DV_VCF"
echo "  Phased VCF:      $PHASED_VCF_GZ"
echo "  Haplotagged BAM: $HAPLOTAGGED_BAM"
echo "  Timing log:      $TIMING_LOG"
