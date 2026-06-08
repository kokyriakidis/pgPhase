#!/usr/bin/env bash
#
# bench_whatshap.sh — Variant-call with DeepVariant, phase with WhatsHap,
#                     and produce a haplotagged BAM.
#
# Pipeline:
#   1. DeepVariant (via Docker) — call variants from aligned BAM
#   2. whatshap phase           — phase heterozygous variants using reads
#   3. whatshap haplotag        — tag each read with HP/PS
#
# Environment setup (run once):
#   micromamba create -f envs/whatshap.yaml -y
#   micromamba activate bench-whatshap
#
# Prerequisites:
#   - docker (for DeepVariant)
#   - micromamba environment bench-whatshap (whatshap, samtools, bgzip, tabix)
#
# Usage:
#   ./scripts/bench_whatshap.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --outdir results/whatshap \
#       --threads 16

set -euo pipefail

REF=""
BAM=""
VCF=""
OUTDIR=""
THREADS=16
SAMPLE=""
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
  --sample NAME    Sample name in VCF [auto-detected]
  --threads INT    Threads [16]
  --no-indels      Skip indel phasing
  --dv-version STR DeepVariant Docker tag [1.10.0]
  --dv-model STR   DeepVariant model type [PACBIO]
                   One of: WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA
  -h, --help       Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)        REF="$2";        shift 2 ;;
        --bam)        BAM="$2";        shift 2 ;;
        --vcf)        VCF="$2";        shift 2 ;;
        --outdir)     OUTDIR="$2";     shift 2 ;;
        --sample)     SAMPLE="$2";     shift 2 ;;
        --threads)    THREADS="$2";    shift 2 ;;
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
for cmd in whatshap samtools bgzip tabix; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

REF_ABS="$(realpath "$REF")"
BAM_ABS="$(realpath "$BAM")"
REF_DIR="$(dirname "$REF_ABS")"
BAM_DIR="$(dirname "$BAM_ABS")"

DV_VCF="$OUTDIR/deepvariant.vcf.gz"
PHASED_VCF="$OUTDIR/phased.vcf.gz"
HAPLOTAGGED_BAM="$OUTDIR/haplotagged.bam"
STATS_TSV="$OUTDIR/stats.tsv"
TIMING_LOG="$OUTDIR/timing.log"

SAMPLE_ARGS=""
[[ -n "$SAMPLE" ]] && SAMPLE_ARGS="--sample $SAMPLE"

INDEL_FLAG=""
[[ "$INDELS" == true ]] && INDEL_FLAG="--indels"

# ── Step 1: Call variants with DeepVariant ────────────────────────────
if [[ -n "$VCF" ]]; then
    echo "[1/4] Using provided VCF: $VCF"
    DV_VCF="$VCF"
elif [[ ! -f "$DV_VCF" ]]; then
    echo "[1/4] Calling variants with DeepVariant ${DV_VERSION} (model: ${DV_MODEL}) ..."
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
    echo "[1/4] DeepVariant VCF exists, skipping."
fi

# ── Step 2: Phase variants with WhatsHap ──────────────────────────────
if [[ ! -f "$PHASED_VCF" ]]; then
    echo "[2/4] Phasing variants with WhatsHap ..."
    { time whatshap phase \
        --reference "$REF" \
        $INDEL_FLAG \
        $SAMPLE_ARGS \
        -o "$OUTDIR/phased.vcf" \
        "$DV_VCF" "$BAM" ; } 2>> "$TIMING_LOG"

    bgzip -c "$OUTDIR/phased.vcf" > "$PHASED_VCF"
    tabix -p vcf "$PHASED_VCF"
    rm -f "$OUTDIR/phased.vcf"
else
    echo "[2/4] Phased VCF exists, skipping."
fi

# ── Step 3: Haplotag reads ────────────────────────────────────────────
if [[ ! -f "$HAPLOTAGGED_BAM" ]]; then
    echo "[3/4] Haplotagging reads ..."
    { time whatshap haplotag \
        --reference "$REF" \
        $SAMPLE_ARGS \
        --output-threads "$THREADS" \
        -o "$HAPLOTAGGED_BAM" \
        "$PHASED_VCF" "$BAM" ; } 2>> "$TIMING_LOG"

    samtools index -@ "$THREADS" "$HAPLOTAGGED_BAM"
else
    echo "[3/4] Haplotagged BAM exists, skipping."
fi

# ── Step 4: Phasing statistics ────────────────────────────────────────
echo "[4/4] Computing phasing statistics ..."
whatshap stats --tsv="$STATS_TSV" "$PHASED_VCF"

echo ""
echo "Done. Outputs in $OUTDIR/"
echo "  DeepVariant VCF: $DV_VCF"
echo "  Phased VCF:      $PHASED_VCF"
echo "  Haplotagged BAM: $HAPLOTAGGED_BAM"
echo "  Stats TSV:       $STATS_TSV"
echo "  Timing log:      $TIMING_LOG"
