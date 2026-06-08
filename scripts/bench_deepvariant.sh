#!/usr/bin/env bash
#
# bench_deepvariant.sh — Call variants with DeepVariant (via Docker).
#
# Produces a VCF that can be passed to the phasing benchmark scripts
# via --vcf, so DeepVariant only needs to run once.
#
# Prerequisites:
#   - docker
#
# Usage:
#   ./scripts/bench_deepvariant.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --outdir results/deepvariant \
#       --threads 16
#
# Then pass the output to each phaser:
#   ./scripts/bench_whatshap.sh  --vcf results/deepvariant/deepvariant.vcf.gz ...
#   ./scripts/bench_hiphase.sh   --vcf results/deepvariant/deepvariant.vcf.gz ...
#   ./scripts/bench_longphase.sh --vcf results/deepvariant/deepvariant.vcf.gz ...

set -euo pipefail

REF=""
BAM=""
OUTDIR=""
THREADS=16
DV_VERSION="1.10.0"
DV_MODEL="PACBIO"
REGIONS=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --ref FILE       Reference FASTA (indexed)
  --bam FILE       Aligned BAM/CRAM (indexed)
  --outdir DIR     Output directory

Optional:
  --threads INT    Threads / num_shards [16]
  --dv-version STR DeepVariant Docker tag [1.10.0]
  --dv-model STR   Model type [PACBIO]
                   One of: WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA
  --regions STR    Restrict to region (e.g. chr20)
  -h, --help       Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)        REF="$2";        shift 2 ;;
        --bam)        BAM="$2";        shift 2 ;;
        --outdir)     OUTDIR="$2";     shift 2 ;;
        --threads)    THREADS="$2";    shift 2 ;;
        --dv-version) DV_VERSION="$2"; shift 2 ;;
        --dv-model)   DV_MODEL="$2";   shift 2 ;;
        --regions)    REGIONS="$2";    shift 2 ;;
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
command -v docker >/dev/null 2>&1 || { echo "Error: docker not found in PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"

REF_ABS="$(realpath "$REF")"
BAM_ABS="$(realpath "$BAM")"
REF_DIR="$(dirname "$REF_ABS")"
BAM_DIR="$(dirname "$BAM_ABS")"
OUTDIR_ABS="$(realpath "$OUTDIR")"

DV_VCF="$OUTDIR/deepvariant.vcf.gz"
DV_GVCF="$OUTDIR/deepvariant.g.vcf.gz"
TIMING_LOG="$OUTDIR/timing.log"

REGION_ARGS=""
[[ -n "$REGIONS" ]] && REGION_ARGS="--regions=$REGIONS"

if [[ ! -f "$DV_VCF" ]]; then
    echo "Calling variants with DeepVariant ${DV_VERSION} (model: ${DV_MODEL}) ..."
    [[ -n "$REGIONS" ]] && echo "  Region: $REGIONS"
    echo "  Threads: $THREADS"

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
        --output_gvcf="/output/deepvariant.g.vcf.gz" \
        --num_shards="$THREADS" \
        $REGION_ARGS ; } 2> >(tee "$TIMING_LOG" >&2)
else
    echo "DeepVariant VCF already exists: $DV_VCF"
fi

echo ""
echo "Done. Outputs in $OUTDIR/"
echo "  VCF:        $DV_VCF"
echo "  gVCF:       $DV_GVCF"
echo "  Timing log: $TIMING_LOG"
echo ""
echo "Pass to phasing scripts:"
echo "  ./scripts/bench_whatshap.sh  --vcf $DV_VCF --ref $REF --bam $BAM --outdir results/whatshap"
echo "  ./scripts/bench_hiphase.sh   --vcf $DV_VCF --ref $REF --bam $BAM --outdir results/hiphase"
echo "  ./scripts/bench_longphase.sh --vcf $DV_VCF --ref $REF --bam $BAM --outdir results/longphase"
