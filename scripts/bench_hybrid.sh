#!/usr/bin/env bash
#
# bench_hybrid.sh — Run pgphase in hybrid (BAM+graph) mode and, optionally,
# the BAM-only and graph-only pipelines on the same inputs for A/B comparison.
#
# This is the Phase 0 measurement harness for the hybrid-model work: it produces
# directly comparable phased BAM/VCF outputs from each pipeline so their phasing
# metrics can be evaluated against the same truth set.
#
# pgphase performs variant calling and phasing in a single step. Hybrid mode
# augments BAM variant calling with graph snarl site observations, so it needs
# a matched BAM + graph-sites VCF + coordinate-indexed GAF on the same contig.
#
# Prerequisites:
#   - pgphase (built from this repo)
#   - samtools (for BAM indexing)
#
# Usage:
#   ./scripts/bench_hybrid.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --sites SITES.vcf.gz \
#       --gaf READS.coord.gaf.gz \
#       --outdir results/hybrid \
#       --threads 8
#
# Add --compare to also run the BAM-only and graph-only pipelines into
# sibling subdirectories (bam/ and graph/) for head-to-head comparison.

set -euo pipefail

readonly SCRIPT_NAME="$(basename "$0")"

REF=""
BAM=""
SITES=""
GAF=""
OUTDIR=""
THREADS=8
PLATFORM="hifi"   # hifi or ont
PGPHASE="./pgphase"
COMPARE=0

die() {
    echo "Error: $*" >&2
    exit 1
}

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options]

Required:
  --ref FILE       Reference FASTA (indexed)
  --bam FILE       Aligned BAM/CRAM (indexed)
  --sites FILE     Graph sites VCF (bgzipped + tabix-indexed)
  --gaf FILE       Coordinate-indexed GAF (bgzipped + tabix-indexed)
  --outdir DIR     Output directory

Optional:
  --compare        Also run BAM-only and graph-only pipelines for A/B
  --ont            ONT mode [default: HiFi]
  --pgphase PATH   Path to pgphase binary [./pgphase]
  --threads INT    Threads [8]
  -h, --help       Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)     REF="$2";     shift 2 ;;
        --bam)     BAM="$2";     shift 2 ;;
        --sites)   SITES="$2";   shift 2 ;;
        --gaf)     GAF="$2";     shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --compare) COMPARE=1;    shift ;;
        --ont)     PLATFORM="ont"; shift ;;
        --pgphase) PGPHASE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *)         echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -n "${REF}" ]]    || die "--ref required"
[[ -n "${BAM}" ]]    || die "--bam required"
[[ -n "${SITES}" ]]  || die "--sites required"
[[ -n "${GAF}" ]]    || die "--gaf required"
[[ -n "${OUTDIR}" ]] || die "--outdir required"

for f in "${REF}" "${BAM}" "${SITES}" "${GAF}"; do
    [[ -f "${f}" ]] || die "file not found: ${f}"
done

# pgphase requires an indexed BAM/CRAM for region chunking.  BAM indexes are
# generated artifacts (not committed), so build one if it is missing.
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" && ! -f "${BAM}.csi" ]]; then
    command -v samtools >/dev/null 2>&1 || die "samtools required to index ${BAM}"
    echo "Indexing ${BAM} ..."
    samtools index "${BAM}"
fi
[[ -x "${PGPHASE}" ]] || [[ -f "${PGPHASE}" ]] || command -v "${PGPHASE}" >/dev/null 2>&1 \
    || die "pgphase not found: ${PGPHASE}"

readonly PLATFORM_FLAG="--${PLATFORM}"

# run_pipeline <subdir> <command> <bam_flag> [extra args...]
# bam_flag differs by pipeline: graph mode uses --phased-bam-out, while the
# bam and hybrid modes use --out-bam.
run_pipeline() {
    local subdir="$1"; shift
    local cmd="$1"; shift
    local bam_flag="$1"; shift
    local dir="${OUTDIR}/${subdir}"
    mkdir -p "${dir}"

    local phased_bam="${dir}/phased.bam"
    local phased_vcf="${dir}/phased.vcf"
    local candidates="${dir}/candidates.tsv"
    local timing="${dir}/timing.log"

    if [[ -f "${phased_bam}" ]]; then
        echo "[${subdir}] phased BAM exists, skipping."
        return 0
    fi

    echo "[${subdir}] Running pgphase ${cmd} ..."
    local start end
    start=$(date +%s)
    "${PGPHASE}" "${cmd}" \
        --ref "${REF}" \
        "${PLATFORM_FLAG}" \
        "${bam_flag}" "${phased_bam}" \
        --phased-vcf-out "${phased_vcf}" \
        -o "${candidates}" \
        -t "${THREADS}" \
        "$@" 2> >(tee "${timing}" >&2)
    end=$(date +%s)
    echo "[${subdir}] wall time: $((end - start))s" | tee -a "${timing}"
}

mkdir -p "${OUTDIR}"

# ── Hybrid pipeline ──────────────────────────────────────────────────
run_pipeline "hybrid" "collect-hybrid-variation" "--out-bam" \
    --bam "${BAM}" --graph-sites "${SITES}" --gaf "${GAF}"

# ── Optional A/B: BAM-only and graph-only ────────────────────────────
if [[ "${COMPARE}" -eq 1 ]]; then
    run_pipeline "bam" "collect-bam-variation" "--out-bam" --bam "${BAM}"
    run_pipeline "graph" "collect-graph-variation" "--phased-bam-out" \
        --gaf "${GAF}" --sites "${SITES}"
fi

echo ""
echo "Done. Outputs in ${OUTDIR}/"
