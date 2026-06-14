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
#
# Add --gapfill [N] (implies --compare) to also emit hybrid_gapfill/phased.bam:
# the hybrid phased core with reads other pipelines phased but hybrid dropped,
# added back as disjoint phase sets. (All pipelines phase the same vg-giraffe
# alignment, so this is a cross-pipeline label transfer, not a re-alignment.)
#   N=1 (default): + BAM-pipeline reads           -> max accuracy gain vs BAM
#   N=2:           + BAM-pipeline + graph reads    -> max coverage (full union)
# This recovers the hard reads the hybrid skip_noisy_kmeans default drops, and
# Pareto-beats the BAM pipeline (more reads, lower error, higher auN). See
# CHECKPOINT.md "hybrid-core + gap-fill".
#
# Note: the N=1 (BAM-pipeline) gap-fill is now built into the binary as
# `collect-hybrid-variation --gap-fill` (no post-process / second pipeline run
# needed) and is the only supported gap-fill: it recovers reads present in the
# projected BAM. N=2 (graph-source) is a deliberate non-goal for the binary and
# stays here for offline analysis only; see CHECKPOINT.md "hybrid-core + gap-fill".

set -euo pipefail

readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REF=""
BAM=""
SITES=""
GAF=""
OUTDIR=""
THREADS=8
PLATFORM="hifi"   # hifi or ont
PGPHASE="./pgphase"
COMPARE=0
GAPFILL=0
GAPFILL_SOURCES=1
readonly GAPFILL_PS_OFFSET_BAM=1000000000
readonly GAPFILL_PS_OFFSET_GRAPH=2000000000

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
  --gapfill [N]    Emit hybrid_gapfill/phased.bam (implies --compare).
                   N=1 (default): + BAM-pipeline reads (max accuracy vs BAM).
                   N=2: + BAM + graph reads (max coverage, full union).
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
        --gapfill)
            GAPFILL=1; COMPARE=1
            if [[ "${2:-}" =~ ^[12]$ ]]; then GAPFILL_SOURCES="$2"; shift 2
            else shift; fi
            ;;
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

# ── Optional additive gap-fill: hybrid core + other-pipeline reads ────
# Stamp the reads other pipelines phased but hybrid dropped onto the hybrid
# phased BAM, each in its own disjoint phase-set namespace (no re-stitch).
# Needs the bam (and, for N=2, graph) pipeline output, which --gapfill forces
# by also setting --compare. Chained: source 2 fills onto source 1's result.
if [[ "${GAPFILL}" -eq 1 ]]; then
    gapfill_dir="${OUTDIR}/hybrid_gapfill"
    gapfill_bam="${gapfill_dir}/phased.bam"
    if [[ -f "${gapfill_bam}" ]]; then
        echo "[hybrid_gapfill] phased BAM exists, skipping."
    else
        command -v samtools >/dev/null 2>&1 || die "samtools required for --gapfill"
        mkdir -p "${gapfill_dir}"
        echo "[hybrid_gapfill] Building hybrid core + ${GAPFILL_SOURCES}-source gap-fill ..."
        # Source 1: BAM-pipeline reads (disjoint PS namespace).
        python3 "${SCRIPT_DIR}/gapfill.py" \
            "${OUTDIR}/hybrid/phased.bam" \
            "${OUTDIR}/bam/phased.bam" \
            "${gapfill_bam}" \
            samtools "${GAPFILL_PS_OFFSET_BAM}"
        # Source 2 (optional): graph-pipeline reads (second disjoint namespace).
        if [[ "${GAPFILL_SOURCES}" -eq 2 ]]; then
            stage1="${gapfill_dir}/stage1.bam"
            mv "${gapfill_bam}" "${stage1}"
            python3 "${SCRIPT_DIR}/gapfill.py" \
                "${stage1}" \
                "${OUTDIR}/graph/phased.bam" \
                "${gapfill_bam}" \
                samtools "${GAPFILL_PS_OFFSET_GRAPH}"
            rm -f "${stage1}"
        fi
        samtools index "${gapfill_bam}"
    fi
fi

echo ""
echo "Done. Outputs in ${OUTDIR}/"
