#!/usr/bin/env bash
#
# bench_pgphase.sh — Phase HG002 HiFi reads with pgphase and produce a phased BAM.
#
# pgphase performs variant calling and phasing in a single step — no
# pre-called VCF is needed. It supports two pipelines:
#   - BAM pipeline: reads aligned to a linear reference
#   - Graph pipeline: reads aligned to a pangenome graph (coordinate-indexed GAF)
#
# Prerequisites:
#   - pgphase (built from this repo)
#   - samtools (for BAM indexing)
#
# Usage (BAM pipeline):
#   ./scripts/bench_pgphase.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --outdir results/pgphase_bam \
#       --threads 16
#
# Usage (graph pipeline):
#   ./scripts/bench_pgphase.sh \
#       --ref REF.fa \
#       --gaf READS.coord.gaf.gz \
#       --sites SITES.vcf.gz \
#       --outdir results/pgphase_graph \
#       --threads 16

set -euo pipefail

REF=""
BAM=""
GAF=""
SITES=""
OUTDIR=""
THREADS=16
PLATFORM="hifi"  # hifi or ont
PGPHASE="./pgphase"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --ref FILE       Reference FASTA (indexed)
  --outdir DIR     Output directory

BAM pipeline (provide --bam):
  --bam FILE       Aligned BAM/CRAM (indexed)

Graph pipeline (provide --gaf + --sites):
  --gaf FILE       Coordinate-indexed GAF (bgzipped + tabix)
  --sites FILE     Sites VCF (bgzipped + tabix-indexed)

Optional:
  --ont            ONT mode [default: HiFi]
  --pgphase PATH   Path to pgphase binary [./pgphase]
  --threads INT    Threads [16]
  -h, --help       Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)       REF="$2";      shift 2 ;;
        --bam)       BAM="$2";      shift 2 ;;
        --gaf)       GAF="$2";      shift 2 ;;
        --sites)     SITES="$2";    shift 2 ;;
        --outdir)    OUTDIR="$2";   shift 2 ;;
        --threads)   THREADS="$2";  shift 2 ;;
        --ont)       PLATFORM="ont"; shift ;;
        --pgphase)   PGPHASE="$2";  shift 2 ;;
        -h|--help)   usage ;;
        *)           echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "$REF" ]]    && { echo "Error: --ref required" >&2; usage; }
[[ -z "$OUTDIR" ]] && { echo "Error: --outdir required" >&2; usage; }
[[ -z "$BAM" ]] && [[ -z "$GAF" ]] && { echo "Error: provide --bam or --gaf" >&2; usage; }

[[ -f "$REF" ]] || { echo "Error: file not found: $REF" >&2; exit 1; }
[[ -x "$PGPHASE" ]] || [[ -f "$PGPHASE" ]] || command -v "$PGPHASE" >/dev/null 2>&1 || \
    { echo "Error: pgphase not found: $PGPHASE" >&2; exit 1; }

mkdir -p "$OUTDIR"

PHASED_BAM="$OUTDIR/phased.bam"
PHASED_VCF="$OUTDIR/phased.vcf"
CANDIDATES_TSV="$OUTDIR/candidates.tsv"
TIMING_LOG="$OUTDIR/timing.log"

PLATFORM_FLAG="--$PLATFORM"

if [[ -n "$BAM" ]]; then
    # ── BAM pipeline ──────────────────────────────────────────────────
    [[ -f "$BAM" ]] || { echo "Error: file not found: $BAM" >&2; exit 1; }

    if [[ ! -f "$PHASED_BAM" ]]; then
        echo "[1/1] Running pgphase collect-bam-variation ..."
        { time "$PGPHASE" collect-bam-variation \
            --ref "$REF" \
            --bam "$BAM" \
            $PLATFORM_FLAG \
            --out-bam "$PHASED_BAM" \
            --phased-vcf-out "$PHASED_VCF" \
            -o "$CANDIDATES_TSV" \
            -t "$THREADS" ; } 2> >(tee "$TIMING_LOG" >&2)
    else
        echo "[1/1] Phased BAM exists, skipping."
    fi

elif [[ -n "$GAF" ]]; then
    # ── Graph pipeline ────────────────────────────────────────────────
    [[ -f "$GAF" ]]   || { echo "Error: file not found: $GAF" >&2; exit 1; }
    [[ -n "$SITES" ]] || { echo "Error: --sites required for graph pipeline" >&2; usage; }
    [[ -f "$SITES" ]] || { echo "Error: file not found: $SITES" >&2; exit 1; }

    if [[ ! -f "$PHASED_BAM" ]]; then
        echo "[1/1] Running pgphase collect-graph-variation ..."
        { time "$PGPHASE" collect-graph-variation \
            --ref "$REF" \
            --gaf "$GAF" \
            --sites "$SITES" \
            $PLATFORM_FLAG \
            --phased-bam-out "$PHASED_BAM" \
            --phased-vcf-out "$PHASED_VCF" \
            -o "$CANDIDATES_TSV" \
            -t "$THREADS" ; } 2> >(tee "$TIMING_LOG" >&2)
    else
        echo "[1/1] Phased BAM exists, skipping."
    fi
fi

echo ""
echo "Done. Outputs in $OUTDIR/"
echo "  Phased BAM:    $PHASED_BAM"
echo "  Phased VCF:    $PHASED_VCF"
echo "  Candidates:    $CANDIDATES_TSV"
echo "  Timing log:    $TIMING_LOG"
