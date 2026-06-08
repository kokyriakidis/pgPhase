#!/usr/bin/env bash
#
# bench_compare.sh — Compare phasing accuracy across all benchmarked tools.
#
# Runs the diplinator-based evaluation on each tool's haplotagged BAM,
# runs whatshap compare on each competitor's phased VCF (not pgphase),
# collects timing, and produces a summary table.
#
# Prerequisites:
#   - All bench_*.sh scripts have been run and produced output
#   - diplinator, minimap2, samtools, whatshap, python3
#   - HG002 diploid assembly (maternal + paternal FASTAs)
#   - Reads FASTQ (for diplinator alignment)
#
# Usage:
#   ./scripts/bench_compare.sh \
#       --reads HG002_hifi.fq.gz \
#       --mat-ref HG002_mat.fa \
#       --pat-ref HG002_pat.fa \
#       --truth-vcf HG002_truth_phased.vcf.gz \
#       --outdir results/comparison \
#       --threads 16 \
#       --pgphase-dir results/pgphase \
#       --whatshap-dir results/whatshap \
#       --hiphase-dir results/hiphase \
#       --longphase-dir results/longphase

set -euo pipefail

READS=""
MAT_REF=""
PAT_REF=""
TRUTH_VCF=""
OUTDIR=""
THREADS=16
CHROMS=""
PGPHASE_DIR=""
WHATSHAP_DIR=""
HIPHASE_DIR=""
LONGPHASE_DIR=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --reads FILE         FASTQ/FASTA with read sequences
  --mat-ref FILE       Maternal haplotype assembly FASTA
  --pat-ref FILE       Paternal haplotype assembly FASTA
  --outdir DIR         Output directory for comparison results

Tool result directories (at least one required):
  --pgphase-dir DIR    pgphase output directory (contains phased.bam)
  --whatshap-dir DIR   WhatsHap output directory (contains haplotagged.bam, phased.vcf.gz)
  --hiphase-dir DIR    HiPhase output directory (contains haplotagged.bam, phased.vcf.gz)
  --longphase-dir DIR  LongPhase output directory (contains haplotagged.bam, phased.vcf.gz)

Optional:
  --truth-vcf FILE     Truth phased VCF for whatshap compare (competitors only)
  --chroms STR         Comma-separated chromosome list [all]
  --threads INT        Threads [16]
  -h, --help           Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)         READS="$2";         shift 2 ;;
        --mat-ref)       MAT_REF="$2";       shift 2 ;;
        --pat-ref)       PAT_REF="$2";       shift 2 ;;
        --truth-vcf)     TRUTH_VCF="$2";     shift 2 ;;
        --outdir)        OUTDIR="$2";        shift 2 ;;
        --threads)       THREADS="$2";       shift 2 ;;
        --chroms)        CHROMS="$2";        shift 2 ;;
        --pgphase-dir)   PGPHASE_DIR="$2";   shift 2 ;;
        --whatshap-dir)  WHATSHAP_DIR="$2";   shift 2 ;;
        --hiphase-dir)   HIPHASE_DIR="$2";   shift 2 ;;
        --longphase-dir) LONGPHASE_DIR="$2"; shift 2 ;;
        -h|--help)       usage ;;
        *)               echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "$READS" ]]   && { echo "Error: --reads required" >&2; usage; }
[[ -z "$MAT_REF" ]] && { echo "Error: --mat-ref required" >&2; usage; }
[[ -z "$PAT_REF" ]] && { echo "Error: --pat-ref required" >&2; usage; }
[[ -z "$OUTDIR" ]]  && { echo "Error: --outdir required" >&2; usage; }

for f in "$READS" "$MAT_REF" "$PAT_REF"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done

# At least one tool directory must be provided.
[[ -z "$PGPHASE_DIR" ]] && [[ -z "$WHATSHAP_DIR" ]] && \
[[ -z "$HIPHASE_DIR" ]] && [[ -z "$LONGPHASE_DIR" ]] && \
    { echo "Error: provide at least one --*-dir" >&2; usage; }

mkdir -p "$OUTDIR"

SUMMARY="$OUTDIR/summary.tsv"

# ── Helper: find the phased BAM in a tool's output directory ──────────
find_phased_bam() {
    local dir="$1"
    local tool="$2"
    if [[ "$tool" == "pgphase" ]]; then
        # pgphase uses phased.bam
        for name in phased.bam; do
            [[ -f "$dir/$name" ]] && { echo "$dir/$name"; return 0; }
        done
    else
        # Competitors use haplotagged.bam
        for name in haplotagged.bam; do
            [[ -f "$dir/$name" ]] && { echo "$dir/$name"; return 0; }
        done
    fi
    echo ""
    return 1
}

# ── Part 1: Diplinator-based read-level evaluation (all tools) ────────
echo "================================================================"
echo " Part 1: Read-level phasing accuracy (diplinator)"
echo "================================================================"

declare -A TOOL_DIRS
[[ -n "$PGPHASE_DIR" ]]   && TOOL_DIRS[pgphase]="$PGPHASE_DIR"
[[ -n "$WHATSHAP_DIR" ]]   && TOOL_DIRS[whatshap]="$WHATSHAP_DIR"
[[ -n "$HIPHASE_DIR" ]]    && TOOL_DIRS[hiphase]="$HIPHASE_DIR"
[[ -n "$LONGPHASE_DIR" ]]  && TOOL_DIRS[longphase]="$LONGPHASE_DIR"

for tool in pgphase whatshap hiphase longphase; do
    [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
    dir="${TOOL_DIRS[$tool]}"
    bam=$(find_phased_bam "$dir" "$tool") || true

    if [[ -z "$bam" ]]; then
        echo "  ⚠️  $tool: no phased BAM found in $dir, skipping."
        continue
    fi

    eval_dir="$OUTDIR/${tool}_eval"
    if [[ -f "$eval_dir/phase_accuracy.json" ]]; then
        echo "  $tool: evaluation exists, skipping."
        continue
    fi

    echo ""
    echo "── Evaluating $tool ──"
    CHROMS_ARG=""
    [[ -n "$CHROMS" ]] && CHROMS_ARG="--chroms $CHROMS"

    "$SCRIPT_DIR/evaluate_phase_accuracy.sh" \
        --phased-bam "$bam" \
        --reads "$READS" \
        --mat-ref "$MAT_REF" \
        --pat-ref "$PAT_REF" \
        --outdir "$eval_dir" \
        --threads "$THREADS" \
        $CHROMS_ARG
done

# ── Part 2: whatshap compare on phased VCFs (competitors only) ────────
if [[ -n "$TRUTH_VCF" ]]; then
    echo ""
    echo "================================================================"
    echo " Part 2: VCF-level phasing accuracy (whatshap compare)"
    echo "================================================================"

    [[ -f "$TRUTH_VCF" ]] || { echo "Error: truth VCF not found: $TRUTH_VCF" >&2; exit 1; }
    command -v whatshap >/dev/null 2>&1 || { echo "Warning: whatshap not found, skipping VCF comparison." >&2; }

    for tool in whatshap hiphase longphase; do
        [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
        dir="${TOOL_DIRS[$tool]}"
        phased_vcf="$dir/phased.vcf.gz"

        if [[ ! -f "$phased_vcf" ]]; then
            echo "  ⚠️  $tool: no phased VCF found in $dir, skipping."
            continue
        fi

        compare_out="$OUTDIR/${tool}_whatshap_compare.tsv"
        if [[ -f "$compare_out" ]]; then
            echo "  $tool: whatshap compare exists, skipping."
            continue
        fi

        echo ""
        echo "── whatshap compare: $tool ──"
        whatshap compare \
            --tsv-pairwise "$compare_out" \
            --names truth,$tool \
            "$TRUTH_VCF" "$phased_vcf"
    done
fi

# ── Part 3: Collect timing ────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Part 3: Timing summary"
echo "================================================================"

{
    printf "tool\twall_clock\n"
    for tool in pgphase whatshap hiphase longphase; do
        [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
        dir="${TOOL_DIRS[$tool]}"
        timing="$dir/timing.log"
        if [[ -f "$timing" ]]; then
            # Extract the last "real" line from the timing log.
            wall=$(grep -oP '(?<=real\s)\S+' "$timing" | tail -1 || echo "N/A")
            printf "%s\t%s\n" "$tool" "$wall"
        else
            printf "%s\tN/A\n" "$tool"
        fi
    done
} | tee "$OUTDIR/timing_summary.tsv"

# ── Part 4: Aggregate summary table ──────────────────────────────────
echo ""
echo "================================================================"
echo " Part 4: Summary"
echo "================================================================"

{
    printf "tool\tswitch_rate\thamming_rate\tphased_reads\tphase_block_n50\twall_clock\n"
    for tool in pgphase whatshap hiphase longphase; do
        [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
        eval_dir="$OUTDIR/${tool}_eval"
        json="$eval_dir/phase_accuracy.json"
        timing="${TOOL_DIRS[$tool]}/timing.log"

        switch="N/A"; hamming="N/A"; reads="N/A"; n50="N/A"; wall="N/A"

        if [[ -f "$json" ]]; then
            switch=$(python3 -c "import json; d=json.load(open('$json')); print(d.get('switch_error_rate','N/A'))" 2>/dev/null || echo "N/A")
            hamming=$(python3 -c "import json; d=json.load(open('$json')); print(d.get('hamming_error_rate','N/A'))" 2>/dev/null || echo "N/A")
            reads=$(python3 -c "import json; d=json.load(open('$json')); print(d.get('total_phased_reads','N/A'))" 2>/dev/null || echo "N/A")
            n50=$(python3 -c "import json; d=json.load(open('$json')); print(d.get('phase_block_n50','N/A'))" 2>/dev/null || echo "N/A")
        fi

        if [[ -f "$timing" ]]; then
            wall=$(grep -oP '(?<=real\s)\S+' "$timing" | tail -1 || echo "N/A")
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tool" "$switch" "$hamming" "$reads" "$n50" "$wall"
    done
} | tee "$SUMMARY"

echo ""
echo "Done. Full results in $OUTDIR/"
echo "  Summary table:    $SUMMARY"
echo "  Timing:           $OUTDIR/timing_summary.tsv"
echo "  Per-tool eval:    $OUTDIR/<tool>_eval/"
[[ -n "$TRUTH_VCF" ]] && echo "  VCF comparisons:  $OUTDIR/<tool>_whatshap_compare.tsv"
