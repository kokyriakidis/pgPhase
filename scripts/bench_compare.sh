#!/usr/bin/env bash
#
# bench_compare.sh — Compare read-level phasing accuracy across tools.
#
# Primary evaluation: diplinator-based read-level accuracy on each
# tool's haplotagged BAM (all tools, including pgphase).
# Optional: whatshap compare on phased VCFs (all tools, requires
# --truth-vcf), gene phasing completeness (requires --gene-bed).
#
# Prerequisites:
#   - All bench_*.sh scripts have been run and produced output
#   - diplinator, minimap2, samtools, python3
#   - HG002 diploid assembly (maternal + paternal FASTAs)
#   - Reads FASTQ (for diplinator alignment)
#   - whatshap (optional, only for VCF-level comparison)
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
GENE_BED=""
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
  --truth-vcf FILE     Truth phased VCF for VCF-level comparison (competitors only)
  --gene-bed FILE      Gene BED (chrom,start,end,name) for % fully phased genes
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
        --gene-bed)      GENE_BED="$2";      shift 2 ;;
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
        echo "  WARNING: $tool: no phased BAM found in $dir, skipping."
        continue
    fi

    eval_dir="$OUTDIR/${tool}_eval"
    if [[ -f "$eval_dir/summary.json" ]]; then
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

# ── Helper: find the phased VCF in a tool's output directory ──────────
find_phased_vcf() {
    local dir="$1"
    for name in phased.vcf.gz phased.vcf; do
        [[ -f "$dir/$name" ]] && { echo "$dir/$name"; return 0; }
    done
    echo ""
    return 1
}

# Helper: extract a field from a JSON file
json_field() {
    local json="$1" field="$2"
    python3 -c "import json; d=json.load(open('$json')); print(d.get('$field','N/A'))" 2>/dev/null || echo "N/A"
}

# ── Part 2: VCF-level metrics ────────────────────────────────────────
echo ""
echo "================================================================"
echo " Part 2: VCF-level metrics"
echo "================================================================"

# 2a: Phased variant count + gene completeness per tool
for tool in pgphase whatshap hiphase longphase; do
    [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
    dir="${TOOL_DIRS[$tool]}"
    vcf=$(find_phased_vcf "$dir") || true

    if [[ -z "$vcf" ]]; then
        echo "  WARNING: $tool: no phased VCF found in $dir, skipping VCF metrics."
        continue
    fi

    gene_json="$OUTDIR/${tool}_gene_phasing.json"
    if [[ -n "$GENE_BED" ]] && [[ -f "$GENE_BED" ]] && [[ ! -f "$gene_json" ]]; then
        echo ""
        echo "── Gene phasing completeness: $tool ──"
        python3 "$SCRIPT_DIR/gene_phasing_completeness.py" \
            "$vcf" "$GENE_BED" > "$gene_json"
    elif [[ -f "$gene_json" ]]; then
        echo "  $tool: gene phasing exists, skipping."
    fi
done

# 2b: whatshap compare on phased VCFs (requires --truth-vcf)
if [[ -n "$TRUTH_VCF" ]]; then
    echo ""
    echo "── VCF-level phasing accuracy (whatshap compare) ──"

    [[ -f "$TRUTH_VCF" ]] || { echo "Error: truth VCF not found: $TRUTH_VCF" >&2; exit 1; }
    if ! command -v whatshap >/dev/null 2>&1; then
        echo "  Warning: whatshap not found, skipping VCF comparison." >&2
    else
        for tool in pgphase whatshap hiphase longphase; do
            [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
            dir="${TOOL_DIRS[$tool]}"
            vcf=$(find_phased_vcf "$dir") || true

            if [[ -z "$vcf" ]]; then
                echo "  WARNING: $tool: no phased VCF, skipping whatshap compare."
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
                "$TRUTH_VCF" "$vcf"
        done
    fi
fi

# ── Part 3: Summary table ────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Part 3: Summary"
echo "================================================================"

{
    printf "tool\taccuracy\thamming_err\tswitch_err\tswitches\tflips\tswitchflips\t"
    printf "phased_reads\tfrac_phased\tng50_bp\t"
    printf "phased_het_vars\tpct_genes_phased\twall_clock\n"
    for tool in pgphase whatshap hiphase longphase; do
        [[ -z "${TOOL_DIRS[$tool]+x}" ]] && continue
        eval_json="$OUTDIR/${tool}_eval/summary.json"
        gene_json="$OUTDIR/${tool}_gene_phasing.json"
        timing="${TOOL_DIRS[$tool]}/timing.log"

        acc="N/A"; hamming="N/A"; switch_rate="N/A"
        switches="N/A"; flips="N/A"; switchflips="N/A"
        reads="N/A"; frac="N/A"; ng50="N/A"
        phased_vars="N/A"; pct_genes="N/A"; wall="N/A"

        if [[ -f "$eval_json" ]]; then
            acc=$(json_field "$eval_json" overall_accuracy)
            hamming=$(json_field "$eval_json" hamming_error_rate)
            switch_rate=$(json_field "$eval_json" switch_error_rate)
            switches=$(json_field "$eval_json" switch_errors)
            flips=$(json_field "$eval_json" flip_errors)
            switchflips=$(json_field "$eval_json" switchflip_errors)
            reads=$(json_field "$eval_json" total_phased_reads)
            frac=$(json_field "$eval_json" fraction_phased)
            ng50=$(json_field "$eval_json" phase_block_ng50_bp)
        fi

        if [[ -f "$gene_json" ]]; then
            phased_vars=$(json_field "$gene_json" phased_het_variants)
            pct_genes=$(json_field "$gene_json" pct_genes_fully_phased)
        fi

        if [[ -f "$timing" ]]; then
            wall=$(grep -oP '(?<=real\s)\S+' "$timing" | tail -1 || echo "N/A")
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$tool" "$acc" "$hamming" "$switch_rate" \
            "$switches" "$flips" "$switchflips" \
            "$reads" "$frac" "$ng50" \
            "$phased_vars" "$pct_genes" "$wall"
    done
} | tee "$SUMMARY"

echo ""
echo "Done. Full results in $OUTDIR/"
echo "  Summary table:    $SUMMARY"
echo "  Per-tool eval:    $OUTDIR/<tool>_eval/"
[[ -n "$GENE_BED" ]] && echo "  Gene phasing:     $OUTDIR/<tool>_gene_phasing.json"
[[ -n "$TRUTH_VCF" ]] && echo "  VCF comparisons:  $OUTDIR/<tool>_whatshap_compare.tsv"
