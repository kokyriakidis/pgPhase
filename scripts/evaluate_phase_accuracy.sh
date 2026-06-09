#!/usr/bin/env bash
#
# evaluate_phase_accuracy.sh
#
# Evaluate pgphase phase sets against the HG002 T2T diploid assembly
# using diplinator for haplotype-aware read assignment.
#
# Diplinator aligns reads to each haplotype separately and picks the
# best haplotype per read, avoiding the MAPQ penalty that a combined
# diploid reference causes.
#
# Requirements: samtools >= 1.17, minimap2 >= 2.26, diplinator,
#               python3 (+ matplotlib optional)
#
# ── Step 0: Generate inputs (uncomment and adjust paths) ─────────────
#
# # BAM pipeline (HiFi reads aligned to linear reference):
# ./pgphase collect-bam-variation \
#     --ref GRCh38.fa \
#     --bam HG002_hifi.bam \
#     --hifi \
#     --phased-vcf-out phased.vcf \
#     --out-bam phased.bam \
#     -o candidates.tsv \
#     -t 16
#
# # Graph pipeline (HiFi reads aligned to pangenome graph):
# # 1. Build site catalog from GBZ
# ./pgphase build-snarl-catalog \
#     --ref-sample CHM13 \
#     -o chr20.sites.vcf.gz \
#     chr20.gbz
#
# # 2. Build databases
# third_party/gbz-base/target/release/gbz2db chr20.gbz chr20.gbz.db
# third_party/gbz-base/target/release/gaf2db \
#     -r chr20.gbz -o HG002.gaf.db --overwrite HG002.gaf
#
# # 3. Run graph-based phasing
# ./pgphase collect-graph-variation \
#     --ref ref.fa \
#     --sites chr20.sites.vcf.gz \
#     --gbz-db chr20.gbz.db \
#     --gaf-db HG002.gaf.db \
#     --sample CHM13 \
#     --phased-vcf-out phased.vcf \
#     --phased-bam-out phased.bam \
#     -o candidates.tsv \
#     -t 16
#
# ── Then run this script ─────────────────────────────────────────────
#
# Usage:
#   ./evaluate_phase_accuracy.sh \
#       --phased-bam pgphase_phased.bam \
#       --reads      reads.fq.gz \
#       --mat-ref    HG002_T2T_maternal.fa \
#       --pat-ref    HG002_T2T_paternal.fa \
#       --outdir     eval_results \
#       --chroms     "chr1,chr20,chr22" \
#       --threads    16

set -euo pipefail

THREADS=16
OUTDIR=""
MIN_MAPQ=0
MIN_READS_PER_PS=5
MIN_HAPQ=0
MINIMAP2="minimap2"
DIPLINATOR="diplinator"
SAMTOOLS="samtools"
EXCLUDE_BED=""
CENSAT_BED=""
SEGDUP_BED=""
SITES_VCF=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --phased-bam  FILE    Phased uBAM from pgphase (HP/PS tags)
  --reads       FILE    FASTQ/FASTA with read sequences (.fq, .fq.gz, .fa, .fa.gz)
  --mat-ref     FILE    Maternal haplotype assembly FASTA
  --pat-ref     FILE    Paternal haplotype assembly FASTA
  --outdir      DIR     Output directory

Optional:
  --truth-bam   FILE    Pre-built diplinator truth BAM (from build_truth_bam.sh).
                        Skips alignment and diplinator steps.
  --chroms      STR     Comma-separated chromosome list [all]
  --threads     INT     Threads for minimap2/samtools/diplinator [16]
  --min-mapq    INT     Min MAPQ for truth alignment [0]
  --min-hapq    INT     Min HapQ from diplinator [0]
  --min-reads   INT     Min reads per phase set to evaluate [5]
  --exclude-bed FILE   BED of difficult regions to exclude from evaluation
  --censat-bed  FILE   Centromere/satellite BED for per-category accuracy
  --segdup-bed  FILE   Segmental duplication BED for per-category accuracy
  --sites-vcf   FILE   Sites VCF from graph decomposition (to detect unphaseable blocks)
  --samtools    PATH    Path to samtools executable [samtools]
  --minimap2    PATH    Path to minimap2 executable [minimap2]
  --diplinator  PATH    Path to diplinator executable [diplinator]
  -h, --help            Show this help
EOF
    exit 1
}

PHASED_BAM=""
TRUTH_BAM=""
READS=""
MAT_REF=""
PAT_REF=""
CHROMS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phased-bam)  PHASED_BAM="$2"; shift 2 ;;
        --truth-bam)   TRUTH_BAM="$2";  shift 2 ;;
        --reads)       READS="$2";      shift 2 ;;
        --mat-ref)     MAT_REF="$2";    shift 2 ;;
        --pat-ref)     PAT_REF="$2";    shift 2 ;;
        --chroms)      CHROMS="$2";     shift 2 ;;
        --threads)     THREADS="$2";    shift 2 ;;
        --min-mapq)    MIN_MAPQ="$2";   shift 2 ;;
        --min-hapq)    MIN_HAPQ="$2";   shift 2 ;;
        --min-reads)   MIN_READS_PER_PS="$2"; shift 2 ;;
        --exclude-bed) EXCLUDE_BED="$2"; shift 2 ;;
        --censat-bed)  CENSAT_BED="$2"; shift 2 ;;
        --segdup-bed)  SEGDUP_BED="$2"; shift 2 ;;
        --sites-vcf)   SITES_VCF="$2";  shift 2 ;;
        --samtools)    SAMTOOLS="$2";   shift 2 ;;
        --minimap2)    MINIMAP2="$2";   shift 2 ;;
        --diplinator)  DIPLINATOR="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2";     shift 2 ;;
        -h|--help)     usage ;;
        *)             echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "$PHASED_BAM" ]] && { echo "Error: --phased-bam required" >&2; usage; }
[[ -z "$OUTDIR" ]]     && { echo "Error: --outdir required" >&2; usage; }
[[ -f "$PHASED_BAM" ]] || { echo "Error: file not found: $PHASED_BAM" >&2; exit 1; }

if [[ -n "$TRUTH_BAM" ]]; then
    # Pre-built truth BAM: skip alignment and diplinator steps
    [[ -f "$TRUTH_BAM" ]] || { echo "Error: truth BAM not found: $TRUTH_BAM" >&2; exit 1; }
    for cmd in "$SAMTOOLS" python3; do
        command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH" >&2; exit 1; }
    done
else
    # Full pipeline: need reads, assemblies, minimap2, diplinator
    [[ -z "$READS" ]]   && { echo "Error: --reads required (or use --truth-bam)" >&2; usage; }
    [[ -z "$MAT_REF" ]] && { echo "Error: --mat-ref required (or use --truth-bam)" >&2; usage; }
    [[ -z "$PAT_REF" ]] && { echo "Error: --pat-ref required (or use --truth-bam)" >&2; usage; }
    for f in "$READS" "$MAT_REF" "$PAT_REF"; do
        [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
    done
    for cmd in "$SAMTOOLS" "$MINIMAP2" "$DIPLINATOR" python3 paftools.js; do
        if [[ "$cmd" == */* ]]; then
            [[ -x "$cmd" ]] || { echo "Error: $cmd not found or not executable" >&2; exit 1; }
        else
            command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH" >&2; exit 1; }
        fi
    done
fi

mkdir -p "$OUTDIR"

# ── Steps 1-4: Build truth BAM (or use pre-built) ────────────────────
if [[ -n "$TRUTH_BAM" ]]; then
    MERGED_BAM="$TRUTH_BAM"
    echo "Using pre-built truth BAM: $MERGED_BAM"
else
    MERGED_BAM="$OUTDIR/diplinator_merged.bam"
    if [[ ! -f "$MERGED_BAM" ]]; then
        echo "Building truth BAM (alignment + diplinator) ..."
        "$(dirname "$0")/build_truth_bam.sh" \
            --reads "$READS" \
            --mat-ref "$MAT_REF" \
            --pat-ref "$PAT_REF" \
            --outdir "$OUTDIR" \
            --threads "$THREADS" \
            --samtools "$SAMTOOLS" \
            --minimap2 "$MINIMAP2" \
            --diplinator "$DIPLINATOR"
    else
        echo "Truth BAM exists: $MERGED_BAM"
    fi
fi

# ── Step 5: Inject HP/PS/YC tags into merged BAM ─────────────────────
TAGGED_BAM="$OUTDIR/diplinator_merged.tagged.bam"

if [[ ! -f "$TAGGED_BAM" ]]; then
    echo "[5/6] Injecting HP/PS/YC tags ..."
    python3 "$(dirname "$0")/inject_hp_tags.py" \
        "$PHASED_BAM" "$MERGED_BAM" "$TAGGED_BAM" "$SAMTOOLS"
    "$SAMTOOLS" index -@ "$THREADS" "$TAGGED_BAM"
else
    echo "[5/6] Tagged BAM exists, skipping."
fi

# ── step 7: evaluate + plot ───────────────────────────────────────────
echo "[6/6] Evaluating phase concordance ..."
python3 "$(dirname "$0")/evaluate_phase_accuracy.py" \
    "$PHASED_BAM" "$TAGGED_BAM" "$MIN_MAPQ" "$MIN_HAPQ" \
    "$MIN_READS_PER_PS" "$CHROMS" "$OUTDIR" "$SAMTOOLS" "$EXCLUDE_BED" \
    "$SITES_VCF" "$CENSAT_BED" "$SEGDUP_BED"

echo "      Generating plots ..."
python3 "$(dirname "$0")/plot_phase_accuracy.py" "$OUTDIR" 2>/dev/null \
    && echo "  Plots saved to $OUTDIR/" \
    || echo "  Skipped (matplotlib not available)."

echo ""
echo "Done. See $OUTDIR/ for results."
