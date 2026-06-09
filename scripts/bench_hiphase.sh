#!/usr/bin/env bash
#
# bench_hiphase.sh — Variant-call with DeepVariant, pbsv, and TRGT,
#                    then jointly phase with HiPhase.
#
# Pipeline:
#   1. DeepVariant — small variant calling (SNVs + indels)
#   2. pbsv       — structural variant calling (optional, requires --trf-bed)
#   3. TRGT       — tandem repeat genotyping (optional, requires --trgt-repeats)
#   4. hiphase    — joint phasing + haplotagging
#
# Modes:
#   "HiPhase (no SV)" — only DeepVariant VCF (default)
#   "HiPhase"         — DeepVariant + pbsv + TRGT with global re-alignment
#
# Environment setup (run once):
#   micromamba create -f envs/hiphase.yaml -y
#   micromamba activate bench-hiphase
#
# Usage:
#   # HiPhase (no SV) — small variants only
#   ./scripts/bench_hiphase.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --outdir results/hiphase_nosv \
#       --threads 16
#
#   # HiPhase — joint phasing with SVs and TRs
#   ./scripts/bench_hiphase.sh \
#       --ref REF.fa \
#       --bam ALIGNED.bam \
#       --trf-bed human_GRCh38_no_alt_analysis_set.trf.bed \
#       --trgt-repeats adotto_repeats.hg38.bed \
#       --outdir results/hiphase \
#       --threads 16
#
# Annotation files for "HiPhase" mode:
#
#   --trf-bed: Tandem repeat annotation for pbsv discover. Improves SV
#     sensitivity in repetitive regions. Download the GRCh38 version from
#     the pbsv repo:
#     https://github.com/PacificBiosciences/pbsv/raw/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
#
#   --trgt-repeats: Genome-wide repeat catalog for TRGT genotyping.
#     The "Adotto" catalog used in the HiPhase paper (Holt et al. 2024)
#     is available at:
#     https://zenodo.org/records/8329210
#     Download: adotto_TRregions_v1.2_GRCh38.bed.gz (gunzip before use)

set -euo pipefail

REF=""
BAM=""
VCF=""
SV_VCF=""
TR_VCF=""
TRF_BED=""
TRGT_REPEATS=""
KARYOTYPE="XX"
OUTDIR=""
THREADS=16
DV_MODEL="PACBIO"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --ref FILE           Reference FASTA (indexed)
  --bam FILE           Aligned BAM/CRAM (indexed)
  --outdir DIR         Output directory

Optional (small variants):
  --vcf FILE           Pre-called small variant VCF (skip DeepVariant)
  --dv-model STR       DeepVariant model type [PACBIO]

Optional (structural variants — enables "HiPhase" mode):
  --sv-vcf FILE        Pre-called SV VCF (skip pbsv)
  --trf-bed FILE       Tandem repeat annotation BED for pbsv discover

Optional (tandem repeats — enables "HiPhase" mode):
  --tr-vcf FILE        Pre-called TR VCF (skip TRGT)
  --trgt-repeats FILE  TRGT repeat catalog BED

Optional:
  --karyotype STR      Sample karyotype for TRGT [XX]
  --threads INT        Threads [16]
  -h, --help           Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)           REF="$2";           shift 2 ;;
        --bam)           BAM="$2";           shift 2 ;;
        --vcf)           VCF="$2";           shift 2 ;;
        --sv-vcf)        SV_VCF="$2";        shift 2 ;;
        --tr-vcf)        TR_VCF="$2";        shift 2 ;;
        --trf-bed)       TRF_BED="$2";       shift 2 ;;
        --trgt-repeats)  TRGT_REPEATS="$2";  shift 2 ;;
        --karyotype)     KARYOTYPE="$2";     shift 2 ;;
        --outdir)        OUTDIR="$2";        shift 2 ;;
        --threads)       THREADS="$2";       shift 2 ;;
        --dv-model)      DV_MODEL="$2";      shift 2 ;;
        -h|--help)       usage ;;
        *)               echo "Unknown option: $1" >&2; usage ;;
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

# Determine mode
RUN_PBSV=false
RUN_TRGT=false
[[ -n "$TRF_BED" || -n "$SV_VCF" ]] && RUN_PBSV=true
[[ -n "$TRGT_REPEATS" || -n "$TR_VCF" ]] && RUN_TRGT=true

if $RUN_PBSV || $RUN_TRGT; then
    echo "Mode: HiPhase (joint phasing with SVs/TRs)"
else
    echo "Mode: HiPhase (no SV)"
fi

# Count total steps
TOTAL_STEPS=2
$RUN_PBSV && TOTAL_STEPS=$((TOTAL_STEPS + 1))
$RUN_TRGT && TOTAL_STEPS=$((TOTAL_STEPS + 1))
STEP=0

next_step() { STEP=$((STEP + 1)); }

# ── Step: Call variants with DeepVariant ──────────────────────────────
next_step
if [[ -n "$VCF" ]]; then
    echo "[$STEP/$TOTAL_STEPS] Using provided VCF: $VCF"
    DV_VCF="$VCF"
elif [[ ! -f "$DV_VCF" ]]; then
    echo "[$STEP/$TOTAL_STEPS] Calling variants with DeepVariant (model: ${DV_MODEL}) ..."
    command -v run_deepvariant >/dev/null 2>&1 || \
        { echo "Error: run_deepvariant not found. Activate bench-hiphase env." >&2; exit 1; }
    { time run_deepvariant \
        --model_type="${DV_MODEL}" \
        --ref="$REF" \
        --reads="$BAM" \
        --output_vcf="$DV_VCF" \
        --num_shards="$THREADS" ; } 2> >(tee "$TIMING_LOG" >&2)
else
    echo "[$STEP/$TOTAL_STEPS] DeepVariant VCF exists, skipping."
fi

# ── Step: Call SVs with pbsv ──────────────────────────────────────────
PBSV_VCF="$OUTDIR/pbsv.vcf.gz"
if $RUN_PBSV; then
    next_step
    if [[ -n "$SV_VCF" ]]; then
        echo "[$STEP/$TOTAL_STEPS] Using provided SV VCF: $SV_VCF"
        PBSV_VCF="$SV_VCF"
    elif [[ ! -f "$PBSV_VCF" ]]; then
        echo "[$STEP/$TOTAL_STEPS] Calling structural variants with pbsv ..."
        command -v pbsv >/dev/null 2>&1 || \
            { echo "Error: pbsv not found in PATH" >&2; exit 1; }

        SVSIG="$OUTDIR/pbsv.svsig.gz"
        TRF_ARGS=""
        [[ -n "$TRF_BED" ]] && TRF_ARGS="--tandem-repeats $TRF_BED"

        if [[ ! -f "$SVSIG" ]]; then
            { time pbsv discover \
                $TRF_ARGS \
                "$BAM" "$SVSIG" ; } 2>> "$TIMING_LOG"
        fi

        { time pbsv call --ccs \
            -j "$THREADS" \
            "$REF" "$SVSIG" "$OUTDIR/pbsv.vcf" ; } 2>> "$TIMING_LOG"

        bgzip -c "$OUTDIR/pbsv.vcf" > "$PBSV_VCF"
        tabix -p vcf "$PBSV_VCF"
        rm -f "$OUTDIR/pbsv.vcf"
    else
        echo "[$STEP/$TOTAL_STEPS] pbsv VCF exists, skipping."
    fi
fi

# ── Step: Genotype tandem repeats with TRGT ───────────────────────────
TRGT_VCF="$OUTDIR/trgt.vcf.gz"
if $RUN_TRGT; then
    next_step
    if [[ -n "$TR_VCF" ]]; then
        echo "[$STEP/$TOTAL_STEPS] Using provided TR VCF: $TR_VCF"
        TRGT_VCF="$TR_VCF"
    elif [[ ! -f "$TRGT_VCF" ]]; then
        echo "[$STEP/$TOTAL_STEPS] Genotyping tandem repeats with TRGT ..."
        command -v trgt >/dev/null 2>&1 || \
            { echo "Error: trgt not found in PATH" >&2; exit 1; }

        TRGT_PREFIX="$OUTDIR/trgt"

        { time trgt genotype \
            --genome "$REF" \
            --reads "$BAM" \
            --repeats "$TRGT_REPEATS" \
            --output-prefix "$TRGT_PREFIX" \
            --karyotype "$KARYOTYPE" \
            --threads "$THREADS" ; } 2>> "$TIMING_LOG"

        # TRGT outputs unsorted VCF; sort and index
        bcftools sort -Oz -o "$TRGT_VCF" "$TRGT_PREFIX.vcf.gz"
        tabix -p vcf "$TRGT_VCF"
    else
        echo "[$STEP/$TOTAL_STEPS] TRGT VCF exists, skipping."
    fi
fi

# ── Step: Phase + haplotag with HiPhase ───────────────────────────────
next_step

# Build hiphase arguments for SV/TR VCFs
EXTRA_VCF_ARGS=""
EXTRA_OUT_ARGS=""
GLOBAL_REALIGN_ARGS=""

if $RUN_PBSV; then
    EXTRA_VCF_ARGS="$EXTRA_VCF_ARGS --vcf $PBSV_VCF"
    EXTRA_OUT_ARGS="$EXTRA_OUT_ARGS --output-vcf $OUTDIR/phased_sv.vcf.gz"
fi
if $RUN_TRGT; then
    EXTRA_VCF_ARGS="$EXTRA_VCF_ARGS --vcf $TRGT_VCF"
    EXTRA_OUT_ARGS="$EXTRA_OUT_ARGS --output-vcf $OUTDIR/phased_tr.vcf.gz"
fi
if $RUN_PBSV || $RUN_TRGT; then
    GLOBAL_REALIGN_ARGS="--global-realignment-cputime 300"
fi

if [[ ! -f "$HAPLOTAGGED_BAM" ]]; then
    echo "[$STEP/$TOTAL_STEPS] Running HiPhase (phase + haplotag) ..."
    { time hiphase \
        --reference "$REF" \
        --bam "$BAM" \
        --vcf "$DV_VCF" \
        $EXTRA_VCF_ARGS \
        $GLOBAL_REALIGN_ARGS \
        --output-bam "$HAPLOTAGGED_BAM" \
        --output-vcf "$PHASED_VCF" \
        $EXTRA_OUT_ARGS \
        --stats-file "$STATS_TSV" \
        --blocks-file "$BLOCKS_TSV" \
        --summary-file "$SUMMARY_TSV" \
        --threads "$THREADS" ; } 2>> "$TIMING_LOG"

    samtools index -@ "$THREADS" "$HAPLOTAGGED_BAM"
else
    echo "[$STEP/$TOTAL_STEPS] Haplotagged BAM exists, skipping."
fi

echo ""
echo "Done. Outputs in $OUTDIR/"
echo "  DeepVariant VCF: $DV_VCF"
$RUN_PBSV && echo "  pbsv VCF:        $PBSV_VCF"
$RUN_TRGT && echo "  TRGT VCF:        $TRGT_VCF"
echo "  Phased VCF:      $PHASED_VCF"
echo "  Haplotagged BAM: $HAPLOTAGGED_BAM"
echo "  Stats:           $STATS_TSV"
echo "  Blocks:          $BLOCKS_TSV"
echo "  Summary:         $SUMMARY_TSV"
echo "  Timing log:      $TIMING_LOG"
