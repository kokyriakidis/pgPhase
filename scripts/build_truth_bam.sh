#!/usr/bin/env bash
#
# build_truth_bam.sh — Build a diplinator truth BAM from reads and
# diploid assembly haplotypes.
#
# Aligns reads to each haplotype separately, runs diplinator to pick
# the best haplotype per read, and merges the results into a single
# BAM with HO:Z:MAT/PAT tags.
#
# This BAM is tool-independent and can be reused across all phasing
# tools via evaluate_phase_accuracy.sh --truth-bam.
#
# Requirements: samtools >= 1.17, minimap2 >= 2.26, diplinator, python3
#
# Usage:
#   ./scripts/build_truth_bam.sh \
#       --reads reads.fq.gz \
#       --mat-ref HG002.mat.fa \
#       --pat-ref HG002.pat.fa \
#       --outdir truth_bam \
#       --threads 16

set -euo pipefail

THREADS=16
OUTDIR=""
MINIMAP2="minimap2"
DIPLINATOR="diplinator"
SAMTOOLS="samtools"
READS=""
MAT_REF=""
PAT_REF=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --reads       FILE    FASTQ/FASTA with read sequences
  --mat-ref     FILE    Maternal haplotype assembly FASTA
  --pat-ref     FILE    Paternal haplotype assembly FASTA
  --outdir      DIR     Output directory

Optional:
  --threads     INT     Threads [16]
  --samtools    PATH    Path to samtools [samtools]
  --minimap2    PATH    Path to minimap2 [minimap2]
  --diplinator  PATH    Path to diplinator [diplinator]
  -h, --help            Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)       READS="$2";       shift 2 ;;
        --mat-ref)     MAT_REF="$2";     shift 2 ;;
        --pat-ref)     PAT_REF="$2";     shift 2 ;;
        --outdir)      OUTDIR="$2";      shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        --samtools)    SAMTOOLS="$2";    shift 2 ;;
        --minimap2)    MINIMAP2="$2";    shift 2 ;;
        --diplinator)  DIPLINATOR="$2";  shift 2 ;;
        -h|--help)     usage ;;
        *)             echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "$READS" ]]   && { echo "Error: --reads required" >&2; usage; }
[[ -z "$MAT_REF" ]] && { echo "Error: --mat-ref required" >&2; usage; }
[[ -z "$PAT_REF" ]] && { echo "Error: --pat-ref required" >&2; usage; }
[[ -z "$OUTDIR" ]]  && { echo "Error: --outdir required" >&2; usage; }
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

mkdir -p "$OUTDIR"

MERGED_BAM="$OUTDIR/diplinator_merged.bam"

if [[ -f "$MERGED_BAM" ]]; then
    echo "Truth BAM already exists: $MERGED_BAM"
    exit 0
fi

# ── Step 1: Index haplotype references ────────────────────────────────
echo "[1/4] Checking inputs and indices ..."
for ref in "$MAT_REF" "$PAT_REF"; do
    if [[ ! -f "${ref}.fai" ]]; then
        "$SAMTOOLS" faidx "$ref"
    fi
done
echo "  MAT: $(wc -l < "${MAT_REF}.fai") contigs"
echo "  PAT: $(wc -l < "${PAT_REF}.fai") contigs"
echo "  Reads: $READS"

# ── Step 2: Align to each haplotype separately ───────────────────────
MAT_SAM="$OUTDIR/vs_mat.sam"
PAT_SAM="$OUTDIR/vs_pat.sam"

if [[ ! -f "$MAT_SAM" ]]; then
    echo "[2/4] Mapping to maternal haplotype ..."
    "$MINIMAP2" -Y --cs --eqx -ax lr:hqae -e 100 --secondary=no \
        -t "$THREADS" "$MAT_REF" "$READS" \
        -o "$MAT_SAM"
else
    echo "[2/4] Maternal SAM exists, skipping."
fi

if [[ ! -f "$PAT_SAM" ]]; then
    echo "      Mapping to paternal haplotype ..."
    "$MINIMAP2" -Y --cs --eqx -ax lr:hqae -e 100 --secondary=no \
        -t "$THREADS" "$PAT_REF" "$READS" \
        -o "$PAT_SAM"
else
    echo "      Paternal SAM exists, skipping."
fi

# ── Step 3: Run diplinator ────────────────────────────────────────────
DIPLO_MAT="$OUTDIR/diplinator_mat.sam"
DIPLO_PAT="$OUTDIR/diplinator_pat.sam"

if [[ ! -f "$DIPLO_MAT" ]] || [[ ! -f "$DIPLO_PAT" ]]; then
    echo "[3/4] Running diplinator ..."
    MAT_ABS="$(realpath "$MAT_SAM")"
    PAT_ABS="$(realpath "$PAT_SAM")"
    if [[ "$DIPLINATOR" == */* ]]; then
        DIPLO_ABS="$(realpath "$DIPLINATOR")"
    else
        DIPLO_ABS="$(command -v "$DIPLINATOR")"
    fi
    (cd "$OUTDIR" && "$DIPLO_ABS" -1 mat -2 pat -t "$THREADS" "$MAT_ABS" "$PAT_ABS")
else
    echo "[3/4] Diplinator output exists, skipping."
fi

# ── Step 4: Merge with haplotype-of-origin tag ───────────────────────
echo "[4/4] Merging diplinator outputs ..."

MAT_HDR="$OUTDIR/_hdr_mat.sam"
PAT_HDR="$OUTDIR/_hdr_pat.sam"
MERGED_HDR="$OUTDIR/_hdr_merged.sam"
"$SAMTOOLS" view -H "$DIPLO_MAT" > "$MAT_HDR"
"$SAMTOOLS" view -H "$DIPLO_PAT" > "$PAT_HDR"

python3 -c "
import sys
sq = {}   # SN -> (LN, full_line)
sq_order = []
other = []
for path in sys.argv[1:]:
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('@SQ'):
                parts = line.split('\t')
                sn = ln = None
                for p in parts:
                    if p.startswith('SN:'): sn = p[3:]
                    if p.startswith('LN:'): ln = int(p[3:])
                if sn is None: continue
                if sn not in sq:
                    sq[sn] = (ln or 0, line)
                    sq_order.append(sn)
                elif (ln or 0) > sq[sn][0]:
                    sq[sn] = (ln or 0, line)
            elif line.startswith('@'):
                if line not in other:
                    other.append(line)
for o in other:
    print(o)
for sn in sq_order:
    print(sq[sn][1])
" "$MAT_HDR" "$PAT_HDR" > "$MERGED_HDR"

{
    cat "$MERGED_HDR"
    "$SAMTOOLS" view "$DIPLO_MAT" | awk 'BEGIN{OFS="\t"}{print $0, "HO:Z:MAT"}'
    "$SAMTOOLS" view "$DIPLO_PAT" | awk 'BEGIN{OFS="\t"}{print $0, "HO:Z:PAT"}'
} | "$SAMTOOLS" sort -@ "$THREADS" -o "$MERGED_BAM"
rm -f "$MAT_HDR" "$PAT_HDR" "$MERGED_HDR"
"$SAMTOOLS" index -@ "$THREADS" "$MERGED_BAM"
echo "  $("$SAMTOOLS" view -c "$MERGED_BAM") reads in merged BAM."

echo ""
echo "Done. Truth BAM: $MERGED_BAM"
echo "Use with: evaluate_phase_accuracy.sh --truth-bam $MERGED_BAM ..."
