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
  --chroms      STR     Comma-separated chromosome list [all]
  --threads     INT     Threads for minimap2/samtools/diplinator [16]
  --min-mapq    INT     Min MAPQ for truth alignment [0]
  --min-hapq    INT     Min HapQ from diplinator [0]
  --min-reads   INT     Min reads per phase set to evaluate [5]
  --exclude-bed FILE   BED of difficult regions to exclude from evaluation
  --samtools    PATH    Path to samtools executable [samtools]
  --minimap2    PATH    Path to minimap2 executable [minimap2]
  --diplinator  PATH    Path to diplinator executable [diplinator]
  -h, --help            Show this help
EOF
    exit 1
}

PHASED_BAM=""
READS=""
MAT_REF=""
PAT_REF=""
CHROMS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phased-bam)  PHASED_BAM="$2"; shift 2 ;;
        --reads)       READS="$2";      shift 2 ;;
        --mat-ref)     MAT_REF="$2";    shift 2 ;;
        --pat-ref)     PAT_REF="$2";    shift 2 ;;
        --chroms)      CHROMS="$2";     shift 2 ;;
        --threads)     THREADS="$2";    shift 2 ;;
        --min-mapq)    MIN_MAPQ="$2";   shift 2 ;;
        --min-hapq)    MIN_HAPQ="$2";   shift 2 ;;
        --min-reads)   MIN_READS_PER_PS="$2"; shift 2 ;;
        --exclude-bed) EXCLUDE_BED="$2"; shift 2 ;;
        --samtools)    SAMTOOLS="$2";   shift 2 ;;
        --minimap2)    MINIMAP2="$2";   shift 2 ;;
        --diplinator)  DIPLINATOR="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2";     shift 2 ;;
        -h|--help)     usage ;;
        *)             echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "$PHASED_BAM" ]] && { echo "Error: --phased-bam required" >&2; usage; }
[[ -z "$READS" ]]      && { echo "Error: --reads required" >&2; usage; }
[[ -z "$MAT_REF" ]]    && { echo "Error: --mat-ref required" >&2; usage; }
[[ -z "$PAT_REF" ]]    && { echo "Error: --pat-ref required" >&2; usage; }
[[ -z "$OUTDIR" ]]     && { echo "Error: --outdir required" >&2; usage; }
for f in "$PHASED_BAM" "$READS" "$MAT_REF" "$PAT_REF"; do
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

# ── step 1: index haplotype references ────────────────────────────────
echo "[1/6] Checking inputs and indices ..."
for ref in "$MAT_REF" "$PAT_REF"; do
    if [[ ! -f "${ref}.fai" ]]; then
        "$SAMTOOLS" faidx "$ref"
    fi
done
echo "  MAT: $(wc -l < "${MAT_REF}.fai") contigs"
echo "  PAT: $(wc -l < "${PAT_REF}.fai") contigs"
echo "  Reads: $READS"
echo "  Phased uBAM: $PHASED_BAM"

# ── step 3: align to each haplotype separately ───────────────────────
# Diplinator requires name-sorted input.  Write SAM directly (minimap2
# outputs in query-name order) so diplinator sees the correct format.
# Also produce BAM and PAF for downstream use.
MAT_SAM="$OUTDIR/vs_mat.sam"
PAT_SAM="$OUTDIR/vs_pat.sam"
MAT_BAM="$OUTDIR/vs_mat.bam"
PAT_BAM="$OUTDIR/vs_pat.bam"
MAT_PAF="$OUTDIR/vs_mat.paf"
PAT_PAF="$OUTDIR/vs_pat.paf"

if [[ ! -f "$MAT_SAM" ]]; then
    echo "[2/6] Mapping to maternal haplotype ..."
    "$MINIMAP2" -Y --cs --eqx -ax lr:hqae -e 100 --secondary=no \
        -t "$THREADS" "$MAT_REF" "$READS" \
        -o "$MAT_SAM"
else
    echo "[2/6] Maternal SAM exists, skipping."
fi

if [[ ! -f "$PAT_SAM" ]]; then
    echo "      Mapping to paternal haplotype ..."
    "$MINIMAP2" -Y --cs --eqx -ax lr:hqae -e 100 --secondary=no \
        -t "$THREADS" "$PAT_REF" "$READS" \
        -o "$PAT_SAM"
else
    echo "      Paternal SAM exists, skipping."
fi

# Convert SAM → BAM (coordinate-sorted, for downstream merge/tagging/IGV)
if [[ ! -f "$MAT_BAM" ]]; then
    echo "      Converting maternal SAM to sorted BAM ..."
    "$SAMTOOLS" sort -@ "$THREADS" -o "$MAT_BAM" "$MAT_SAM"
    "$SAMTOOLS" index -@ "$THREADS" "$MAT_BAM"
else
    echo "      Maternal BAM exists, skipping."
fi

if [[ ! -f "$PAT_BAM" ]]; then
    echo "      Converting paternal SAM to sorted BAM ..."
    "$SAMTOOLS" sort -@ "$THREADS" -o "$PAT_BAM" "$PAT_SAM"
    "$SAMTOOLS" index -@ "$THREADS" "$PAT_BAM"
else
    echo "      Paternal BAM exists, skipping."
fi

# Convert SAM → PAF (paftools.js sam2paf preserves ms:i: tag;
# diplinator can use it via --paf --ms).
if [[ ! -f "$MAT_PAF" ]]; then
    echo "      Converting maternal SAM to PAF ..."
    paftools.js sam2paf "$MAT_SAM" > "$MAT_PAF"
else
    echo "      Maternal PAF exists, skipping."
fi

if [[ ! -f "$PAT_PAF" ]]; then
    echo "      Converting paternal SAM to PAF ..."
    paftools.js sam2paf "$PAT_SAM" > "$PAT_PAF"
else
    echo "      Paternal PAF exists, skipping."
fi

# ── step 4: run diplinator ───────────────────────────────────────────
DIPLO_MAT="$OUTDIR/diplinator_mat.sam"
DIPLO_PAT="$OUTDIR/diplinator_pat.sam"

if [[ ! -f "$DIPLO_MAT" ]] || [[ ! -f "$DIPLO_PAT" ]]; then
    echo "[3/6] Running diplinator ..."
    # Diplinator writes output files in the current directory as
    # diplinator_<label>.sam, so run from OUTDIR.
    MAT_ABS="$(realpath "$MAT_SAM")"
    PAT_ABS="$(realpath "$PAT_SAM")"
    # Resolve diplinator to an absolute path so it works after cd.
    if [[ "$DIPLINATOR" == */* ]]; then
        DIPLO_ABS="$(realpath "$DIPLINATOR")"
    else
        DIPLO_ABS="$(command -v "$DIPLINATOR")"
    fi
    (cd "$OUTDIR" && "$DIPLO_ABS" -1 mat -2 pat -t "$THREADS" "$MAT_ABS" "$PAT_ABS")
else
    echo "[3/6] Diplinator output exists, skipping."
fi

# ── step 5: merge diplinator outputs with haplotype-of-origin tag ────
MERGED_BAM="$OUTDIR/diplinator_merged.bam"

if [[ ! -f "$MERGED_BAM" ]]; then
    echo "[4/6] Merging diplinator outputs ..."
    # Tag each read with HO:Z:MAT or HO:Z:PAT so the evaluation script
    # knows which haplotype diplinator assigned it to.
    # Build a merged header.  Both haplotype FASTAs may share contig
    # names (e.g. chr1) but differ in length.  Keep the longer LN for
    # each SN so that coordinates from either haplotype are valid.
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
        # Maternal reads
        "$SAMTOOLS" view "$DIPLO_MAT" | awk 'BEGIN{OFS="\t"}{print $0, "HO:Z:MAT"}'
        # Paternal reads
        "$SAMTOOLS" view "$DIPLO_PAT" | awk 'BEGIN{OFS="\t"}{print $0, "HO:Z:PAT"}'
    } | "$SAMTOOLS" sort -@ "$THREADS" -o "$MERGED_BAM"
    rm -f "$MAT_HDR" "$PAT_HDR" "$MERGED_HDR"
    "$SAMTOOLS" index -@ "$THREADS" "$MERGED_BAM"
    echo "  $("$SAMTOOLS" view -c "$MERGED_BAM") reads in merged BAM."
else
    echo "[4/6] Merged BAM exists, skipping."
fi

# ── step 6: inject HP/PS/YC tags into merged BAM ─────────────────────
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
    "$MIN_READS_PER_PS" "$CHROMS" "$OUTDIR" "$SAMTOOLS" "$EXCLUDE_BED"

echo "      Generating plots ..."
python3 "$(dirname "$0")/plot_phase_accuracy.py" "$OUTDIR" 2>/dev/null \
    && echo "  Plots saved to $OUTDIR/" \
    || echo "  Skipped (matplotlib not available)."

echo ""
echo "Done. See $OUTDIR/ for results."
