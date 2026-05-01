#!/usr/bin/env bash
# verify_refine_bam_parity.sh
#
# Quick regression check for phased BAM parity against longcallD in:
#   1) normal phased-output mode
#   2) --refine-aln phased-output mode
#
# Usage:
#   ./scripts/verify_refine_bam_parity.sh [REGION] [--hifi|--ont]
#
# Defaults:
#   REGION  CHM13#0#chr20:15200000-15205000 (HG002 HiFi chr20 slice; FA+BAM under test_data/chr20_quick/)
#   mode    --hifi
#
# Environment overrides:
#   PGPHASE      path to pgphase binary
#   LONGCALLD    path to longcallD binary
#   FA           reference FASTA
#   BAM          input BAM/CRAM (if unset, derived from mode)
#   OUTDIR       output directory for artifacts
#
# Exit codes:
#   0 parity pass (normal + refine BAM identical to longcallD after sort)
#   1 parity fail (normal or refine diff non-zero)
#   2 setup/runtime error

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
LONGCALLD="${LONGCALLD:-}"
REGION="${1:-CHM13#0#chr20:15200000-15205000}"
MODE="${2:---hifi}"
OUTDIR="${OUTDIR:-$ROOT/test_data/parity_compare/verify_refine}"

if [[ -z "${FA:-}" ]]; then
    if [[ "$MODE" == "--ont" ]]; then
        FA="$ROOT/test_data/chr11_2M.fa"
    else
        FA="$ROOT/test_data/chr20_quick/CHM13_chr20_only.fa"
    fi
fi

if [[ -z "$LONGCALLD" ]]; then
    if [[ -x "/tmp/longcalld/longcallD-v0.0.11_x64-linux/longcallD" ]]; then
        LONGCALLD="/tmp/longcalld/longcallD-v0.0.11_x64-linux/longcallD"
    else
        LONGCALLD="$ROOT/../longcallD/bin/longcallD"
    fi
fi

if [[ -z "${BAM:-}" ]]; then
    if [[ "$MODE" == "--ont" ]]; then
        BAM="$ROOT/test_data/HG002_chr11_ont_test.bam"
    else
        BAM="$ROOT/test_data/chr20_quick/HG002_CHM13_chr20_15000001_15500000.bam"
    fi
fi

fail_setup() { echo "error: $*" >&2; exit 2; }

command -v samtools >/dev/null 2>&1 || fail_setup "samtools not found in PATH"
[[ -x "$PGPHASE" ]] || fail_setup "pgphase binary not found: $PGPHASE (run: make)"
[[ -x "$LONGCALLD" ]] || fail_setup "longcallD binary not found: $LONGCALLD"
[[ -f "$FA" ]] || fail_setup "reference FASTA not found: $FA"
[[ -f "$BAM" ]] || fail_setup "input BAM/CRAM not found: $BAM"

mkdir -p "$OUTDIR"

echo "region : $REGION"
echo "mode   : $MODE"
echo "FA     : $FA"
echo "BAM    : $BAM"
echo "outdir : $OUTDIR"
echo ""

# Run pgphase
"$PGPHASE" collect-bam-variation -t 1 --include-filtered "$MODE" \
    -b "$OUTDIR/pgphase.bam" \
    "$FA" "$BAM" "$REGION" -o "$OUTDIR/pgphase.tsv"

"$PGPHASE" collect-bam-variation -t 1 --include-filtered "$MODE" --refine-aln \
    -b "$OUTDIR/pgphase_refine.bam" \
    "$FA" "$BAM" "$REGION" -o "$OUTDIR/pgphase_refine.tsv"

# Run longcallD
"$LONGCALLD" call -t1 "$MODE" \
    "$FA" "$BAM" -b "$OUTDIR/longcalld.bam" "$REGION" > "$OUTDIR/longcalld.vcf"

"$LONGCALLD" call -t1 "$MODE" --refine-aln \
    "$FA" "$BAM" -b "$OUTDIR/longcalld_refine.bam" "$REGION" > "$OUTDIR/longcalld_refine.vcf"

# Normalize order for stable comparison
samtools view "$OUTDIR/longcalld.bam" | LC_ALL=C sort > "$OUTDIR/longcalld.sam.sorted"
samtools view "$OUTDIR/pgphase.bam" | LC_ALL=C sort > "$OUTDIR/pgphase.sam.sorted"
samtools view "$OUTDIR/longcalld_refine.bam" | LC_ALL=C sort > "$OUTDIR/longcalld_refine.sam.sorted"
samtools view "$OUTDIR/pgphase_refine.bam" | LC_ALL=C sort > "$OUTDIR/pgphase_refine.sam.sorted"

diff -u "$OUTDIR/longcalld.sam.sorted" "$OUTDIR/pgphase.sam.sorted" > "$OUTDIR/diff.normal.txt" || true
diff -u "$OUTDIR/longcalld_refine.sam.sorted" "$OUTDIR/pgphase_refine.sam.sorted" > "$OUTDIR/diff.refine.txt" || true

normal_diff_lines=$(wc -l < "$OUTDIR/diff.normal.txt")
refine_diff_lines=$(wc -l < "$OUTDIR/diff.refine.txt")

echo "normal diff lines : $normal_diff_lines"
echo "refine diff lines : $refine_diff_lines"
echo ""

if [[ "$normal_diff_lines" -ne 0 || "$refine_diff_lines" -ne 0 ]]; then
    echo "FAIL: BAM parity mismatch detected."
    echo "See:"
    echo "  $OUTDIR/diff.normal.txt"
    echo "  $OUTDIR/diff.refine.txt"
    exit 1
fi

echo "OK: normal and refine BAM outputs match longcallD."
