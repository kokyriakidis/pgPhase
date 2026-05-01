#!/usr/bin/env bash
# Split CHM13#0#chr20 into ~CHUNK bases (default 20 Mb) windows and run
# compare_vcf_parity.sh + compare_phased_vcf.sh on each (same FA/BAM/MODE).
#
# Env:
#   ROOT        repo root (default: parent of scripts/)
#   FA          reference FASTA with .fai (default: pgbam-experiments CHM13 chr20 FA)
#   BAM         reads aligned to CHM13 (default: pgbam HG002.chr20.bam)
#   CONTIG      contig name (default: CHM13#0#chr20)
#   CHUNK       window size in bp (default: 20000000)
#   MODE        --hifi | --ont   (default: --hifi)
#   PGPHASE     LONGCALLD        optional binary paths
#
# Example:
#   FA=.../chm13v2.0.chr20.renamed.fa BAM=.../HG002.chr20.bam bash scripts/run_chm13_chr20_20mb_parity.sh
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FA="${FA:-/home/kokyriakidis/Downloads/pgbam-experiments/pgphase_ont_compare/chm13v2.0.chr20.renamed.fa}"
BAM="${BAM:-/home/kokyriakidis/Downloads/pgbam-experiments/HG002.chr20.bam}"
CONTIG="${CONTIG:-CHM13#0#chr20}"
CHUNK="${CHUNK:-20000000}"
MODE="${MODE:---hifi}"

fail() { echo "error: $*" >&2; exit 2; }

[[ -f "$FA" ]]  || fail "FASTA not found: $FA"
[[ -f "${FA}.fai" ]] || fail "FASTA index missing: ${FA}.fai"
[[ -f "$BAM" ]] || fail "BAM not found: $BAM"

export FA BAM
export PGPHASE="${PGPHASE:-$ROOT/pgphase}"
export LONGCALLD="${LONGCALLD:-$ROOT/../longcallD/bin/longcallD}"

[[ -x "$PGPHASE" ]] || fail "pgphase not executable: $PGPHASE"
[[ -x "$LONGCALLD" ]] || fail "longcallD not executable: $LONGCALLD"

CHRLEN="$(awk -v c="$CONTIG" '$1 == c { print $2; exit }' "${FA}.fai")"
[[ -n "$CHRLEN" ]] || fail "contig $CONTIG not in ${FA}.fai"

echo "FA=$FA"
echo "BAM=$BAM"
echo "CONTIG=$CONTIG length=$CHRLEN chunk=$CHUNK mode=$MODE"
echo ""

failures=0
start=1
while [[ "$start" -le "$CHRLEN" ]]; do
  end=$(( start + CHUNK - 1 ))
  [[ "$end" -gt "$CHRLEN" ]] && end="$CHRLEN"
  region="${CONTIG}:${start}-${end}"

  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "WINDOW ${start}-${end} (${region})"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  if ! bash "$ROOT/scripts/compare_vcf_parity.sh" "$region" "$MODE"; then
    echo "FAIL: compare_vcf_parity.sh $region" >&2
    failures=$(( failures + 1 ))
  fi
  if ! bash "$ROOT/scripts/compare_phased_vcf.sh" "$region" "$MODE"; then
    echo "FAIL: compare_phased_vcf.sh $region" >&2
    failures=$(( failures + 1 ))
  fi

  start=$(( end + 1 ))
done

echo ""
if [[ "$failures" -eq 0 ]]; then
  echo "OK: all CHM13 chr20 windows passed (${MODE})"
  exit 0
fi
echo "error: $failures parity step(s) failed" >&2
exit 1
