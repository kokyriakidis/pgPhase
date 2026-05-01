#!/usr/bin/env bash
# Integration test: pgPhase TSV vs longcallD stderr (verbose) on HG002 HiFi chr20 slice.
# Requires: built ./pgphase, samtools (optional, for sanity), longcallD binary.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
LONGCALLD_BIN="${LONGCALLD:-$ROOT/../longcallD/bin/longcallD}"
FA="$ROOT/test_data/chr20_quick/CHM13_chr20_only.fa"
BAM="$ROOT/test_data/chr20_quick/HG002_CHM13_chr20_15000001_15500000.bam"
REGION="CHM13#0#chr20:15000001-15500000"

if [[ ! -x "$PGPHASE" ]]; then
  echo "error: pgphase not found at $PGPHASE (run: make)" >&2
  exit 1
fi
if [[ ! -x "$LONGCALLD_BIN" ]]; then
  echo "skip: longcallD not executable at $LONGCALLD_BIN; set LONGCALLD to your binary" >&2
  exit 0
fi
if [[ ! -f "$FA" || ! -f "$BAM" ]]; then
  echo "error: missing fixture FA=$FA BAM=$BAM" >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

"$PGPHASE" collect-bam-variation -t 1 --include-filtered "$FA" "$BAM" "$REGION" -o "$TMP/pg.tsv" \
  >/dev/null
"$LONGCALLD_BIN" call -t 2 -V 2 --hifi -o "$TMP/lcd.vcf" "$FA" "$BAM" "$REGION" 2>"$TMP/lcd.log"

# final: longcallD's post-compaction CandVarCate block (small set, strict parity on overlaps).
python3 "$ROOT/scripts/compare_candidates.py" \
  --pgphase "$TMP/pg.tsv" \
  --longcalld "$TMP/lcd.log" \
  --chrom chr11 \
  --category-stage final

echo "ok: compare_candidates.py (--category-stage final) passed on $REGION"

python3 "$ROOT/scripts/compare_candidates.py" \
  --pgphase "$TMP/pg.tsv" \
  --longcalld "$TMP/lcd.log" \
  --chrom 'CHM13#0#chr20' \
  --category-stage initial

echo "ok: compare_candidates.py (--category-stage initial, INIT_CAT) passed on $REGION"
