#!/usr/bin/env bash
# verify_refine_regression.sh
#
# Fast local regression check for phased BAM output parity behavior WITHOUT longcallD.
# It runs pgphase on the bundled HG002 HiFi chr20 slice (CHM13) in:
#   1) normal phased BAM mode
#   2) --refine-aln phased BAM mode
# and compares sorted SAM SHA-256 hashes against known-good expected values.
#
# Usage:
#   ./scripts/verify_refine_regression.sh
#
# Exit codes:
#   0  pass
#   1  mismatch (regression)
#   2  setup/runtime error

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
FA="${FA:-$ROOT/test_data/chr20_quick/CHM13_chr20_only.fa}"
BAM="${BAM:-$ROOT/test_data/chr20_quick/HG002_CHM13_chr20_15000001_15500000.bam}"
REGION="${REGION:-CHM13#0#chr20:15200000-15205000}"

# Known-good hashes from strict parity baseline (normal + refine).
EXPECTED_NORMAL_SHA="27761c029a62018d00f4c9c7d045db958362b57fce9ffd1409f3c41000bbee66"
EXPECTED_REFINE_SHA="e45ecbd1b3d5a2e7d31b1929f2941cd18aaa28c650acba23e5f6bfa9568a742a"

fail_setup() { echo "error: $*" >&2; exit 2; }

command -v samtools >/dev/null 2>&1 || fail_setup "samtools not found in PATH"
command -v sha256sum >/dev/null 2>&1 || fail_setup "sha256sum not found in PATH"
[[ -x "$PGPHASE" ]] || fail_setup "pgphase binary not found: $PGPHASE (run: make)"
[[ -f "$FA" ]] || fail_setup "reference FASTA not found: $FA"
[[ -f "$BAM" ]] || fail_setup "input BAM not found: $BAM"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "Running regression fixture:"
echo "  region: $REGION"
echo "  fasta : $FA"
echo "  bam   : $BAM"
echo

"$PGPHASE" collect-bam-variation -t 1 --hifi --include-filtered \
  -b "$TMP/pgphase.bam" \
  "$FA" "$BAM" "$REGION" -o "$TMP/pgphase.tsv" >/dev/null

"$PGPHASE" collect-bam-variation -t 1 --hifi --include-filtered --refine-aln \
  -b "$TMP/pgphase_refine.bam" \
  "$FA" "$BAM" "$REGION" -o "$TMP/pgphase_refine.tsv" >/dev/null

samtools view "$TMP/pgphase.bam" | LC_ALL=C sort > "$TMP/pgphase.sam.sorted"
samtools view "$TMP/pgphase_refine.bam" | LC_ALL=C sort > "$TMP/pgphase_refine.sam.sorted"

NORMAL_SHA="$(sha256sum "$TMP/pgphase.sam.sorted" | awk '{print $1}')"
REFINE_SHA="$(sha256sum "$TMP/pgphase_refine.sam.sorted" | awk '{print $1}')"

echo "normal sha : $NORMAL_SHA"
echo "refine sha : $REFINE_SHA"
echo

ok=1
if [[ "$NORMAL_SHA" != "$EXPECTED_NORMAL_SHA" ]]; then
  echo "MISMATCH normal: expected $EXPECTED_NORMAL_SHA"
  if [[ -f "$ROOT/test_data/parity_compare/pgphase.new4.sam.sorted" ]]; then
    diff -u "$ROOT/test_data/parity_compare/pgphase.new4.sam.sorted" "$TMP/pgphase.sam.sorted" \
      > "$TMP/diff.normal.txt" || true
    echo "--- Offending normal lines (first 200 diff lines) ---"
    awk 'NR<=200 {print}' "$TMP/diff.normal.txt"
    echo "--- end normal diff ---"
  else
    echo "note: baseline file missing: test_data/parity_compare/pgphase.new4.sam.sorted"
  fi
  ok=0
fi
if [[ "$REFINE_SHA" != "$EXPECTED_REFINE_SHA" ]]; then
  echo "MISMATCH refine: expected $EXPECTED_REFINE_SHA"
  if [[ -f "$ROOT/test_data/parity_compare/pgphase_refine.new4.sam.sorted" ]]; then
    diff -u "$ROOT/test_data/parity_compare/pgphase_refine.new4.sam.sorted" "$TMP/pgphase_refine.sam.sorted" \
      > "$TMP/diff.refine.txt" || true
    echo "--- Offending refine lines (first 200 diff lines) ---"
    awk 'NR<=200 {print}' "$TMP/diff.refine.txt"
    echo "--- end refine diff ---"
  else
    echo "note: baseline file missing: test_data/parity_compare/pgphase_refine.new4.sam.sorted"
  fi
  ok=0
fi

if [[ "$ok" -ne 1 ]]; then
  echo "FAIL: regression detected."
  echo "Diff artifacts are under: $TMP"
  trap - EXIT
  exit 1
fi

echo "OK: regression check passed."
