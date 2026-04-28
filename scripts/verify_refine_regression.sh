#!/usr/bin/env bash
# verify_refine_regression.sh
#
# Fast local regression check for phased BAM output parity behavior WITHOUT longcallD.
# It runs pgphase on the bundled HiFi fixture in:
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
FA="${FA:-$ROOT/test_data/chr11_2M.fa}"
BAM="${BAM:-$ROOT/test_data/HG002_chr11_hifi_test.bam}"
REGION="${REGION:-chr11:1255000-1260000}"

# Known-good hashes from strict parity baseline (normal + refine).
EXPECTED_NORMAL_SHA="c8c45960d884029887bd85b8ffbb0836af007b3178d2b2feee66243b1cbc5811"
EXPECTED_REFINE_SHA="02c217b128897cf264d9a6fb66e6d8014a2fac46285645959bfd3e9f97154884"

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
