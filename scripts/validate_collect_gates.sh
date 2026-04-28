#!/usr/bin/env bash
# Parity / regression gates for `collect-bam-variation` (HiFi candidate count + deterministic threads).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
HIFI_FA="${HIFI_FA:-$ROOT/test_data/chr11_2M.fa}"
HIFI_BAM="${HIFI_BAM:-$ROOT/test_data/HG002_chr11_hifi_test.bam}"
ONT_FA="${ONT_FA:-$ROOT/test_data/chr11_2M.fa}"
ONT_BAM="${ONT_BAM:-$ROOT/test_data/HG002_chr11_ont_test.bam}"

# Backward-compatible fallback for older mixed-case fixture paths.
if [[ ! -f "$ONT_FA" && -f "$ROOT/test_data/longCallD/chr11_2M.fa" ]]; then
  ONT_FA="$ROOT/test_data/longCallD/chr11_2M.fa"
fi
if [[ ! -f "$ONT_BAM" && -f "$ROOT/test_data/longCallD/HG002_chr11_ont_test.bam" ]]; then
  ONT_BAM="$ROOT/test_data/longCallD/HG002_chr11_ont_test.bam"
fi
EXPECTED_HIFI="${EXPECTED_HIFI:-8745}"
TMPDIR="${TMPDIR:-/tmp}"

if [[ ! -x "$PGPHASE" ]]; then
  echo "error: $PGPHASE not found or not executable (build with: make)" >&2
  exit 1
fi

if [[ ! -f "$HIFI_FA" || ! -f "$HIFI_BAM" ]]; then
  echo "error: HiFi fixture missing (FA=$HIFI_FA BAM=$HIFI_BAM)" >&2
  exit 1
fi

run_collect() {
  local threads="$1"
  local out="$2"
  local fa="$3"
  local bam="$4"
  local log
  log="$(mktemp)"
  set +e
  "$PGPHASE" collect-bam-variation -t "$threads" --include-filtered "$fa" "$bam" -o "$out" >"$log" 2>&1
  local st=$?
  set -e
  if [[ "$st" -ne 0 ]]; then
    cat "$log" >&2
    rm -f "$log"
    exit "$st"
  fi
  grep -oE 'Collected [0-9]+ candidate variant sites' "$log" | grep -oE '[0-9]+' | head -1
  rm -f "$log"
}

t1="$TMPDIR/pgphase_validate_hifi_t1.tsv"
t4="$TMPDIR/pgphase_validate_hifi_t4.tsv"
count_t1="$(run_collect 1 "$t1" "$HIFI_FA" "$HIFI_BAM")"
count_t4="$(run_collect 4 "$t4" "$HIFI_FA" "$HIFI_BAM")"

if [[ "$count_t1" != "$EXPECTED_HIFI" ]]; then
  echo "error: HiFi candidate count t1=$count_t1 expected $EXPECTED_HIFI" >&2
  exit 1
fi
if [[ "$count_t4" != "$EXPECTED_HIFI" ]]; then
  echo "error: HiFi candidate count t4=$count_t4 expected $EXPECTED_HIFI" >&2
  exit 1
fi

hash_t1="$(sha256sum < "$t1" | awk '{print $1}')"
hash_t4="$(sha256sum < "$t4" | awk '{print $1}')"
if [[ "$hash_t1" != "$hash_t4" ]]; then
  echo "error: deterministic output mismatch (t1 sha256 $hash_t1 vs t4 $hash_t4)" >&2
  exit 1
fi

echo "ok: HiFi Collected $count_t1 candidate sites; t1 vs t4 output identical ($hash_t1)"

if [[ -f "$ONT_FA" && -f "$ONT_BAM" ]]; then
  ont_out="$TMPDIR/pgphase_validate_ont.tsv"
  count_ont="$(run_collect 1 "$ont_out" "$ONT_FA" "$ONT_BAM")"
  echo "ok: ONT fixture ran ($count_ont sites) -> $ont_out"
else
  echo "skip: ONT fixture not present (FA=$ONT_FA BAM=$ONT_BAM)"
fi
