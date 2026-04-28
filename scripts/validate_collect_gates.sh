#!/usr/bin/env bash
# Parity / regression gates for `collect-bam-variation` using fixture golden outputs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
HIFI_FA="${HIFI_FA:-$ROOT/test_data/chr11_2M.fa}"
HIFI_BAM="${HIFI_BAM:-$ROOT/test_data/HG002_chr11_hifi_test.bam}"
ONT_FA="${ONT_FA:-$ROOT/test_data/chr11_2M.fa}"
ONT_BAM="${ONT_BAM:-$ROOT/test_data/HG002_chr11_ont_test.bam}"
HIFI_EXPECTED_TSV="${HIFI_EXPECTED_TSV:-$ROOT/test_data/expected/hifi_collect_expected.tsv}"
ONT_EXPECTED_TSV="${ONT_EXPECTED_TSV:-$ROOT/test_data/expected/ont_collect_expected.tsv}"
HIFI_EXPECTED_VCF="${HIFI_EXPECTED_VCF:-$ROOT/test_data/expected/hifi_collect_expected.phased.vcf}"
ONT_EXPECTED_VCF="${ONT_EXPECTED_VCF:-$ROOT/test_data/expected/ont_collect_expected.phased.vcf}"

# Backward-compatible fallback for older mixed-case fixture paths.
if [[ ! -f "$ONT_FA" && -f "$ROOT/test_data/longCallD/chr11_2M.fa" ]]; then
  ONT_FA="$ROOT/test_data/longCallD/chr11_2M.fa"
fi
if [[ ! -f "$ONT_BAM" && -f "$ROOT/test_data/longCallD/HG002_chr11_ont_test.bam" ]]; then
  ONT_BAM="$ROOT/test_data/longCallD/HG002_chr11_ont_test.bam"
fi
TMPDIR="${TMPDIR:-/tmp}"

if [[ ! -x "$PGPHASE" ]]; then
  echo "error: $PGPHASE not found or not executable (build with: make)" >&2
  exit 1
fi

if [[ ! -f "$HIFI_FA" || ! -f "$HIFI_BAM" ]]; then
  echo "error: HiFi fixture missing (FA=$HIFI_FA BAM=$HIFI_BAM)" >&2
  exit 1
fi
if [[ ! -f "$HIFI_EXPECTED_TSV" ]]; then
  echo "error: HiFi expected TSV missing: $HIFI_EXPECTED_TSV" >&2
  exit 1
fi
if [[ ! -f "$HIFI_EXPECTED_VCF" ]]; then
  echo "error: HiFi expected phased VCF missing: $HIFI_EXPECTED_VCF" >&2
  exit 1
fi

hash_file() {
  local f="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum < "$f" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$f" | awk '{print $1}'
  else
    echo "error: no sha256 tool found (need sha256sum or shasum)" >&2
    exit 1
  fi
}

run_collect() {
  local threads="$1"
  local out_tsv="$2"
  local out_vcf="$3"
  local fa="$4"
  local bam="$5"
  local mode="$6"
  local log
  log="$(mktemp)"
  set +e
  "$PGPHASE" collect-bam-variation -t "$threads" "$mode" --include-filtered "$fa" "$bam" \
    -o "$out_tsv" --phased-vcf-output "$out_vcf" >"$log" 2>&1
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

t1_tsv="$TMPDIR/pgphase_validate_hifi_t1.tsv"
t1_vcf="$TMPDIR/pgphase_validate_hifi_t1.phased.vcf"
t4_tsv="$TMPDIR/pgphase_validate_hifi_t4.tsv"
t4_vcf="$TMPDIR/pgphase_validate_hifi_t4.phased.vcf"
count_t1="$(run_collect 1 "$t1_tsv" "$t1_vcf" "$HIFI_FA" "$HIFI_BAM" --hifi)"
count_t4="$(run_collect 4 "$t4_tsv" "$t4_vcf" "$HIFI_FA" "$HIFI_BAM" --hifi)"

hash_exp_hifi_tsv="$(hash_file "$HIFI_EXPECTED_TSV")"
hash_t1_tsv="$(hash_file "$t1_tsv")"
hash_t4_tsv="$(hash_file "$t4_tsv")"
hash_exp_hifi_vcf="$(hash_file "$HIFI_EXPECTED_VCF")"
hash_t1_vcf="$(hash_file "$t1_vcf")"
hash_t4_vcf="$(hash_file "$t4_vcf")"

if [[ "$hash_t1_tsv" != "$hash_exp_hifi_tsv" ]]; then
  echo "error: HiFi TSV golden mismatch (t1 sha256 $hash_t1_tsv expected $hash_exp_hifi_tsv)" >&2
  echo "       expected: $HIFI_EXPECTED_TSV" >&2
  echo "       actual:   $t1_tsv" >&2
  exit 1
fi
if [[ "$hash_t1_tsv" != "$hash_t4_tsv" ]]; then
  echo "error: deterministic TSV mismatch (t1 sha256 $hash_t1_tsv vs t4 $hash_t4_tsv)" >&2
  exit 1
fi

if [[ "$hash_t1_vcf" != "$hash_exp_hifi_vcf" ]]; then
  echo "error: HiFi phased VCF golden mismatch (t1 sha256 $hash_t1_vcf expected $hash_exp_hifi_vcf)" >&2
  echo "       expected: $HIFI_EXPECTED_VCF" >&2
  echo "       actual:   $t1_vcf" >&2
  exit 1
fi
if [[ "$hash_t1_vcf" != "$hash_t4_vcf" ]]; then
  echo "error: deterministic phased VCF mismatch (t1 sha256 $hash_t1_vcf vs t4 $hash_t4_vcf)" >&2
  exit 1
fi

echo "ok: HiFi TSV golden match ($hash_t1_tsv), t1 vs t4 deterministic; sites=$count_t1"
echo "ok: HiFi phased VCF golden match ($hash_t1_vcf), t1 vs t4 deterministic"

if [[ -f "$ONT_FA" && -f "$ONT_BAM" ]]; then
  ont_tsv="$TMPDIR/pgphase_validate_ont.tsv"
  ont_vcf="$TMPDIR/pgphase_validate_ont.phased.vcf"
  count_ont="$(run_collect 1 "$ont_tsv" "$ont_vcf" "$ONT_FA" "$ONT_BAM" --ont)"
  if [[ -f "$ONT_EXPECTED_TSV" ]]; then
    hash_exp_ont_tsv="$(hash_file "$ONT_EXPECTED_TSV")"
    hash_ont_tsv="$(hash_file "$ont_tsv")"
    if [[ "$hash_ont_tsv" != "$hash_exp_ont_tsv" ]]; then
      echo "error: ONT TSV golden mismatch (sha256 $hash_ont_tsv expected $hash_exp_ont_tsv)" >&2
      echo "       expected: $ONT_EXPECTED_TSV" >&2
      echo "       actual:   $ont_tsv" >&2
      exit 1
    fi
    echo "ok: ONT TSV golden match ($hash_ont_tsv); sites=$count_ont"
  else
    echo "ok: ONT fixture ran ($count_ont sites) -> $ont_tsv"
    echo "warn: ONT expected TSV not found: $ONT_EXPECTED_TSV (golden check skipped)"
  fi
  if [[ -f "$ONT_EXPECTED_VCF" ]]; then
    hash_exp_ont_vcf="$(hash_file "$ONT_EXPECTED_VCF")"
    hash_ont_vcf="$(hash_file "$ont_vcf")"
    if [[ "$hash_ont_vcf" != "$hash_exp_ont_vcf" ]]; then
      echo "error: ONT phased VCF golden mismatch (sha256 $hash_ont_vcf expected $hash_exp_ont_vcf)" >&2
      echo "       expected: $ONT_EXPECTED_VCF" >&2
      echo "       actual:   $ont_vcf" >&2
      exit 1
    fi
    echo "ok: ONT phased VCF golden match ($hash_ont_vcf)"
  else
    echo "warn: ONT expected phased VCF not found: $ONT_EXPECTED_VCF (golden check skipped)"
  fi
else
  echo "skip: ONT fixture not present (FA=$ONT_FA BAM=$ONT_BAM)"
fi
