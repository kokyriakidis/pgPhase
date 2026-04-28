#!/usr/bin/env bash
# compare_phased_vcf.sh — compare pgPhase --phased-vcf-output against longcallD's phased VCF.
#
# Usage:
#   ./scripts/compare_phased_vcf.sh [REGION] [--hifi|--ont]
#
# Defaults:
#   REGION  chr11:1230000-1260000
#   mode    --hifi  (HiFi BAM)
#
# Environment overrides:
#   PGPHASE     path to pgphase binary    (default: ./pgphase)
#   LONGCALLD   path to longcallD binary  (default: ../longcallD/bin/longcallD)
#   FA          path to reference FASTA   (default: test_data/chr11_2M.fa)
#   BAM         path to BAM               (default: derived from mode)
#   KEEP_TMP    set to 1 to keep temp dir for inspection
#
# Exit codes:
#   0  all phased het sites match GT and PS between both tools
#   1  one or more GT/PS mismatches detected
#   2  setup error (missing binary / file)
#
# Duplicate POS: multiple phased hets can share the same VCF POS (e.g. nested DELs).
# Keys use POS + SVLEN from INFO so both records are compared.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
LONGCALLD_BIN="${LONGCALLD:-$ROOT/../longcallD/bin/longcallD}"
FA="${FA:-$ROOT/test_data/chr11_2M.fa}"
REGION="${1:-chr11:1230000-1260000}"
MODE="${2:---hifi}"

# Derive BAM from mode if not set explicitly
if [[ -z "${BAM:-}" ]]; then
    if [[ "$MODE" == "--ont" ]]; then
        BAM="$ROOT/test_data/HG002_chr11_ont_test.bam"
    else
        BAM="$ROOT/test_data/HG002_chr11_hifi_test.bam"
    fi
fi

# ── Validation ───────────────────────────────────────────────────────────────

fail() { echo "error: $*" >&2; exit 2; }

[[ -x "$PGPHASE" ]]      || fail "pgphase not found at $PGPHASE (run: make)"
[[ -x "$LONGCALLD_BIN" ]] || fail "longcallD not found at $LONGCALLD_BIN"
[[ -f "$FA" ]]            || fail "reference FASTA not found: $FA"
[[ -f "$BAM" ]]           || fail "BAM not found: $BAM"

TMP="$(mktemp -d)"
if [[ "${KEEP_TMP:-0}" == "1" ]]; then
    echo "temp dir: $TMP"
else
    trap 'rm -rf "$TMP"' EXIT
fi

# ── Run both tools ────────────────────────────────────────────────────────────

echo "region : $REGION"
echo "mode   : $MODE"
echo "BAM    : $BAM"
echo ""

"$PGPHASE" collect-bam-variation -t 1 --include-filtered \
    ${MODE:---hifi} \
    --phased-vcf-output "$TMP/pg_phased.vcf" \
    "$FA" "$BAM" "$REGION" -o "$TMP/pg.tsv" 2>/dev/null

"$LONGCALLD_BIN" call -t 1 "$MODE" \
    "$FA" "$BAM" "$REGION" > "$TMP/lcd.vcf" 2>/dev/null

# ── Extract phased het sites: POS SVLEN_key GT PS ─────────────────────────────

# SVLEN from INFO (signed for DEL). If missing (SNPs / small indels),
# use ALT to disambiguate duplicate POS records (e.g. multi-allelic INS at same POS).
extract_pg() {
    grep -v '^#' "$1" | awk '
        function svkey(info, alt) {
            if (match(info, /SVLEN=-?[0-9]+/))
                return substr(info, RSTART+6, RLENGTH-6);
            return "ALT=" alt;
        }
        $NF ~ /[|]/ {
            split($NF, f, ":");
            gt = f[1]; ps = f[length(f)];
            if (gt == "1|1") next;
            print $2, svkey($8, $5), gt, ps
        }
    '
}

extract_lcd() {
    grep -v '^#' "$1" | awk '
        function svkey(info, alt) {
            if (match(info, /SVLEN=-?[0-9]+/))
                return substr(info, RSTART+6, RLENGTH-6);
            return "ALT=" alt;
        }
        $NF ~ /[|]/ {
            split($NF, f, ":");
            gt = f[1];
            if (gt == "1|1") next;
            ps = f[length(f)];
            print $2, svkey($8, $5), gt, ps
        }
    '
}

extract_pg  "$TMP/pg_phased.vcf" | sort -k1,1n -k2,2 > "$TMP/pg_het.tsv"
extract_lcd "$TMP/lcd.vcf"       | sort -k1,1n -k2,2 > "$TMP/lcd_het.tsv"

pg_n=$(wc -l  < "$TMP/pg_het.tsv")
lcd_n=$(wc -l < "$TMP/lcd_het.tsv")
echo "pgPhase phased het sites : $pg_n"
echo "longcallD phased het sites: $lcd_n"
echo ""

# ── Compare shared sites ──────────────────────────────────────────────────────

n_match=0
n_diff=0
n_only_pg=0
n_only_lcd=0

# Build lookup maps: key = POS_SVLEN (handles duplicate POS)
declare -A pg_map
declare -A lcd_map
while read -r pos svlen gt ps; do
    pg_map["${pos}_${svlen}"]="$gt:$ps"
done < "$TMP/pg_het.tsv"

while read -r pos svlen gt ps; do
    lcd_map["${pos}_${svlen}"]="$gt:$ps"
done < "$TMP/lcd_het.tsv"

# Check all pgPhase sites against longcallD
for k in "${!pg_map[@]}"; do
    pg_val="${pg_map[$k]}"
    if [[ -n "${lcd_map[$k]+set}" ]]; then
        lcd_val="${lcd_map[$k]}"
        if [[ "$pg_val" == "$lcd_val" ]]; then
            n_match=$(( n_match + 1 ))
        else
            echo "DIFF     key=$k  pg=$pg_val  lcd=$lcd_val"
            n_diff=$(( n_diff + 1 ))
        fi
    else
        echo "ONLY_PG  key=$k  pg=$pg_val"
        n_only_pg=$(( n_only_pg + 1 ))
    fi
done

for k in "${!lcd_map[@]}"; do
    if [[ -z "${pg_map[$k]+set}" ]]; then
        lcd_val="${lcd_map[$k]}"
        echo "ONLY_LCD key=$k  lcd=$lcd_val"
        n_only_lcd=$(( n_only_lcd + 1 ))
    fi
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "─────────────────────────────────────"
echo "MATCH    : $n_match"
echo "DIFF     : $n_diff     (GT or PS mismatch — failures)"
echo "ONLY_PG  : $n_only_pg  (phased by pgPhase, not by longcallD)"
echo "ONLY_LCD : $n_only_lcd (phased by longcallD, missing or different in pgPhase)"
echo "─────────────────────────────────────"

if (( n_diff > 0 )); then
    echo "FAIL: $n_diff GT/PS mismatch(es) detected"
    exit 1
fi

echo "OK: all shared het sites match GT and PS"
