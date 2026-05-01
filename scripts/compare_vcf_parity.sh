#!/usr/bin/env bash
# compare_vcf_parity.sh — compare pgPhase phased VCF against longcallD call VCF.
# Rows keyed by CHROM, POS, REF, ALT; payload includes a derived variant "VCIGAR"
# (allele-shape digest: X/I/D for concrete REF/ALT, or S:… for symbolic ALT using INFO),
# plus FILTER / INFO subset / GT:PS.
#
# Usage:
#   ./scripts/compare_vcf_parity.sh [REGION] [--hifi|--ont]
#
# Defaults:
#   REGION  CHM13#0#chr20:15000001-15500000 (HG002 HiFi chr20 slice; FA+BAM under test_data/chr20_quick/)
#   mode    --hifi
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
LONGCALLD_BIN="${LONGCALLD:-$ROOT/../longcallD/bin/longcallD}"
REGION="${1:-CHM13#0#chr20:15000001-15500000}"
MODE="${2:---hifi}"

if [[ -z "${FA:-}" ]]; then
    if [[ "$MODE" == "--ont" ]]; then
        FA="$ROOT/test_data/chr11_2M.fa"
    else
        FA="$ROOT/test_data/chr20_quick/CHM13_chr20_only.fa"
    fi
fi

if [[ -z "${BAM:-}" ]]; then
    if [[ "$MODE" == "--ont" ]]; then
        BAM="$ROOT/test_data/HG002_chr11_ont_test.bam"
    else
        BAM="$ROOT/test_data/chr20_quick/HG002_CHM13_chr20_15000001_15500000.bam"
    fi
fi

fail() { echo "error: $*" >&2; exit 2; }
[[ -x "$PGPHASE" ]] || fail "pgphase not found: $PGPHASE"
[[ -x "$LONGCALLD_BIN" ]] || fail "longcallD not found: $LONGCALLD_BIN"
[[ -f "$FA" ]] || fail "reference FASTA not found: $FA"
[[ -f "$BAM" ]] || fail "BAM not found: $BAM"

TMP="$(mktemp -d)"
if [[ "${KEEP_TMP:-0}" == "1" ]]; then
    echo "temp dir: $TMP"
else
    trap 'rm -rf "$TMP"' EXIT
fi

echo "region : $REGION"
echo "mode   : $MODE"
echo "BAM    : $BAM"
echo ""

"$PGPHASE" collect-bam-variation -t 1 --include-filtered \
    "$MODE" \
    --phased-vcf-output "$TMP/pg.vcf" \
    "$FA" "$BAM" "$REGION" -o "$TMP/pg.tsv" >/dev/null 2>&1

"$LONGCALLD_BIN" call -t 1 "$MODE" \
    "$FA" "$BAM" "$REGION" >"$TMP/lcd.vcf" 2>/dev/null

awk '!/^#/ {print $1"\t"$2"\t"$4"\t"$5}' "$TMP/pg.vcf" | sort -u > "$TMP/pg.keys"
awk '!/^#/ {print $1"\t"$2"\t"$4"\t"$5}' "$TMP/lcd.vcf" | sort -u > "$TMP/lcd.keys"
comm -23 "$TMP/pg.keys" "$TMP/lcd.keys" > "$TMP/only_pg.keys"
comm -13 "$TMP/pg.keys" "$TMP/lcd.keys" > "$TMP/only_lcd.keys"

n_only_pg="$(wc -l < "$TMP/only_pg.keys")"
n_only_lcd="$(wc -l < "$TMP/only_lcd.keys")"
echo "pg keys : $(wc -l < "$TMP/pg.keys")"
echo "lcd keys: $(wc -l < "$TMP/lcd.keys")"
echo "only_pg : $n_only_pg"
echo "only_lcd: $n_only_lcd"

if [[ "$n_only_pg" -gt 0 ]]; then
    echo ""
    echo "ONLY_PG examples:"
    head -20 "$TMP/only_pg.keys"
fi
if [[ "$n_only_lcd" -gt 0 ]]; then
    echo ""
    echo "ONLY_LCD examples:"
    head -20 "$TMP/only_lcd.keys"
fi

extract_cmp_table() {
    local in_vcf="$1"
    local out_tsv="$2"
    awk '
        function info_get(info, key,    n,i,kv,a) {
            n = split(info, kv, ";");
            for (i = 1; i <= n; ++i) {
                split(kv[i], a, "=");
                if (a[1] == key) return (length(a[2]) ? a[2] : ".");
            }
            return ".";
        }
        # Allele-shape token comparable between callers (VCF has no SAM CIGAR in INFO).
        # Concrete REF/ALT: X=rslen==altlen (SNV/MNV), I/D from lengths.
        # Symbolic ALT: S:SVTYPE:SVLEN:END:rflen using INFO + REF column length as fallback.
        function variant_cigar(pos_s, ref, alt, info,    pos, rlen, alen, op, endv, svlen, st, sym) {
            pos = pos_s + 0;
            sym = (alt ~ /^<[^>]+>$/);
            if (!sym) {
                rlen = length(ref);
                alen = length(alt);
                if (rlen > alen) op = "D";
                else if (rlen < alen) op = "I";
                else op = "X";
                return op ":" rlen ":" alen;
            }
            endv = info_get(info, "END");
            svlen = info_get(info, "SVLEN");
            st = info_get(info, "SVTYPE");
            rlen = length(ref);
            if (endv != "." && endv != "" && endv ~ /^[0-9]+$/) {
                op = int(endv) - pos + 1;
                if (op >= 1) rlen = op;
            }
            return "S:" st ":" svlen ":" endv ":" length(ref);
        }
        BEGIN { OFS = "\t"; }
        /^#/ { next; }
        {
            chrom = $1;
            pos = $2;
            ref = $4;
            alt = $5;
            k = chrom "\t" pos "\t" ref "\t" alt;
            vcig = variant_cigar(pos, ref, alt, $8);
            gt = ".";
            ps = ".";
            nfmt = split($9, fmtk, ":");
            nval = split($10, fmtv, ":");
            for (i = 1; i <= nfmt; ++i) {
                if (fmtk[i] == "GT" && i <= nval && length(fmtv[i])) gt = fmtv[i];
                else if (fmtk[i] == "PS" && i <= nval && length(fmtv[i])) ps = fmtv[i];
            }
            clean = (index(";" $8 ";", ";CLEAN;") > 0 ? "1" : "0");
            print chrom, pos, ref, alt, vcig, $7, info_get($8, "END"), info_get($8, "SVTYPE"),
                  info_get($8, "SVLEN"), clean, gt, ps;
        }
    ' "$in_vcf" | sort -u > "$out_tsv"
}

extract_cmp_table "$TMP/pg.vcf" "$TMP/pg.cmp.tsv"
extract_cmp_table "$TMP/lcd.vcf" "$TMP/lcd.cmp.tsv"

set +e
awk -F'\t' '
    BEGIN { mism=0; }
    NR==FNR {
        k = $1 "\t" $2 "\t" $3 "\t" $4;
        lcd[k] = $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12;
        next;
    }
    {
        key = $1 "\t" $2 "\t" $3 "\t" $4;
        if (!(key in lcd)) next;
        pg = $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12;
        lc = lcd[key];
        if (pg != lc) {
            ++mism;
            if (mism <= 30) {
                print "DIFF key=" key " (CHROM POS REF ALT)";
                print "  pg : VCIGAR FILTER END SVTYPE SVLEN CLEAN GT PS = " pg;
                print "  lcd: VCIGAR FILTER END SVTYPE SVLEN CLEAN GT PS = " lc;
            }
        }
    }
    END {
        print "shared_diff: " mism;
        exit(mism > 0 ? 1 : 0);
    }
' "$TMP/lcd.cmp.tsv" "$TMP/pg.cmp.tsv"
shared_diff_exit=$?
set -e

if [[ "$n_only_pg" -gt 0 || "$n_only_lcd" -gt 0 ]]; then
    exit 1
fi
if [[ "$shared_diff_exit" -ne 0 ]]; then
    exit 1
fi
echo ""
echo "OK: VCF parity matched (POS REF ALT + VCIGAR + FILTER/INFO subset + GT:PS)"
