#!/usr/bin/env bash
#SBATCH --job-name=cat_v21_allchr
#SBATCH --partition=batch_cpu
#SBATCH --account=dsricg
#SBATCH --qos=1d
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --array=1-24
#SBATCH --output=/sc1/groups/sbx/workspace/kyriakik/pgphase_eval/hg003_vg/logs/cat_v21_allchr_%A_%a.out
#SBATCH --error=/sc1/groups/sbx/workspace/kyriakik/pgphase_eval/hg003_vg/logs/cat_v21_allchr_%A_%a.out
#
# Build per-chromosome snarl-catalog sites VCFs from the HPRC v2.1 GBZ.
# One SLURM array task per chromosome (chr1..22, X, Y) -- each is an independent
# job with its own resources, so they run in parallel and a failure on one
# chromosome does not affect the others.
#
# Each task runs `pgphase build-snarl-catalog --contig <chr>` directly on the
# full GBZ (vg deconstruct -p GRCh38#0#<chr> extracts the contig itself), writes
# to node-local scratch, then stages the small result to $HOME (the shared
# /sc1/groups quota is volatile and truncates long writes), and best-effort
# mirrors to the shared FS. A striped copy (GRCh38#0#<chr> -> <chr>) is also
# produced for the linear-reference phasing step.
#
# Resubmit a single chromosome after a failure, e.g.:
#   sbatch --array=1 build_catalog_allchr_v21.sbatch     # chr1 only
#   sbatch --array=23 build_catalog_allchr_v21.sbatch    # chrX only

set -euo pipefail

# --- chromosome for this array task --------------------------------------
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "$IDX" -lt 0 ] || [ "$IDX" -ge "${#CHROMS[@]}" ]; then
  echo "ERROR: array task $SLURM_ARRAY_TASK_ID out of range 1..${#CHROMS[@]}" >&2
  exit 1
fi
CHR="${CHROMS[$IDX]}"
echo "=== array task $SLURM_ARRAY_TASK_ID -> $CHR ==="

# --- paths ----------------------------------------------------------------
B=/sc1/groups/sbx/workspace/kyriakik/pgphase_eval
export MAMBA_ROOT_PREFIX="$B/micromamba"
MM="$B/bin/micromamba"; ENV="$B/env/dv_eval"
VG="$B/bin/vg"; PGPHASE="$B/pgphase/pgphase"
D="$B/hg003_vg"; G="$D/graph_v21"; OUT="$D/catalog_v21"
mkdir -p "$OUT"
GBZ="$G/hprc-v2.1-mc-grch38.gbz"
SNARLS="$G/hprc-v2.1-mc-grch38.snarls"

SITES="$OUT/${CHR}.sites.vcf.gz"
STRIPED="$OUT/${CHR}.sites.striped.vcf.gz"

# --- node-local scratch (avoid long writes to contended Lustre) -----------
SCRATCH="${TMPDIR:-/tmp}/cat_v21_${CHR}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH"
SITES_LOCAL="$SCRATCH/${CHR}.sites.vcf.gz"
STRIPED_LOCAL="$SCRATCH/${CHR}.sites.striped.vcf.gz"
trap 'rm -rf "$SCRATCH"' EXIT

# --- $HOME staging (separate, uncontended quota) --------------------------
HOME_STAGE="$HOME/pgphase_catalog_v21"
mkdir -p "$HOME_STAGE"
SITES_HOME="$HOME_STAGE/${CHR}.sites.vcf.gz"
STRIPED_HOME="$HOME_STAGE/${CHR}.sites.striped.vcf.gz"

# build-snarl-catalog shells out to vg deconstruct/bgzip/tabix; put them on PATH.
export PATH="$(dirname "$VG"):$PATH"
ln -sf "$VG" "$(dirname "$VG")/vg" 2>/dev/null || true
run() { "$MM" run -p "$ENV" bash -c "$1"; }

# --- build sites VCF ------------------------------------------------------
if [ -s "$SITES_HOME" ]; then
  echo "sites VCF already staged in HOME, skipping build: $SITES_HOME"
else
  echo "=== build-snarl-catalog (GRCh38 $CHR, v2.1 snarls) -> node-local scratch ==="
  run "set -e; export PATH=$(dirname "$VG"):\$PATH; \
       $PGPHASE build-snarl-catalog \
         --ref-sample GRCh38 --contig $CHR \
         --snarls $SNARLS --threads 16 \
         --vg-bin $VG \
         -o $SITES_LOCAL \
         $GBZ"
  echo "=== stage sites: scratch -> \$HOME (quota-safe) ==="
  cp -f "$SITES_LOCAL" "$SITES_HOME"
  echo "staged in HOME: $(ls -lh "$SITES_HOME" | awk '{print $5, $9}')"
fi

# Best-effort mirror to the shared FS (may fail under quota; not fatal).
if cp -f "$SITES_HOME" "$SITES" 2>/dev/null; then
  echo "mirrored to shared FS: $SITES"
else
  echo "WARN: shared-FS copy blocked by quota; using HOME copy." >&2
  rm -f "$SITES" 2>/dev/null || true
fi

# --- validate -------------------------------------------------------------
echo "=== sites VCF ($CHR) ==="; ls -lh "$SITES_HOME"* 2>&1
NREC=$(run "zcat $SITES_HOME | grep -vc '^#'" || echo 0)
echo "records: $NREC"
run "echo ---first3---; zcat $SITES_HOME | grep -v ^## | head -3 | cut -f1-5" || true
# Guard: vg deconstruct can be OOM-killed mid-stream leaving a header-only VCF
# while bgzip still exits 0. Fail loudly so an empty catalog is never used.
if [ "$NREC" -lt 1000 ]; then
  echo "ERROR: $CHR catalog has only $NREC records (expected many thousands)." \
       "Likely OOM in vg deconstruct -- resubmit this chromosome with more --mem." >&2
  exit 1
fi
echo "CATALOG_V21_${CHR}_DONE"

# --- striped copy (GRCh38#0#<chr> -> <chr>) -------------------------------
if [ -s "$STRIPED_HOME" ]; then
  echo "striped VCF already staged in HOME, skipping: $STRIPED_HOME"
else
  echo "=== stripe contig name -> $CHR (scratch) ==="
  run "set -e; zcat $SITES_HOME | sed -e 's/GRCh38#0#${CHR}\b/${CHR}/g' | bgzip > $STRIPED_LOCAL; tabix -p vcf $STRIPED_LOCAL"
  echo "=== stage striped -> \$HOME ==="
  cp -f "$STRIPED_LOCAL" "$STRIPED_HOME"
  cp -f "$STRIPED_LOCAL.tbi" "$STRIPED_HOME.tbi"
  echo "staged striped in HOME: $(ls -lh "$STRIPED_HOME" | awk '{print $5, $9}')"
fi

# Best-effort mirror striped to shared FS.
if cp -f "$STRIPED_HOME" "$STRIPED" 2>/dev/null && cp -f "$STRIPED_HOME.tbi" "$STRIPED.tbi" 2>/dev/null; then
  echo "mirrored striped to shared FS: $STRIPED"
else
  echo "WARN: shared-FS striped copy blocked by quota; HOME copy at $STRIPED_HOME" >&2
  rm -f "$STRIPED" "$STRIPED.tbi" 2>/dev/null || true
fi
echo "=== striped sites (HOME, $CHR) ==="; ls -lh "$STRIPED_HOME"* 2>&1
echo "CATALOG_V21_${CHR}_STRIPED_DONE"
