#!/usr/bin/env bash
#
# run_eval.sh — Does pgphase read phasing help DeepVariant call variants?
#
# Runs DeepVariant in two arms on the same reads and benchmarks both against
# the GIAB HG002 v4.2.1 truth set with hap.py:
#
#   Arm A (dv_default): DeepVariant on the RAW BAM. DeepVariant runs its own
#                       internal read phasing (PACBIO default phase_reads=true),
#                       then calls variants. This is the stock DeepVariant
#                       PacBio result and the baseline.
#
#   Arm B (dv_our_hp):  DeepVariant on the pgphase PHASED BAM (carries HP/PS
#                       tags). DeepVariant's internal phasing is disabled
#                       (phase_reads=false) but it still sorts the pileup by
#                       haplotype, which makes make_examples parse the HP tag
#                       from OUR reads (sort_by_haplotypes=true). This isolates
#                       the effect of substituting pgphase phasing for
#                       DeepVariant's own.
#
# If arm B's SNP/INDEL F1 exceeds arm A's, pgphase phasing improves DeepVariant
# variant calling. The comparison is reported by summarize_happy.py.
#
# Flag basis (DeepVariant r1.10): from v1.10.0 the per-model make_examples args
# (incl. the PACBIO default phase_reads=true) live in the model's
# example_info.json, so arm A needs no flags. make_examples_options.py documents
# the HP-tag contract used by arm B: "HP is parsed if --phase_reads=False and
# --sort_by_haplotypes ... are set." parse_sam_aux_fields is deprecated in 1.10
# (AUX parsing is now auto-controlled by the requested channels/flags), so it is
# not set here. Arm B only overrides phase_reads=false,sort_by_haplotypes=true.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
readonly DV_VERSION="1.10.0"
readonly HAPPY_VERSION="v0.3.12"
readonly MODEL_TYPE="PACBIO"

RAW_BAM=""
OUR_BAM=""
OUTDIR=""
RESOURCE_DIR="${SCRIPT_DIR}/resources"
REGIONS="chr20"
THREADS="$(nproc)"
ENGINE=""

die() { echo "Error: $*" >&2; exit 1; }

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --raw-bam FILE      Aligned, indexed BAM WITHOUT pgphase HP tags (arm A input).
  --our-bam FILE      pgphase phased BAM WITH HP/PS tags (arm B input).
                      Both must be the same reads/alignment; only phasing differs.

Optional:
  --outdir DIR        Output directory [${SCRIPT_DIR}/results]
  --resource-dir DIR  Resources from fetch_resources.sh [${RESOURCE_DIR}]
  --regions STR       Region(s) to evaluate [${REGIONS}]
  --threads INT       Threads / shards [nproc]
  --engine NAME       singularity | apptainer | docker [auto-detected]
  -h, --help          Show this help
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-bam)      RAW_BAM="$2";      shift 2 ;;
        --our-bam)      OUR_BAM="$2";      shift 2 ;;
        --outdir)       OUTDIR="$2";       shift 2 ;;
        --resource-dir) RESOURCE_DIR="$2"; shift 2 ;;
        --regions)      REGIONS="$2";      shift 2 ;;
        --threads)      THREADS="$2";      shift 2 ;;
        --engine)       ENGINE="$2";       shift 2 ;;
        -h|--help)      usage 0 ;;
        *)              echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

[[ -n "${RAW_BAM}" ]] || { echo "Error: --raw-bam required" >&2; usage 1; }
[[ -n "${OUR_BAM}" ]] || { echo "Error: --our-bam required" >&2; usage 1; }
[[ -z "${OUTDIR}"  ]] && OUTDIR="${SCRIPT_DIR}/results"

for f in "${RAW_BAM}" "${OUR_BAM}"; do
    [[ -f "${f}" ]] || die "BAM not found: ${f}"
    [[ -f "${f}.bai" || -f "${f%.bam}.bai" || -f "${f}.csi" ]] \
        || die "BAM index not found for ${f} (run: samtools index ${f})"
done

readonly REF_FASTA="${RESOURCE_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
readonly TRUTH_VCF="${RESOURCE_DIR}/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
readonly TRUTH_BED="${RESOURCE_DIR}/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
for f in "${REF_FASTA}" "${TRUTH_VCF}" "${TRUTH_BED}"; do
    [[ -f "${f}" ]] || die "resource missing: ${f} (run fetch_resources.sh first)"
done

detect_engine() {
    if [[ -n "${ENGINE}" ]]; then echo "${ENGINE}"; return; fi
    if command -v singularity >/dev/null 2>&1; then echo singularity; return; fi
    if command -v apptainer   >/dev/null 2>&1; then echo apptainer;   return; fi
    if command -v docker      >/dev/null 2>&1; then echo docker;      return; fi
    die "no container engine found; pass --engine"
}
CONTAINER_ENGINE="$(detect_engine)"
readonly CONTAINER_ENGINE
readonly IMG_DIR="${RESOURCE_DIR}/images"
readonly DV_SIF="${IMG_DIR}/deepvariant_${DV_VERSION}.sif"
readonly HAPPY_SIF="${IMG_DIR}/happy_${HAPPY_VERSION}.sif"

mkdir -p "${OUTDIR}"

# run_in_container IMAGE_KIND -- command...
# Binds the host filesystem root so absolute paths resolve identically inside
# and outside the container (the trick the DeepVariant FAQ recommends). hap.py
# and DeepVariant both run this way; IMAGE_KIND selects which image.
run_in_container() {
    local kind="$1"; shift
    [[ "$1" == "--" ]] && shift
    case "${CONTAINER_ENGINE}" in
        singularity|apptainer)
            local sif
            case "${kind}" in
                dv)    sif="${DV_SIF}" ;;
                happy) sif="${HAPPY_SIF}" ;;
            esac
            "${CONTAINER_ENGINE}" exec --bind / --pwd "${PWD}" "${sif}" "$@"
            ;;
        docker)
            local img
            case "${kind}" in
                dv)    img="google/deepvariant:${DV_VERSION}" ;;
                happy) img="jmcdani20/hap.py:${HAPPY_VERSION}" ;;
            esac
            docker run --rm -v /:/ -w "${PWD}" \
                -u "$(id -u):$(id -g)" "${img}" "$@"
            ;;
    esac
}

run_deepvariant_arm() {
    # run_deepvariant_arm ARM_NAME INPUT_BAM [EXTRA_MAKE_EXAMPLES_ARGS]
    local arm="$1" bam="$2" extra_args="${3:-}"
    local arm_dir="${OUTDIR}/${arm}"
    local out_vcf="${arm_dir}/deepvariant.vcf.gz"
    mkdir -p "${arm_dir}"

    if [[ -s "${out_vcf}" ]]; then
        echo "[${arm}] DeepVariant VCF exists, skipping call."
        return
    fi

    echo "[${arm}] Running DeepVariant (${MODEL_TYPE}) on $(basename "${bam}") ..."
    [[ -n "${extra_args}" ]] && echo "[${arm}]   make_examples_extra_args=${extra_args}"

    local -a cmd=(
        /opt/deepvariant/bin/run_deepvariant
        --model_type "${MODEL_TYPE}"
        --ref "${REF_FASTA}"
        --reads "${bam}"
        --output_vcf "${out_vcf}"
        --output_gvcf "${arm_dir}/deepvariant.g.vcf.gz"
        --num_shards "${THREADS}"
        --regions "${REGIONS}"
        --intermediate_results_dir "${arm_dir}/intermediate"
    )
    [[ -n "${extra_args}" ]] && cmd+=(--make_examples_extra_args "${extra_args}")

    { time run_in_container dv -- "${cmd[@]}" ; } \
        2> >(tee "${arm_dir}/timing.log" >&2)
}

run_happy_arm() {
    # run_happy_arm ARM_NAME — benchmark that arm's VCF against GIAB truth.
    local arm="$1"
    local arm_dir="${OUTDIR}/${arm}"
    local query_vcf="${arm_dir}/deepvariant.vcf.gz"
    local happy_prefix="${arm_dir}/happy"

    if [[ -s "${happy_prefix}.summary.csv" ]]; then
        echo "[${arm}] hap.py summary exists, skipping."
        return
    fi
    echo "[${arm}] Benchmarking with hap.py (vcfeval engine) ..."
    run_in_container happy -- /opt/hap.py/bin/hap.py \
        --threads "${THREADS}" \
        -r "${REF_FASTA}" \
        -f "${TRUTH_BED}" \
        -o "${happy_prefix}" \
        --engine vcfeval \
        --pass-only \
        -l "${REGIONS}" \
        "${TRUTH_VCF}" \
        "${query_vcf}"
}

echo "=== Engine: ${CONTAINER_ENGINE} | Region: ${REGIONS} | Threads: ${THREADS} ==="

# Arm A: stock DeepVariant on the raw BAM (DeepVariant phases internally).
run_deepvariant_arm "dv_default" "${RAW_BAM}"

# Arm B: DeepVariant on the pgphase HP-tagged BAM; disable DV's own phasing so
# it uses OUR HP tags (sort_by_haplotypes stays on; add_hp_channel is default).
run_deepvariant_arm "dv_our_hp" "${OUR_BAM}" "phase_reads=false,sort_by_haplotypes=true"

run_happy_arm "dv_default"
run_happy_arm "dv_our_hp"

echo ""
echo "=== Summary: pgphase phasing effect on DeepVariant ==="
python3 "${SCRIPT_DIR}/summarize_happy.py" \
    --baseline "${OUTDIR}/dv_default/happy.summary.csv" \
    --treatment "${OUTDIR}/dv_our_hp/happy.summary.csv" \
    --baseline-label "DV default (self-phased)" \
    --treatment-label "DV + pgphase HP" \
    | tee "${OUTDIR}/comparison.txt"

echo ""
echo "Outputs in: ${OUTDIR}/"
echo "  dv_default/happy.summary.csv   (baseline)"
echo "  dv_our_hp/happy.summary.csv    (pgphase HP)"
echo "  comparison.txt                 (side-by-side)"
