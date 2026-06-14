#!/usr/bin/env bash
#
# fetch_resources.sh — Download everything the DeepVariant HP-tag evaluation
# needs on the cluster: GRCh38 reference, GIAB HG002 v4.2.1 truth set, and the
# DeepVariant + hap.py container images.
#
# Run once on the cluster before run_eval.sh. Re-running skips files that
# already exist, so it is safe to resume after an interrupted download.
#
# Requirements: curl, samtools, and a container runtime (singularity/apptainer
# or docker). Container engine is auto-detected; override with --engine.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
readonly DV_VERSION="1.6.1"
readonly HAPPY_VERSION="v0.3.12"

RESOURCE_DIR="${SCRIPT_DIR}/resources"
ENGINE=""

readonly REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
readonly GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
readonly TRUTH_VCF_NAME="HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
readonly TRUTH_BED_NAME="HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

die() { echo "Error: $*" >&2; exit 1; }

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --resource-dir DIR  Where to store downloads [${RESOURCE_DIR}]
  --engine NAME       Container engine: singularity | apptainer | docker
                      [auto-detected]
  -h, --help          Show this help

Downloads:
  reference/GRCh38_no_alt_analysis_set.fasta(.fai)
  benchmark/${TRUTH_VCF_NAME}(.tbi)
  benchmark/${TRUTH_BED_NAME}
  container images: google/deepvariant:${DV_VERSION}, jmcdani20/hap.py:${HAPPY_VERSION}
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --resource-dir) RESOURCE_DIR="$2"; shift 2 ;;
        --engine)       ENGINE="$2";       shift 2 ;;
        -h|--help)      usage 0 ;;
        *)              echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

detect_engine() {
    if [[ -n "${ENGINE}" ]]; then echo "${ENGINE}"; return; fi
    if command -v singularity >/dev/null 2>&1; then echo singularity; return; fi
    if command -v apptainer   >/dev/null 2>&1; then echo apptainer;   return; fi
    if command -v docker      >/dev/null 2>&1; then echo docker;      return; fi
    die "no container engine found (need singularity, apptainer, or docker); pass --engine"
}

command -v curl     >/dev/null 2>&1 || die "curl not found"
command -v samtools >/dev/null 2>&1 || die "samtools not found"

readonly REF_DIR="${RESOURCE_DIR}/reference"
readonly BENCH_DIR="${RESOURCE_DIR}/benchmark"
readonly IMG_DIR="${RESOURCE_DIR}/images"
mkdir -p "${REF_DIR}" "${BENCH_DIR}" "${IMG_DIR}"

fetch() {
    # fetch URL DEST — download to a temp file then move into place, so a
    # half-finished download is never mistaken for a completed one on resume.
    local url="$1" dest="$2"
    if [[ -s "${dest}" ]]; then
        echo "  exists: $(basename "${dest}")"
        return
    fi
    echo "  downloading: $(basename "${dest}")"
    curl -fL --retry 3 --retry-delay 5 -o "${dest}.part" "${url}"
    mv "${dest}.part" "${dest}"
}

echo "[1/4] Reference (GRCh38 no-alt analysis set)"
readonly REF_FASTA="${REF_DIR}/GRCh38_no_alt_analysis_set.fasta"
if [[ ! -s "${REF_FASTA}" ]]; then
    fetch "${REF_URL}" "${REF_FASTA}.gz"
    echo "  decompressing reference ..."
    gunzip -c "${REF_FASTA}.gz" > "${REF_FASTA}"
    rm -f "${REF_FASTA}.gz"
else
    echo "  exists: $(basename "${REF_FASTA}")"
fi
[[ -s "${REF_FASTA}.fai" ]] || { echo "  indexing reference ..."; samtools faidx "${REF_FASTA}"; }

echo "[2/4] GIAB HG002 v4.2.1 truth set"
fetch "${GIAB_BASE}/${TRUTH_VCF_NAME}"      "${BENCH_DIR}/${TRUTH_VCF_NAME}"
fetch "${GIAB_BASE}/${TRUTH_VCF_NAME}.tbi"  "${BENCH_DIR}/${TRUTH_VCF_NAME}.tbi"
fetch "${GIAB_BASE}/${TRUTH_BED_NAME}"      "${BENCH_DIR}/${TRUTH_BED_NAME}"

CONTAINER_ENGINE="$(detect_engine)"
readonly CONTAINER_ENGINE
echo "[3/4] DeepVariant image (engine: ${CONTAINER_ENGINE}, version: ${DV_VERSION})"
case "${CONTAINER_ENGINE}" in
    singularity|apptainer)
        readonly DV_SIF="${IMG_DIR}/deepvariant_${DV_VERSION}.sif"
        if [[ -s "${DV_SIF}" ]]; then
            echo "  exists: $(basename "${DV_SIF}")"
        else
            echo "  pulling DeepVariant SIF ..."
            "${CONTAINER_ENGINE}" pull "${DV_SIF}" "docker://google/deepvariant:${DV_VERSION}"
        fi
        ;;
    docker)
        echo "  pulling DeepVariant docker image ..."
        docker pull "google/deepvariant:${DV_VERSION}"
        ;;
esac

echo "[4/4] hap.py image (engine: ${CONTAINER_ENGINE}, version: ${HAPPY_VERSION})"
case "${CONTAINER_ENGINE}" in
    singularity|apptainer)
        readonly HAPPY_SIF="${IMG_DIR}/happy_${HAPPY_VERSION}.sif"
        if [[ -s "${HAPPY_SIF}" ]]; then
            echo "  exists: $(basename "${HAPPY_SIF}")"
        else
            echo "  pulling hap.py SIF ..."
            "${CONTAINER_ENGINE}" pull "${HAPPY_SIF}" "docker://jmcdani20/hap.py:${HAPPY_VERSION}"
        fi
        ;;
    docker)
        echo "  pulling hap.py docker image ..."
        docker pull "jmcdani20/hap.py:${HAPPY_VERSION}"
        ;;
esac

cat <<EOF

Done. Resources in: ${RESOURCE_DIR}
  reference: ${REF_FASTA}
  truth VCF: ${BENCH_DIR}/${TRUTH_VCF_NAME}
  truth BED: ${BENCH_DIR}/${TRUTH_BED_NAME}
  engine:    ${CONTAINER_ENGINE}

Next:
  ./run_eval.sh --our-bam PGPHASE_PHASED.bam --raw-bam RAW.bam --outdir results
EOF
