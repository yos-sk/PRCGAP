#!/bin/bash
# CBS segmentation for one haplotype, then clip segments that span assembly
# gaps. Ploidy is estimated automatically unless PLOIDY_ARG is non-empty.
#
# Required positional args:
#   $1  TUMOR       tumor sample name
#   $2  HAP         hap1 | hap2
#   $3  CN_TSV      copynumber TSV from copynumber_depth.sh
#   $4  REF_TABLE   published reference table for this haplotype
#   $5  OUT_DIR     destination dir
#   $6  PLOIDY_ARG  manual ploidy override, or empty to auto-estimate
#   $7  BIN_WIDTH   bin width for the depth-ratio mode/ploidy estimate
#   $8  SCRIPT_DIR  absolute path to workflow/scripts

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
HAP=$2
CN_TSV=$3
REF_TABLE=$4
OUT_DIR=$5
PLOIDY_ARG=${6:-}
BIN_WIDTH=${7:-0.05}
SCRIPT_DIR=$8

mkdir -p "${OUT_DIR}"

CBS="${OUT_DIR}/${TUMOR}.${HAP}.cbs.txt"

if [ -n "${PLOIDY_ARG}" ]; then
    PLOIDY_FLAGS=(-p "${PLOIDY_ARG}")
else
    # cbs.R writes the resolved ploidy to --ploidy-out either way.
    PLOIDY_FLAGS=(-a)
fi

Rscript "${SCRIPT_DIR}"/copynumber/cbs.R \
    -i "${CN_TSV}" \
    -s "${TUMOR}" \
    -o "${CBS}" \
    "${PLOIDY_FLAGS[@]}" \
    -w "${BIN_WIDTH}" \
    --ploidy-out "${OUT_DIR}/${TUMOR}.${HAP}.ploidy"

# split_gaps.py derives gap intervals from the ref.table (replaying
# copynumber_window.py's index loop), then clips CBS segments spanning gaps.
python3 "${SCRIPT_DIR}"/copynumber/split_gaps.py \
    "${REF_TABLE}" \
    "${CBS}" \
    > "${OUT_DIR}/${TUMOR}.${HAP}.cbs.split.txt"

echo "[copynumber_segment] done: ${TUMOR} ${HAP}"
