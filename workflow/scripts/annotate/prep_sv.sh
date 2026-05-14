#!/bin/bash
# Filter raw nanomonsv result and extract breakpoint BED.
# Runs inside the point_mutation_postprocess singularity container
# (python3 + grep + awk are available).
#
# Required positional args:
#   $1  NANOMONSV    *.nanomonsv.new_result.sv_typed.insert_classified.txt
#   $2  OUTPUT_DIR   workspace dir to place outputs
#   $3  SAMPLE       tumor sample name (used as filename prefix)
#   $4  SCRIPT_DIR   absolute path to workflow/scripts/annotate

set -xv
set -o errexit
set -o nounset
set -o pipefail

NANOMONSV=$1
OUTPUT_DIR=$2
SAMPLE=$3
SCRIPT_DIR=$4

mkdir -p "${OUTPUT_DIR}"

FILT_TXT="${OUTPUT_DIR}/${SAMPLE}.filt.txt"
PASS_TXT="${OUTPUT_DIR}/${SAMPLE}.filt.pass.txt"
BP_BED="${OUTPUT_DIR}/${SAMPLE}.bp.bed"

python3 "${SCRIPT_DIR}/filter_PRCGAP.py" -i "${NANOMONSV}" -o "${FILT_TXT}"

head -n 1 "${FILT_TXT}" > "${PASS_TXT}"
grep "PASS" "${FILT_TXT}" >> "${PASS_TXT}" || true

awk 'NR!=1 { print $1 "\t" $2-1 "\t" $2 "\t" $8 "\n" $4 "\t" $5-1 "\t" $5 "\t" $8 }' \
    "${PASS_TXT}" > "${BP_BED}"

echo "[prep_sv] done: ${SAMPLE}"
