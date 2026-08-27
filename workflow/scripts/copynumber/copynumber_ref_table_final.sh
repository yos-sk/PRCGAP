#!/bin/bash
# Publish the two per-haplotype reference tables. For a male sample the sex
# chromosomes have to be reconciled across haplotypes (chrX on one, chrY on the
# other), which is the one step that needs both tables at once.
#
# Required positional args:
#   $1  SEX         female | male
#   $2  RAW_HAP1    hap1 table from copynumber_ref_table.sh
#   $3  RAW_HAP2    hap2 table
#   $4  OUT_HAP1    published hap1 table
#   $5  OUT_HAP2    published hap2 table
#   $6  SCRIPT_DIR  absolute path to workflow/scripts
#   $7  PAF_HAP1    hap1 PAF from copynumber_ref_table.sh (male only)
#   $8  PAF_HAP2    hap2 PAF (male only)

set -xv
set -o errexit
set -o nounset
set -o pipefail

SEX=$1
RAW_HAP1=$2
RAW_HAP2=$3
OUT_HAP1=$4
OUT_HAP2=$5
SCRIPT_DIR=$6
PAF_HAP1=${7:-}
PAF_HAP2=${8:-}

mkdir -p "$(dirname "${OUT_HAP1}")" "$(dirname "${OUT_HAP2}")"

cp "${RAW_HAP1}" "${OUT_HAP1}"
cp "${RAW_HAP2}" "${OUT_HAP2}"

if [ "${SEX}" = "male" ]; then
    # With the PAFs the overlap between the chrX and chrY records is settled by
    # the alignment actually present in the disputed interval instead of by
    # envelope size; without them the script falls back to the larger span.
    PAF_ARGS=()
    if [ -n "${PAF_HAP1}" ] && [ -n "${PAF_HAP2}" ]; then
        PAF_ARGS=(--paf "${PAF_HAP1}" "${PAF_HAP2}")
    fi
    python3 "${SCRIPT_DIR}"/copynumber/postprocess_sex_chrom.py \
        --hap1 "${OUT_HAP1}" \
        --hap2 "${OUT_HAP2}" \
        --out1 "${OUT_HAP1}" \
        --out2 "${OUT_HAP2}" \
        "${PAF_ARGS[@]}"
fi

echo "[copynumber_ref_table_final] done (${SEX})"
