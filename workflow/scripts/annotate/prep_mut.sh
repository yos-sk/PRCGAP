#!/bin/bash
# Filter mutation_postprocess's `*.haplotyped.bed` into SNV-only or
# INDEL-only rows. Output remains in the original 14-column,
# header-less format expected by downstream `coordconv` (SNV) /
# `bed2vcf.py` (INDEL) / `add_lift_coords.py`.
#
# Required positional args:
#   $1  HAPLOTYPED_BED   *.haplotyped.bed (mutation_postprocess haplotype output)
#   $2  OUTPUT_DIR       workspace dir
#   $3  SAMPLE           tumor sample
#   $4  TOOL             clairs | deepsomatic
#   $5  MODE             snv | indel

set -xv
set -o errexit
set -o nounset
set -o pipefail

HAPLOTYPED_BED=$1
OUTPUT_DIR=$2
SAMPLE=$3
TOOL=$4
MODE=$5

mkdir -p "${OUTPUT_DIR}"

PREP_BED="${OUTPUT_DIR}/${SAMPLE}.${TOOL}.${MODE}.bed"

if [ "${MODE}" = "snv" ]; then
    awk -F'\t' 'length($4)==length($5)' "${HAPLOTYPED_BED}" > "${PREP_BED}"
else
    awk -F'\t' 'length($4)!=length($5)' "${HAPLOTYPED_BED}" > "${PREP_BED}"
fi

echo "[prep_mut] done: ${SAMPLE}/${TOOL}/${MODE}"
