#!/bin/bash
# Parse VCF script
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

SCRIPTS_DIR=$1
INPUT_VCF=$2
REFERENCE_FA=$3
OUTPUT_BED=$4
OUTPUT_FASTA=$5
TOOL=$6
OUTPUT_DIR=$7

mkdir -p ${OUTPUT_DIR}/realign
python3 ${SCRIPTS_DIR}/mutation_postprocess/parse_vcf.py \
    -i ${INPUT_VCF} \
    -f ${REFERENCE_FA} \
    -o ${OUTPUT_BED} \
    -g ${OUTPUT_FASTA} \
    -t ${TOOL}

echo ${?}
