#!/bin/bash
# Split BED file for parallel pileup processing
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

PARSED_BED=$1
OUTPUT_DIR=$2
PILEUP_TASKS=$3

mkdir -p ${OUTPUT_DIR}/pileup/workspace
split -n l/16 ${PARSED_BED} ${OUTPUT_DIR}/pileup/workspace/input_
ls ${OUTPUT_DIR}/pileup/workspace/input_* | while read f; do
    echo -e "${f}\t${f}.pileup"
done > ${PILEUP_TASKS}

echo ${?}
