#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
INPUT_BAM=$2
OUTPUT_DIR=$3
OUTPUT_PREFIX=${OUTPUT_DIR}/${SAMPLE}

mkdir -p ${OUTPUT_DIR}

nanomonsv parse \
    ${INPUT_BAM} \
    ${OUTPUT_PREFIX} \
    --split_alignment_check_margin 100 
    
echo ${?}
