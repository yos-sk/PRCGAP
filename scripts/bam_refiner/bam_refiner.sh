#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP1_ASSEMBLY=$2
HAP2_ASSEMBLY=$3
INPUT=$4
OUTPUT_DIR=$5
THREAD=$6
DATA=$7

ext="${INPUT##*.}"
if [ $ext = "bam" ]; then
    bash run_refine.sh \
        -d \
        -b ${INPUT} \
        -h ${HAP1_ASSEMBLY} \
        -i ${HAP2_ASSEMBLY} \
        -m "--MD" \
        -o ${OUTPUT_DIR} \
        -p \
        -s ${SAMPLE} \
        -t 16 \
        -u ${DATA}
else
    bash run_refine.sh \
        -d \
        -f ${INPUT} \
        -h ${HAP1_ASSEMBLY} \
        -i ${HAP2_ASSEMBLY} \
        -m "--MD" \
        -o ${OUTPUT_DIR} \
        -p \
        -s ${SAMPLE} \
        -t 16 \
        -u ${DATA}
fi

echo ${?}
