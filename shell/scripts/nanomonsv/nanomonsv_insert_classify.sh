#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SV_FILE=$1
OUTPUT_DIR=$2
OUTPUT_FILE=$3
HAP1_ASSEMBLY=$4
HAP2_ASSEMBLY=$5
GTF_FILE=$6
LINE1_BED=$7

if [ ! -e ${OUTPUT_DIR}/reference_hap1.fa ]; then
    cp ${HAP1_ASSEMBLY} ${OUTPUT_DIR}/reference_hap1.fa
fi

if [ ! -e ${OUTPUT_DIR}/reference_hap2.fa ]; then
    cp ${HAP2_ASSEMBLY} ${OUTPUT_DIR}/reference_hap2.fa
fi

if [ ! -e ${OUTPUT_DIR}/reference_hap1.fa.sa ]; then
    bwa index ${OUTPUT_DIR}/reference_hap1.fa
fi

if [ ! -e ${OUTPUT_DIR}/reference_hap2.fa.sa ]; then
    bwa index ${OUTPUT_DIR}/reference_hap2.fa
fi

python3 ./scripts/nanomonsv/insert_classify/insert_classify.py \
    ${SV_FILE} \
    ${OUTPUT_DIR}/${OUTPUT_FILE} \
    ${OUTPUT_DIR}/reference_hap1.fa \
    ${OUTPUT_DIR}/reference_hap2.fa \
    ${GTF_FILE} \
    ${LINE1_BED}

echo ${?}
