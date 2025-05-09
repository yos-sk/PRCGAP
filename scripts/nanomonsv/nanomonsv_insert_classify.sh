#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SV_FILE=$1
OUTPUT_DIR=$2
OUTPUT_FILE=$3
GTF_FILE=$4
LINE1_BED=$5

if [ ! -e ${OUTPUT_DIR}/reference.fa.sa ]; then
    bwa index ${OUTPUT_DIR}/reference.fa
fi

nanomonsv insert_classify \
    ${SV_FILE} \
    ${OUTPUT_DIR}/${OUTPUT_FILE} \
    ${OUTPUT_DIR}/reference.fa \
    ${GTF_FILE} \
    ${LINE1_BED}

echo ${?}
