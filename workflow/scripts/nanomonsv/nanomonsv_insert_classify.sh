#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SV_FILE=$1
OUTPUT_DIR=$2
OUTPUT_FILE=$3
ASSEMBLY_HAP1=$4
ASSEMBLY_HAP2=$5
GTF_FILE=$6
LINE1_BED=$7
SCRIPT_DIR=$8

mkdir -p ${OUTPUT_DIR}

python3 ${SCRIPT_DIR}/nanomonsv/insert_classify/insert_classify.py \
    ${SV_FILE} \
    ${OUTPUT_DIR}/${OUTPUT_FILE} \
    ${ASSEMBLY_HAP1} \
    ${ASSEMBLY_HAP2} \
    ${GTF_FILE} \
    ${LINE1_BED}

echo ${?}
