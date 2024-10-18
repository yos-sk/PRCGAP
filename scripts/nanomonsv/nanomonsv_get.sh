#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7
DATA=$8
SIMPLE_REPEAT=$9

TUMOR_PREFIX=${OUTPUT_DIR}/${TUMOR}
NORMAL_PREFIX=${OUTPUT_DIR}/${NORMAL}

<<_
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa

if [ $DATA = "hifi" ];then
    nanomonsv get \
        ${TUMOR_PREFIX} \
        ${TUMOR_BAM} \
        ${OUTPUT_DIR}/reference.fa \
        --control_prefix ${NORMAL_PREFIX} \
        --control_bam ${NORMAL_BAM} \
        --processes 8 \
        --single_bnd \
        --use_racon \
        --max_memory_minimap2 16 \
        --qv25
else
    nanomonsv get \
        ${TUMOR_PREFIX} \
        ${TUMOR_BAM} \
        ${OUTPUT_DIR}/reference.fa \
        --control_prefix ${NORMAL_PREFIX} \
        --control_bam ${NORMAL_BAM} \
        --processes 8 \
        --single_bnd \
        --use_racon \
        --max_memory_minimap2 16 \
        --qv15
fi
_
python3 /nanomonsv/misc/add_simple_repeat.py \
    ${OUTPUT_DIR}/${TUMOR}.nanomonsv.result.txt \
    ${OUTPUT_DIR}/${TUMOR}.nanomonsv.result.filt.txt \
    ${SIMPLE_REPEAT}

echo ${?}
