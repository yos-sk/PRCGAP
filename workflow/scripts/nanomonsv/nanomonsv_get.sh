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
THREAD=${10:-8}
SCRIPTS_DIR=${11}

TUMOR_PREFIX=${OUTPUT_DIR}/${TUMOR}
NORMAL_PREFIX=${OUTPUT_DIR}/${NORMAL}


cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa

if [ $DATA = "hifi" ];then
    nanomonsv get \
        ${TUMOR_PREFIX} \
        ${TUMOR_BAM} \
        ${OUTPUT_DIR}/reference.fa \
        --control_prefix ${NORMAL_PREFIX} \
        --control_bam ${NORMAL_BAM} \
        --processes ${THREAD} \
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
        --processes ${THREAD} \
        --single_bnd \
        --use_racon \
        --max_memory_minimap2 16 \
        --qv15
fi

# The tandem repeat flag only gates the test; whether a call is a copy-number
# change of the array is decided from its inserted sequence, so score that first.
# reference.fa above is the two haplotypes concatenated, which is what the
# breakpoint contigs are named against.
SCORE_ARGS=()
if [ -n "${SIMPLE_REPEAT}" ]; then
    python3 ${SCRIPTS_DIR}/nanomonsv/simple_repeat_score.py \
        ${OUTPUT_DIR}/${TUMOR}.nanomonsv.result.txt \
        ${OUTPUT_DIR}/${TUMOR}.simple_repeat_scores.tsv \
        ${OUTPUT_DIR}/reference.fa
    SCORE_ARGS=(--scores ${OUTPUT_DIR}/${TUMOR}.simple_repeat_scores.tsv)
fi

python3 ${SCRIPTS_DIR}/nanomonsv/add_simple_repeat.py \
    ${OUTPUT_DIR}/${TUMOR}.nanomonsv.result.txt \
    ${OUTPUT_DIR}/${TUMOR}.nanomonsv.result.filt.txt \
    ${SIMPLE_REPEAT} \
    "${SCORE_ARGS[@]}"

echo ${?}
