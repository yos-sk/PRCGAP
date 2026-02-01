#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
INPUT_BAM=$2
WORK_DIR=$3
OUTPUT_DIR=$4
ASSEMBLY_HAP1=$5
ASSEMBLY_HAP2=$6
TYPE=$7
THREAD=${8:-8}

mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}


cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${WORK_DIR}/reference.fa
samtools faidx ${WORK_DIR}/reference.fa
BASE_NAME=`basename ${INPUT_BAM}`
BAM_FILE=${WORK_DIR}/${BASE_NAME%.bam}.filtered.bam
samtools view -@ ${THREAD} -F 2308 -Shb ${INPUT_BAM} > ${BAM_FILE}
samtools index -@ ${THREAD} ${BAM_FILE}


if [ ${TYPE} = "HiFi" ]; then
    aligned_bam_to_cpg_scores \
        --bam ${BAM_FILE} \
        --output-prefix ${OUTPUT_DIR}/${SAMPLE}_methylation \
        --model /opt/pb-CpG-tools/models/pileup_calling_model.v1.tflite \
        --threads ${THREAD} \
        --pileup-mode model
    find ${OUTPUT_DIR} -type f -name "${SAMPLE}_methylation*.bed" | while read f; do
        bgzip -f ${f}
        tabix -p bed ${f}.gz
    done
else
    modkit pileup \
        --cpg \
        -t ${THREAD} \
        --ref ${WORK_DIR}/reference.fa \
        --combine-strands \
        ${BAM_FILE} \
        ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
    bgzip -f ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
    tabix -p bed ${OUTPUT_DIR}/${SAMPLE}_methylation.bed.gz
fi

rm ${BAM_FILE}*

echo ${?}
