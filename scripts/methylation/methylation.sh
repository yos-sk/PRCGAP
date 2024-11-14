#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
INPUT_BAM=$2
INPUT_MODBAM=$3
WORK_DIR=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7
MODE=$8

mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}


cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${WORK_DIR}/reference.fa
samtools faidx ${WORK_DIR}/reference.fa


methylation_utils modbam-utils \
    -i ${INPUT_BAM} \
    -r ${INPUT_MODBAM} \
    -o ${INPUT_BAM%.bam}_methylation_tagged.bam \
    -t 16 \
    > ${WORK_DIR}/${SAMPLE}_methylation_tags.txt

samtools index -@ 16 ${INPUT_BAM%.bam}_methylation_tagged.bam 


if [ ${MODE} = "HiFi" ]; then
    #modkit update-tags \
    #    -m implicit \
    #    -t 16 \
    #    ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged.bam \
    #    ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged_update.bam

    #samtools index ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged_update.bam
    #rm ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged.bam
    #rm ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged.bam.bai

    #modkit pileup \
    #    --cpg \
    #    -t 16 \
    #    --ref ${WORK_DIR}/reference.fa \
    #    --filter-threshold 0.8 \
    #    --with-header \
    #    ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged_update.bam \
    #    ${OUTPUT_DIR}/${SAMPLE}_methylation.bed

    ~/bin/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
        --bam ${INPUT_BAM%.bam}_methylation_tagged.bam \
        --output-prefix ${OUTPUT_DIR}/${SAMPLE}_methylation \
        --model ~/bin/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
        --threads 16 \
        --pileup-mode model
    gzip -f ${OUTPUT_DIR}/${SAMPLE}_methylation*.bed
else
    modkit pileup \
        --cpg \
        -t 16 \
        --ref ${WORK_DIR}/reference.fa \
        --with-header \
        --combine-strands \
        ${INPUT_BAM%.bam}_methylation_tagged.bam \
        ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
    gzip -f ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
fi

echo ${?}
