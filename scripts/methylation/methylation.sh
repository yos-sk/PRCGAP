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

mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}


cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${WORK_DIR}/reference.fa
samtools faidx ${WORK_DIR}/reference.fa
samtools view -@ 16 -F 2308 -Shb ${INPUT_BAM} > ${INPUT_BAM%.bam}.filtered.bam
samtools index -@ 16 ${INPUT_BAM%.bam}.filtered.bam


if [ ${TYPE} = "HiFi" ]; then
    ~/bin/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
        --bam ${INPUT_BAM%.bam}.filtered.bam \
        --output-prefix ${OUTPUT_DIR}/${SAMPLE}_methylation \
        --model ~/bin/pb-CpG-tools/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
        --threads 16 \
        --pileup-mode model
    find ${OUTPUT_DIR} -type f -name "${SAMPLE}_methylation*.bed" | while read f; do
        bgzip -f ${f} 
        tabix -p bed ${f}.gz 
    done
else
    modkit pileup \
        --cpg \
        -t 16 \
        --ref ${WORK_DIR}/reference.fa \
        --combine-strands \
        ${INPUT_BAM%.bam}.filtered.bam \
        ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
    bgzip -f ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
    tabix -p bed ${OUTPUT_DIR}/${SAMPLE}_methylation.bed.gz
fi

rm ${INPUT_BAM%.bam}.filtered.bam

echo ${?}
