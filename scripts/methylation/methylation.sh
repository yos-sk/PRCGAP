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

mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}

cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${WORK_DIR}/reference.fa
samtools faidx ${WORK_DIR}/reference.fa

methylation_utils modbam-utils \
    -i ${INPUT_BAM} \
    -r ${INPUT_MODBAM} \
    -o ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged.bam \
    -t 16 \
    > ${WORK_DIR}/${SAMPLE}_methylation_tags.txt

samtools index ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged.bam

modkit pileup \
    --cpg \
    -t 16 \
    --ref ${WORK_DIR}/reference.fa \
    --with-header \
    ${OUTPUT_DIR}/${SAMPLE}_methylation_tagged.bam \
    ${OUTPUT_DIR}/${SAMPLE}_methylation.bed

gzip -f ${OUTPUT_DIR}/${SAMPLE}_methylation.bed

echo ${?}
