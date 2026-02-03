#!/bin/bash
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

REALIGN_FASTA=$1
OUTPUT_DIR=$2
THREADS=$3

bwa index ${REALIGN_FASTA}
bwa mem -a -k 50 -c 1000000 -t ${THREADS} ${REALIGN_FASTA} ${REALIGN_FASTA} | \
    samtools view -Shb > ${OUTPUT_DIR}/realign/realign.unsorted
samtools sort -@ ${THREADS} ${OUTPUT_DIR}/realign/realign.unsorted > ${OUTPUT_DIR}/realign/realign.bam
samtools index -@ ${THREADS} ${OUTPUT_DIR}/realign/realign.bam

rm ${OUTPUT_DIR}/realign/realign.unsorted

echo ${?}
