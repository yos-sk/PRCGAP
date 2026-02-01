#!/bin/bash
# Standalone realignment script (without SLURM)
# Used by Snakemake for mutation postprocessing

set -euo pipefail

INPUT_VCF=$1
HAP1_FASTA=$2
HAP2_FASTA=$3
OUTPUT_DIR=$4
TOOL=$5

bwa index ${OUTPUT_DIR}/realign/realign_ref.fasta
bwa mem -a -k 50 -c 1000000 -t 16 ${OUTPUT_DIR}/realign/realign_ref.fasta ${OUTPUT_DIR}/realign/realign_ref.fasta | \
samtools view -Shb > ${OUTPUT_DIR}/realign/realign.unsorted
samtools sort -@ 16 ${OUTPUT_DIR}/realign/realign.unsorted > ${OUTPUT_DIR}/realign/realign.bam
samtools index -@ 16 ${OUTPUT_DIR}/realign/realign.bam

rm ${OUTPUT_DIR}/realign/realign.unsorted
