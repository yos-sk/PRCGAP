#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu=8G
#SBATCH -p mjobs,rjobs

set -xv
set -o errexit
set -o nounset
set -o pipefail

INPUT_VCF=$1
HAP1_FASTA=$2
HAP2_FASTA=$3
OUTPUT_DIR=$4
TOOL=$5
IMAGE=$6

# 1. Realignment by bwa mem 
singularity exec ${IMAGE} bwa index ${OUTPUT_DIR}/realign/realign_ref.fasta
singularity exec ${IMAGE} bwa mem -a -k 50 -c 1000000 -t 16 ${OUTPUT_DIR}/realign/realign_ref.fasta ${OUTPUT_DIR}/realign/realign_ref.fasta | singularity exec ${IMAGE} samtools view -Shb > ${OUTPUT_DIR}/realign/realign.unsorted
singularity exec ${IMAGE} samtools sort -@ 16 ${OUTPUT_DIR}/realign/realign.unsorted > ${OUTPUT_DIR}/realign/realign.bam
singularity exec ${IMAGE} samtools index -@ 16 ${OUTPUT_DIR}/realign/realign.bam

rm ${OUTPUT_DIR}/realign/realign.unsorted

echo ${?}

