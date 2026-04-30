#!/bin/bash
#SBATCH --mem-per-cpu=16G
#SBATCH -c 8
#SBATCH -p mjobs,rjobs

SV_FILE=$1
OUTPUT_DIR=$2
OUTPUT_FILE=$3
HAP1_ASSEMBLY=$4
HAP2_ASSEMBLY=$5
GTF_FILE=$6
LINE1_BED=$7

singularity exec ./images/nanomonsv-v0.8.0.sif \
    bash ./scripts/nanomonsv/nanomonsv_insert_classify.sh \
        ${SV_FILE} \
        ${OUTPUT_DIR} \
        ${OUTPUT_FILE} \
        ${HAP1_ASSEMBLY} \
        ${HAP2_ASSEMBLY} \
        ${GTF_FILE} \
        ${LINE1_BED} 
