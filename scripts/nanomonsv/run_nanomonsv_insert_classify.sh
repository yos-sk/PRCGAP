#!/bin/bash
#SBATCH --mem-per-cpu=16G
#SBATCH -c 8
#SBATCH -p mjobs,rjobs

SV_FILE=$1
OUTPUT_DIR=$2
OUTPUT_FILE=$3
GTF_FILE=$4
LINE1_BED=$5

singularity exec ./images/nanomonsv-v0.8.0.sif \
    bash ./scripts/nanomonsv/nanomonsv_insert_classify.sh \
        ${SV_FILE} \
        ${OUTPUT_DIR} \
        ${OUTPUT_FILE} \
        ${GTF_FILE} \
        ${LINE1_BED} \
