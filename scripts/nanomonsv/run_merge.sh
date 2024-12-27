#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=30G
#SBATCH -p mjobs,rjobs

NANOMONSV_RESULT_1=$1
NANOMONSV_RESULT_2=$2
OUTPUT_FILE=$3

singularity exec ./images/nanomonsv_postprocess-v0.2.5.sif \
    bash ./scripts/nanomonsv/merge.sh \
        ${NANOMONSV_RESULT_1} \
        ${NANOMONSV_RESULT_2} \
        ${OUTPUT_FILE}

