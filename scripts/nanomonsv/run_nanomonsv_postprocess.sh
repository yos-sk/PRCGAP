#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=30G
#SBATCH -p mjobs,rjobs

SAMPLE=$1
OUTPUT_DIR=$2
SIMPLE_REPEAT=$3
BAM_FILE=$4

singularity exec ./images/nannomonsv_postprocess-0.2.3.sif \
    bash ./scripts/nanomonsv/nanomonsv_postprocess.sh \
        ${SAMPLE} \
        ${OUTPUT_DIR} \
        ${SIMPLE_REPEAT} \
        ${BAM_FILE} \

