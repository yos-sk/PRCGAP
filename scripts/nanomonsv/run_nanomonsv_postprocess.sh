#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=30G
#SBATCH -p mjobs,rjobs

SAMPLE=$1
OUTPUT_DIR=$2
BAM_FILE=$3

singularity exec \
    ./images/nanomonsv_postprocess-v0.2.4b.sif \
    bash ./scripts/nanomonsv/nanomonsv_postprocess.sh \
        ${SAMPLE} \
        ${OUTPUT_DIR} \
        ${BAM_FILE} 

