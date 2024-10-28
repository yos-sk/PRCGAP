#!/bin/bash
#SBATCH --mem-per-cpu=30G
#SBATCH -c 1
#SBATCH -p mjobs,rjobs

NANOMONSV_RESULT=$1
SUPPORT_READ_FILE=$2
OUTPUT_PREFIX=$3

singularity exec \
    ./images/nanomonsv-devel.sif \
    bash ./scripts/nanomonsv/nanomonsv_connect.sh \
        ${NANOMONSV_RESULT} \
        ${SUPPORT_READ_FILE} \
        ${OUTPUT_PREFIX}
