#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -c 1
#SBATCH --mem-per-cpu=30G

CLAIRS_RESULTS=$1
REFERENCE=$2
WORK_DIR=$3
OUTPUT_DIR=$4

singularity exec ./images/mutation_postprocess-v0.1.1.sif \
    bash ./scripts/clairs/postprocess.sh \
        ${CLAIRS_RESULTS} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${WORK_DIR} \
        ${OUTPUT_DIR} 

