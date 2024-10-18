#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -c 1
#SBATCH --mem-per-cpu=30G

CLAIRS_RESULTS=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
WORK_DIR=$4
OUTPUT_DIR=$5

ASSEMBLY_DIR=`dirname ${ASSEMBLY_HAP1}`

singularity exec \
    --bind ${ASSEMBLY_DIR}/ \
    ./images/mutation_postprocess-v0.1.1.sif \
    bash ./scripts/clairs/postprocess.sh \
        ${CLAIRS_RESULTS} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${WORK_DIR} \
        ${OUTPUT_DIR} 

