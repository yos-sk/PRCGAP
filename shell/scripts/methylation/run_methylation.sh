#!/bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
#SBATCH -p mjobs,rjobs

SAMPLE=$1
INPUT_BAM=$2
WORK_DIR=$3
OUTPUT_DIR=$4
ASSEMBLY_HAP1=$5
ASSEMBLY_HAP2=$6
TYPE=$7

ASSEMBLY_DIR=`dirname ${ASSEMBLY_HAP1}`

singularity exec \
    --bind /lustre1/,${ASSEMBLY_DIR}/ \
    ./images/methylation_utils-v0.1.2.sif \
    bash ./scripts/methylation/methylation.sh \
        ${SAMPLE} \
        ${INPUT_BAM} \
        ${WORK_DIR} \
        ${OUTPUT_DIR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${TYPE}
