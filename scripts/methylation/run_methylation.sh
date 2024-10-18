#!/bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
#SBATCH -p mjobs,rjobs

SAMPLE=$1
INPUT_BAM=$2
INPUT_MODBAM=$3
WORK_DIR=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7

INPUT_MODBAM_DIR=`dirname ${INPUT_MODBAM}`
ASSEMBLY_DIR=`dirname ${ASSEMBLY_HAP1}`

singularity exec \
    --bind ${INPUT_MODBAM_DIR}/,${ASSEMBLY_DIR}/ \
    ./images/methylation_utils-v0.1.1.sif \
    bash ./scripts/methylation/methylation.sh \
        ${SAMPLE} \
        ${INPUT_BAM} \
        ${INPUT_MODBAM} \
        ${WORK_DIR} \
        ${OUTPUT_DIR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2}
