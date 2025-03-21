#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -c 16
#SBATCH --mem-per-cpu=8G

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7

ASSEMBLY_DIR=`dirname ${ASSEMBLY_HAP1}`

singularity exec \
    --bind ${ASSEMBLY_DIR}/ \
    ./images/DeepSomatic-v1.8.0.sif \
    bash ./scripts/deepsomatic/deepsomatic.sh \
        ${TUMOR} \
        ${NORMAL} \
        ${TUMOR_BAM} \
        ${NORMAL_BAM} \
        ${OUTPUT_DIR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} 

