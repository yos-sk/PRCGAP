#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -c 16
#SBATCH --mem-per-cpu=8G


TUMOR_BAM=$1
CONTROL_BAM=$2
ASSEMBLY_HAP1=$3
ASSEMBLY_HAP2=$4
OUTPUT_DIR=$5

singularity exec ./images/clairs-latest.sif \
    bash ./scripts/clairs/clairs.sh \
        ${TUMOR_BAM} \
        ${CONTROL_BAM} \
        ${OUTPUT_DIR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2}

