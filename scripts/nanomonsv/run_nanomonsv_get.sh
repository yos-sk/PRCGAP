#!/bin/bash
#SBATCH --mem-per-cpu=30G
#SBATCH -c 8
#SBATCH -p mjobs,rjobs

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
OUTPUT_DIR=$5
ASSENBLY_HAP1=$6
ASSENBLY_HAP2=$7
DATA=$8

singularity exec nanomonsv-latest.sif \
    bash nanomonsv_get.sh \
        ${TUMOR} \
        ${NORMAL} \
        ${TUMOR_BAM} \
        ${NORMAL_BAM} \
        ${OUTPUT_DIR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${DATA}
