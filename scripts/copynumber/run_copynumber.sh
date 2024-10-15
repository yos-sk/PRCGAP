#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu=8G
#SBATCH -p mjobs,rjobs

TUMOR=$1
NORMAL=$2
ASSEMBLY_HAP1=$3
ASSEMBLY_HAP2=$4
TUMOR_BAM=$5
NORMAL_BAM=$6
REFERENCE=$7
WORK_DIR=$8
OUTPUT_DIR=$9
DNANN_MODEL=${10}

singularity exec ../CN_utils-latest.sif \
   copynumber.sh \
        ${TUMOR} \
        ${NORMAL} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${TUMOR_BAM} \
        ${NORMAL_BAM} \
        ${REFERECNE} \
        ${WORK_DIR} \
        ${OUTPUT_DIR} \
        /tools/models/attcc-alpha.knm 
