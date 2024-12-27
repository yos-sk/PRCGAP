#!/bin/bash
#SBATCH --mem-per-cpu=30G
#SBATCH -c 8
#SBATCH -p mjobs,rjobs

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7
DATA=$8
SIMPLE_REPEAT=$9

ASSEMBLY_DIR=`dirname ${ASSEMBLY_HAP1}`
SIMPLE_REPEAT_DIR=`dirname ${SIMPLE_REPEAT}`

singularity exec \
    --bind /lustre1/,${ASSEMBLY_DIR}/,${SIMPLE_REPEAT_DIR}/ \
    ./images/nanomonsv-devel.sif \
    bash ./scripts/nanomonsv/nanomonsv_get.sh \
        ${TUMOR} \
        ${NORMAL} \
        ${TUMOR_BAM} \
        ${NORMAL_BAM} \
        ${OUTPUT_DIR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${DATA} \
        ${SIMPLE_REPEAT}
