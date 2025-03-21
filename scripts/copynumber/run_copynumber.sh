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
HAP1_SATELLITE=${10}
HAP2_SATELLITE=${11}
SEX=${12}

ASSEMBLY_DIR=`dirname ${ASSEMBLY_HAP1}`
REFERENCE_DIR=`dirname ${REFERENCE}`

#singularity exec \
#    --bind ${ASSEMBLY_DIR}/,${REFERENCE_DIR}/ \
#    ./images/CN_utils-v0.1.3.sif \
bash ./scripts/copynumber//copynumber.sh \
    ${TUMOR} \
    ${NORMAL} \
    ${ASSEMBLY_HAP1} \
    ${ASSEMBLY_HAP2} \
    ${TUMOR_BAM} \
    ${NORMAL_BAM} \
    ${REFERENCE} \
    ${WORK_DIR} \
    ${OUTPUT_DIR} \
    ${HAP1_SATELLITE} \
    ${HAP2_SATELLITE} \
    ${SEX}
