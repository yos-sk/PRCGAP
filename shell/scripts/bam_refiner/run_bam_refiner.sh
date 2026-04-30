#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu=8G
#SBATCH -p rjobs,mjobs

SAMPLE=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
INPUT=$4
OUTPUT_DIR=$5
DATA=$6

echo ${SUMBIT_DIR}
echo ${SUBMIT_CMD}
echo ${SCRIPT_RELATIVE_PATH}

FASTQ_DIR=`dirname ${FASTQ}`

singularity exec \
    --bind /lustre1/,${FASTQ_DIR} \
    ./images/bam_refiner_v0.3.6b.sif \
    /bin/bash \
    ./scripts/bam_refiner/bam_refiner.sh \
    ${SAMPLE} \
    ${ASSEMBLY_HAP1} \
    ${ASSEMBLY_HAP2} \
    ${INPUT} \
    ${OUTPUT_DIR} \
    16 \
    ${DATA}
