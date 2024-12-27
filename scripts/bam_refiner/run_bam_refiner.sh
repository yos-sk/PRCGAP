#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu=8G
#SBATCH -p rjobs,mjobs

SAMPLE=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
FASTQ=$4
WORK_DIR=$5
OUTPUT_DIR=$6
DATA=$7

FASTQ_DIR=`dirname ${FASTQ}`

singularity exec \
    --bind /lustre1/,${FASTQ_DIR} \
    ./images/bam_refiner_v0.3.3.sif \
    /bin/bash \
    ./scripts/bam_refiner/bam_refiner.sh \
    $SAMPLE \
    $ASSEMBLY_HAP1 \
    $ASSEMBLY_HAP2 \
    $FASTQ \
    true \
    $WORK_DIR \
    $OUTPUT_DIR \
    16 \
    $DATA
