#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=200G
#SBATCH -p mjobs,rjobs

set -xv
set -o errexit
set -o nounset
set -o pipefail

INPUT_BAM=$1
TASK_FILE=$2
INPUT_FASTA=$3
IMAGE=$4

TASK=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${TASK_FILE})

INPUT_BED=$(echo ${TASK} | awk '{print $1}')
OUTPUT_FILE=$(echo ${TASK} | awk '{print $2}')
singularity exec ${IMAGE} samtools mpileup -l ${INPUT_BED} -f ${INPUT_FASTA} --output-QNAME ${INPUT_BAM} | \
    awk '{print $1 "\t" $2 -1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > ${OUTPUT_FILE}

echo ${?}
