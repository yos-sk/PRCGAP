#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=150G
#SBATCH -J pileup
#SBATCH -e log/pileup_%j.err
#SBATCH -o log/pileup_%j.out
#SBATCH -p mjobs,rjobs

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2
IMAGE=$3

cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
singularity exec ${IMAGE} bgzip -@ 16 -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
singularity exec ${IMAGE} tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz

echo ${?}
