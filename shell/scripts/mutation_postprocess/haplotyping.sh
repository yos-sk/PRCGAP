#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=30G
#SBATCH -p mjobs,rjobs

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2
INPUT_BAM=$3
KMER_FILE=$4
IMAGE=$5

# 3. haplotyping 
singularity exec ${IMAGE} mutation_postprocess haplotype \
    -i ${OUTPUT_DIR}/realign/parsed_vcf.bed \
    -b ${INPUT_BAM} \
    -r ${OUTPUT_DIR}/realign/realign.bam \
    -p ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz \
    -k ${KMER_FILE} \
    -d 0.98 \
1>${OUTPUT_DIR}/${SAMPLE}.haplotyped.bed 2>${OUTPUT_DIR}/haplotyping.log
    
echo ${?}
