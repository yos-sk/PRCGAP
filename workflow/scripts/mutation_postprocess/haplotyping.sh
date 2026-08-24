#!/bin/bash
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2
INPUT_BAM=$3
KMER_RATIO=$4
HAP1_LIST=$5
HAP2_LIST=$6

mutation_postprocess haplotype \
    -i ${OUTPUT_DIR}/realign/parsed_vcf.bed \
    -b ${INPUT_BAM} \
    -r ${OUTPUT_DIR}/realign/realign.bam \
    -p ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz \
    -d 0.98 \
    -k ${KMER_RATIO} \
    --hap1-list ${HAP1_LIST} \
    --hap2-list ${HAP2_LIST} \
    1>${OUTPUT_DIR}/${SAMPLE}.haplotyped.bed 2>${OUTPUT_DIR}/haplotyping.log

echo ${?}