#!/bin/bash
# Standalone haplotyping script (without SLURM)
# Used by Snakemake for mutation postprocessing

set -euo pipefail

SAMPLE=$1
OUTPUT_DIR=$2
INPUT_BAM=$3

mutation_postprocess haplotype \
    -i ${OUTPUT_DIR}/realign/parsed_vcf.bed \
    -b ${INPUT_BAM} \
    -r ${OUTPUT_DIR}/realign/realign.bam \
    -p ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz \
    -d 0.98 \
    1>${OUTPUT_DIR}/${SAMPLE}.haplotyped.bed 2>${OUTPUT_DIR}/haplotyping.log
