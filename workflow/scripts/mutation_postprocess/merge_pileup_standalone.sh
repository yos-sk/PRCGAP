#!/bin/bash
# Standalone merge pileup script (without SLURM)
# Used by Snakemake for mutation postprocessing

set -euo pipefail

SAMPLE=$1
OUTPUT_DIR=$2

cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
bgzip -@ 16 -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz
