#!/bin/bash
# Combine the two per-haplotype tandem repeat BEDs into the tabix-indexed file
# nanomonsv get reads when filtering indel-like SVs.
# Runs inside the annotation singularity container (bgzip/tabix).
#
# Required positional args:
#   $1  SAMPLE      sample name (filename prefix)
#   $2  HAP1_BED    hap1 simple repeat BED
#   $3  HAP2_BED    hap2 simple repeat BED
#   $4  OUTPUT_DIR  destination dir

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP1_BED=$2
HAP2_BED=$3
OUTPUT_DIR=$4

mkdir -p "${OUTPUT_DIR}"

BED="${OUTPUT_DIR}/${SAMPLE}.simple_repeats.bed"
export LC_ALL=C
# The haplotypes have disjoint contig names, so this is a concatenation; the
# sort is what tabix needs, not a merge across haplotypes.
cat "${HAP1_BED}" "${HAP2_BED}" | sort -k1,1 -k2,2n > "${BED}"

bgzip -f "${BED}"
tabix -f -p bed "${BED}.gz"

echo "[simple_repeat_merge] done: ${SAMPLE} ($(zcat "${BED}.gz" | wc -l) regions)"
