#!/bin/bash
# Combine the two per-haplotype LINE-1 BEDs into the tabix-indexed file
# nanomonsv insert_classify reads as its LINE-1 source database.
# Runs inside the annotation singularity container (bgzip/tabix).
#
# Required positional args:
#   $1  SAMPLE      sample name (filename prefix)
#   $2  HAP1_BED    hap1 LINE-1 BED
#   $3  HAP2_BED    hap2 LINE-1 BED
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

BED="${OUTPUT_DIR}/${SAMPLE}.LINE1.bed"
export LC_ALL=C
cat "${HAP1_BED}" "${HAP2_BED}" | sort -k1,1 -k2,2n > "${BED}"

bgzip -f "${BED}"
tabix -f -p bed "${BED}.gz"

echo "[line1_merge] done: ${SAMPLE} ($(zcat "${BED}.gz" | wc -l) elements)"
