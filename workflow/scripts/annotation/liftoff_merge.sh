#!/bin/bash
# Combine the two per-haplotype liftoff GFFs into the workflow's gene
# annotation: a sorted, tabix-indexed GFF plus the GTF nanomonsv
# insert_classify consumes.
# Runs inside the liftoff singularity container (gffread + bgzip/tabix).
#
# Required positional args:
#   $1  SAMPLE      sample name (filename prefix)
#   $2  HAP1_GFF    hap1 liftoff GFF
#   $3  HAP2_GFF    hap2 liftoff GFF
#   $4  OUTPUT_DIR  destination dir

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP1_GFF=$2
HAP2_GFF=$3
OUTPUT_DIR=$4

mkdir -p "${OUTPUT_DIR}"

GFF="${OUTPUT_DIR}/${SAMPLE}.liftoff.gff"
cat "${HAP1_GFF}" "${HAP2_GFF}" \
    | grep -v "^#" \
    | sort -k1,1 -k4,4n -k5,5n -t$'\t' > "${GFF}"

# gffread must read the plain GFF, so convert before bgzipping it.
gffread "${GFF}" -T -o "${OUTPUT_DIR}/${SAMPLE}.liftoff.gtf"
bgzip -f "${OUTPUT_DIR}/${SAMPLE}.liftoff.gtf"

bgzip -f "${GFF}"
tabix -f -p gff "${GFF}.gz"

echo "[liftoff_merge] done: ${SAMPLE}"
