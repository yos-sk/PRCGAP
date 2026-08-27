#!/bin/bash
# dna-brnn alpha-satellite annotation for one haplotype.
# Runs inside the dna_nn singularity container (dna-brnn + bgzip/tabix).
#
# Required positional args:
#   $1  SAMPLE       sample name (filename prefix)
#   $2  HAP          hap1 | hap2
#   $3  FASTA        that haplotype's assembly FASTA
#   $4  OUTPUT_DIR   destination dir
#   $5  THREADS      dna-brnn threads

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP=$2
FASTA=$3
OUTPUT_DIR=$4
THREADS=${5:-8}

# Shipped with the dna-nn release the container builds from.
MODEL=/opt/dna-nn-0.1/models/attcc-alpha.knm

mkdir -p "${OUTPUT_DIR}"

BED="${OUTPUT_DIR}/${SAMPLE}.${HAP}_dna-brnn.bed"
dna-brnn -Ai "${MODEL}" -t"${THREADS}" "${FASTA}" \
    | sort -k1,1 -k2,2n > "${BED}"
bgzip -f "${BED}"
tabix -f -p bed "${BED}.gz"

echo "[dna_brnn] done: ${SAMPLE} ${HAP}"
