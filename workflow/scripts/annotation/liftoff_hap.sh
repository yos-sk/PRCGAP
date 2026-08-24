#!/bin/bash
# Lift the staged GRCh38 gene annotation onto one haplotype.
# Runs inside the liftoff singularity container (liftoff + minimap2).
#
# Required positional args:
#   $1  SAMPLE       sample name (filename prefix)
#   $2  HAP          hap1 | hap2
#   $3  FASTA        that haplotype's assembly FASTA
#   $4  GRCH38       staged GRCh38 FASTA (liftoff_reference.sh)
#   $5  GRCH38_GTF   staged GRCh38 GTF
#   $6  OUT_GFF      per-haplotype liftoff GFF to write
#   $7  THREADS      liftoff -p

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP=$2
FASTA=$3
GRCH38=$4
GRCH38_GTF=$5
OUT_GFF=$6
THREADS=${7:-8}

# liftoff caches its GTF database and intermediate minimap2 output under -dir;
# a stale cache from a previous (possibly different) reference silently wins, so
# start each haplotype from an empty scratch dir of its own.
WORK_DIR="$(dirname "${OUT_GFF}")/${HAP}"
rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"

liftoff \
    -g "${GRCH38_GTF}" \
    -o "${OUT_GFF}" \
    -u "${WORK_DIR}/unmapped_features.txt" \
    -dir "${WORK_DIR}" \
    -p "${THREADS}" \
    -m /usr/local/bin/minimap2 \
    "${FASTA}" \
    "${GRCH38}"

echo "[liftoff_hap] done: ${SAMPLE} ${HAP}"
