#!/bin/bash

# Stage the CHM13 reference for the copy-number reference tables: female
# samples align against CHM13 without chrY, male samples against the full
# reference. Also writes the chromosome-length tables used by the plot.
#
# Usage: prepare_reference.sh REFERENCE OUT_DIR   (run inside the copynumber image)

set -xv
set -o errexit
set -o nounset
set -o pipefail

REFERENCE=$1
OUT_DIR=$2

mkdir -p ${OUT_DIR}

awk '/^>/ {p = ($0 !~ /^>chrY/)} p' ${REFERENCE} > ${OUT_DIR}/chm13_female.fa.tmp
mv ${OUT_DIR}/chm13_female.fa.tmp ${OUT_DIR}/chm13_female.fa

python3 ./scripts/copynumber/chromosome_length.py ${OUT_DIR}/chm13_female.fa > ${OUT_DIR}/chm13_female_chrom_length.txt
python3 ./scripts/copynumber/chromosome_length.py ${REFERENCE} > ${OUT_DIR}/chm13_male_chrom_length.txt

echo "[prepare_reference] done"
