#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR_BAM=$1
CONTROL_BAM=$2
OUTPUT_DIR=$3
ASSEMBLY_HAP1=$4
ASSEMBLY_HAP2=$5

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa
samtools faidx ${OUTPUT_DIR}/reference.fa

/opt/bin/ClairS/run_clairs \
    -T ${TUMOR_BAM} \
    -N ${CONTROL_BAM} \
    -R ${OUTPUT_DIR}/reference.fa \
    -o ${OUTPUT_DIR} \
    -t 16 \
    -p hifi_sequel2 \
    --include_all_ctgs 

echo ${?}
