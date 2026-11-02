#!/bin/bash

# Reference tables of one sample:
#   minimap2 -cx asm5 --secondary=no of CHM13 to the unmasked assembly, seqtk telo,
#   make_reference_table.py --censat --telo, then postprocess_sex_chrom.py --paf for males
# The alignment and the telomere calls depend only on the assembly and CHM13, so an existing
# workspace <hap>_ref.paf and <hap>.telo.bed are reused.
#
# Usage: ref_table.sh NORMAL ASSEMBLY_HAP1 ASSEMBLY_HAP2 SEX REFERENCE_DIR FULL_REFERENCE CHM13_CENSAT
#            WORK_DIR OUTPUT_DIR THREADS   (run inside the copynumber image)
# Output: OUTPUT_DIR/<normal>.<hap>.ref.table (+ .ref.table.log, the make_reference_table.py log)

set -xv
set -o errexit
set -o nounset
set -o pipefail

NORMAL=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
SEX=$4
REFERENCE_DIR=$5
FULL_REFERENCE=$6
CHM13_CENSAT=$7
WORK_DIR=$8
OUTPUT_DIR=$9
THREADS=${10:-16}

if [ ${SEX} = "female" ]; then
    REFERENCE=${REFERENCE_DIR}/chm13_female.fa
else
    REFERENCE=${FULL_REFERENCE}
fi

mkdir -p ${WORK_DIR} ${OUTPUT_DIR}

for hap in hap1 hap2; do
    if [ ${hap} = "hap1" ]; then ASSEMBLY=${ASSEMBLY_HAP1}; else ASSEMBLY=${ASSEMBLY_HAP2}; fi
    if [ ! -s ${WORK_DIR}/${hap}_ref.paf ]; then
        minimap2 -cx asm5 --secondary=no -t ${THREADS} ${ASSEMBLY} ${REFERENCE} > ${WORK_DIR}/${hap}_ref.paf
    fi
    if [ ! -s ${WORK_DIR}/${hap}.telo.bed ]; then
        seqtk telo ${ASSEMBLY} > ${WORK_DIR}/${hap}.telo.bed
    fi
    python3 ./scripts/copynumber/make_reference_table.py \
        -i ${WORK_DIR}/${hap}_ref.paf \
        --censat ${CHM13_CENSAT} \
        --telo ${WORK_DIR}/${hap}.telo.bed \
        > ${WORK_DIR}/ref_table_${hap}.raw \
        2> ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table.log
    cp ${WORK_DIR}/ref_table_${hap}.raw ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table
done

if [ ${SEX} = "male" ]; then
    python3 ./scripts/copynumber/postprocess_sex_chrom.py \
        --hap1 ${OUTPUT_DIR}/${NORMAL}.hap1.ref.table \
        --hap2 ${OUTPUT_DIR}/${NORMAL}.hap2.ref.table \
        --out1 ${OUTPUT_DIR}/${NORMAL}.hap1.ref.table \
        --out2 ${OUTPUT_DIR}/${NORMAL}.hap2.ref.table \
        --paf ${WORK_DIR}/hap1_ref.paf ${WORK_DIR}/hap2_ref.paf
fi

echo "[ref_table] done: ${NORMAL}"
