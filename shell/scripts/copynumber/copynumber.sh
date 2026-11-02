#!/bin/bash

# Copy-number analysis of one tumor on the reference tables of ref_table.sh:
#   1. windowed tumor/normal depth: samtools depth -a -Q 40 on the reference-table intervals,
#      copynumber_window.py
#   2. cbs.R (auto ploidy unless PLOIDY_HAP1/2 are given), split_gaps.py
#   3. plot_copy_number.R
#
# Usage: copynumber.sh TUMOR NORMAL LABEL TUMOR_BAM NORMAL_BAM SEX REFERENCE_DIR CHM13_CENSAT CENSAT_BED
#            HAP1_SATELLITE HAP2_SATELLITE BIN_WIDTH PLOIDY_HAP1 PLOIDY_HAP2 HAP1_LABEL HAP2_LABEL
#            WORK_DIR OUTPUT_DIR THREADS   (run inside the copynumber image)
# PLOIDY_HAP1 / PLOIDY_HAP2: "-" for automatic estimation (both must be given to override)
# CENSAT_BED: cenSat annotation of the assembly for the plot; "-" to use HAP1_SATELLITE + HAP2_SATELLITE

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
LABEL=$3
TUMOR_BAM=$4
NORMAL_BAM=$5
SEX=$6
REFERENCE_DIR=$7
CHM13_CENSAT=$8
CENSAT_BED=$9
HAP1_SATELLITE=${10}
HAP2_SATELLITE=${11}
BIN_WIDTH=${12}
PLOIDY_HAP1=${13}
PLOIDY_HAP2=${14}
HAP1_LABEL=${15}
HAP2_LABEL=${16}
WORK_DIR=${17}
OUTPUT_DIR=${18}
THREADS=${19:-16}

CHM13_LENGTHS=${REFERENCE_DIR}/chm13_${SEX}_chrom_length.txt

mkdir -p ${WORK_DIR}

for hap in hap1 hap2; do
    [ -s ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table ] || { echo "ERROR: run ref_table.sh first: missing ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table" >&2; exit 1; }
done

# 1. Windowed depth
for hap in hap1 hap2; do
    REF_TABLE=${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table
    awk '{print $1 "\t" $2 "\t" $3}' ${REF_TABLE} > ${WORK_DIR}/${NORMAL}.${hap}.ref.bed
    for entry in ${TUMOR}:${TUMOR_BAM} ${NORMAL}:${NORMAL_BAM}; do
        SAMPLE=${entry%%:*}
        BAM=${entry#*:}
        samtools depth -@ ${THREADS} -a -b ${WORK_DIR}/${NORMAL}.${hap}.ref.bed ${BAM} -Q 40 \
            | awk '{print $1 "\t" $2 - 1 "\t" $2 "\t" $3}' \
            | bgzip -@ ${THREADS} > ${WORK_DIR}/${SAMPLE}.${hap}.depth.bed.gz
        tabix -f -p bed ${WORK_DIR}/${SAMPLE}.${hap}.depth.bed.gz
    done
    python3 ./scripts/copynumber/copynumber_window.py \
        -t ${WORK_DIR}/${TUMOR}.${hap}.depth.bed.gz \
        -n ${WORK_DIR}/${NORMAL}.${hap}.depth.bed.gz \
        -r ${REF_TABLE} > ${OUTPUT_DIR}/${TUMOR}.${hap}.copynumber.tsv
done

# 2. CBS and gap splitting
for hap in hap1 hap2; do
    if [ ${PLOIDY_HAP1} != "-" ] && [ ${PLOIDY_HAP2} != "-" ]; then
        if [ ${hap} = "hap1" ]; then PLOIDY_FLAGS="-p ${PLOIDY_HAP1}"; else PLOIDY_FLAGS="-p ${PLOIDY_HAP2}"; fi
    else
        PLOIDY_FLAGS="-a"
    fi
    Rscript ./scripts/copynumber/cbs.R \
        -i ${OUTPUT_DIR}/${TUMOR}.${hap}.copynumber.tsv \
        -s ${LABEL} \
        -o ${OUTPUT_DIR}/${TUMOR}.${hap}.cbs.txt \
        ${PLOIDY_FLAGS} \
        -w ${BIN_WIDTH} \
        --ploidy-out ${OUTPUT_DIR}/${TUMOR}.${hap}.ploidy 2>&1 | tee ${OUTPUT_DIR}/${TUMOR}.${hap}.cbs.log
    python3 ./scripts/copynumber/split_gaps.py \
        ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table \
        ${OUTPUT_DIR}/${TUMOR}.${hap}.cbs.txt \
        > ${OUTPUT_DIR}/${TUMOR}.${hap}.cbs.split.txt
done

# 3. Plot
if [ ${CENSAT_BED} != "-" ]; then
    ANNOTATION=${CENSAT_BED}
else
    zcat ${HAP1_SATELLITE} ${HAP2_SATELLITE} | bgzip -c > ${WORK_DIR}/${NORMAL}.satellite.annotation.bed.gz
    ANNOTATION=${WORK_DIR}/${NORMAL}.satellite.annotation.bed.gz
fi

Rscript ./scripts/copynumber/plot_copy_number.R \
    -i ${OUTPUT_DIR}/${TUMOR}.hap1.copynumber.tsv \
    -j ${OUTPUT_DIR}/${TUMOR}.hap2.copynumber.tsv \
    -k ${OUTPUT_DIR}/${TUMOR}.hap1.cbs.split.txt \
    -l ${OUTPUT_DIR}/${TUMOR}.hap2.cbs.split.txt \
    -o ${OUTPUT_DIR}/${TUMOR}.copynumber.png \
    -s ${LABEL} \
    -p $(< ${OUTPUT_DIR}/${TUMOR}.hap1.ploidy) \
    -t $(< ${OUTPUT_DIR}/${TUMOR}.hap2.ploidy) \
    -q ${OUTPUT_DIR}/${NORMAL}.hap1.ref.table \
    -r ${OUTPUT_DIR}/${NORMAL}.hap2.ref.table \
    -a ${ANNOTATION} \
    -w ${BIN_WIDTH} \
    --hap1_label ${HAP1_LABEL} \
    --hap2_label ${HAP2_LABEL} \
    --plot-sex-chrom TRUE \
    --chm13_censat ${CHM13_CENSAT} \
    --chm13_lengths ${CHM13_LENGTHS}

echo "[copynumber] done: ${TUMOR}"
