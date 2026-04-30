#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP1_CONTIG=$2
HAP2_CONTIG=$3
OUTPUT_DIR=$4
THREAD=${5:-8}

# Step 2 + Step 3 of run_refine.sh — depend only on the assembly,
# so they are computed once per sample and shared across HiFi/ONT.

WORK_DIR=${OUTPUT_DIR}/workspace
mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}

# Step 2: Extract haplotype-specific unique k-mers
mkdir -p ${WORK_DIR}/meryl
meryl count k=21 threads=${THREAD} ${HAP1_CONTIG} output ${WORK_DIR}/meryl/hap1.meryl
meryl count k=21 threads=${THREAD} ${HAP2_CONTIG} output ${WORK_DIR}/meryl/hap2.meryl

meryl difference ${WORK_DIR}/meryl/hap1.meryl ${WORK_DIR}/meryl/hap2.meryl output ${WORK_DIR}/meryl/hap1.uniq.meryl
meryl difference ${WORK_DIR}/meryl/hap2.meryl ${WORK_DIR}/meryl/hap1.meryl output ${WORK_DIR}/meryl/hap2.uniq.meryl

meryl print threads=${THREAD} ${WORK_DIR}/meryl/hap1.uniq.meryl > ${WORK_DIR}/meryl/hap1.uniq.tsv
meryl print threads=${THREAD} ${WORK_DIR}/meryl/hap2.uniq.meryl > ${WORK_DIR}/meryl/hap2.uniq.tsv

for hap in hap1 hap2
do
    awk '{if ($2 == 1) print}' ${WORK_DIR}/meryl/${hap}.uniq.tsv > ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv
    gzip -f ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv
    gzip -f ${WORK_DIR}/meryl/${hap}.uniq.tsv
done

bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/hap1.cnt.uniq.tsv.gz \
    -f ${HAP1_CONTIG} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed.gz

bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/hap2.cnt.uniq.tsv.gz \
    -f ${HAP2_CONTIG} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed.gz

# Step 3: List contig names per haplotype
grep ">" ${HAP1_CONTIG} | sed 's/>//' > ${OUTPUT_DIR}/hap1_list.txt
gzip -f ${OUTPUT_DIR}/hap1_list.txt
grep ">" ${HAP2_CONTIG} | sed 's/>//' > ${OUTPUT_DIR}/hap2_list.txt
gzip -f ${OUTPUT_DIR}/hap2_list.txt

rm -rf ${WORK_DIR}

echo ${?}
