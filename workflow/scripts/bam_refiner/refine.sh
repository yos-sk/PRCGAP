#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP1_CONTIG=$2
HAP2_CONTIG=$3
INPUT=$4              # Single file or comma-separated list of files (BAM or FASTQ)
OUTPUT_DIR=$5
THREAD=$6
DATA=$7               # hifi or ont
KMER_DIR=$8           # Directory with precomputed k-mer outputs from kmer_extract.sh

# Steps 1 + 4 of run_refine.sh: mapping + refine using precomputed k-mer files.

mkdir -p ${OUTPUT_DIR}
WORK_DIR=${OUTPUT_DIR}/workspace
mkdir -p ${WORK_DIR}

# Concatenate haplotype contigs to form mapping reference
cat ${HAP1_CONTIG} ${HAP2_CONTIG} > ${WORK_DIR}/reference.fa

# ---------- Resolve input: support single or comma-separated list ----------
IFS=',' read -ra INPUT_FILES <<< "${INPUT}"

if [ ${#INPUT_FILES[@]} -eq 1 ]; then
    SINGLE_INPUT="${INPUT_FILES[0]}"
    EXT="${SINGLE_INPUT##*.}"
    if [ "${EXT}" = "bam" ]; then
        INPUT_BAM="${SINGLE_INPUT}"
        INPUT_FASTQ=""
    else
        INPUT_BAM=""
        INPUT_FASTQ="${SINGLE_INPUT}"
    fi
else
    FIRST_EXT="${INPUT_FILES[0]##*.}"
    MERGED_DIR="${WORK_DIR}/merged_input"
    mkdir -p ${MERGED_DIR}
    if [ "${FIRST_EXT}" = "bam" ]; then
        INPUT_BAM="${MERGED_DIR}/${SAMPLE}_merged.bam"
        samtools merge -@ ${THREAD} -f ${INPUT_BAM} ${INPUT_FILES[@]}
        samtools index -@ ${THREAD} ${INPUT_BAM}
        INPUT_FASTQ=""
    else
        INPUT_FASTQ="${MERGED_DIR}/${SAMPLE}_merged.fastq.gz"
        cat ${INPUT_FILES[@]} > ${INPUT_FASTQ}
        INPUT_BAM=""
    fi
fi

# ---------- Step 1: Mapping ----------
OUTPUT_BAM_PREFIX=${WORK_DIR}/${SAMPLE}
if [ "${DATA}" = "hifi" ]; then
    MM2_PRESET="asm5"
else
    MM2_PRESET="asm10"
fi

if [ -n "${INPUT_FASTQ}" ]; then
    minimap2 -t ${THREAD} -ax ${MM2_PRESET} --MD -y ${WORK_DIR}/reference.fa ${INPUT_FASTQ} \
        | samtools view -@ ${THREAD} -Shb - > ${OUTPUT_BAM_PREFIX}.unsorted
else
    samtools fastq -@ ${THREAD} -TMM,ML ${INPUT_BAM} \
        | minimap2 -t ${THREAD} -ax ${MM2_PRESET} --MD -y ${WORK_DIR}/reference.fa - \
        | samtools view -@ ${THREAD} -Shb - > ${OUTPUT_BAM_PREFIX}.unsorted
fi

samtools sort -@ ${THREAD} -m 2G -n ${OUTPUT_BAM_PREFIX}.unsorted -o ${OUTPUT_BAM_PREFIX}.bam
rm ${OUTPUT_BAM_PREFIX}.unsorted

# ---------- Step 4: Refine BAM using precomputed k-mer files ----------
mkdir -p ${WORK_DIR}/split
SIZE=$(split_bam size --input-file ${OUTPUT_BAM_PREFIX}.bam)
split_bam split \
    --input-file ${OUTPUT_BAM_PREFIX}.bam \
    --output-dir ${WORK_DIR}/split \
    --input-size ${SIZE} \
    --num-split ${THREAD}

for i in $(seq 0 $(( ${THREAD} - 1 ))); do
    bam_refiner refine \
        --input-bam ${WORK_DIR}/split/${i}.bam \
        --output-bam ${WORK_DIR}/split/${i}.refined.bam \
        --hap1-tabix ${KMER_DIR}/hap1_cnt_kmerposition.bed.gz \
        --hap2-tabix ${KMER_DIR}/hap2_cnt_kmerposition.bed.gz \
        --hap1-list ${KMER_DIR}/hap1_list.txt.gz \
        --hap2-list ${KMER_DIR}/hap2_list.txt.gz \
        --kmer-size 21 \
        1>${WORK_DIR}/split/${i}.bam_refiner.tsv 2>${WORK_DIR}/split/${i}.bam_refiner.log &
done
wait

cat ${WORK_DIR}/split/*.bam_refiner.tsv > ${OUTPUT_DIR}/bam_refiner_result.tsv
cat ${WORK_DIR}/split/*.bam_refiner.log > ${OUTPUT_DIR}/bam_refiner.log
samtools merge -@ ${THREAD} -f -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam ${WORK_DIR}/split/*.refined.bam

samtools sort -@ ${THREAD} -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam
samtools index -@ ${THREAD} ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam
rm ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam

gzip -f ${OUTPUT_DIR}/bam_refiner_result.tsv
gzip -f ${OUTPUT_DIR}/bam_refiner.log

bam_refiner kmer-ratio \
    ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam \
    --threads ${THREAD} \
> ${OUTPUT_DIR}/${SAMPLE}_kmer_ratio.txt

rm -rf ${WORK_DIR}

echo ${?}
