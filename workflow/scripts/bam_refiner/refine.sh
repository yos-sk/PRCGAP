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
# bam_refiner parallelises internally, so the read-name-sorted BAM no longer has
# to be split into ${THREAD} files and merged back: one call replaces the
# split_bam + for-loop + samtools merge trio (and the extra copy of the BAM on
# disk that came with it).
bam_refiner refine \
    --input-bam ${OUTPUT_BAM_PREFIX}.bam \
    --output-bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam \
    --hap1-tabix ${KMER_DIR}/hap1_cnt_kmerposition.bed.gz \
    --hap2-tabix ${KMER_DIR}/hap2_cnt_kmerposition.bed.gz \
    --hap1-list ${KMER_DIR}/hap1_list.txt.gz \
    --hap2-list ${KMER_DIR}/hap2_list.txt.gz \
    --kmer-size 21 \
    --threads ${THREAD} \
    1>${OUTPUT_DIR}/bam_refiner_result.tsv 2>${OUTPUT_DIR}/bam_refiner.log

samtools sort -@ ${THREAD} -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam
samtools index -@ ${THREAD} ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam
rm ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam

gzip -f ${OUTPUT_DIR}/bam_refiner_result.tsv
gzip -f ${OUTPUT_DIR}/bam_refiner.log

# Reads spanning no haplotype-specific locus have no ratio to report, so they get
# this prior instead of a bogus 1.0. The values are picked to land just above the
# 0.6 Kmer_ratio cutoff applied downstream, which is what keeps those reads in play;
# leaving the option off would fall back to 0.5 and put every one of them on the
# discarded side. ONT is noisier, so its no-evidence reads get less credit.
if [ "${DATA}" = "hifi" ]; then
    PRIOR_MEAN=0.8
else
    PRIOR_MEAN=0.6
fi

bam_refiner kmer-ratio \
    ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam \
    --threads ${THREAD} \
    --prior-mean ${PRIOR_MEAN} \
> ${OUTPUT_DIR}/${SAMPLE}_kmer_ratio.txt

rm -rf ${WORK_DIR}

echo ${?}
