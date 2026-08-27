#!/bin/bash
# Build the masked CHM13 / GRCh38 references the chain-file step aligns
# against. Satellite and centromere arrays align promiscuously, so masking
# them keeps minimap2 asm5 from emitting chains through them.
# Runs inside the chain_files singularity container (bedtools + samtools +
# python3).
#
# Emits, in OUTPUT_DIR, the pair the sample's sex selects (+ .fai for each):
#   DROP_CHRY=true   chm13.masked_noY.fa  GRCh38.masked_noY.fa
#   DROP_CHRY=false  chm13.masked.fa      GRCh38.masked.fa
# Only one pair is ever read; building both leaves a full extra copy of each
# reference on disk.
#
# Required positional args:
#   $1  CHM13_FASTA        CHM13v2.0 FASTA
#   $2  CHM13_CENSAT       CHM13 cenSat BED (plain or .gz)
#   $3  GRCH38_FASTA       GRCh38 FASTA
#   $4  GRCH38_CENTROMERES UCSC hg38 centromeres.txt(.gz)
#   $5  GRCH38_EXCLUSIONS  GRC exclusion regions BED
#   $6  OUTPUT_DIR         destination dir for the masked FASTAs
#   $7  WORK_DIR           scratch dir for the mask BEDs
#   $8  SCRIPT_DIR         absolute path to workflow/scripts/annotation
#   $9  DROP_CHRY          true (female) | false (male)

set -xv
set -o errexit
set -o nounset
set -o pipefail

CHM13_FASTA=$1
CHM13_CENSAT=$2
GRCH38_FASTA=$3
GRCH38_CENTROMERES=$4
GRCH38_EXCLUSIONS=$5
OUTPUT_DIR=$6
WORK_DIR=$7
SCRIPT_DIR=$8
DROP_CHRY=${9:-true}

mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

# ---- CHM13: mask HOR / monomeric / HSat arrays from the cenSat BED ----
zcat -f "${CHM13_CENSAT}" \
    | grep -e hor -e mon -e hsat \
    | awk 'NR != 1 {print $1 "\t" $2 "\t" $3}' > "${WORK_DIR}/chm13_mask_regions.bed"

if [ "${DROP_CHRY}" = "true" ]; then
    bedtools maskfasta \
        -fi "${CHM13_FASTA}" \
        -fo "${WORK_DIR}/chm13.masked.withY.fa" \
        -bed "${WORK_DIR}/chm13_mask_regions.bed"
    awk '/^>/ {p = ($0 !~ /^>chrY/)} p' "${WORK_DIR}/chm13.masked.withY.fa" \
        > "${OUTPUT_DIR}/chm13.masked_noY.fa"
    rm -f "${WORK_DIR}/chm13.masked.withY.fa"
    samtools faidx "${OUTPUT_DIR}/chm13.masked_noY.fa"
else
    bedtools maskfasta \
        -fi "${CHM13_FASTA}" \
        -fo "${OUTPUT_DIR}/chm13.masked.fa" \
        -bed "${WORK_DIR}/chm13_mask_regions.bed"
    samtools faidx "${OUTPUT_DIR}/chm13.masked.fa"
fi

# ---- GRCh38: mask centromeres + GRC exclusion regions, drop non-chromosomes ----
# centromeres.txt is a UCSC table dump: bin, chrom, chromStart, chromEnd, name.
zcat -f "${GRCH38_CENTROMERES}" \
    | awk '{print $2 "\t" $3 - 1 "\t" $4}' > "${WORK_DIR}/GRCh38_mask_regions.bed"
cat "${GRCH38_EXCLUSIONS}" >> "${WORK_DIR}/GRCh38_mask_regions.bed"

python3 "${SCRIPT_DIR}/remove_unlocalized_GRCh38.py" "${GRCH38_FASTA}" \
    > "${WORK_DIR}/GRCh38_removed_unlocalized.fa"

if [ "${DROP_CHRY}" = "true" ]; then
    bedtools maskfasta \
        -fi "${WORK_DIR}/GRCh38_removed_unlocalized.fa" \
        -fo "${WORK_DIR}/GRCh38.masked.withY.fa" \
        -bed "${WORK_DIR}/GRCh38_mask_regions.bed"
    awk '/^>/ {p = ($0 !~ /^>chrY/)} p' "${WORK_DIR}/GRCh38.masked.withY.fa" \
        > "${OUTPUT_DIR}/GRCh38.masked_noY.fa"
    rm -f "${WORK_DIR}/GRCh38.masked.withY.fa"
    samtools faidx "${OUTPUT_DIR}/GRCh38.masked_noY.fa"
else
    bedtools maskfasta \
        -fi "${WORK_DIR}/GRCh38_removed_unlocalized.fa" \
        -fo "${OUTPUT_DIR}/GRCh38.masked.fa" \
        -bed "${WORK_DIR}/GRCh38_mask_regions.bed"
    samtools faidx "${OUTPUT_DIR}/GRCh38.masked.fa"
fi

echo "[prepare_mask_regions] done"
