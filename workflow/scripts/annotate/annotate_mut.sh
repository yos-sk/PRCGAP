#!/bin/bash
# SNV / INDEL annotation chain (per tumor, per tool ∈ {clairs, deepsomatic},
# per mode ∈ {snv, indel}).
#
# Runs inside the point_mutation_postprocess singularity container.
#
# Mirrors snv/annotate_snv.sh + indel/annotate_indel.sh from the
# reference scripts. Coordinate liftover is bundled in by the upstream
# coordconv (SNV) or transanno (INDEL) rules; the GRCh38/CHM13 *variant
# comparison* steps (compare_current_reference / filter_diff_references
# / check_homozygous) are intentionally skipped per proposal step 3.
#
# Required positional args:
#   $1  SAMPLE       tumor sample name
#   $2  TOOL         clairs | deepsomatic
#   $3  MODE         snv | indel
#   $4  PREP_BED     14-col haplotyped.bed filtered to SNV or INDEL rows
#   $5  OTHER_VCF    other tool's vcf for the "other" cross-check step
#   $6  HAP1_FA
#   $7  HAP2_FA
#   $8  OUTPUT_DIR
#   $9  WORK_DIR
#   $10 SCRIPT_DIR
#
# Optional positional args (pass "" to skip):
#   $11 LIFTOFF_GFF  tabix-indexed liftoff GFF (.gff.gz + .tbi)
#   $12 CGC_TSV
#   $13 CMRG_TSV
#   $14 GENCODE_BED
#   $15 RMSK_BED
#   $16 CENSAT_BED
#   $17 SEGDUP_BED
#   $18 MISASSEMBLY_HAP1_BED
#   $19 MISASSEMBLY_HAP2_BED
#   $20 GNOMAD_VCF
#   $21 GRCH38_LIFT   coordconv output (SNV) or transanno VCF (INDEL)
#   $22 CHM13_LIFT    coordconv output (SNV) or transanno VCF (INDEL)

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
TOOL=$2
MODE=$3
PREP_BED=$4
OTHER_VCF=$5
HAP1_FA=$6
HAP2_FA=$7
OUTPUT_DIR=$8
WORK_DIR=$9
SCRIPT_DIR=${10}

LIFTOFF_GFF=${11:-}
CGC_TSV=${12:-}
CMRG_TSV=${13:-}
GENCODE_BED=${14:-}
RMSK_BED=${15:-}
CENSAT_BED=${16:-}
SEGDUP_BED=${17:-}
MISASSEMBLY_HAP1_BED=${18:-}
MISASSEMBLY_HAP2_BED=${19:-}
GNOMAD_VCF=${20:-}
GRCH38_LIFT=${21:-}
CHM13_LIFT=${22:-}

mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

ADD="python3 ${SCRIPT_DIR}/add_annotation_mut.py"

PREFIX="${SAMPLE}.${TOOL}.${MODE}"

# 0. Add header + GRCh38/CHM13 lifted-coordinate columns. This is the
#    coordinate-annotation portion of compare_current_reference.py
#    (the GRCh38/CHM13 variant comparison part is deferred).
CUR="${WORK_DIR}/${PREFIX}.lifted.txt"
python3 "${SCRIPT_DIR}/add_lift_coords.py" \
    --mode "${MODE}" \
    -i "${PREP_BED}" \
    -o "${CUR}" \
    --grch38 "${GRCH38_LIFT}" \
    --chm13 "${CHM13_LIFT}"

# 1. gene
if [ -n "${LIFTOFF_GFF}" ] && [ -n "${CGC_TSV}" ] && [ -n "${CMRG_TSV}" ] && [ -n "${GENCODE_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.gene.txt"
    ${ADD} gene -i "${CUR}" -o "${NEXT}" -g "${LIFTOFF_GFF}" -c "${CGC_TSV}" -m "${CMRG_TSV}" -t "${GENCODE_BED}"
    CUR="${NEXT}"
fi

# 2. RepeatMasker
if [ -n "${RMSK_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.rmsk.txt"
    ${ADD} rmsk -i "${CUR}" -o "${NEXT}" -b "${RMSK_BED}"
    CUR="${NEXT}"
fi

# 3. Contig size
NEXT="${WORK_DIR}/${PREFIX}.size.txt"
${ADD} size -i "${CUR}" -o "${NEXT}" -f "${HAP1_FA}" -g "${HAP2_FA}"
CUR="${NEXT}"

# 4. misassembly (optional)
if [ -n "${MISASSEMBLY_HAP1_BED}" ] && [ -n "${MISASSEMBLY_HAP2_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.misassembly.txt"
    ${ADD} misassembly -i "${CUR}" -o "${NEXT}" -b "${MISASSEMBLY_HAP1_BED}" -c "${MISASSEMBLY_HAP2_BED}"
    CUR="${NEXT}"
fi

# 5. centromere
if [ -n "${CENSAT_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.cen.txt"
    ${ADD} cen -i "${CUR}" -o "${NEXT}" -s "${CENSAT_BED}"
    CUR="${NEXT}"
fi

# 6. segdup
if [ -n "${SEGDUP_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.segdup.txt"
    ${ADD} segdup -i "${CUR}" -o "${NEXT}" -s "${SEGDUP_BED}"
    CUR="${NEXT}"
fi

# 7. other (cross-check with the other tool's vcf)
NEXT="${WORK_DIR}/${PREFIX}.other.txt"
${ADD} other -i "${CUR}" -o "${NEXT}" -j "${OTHER_VCF}"
CUR="${NEXT}"

# 8. gnomAD (final step in reference snv/indel scripts). Skipped if no
#    gnomAD VCF configured or if no GRCh38 chain was configured (the
#    annotation looks up GRCh38_pos which won't exist without liftover).
if [ -n "${GNOMAD_VCF}" ] && [ -n "${GRCH38_LIFT}" ]; then
    NEXT="${OUTPUT_DIR}/${PREFIX}.annotated.txt"
    ${ADD} gnomad -i "${CUR}" -o "${NEXT}" -k "${GNOMAD_VCF}"
else
    cp "${CUR}" "${OUTPUT_DIR}/${PREFIX}.annotated.txt"
fi

echo "[annotate_mut] done: ${SAMPLE}/${TOOL}/${MODE}"
