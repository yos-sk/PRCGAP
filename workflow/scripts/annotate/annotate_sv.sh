#!/bin/bash
# SV annotation chain (per tumor, per seqtype) for the PRCGAP workflow.
#
# Designed to be run inside the point_mutation_postprocess singularity
# container (snakemake's `singularity:` directive on the rule wraps the
# shell automatically). The coordconv step has been factored out into
# its own rule, so this script does NOT invoke coordconv itself — it
# accepts pre-computed coordconv BEDs as args (or empty strings to skip
# the corresponding conv annotation step).
#
# Required positional args:
#   $1  SAMPLE         tumor sample name
#   $2  PASS_TXT       PASS-filtered nanomonsv TSV (from prep_sv rule)
#   $3  NANOMONSV_OTHER  the OTHER seqtype's insert_classified.txt
#   $4  SUPPORT_READS  this seqtype's *.nanomonsv.supporting_read.txt
#   $5  KMER_RATIO     this seqtype's *_kmer_ratio.txt
#   $6  HAP1_FA        haplotype 1 assembly fasta
#   $7  HAP2_FA        haplotype 2 assembly fasta
#   $8  OUTPUT_DIR     directory for per-seqtype annotation outputs
#   $9  WORK_DIR       intermediate dir
#   $10 SCRIPT_DIR     absolute path to workflow/scripts/annotate
#
# Optional positional args (pass "" to skip the corresponding step):
#   $11 LIFTOFF_BED           tabix-indexed liftoff gene BED (.bed.gz + .tbi)
#                              (derived from the user-supplied liftoff GFF
#                              by the upstream `gff_to_bed` rule)
#   $12 CGC_TSV               cancer_gene_census tsv
#   $13 RMSK_BED              tabix-indexed RepeatMasker BED
#   $14 CENSAT_BED            tabix-indexed centromere BED
#   $15 SEGDUP_BED            tabix-indexed segdup BED
#   $16 MISASSEMBLY_HAP1_BED  hap1 misassembly BED
#   $17 MISASSEMBLY_HAP2_BED  hap2 misassembly BED
#   $18 GRCH38_BED            pre-computed coordconv-to-GRCh38 BED
#   $19 CHM13_BED             pre-computed coordconv-to-chm13 BED
#   $20 GNOMAD_BED            tabix-indexed gnomAD SV bed (requires GRCH38_BED)

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
PASS_TXT=$2
NANOMONSV_OTHER=$3
SUPPORT_READS=$4
KMER_RATIO=$5
HAP1_FA=$6
HAP2_FA=$7
OUTPUT_DIR=$8
WORK_DIR=$9
SCRIPT_DIR=${10}

LIFTOFF_BED=${11:-}
CGC_TSV=${12:-}
RMSK_BED=${13:-}
CENSAT_BED=${14:-}
SEGDUP_BED=${15:-}
MISASSEMBLY_HAP1_BED=${16:-}
MISASSEMBLY_HAP2_BED=${17:-}
GRCH38_BED=${18:-}
CHM13_BED=${19:-}
GNOMAD_BED=${20:-}

mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

ADD="python3 ${SCRIPT_DIR}/add_annotation_sv.py"
PREFIX="${SAMPLE}.PRCGAP"

CUR="${PASS_TXT}"
NEXT=""

# 1. gene annotation — pass the tabix-indexed liftoff BED (derived by
#    upstream `gff_to_bed` rule), matching reference
#    SV_analysis/scripts/PRCGAP/annotate_PRCGAP.sh:148 which passes
#    `liftoff.bed.gz` directly to `add_annotation.py gene -b`.
if [ -n "${LIFTOFF_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.gene.txt"
    ${ADD} gene -i "${CUR}" -o "${NEXT}" -b "${LIFTOFF_BED}" ${CGC_TSV:+-c "${CGC_TSV}"}
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

# 4/5. liftover-able to GRCh38 / chm13 (conv) — uses pre-computed coordconv BEDs
if [ -n "${GRCH38_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.conv_GRCh38.txt"
    ${ADD} conv -i "${CUR}" -o "${NEXT}" -l "${GRCH38_BED}" --feature liftover_GRCh38
    CUR="${NEXT}"
fi
if [ -n "${CHM13_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.conv_chm13.txt"
    ${ADD} conv -i "${CUR}" -o "${NEXT}" -l "${CHM13_BED}" --feature liftover_chm13
    CUR="${NEXT}"
fi

# 6. kmer ratio
NEXT="${WORK_DIR}/${PREFIX}.kmer.txt"
${ADD} kmer -i "${CUR}" -t "${SUPPORT_READS}" -m "${KMER_RATIO}" -o "${NEXT}"
CUR="${NEXT}"

# 7. centromere
if [ -n "${CENSAT_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.cen.txt"
    ${ADD} cen -i "${CUR}" -o "${NEXT}" -s "${CENSAT_BED}"
    CUR="${NEXT}"
fi

# 8. segdup
if [ -n "${SEGDUP_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.segdup.txt"
    ${ADD} segdup -i "${CUR}" -o "${NEXT}" -s "${SEGDUP_BED}"
    CUR="${NEXT}"
fi

# 9. nanomonsv other (HiFi vs ONT cross-check; skipped when the other seqtype
#    is unavailable, i.e. HiFi-only / ONT-only samples)
if [ -n "${NANOMONSV_OTHER}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.other.txt"
    ${ADD} other -i "${CUR}" -o "${NEXT}" -j "${NANOMONSV_OTHER}"
    CUR="${NEXT}"
fi

# 10. misassembly (optional)
if [ -n "${MISASSEMBLY_HAP1_BED}" ] && [ -n "${MISASSEMBLY_HAP2_BED}" ]; then
    NEXT="${WORK_DIR}/${PREFIX}.misassembly.txt"
    ${ADD} misassembly -i "${CUR}" -o "${NEXT}" -b "${MISASSEMBLY_HAP1_BED}" -c "${MISASSEMBLY_HAP2_BED}"
    CUR="${NEXT}"
fi

# Stable pre-gnomAD result name (consumed by reclassify_sv)
PRE_GNOMAD="${OUTPUT_DIR}/${PREFIX}.nanomonsv_results.annotated.txt"
cp "${CUR}" "${PRE_GNOMAD}"

# 11. gnomAD (optional, requires GRCh38 coordconv output)
if [ -n "${GNOMAD_BED}" ] && [ -n "${GRCH38_BED}" ]; then
    GNOMAD_OUT="${OUTPUT_DIR}/${PREFIX}.nanomonsv_results.annotated_gnomad.txt"
    ${ADD} gnomad -i "${PRE_GNOMAD}" -o "${GNOMAD_OUT}" -g "${GNOMAD_BED}"
fi

echo "[annotate_sv] done: ${SAMPLE}"
