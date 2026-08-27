#!/bin/bash
# Cut the genome-wide references down to chr20 for the HG008 chr20 test case.
#
# The assembly is chr20 only, so aligning it against 3 Gb of reference is waste:
# the minimap2 index dominates chain_files_hap and liftoff_hap, and restricting
# the reference makes that index roughly 1/50 the size.
#
# Reads from the repo-root resource/reference/, which resource/scripts/
# download_reference.sh fills once for every run in the repository. Writes into
# this test case's resources/reference/ under names carrying a _chr20 suffix, so
# a chr20 cut is never mistaken for the genome-wide original it came from.
# test_configure.sh points at the suffixed names.
#
# The GTF has to be cut as well: with a chr20-only FASTA, liftoff would still try
# to lift every other chromosome's genes and find no sequence for them.
#
# Usage: bash extract_chr20_reference.sh [SRC_DIR]
#   SRC_DIR  genome-wide references (default: ../../../../resource/reference)

set -o errexit
set -o nounset
set -o pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
SRC="${1:-$(cd "${HERE}/../../../../resource/reference" && pwd)}"
DST="$(cd "${HERE}/.." && pwd)/reference"
mkdir -p "${DST}"

SRC="$(cd "${SRC}" && pwd -P)"
if [ "$(cd "${DST}" && pwd -P)" = "${SRC}" ]; then
    echo "error: source and destination are the same directory (${SRC})" >&2
    exit 1
fi

echo "### source ${SRC}"
echo "### dest   ${DST}"

# ---- CHM13 chr20 ----
samtools faidx "${SRC}/chm13v2.0_maskedY_rCRS.fa" chr20 \
    > "${DST}/chm13v2.0_maskedY_rCRS_chr20.fa"
samtools faidx "${DST}/chm13v2.0_maskedY_rCRS_chr20.fa"

# ---- GRCh38 chr20 ----
samtools faidx "${SRC}/GRCh38.d1.vd1.fa" chr20 \
    > "${DST}/GRCh38.d1.vd1_chr20.fa"
samtools faidx "${DST}/GRCh38.d1.vd1_chr20.fa"

# ---- GRCh38 Ensembl GTF, chr20 ----
awk -F'\t' '$1 == "chr20"' "${SRC}/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf" \
    > "${DST}/Homo_sapiens.GRCh38.Ensembl.112.chr.format_chr20.gtf"

# ---- BEDs: chr20 rows only ----
# centromeres.txt.gz is a UCSC table (bin, chrom, start, end, name), so chrom is
# field 2, not 1.
zcat "${SRC}/centromeres.txt.gz" | awk -F'\t' '$2 == "chr20"' | gzip -c \
    > "${DST}/centromeres_chr20.txt.gz"
awk -F'\t' '$1 == "chr20"' "${SRC}/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed" \
    > "${DST}/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2_chr20.bed"
zcat "${SRC}/chm13v2.0_censat_v2.1.bed.gz" | awk -F'\t' '$1 == "chr20"' \
    > "${DST}/chm13v2.0_censat_v2.1_chr20.bed"
bgzip -f "${DST}/chm13v2.0_censat_v2.1_chr20.bed"
tabix -f -p bed "${DST}/chm13v2.0_censat_v2.1_chr20.bed.gz"

echo
for f in chm13v2.0_maskedY_rCRS_chr20.fa GRCh38.d1.vd1_chr20.fa \
         Homo_sapiens.GRCh38.Ensembl.112.chr.format_chr20.gtf \
         centromeres_chr20.txt.gz \
         GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2_chr20.bed \
         chm13v2.0_censat_v2.1_chr20.bed.gz; do
    printf "  %-54s %s\n" "$f" "$(du -h "${DST}/${f}" | cut -f1)"
done
echo "[extract_chr20_reference] done"
