#!/bin/bash
# Cut the genome-wide references down to chr20 for the HG008 chr20 test case.
#
# The assembly is chr20 only, so aligning it against 3 Gb of reference is waste:
# chain_files_hap ran over an hour per haplotype/reference and liftoff_hap 41
# minutes, both dominated by the minimap2 index. Restricting the reference makes
# the index ~1/50.
#
# Writes into resources/reference/, which test_configure_{full,min}.sh point at
# via --chm13-fasta / --grch38-fasta / --grch38-gtf. The genome-wide originals
# stay in the repo-root resource/reference/ and are read from there.
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

# The destination used to hold symlinks INTO the source, so writing chr20 over a
# name whose symlink still pointed at the genome-wide original truncated the
# original to 0 bytes. Refuse to run unless the two trees are distinct, and read
# through to the real files.
SRC="$(cd "${SRC}" && pwd -P)"
if [ "$(cd "${DST}" && pwd -P)" = "${SRC}" ]; then
    echo "error: source and destination are the same directory (${SRC})" >&2
    exit 1
fi
for f in chm13v2.0_maskedY_rCRS.fa GRCh38.d1.vd1.fa \
         Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf centromeres.txt.gz \
         GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed chm13v2.0_censat_v2.1.bed.gz; do
    if [ -e "${DST}/${f}" ] && [ "$(readlink -f "${DST}/${f}")" = "$(readlink -f "${SRC}/${f}")" ]; then
        echo "error: ${DST}/${f} resolves to the source file; remove it first" >&2
        exit 1
    fi
done

echo "### source ${SRC}"
echo "### dest   ${DST}"

# ---- CHM13 chr20 ----
samtools faidx "${SRC}/chm13v2.0_maskedY_rCRS.fa" chr20 \
    > "${DST}/chm13v2.0_maskedY_rCRS.fa"
samtools faidx "${DST}/chm13v2.0_maskedY_rCRS.fa"

# ---- GRCh38 chr20 ----
samtools faidx "${SRC}/GRCh38.d1.vd1.fa" chr20 \
    > "${DST}/GRCh38.d1.vd1.fa"
samtools faidx "${DST}/GRCh38.d1.vd1.fa"

# ---- GRCh38 Ensembl GTF, chr20 ----
awk -F'\t' '$1 == "chr20"' "${SRC}/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf" \
    > "${DST}/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf"

# ---- BEDs: chr20 rows only ----
# centromeres.txt.gz is a UCSC table (bin, chrom, start, end, name), so chrom is
# field 2, not 1.
zcat "${SRC}/centromeres.txt.gz" | awk -F'\t' '$2 == "chr20"' | gzip -c \
    > "${DST}/centromeres.txt.gz"
awk -F'\t' '$1 == "chr20"' "${SRC}/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed" \
    > "${DST}/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed"
zcat "${SRC}/chm13v2.0_censat_v2.1.bed.gz" | awk -F'\t' '$1 == "chr20"' \
    > "${DST}/chm13v2.0_censat_v2.1.bed"
bgzip -f "${DST}/chm13v2.0_censat_v2.1.bed"
tabix -f -p bed "${DST}/chm13v2.0_censat_v2.1.bed.gz"

echo
for f in chm13v2.0_maskedY_rCRS.fa GRCh38.d1.vd1.fa \
         Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf \
         centromeres.txt.gz GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed \
         chm13v2.0_censat_v2.1.bed.gz; do
    printf "  %-48s %s\n" "$f" "$(du -h "${DST}/${f}" | cut -f1)"
done
echo "[extract_chr20_reference] done"
