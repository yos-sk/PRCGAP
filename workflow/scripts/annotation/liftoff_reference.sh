#!/bin/bash
# Stage the GRCh38 FASTA + GTF that the per-haplotype liftoff jobs share.
#
# A female assembly has no chrY, so chrY is dropped from both — otherwise
# liftoff spends the whole run trying to place chrY features. For a male sample
# the inputs are used as-is and only symlinked, so nothing is copied.
# Runs inside the liftoff singularity container.
#
# Required positional args:
#   $1  SEX          female | male
#   $2  GRCH38       GRCh38 reference FASTA
#   $3  GRCH38_GTF   GRCh38 GTF with chr* contig names (plain .gtf)
#   $4  OUT_FASTA    staged FASTA path
#   $5  OUT_GTF      staged GTF path

set -xv
set -o errexit
set -o nounset
set -o pipefail

SEX=$1
GRCH38=$2
GRCH38_GTF=$3
OUT_FASTA=$4
OUT_GTF=$5

mkdir -p "$(dirname "${OUT_FASTA}")"

if [ "${SEX}" = "female" ]; then
    awk '/^>/ {p = ($0 !~ /^>chrY/)} p' "${GRCH38}" > "${OUT_FASTA}"
    grep -v "chrY" "${GRCH38_GTF}" > "${OUT_GTF}"
else
    ln -sf "$(readlink -f "${GRCH38}")" "${OUT_FASTA}"
    ln -sf "$(readlink -f "${GRCH38_GTF}")" "${OUT_GTF}"
fi

echo "[liftoff_reference] done (${SEX})"
