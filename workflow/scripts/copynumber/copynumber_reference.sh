#!/bin/bash
# Stage the CHM13 reference the copy-number steps align against, plus the
# chromosome-length table the plot's terminal margin needs. Shared by every
# tumor and haplotype, so the chrY filter runs once per workflow.
#
# Required positional args:
#   $1  SEX          female | male ("female" drops chrY)
#   $2  REFERENCE    CHM13 reference FASTA
#   $3  OUT_FASTA    staged reference path
#   $4  OUT_LENGTHS  chromosome-length table path
#   $5  SCRIPT_DIR   absolute path to workflow/scripts

set -xv
set -o errexit
set -o nounset
set -o pipefail

SEX=$1
REFERENCE=$2
OUT_FASTA=$3
OUT_LENGTHS=$4
SCRIPT_DIR=$5

mkdir -p "$(dirname "${OUT_FASTA}")"

if [ "${SEX}" = "female" ]; then
    awk '/^>/ {p = ($0 !~ /^>chrY/)} p' "${REFERENCE}" > "${OUT_FASTA}"
else
    ln -sf "$(readlink -f "${REFERENCE}")" "${OUT_FASTA}"
fi

# chr/start/end with a header; chrM is skipped.
python3 "${SCRIPT_DIR}"/copynumber/chromosome_length.py "${OUT_FASTA}" > "${OUT_LENGTHS}"

echo "[copynumber_reference] done (${SEX})"
