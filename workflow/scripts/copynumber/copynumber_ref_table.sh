#!/bin/bash
# Contig ↔ CHM13 correspondence table for one haplotype: align the reference to
# the assembly and reduce the primary alignments to a table.
#
# The assembly is not masked. Masking with the dna-brnn satellite BED changed no
# chromosome or strand assignment and removed 32.8 Mb of span, all acrocentric
# (workspace/reference_table/README.md §2).
#
# Required positional args:
#   $1  HAP        hap1 | hap2 (labels the scratch files only)
#   $2  ASSEMBLY   that haplotype's assembly FASTA
#   $3  REFERENCE  staged CHM13 reference (copynumber_reference.sh)
#   $4  CENSAT     CHM13 cenSat BED, or empty
#   $5  WORK_DIR   scratch dir
#   $6  OUT_TABLE  reference table to write
#   $7  SCRIPT_DIR absolute path to workflow/scripts
#   $8  THREADS    minimap2 -t

set -xv
set -o errexit
set -o nounset
set -o pipefail

HAP=$1
ASSEMBLY=$2
REFERENCE=$3
CENSAT=${4:-}
WORK_DIR=$5
OUT_TABLE=$6
SCRIPT_DIR=$7
THREADS=${8:-8}

mkdir -p "${WORK_DIR}" "$(dirname "${OUT_TABLE}")"

# Only the primary contig↔chromosome mapping is wanted. --secondary=no rather
# than dropping tp:A:S afterwards: that grep leaves the secondary inversions
# (tp:A:i, lowercase) behind, and only --min-mapq keeps them out of the table.
minimap2 -cx asm5 --secondary=no -t "${THREADS}" \
    "${ASSEMBLY}" "${REFERENCE}" \
    > "${WORK_DIR}/${HAP}_ref.paf"

# Telomeric repeats at a contig end mark a chromosome terminus, which is how the
# acrocentric short arms get their span back.
TELO_ARGS=()
if command -v seqtk > /dev/null 2>&1; then
    seqtk telo "${ASSEMBLY}" > "${WORK_DIR}/${HAP}.telo.bed"
    TELO_ARGS=(--telo "${WORK_DIR}/${HAP}.telo.bed")
else
    echo "[copynumber_ref_table] seqtk not found; skipping the telomere step" >&2
fi

CENSAT_ARGS=()
if [ -n "${CENSAT}" ]; then
    CENSAT_ARGS=(--censat "${CENSAT}")
fi

python3 "${SCRIPT_DIR}"/copynumber/make_reference_table.py \
    -i "${WORK_DIR}/${HAP}_ref.paf" \
    "${CENSAT_ARGS[@]}" \
    "${TELO_ARGS[@]}" \
    > "${OUT_TABLE}"

echo "[copynumber_ref_table] done: ${HAP}"
