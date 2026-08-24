#!/bin/bash
# Short-period tandem repeat regions for one haplotype, with ULTRA.
#
# nanomonsv uses this to drop indel-like SVs whose breakpoints both sit in the
# same tandem repeat, so what matters is the repeat unit itself, not long
# satellite arrays: period is capped at 10, matching the MaxPeriod RepeatMasker
# gives TRF in its Simple_repeat stage. The copy-number cut then keeps only
# arrays with enough repeats for an indel to be attributable to the unit.
#
# Runs inside the annotation singularity container (ultra + bgzip/tabix).
#
# Required positional args:
#   $1  SAMPLE       sample name (filename prefix)
#   $2  HAP          hap1 | hap2
#   $3  FASTA        that haplotype's assembly FASTA
#   $4  OUTPUT_DIR   destination dir
#   $5  THREADS      ultra -t
#   $6  MIN_COPIES   minimum #copies to keep (default 4)

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP=$2
FASTA=$3
OUTPUT_DIR=$4
THREADS=${5:-8}
MIN_COPIES=${6:-4}

WORK_DIR="${OUTPUT_DIR}/workspace"
mkdir -p "${WORK_DIR}"

TSV="${WORK_DIR}/${SAMPLE}.${HAP}.ultra_p10.tsv"
BED="${OUTPUT_DIR}/${SAMPLE}.${HAP}.simple_repeats.bed"

# --read_all: annotate every contig, not a sample of them.
# -i 3 -d 3: insertion/deletion penalties as used in the benchmark.
# --tsv -c: the TSV carries Score and #copies, which the BED does not.
ultra --read_all -p 10 -t "${THREADS}" -i 3 -d 3 --tsv -c \
    -o "${TSV}" "${FASTA}"

# TSV columns: SeqID Start End Period Score Consensus #copies ...
# Merge overlapping/abutting intervals in one pass; bedtools is not in this
# container and a merge over sorted intervals does not need it.
export LC_ALL=C
awk -F'\t' -v OFS='\t' -v m="${MIN_COPIES}" 'FNR > 1 && $7 >= m {print $1, $2, $3}' "${TSV}" \
    | sort -k1,1 -k2,2n \
    | awk -F'\t' -v OFS='\t' '
        NR == 1 { c = $1; s = $2; e = $3; next }
        $1 == c && $2 <= e { if ($3 > e) e = $3; next }
        { print c, s, e; c = $1; s = $2; e = $3 }
        END { if (NR) print c, s, e }' > "${BED}"

echo "[simple_repeat_hap] done: ${SAMPLE} ${HAP} ($(wc -l < "${BED}") regions)"
