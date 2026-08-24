#!/bin/bash
# Full-length young LINE-1 (L1HS / L1PA2-L1PA5) for one haplotype.
#
# Two passes. The first finds candidate regions with L1.3, which is 95.5-99%
# identical to every young subfamily so one query is enough. The second lets
# Dfam's subunit models set the boundary and the subfamily: the 3'end model
# names the element -- that is how RepeatMasker itself decides -- and the 5'end
# model marks where it starts, which the whole-length query cannot do because
# blastn does not reach the 5'UTR of the older subfamilies.
#
# Runs inside the nanomonsv singularity container (blastn + bedtools + python3).
#
# Required positional args:
#   $1  SAMPLE       sample name (filename prefix)
#   $2  HAP          hap1 | hap2
#   $3  FASTA        that haplotype's assembly FASTA
#   $4  OUTPUT_DIR   destination dir
#   $5  THREADS      blastn -num_threads
#   $6  SCRIPT_DIR   absolute path to workflow/scripts
#   $7  RESOURCE_DIR absolute path to workflow/resources/line1

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP=$2
FASTA=$3
OUTPUT_DIR=$4
THREADS=${5:-8}
SCRIPT_DIR=$6
RESOURCE_DIR=$7

WORK_DIR="${OUTPUT_DIR}/workspace"
mkdir -p "${WORK_DIR}"

PREFIX="${WORK_DIR}/${SAMPLE}.${HAP}"
QUERY="${RESOURCE_DIR}/L1.3.fa"
QLEN=$(grep -v "^>" "${QUERY}" | tr -d "\n" | wc -c)

# makeblastdb builds an LMDB-backed v5 database through mmap, which is not
# reliable on Lustre, so keep both databases on node-local disk.
LOCAL_DIR=$(mktemp -d "${TMPDIR:-/tmp}/line1_${SAMPLE}_${HAP}_XXXXXX")
trap 'rm -rf "${LOCAL_DIR}"' EXIT

# bedtools writes a .fai next to the FASTA it is given, so give it a symlink in
# the workspace rather than the assembly itself, which may sit somewhere
# read-only.
ASM_LINK="${WORK_DIR}/${SAMPLE}.${HAP}.asm.fa"
ln -sf "$(readlink -f "${FASTA}")" "${ASM_LINK}"

# ---------- pass 1: candidate regions ----------
# -parse_seqids so blastdbcmd can report contig names, which saves needing a
# .fai just to clamp the padding to contig ends.
makeblastdb -in "${FASTA}" -dbtype nucl -parse_seqids -out "${LOCAL_DIR}/asm" > /dev/null

blastn -task dc-megablast -query "${QUERY}" -db "${LOCAL_DIR}/asm" \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
    -max_target_seqs 1000000 -evalue 1e-20 -num_threads "${THREADS}" \
    > "${PREFIX}.detect.tsv"

python3 "${SCRIPT_DIR}"/annotation/l1_chain.py \
    --input "${PREFIX}.detect.tsv" --query-len "${QLEN}" \
    --out "${PREFIX}.candidates.bed"

# Pad so the subunit models can reach past where the chain stopped, and name
# each region so the FASTA header identifies it. bedtools appends "(strand)",
# which l1_refine.py strips.
blastdbcmd -db "${LOCAL_DIR}/asm" -entry all -outfmt "%a	%l" > "${PREFIX}.contig_len.tsv"
awk -F'\t' -v OFS='\t' -v p=1000 'NR == FNR { len[$1] = $2; next }
    { s = $2 - p; if (s < 0) s = 0
      e = $3 + p; if (e > len[$1]) e = len[$1]
      print $1, s, e, $1 ":" s + 1 "-" e, 0, $4 }' \
    "${PREFIX}.contig_len.tsv" "${PREFIX}.candidates.bed" > "${PREFIX}.padded.bed"

# ---------- pass 2: boundary + subfamily ----------
if [ -s "${PREFIX}.padded.bed" ]; then
    # -s reverse-complements minus-strand regions, so inside a locus the element
    # always reads 5'->3' and the subunit models land at predictable ends.
    bedtools getfasta -fi "${ASM_LINK}" -bed "${PREFIX}.padded.bed" -s -nameOnly \
        > "${PREFIX}.loci.fa"

    makeblastdb -in "${PREFIX}.loci.fa" -dbtype nucl -out "${LOCAL_DIR}/loci" > /dev/null
    for MODELS in 3end 5end; do
        blastn -task blastn -word_size 11 -reward 1 -penalty -1 \
            -gapopen 2 -gapextend 1 -evalue 1e-10 -num_threads "${THREADS}" \
            -query "${RESOURCE_DIR}/l1_${MODELS}.fa" -db "${LOCAL_DIR}/loci" \
            -max_target_seqs 100000 \
            -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore" \
            > "${PREFIX}.${MODELS}.tsv"
    done
else
    : > "${PREFIX}.loci.fa"
    : > "${PREFIX}.3end.tsv"
    : > "${PREFIX}.5end.tsv"
fi

python3 "${SCRIPT_DIR}"/annotation/l1_refine.py \
    --loci "${PREFIX}.padded.bed" \
    --three "${PREFIX}.3end.tsv" --five "${PREFIX}.5end.tsv" \
    --out "${OUTPUT_DIR}/${SAMPLE}.${HAP}.LINE1.bed"

echo "[line1_hap] done: ${SAMPLE} ${HAP}"
