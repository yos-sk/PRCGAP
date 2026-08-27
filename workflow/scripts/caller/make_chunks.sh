#!/bin/bash
# Group a reference's contigs into caller chunks. Contigs are never split: every
# chunk BED covers whole contigs, 0..length. Contigs >= MIN_BP get a chunk of
# their own (largest first, so the long jobs are queued first) and all remaining
# small/unplaced contigs are bundled into one final chunk, which keeps the
# per-job fixed cost (model load) off dozens of tiny contigs. MIN_BP=0 gives
# strictly one chunk per contig.
#
# The critical path of the scattered run is the largest single contig, so
# packing beyond "one big contig per chunk" cannot shorten it.
#
# Per chunk NNN emits:
#   NNN.bed   BED covering the whole contig(s)   -> DeepSomatic --regions
#   NNN.ctg   comma-separated contig list        -> ClairS -c
# plus manifest.tsv: chunk_id, n_contig, total_bp, contigs

set -xv
set -o errexit
set -o nounset
set -o pipefail

FAI=$1
CHUNK_DIR=$2
MIN_BP=${3:-1000000}

mkdir -p ${CHUNK_DIR}
rm -f ${CHUNK_DIR}/*.bed ${CHUNK_DIR}/*.ctg ${CHUNK_DIR}/manifest.tsv
touch ${CHUNK_DIR}/manifest.tsv

i=0
while IFS=$'\t' read -r ctg len; do
    i=$((i+1))
    id=$(printf "%03d" "${i}")
    printf "%s\t0\t%s\n" "${ctg}" "${len}" > "${CHUNK_DIR}/${id}.bed"
    printf "%s\n" "${ctg}" > "${CHUNK_DIR}/${id}.ctg"
    printf "%s\t1\t%s\t%s\n" "${id}" "${len}" "${ctg}" >> "${CHUNK_DIR}/manifest.tsv"
done < <(sort -k2,2nr "${FAI}" | awk -F'\t' -v m="${MIN_BP}" '$2>=m {print $1"\t"$2}')

n_large=$(sort -k2,2nr "${FAI}" | awk -F'\t' -v m="${MIN_BP}" '$2>=m' | wc -l)
id=$(printf "%03d" $((n_large+1)))
sort -k2,2nr "${FAI}" | awk -F'\t' -v m="${MIN_BP}" '$2<m {print $1"\t0\t"$2}' > "${CHUNK_DIR}/${id}.bed"
if [ -s "${CHUNK_DIR}/${id}.bed" ]; then
    cut -f1 "${CHUNK_DIR}/${id}.bed" | paste -sd, - > "${CHUNK_DIR}/${id}.ctg"
    awk -F'\t' -v id="${id}" '{n++; s+=$3; c=(c=="")?$1:c","$1} END{print id"\t"n"\t"s"\t"c}' \
        "${CHUNK_DIR}/${id}.bed" >> "${CHUNK_DIR}/manifest.tsv"
else
    rm -f "${CHUNK_DIR}/${id}.bed"
fi

sort -k1,1 -o "${CHUNK_DIR}/manifest.tsv" "${CHUNK_DIR}/manifest.tsv"

if [ ! -s "${CHUNK_DIR}/manifest.tsv" ]; then
    echo "caller_make_chunks.sh: no contigs found in ${FAI}" >&2
    exit 1
fi

echo "chunks: $(wc -l < "${CHUNK_DIR}/manifest.tsv")"
