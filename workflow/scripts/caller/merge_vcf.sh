#!/bin/bash
# Merge per-chunk caller VCFs into one sorted, bgzipped, tabix-indexed VCF.
#
# The chunks cover disjoint contig sets, so merging is a concatenation plus a
# sort into reference order. Two details make it worth a script:
#
#   - ClairS writes only the contigs it called on into each chunk header (1 line
#     for a single-contig chunk), so the merged ##contig block is regenerated
#     from the .fai instead of copied from a chunk. That also guarantees the
#     header matches the reference exactly. DeepSomatic already writes all
#     contigs, and regenerating is a no-op for it.
#   - Records are ordered by the .fai contig order, not by chunk order.
#
# '_HASH_' is restored to '#' when the reference was sanitised by
# caller_scatter_setup.sh.

set -xv
set -o errexit
set -o nounset
set -o pipefail

REFERENCE_FAI=$1
SANITIZED_FLAG=$2
OUTPUT_VCF=$3
SORT_MEM_MB=$4
shift 4
INPUTS=("$@")

if [ ${#INPUTS[@]} -eq 0 ]; then
    echo "merge_caller_vcf.sh: no input VCFs given" >&2
    exit 1
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/merge_caller_vcf.XXXXXX")
trap 'rm -rf "${WORK}"' EXIT

# Print just the header of a VCF. The awk stops at the end of the header, which
# makes zcat die of SIGPIPE, so pipefail has to be off for the duration or
# errexit would abort on a successful read.
vcf_header() {
    set +o pipefail
    zcat "$1" | awk '/^##/{print; next} /^#CHROM/{print; exit}'
    set -o pipefail
}

# ---- header ----------------------------------------------------------------
# ##fileformat must stay on line 1; collect the remaining non-contig meta lines
# from every chunk (first-seen order) so a FILTER/INFO only emitted by one chunk
# is not lost.
vcf_header "${INPUTS[0]}" | awk '/^##fileformat/{print; exit}' > "${WORK}/header"

for f in "${INPUTS[@]}"; do
    vcf_header "${f}" | awk '/^##/ && !/^##fileformat/ && !/^##contig=/'
done | awk '!seen[$0]++' >> "${WORK}/header"

awk -F'\t' '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${REFERENCE_FAI}" >> "${WORK}/header"

# #CHROM must agree across chunks; a mismatch means differing samples.
vcf_header "${INPUTS[0]}" | grep '^#CHROM' > "${WORK}/chrom"
for f in "${INPUTS[@]}"; do
    vcf_header "${f}" | grep '^#CHROM' > "${WORK}/chrom.this"
    if ! cmp -s "${WORK}/chrom" "${WORK}/chrom.this"; then
        echo "merge_caller_vcf.sh: #CHROM line differs in ${f}" >&2
        exit 1
    fi
done
cat "${WORK}/chrom" >> "${WORK}/header"

# ---- records ---------------------------------------------------------------
for f in "${INPUTS[@]}"; do
    zcat "${f}" | awk 'p{print} /^#CHROM/{p=1}'
done > "${WORK}/records"

# Decorate with the .fai contig ordinal, sort by (contig order, POS), undecorate.
awk -F'\t' -v OFS='\t' '
    NR==FNR { ord[$1]=FNR; next }
    {
        if (!($1 in ord)) { print "merge_caller_vcf.sh: contig " $1 " not in fai" > "/dev/stderr"; exit 1 }
        print ord[$1], $2, $0
    }' "${REFERENCE_FAI}" "${WORK}/records" \
    | sort -S "${SORT_MEM_MB}M" -k1,1n -k2,2n \
    | cut -f3- > "${WORK}/records.sorted"

# ---- emit ------------------------------------------------------------------
if [ "$(cat "${SANITIZED_FLAG}")" = "1" ]; then
    cat "${WORK}/header" "${WORK}/records.sorted" | sed 's/_HASH_/#/g' | bgzip -c > "${OUTPUT_VCF}"
else
    cat "${WORK}/header" "${WORK}/records.sorted" | bgzip -c > "${OUTPUT_VCF}"
fi

tabix -f -p vcf "${OUTPUT_VCF}"

echo ${?}
