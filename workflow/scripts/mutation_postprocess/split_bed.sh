#!/bin/bash
# Split BED file for parallel pileup processing.
# Used by Snakemake for mutation postprocessing.
#
# Phase A change (2026-07): instead of `split -n l/16` (which scatters positions
# across the whole genome so every chunk must scan the BAM end-to-end), we now
# split BY CONTIG and pack contigs into balanced chunks (Longest-Processing-Time
# bin packing on per-contig variant counts). Each chunk therefore covers a small
# set of contigs, letting pileup.sh use `samtools mpileup -r <contig>` (index
# seek + early termination) instead of a full-BAM scan.

set -xv
set -o errexit
set -o nounset
set -o pipefail

PARSED_BED=$1
OUTPUT_DIR=$2
PILEUP_TASKS=$3
NUM_CHUNKS=${4:-16}   # target maximum number of chunks (upper bound)

WORK="${OUTPUT_DIR}/pileup/workspace"
mkdir -p "${WORK}"
# clean any stale chunk files from a previous run
find "${WORK}" -maxdepth 1 -name 'input_*' -delete 2>/dev/null || true

# 1) per-contig variant counts -> LPT-balance contigs into <= NUM_CHUNKS bins.
#    Emit "contig<TAB>bin" (bin = 1..NUM_CHUNKS). Each contig stays whole.
awk -v n="${NUM_CHUNKS}" '
    { cnt[$1]++; if(!($1 in seen)){seen[$1]=1; order[++m]=$1} }
    END{
        for(i=1;i<=m;i++) c[i]=order[i]
        # sort contigs by descending variant count (small m; simple selection)
        for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) if(cnt[c[j]]>cnt[c[i]]){t=c[i];c[i]=c[j];c[j]=t}
        if(n>m) n=m
        for(b=1;b<=n;b++) load[b]=0
        for(i=1;i<=m;i++){
            mb=1; for(b=2;b<=n;b++) if(load[b]<load[mb]) mb=b
            bin[c[i]]=mb; load[mb]+=cnt[c[i]]
        }
        for(i=1;i<=m;i++) print c[i]"\t"bin[c[i]]
    }' "${PARSED_BED}" > "${WORK}/contig_bin.map"

# 2) route each variant row to its bin sub-bed (sorted by contig,pos within bin)
sort -k1,1 -k2,2n "${PARSED_BED}" > "${WORK}/parsed_sorted.bed"
awk -v w="${WORK}" '
    NR==FNR{ bin[$1]=$2; next }
    { print > (w"/input_" sprintf("%02d", bin[$1])) }' \
    "${WORK}/contig_bin.map" "${WORK}/parsed_sorted.bed"

# 3) build the task list from the chunks that actually got rows
: > "${PILEUP_TASKS}"
for f in "${WORK}"/input_*; do
    [ -e "${f}" ] || continue
    case "${f}" in *.pileup) continue;; esac
    echo -e "${f}\t${f}.pileup" >> "${PILEUP_TASKS}"
done

echo ${?}
