#!/bin/bash
# Run ClairS on one contig chunk. The reference and the (possibly reheadered)
# BAMs come from caller_scatter_setup.sh, so no '#' handling happens here.

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR_BAM=$1
NORMAL_BAM=$2
REFERENCE_FA=$3
CTG_FILE=$4
OUTPUT_DIR=$5
THREAD=${6:-8}
# ClairS --platform (set per seqtype by the workflow; no default), e.g.
# hifi_sequel2, hifi_revio, ont_r9_guppy, ont_r10_guppy,
# ont_r10_dorado_sup_5khz_ssrs, ont_r10_dorado_sup_4khz, ...
MODEL=$7

mkdir -p ${OUTPUT_DIR}

/opt/bin/run_clairs \
    -T ${TUMOR_BAM} \
    -N ${NORMAL_BAM} \
    -R ${REFERENCE_FA} \
    -o ${OUTPUT_DIR} \
    -t ${THREAD} \
    -p ${MODEL} \
    -c "$(cat ${CTG_FILE})" \
    --include_all_ctgs \
    --enable_indel_calling \
    --remove_intermediate_dir

echo ${?}
