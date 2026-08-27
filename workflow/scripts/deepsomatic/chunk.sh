#!/bin/bash
# Run DeepSomatic on one contig chunk. The reference and BAMs come from
# caller_scatter_setup.sh.

set -xv
set -o errexit
set -o nounset
set -o pipefail

# Cap glibc malloc arenas. On many-core hosts (e.g. 192 cores) the default
# (cores * 8) arenas reserve hundreds of GB of virtual memory upfront, which
# trips SGE's s_vmem limit even though resident memory is small.
export MALLOC_ARENA_MAX=4

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
REFERENCE_FA=$5
REGIONS_BED=$6
OUTPUT_DIR=$7
THREAD=${8:-8}
MODEL_TYPE=${9:-PACBIO}
# postprocess_variants' own worker count. Its default is the node's online core
# count and it ignores both --num_shards and the scheduler allocation, so on a
# many-core host it forks one worker per host core onto our allocated cores.
POSTPROCESS_CPUS=${10:-1}

mkdir -p ${OUTPUT_DIR}
INTERMEDIATE=${TMPDIR:-${OUTPUT_DIR}}/intermediate_results_dir_$$

run_deepsomatic \
    --model_type=${MODEL_TYPE} \
    --ref=${REFERENCE_FA} \
    --reads_normal=${NORMAL_BAM} \
    --reads_tumor=${TUMOR_BAM} \
    --output_vcf=${OUTPUT_DIR}/output.vcf.gz \
    --sample_name_tumor=${TUMOR} \
    --sample_name_normal=${NORMAL} \
    --num_shards=${THREAD} \
    --postprocess_variants_extra_args="cpus=${POSTPROCESS_CPUS}" \
    --regions=${REGIONS_BED} \
    --logging_dir=${OUTPUT_DIR}/logs \
    --intermediate_results_dir=${INTERMEDIATE}

# The examples are tens of GB per chunk and are not needed downstream.
rm -rf ${INTERMEDIATE}

echo ${?}
