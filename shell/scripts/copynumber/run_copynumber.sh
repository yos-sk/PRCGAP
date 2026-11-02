#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu=4G
#SBATCH -p mjobs,rjobs

# Copy-number analysis of one tumor (copynumber image):
#   1. prepare_reference.sh : CHM13 for the reference tables (female: without chrY) + chromosome lengths
#   2. ref_table.sh         : CHM13 -> assembly alignment, reference tables per haplotype
#   3. copynumber.sh        : windowed tumor/normal depth, CBS, plot
#
# PLOIDY_HAP1 / PLOIDY_HAP2: "-" for automatic estimation (both must be given to override)
# CENSAT_BED: cenSat annotation of the assembly for the plot; "-" to use HAP1_SATELLITE + HAP2_SATELLITE

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
LABEL=$3
ASSEMBLY_HAP1=$4
ASSEMBLY_HAP2=$5
TUMOR_BAM=$6
NORMAL_BAM=$7
REFERENCE=$8
WORK_DIR=$9
OUTPUT_DIR=${10}
SEX=${11}
CHM13_CENSAT=${12}
CENSAT_BED=${13}
HAP1_SATELLITE=${14}
HAP2_SATELLITE=${15}
BIN_WIDTH=${16}
PLOIDY_HAP1=${17}
PLOIDY_HAP2=${18}
HAP1_LABEL=${19}
HAP2_LABEL=${20}
THREADS=${SLURM_CPUS_PER_TASK:-16}

IMAGE=./images/copynumber.sif

BIND=$(for f in ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} ${REFERENCE} ${CHM13_CENSAT} ${CENSAT_BED} ${HAP1_SATELLITE} ${HAP2_SATELLITE}; do
    if [ ${f} != "-" ]; then dirname ${f}; dirname $(readlink -f ${f}); fi
done | sort -u | paste -sd, -)

# The image sets WORKDIR /opt: keep the current directory explicitly.
run() { singularity exec --bind ${BIND} --pwd "$(pwd -P)" ${IMAGE} "$@"; }

mkdir -p ${WORK_DIR} ${OUTPUT_DIR}

run bash ./scripts/copynumber/prepare_reference.sh \
    ${REFERENCE} \
    ${WORK_DIR}/reference

run bash ./scripts/copynumber/ref_table.sh \
    ${NORMAL} \
    ${ASSEMBLY_HAP1} \
    ${ASSEMBLY_HAP2} \
    ${SEX} \
    ${WORK_DIR}/reference \
    ${REFERENCE} \
    ${CHM13_CENSAT} \
    ${WORK_DIR} \
    ${OUTPUT_DIR} \
    ${THREADS}

run bash ./scripts/copynumber/copynumber.sh \
    ${TUMOR} \
    ${NORMAL} \
    ${LABEL} \
    ${TUMOR_BAM} \
    ${NORMAL_BAM} \
    ${SEX} \
    ${WORK_DIR}/reference \
    ${CHM13_CENSAT} \
    ${CENSAT_BED} \
    ${HAP1_SATELLITE} \
    ${HAP2_SATELLITE} \
    ${BIN_WIDTH} \
    ${PLOIDY_HAP1} \
    ${PLOIDY_HAP2} \
    ${HAP1_LABEL} \
    ${HAP2_LABEL} \
    ${WORK_DIR} \
    ${OUTPUT_DIR} \
    ${THREADS}
