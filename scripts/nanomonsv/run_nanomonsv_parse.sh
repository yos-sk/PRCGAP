#!/bin/bash
#SBATCH --mem-per-cpu=30G
#SBATCH -c 8
#SBATCH -p rjobs

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
INPUT_BAM=$2
OUTPUT_DIR=$3

singularity exec ./images/nanomonsv-devel.sif \
    bash ./scripts/nanomonsv/nanomonsv_parse.sh \
    ${SAMPLE} \
    ${INPUT_BAM} \
    ${OUTPUT_DIR}
    
