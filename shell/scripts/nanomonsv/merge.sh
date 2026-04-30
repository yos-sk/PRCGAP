#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

NANOMONSV_RESULT_1=$1
NANOMONSV_RESULT_2=$2
OUTPUT_FILE=$3

nanomonsv_postprocess merge \
    -i ${NANOMONSV_RESULT_1} \
    -j ${NANOMONSV_RESULT_2} \
    -o ${OUTPUT_FILE} \
    -s 90

echo ${?}
