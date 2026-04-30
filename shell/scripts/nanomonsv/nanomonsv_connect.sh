#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

NANOMONSV_RESULT=$1
SUPPORT_READ_FILE=$2
OUTPUT_PREFIX=$3

python3 ./scripts/nanomonsv/connect.py \
    ${NANOMONSV_RESULT} \
    ${SUPPORT_READ_FILE} \
    ${OUTPUT_PREFIX} 
    
echo ${?}
