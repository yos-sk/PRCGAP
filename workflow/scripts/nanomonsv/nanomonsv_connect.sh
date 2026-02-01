#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

NANOMONSV_RESULT=$1
SUPPORT_READ_FILE=$2
OUTPUT_PREFIX=$3
SCRIPT_DIR=$4

python3 ${SCRIPT_DIR}/nanomonsv/connect.py \
    ${NANOMONSV_RESULT} \
    ${SUPPORT_READ_FILE} \
    ${OUTPUT_PREFIX} 
    
echo ${?}
