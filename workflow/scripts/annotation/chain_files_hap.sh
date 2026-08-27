#!/bin/bash
# Build the assembly → reference chain for one (haplotype, reference) pair.
# Runs inside the chain_files singularity container (minimap2 + transanno +
# chaintools + rustybam + paf2chain).
#
# Pipeline, following the assembly_workflow implementation: minimap2 asm5 with
# the assembly as target and the masked reference as query, converted to a
# chain, split into single-block chains, re-expressed as PAF, broken/trimmed so
# no block exceeds 10 kb or overlaps its neighbours, converted back to a chain,
# and finally inverted so the chain reads assembly → reference.
#
# Required positional args:
#   $1  SAMPLE           sample name (filename prefix)
#   $2  HAP              hap1 | hap2
#   $3  FASTA            that haplotype's assembly FASTA
#   $4  REFERENCE        GRCh38 | chm13 (label only)
#   $5  REFERENCE_FASTA  masked reference FASTA (prepare_mask_regions.sh)
#   $6  OUT_CHAIN        inverted chain to write
#   $7  THREADS          minimap2 -t

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
HAP=$2
FASTA=$3
REFERENCE=$4
REFERENCE_FASTA=$5
OUT_CHAIN=$6
THREADS=${7:-8}

CHAINTOOLS=/tools/chaintools/chaintools

WORK_DIR="$(dirname "${OUT_CHAIN}")"
mkdir -p "${WORK_DIR}"
PREFIX=${WORK_DIR}/${SAMPLE}_${HAP}.${REFERENCE}

minimap2 -cx asm5 -t "${THREADS}" "${FASTA}" "${REFERENCE_FASTA}" \
    > "${PREFIX}.paf"
transanno minimap2chain "${PREFIX}.paf" --output "${PREFIX}.chain"

python3 "${CHAINTOOLS}/split.py" \
    -c "${PREFIX}.chain" \
    -o "${PREFIX}-split.chain"
python3 "${CHAINTOOLS}/to_paf.py" \
    -c "${PREFIX}-split.chain" \
    -t "${REFERENCE_FASTA}" \
    -q "${FASTA}" \
    -o "${PREFIX}-split.paf"

rb break-paf --max-size 10000 < "${PREFIX}-split.paf" \
    | rb trim-paf -r | rb invert | rb trim-paf -r | rb invert \
    > "${PREFIX}-out.paf"
paf2chain -i "${PREFIX}-out.paf" > "${PREFIX}-out.chain"
python3 "${CHAINTOOLS}/invert.py" \
    -c "${PREFIX}-out.chain" \
    -o "${OUT_CHAIN}"

echo "[chain_files_hap] done: ${SAMPLE} ${HAP} → ${REFERENCE}"
