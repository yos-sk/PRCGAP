#!/bin/bash
#
# Pull all singularity images required by the PRCGAP workflow into ./images/
# (or the directory passed as $1).
#
# After this script finishes, setup_workflow.py picks up the images by their
# default paths (images/<tool>.sif) so users do not need to pass --*-image
# flags on the command line.
#
# Image URLs below are placeholders — replace with actual registry locations
# (docker://quay.io/..., docker://ghcr.io/..., etc.) once images are published.
# Until then, build from Dockerfile/<tool>/Dockerfile via singularity build.

set -euo pipefail

# tool name -> registry URL (replace once images are hosted)
declare -A IMAGES=(
    [bam_refiner]="docker://yosakam2/bam_refiner:0.3.6"
    [methylation]="docker://yosakam2/methylation:0.1.0"
    [copynumber]="docker://yosakam2/copynumber:0.1.0"
    [nanomonsv]="docker://friend1ws/nanomonsv:v0.8.0"
    [nanomonsv_postprocess]="docker://yosakam2/nanomonsv_postprocess:0.2.5"
    [clairs]="docker://yosakam2/clairs:0.4.0"
    [deepsomatic]="docker://yosakam2/deepsomatic:1.8.0"
    [point_mutation_postprocess]="docker://yosakam2/mutation_postprocess:0.1.2"
    [annotation]="docker://yosakam2/annotation:0.1.0"
)

for tool in "${!IMAGES[@]}"; do
    sif="./${tool}.sif"
    url="${IMAGES[$tool]}"
    if [ -f "${sif}" ]; then
        echo "[skip] ${sif} already exists"
        continue
    fi
    echo "[pull] ${tool} -> ${sif}"
    singularity pull "${sif}" "${url}"
done

echo "Done. Images placed in ./"