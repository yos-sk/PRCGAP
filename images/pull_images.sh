#!/bin/bash
# Pull every singularity image the PRCGAP workflow needs into this directory.
#
# File names follow `<image-key>.sif`, where <image-key> matches the keys used
# by setup_workflow.py / the workflow rules. With every image staged under the
# same dir using these names, `setup_workflow.py --images-dir <dir>` is enough
# — no per-image override flags needed.
#
# Usage (run from this images/ directory):
#   bash images/pull_images.sh                       # pull missing images
#   bash images/pull_images.sh --force               # re-pull everything (even if .sif exists)
#   bash images/pull_images.sh clairs deepsomatic    # pull just these modules
#   bash images/pull_images.sh --force clairs        # re-pull a single module
#
# Notes:
# - Requires `singularity` (or `apptainer`) on PATH.
# - Image URLs below are placeholders for unpublished images — replace once they
#   are hosted. Until then, build from Dockerfile/<tool>/Dockerfile.

set -euo pipefail

# (module-key, docker URI). Each entry produces <key>.sif in this directory.
# Stored as a plain array of "key uri" strings (not declare -A) so the script
# also runs under the bash 3.2 shipped on macOS.
declare -a IMAGES=(
    "bam_refiner                 docker://yosakam2/bam_refiner:v0.4.0"
    "methylation                 docker://yosakam2/methylation:v0.1.0"
    "copynumber                  docker://yosakam2/copynumber:v0.3.0"
    "nanomonsv                   docker://yosakam2/nanomonsv:v0.8.0"
    "nanomonsv_postprocess       docker://yosakam2/nanomonsv_postprocess:v0.2.6"
    "clairs                      docker://yosakam2/clairs:v0.4.0"
    "deepsomatic                 docker://yosakam2/deepsomatic:v1.8.0"
    "point_mutation_postprocess  docker://yosakam2/mutation_postprocess:v0.1.3"
    "annotation                  docker://yosakam2/annotation:v0.2"
    # In-workflow annotation generation (setup_workflow.py --run-dna-brnn /
    # --run-liftoff / --run-chain-files). Shared with the assembly_workflow repo.
    "dna_nn                      docker://yosakam2/dna-nn:v0.1"
    "liftoff                     docker://yosakam2/liftoff:1.6.3"
    "chain_files                 docker://yosakam2/chaintools:2a3b47e"
)

FORCE=0
SELECT=()
for arg in "$@"; do
    case "$arg" in
        --force|-f) FORCE=1 ;;
        -h|--help)
            sed -n '1,/^set -euo pipefail/p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        -*) echo "Unknown option: $arg" >&2; exit 1 ;;
        *) SELECT+=("$arg") ;;
    esac
done

# Pick singularity or apptainer (after arg parsing so --help works on hosts
# without either installed).
if command -v singularity >/dev/null 2>&1; then
    SING=singularity
elif command -v apptainer >/dev/null 2>&1; then
    SING=apptainer
else
    echo "Error: neither singularity nor apptainer is on PATH" >&2
    exit 1
fi

# Images are written to the current directory; run this script from images/.
# Resolve from the script's own location so the images land beside it
# whatever the working directory is.
OUT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

pull_one() {
    local key="$1"
    local uri="$2"
    local out="$OUT_DIR/${key}.sif"
    if [[ -f "$out" && $FORCE -eq 0 ]]; then
        echo "[skip ] $key  (already at $out; use --force to re-pull)"
        return 0
    fi
    echo "[pull ] $key  ←  $uri"
    # --force handles the case of a leftover partial download; singularity pull
    # otherwise refuses to overwrite.
    "$SING" pull --force "$out" "$uri"
}

for line in "${IMAGES[@]}"; do
    # Split on whitespace; first token = key, remainder = uri.
    read -r key uri <<<"$line"
    if [[ ${#SELECT[@]} -gt 0 ]]; then
        # When specific modules were requested, skip everything else.
        match=0
        for s in "${SELECT[@]}"; do
            [[ "$key" == "$s" ]] && match=1 && break
        done
        [[ $match -eq 0 ]] && continue
    fi
    pull_one "$key" "$uri"
done

echo
echo "Done. Use this directory with: setup_workflow.py --images-dir $(pwd)"
