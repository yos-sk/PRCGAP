#!/bin/bash
# Rebuild the three FASTA files this directory ships.
#
# They are checked in because they are small and public, so enabling run_line1
# needs no download. This script exists so their provenance stays reproducible:
# run it and diff the output against what is committed.
#
#   L1.3.fa      GenBank L19088, fetched from NCBI (network required)
#   l1_3end.fa   *_3end LINE/L1 subunit models, human lineage
#   l1_5end.fa   *_5end   "
#
# The subunit models are taken from the RepeatMasker library inside the
# nanomonsv container rather than from Dfam directly: RepeatMasker ships the
# same Dfam entries already split into the subunit form the pipeline scores
# against, and its headers carry the repeat class, which Dfam's own FASTA does
# not. That library covers every species, so the clade in each header is used
# to keep the human lineage and drop the rodent-specific L1.
#
# Usage:
#   build_line1_resources.sh [-c <container>] [-o <output dir>] [-l <lib path>]
#
# Defaults to <repo>/images/nanomonsv.sif and this directory.

set -o errexit
set -o nounset
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The script lives in resource/scripts/ and writes the models to resource/line1/.
OUT_DIR="${HERE}/../line1"
mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"
CONTAINER="${REPO}/images/nanomonsv.sif"
RM_LIB=/opt/conda/share/RepeatMasker/Libraries/RepeatMasker.lib
SKIP_FETCH=0

# The clades on the path from the root to Homo sapiens. Everything else in the
# library is another lineage's L1 -- 51 rodent subunit models in Dfam 3.9.
HUMAN_CLADES="root|Eukaryota|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria_mammals|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo_sapiens"

while getopts "c:o:l:nh" opt; do
    case "${opt}" in
        c) CONTAINER=${OPTARG} ;;
        o) OUT_DIR=${OPTARG} ;;
        l) RM_LIB=${OPTARG} ;;
        n) SKIP_FETCH=1 ;;
        h) sed -n '2,22p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) exit 2 ;;
    esac
done

mkdir -p "${OUT_DIR}"

if [ "${SKIP_FETCH}" -eq 0 ]; then
    echo "### L1.3 (GenBank L19088) -> ${OUT_DIR}/L1.3.fa"
    curl -fsS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\
?db=nuccore&id=L19088&rettype=fasta&retmode=text" -o "${OUT_DIR}/L1.3.fa"
    printf "  %s bp\n" "$(grep -v '^>' "${OUT_DIR}/L1.3.fa" | tr -d '\n' | wc -c)"
fi

# Header form: >L1HS_3end#LINE/L1 @Homininae [S:45,55]
# Keep the family name only, so a model name in the blast output reads L1HS_3end.
for TAG in 3end 5end; do
    echo "### ${TAG} models -> ${OUT_DIR}/l1_${TAG}.fa"
    singularity exec "${CONTAINER}" \
        awk -v tag="_${TAG}#LINE/L1" -v clades="${HUMAN_CLADES}" '
            BEGIN { n = split(clades, a, "|"); for (i = 1; i <= n; i++) ok[a[i]] = 1 }
            /^>/ {
                if (keep) print seq
                keep = 0; seq = ""
                if (index($0, tag) > 0) {
                    clade = $0; sub(/.*@/, "", clade); sub(/[ \t].*/, "", clade)
                    if (clade in ok) {
                        keep = 1
                        split(substr($0, 2), f, "#"); print ">" f[1]
                    }
                }
                next
            }
            keep { seq = seq $0 }
            END { if (keep) print seq }' "${RM_LIB}" > "${OUT_DIR}/l1_${TAG}.fa"
    printf "  %s models\n" "$(grep -c '^>' "${OUT_DIR}/l1_${TAG}.fa")"
done

echo "[build_line1_resources] done"
