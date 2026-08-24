#!/bin/bash
# Download the reference files PRCGAP's in-workflow annotation steps need.
#
# Usage: bash download_reference.sh [OUTDIR]
#   OUTDIR   destination directory (default: reference)
#
# Fetches:
#   chm13v2.0_maskedY_rCRS.fa                       --chm13-fasta        (always required)
#   GRCh38.d1.vd1.fa                                --grch38-fasta       (--run-liftoff / --run-chain-files)
#   Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf  --grch38-gtf         (--run-liftoff)
#   chm13v2.0_censat_v2.1.bed                       --chm13-censat       (--run-chain-files)
#   centromeres.txt.gz                              --grch38-centromeres (--run-chain-files)
#   GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed  --grch38-exclusions (--run-chain-files)
#
# The heavier annotation PRCGAP still expects from assembly_workflow (sedef
# segdups, cenSat, RepeatMasker, misassembly) is not covered here.

set -o errexit
set -o nounset
set -o pipefail

OUTDIR="${1:-reference}"
mkdir -p "${OUTDIR}"

# ---- CHM13v2.0 genome FASTA (copynumber, INDEL liftover, chain files) ----
wget -O "${OUTDIR}/chm13v2.0_maskedY_rCRS.fa.gz" \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
gunzip -f "${OUTDIR}/chm13v2.0_maskedY_rCRS.fa.gz"

# ---- GRCh38 genome FASTA (liftoff source, chain files, INDEL liftover) ----
wget --content-disposition -O "${OUTDIR}/GRCh38.d1.vd1.fa.tar.gz" \
  'https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834'
tar xvzf "${OUTDIR}/GRCh38.d1.vd1.fa.tar.gz" -C "${OUTDIR}"

# ---- CHM13 cenSat annotation BED (chain-file masking; copy-number plot) ----
wget -P "${OUTDIR}" \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed

# ---- GRCh38 centromeres (chain-file masking) ----
wget -P "${OUTDIR}" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz

# ---- GRCh38 GRC exclusion regions BED (chain-file masking) ----
wget -P "${OUTDIR}" \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed

# ---- GRCh38 Ensembl 112 GTF for liftoff, renamed to chr* contigs ----
# Liftoff matches contig names between the GTF and the reference FASTA, and the
# GDC GRCh38 FASTA uses chr*; Ensembl ships bare 1/2/.../X/Y/MT.
wget -P "${OUTDIR}" \
  https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.chr.gtf.gz
zgrep -v "#" "${OUTDIR}/Homo_sapiens.GRCh38.112.chr.gtf.gz" \
  | sed 's/^\([0-9]\|X\|Y\|MT\)/chr\1/' \
  | sed 's/^chrMT/chrM/' \
  > "${OUTDIR}/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf"

cat <<EOF

References written to ${OUTDIR}

Pass them to setup_workflow.py, e.g.:

  python setup_workflow.py \\
    --samplesheet samples.tsv \\
    --chm13-fasta ${OUTDIR}/chm13v2.0_maskedY_rCRS.fa \\
    --grch38-fasta ${OUTDIR}/GRCh38.d1.vd1.fa \\
    --grch38-gtf ${OUTDIR}/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf \\
    --chm13-censat ${OUTDIR}/chm13v2.0_censat_v2.1.bed \\
    --grch38-centromeres ${OUTDIR}/centromeres.txt.gz \\
    --grch38-exclusions ${OUTDIR}/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed \\
    --run-dna-brnn --run-liftoff --run-chain-files
EOF
