#!/bin/bash

mkdir -p ../reference
cd ../reference
# T2T-CHM13
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
gunzip chm13v2.0_maskedY_rCRS.fa.gz
samtools faidx chm13v2.0_maskedY_rCRS.fa

# CHM13 cenSat v2.1 (copy-number plot: fill assembly gaps with reference satellite).
# (CHM13 chromosome lengths are derived at plot time from the fasta via
# chromosome_length.py, so no separate length file is staged here.)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed
bgzip -f chm13v2.0_censat_v2.1.bed
tabix -p bed chm13v2.0_censat_v2.1.bed.gz

# GRCh38
wget --content-disposition 'https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834' -O GRCh38.d1.vd1.fa.tar.gz
tar xvzf GRCh38.d1.vd1.fa.tar.gz
samtools faidx GRCh38.d1.vd1.fa
rm GRCh38.d1.vd1.fa.tar.gz

# The three inputs PRCGAP's own annotation steps need on top of the two FASTAs
# and the cenSat BED above. test_configure.sh only passes them when the
# RUN_ANNOTATION=1 variant is configured; the default variant keeps using the
# pre-extracted chr20 annotation under ../annotation/.

# GRCh38 Ensembl 112 GTF for liftoff, renamed to chr* contigs to match the
# GDC GRCh38 FASTA (Ensembl ships bare 1/2/.../X/Y/MT).
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.chr.gtf.gz
zgrep -v "#" Homo_sapiens.GRCh38.112.chr.gtf.gz \
    | sed 's/^\([0-9]\|X\|Y\|MT\)/chr\1/' \
    | sed 's/^chrMT/chrM/' \
    > Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf

# GRCh38 centromeres + GRC exclusion regions, masked before chain-file alignment.
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed

echo ${?}