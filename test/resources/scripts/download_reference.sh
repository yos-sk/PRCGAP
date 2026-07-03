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

echo ${?}