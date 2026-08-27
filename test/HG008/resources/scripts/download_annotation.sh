#!/bin/bash

# CMRG
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-021-01158-1/MediaObjects/41587_2021_1158_MOESM4_ESM.tsv
python3 extract_cmrg_gene.py > ../annotation/cmrg_genes.list
rm 41587_2021_1158_MOESM4_ESM.tsv

# Cancer Gene Census - need to register and login to download the files. Please download manually and place in test/resources/annotation/cancer_gene_census.tsv.

# GnomAD SVs
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz
zcat gnomad.v4.1.sv.sites.bed.gz | sort -k 1,1 -k 2,2n | bgzip -f -c > ../annotation/gnomad.v4.1.sv.sites.bed.gz
tabix -p bed ../annotation/gnomad.v4.1.sv.sites.bed.gz

# GnomAD SNVs and indels
wget -q https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr20.vcf.bgz -O ../annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz
wget -q https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr20.vcf.bgz.tbi -O ../annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz.tbi

# The MANE summary is a reference, not a test fixture: it is downloaded once by
# ../../../../resource/scripts/download_reference.sh and read from there, so the
# GENCODE-derived transcript BED this script used to build is no longer needed.