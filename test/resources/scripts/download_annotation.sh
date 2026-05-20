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

# GENCODE v4.6 with MANE v1.3
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gff3.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.summary.txt.gz
python3 proc_gencode_bed_mane_chr20.py gencode.v46.basic.annotation.gff3.gz MANE.GRCh38.v1.3.summary.txt.gz ../annotation/gencode.v46.basic.annotation.chr20.mane
gzip ../annotation/gencode.v46.basic.annotation.chr20.mane.transcript.bed
rm gencode.v46.basic.annotation.gff3.gz MANE.GRCh38.v1.3.summary.txt.gz