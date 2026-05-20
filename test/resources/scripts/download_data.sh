#!/bin/bash

set -ex 

mkdir -p ../reads/ont
mkdir -p ../reads/hifi

BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab"

# HG008T ONT
samtools view -Shb -F 2308 \
    "$BASE/Northeastern-ONT-UL-20241216/HG008-T_CHM13v2.0_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.bam" \
    chr20 \
> ../reads/ont/HG008T.chr20.ont.bam

# HG008N ONT
samtools view -Shb -F 2308 \
    "$BASE/Northeastern_ONT-std_20240422/HG008-N-P_CHM13v2.0_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x.bam" \
    chr20 \
> ../reads/ont/HG008N.chr20.ont.bam

# HG008T HiFi
samtools view -Shb -F 2308 \
    "$BASE/PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x_CHM13v2.0.bam" \
    chr20 \
> ../reads/hifi/HG008T.chr20.hifi.bam

# HG008N HiFi
samtools view -Shb -F 2308 \
    "$BASE/PacBio_Revio_20240125/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_CHM13v2.0.bam" \
    chr20 \
> ../reads/hifi/HG008N.chr20.hifi.bam

mkdir -p ../asm
BASE_ASM="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Verkko_assemblies_05162024/HG008_N_asm_hifiherrohic_verkko2.2_20250218"
wget "$BASE_ASM/HG008N_verkko-assembly_20250218.haplotype1.fasta.gz" -O ../asm/HG008N.hap1.fa.gz && gunzip ../asm/HG008N.hap1.fa.gz
wget "$BASE_ASM/HG008N_verkko-assembly_20250218.haplotype2.fasta.gz" -O ../asm/HG008N.hap2.fa.gz && gunzip ../asm/HG008N.hap2.fa.gz

samtools faidx ../asm/HG008N.hap1.fa
samtools faidx ../asm/HG008N.hap2.fa
samtools faidx ../asm/HG008N.hap1.fa haplotype1-0000020 > ../asm/HG008N.hap1.chr20.fa
samtools faidx ../asm/HG008N.hap2.fa haplotype2-0000110 > ../asm/HG008N.hap2.chr20.fa

echo ${?}