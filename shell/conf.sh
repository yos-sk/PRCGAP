#!/bin/bash

OUTPUT=
TUMOR=
CONTROL=

ASSEMBLER=

SEX=

OUTPUT_PREFIX=

### De novo assembly ###
ASSEMBLY_HAP1=
ASSEMBLY_HAP2=


### Input sequencing data ###
TUMOR_HIFI=
CONTROL_HIFI=

TUMOR_ONT=
CONTROL_ONT=

TUMOR_METHYLATION_HIFI=
CONTROL_METHYLATION_HIFI=


### Rerefence genome and annotation files ###
REFERENCE=

SIMPLE_REPEAT=
HAP1_SATELLITE=
HAP2_SATELLITE=
GTF_FILE=
LINE1_BED=

### Copy number ###
# CHM13 cenSat v2.1 BED(.gz) (reference tables and plot)
CHM13_CENSAT=
# cenSat annotation of the assembly for the plot; empty = HAP1_SATELLITE + HAP2_SATELLITE
CENSAT_BED=
CN_BIN_WIDTH=0.05
# "-" = estimate the tumor ploidy (set both to override)
CN_PLOIDY_HAP1=-
CN_PLOIDY_HAP2=-
CN_HAP1_LABEL=Haplotype1
CN_HAP2_LABEL=Haplotype2

### Containers ###
# Comma-separated paths to bind into the containers (data outside $HOME and the current directory)
BIND_PATHS=

# TRUE or FALSE
BAM_REFINER="TRUE"
METHYLATION="TRUE"
COPYNUMBER="TRUE"
NANOMONSV_PARSE="TRUE"
NANOMONSV_GET="TRUE"
NANOMONSV_POST="TRUE"
NANOMONSV_INSERT_CLASSIFY="TRUE"
NANOMONSV_CONNECT="FALSE"
NANOMONSV_MERGE="FALSE"
CLAIRS="FALSE"
DEEPSOMATIC="TRUE"
CLAIRS_POST="FALSE"
DEEPSOMATIC_POST="TRUE"
