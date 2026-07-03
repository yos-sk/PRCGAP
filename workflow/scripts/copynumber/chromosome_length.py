#!/usr/bin/env python3
# Emit per-chromosome lengths (chr<TAB>start<TAB>end, with header; chrM skipped)
# from a FASTA. Used by copynumber.sh at plot time to build the CHM13 chromosome
# lengths passed to plot_copy_number.R --chm13_lengths (so each chromosome can be
# extended to its CHM13 end / unaligned terminal margin). Ported from
# PRCGAP-paper/reference/scripts/chromosome_length.py (added the missing
# `import gzip` so .gz FASTAs also work).

import sys
import gzip

def main():
    input_fasta = sys.argv[1]
    print("chr\tstart\tend")

    if input_fasta[-3:] == ".gz":
        f = gzip.open(input_fasta, "rt")
    else:
        f = open(input_fasta, "r")

    chrom = ""
    seq_len = 0
    for line in f:
        if line[0] == ">":
            if chrom != "":
                if chrom != "chrM":
                    print(chrom, 1, seq_len, sep="\t")
            chrom = line.rstrip("\n").split()[0][1:]
            seq_len = 0
        else:
            seq_len += len(line.rstrip("\n"))

    if chrom != "chrM":
        print(chrom, 1, seq_len, sep="\t")
    f.close()

if __name__ == "__main__":
    main()
