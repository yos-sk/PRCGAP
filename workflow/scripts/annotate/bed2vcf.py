#!/usr/bin/env python3
"""Convert an INDEL-filtered `*.haplotyped.bed` row into a VCF record by
joining against the original (per-tool) somatic VCF on (contig, pos, ref,
alt). The output VCF inherits the original VCF header but prepends five
INFO descriptions so that the per-record INFO column can carry
`haplotype/contig/pos/ref/alt` through downstream `transanno liftvcf`.

This is functionally the reference `indel/bed2vcf.py` (PRCGAP-paper
analysis/point_mutation_analysis/scripts/annotation/indel/bed2vcf.py).
"""

import argparse
import gzip
import sys

import pysam


ADDED_INFO_HEADER = (
    '##INFO=<ID=haplotype,Number=1,Type=String,Description="Haplotype of variant allele.">\n'
    '##INFO=<ID=contig,Number=1,Type=String,Description="contig.">\n'
    '##INFO=<ID=pos,Number=1,Type=String,Description="position.">\n'
    '##INFO=<ID=ref,Number=1,Type=String,Description="reference allele.">\n'
    '##INFO=<ID=alt,Number=1,Type=String,Description="alternative allele.">'
)


def convert_main(args):
    with open(args.output_vcf_file, "w") as w:
        prev_header = ""
        with gzip.open(args.input_vcf_file, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    break
                items = line.rstrip("\n").split("=")
                if prev_header != items[0]:
                    if prev_header == "##INFO":
                        print(ADDED_INFO_HEADER, file=w)
                prev_header = items[0]
                print(line.rstrip("\n"), file=w)

        with open(args.input_bed_file, "r") as f, pysam.TabixFile(args.input_vcf_file) as tbx:
            for line in f:
                items = line.rstrip("\n").split("\t")
                contig = items[0]
                start = int(items[1])
                end = int(items[2])
                ref = items[3]
                alt = items[4]
                haplotype = items[12]

                try:
                    records = tbx.fetch(contig, start, end)
                except Exception:
                    print(
                        f"Could not fetch the region: {contig}:{start}-{end}, "
                        f"Ref: {ref}, Alt: {alt}",
                        file=sys.stderr,
                    )
                    sys.exit(1)

                for row in records:
                    r_items = str(row).split("\t")
                    if ref == r_items[3] and alt == r_items[4]:
                        info = (
                            f"haplotype={haplotype};contig={contig};"
                            f"pos={end};ref={ref};alt={alt}"
                        )
                        print(
                            "\t".join(r_items[:7]),
                            info,
                            "\t".join(r_items[8:]),
                            sep="\t",
                            file=w,
                        )
                        break


def arg_parser():
    parser = argparse.ArgumentParser(
        prog="bed2vcf.py",
        description="Convert *.haplotyped.bed INDEL rows into a VCF for transanno liftvcf",
    )
    parser.add_argument("--input_bed_file", "-i", required=True,
                        help="14-column haplotyped.bed (INDEL rows only)")
    parser.add_argument("--input_vcf_file", "-v", required=True,
                        help="Original per-tool tabix-indexed somatic VCF (clairs indel.vcf.gz or deepsomatic output.vcf.gz)")
    parser.add_argument("--output_vcf_file", "-o", required=True,
                        help="Output VCF file")
    return parser


if __name__ == "__main__":
    parser = arg_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    convert_main(args)
