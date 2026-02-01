#! /usr/bin/env python3

import sys
import gzip
import pysam
import argparse

def extract_seq(input_fasta: str, input_bed: str, output_fasta: str) -> None:
    fasta_reader = pysam.FastaFile(input_fasta)

    with open(input_bed, "r") as f, open(output_fasta, "w") as w:
        for line in f:
            items = line.rstrip("\n").split("\t")
            chrom = items[0]
            start = int(items[1]) 
            end = int(items[2])
            mut_id = items[11]
            ref_len = fasta_reader.get_reference_length(chrom)

            bstart, bend = max(start - 100, 1), min(end + 100, ref_len)
            #print(chrom, bstart, bend, file=sys.stderr)

            seq = fasta_reader.fetch(chrom, bstart, bend)

            print(f'>{mut_id}', file=w)
            print(seq, file=w)

    fasta_reader.close()

def convert_bed(input_vcf: str, output_bed: str, tool: str) -> None:
    index = 0
    with gzip.open(input_vcf, "rt") as f, open(output_bed, "w") as w:
        for line in f:
            if line[0] == "#": continue
            items = line.rstrip("\n").split("\t")
            chrom = items[0]
            pos = int(items[1])
            ref = items[3]
            alt = items[4]
            qual = items[5]
            filt = items[6]
            
            if tool == "ClairS":
                tumor_vaf = items[9].split(":")[3]
                tumor_cov = items[9].split(":")[2]
                normal_vaf = items[9].split(":")[5]
                normal_cov = items[9].split(":")[6]
            elif tool == "DeepSomatic":
                tumor_vaf = items[9].split(":")[4]
                if tumor_vaf == "1":
                    tumor_vaf = "1.0"
                tumor_cov = items[9].split(":")[2]
                normal_vaf = 0.0
                normal_cov = 0
            
            if filt == "PASS":
                if len(alt.split(",")) > 1:
                    l_tumor_vaf = tumor_vaf.split(",")
                    l_alt = alt.split(",")
                    for i in range(0, len(l_alt)):
                        print(chrom, pos-1, pos-1+len(ref), ref, l_alt[i], qual, filt, l_tumor_vaf[i], tumor_cov, normal_vaf, normal_cov, "ID_" + str(index), sep="\t", file=w)
                        index += 1
                else:
                    print(chrom, pos-1, pos-1+len(ref), ref, alt, qual, filt, tumor_vaf, tumor_cov, normal_vaf, normal_cov, "ID_" + str(index), sep="\t", file=w)
                    index += 1

def parse_vcf_main(args):

    convert_bed(args.input_vcf, args.output_bed, args.tool)
    extract_seq(args.input_fasta, args.output_bed, args.output_fasta)
    

def arg_parser():
    parser = argparse.ArgumentParser(prog = "parse_vcf",
                description = "Parse VCF file from ClairS/DeepSomatic")

    parser.add_argument("--input_vcf", "-i",  type = str,
            help = "Path to the input VCF file")

    parser.add_argument("--input_fasta", "-f",  type = str,
            help = "Path to the input FASTA file")
    
    parser.add_argument("--output_bed", "-o", type = str,
            help = "Path to the output BED file")

    parser.add_argument("--output_fasta", "-g", type = str,
            help = "Path to the output FASTA file")

    parser.add_argument("--tool", "-t", type = str,
            help = "ClairS/DeepSomatic")
    
    parser.set_defaults(func = parse_vcf_main)

    return parser

def main():
    parser = arg_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args.func(args)

if __name__ == "__main__":
    main() 
