#!/usr/bin/env python3

import sys
import pysam
import argparse

def sum_depth_window(args, window_size=50000):
    tbx_file_tumor = pysam.TabixFile(args.tumor)
    tbx_file_normal = pysam.TabixFile(args.normal)

    diff = 0.0
    prev_chrom = ""
    prev_end = 0
    prev_index = 0
    with open(args.reftable, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            p_chrom = items[0]
            p_start = int(items[1])
            p_end = int(items[2])
            strand = items[3]
            c_chrom = items[4]
            c_start = int(items[5])
            c_end = int(items[6])

            if prev_chrom != c_chrom:
                diff = c_start / window_size
                prev_index = 0
            else:
                diff += (c_start - prev_end) / window_size

            prev_chrom = c_chrom
            prev_end = c_end

            index = 0
            if strand == "+":
                for start in range(p_start, p_end, window_size):
                    tumor_depth = 0
                    normal_depth = 0
                    try: 
                        for record in tbx_file_tumor.fetch(p_chrom, start, min(start + window_size, p_end)):
                            items = record.rstrip('\n').split('\t')
                            tumor_depth += int(items[3])
                    except:
                        print(f'Could not fetch by tumor tabix:{p_chrom}:{start}-{min(start + window_size, p_end)}', file = sys.stderr)
                        
                    try: 
                        for record in tbx_file_normal.fetch(p_chrom, start, min(start + window_size, p_end)):
                            items = record.rstrip('\n').split('\t')
                            normal_depth += int(items[3])
                    except:
                        print(f'Could not fetch by tumor tabix:{p_chrom}:{start}-{min(start + window_size, p_end)}', file = sys.stderr)

                    print(c_chrom, p_chrom, diff + prev_index + index, tumor_depth, normal_depth, sep="\t") 
                    index += 1
            else:
                for end in range(p_end,  p_start, -window_size):
                    tumor_depth = 0
                    normal_depth = 0
                    try: 
                        for record in tbx_file_tumor.fetch(p_chrom, max(p_start, end - window_size), end):
                            items = record.rstrip('\n').split('\t')
                            tumor_depth += int(items[3])
                    except:
                        print(f'Could not fetch by tumor tabix:{p_chrom}:{max(end - window_size, p_start)}-{end}', file = sys.stderr)

                    try: 
                        for record in tbx_file_normal.fetch(p_chrom, max(p_start, end - window_size), end):
                            items = record.rstrip('\n').split('\t')
                            normal_depth += int(items[3])
                    except:
                        print(f'Could not fetch by tumor tabix:{p_chrom}:{max(end - window_size, p_start)}-{end}', file = sys.stderr)

                    print(c_chrom, p_chrom, diff + prev_index + index, tumor_depth, normal_depth, sep="\t") 
                    index += 1
            prev_index += index

def create_parser():
    parser = argparse.ArgumentParser(prog = "copynumber_window.py")
    parser.add_argument('-t', '--tumor', help="The PATH to the tumor depth file.")
    parser.add_argument('-n', '--normal', help="The PATH to the normal depth file.")
    parser.add_argument('-r', '--reftable', help="The PATH to the reference table.")

    args = parser.parse_args()

    return args

def main():
    args = create_parser()

    sum_depth_window(args)

if __name__ == "__main__":
    main()

