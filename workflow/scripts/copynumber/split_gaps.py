#!/usr/bin/env python3

import sys
import numpy as np

def split_gaps(input_file: str, cbs_result: str) -> None:
    gap_db = dict()
    prev_contig = ""
    prev_index = 0
    with open(input_file, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            contig = items[1]
            index = float(items[2])

            if prev_contig == "" and prev_index == 0:
                prev_contig = contig
                prev_index = index
                continue

            if index - prev_index > 1000:
                if items[0] in gap_db:
                    gap_db[items[0]].append((prev_index, index))
                else:
                    gap_db[items[0]] = [(prev_index, index)]
            
            prev_contig = contig
            prev_index = index

    with open(cbs_result, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            chrom = "chr" + items[1]
            start = float(items[2])
            end = float(items[3])

            flag = False
            if chrom in gap_db:
                for interval in gap_db[chrom]:
                    if start <= interval[0] and interval[1] <= end:
                        print(items[0], items[1], start, interval[0], items[4], items[5], sep="\t")
                        print(items[0], items[1], interval[1], end, items[4], items[5], sep="\t")
                        flag = True
                        break

            if not flag:
                print(line.rstrip("\n"))
                

def main():
    input_file = sys.argv[1]
    cbs_result = sys.argv[2]
    split_gaps(input_file, cbs_result)

if __name__ == "__main__":
    main()

