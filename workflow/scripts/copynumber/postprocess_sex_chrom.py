#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post-process reference tables for male samples.
Consolidate chrX to one haplotype and chrY to the other.
"""

import sys
import csv
import argparse


def load_table(path):
    records = []
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            records.append(row)
    return records


def ref_span(records, chrom):
    return sum(int(r[6]) - int(r[5]) for r in records if r[4] == chrom)


def remove_overlapping(records):
    """Remove records that overlap with a larger record on the same chrom."""
    sorted_records = sorted(records, key=lambda x: (x[4], int(x[5])))
    keep = [True] * len(sorted_records)
    for i in range(len(sorted_records)):
        if not keep[i]:
            continue
        chrom_i = sorted_records[i][4]
        start_i = int(sorted_records[i][5])
        end_i = int(sorted_records[i][6])
        span_i = end_i - start_i
        for j in range(i + 1, len(sorted_records)):
            if sorted_records[j][4] != chrom_i:
                break
            start_j = int(sorted_records[j][5])
            end_j = int(sorted_records[j][6])
            if start_j >= end_i:
                break
            # Records overlap; remove the smaller one
            span_j = end_j - start_j
            if span_i >= span_j:
                keep[j] = False
            else:
                keep[i] = False
                break
    return [r for r, k in zip(sorted_records, keep) if k]


def write_table(records, path):
    # Apply overlap removal only to sex chromosomes
    sex_chrom = [r for r in records if r[4] in ('chrX', 'chrY')]
    autosome = [r for r in records if r[4] not in ('chrX', 'chrY')]
    sex_chrom = remove_overlapping(sex_chrom)
    records = autosome + sex_chrom
    with open(path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for r in sorted(records, key=lambda x: (x[4], int(x[5]))):
            writer.writerow(r)


def main():
    parser = argparse.ArgumentParser(prog="postprocess_sex_chrom.py")
    parser.add_argument('--hap1', required=True, help="Path to hap1 ref.table")
    parser.add_argument('--hap2', required=True, help="Path to hap2 ref.table")
    parser.add_argument('--out1', required=True, help="Output path for hap1 ref.table")
    parser.add_argument('--out2', required=True, help="Output path for hap2 ref.table")
    args = parser.parse_args()

    hap1 = load_table(args.hap1)
    hap2 = load_table(args.hap2)

    # Determine which haplotype has chrX
    chrX_span_1 = ref_span(hap1, 'chrX')
    chrX_span_2 = ref_span(hap2, 'chrX')

    if chrX_span_1 >= chrX_span_2:
        chrX_hap, chrY_hap = hap1, hap2
        chrX_label, chrY_label = 'hap1', 'hap2'
    else:
        chrX_hap, chrY_hap = hap2, hap1
        chrX_label, chrY_label = 'hap2', 'hap1'

    print(f"chrX assigned to {chrX_label} (span: {max(chrX_span_1, chrX_span_2):,} bp)", file=sys.stderr)
    print(f"chrY assigned to {chrY_label}", file=sys.stderr)

    # Move chrY from chrX_hap to chrY_hap, move chrX from chrY_hap to chrX_hap
    chrY_from_chrX_hap = [r for r in chrX_hap if r[4] == 'chrY']
    chrX_from_chrY_hap = [r for r in chrY_hap if r[4] == 'chrX']

    new_chrX_hap = [r for r in chrX_hap if r[4] != 'chrY'] + chrX_from_chrY_hap
    new_chrY_hap = [r for r in chrY_hap if r[4] != 'chrX'] + chrY_from_chrX_hap

    if chrX_span_1 >= chrX_span_2:
        write_table(new_chrX_hap, args.out1)
        write_table(new_chrY_hap, args.out2)
    else:
        write_table(new_chrY_hap, args.out1)
        write_table(new_chrX_hap, args.out2)

    chrY_total = ref_span(new_chrY_hap, 'chrY')
    print(f"chrY total span in {chrY_label}: {chrY_total:,} bp", file=sys.stderr)


if __name__ == '__main__':
    main()
