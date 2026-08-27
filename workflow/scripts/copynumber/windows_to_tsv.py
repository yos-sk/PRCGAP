#!/usr/bin/env python3
"""mosdepth --by region means → the copynumber.tsv copynumber_window.py emits.

Columns match the existing format exactly:
    c_chrom, p_chrom, window_index, tumor_depth, normal_depth

cbs.R applies absolute thresholds to the depth columns (normal_depth > 100000,
tumor_depth > 10000 over a 50 kb window), so they have to stay SUMS of per-base
depth, not means. mosdepth reports a region mean rounded to 2 decimals, so the
sum is recovered as round(mean * window_length) — quantized in steps of
0.005 * window_length (250 for a full 50 kb window).
"""

import argparse
import gzip

WINDOW_SIZE = 50000


def load_means(path):
    """(chrom, start, end) → mean depth, from a mosdepth regions.bed.gz."""
    means = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            # A 4th column is the region name when the input BED carried one;
            # the mean is always the last field.
            means[(items[0], int(items[1]), int(items[2]))] = float(items[-1])
    return means


def window_sum(means, chrom, start, end, missing):
    mean = means.get((chrom, start, end))
    if mean is None:
        missing.append((chrom, start, end))
        return 0
    return int(round(mean * (end - start)))


def main():
    parser = argparse.ArgumentParser(prog="windows_to_tsv.py")
    parser.add_argument("-t", "--tumor", required=True,
                        help="tumor mosdepth regions.bed.gz")
    parser.add_argument("-n", "--normal", required=True,
                        help="normal mosdepth regions.bed.gz")
    parser.add_argument("-r", "--reftable", required=True)
    args = parser.parse_args()

    tumor = load_means(args.tumor)
    normal = load_means(args.normal)
    missing = []

    diff = 0.0
    prev_chrom = ""
    prev_end = 0
    prev_index = 0
    with open(args.reftable) as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            p_chrom, p_start, p_end = items[0], int(items[1]), int(items[2])
            strand = items[3]
            c_chrom, c_start, c_end = items[4], int(items[5]), int(items[6])

            if prev_chrom != c_chrom:
                diff = c_start / WINDOW_SIZE
                prev_index = 0
            else:
                diff += (c_start - prev_end) / WINDOW_SIZE

            prev_chrom = c_chrom
            prev_end = c_end

            if strand == "+":
                bounds = [(s, min(s + WINDOW_SIZE, p_end))
                          for s in range(p_start, p_end, WINDOW_SIZE)]
            else:
                bounds = [(max(p_start, e - WINDOW_SIZE), e)
                          for e in range(p_end, p_start, -WINDOW_SIZE)]

            for index, (w_start, w_end) in enumerate(bounds):
                print(c_chrom, p_chrom, diff + prev_index + index,
                      window_sum(tumor, p_chrom, w_start, w_end, missing),
                      window_sum(normal, p_chrom, w_start, w_end, missing),
                      sep="\t")
            prev_index += len(bounds)

    if missing:
        import sys
        print("windows_to_tsv.py: {} windows absent from the mosdepth output, "
              "e.g. {}".format(len(missing), missing[:3]), file=sys.stderr)


if __name__ == "__main__":
    main()
