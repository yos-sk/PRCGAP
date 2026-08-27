#!/usr/bin/env python3
"""ref.table → the 50 kb window BED that mosdepth --by should measure.

Emits exactly the windows copynumber_window.py iterates over, coordinate-sorted
so mosdepth can stream them. Window ORDER and the fractional window index are
not encoded here: windows_to_tsv.py re-walks the same ref.table to recover both,
so nothing depends on how mosdepth orders or names its output regions.
"""

import argparse

WINDOW_SIZE = 50000


def windows_of(p_start, p_end, strand):
    """The windows tiling [p_start, p_end), in copynumber_window.py's order."""
    if strand == "+":
        return [(s, min(s + WINDOW_SIZE, p_end))
                for s in range(p_start, p_end, WINDOW_SIZE)]
    return [(max(p_start, e - WINDOW_SIZE), e)
            for e in range(p_end, p_start, -WINDOW_SIZE)]


def main():
    parser = argparse.ArgumentParser(prog="make_windows.py")
    parser.add_argument("reftable")
    args = parser.parse_args()

    rows = []
    with open(args.reftable) as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            p_chrom, p_start, p_end = items[0], int(items[1]), int(items[2])
            for w_start, w_end in windows_of(p_start, p_end, items[3]):
                rows.append((p_chrom, w_start, w_end))

    for p_chrom, w_start, w_end in sorted(rows):
        print(p_chrom, w_start, w_end, sep="\t")


if __name__ == "__main__":
    main()
