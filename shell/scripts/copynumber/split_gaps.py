#!/usr/bin/env python3

# Clip CBS segments that span assembly gaps.
#
# The copy-number index axis is built by copynumber_window.py: per c_chrom the
# index advances +1 per WINDOW_SIZE window inside each contig, and jumps by
# (c_start - prev_c_end)/WINDOW_SIZE between consecutive contigs of the same
# c_chrom. A large jump is a gap region with no assembled sequence. Rather than
# re-detecting those jumps from the (already binned) copynumber.tsv, we replay
# copynumber_window.py's index loop directly on the ref.table so the gap
# coordinates come from the authoritative source.

import argparse

WINDOW_SIZE = 50000


def build_gaps(ref_table: str, min_gap_bins: float) -> dict:
    """Replay copynumber_window.py's index loop and collect gap intervals
    (gap_lo, gap_hi) on the index axis, keyed by c_chrom (e.g. "chr1")."""
    gap_db = dict()
    diff = 0.0
    prev_chrom = ""
    prev_end = 0
    prev_index = 0
    with open(ref_table, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            p_start = int(items[1])
            p_end = int(items[2])
            c_chrom = items[4]
            c_start = int(items[5])
            c_end = int(items[6])

            if prev_chrom != c_chrom:
                diff = c_start / WINDOW_SIZE
                prev_index = 0
            else:
                # Gap between the previous contig and this one on the index axis:
                # gap_lo = last window of the previous contig,
                # gap_hi = first window of this contig.
                gap_lo = diff + prev_index - 1
                diff += (c_start - prev_end) / WINDOW_SIZE
                gap_hi = diff + prev_index
                if gap_hi - gap_lo > min_gap_bins:
                    gap_db.setdefault(c_chrom, []).append((gap_lo, gap_hi))

            prev_chrom = c_chrom
            prev_end = c_end

            n_windows = len(range(p_start, p_end, WINDOW_SIZE))
            prev_index += n_windows

    for chrom in gap_db:
        gap_db[chrom].sort()
    return gap_db


def subtract_gaps(start: float, end: float, gaps: list) -> list:
    """Remove every gap interval from [start, end], returning the remaining
    pieces. Handles partial overlap, multiple gaps per segment, and drops
    pieces (or whole segments) that fall entirely inside a gap."""
    pieces = [(start, end)]
    for g_lo, g_hi in gaps:
        new_pieces = []
        for s, e in pieces:
            if e <= g_lo or s >= g_hi:
                new_pieces.append((s, e))
                continue
            if s < g_lo:
                new_pieces.append((s, g_lo))
            if g_hi < e:
                new_pieces.append((g_hi, e))
        pieces = new_pieces
    return pieces


def split_gaps(ref_table: str, cbs_result: str, min_gap_bins: float) -> None:
    gap_db = build_gaps(ref_table, min_gap_bins)

    with open(cbs_result, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            chrom = "chr" + items[1]
            start = float(items[2])
            end = float(items[3])

            if chrom not in gap_db:
                print(line.rstrip("\n"))
                continue

            pieces = subtract_gaps(start, end, gap_db[chrom])
            if pieces == [(start, end)]:
                print(line.rstrip("\n"))
                continue

            for piece_start, piece_end in pieces:
                print(items[0], items[1], piece_start, piece_end,
                      items[4], items[5], sep="\t")


def create_parser():
    parser = argparse.ArgumentParser(prog="split_gaps.py")
    parser.add_argument("ref_table", help="The PATH to the reference table.")
    parser.add_argument("cbs_result", help="The PATH to the CBS result.")
    parser.add_argument("--min-gap-bins", type=float, default=5,
                        help="Register a gap when its width on the index axis "
                             "exceeds this many bins (default: 5). Matches the "
                             "GAP_THRESHOLD used by the cen/sat annotation in "
                             "plot_copy_number.R so CBS clipping and the "
                             "annotation share the same gap regions.")
    return parser.parse_args()


def main():
    args = create_parser()
    split_gaps(args.ref_table, args.cbs_result, args.min_gap_bins)


if __name__ == "__main__":
    main()
