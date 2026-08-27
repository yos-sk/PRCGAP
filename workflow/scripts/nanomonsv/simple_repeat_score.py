#!/usr/bin/env python3
"""Score each SV's inserted sequence against the reference flank and against itself.

Flagging every call inside a tandem repeat removes true variants that merely sit
in one. What should be removed is a change in the array's own copy number, and an
expansion of the local array is by definition a copy of sequence next to the
breakpoint -- so the inserted sequence is aligned back to the flank rather than
compared with a consensus motif. The motif route does not work: one ULTRA
interval carries a single consensus even where it spans kilobases of
heterogeneous array, so the consensus rarely matches the unit at the breakpoint.

Two independent scores per call:

flank_coverage / flank_aln_identity
    fraction of the inserted sequence matched by a local alignment to one flank,
    and the identity of that alignment. Local rather than end-to-end because
    diverged VNTR copies only line up over part of their length.

    Same strand only. An expansion of an array reading (T)n on the plus strand
    inserts more (T)n, not (A)n; allowing the reverse complement only invents
    matches, and is how a poly-A tail on the opposite strand comes to look like
    the local array.

    The gap penalties have to be strict. At gap open 4 / extend 2 against a +5
    match a gap costs barely more than the match it buys, so in A/T-rich sequence
    the aligner strings matches together through gaps and unrelated sequence
    reaches coverage 0.85 at identity 0.50. Under the defaults here that same
    sequence scores 0.01-0.12 while a real extra copy stays above 0.6.

self_repeat
    the largest fraction of positions equal to the base one period ahead, over
    all periods. Catches an inserted sequence that is a tandem repeat in its own
    right even when the array it came from is not in the flank window, and unlike
    ULTRA it still works on a 36 bp insertion. Also reported for the stretch of
    flank the insertion aligned to (flank_self_repeat): the flank criterion only
    means "one more copy of the array" if what it matched is an array, which is
    what keeps a full-length mobile element that landed beside a homologous copy
    from being filtered.

Derived from workspace/repeat_bench Experiment 10; see that README for the IGV
verification behind the penalties and thresholds.
"""

import argparse
import csv
import sys

import parasail
import pysam

try:                                  # only used for inserted sequences too long
    import edlib                      # to hold a traceback matrix; ~1 call in 2000
except ImportError:                   # falls back to parasail statistics
    edlib = None


def _long_identity(query, target, matrix, gap_open, gap_extend):
    """Coverage/identity for a query too long for a traceback alignment."""
    if not query or not target:
        return 0.0, 0.0, ""
    if edlib is not None:
        r = edlib.align(query, target, mode="HW", task="distance")
        ident = max(0.0, 1.0 - r["editDistance"] / len(query))
        return ident, ident, target
    # parasail's stats variants report matches without building the traceback,
    # so memory stays O(query) instead of the O(query x target) a traceback
    # matrix would need. The _sat variant re-runs at a wider integer width if
    # the score saturates, which a long near-identical insertion could do.
    r = parasail.sw_stats_striped_sat(query, target, gap_open, gap_extend, matrix)
    cov = r.matches / len(query)
    ident = r.matches / r.length if r.length else 0.0
    block = target[max(r.end_ref - r.length + 1, 0):r.end_ref + 1]
    return cov, ident, block


def sw_coverage(query, target, matrix, gap_open, gap_extend):
    """Query coverage, identity and the matched stretch of target."""
    if not query or not target:
        return 0.0, 0.0, ""
    r = parasail.sw_trace_striped_16(query, target, gap_open, gap_extend, matrix)
    tb = r.traceback
    matches, columns = tb.comp.count("|"), len(tb.comp)
    block = len(tb.ref.replace("-", ""))
    return (matches / len(query), (matches / columns if columns else 0.0),
            target[max(r.end_ref - block + 1, 0):r.end_ref + 1])


def self_repeat(seq, max_period=200):
    """Largest fraction of positions equal to the base one period ahead."""
    best = (0.0, 0)
    for p in range(1, min(max_period, len(seq) // 2) + 1):
        m = sum(1 for i in range(len(seq) - p) if seq[i] == seq[i + p])
        f = m / (len(seq) - p)
        if f > best[0]:
            best = (f, p)
    return best


def _fetch(fa, chrom, start, end):
    ctg_len = fa.get_reference_length(chrom)
    clip = lambda x: min(max(x, 0), ctg_len)
    return fa.fetch(chrom, clip(start), clip(end)).upper()


def flank_support(fa, chrom1, pos1, chrom2, pos2, seq,
                  max_sw_len, matrix, gap_open, gap_extend):
    # Each breakpoint brings its own flank, so an inter-contig call is handled
    # the same way as a local one.
    pad = 2 * len(seq) + 100
    left = _fetch(fa, chrom1, pos1 - pad, pos1)
    right = _fetch(fa, chrom2, pos2 - 1, pos2 - 1 + pad)

    best = (0.0, "-", 0.0, "")
    for side, target in (("left", left), ("right", right)):
        if len(seq) > max_sw_len:
            cov, ident, block = _long_identity(seq, target, matrix,
                                               gap_open, gap_extend)
        else:
            cov, ident, block = sw_coverage(seq, target, matrix,
                                            gap_open, gap_extend)
        if cov > best[0]:
            best = (cov, side, ident, block)
    return best


class _Assembly:
    """Both haplotype FASTAs behind one contig -> sequence lookup."""

    def __init__(self, paths):
        self._fa = [pysam.FastaFile(p) for p in paths if p]
        self._of = {}
        for fa in self._fa:
            for ref in fa.references:
                self._of.setdefault(ref, fa)

    def get_reference_length(self, chrom):
        return self._of[chrom].get_reference_length(chrom)

    def fetch(self, chrom, start, end):
        return self._of[chrom].fetch(chrom, start, end)


def main(argv=None):
    p = argparse.ArgumentParser(prog="simple_repeat_score")
    p.add_argument("sv_file", help="nanomonsv result table")
    p.add_argument("output_file")
    p.add_argument("assembly", nargs="+", help="haplotype FASTA(s)")
    p.add_argument("--max-sw-len", type=int, default=20000,
                   help="fall back to a statistics-only alignment above this "
                        "inserted length (traceback memory)")
    p.add_argument("--match", type=int, default=2)
    p.add_argument("--mismatch", type=int, default=-6)
    p.add_argument("--gap-open", type=int, default=10)
    p.add_argument("--gap-extend", type=int, default=1)
    a = p.parse_args(argv)

    fa = _Assembly(a.assembly)
    matrix = parasail.matrix_create("ACGT", a.match, a.mismatch)

    with open(a.sv_file) as fh, open(a.output_file, "w") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["SV_ID", "ins_len", "flank_coverage", "best_flank",
                    "flank_aln_identity", "self_repeat", "self_period",
                    "flank_self_repeat"])
        for F in csv.DictReader(fh, delimiter="\t"):
            seq = F["Inserted_Seq"]
            if seq in ("", "---"):
                continue
            seq = seq.upper()
            cov, side, ident, block = flank_support(
                fa, F["Chr_1"], int(F["Pos_1"]), F["Chr_2"], int(F["Pos_2"]),
                seq, a.max_sw_len, matrix, a.gap_open, a.gap_extend)
            srep, sper = self_repeat(seq)
            fsrep = self_repeat(block)[0] if len(block) >= 4 else 0.0
            w.writerow([F["SV_ID"], len(seq), f"{cov:.3f}", side, f"{ident:.3f}",
                        f"{srep:.3f}", sper, f"{fsrep:.3f}"])

    if edlib is None:
        print("[simple_repeat_score] edlib not available; inserted sequences "
              "over {:,} bp scored with parasail statistics".format(a.max_sw_len),
              file=sys.stderr)


if __name__ == "__main__":
    main()
