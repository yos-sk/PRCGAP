#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post-process reference tables for male samples.
Consolidate chrX to one haplotype and chrY to the other.

Overlapping sex-chromosome records are resolved by real alignment support
inside the disputed interval, read from the rmsec PAFs (--paf), not by
envelope size: a chain stretched by satellite matches can claim most of chrY
while holding almost no alignment over the stretch it disputes with the true
contig, so preferring the larger span there deletes the real contig. The
loser is dropped only when it lies inside the winner; otherwise it is clipped
back to its own supported extent. Without --paf the larger-span rule is used.
"""

import sys
import csv
import argparse
from collections import defaultdict


def load_table(path):
    records = []
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            records.append(row)
    return records


def load_paf(paths, min_mapq):
    """(chrom, contig) -> records. PAF query is chm13, target is the contig,
    matching the minimap2 orientation the ref.tables were built from."""
    idx = defaultdict(list)
    for path in paths:
        with open(path) as f:
            for line in f:
                a = line.split('\t')
                if len(a) < 12 or int(a[11]) < min_mapq:
                    continue
                idx[(a[0], a[5])].append((int(a[2]), int(a[3]), int(a[7]), int(a[8]),
                                          int(a[9]), a[4]))
    return idx


def ref_span(records, chrom):
    return sum(int(r[6]) - int(r[5]) for r in records if r[4] == chrom)


def row_support(idx, row, lo, hi):
    """Matched bases of the row's contig aligned to its chromosome inside
    [lo, hi), records clipped proportionally and restricted to the row's
    contig span."""
    cs, ce = int(row[1]), int(row[2])
    tot = 0.0
    for ts, te, qs, qe, m, _ in idx.get((row[4], row[0]), []):
        if qe <= cs or qs >= ce:
            continue
        ov = min(te, hi) - max(ts, lo)
        if ov > 0:
            tot += m * ov / (te - ts)
    return tot


def clip_row_refside(row, keep_lo, keep_hi, idx):
    """Clip the row's reference claim to [keep_lo, keep_hi). The clipped side
    is recomputed from the row's own supporting records so no empty envelope
    survives up to the cut; the untouched side keeps its original (possibly
    telomere-snapped) coordinates. Returns a new record or None."""
    contig, strand, chrom = row[0], row[3], row[4]
    cs, ce, rs, re = int(row[1]), int(row[2]), int(row[5]), int(row[6])
    frs = fre = fcs = fce = None
    total = 0
    for ts, te, qs, qe, m, st in idx.get((chrom, contig), []):
        if qe <= cs or qs >= ce:
            continue
        lo2, hi2 = max(ts, keep_lo), min(te, keep_hi)
        if hi2 <= lo2:
            continue
        span_t = te - ts
        q_len = qe - qs
        cut_l = round((lo2 - ts) / span_t * q_len)
        cut_r = round((te - hi2) / span_t * q_len)
        q_lo = qs + (cut_l if st == "+" else cut_r)
        q_hi = qe - (cut_r if st == "+" else cut_l)
        total += round(m * (hi2 - lo2) / span_t)
        frs = lo2 if frs is None else min(frs, lo2)
        fre = hi2 if fre is None else max(fre, hi2)
        fcs = q_lo if fcs is None else min(fcs, q_lo)
        fce = q_hi if fce is None else max(fce, q_hi)
    if frs is None:
        return None
    if keep_lo <= rs:  # start side untouched
        frs = rs
        if strand == "+":
            fcs = cs
        else:
            fce = ce
    if keep_hi >= re:  # end side untouched
        fre = re
        if strand == "+":
            fce = ce
        else:
            fcs = cs
    if fre <= frs or fce <= fcs:
        return None
    out = list(row)
    out[1], out[2], out[5], out[6] = str(fcs), str(fce), str(frs), str(fre)
    if len(out) > 7:
        out[7] = str(total)
    return out


def _spans_physically(row, lo, hi, idx, lo_ratio=0.5, hi_ratio=2.0):
    """Does the row's contig physically continue across [lo, hi)? True when
    its alignment anchors on both sides sit a contig distance apart
    comparable to the reference distance: a length-conserved, diverged
    segment (assembly-graph branch), not a hole. Returns the ratio or None."""
    rs, re = int(row[5]), int(row[6])
    left = clip_row_refside(row, rs, lo, idx)
    right = clip_row_refside(row, hi, re, idx)
    if left is None or right is None:
        return None
    ref_gap = int(right[5]) - int(left[6])
    if row[3] == "+":
        contig_gap = int(right[1]) - int(left[2])
    else:
        contig_gap = int(left[1]) - int(right[2])
    if ref_gap <= 0 or contig_gap <= 0:
        return None
    ratio = contig_gap / ref_gap
    return ratio if lo_ratio <= ratio <= hi_ratio else None


def _piece_matched(p):
    return int(p[7]) if len(p) > 7 else int(p[2]) - int(p[1])


def _resolve_loser(loser, lo, hi, idx):
    """The loser's claim on [lo, hi) is removed. Inside the winner -> gone
    (0 pieces); partial overlap -> its non-overlapping side (1); containing
    the winner -> both flanks (2). Two flanks whose CONTIG ranges mostly
    overlap are the same contig segment multi-mapped to two reference
    locations (a satellite contig end), which breaks contig linearity --
    such an end can back both a small spurious claim and the real row, so
    only the placement with more aligned bases survives.
    Returns the surviving pieces."""
    l_rs, l_re = int(loser[5]), int(loser[6])
    pieces = []
    if l_rs < lo:
        p = clip_row_refside(loser, l_rs, lo, idx)
        if p is not None:
            pieces.append(p)
    if l_re > hi:
        p = clip_row_refside(loser, hi, l_re, idx)
        if p is not None:
            pieces.append(p)
    if len(pieces) == 2:
        a, b = pieces
        ov = min(int(a[2]), int(b[2])) - max(int(a[1]), int(b[1]))
        min_span = min(int(a[2]) - int(a[1]), int(b[2]) - int(b[1]))
        if min_span > 0 and ov > 0.5 * min_span:
            keep = max(pieces, key=_piece_matched)
            drop = a if keep is b else b
            print(f"[contig-linearity] {loser[0]}: flanks share contig "
                  f"{max(int(a[1]), int(b[1])):,}-{min(int(a[2]), int(b[2])):,}; "
                  f"dropping ref {drop[5]}-{drop[6]} (matched {_piece_matched(drop):,}) "
                  f"in favour of ref {keep[5]}-{keep[6]} (matched {_piece_matched(keep):,})",
                  file=sys.stderr)
            pieces = [keep]
    return pieces


def resolve_overlapping_support(records, idx):
    """Same-chromosome overlaps go to the record with more real alignment
    inside the disputed interval (ties: larger span, then contig name).
    Clipped pieces re-enter the queue and are re-checked."""
    by_chrom = defaultdict(list)
    for r in records:
        by_chrom[r[4]].append(r)

    out = []
    for chrom in sorted(by_chrom):
        queue = sorted(by_chrom[chrom],
                       key=lambda x: (-(int(x[6]) - int(x[5])), x[0]))
        kept = []
        while queue:
            r = queue.pop(0)
            alive = True
            i = 0
            while i < len(kept):
                k = kept[i]
                lo = max(int(k[5]), int(r[5]))
                hi = min(int(k[6]), int(r[6]))
                if hi <= lo:
                    i += 1
                    continue
                s_k = row_support(idx, k, lo, hi)
                s_r = row_support(idx, r, lo, hi)
                if (s_r, int(r[6]) - int(r[5]), r[0]) > (s_k, int(k[6]) - int(k[5]), k[0]):
                    winner, loser = r, k
                else:
                    winner, loser = k, r
                # A loser physically continuous across the interval carries a
                # diverged branch of the locus: keep it whole, drop the copy.
                filled = None
                if int(loser[5]) < lo and int(loser[6]) > hi:
                    filled = _spans_physically(loser, lo, hi, idx)
                    if filled is not None:
                        winner, loser = loser, winner
                pieces = _resolve_loser(loser, lo, hi, idx)
                outcome = (f"{loser[0]} dropped" if not pieces else
                           f"{loser[0]} clipped to "
                           + ", ".join(f"{p[5]}-{p[6]}" for p in pieces))
                how = (f"physically continuous (ratio {filled:.2f})" if filled is not None
                       else f"support {(s_r if winner is r else s_k)/1e6:.3f} Mb beats "
                            f"{loser[0]} {(s_k if winner is r else s_r)/1e6:.3f} Mb")
                print(f"[overlap] {chrom} {lo:,}-{hi:,}: {winner[0]} {how} -> {outcome}",
                      file=sys.stderr)
                if loser is k:
                    kept[i:i + 1] = pieces
                    if not pieces:
                        continue
                    i += len(pieces)
                else:
                    queue.extend(pieces)
                    alive = False
                    break
            if alive:
                kept.append(r)
        out.extend(kept)
    return out


def remove_overlapping(records):
    """Legacy rule (no --paf): remove records that overlap with a larger
    record on the same chrom."""
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


def write_table(records, path, idx):
    # Apply overlap resolution only to sex chromosomes
    sex_chrom = [r for r in records if r[4] in ('chrX', 'chrY')]
    autosome = [r for r in records if r[4] not in ('chrX', 'chrY')]
    if idx is not None:
        sex_chrom = resolve_overlapping_support(sex_chrom, idx)
    else:
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
    parser.add_argument('--paf', nargs='*', default=None,
                        help="rmsec PAFs of both haps (minimap2 -cx asm5 <assembly> "
                             "<chm13>, secondaries removed); enables support-based "
                             "overlap resolution")
    parser.add_argument('--min-mapq', type=int, default=30)
    args = parser.parse_args()

    hap1 = load_table(args.hap1)
    hap2 = load_table(args.hap2)
    idx = load_paf(args.paf, args.min_mapq) if args.paf else None

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
        write_table(new_chrX_hap, args.out1, idx)
        write_table(new_chrY_hap, args.out2, idx)
    else:
        write_table(new_chrY_hap, args.out1, idx)
        write_table(new_chrX_hap, args.out2, idx)

    chrY_total = ref_span(load_table(args.out1 if chrY_label == 'hap1' else args.out2), 'chrY')
    print(f"chrY total span in {chrY_label}: {chrY_total:,} bp", file=sys.stderr)


if __name__ == '__main__':
    main()
