#!/usr/bin/env python3
"""Contig -> chm13 chromosome/strand/coordinate table for copy number.

Input: PAF from `minimap2 -cx asm5 <assembly> <chm13>` with secondaries removed,
so the PAF target is the contig and the PAF query is chm13; parse_paf_gt_as_query
swaps them. The assembly is NOT masked -- masking changed no assignment and cost
32.8 Mb of span.

Output, one row per contig unless it is genuinely chimeric:
  contig, contig_start, contig_end, strand, chrom, ref_start, ref_end, matched
Rows are ordered by (chrom, ref_start), not by contig: copynumber_window.py
accumulates the reference gap only between consecutive rows of one chromosome,
so non-adjacent rows of the same chromosome would both start at window index 0.

Stages:
  1+2 per (chrom, strand), the matched-weighted best chain, gated by a seed
  3   merge chains of one chromosome across gaps up to --merge-gap
  4   drop chains below --min-density, and below --min-ref-density on the
      chm13 side
  5   two chromosomes overlapping in contig coordinates claim one physical
      stretch, so one is wrong: the much denser chain wins, else clip both
  7   choose the contig's chromosome(s). An acrocentric's evidence is its chm13
      span minus the censat part, so short-arm homology to the other four
      acrocentrics falls below --split-min
  9   a telomeric repeat at a contig end means that end is a chromosome
      terminus; snap contig and chm13 sides together
  8   drop a row whose reference interval is contained in another row's;
      a PARTIAL overlap goes to the contig with more real alignment inside
      it, and the loser is clipped back to its own supported extent
"""

import sys
import gzip
import argparse
from collections import defaultdict


def _open_text(path):
    """The workflow passes annotation BEDs bgzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_paf_gt_as_query(path):
    rows = []
    with _open_text(path) as f:
        for line in f:
            f_ = line.rstrip("\n").split("\t")
            if len(f_) < 12:
                continue
            chrom, chrom_len = f_[0], int(f_[1])
            chrom_start, chrom_end = int(f_[2]), int(f_[3])
            strand = f_[4]
            contig, contig_len = f_[5], int(f_[6])
            contig_start, contig_end = int(f_[7]), int(f_[8])
            matched = int(f_[9])
            mapq = int(f_[11])
            rows.append(dict(qname=contig, qlen=contig_len,
                              qstart=contig_start, qend=contig_end,
                              strand=strand, tname=chrom, tlen=chrom_len,
                              tstart=chrom_start, tend=chrom_end,
                              matched=matched, mapq=mapq))
    return rows


def weighted_chain_indices(keys, weights):
    """Weighted longest path i->j (i<j, keys[i]<keys[j]) maximizing
    sum(weights); O(n log n) via a Fenwick tree over compressed keys."""
    n = len(keys)
    if n == 0:
        return []
    sorted_keys = sorted(set(keys))
    rank = {v: i + 1 for i, v in enumerate(sorted_keys)}
    m = len(sorted_keys)
    tree_val = [0] * (m + 1)
    tree_idx = [-1] * (m + 1)

    def update(i, val, idx):
        while i <= m:
            if tree_val[i] < val:
                tree_val[i] = val
                tree_idx[i] = idx
            i += i & (-i)

    def query_prefix_max(i):
        best_val, best_idx = 0, -1
        while i > 0:
            if tree_val[i] > best_val:
                best_val, best_idx = tree_val[i], tree_idx[i]
            i -= i & (-i)
        return best_val, best_idx

    dp = [0] * n
    prev = [-1] * n
    for i in range(n):
        r = rank[keys[i]]
        best_val, best_idx = query_prefix_max(r - 1)
        dp[i] = weights[i] + best_val
        prev[i] = best_idx
        update(r, dp[i], i)

    end = max(range(n), key=lambda i: dp[i])
    chain = []
    k = end
    while k != -1:
        chain.append(k)
        k = prev[k]
    chain.reverse()
    return chain


def has_seed(chain_rows, seed_mapq, seed_alnlen):
    return any(r["mapq"] >= seed_mapq and (r["qend"] - r["qstart"]) >= seed_alnlen
               for r in chain_rows)


def chain_from_members(tname, strand, members, truncated_left=False, truncated_right=False):
    members = sorted(members, key=lambda r: r["qstart"])
    return dict(tname=tname, strand=strand, members=members,
                qstart=members[0]["qstart"], qend=max(r["qend"] for r in members),
                tstart=min(r["tstart"] for r in members), tend=max(r["tend"] for r in members),
                matched=sum(r["matched"] for r in members),
                qmember=sum(r["qend"] - r["qstart"] for r in members),
                tmember=sum(r["tend"] - r["tstart"] for r in members),
                n_anchors=len(members),
                truncated_left=truncated_left, truncated_right=truncated_right)


def _density(matched, envelope, member_total):
    """Capped at 1.0 by dividing by the summed member spans when they exceed the
    envelope. Overlapping members otherwise count the same contig bases twice and
    a satellite pile-up reads as density 4.9, which stage 5 believed."""
    denom = max(envelope, member_total)
    return matched / denom if denom > 0 else 0.0


def density(chain):
    return _density(chain["matched"], chain["qend"] - chain["qstart"], chain["qmember"])


def ref_density(chain):
    """density() taken on the chm13 side. A pericentromeric fragment can pass the
    contig-side cutoff while its chain reaches across tens of Mb of reference,
    because the contig span it divides by is tiny -- a fraction of a Mb of contig
    claiming most of a chromosome arm. Such a row overlaps every other row on
    that chromosome, and copynumber_window.py's
    `diff += (c_start - prev_end) / W` then goes negative and walks the window
    axis backwards."""
    return _density(chain["matched"], chain["tend"] - chain["tstart"], chain["tmember"])


# ---------------------------------------------------------------- stage 1+2
def build_group_chains(rows, min_mapq, seed_mapq, seed_alnlen):
    """Per (tname, strand), the matched-weighted-best chain, gated by seed."""
    groups = defaultdict(list)
    for r in rows:
        if r["mapq"] < min_mapq:
            continue
        groups[(r["tname"], r["strand"])].append(r)

    chains = []
    for (tname, strand), members in groups.items():
        members.sort(key=lambda r: r["qstart"])
        key = [(-r["tstart"] if strand == "-" else r["tstart"]) for r in members]
        weight = [r["matched"] for r in members]
        idx = weighted_chain_indices(key, weight)
        if not idx:
            continue
        chain_rows = [members[i] for i in idx]
        if not has_seed(chain_rows, seed_mapq, seed_alnlen):
            continue
        chains.append(chain_from_members(tname, strand, chain_rows))
    chains.sort(key=lambda c: c["qstart"])
    return chains


# --------------------------------------------------------------------- stage 3
def clip_member_outside(member, lo, hi):
    """Fragment(s) of `member` outside query [lo, hi). Clips rather than
    drops the whole record, since one PAF row can span megabases and a small
    edge overlap shouldn't cost all of it. 0-2 fragments."""
    qs, qe = member["qstart"], member["qend"]
    if qe <= lo or qs >= hi:
        return [member]
    span_q = qe - qs
    ts, te = member["tstart"], member["tend"]
    is_minus = member["strand"] == "-"
    frags = []
    if qs < lo:
        frac = (lo - qs) / span_q
        t_cut = (te - round((te - ts) * frac)) if is_minus else (ts + round((te - ts) * frac))
        left = dict(member, qend=lo, matched=round(member["matched"] * frac))
        left["tstart"], left["tend"] = (t_cut, te) if is_minus else (ts, t_cut)
        frags.append(left)
    if qe > hi:
        frac = (qe - hi) / span_q
        t_cut = (ts + round((te - ts) * frac)) if is_minus else (te - round((te - ts) * frac))
        right = dict(member, qstart=hi, matched=round(member["matched"] * frac))
        right["tstart"], right["tend"] = (ts, t_cut) if is_minus else (t_cut, te)
        frags.append(right)
    return frags


def resolve_cross_contig_overlaps(chains, seed_mapq, seed_alnlen, min_density, win_ratio=0.0):
    """Different chromosomes cannot both be right for the same contig stretch.
    Containment -> keep the better chain. Partial overlap -> the much denser
    chain wins, otherwise clip both at the shared window."""
    chains = sorted(chains, key=lambda c: c["qstart"])
    changed = True
    while changed:
        changed = False
        chains.sort(key=lambda c: c["qstart"])
        containment_removed = False
        for i in range(len(chains)):
            for j in range(i + 1, len(chains)):
                a, b = chains[i], chains[j]
                if a["tname"] == b["tname"]:
                    continue
                if _contains(a, b) or _contains(b, a):
                    loser = b if _containment_winner(a, b) is a else a
                    chains = [c for c in chains if c is not loser]
                    changed = True
                    containment_removed = True
                    break
                lo, hi = max(a["qstart"], b["qstart"]), min(a["qend"], b["qend"])
                if lo >= hi:
                    continue
                # Clipping both sides treats the two chains as equally credible.
                # Overlapping in contig coordinates means one physical stretch
                # is claimed by two chromosomes, so one of them is wrong -- a
                # real chimera puts its arms on DIFFERENT parts of the contig
                # and does not overlap at all.
                #
                # Density, not matched bases, says which is wrong: the homology
                # chain is sparse over the stretch it claims while the real one
                # is solid. Matched bases separate the two far less, and both
                # arms of a real chimera are dense, so a density ratio leaves
                # them alone where a matched-base ratio would drop one.
                if win_ratio > 0:
                    da, db = density(a), density(b)
                    if db > 0 and da >= win_ratio * db:
                        loser = b
                    elif da > 0 and db >= win_ratio * da:
                        loser = a
                    else:
                        loser = None
                    if loser is not None:
                        chains = [c for c in chains if c is not loser]
                        changed = True
                        containment_removed = True
                        break
                a["members"] = [frag for m in a["members"] for frag in clip_member_outside(m, lo, hi)]
                b["members"] = [frag for m in b["members"] for frag in clip_member_outside(m, lo, hi)]
                a["truncated_right"] = True
                b["truncated_left"] = True
                changed = True
            if containment_removed:
                break
        if not changed:
            break
        rebuilt = []
        for c in chains:
            if not c["members"]:
                continue
            new_c = chain_from_members(c["tname"], c["strand"], c["members"],
                                        truncated_left=c["truncated_left"], truncated_right=c["truncated_right"])
            if not has_seed(new_c["members"], seed_mapq, seed_alnlen):
                continue
            if density(new_c) < min_density:
                continue
            rebuilt.append(new_c)
        chains = rebuilt
    return chains


# ------------------------------------------------------------------------ stage 4
def apply_density_cutoff(chains, min_density):
    return [c for c in chains if density(c) >= min_density]


def apply_ref_density_cutoff(chains, min_ref_density):
    """Applied after stage 3's consolidation, since merging across --merge-gap is
    itself a way to stretch the reference envelope."""
    if min_ref_density <= 0:
        return chains, []
    kept = [c for c in chains if ref_density(c) >= min_ref_density]
    dropped = [c for c in chains if ref_density(c) < min_ref_density]
    return kept, dropped


# ------------------------------------------------------------------------ stage 5
def _contains(outer, inner):
    return (outer["qstart"] <= inner["qstart"] and inner["qend"] <= outer["qend"] and
            outer["tstart"] <= inner["tstart"] and inner["tend"] <= outer["tend"])


def _containment_winner(a, b):
    """Stage 3 only. Density first (scale-free span alone lets a wide,
    sparse chain -- e.g. chrX/chrY homology noise -- swallow a much smaller,
    much denser real chain). Once both clear 50% density, prefer more
    matched bases over more span (span is what an inflated envelope, see
    stage 7's docstring, would win on)."""
    da, db = density(a), density(b)
    if da >= 0.5 and db >= 0.5:
        return a if a["matched"] >= b["matched"] else b
    return a if da >= db else b


def _majority_strand(members):
    """Matched-weighted, not span-weighted (span doesn't distinguish real
    backbone from unaligned filler)."""
    plus = sum(r["matched"] for r in members if r["strand"] == "+")
    minus = sum(r["matched"] for r in members if r["strand"] == "-")
    return "+" if plus >= minus else "-"


def consolidate_same_contig(chains, merge_gap):
    """Same-chromosome, opposite-strand chains always merge, never
    pick-a-winner (unlike stage 3): both already agree on the chromosome, so
    a strand disagreement is real assembly structure, not a conflict.
    "Adjacent" requires both query and chm13 gaps within merge_gap (chm13
    alone isn't enough -- two chains far apart on the contig can land close
    on chm13 by chance)."""
    chains = list(chains)
    changed = True
    while changed:
        changed = False
        by_tname = defaultdict(list)
        for c in chains:
            by_tname[c["tname"]].append(c)
        for tname, group in by_tname.items():
            if len(group) < 2:
                continue
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    a, b = group[i], group[j]
                    if a["strand"] == b["strand"]:
                        continue
                    tgap = max(a["tstart"], b["tstart"]) - min(a["tend"], b["tend"])
                    if a["tstart"] <= b["tend"] and b["tstart"] <= a["tend"]:
                        tgap = 0
                    qgap = max(a["qstart"], b["qstart"]) - min(a["qend"], b["qend"])
                    if a["qstart"] <= b["qend"] and b["qstart"] <= a["qend"]:
                        qgap = 0
                    if tgap > merge_gap or qgap > merge_gap:
                        continue
                    merged_members = a["members"] + b["members"]
                    merged = chain_from_members(tname, _majority_strand(merged_members), merged_members,
                                                 truncated_left=a["truncated_left"] or b["truncated_left"],
                                                 truncated_right=a["truncated_right"] or b["truncated_right"])
                    chains = [c for c in chains if c is not a and c is not b] + [merged]
                    changed = True
                    break
                if changed:
                    break
            if changed:
                break
    chains.sort(key=lambda c: c["qstart"])
    return chains


def run_pipeline_one_contig(rows, min_mapq, seed_mapq, seed_alnlen, min_density,
                            merge_gap, win_ratio, min_ref_density):
    chains = build_group_chains(rows, min_mapq, seed_mapq, seed_alnlen)
    chains = apply_density_cutoff(chains, min_density)
    chains = resolve_cross_contig_overlaps(chains, seed_mapq, seed_alnlen, min_density, win_ratio)
    chains = consolidate_same_contig(chains, merge_gap)
    chains, ref_dropped = apply_ref_density_cutoff(chains, min_ref_density)
    return chains, ref_dropped


# ------------------------------------------------------------------------ stage 7
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def load_censat(path):
    d = defaultdict(list)
    with _open_text(path) as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            try:
                d[a[0]].append((int(a[1]), int(a[2])))
            except ValueError:
                continue
    for c in d:
        d[c].sort()
    return d


def censat_overlap(chrom, s, e, censat):
    ints = censat.get(chrom)
    if not ints:
        return 0
    tot = 0
    for a, b in ints:
        if b <= s:
            continue
        if a >= e:
            break
        o = min(e, b) - max(s, a)
        if o > 0:
            tot += o
    return tot


def resolve_contig_chromosomes(chains, censat, split_min):
    """Stages 1-5 still report >=2 chromosomes for ~12% of contigs, nearly all
    acrocentric short-arm homology rather than chimerism. A chromosome counts as
    real when its chm13 span, minus the censat part for acrocentrics, clears
    split_min. None significant -> the largest raw span, so a pure-satellite
    contig is still placed; one -> that chromosome; two or more -> a real
    chimera, one row each."""
    by_chrom = defaultdict(list)
    for c in chains:
        by_chrom[c["tname"]].append(c)

    blocks = {}
    for chrom, group in by_chrom.items():
        members = [m for c in group for m in c["members"]]
        ref_span = sum(m["tend"] - m["tstart"] for m in members)
        cens = sum(censat_overlap(chrom, m["tstart"], m["tend"], censat) for m in members) if chrom in ACRO else 0
        blocks[chrom] = dict(members=members, ref_span=ref_span, eff=max(0, ref_span - cens))

    sig = [c for c, b in blocks.items() if b["eff"] >= split_min]

    def build_row(chrom):
        b = blocks[chrom]
        return chain_from_members(chrom, _majority_strand(b["members"]), b["members"])

    if not sig:
        return [build_row(max(blocks, key=lambda c: blocks[c]["ref_span"]))]
    if len(sig) == 1:
        return [build_row(sig[0])]
    return sorted((build_row(c) for c in sig), key=lambda r: r["qstart"])


# ------------------------------------------------------------------------ stage 9
def load_telo(path):
    """seqtk telo BED: contig, start, end, contig_len (4th column ignored --
    contig lengths come from the PAF, which is what the rest of the table is
    built from)."""
    d = defaultdict(list)
    with _open_text(path) as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            try:
                d[a[0]].append((int(a[1]), int(a[2])))
            except ValueError:
                continue
    return d


def extend_to_telomeres(rows, contig_len, chrom_len, telo, max_gap_frac=0.5):
    """A telomeric repeat at a contig end is evidence independent of the
    alignment that the end IS a chromosome terminus. Stages 1-7 stop at the last
    alignment block, so an acrocentric whose short arm is assembled but
    unalignable ends megabases short and reads as an assembly failure on the
    copy-number plot. Where a row reaches such an end, both sides snap to the
    terminus: contig to 0/contig_len, chm13 to 0/chrom_len.

    Both sides, or neither. copynumber_window.py takes its per-chromosome offset
    from ref_start (`diff = c_start / window_size`) while counting windows on the
    contig, so moving only the contig side would push the whole profile right by
    exactly the windows added. Moving ref_start to 0 as well puts the added
    windows where the reference says the short arm is.

    Strand decides which reference end a contig end maps to: on '-' the contig's
    start is the chromosome's q terminus. Only one row may claim each
    (chromosome, terminus) -- acrocentric short arms carry interstitial telomeric
    repeats, so a 1 Mb fragment can look like a terminus too; ties go to matched
    bases, which separates it from the chromosome-scale contig by three orders of
    magnitude. A row must also be the one that actually reaches that contig end,
    so a chimera's second arm cannot extend across the first, and the gap it
    would bridge must stay under max_gap_frac of its own contig span -- being the
    outermost row is not the same as being near the end: without that test a
    satellite contig whose only row covers a few percent of it extends over the
    whole contig. Legitimate extensions bridge a small fraction of the contig
    span, which is what max_gap_frac bounds.
    Returns (rows, applied, skipped)."""
    by_contig = defaultdict(list)
    for r in rows:
        by_contig[r[0]].append(r)

    def touches(contig, at_end):
        n = contig_len.get(contig)
        ivs = telo.get(contig, [])
        return any(e == n for s, e in ivs) if at_end else any(s == 0 for s, e in ivs)

    cand = defaultdict(list)
    for r in rows:
        contig, cs, ce, strand, chrom = r[0], r[1], r[2], r[3], r[4]
        n = contig_len.get(contig)
        if n is None or chrom not in chrom_len:
            continue
        group = by_contig[contig]
        limit = max_gap_frac * (ce - cs)
        if cs > 0 and cs <= limit and cs == min(x[1] for x in group) and touches(contig, False):
            cand[(chrom, "p" if strand == "+" else "q")].append((r[7], contig, id(r), r, "start"))
        if ce < n and (n - ce) <= limit and ce == max(x[2] for x in group) and touches(contig, True):
            cand[(chrom, "q" if strand == "+" else "p")].append((r[7], contig, id(r), r, "end"))

    applied, skipped = [], []
    for (chrom, side), v in sorted(cand.items()):
        v.sort(key=lambda t: (-t[0], t[1]))
        for _, contig, _, _, which in v[1:]:
            skipped.append((contig, chrom, side, which))
        _, contig, _, r, which = v[0]
        n = contig_len[contig]
        before = list(r)
        if which == "start":
            r[1] = 0
            if r[3] == "+":
                r[5] = 0
            else:
                r[6] = chrom_len[chrom]
        else:
            r[2] = n
            if r[3] == "+":
                r[6] = chrom_len[chrom]
            else:
                r[5] = 0
        applied.append((contig, chrom, side, which, before, list(r)))
    return rows, applied, skipped


# ------------------------------------------------------------------------ stage 8
def build_support_index(paf_rows, min_mapq):
    """(chrom, contig) -> mapq-filtered PAF records, for measuring how much
    real alignment a row has inside a disputed reference interval."""
    idx = defaultdict(list)
    for r in paf_rows:
        if r["mapq"] < min_mapq:
            continue
        idx[(r["tname"], r["qname"])].append(
            (r["tstart"], r["tend"], r["qstart"], r["qend"], r["matched"], r["strand"]))
    return idx


def row_support(idx, row, lo, hi):
    """Matched bases of the row's contig aligned to its chromosome inside
    [lo, hi), records clipped proportionally. Restricted to the row's contig
    span so another row of the same contig does not lend it support."""
    tot = 0.0
    for ts, te, qs, qe, m, _ in idx.get((row[4], row[0]), []):
        if qe <= row[1] or qs >= row[2]:
            continue
        ov = min(te, hi) - max(ts, lo)
        if ov > 0:
            tot += m * ov / (te - ts)
    return tot


def clip_row_refside(row, keep_lo, keep_hi, idx):
    """Clip a row's reference claim to [keep_lo, keep_hi). The clipped side is
    recomputed from the row's own supporting records, so an envelope stretched
    across satellite does not survive as an empty claim reaching exactly to
    the cut -- the first supporting record beyond it can sit several Mb further
    out than the cut itself. The untouched side keeps the original coordinates,
    which stage 9 may have snapped to a terminus. Returns a new row or None
    when nothing supported remains."""
    contig, cs, ce, strand, chrom, rs, re = row[0], row[1], row[2], row[3], row[4], row[5], row[6]
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
    return [contig, fcs, fce, strand, chrom, frs, fre, total]


def resolve_overlaps_between_contigs(rows, support_idx):
    """Two contigs cannot both hold the same stretch of chm13.

    Stages 1-7 decide one contig at a time, so nothing stops a pericentromeric
    fragment from claiming reference already covered by the chromosome-scale
    contig it sits inside. The copy-number window axis is built by walking the
    ref.table in order and accumulating reference gaps, so such a pair lands on
    the same x positions and counts one locus twice.

    Containment -> the contained row is dropped, keeping one row per contig
    stretch and the backbone contigs whole (adjudicating containment by local
    support instead splices satellite fragments into every backbone whose
    centromere is collapsed; sex-chromosome consolidation, where an inflated
    envelope really does swallow the true contig, is handled downstream by
    postprocess_sex_chrom.py). A PARTIAL overlap goes to the contig with more
    real alignment inside it (row_support), NOT to the wider envelope, and
    costs the loser only the disputed interval (dropping outright would cost
    a genuinely chimeric arm its row -- two adjacent contigs on the same
    chromosome can overlap by a tenth of their length without either being a
    duplicate of the other). The loser's new boundary is
    recomputed from its own supporting records so no empty envelope survives
    up to the cut. Ties fall back to span, then contig name, so the choice
    does not depend on input order. Returns (kept, dropped, log)."""
    by_chrom = defaultdict(list)
    for r in rows:
        by_chrom[r[4]].append(r)

    kept_all, dropped, log = [], [], []
    for chrom in sorted(by_chrom):
        # Largest first; contig name breaks ties so the outcome does not
        # depend on input order. A clipped row re-enters the queue and is
        # re-checked against everything already kept.
        queue = sorted(by_chrom[chrom], key=lambda r: (-(r[6] - r[5]), r[0]))
        kept = []
        while queue:
            r = queue.pop(0)
            alive = True
            i = 0
            while i < len(kept):
                k = kept[i]
                lo, hi = max(k[5], r[5]), min(k[6], r[6])
                if hi <= lo:
                    i += 1
                    continue
                if k[5] <= r[5] and r[6] <= k[6]:
                    log.append((chrom, lo, hi, k[0], None, r[0], None, None))
                    dropped.append(r)
                    alive = False
                    break
                if r[5] <= k[5] and k[6] <= r[6]:
                    log.append((chrom, lo, hi, r[0], None, k[0], None, None))
                    dropped.append(k)
                    kept.pop(i)
                    continue
                s_k = row_support(support_idx, k, lo, hi)
                s_r = row_support(support_idx, r, lo, hi)
                if (s_r, r[6] - r[5], r[0]) > (s_k, k[6] - k[5], k[0]):
                    winner, loser = r, k
                else:
                    winner, loser = k, r
                # Partial: the loser keeps only its non-overlapping side.
                if loser[5] < lo:
                    piece = clip_row_refside(loser, loser[5], lo, support_idx)
                else:
                    piece = clip_row_refside(loser, hi, loser[6], support_idx)
                log.append((chrom, lo, hi, winner[0],
                            s_r if winner is r else s_k,
                            loser[0], s_k if winner is r else s_r, piece))
                if loser is k:
                    if piece is None:
                        dropped.append(k)
                        kept.pop(i)
                        continue
                    kept[i] = piece
                    i += 1
                else:
                    if piece is None:
                        dropped.append(r)
                    else:
                        queue.append(piece)
                    alive = False
                    break
            if alive:
                kept.append(r)
        kept_all.extend(kept)
    return kept_all, dropped, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True,
                     help="PAF from minimap2 -cx asm5 <assembly> <chm13>, secondaries removed")
    ap.add_argument("--censat", default=None,
                     help="chm13 censat BED. Without it stage 7 cannot tell an acrocentric "
                          "short-arm match from a real second chromosome")
    ap.add_argument("--telo", default=None,
                     help="seqtk telo BED of the assembly, enabling stage 9")
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--seed-mapq", type=int, default=50)
    ap.add_argument("--seed-alnlen", type=int, default=10000)
    ap.add_argument("--min-density", type=float, default=0.10)
    ap.add_argument("--min-ref-density", type=float, default=0.02,
                     help="stage 4 cutoff on the chm13 side (0 = off)")
    ap.add_argument("--merge-gap", type=int, default=1_000_000)
    ap.add_argument("--split-min", type=int, default=1_000_000,
                     help="stage 7: censat-adjusted span a chromosome needs to count as real")
    ap.add_argument("--overlap-win-density-ratio", type=float, default=2.0,
                     dest="overlap_win_ratio",
                     help="stage 5: density ratio at which the denser chain wins outright "
                          "instead of both being clipped (0 = always clip)")
    ap.add_argument("--telo-max-gap-frac", type=float, default=0.5,
                     help="stage 9: largest gap to a contig end, as a fraction of the row's "
                          "own span, that may be bridged")
    args = ap.parse_args()

    rows = parse_paf_gt_as_query(args.input)
    print(f"[load] {args.input}: {len(rows)} PAF records", file=sys.stderr)

    by_contig = defaultdict(list)
    for r in rows:
        by_contig[r["qname"]].append(r)
    print(f"[load] {len(by_contig)} distinct GT contigs", file=sys.stderr)

    censat = load_censat(args.censat) if args.censat else defaultdict(list)

    out_rows = []
    for contig in sorted(by_contig):
        contig_rows = by_contig[contig]
        chains, ref_dropped = run_pipeline_one_contig(
            contig_rows, args.min_mapq, args.seed_mapq, args.seed_alnlen,
            args.min_density, args.merge_gap, args.overlap_win_ratio,
            args.min_ref_density)
        for c in ref_dropped:
            print(f"[stage4b] dropped {contig} on {c['tname']}: contig {c['qend']-c['qstart']:,} bp "
                  f"claims ref {c['tend']-c['tstart']:,} bp, ref-side density {ref_density(c):.5f} "
                  f"< {args.min_ref_density}", file=sys.stderr)
        if chains:
            chains = resolve_contig_chromosomes(chains, censat, args.split_min)
        for c in chains:
            out_rows.append([contig, c["qstart"], c["qend"], c["strand"],
                             c["tname"], c["tstart"], c["tend"], c["matched"]])

    if args.telo:
        # Lengths both come from the PAF, so stage 9 needs no extra input beyond
        # the BED: chm13 is the PAF query (chromosome length in field 2) and the
        # GT contig is the PAF target (contig length in field 7).
        contig_len = {r["qname"]: r["qlen"] for r in rows}
        chrom_len = {r["tname"]: r["tlen"] for r in rows}
        telo = load_telo(args.telo)
        print(f"[stage9] {args.telo}: telomere calls on {len(telo)} contigs", file=sys.stderr)
        out_rows, applied, skipped = extend_to_telomeres(out_rows, contig_len, chrom_len, telo,
                                                         args.telo_max_gap_frac)
        for contig, chrom, side, which, before, after in applied:
            d_c = (after[2] - after[1]) - (before[2] - before[1])
            d_r = (after[6] - after[5]) - (before[6] - before[5])
            print(f"[stage9] {contig} {chrom} {side}-terminus via contig {which}: "
                  f"contig {before[1]:,}-{before[2]:,} -> {after[1]:,}-{after[2]:,} ({d_c:+,} bp), "
                  f"ref {before[5]:,}-{before[6]:,} -> {after[5]:,}-{after[6]:,} ({d_r:+,} bp)",
                  file=sys.stderr)
        for contig, chrom, side, which in skipped:
            print(f"[stage9] {contig} also claims {chrom} {side}-terminus (contig {which}): "
                  f"skipped, fewer matched bases", file=sys.stderr)

    support_idx = build_support_index(rows, args.min_mapq)
    out_rows, dropped, ov_log = resolve_overlaps_between_contigs(out_rows, support_idx)
    for chrom, lo, hi, w, sw, l, sl, piece in ov_log:
        if sw is None:
            print(f"[stage8] {chrom} {lo:,}-{hi:,} ({(hi-lo)/1e6:.2f} Mb): "
                  f"{l} contained in {w} -> {l} dropped", file=sys.stderr)
            continue
        outcome = (f"{l} dropped" if piece is None else
                   f"{l} clipped to {piece[5]:,}-{piece[6]:,}")
        print(f"[stage8] {chrom} {lo:,}-{hi:,} ({(hi-lo)/1e6:.2f} Mb): "
              f"{w} support {sw/1e6:.3f} Mb beats {l} {sl/1e6:.3f} Mb -> {outcome}",
              file=sys.stderr)

    # Chromosome, then reference position -- NOT contig name. The copy-number
    # window axis accumulates the reference gap only while consecutive rows
    # share a chromosome and restarts otherwise, so two contigs of one
    # chromosome have to be adjacent or their windows both start at index 0 and
    # overlap. Plain string order on the name, matching make_reference_table.py.
    out_rows.sort(key=lambda r: (r[4], r[5]))

    for r in out_rows:
        print("\t".join(str(x) for x in r))
    print(f"[done] {len(out_rows)} ref.table rows for {len(by_contig)} contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
