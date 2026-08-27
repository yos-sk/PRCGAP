#! /usr/bin/env python3
"""Flag SV calls that are a copy-number change of the tandem repeat they sit in.

Being inside a tandem repeat is not itself a reason to drop a call -- that
removes true variants that merely happen to land in one. What should be dropped
is a change in the array's own copy number. So the repeat region only gates the
test, and the decision comes from the inserted sequence:

  - the call loses sequence (nanomonsv's d_ detection class, or no inserted
    sequence at all) -- a contraction of the array, nothing to align;
  - its inserted sequence is a tandem repeat in its own right
    (self_repeat >= --min-self-repeat);
  - its inserted sequence realigns to a breakpoint flank over most of its length
    (coverage >= --min-coverage) AND the flank stretch it matched is itself
    tandem repetitive (>= --min-flank-self-repeat) -- one more copy of the local
    array;
  - its inserted sequence is shorter than --min-length.

The second half of the flank criterion is what separates an array expansion from
a full-length mobile element that landed beside a homologous copy: the element is
just as similar to its flank, but what it matched is another element, not an
array. Requiring the matched flank to be repetitive needs no annotation.

Scores come from simple_repeat_score.py. Without --scores the script keeps its
older behaviour and flags on the region alone.

Thresholds and the alignment penalties behind them were set in
workspace/repeat_bench Experiment 10, which agreed with all 15 IGV calls checked.
"""

import sys, csv
import pysam


def filter_indel_in_simple_repeat(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, simple_repeat_tb, simple_repeat_dist_margin = 30):

    if tchr1 == tchr2 and tdir1 == '+' and tdir2 == '-':
        sv_size = tpos2 - tpos1 + len(tinseq) - 1

        tabix_error_flag = False
        try:
            records = simple_repeat_tb.fetch(tchr1, max(tpos1 - simple_repeat_dist_margin + 1, 0),
                tpos1 + simple_repeat_dist_margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos1 >= int(record[1]) - simple_repeat_dist_margin and \
                    int(tpos2) <= int(record[2]) + simple_repeat_dist_margin:
                    return True

        return False


def load_scores(path):
    """SV_ID -> score row from simple_repeat_score.py."""
    scores = {}
    if path is None:
        return scores
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            scores[row["SV_ID"]] = row
    return scores


def decide(sv_id, in_region, scores, args):
    """(flag, reason) for one call. flag=True means mark it Simple_repeat."""
    if not in_region:
        return False, "not_in_repeat_region"
    if not scores:                       # region-only mode
        return True, "in_repeat_region"
    row = scores.get(sv_id)
    if row is None:                      # no inserted sequence to test
        return True, "no_inserted_seq"
    if sv_id.startswith("d_"):
        return True, "deletion"
    srep = float(row["self_repeat"])
    if srep >= args.min_self_repeat:
        return True, f"self_repeat={srep:.3f}"
    cov = float(row["flank_coverage"])
    fsrep = float(row["flank_self_repeat"])
    if cov >= args.min_coverage and fsrep >= args.min_flank_self_repeat:
        return True, f"flank={cov:.3f} flank_self={fsrep:.3f}"
    if int(row["ins_len"]) < args.min_length:
        return True, f"short={row['ins_len']}"
    return False, f"flank={cov:.3f} flank_self={fsrep:.3f} self_repeat={srep:.3f}"


def post_filter_main(args):

    simple_repeat_tb = pysam.TabixFile(args.simple_repeat_bed) if args.simple_repeat_bed is not None else None
    scores = load_scores(args.scores)

    n_flagged = 0
    reasons = {}
    with open(args.sv_list_file, 'r') as hin, open(args.output_file, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        print('\t'.join(header), file = hout)

        for F in dreader:
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F["Chr_1"], int(F["Pos_1"]), F["Dir_1"], F["Chr_2"], int(F["Pos_2"]), F["Dir_2"], F["Inserted_Seq"]

            if args.simple_repeat_bed is not None:
                in_region = filter_indel_in_simple_repeat(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, simple_repeat_tb)
            else:
                in_region = False

            simple_repeat_flag, reason = decide(F.get("SV_ID", ""), in_region, scores, args)
            if simple_repeat_flag:
                n_flagged += 1
                key = reason.split("=")[0]
                reasons[key] = reasons.get(key, 0) + 1
                if F["Is_Filter"] == "PASS":
                    F["Is_Filter"] = "Simple_repeat"
                else:
                    F["Is_Filter"] = F["Is_Filter"] + ';' + "Simple_repeat"

            print('\t'.join(F.values()), file = hout)

    print("[add_simple_repeat] flagged {} calls: {}".format(
        n_flagged, ", ".join(f"{k}={v}" for k, v in sorted(reasons.items(),
                                                           key=lambda kv: -kv[1]))),
          file=sys.stderr)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog = "nanomonsv_simple_repeat_annot",
        description = "Add simple repeat annotation to the result of nanomonsv")

    parser.add_argument("sv_list_file", type = str,
                        help = "Path to the nanomonsv result file")

    parser.add_argument("output_file", type = str,
                        help = "Path to the output file")

    parser.add_argument("simple_repeat_bed", metavar = "simpleRepeat.bed.gz", type = str, default = None,
                        help = "Path to the tabix indexed simple repeat bed file")

    parser.add_argument("--scores", default = None,
                        help = "score table from simple_repeat_score.py; without "
                               "it the region alone decides")
    parser.add_argument("--min-coverage", type = float, default = 0.60)
    parser.add_argument("--min-flank-self-repeat", type = float, default = 0.48)
    parser.add_argument("--min-self-repeat", type = float, default = 0.70)
    parser.add_argument("--min-length", type = int, default = 50)

    args = parser.parse_args()

    post_filter_main(args)
