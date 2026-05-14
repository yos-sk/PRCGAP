#!/usr/bin/env python3

import argparse
import sys
import csv

def filter_main(args):
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        print('\t'.join(dreader.fieldnames), file=hout)
        for F in dreader:
            if int(F["Checked_Read_Num_Control"]) < args.min_normal_depth:
                if F["Is_Filter"] == "PASS":
                    F["Is_Filter"] = "Low_normal_depth"
                else:
                    F["Is_Filter"] = F["Is_Filter"] + ';' + "Low_normal_depth"

            if int(F["Checked_Read_Num_Tumor"]) < args.min_tumor_depth:
                if F["Is_Filter"] == "PASS":
                    F["Is_Filter"] = "Low_tumor_depth"
                else:
                    F["Is_Filter"] = F["Is_Filter"] + ';' + "Low_tumor_depth"

            if int(F["Supporting_Read_Num_Tumor"]) / int(F["Checked_Read_Num_Tumor"]) < args.min_VAF:
                if F["Is_Filter"] == "PASS":
                    F["Is_Filter"] = "Low_VAF"
                else:
                    F["Is_Filter"] = F["Is_Filter"] + ';' + "Low_VAF"

            if F["SV_ID"][0] == 'd' and F["Dir_1"] == "-" and F["Dir_2"] == "+":
                if F["Is_Filter"] == "PASS":
                    F["Is_Filter"] = "False_dup"
                else:
                    F["Is_Filter"] = F["Is_Filter"] + ';' + "False_dup"

            print('\t'.join(F.values()), file = hout)
    

def arg_parser():
    parser = argparse.ArgumentParser(prog = "nanomonsv_filter",
            description = "Filter nanomonsv results")

    parser.add_argument("--input", "-i",  type = str,
                            help = "Path to the nanomonsv postprocess result file")

    parser.add_argument("--output", "-o", type = str,
            help = "Path to the output file")


    parser.add_argument("--min_VAF", default = 0.05, type = float,
            help = "Minimum required estimated variant allele frequency (default: 0.05)")

    parser.add_argument("--min_normal_depth", default = 10, type = int,
            help = "Minimum required normal depth  (default: 10)")

    parser.add_argument("--min_tumor_depth", default = 10, type = int,
            help = "Minimum required tumor depth  (default: 10)")

    parser.set_defaults(func = filter_main)


    return parser

def main():
    parser = arg_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args.func(args)

if __name__ == "__main__":
    main()
