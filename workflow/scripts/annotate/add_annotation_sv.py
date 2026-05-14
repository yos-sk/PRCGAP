#! /usr/bin/env python3

import sys, csv
import pysam
import argparse

def annotation_ref(args) -> None:
    def open_list_file(sv_list_file: list) -> dict:
        out = dict()
        with open(sv_list_file, 'r') as f:
            for line in f:
                items = line.rstrip('\n').split('\t')
                if items[0] in out:
                    out[items[0]].append(items[1])
                else:
                    out[items[0]] = [items[1]]
        return out
    def list_annotation_feature(sv_id: str, sv_list: dict, nanomonsv_list) -> str:
        if sv_id in sv_list:
            return ",".join([nanomonsv_list[s_id] for s_id in sv_list[sv_id]])
        else:
            return ""
    
    def open_nanomonsv_file(nanomonsv_file: str) -> dict:
        out = dict()
        with open(nanomonsv_file, "r") as f:
            dreader = csv.DictReader(f, delimiter = '\t')
            for F in dreader:
                tsvid = F["SV_ID"]
                tfilter = str(F["Is_Filter"])

                out[tsvid] = tfilter
            return out
    
    feature = args.feature
    sv_list = open_list_file(args.list)
    nanomonsv_results = open_nanomonsv_file(args.nanomonsv)

    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        print('\t'.join(header), feature, feature + "_filter", sep="\t", file = hout)

        for F in dreader:
            tsvid = F["SV_ID"]

            annotation = list_annotation_feature(tsvid, sv_list, nanomonsv_results)
            if annotation != "":
                record = "True\t" + annotation
            else:
                record = "False\t" + "-"

            print('\t'.join(F.values()), record, sep="\t", file = hout)
    

def annotation_gene(args) -> None:
    def open_cgc_file(cgc_file: str) -> set:
        out = set()
        
        with open(cgc_file, 'r') as f:
            for i, line in enumerate(f):
                if i == 0: continue
                items = line.rstrip('\n').split('\t')
                out.add(items[0])
        
        return out
    
    def tbx_annotation(tchr1: str, tpos1: int, annotation_tbx, annotation_dist_margin = 1) -> str:
        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr1, max(tpos1 - annotation_dist_margin + 1, 0), 
                tpos1 + annotation_dist_margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos1 >= int(record[1]) - annotation_dist_margin and \
                    int(tpos1) <= int(record[2]) + annotation_dist_margin:
                    return record[3]

        return ""
    
    annotation_tbx = pysam.TabixFile(args.bed) if args.bed is not None else None
    cgc_gene_set = open_cgc_file(args.cgc) if args.cgc is not None else None
    
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        print('\t'.join(header), "bp1_gene", "bp2_gene", "bp1_cgc", "bp2_cgc", sep="\t", file = hout)
        for F in dreader:
            tchr1, tpos1, tchr2, tpos2 = F["Chr_1"], int(F["Pos_1"]), F["Chr_2"], int(F["Pos_2"])
            bp1_gene = tbx_annotation(tchr1, tpos1, annotation_tbx)
            if bp1_gene != "":
                record = bp1_gene
            else:
                record = "-"
            
            bp2_gene = tbx_annotation(tchr2, tpos2, annotation_tbx)
            if bp2_gene != "":
                record = record + "\t" + bp2_gene
            else:
                record = record + "\t" + "-"
            
            if bp1_gene in cgc_gene_set:
                record = record + "\tTrue"
            else:
                record = record + "\tFalse"
                
            if bp2_gene in cgc_gene_set:
                record = record + "\tTrue"
            else:
                record = record + "\tFalse"
        
            print('\t'.join(F.values()), record, sep="\t", file = hout)
        
def annotation_rmsk(args) -> None:
    def tbx_annotation_feature(tchr1: str, tpos1: int, annotation_tbx, annotation_dist_margin = 1) -> list:
        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr1, max(tpos1 - annotation_dist_margin + 1, 0),
                tpos1 + annotation_dist_margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos1 >= int(record[1]) - annotation_dist_margin and \
                    int(tpos1) <= int(record[2]) + annotation_dist_margin:
                    return [record[-2], record[-1]]

        return []
    
    annotation_tbx = pysam.TabixFile(args.bed) if args.bed is not None else None

    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "bp1_rmsk_class", "bp1_rmsk_type", "bp2_rmsk_class", "bp2_rmsk_type", sep="\t", file = hout)

        for F in dreader:
            tchr1, tpos1, tchr2, tpos2= F["Chr_1"], int(F["Pos_1"]), F["Chr_2"], int(F["Pos_2"])

            annotation_bp1 = tbx_annotation_feature(tchr1, tpos1, annotation_tbx)
            annotation_bp2 = tbx_annotation_feature(tchr2, tpos2, annotation_tbx)

            if annotation_bp1 != []:
                record_bp1 = annotation_bp1[1] + "\t" + annotation_bp1[0]
            else:
                record_bp1 = "-\t-"

            if annotation_bp2 != []:
                record_bp2 = annotation_bp2[1] + "\t" + annotation_bp2[0]
            else:
                record_bp2 = "-\t-"

            print('\t'.join(F.values()), record_bp1, record_bp2, sep="\t", file = hout)

def annotation_size(args) -> None:
    def extract_contig_size(input: str, out: dict) -> None:
        with open(input, "r") as f:
            seq = ""
            rname = ""
            for line in f:
                if line[0] == ">":
                    if rname == "":
                        rname = line.rstrip("\n")[1:]
                    else:
                        out[rname] = len(seq)
                        seq = ""
                        rname = line.rstrip("\n")[1:]
                else:
                    seq += line.rstrip("\n")
            
            out[rname] = len(seq)
    
    contig_size = dict()
    extract_contig_size(args.hap1_fasta, contig_size)
    extract_contig_size(args.hap2_fasta, contig_size)
    
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        print('\t'.join(header), "bp1_contig_size", "bp2_contig_size", sep="\t", file = hout)
        for F in dreader:
            tchr1, tchr2= F["Chr_1"], F["Chr_2"]
            record = str(contig_size[tchr1]) + "\t" + str(contig_size[tchr2])
            
            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_conv(args) -> None:
    
    def open_lift_over_file(lift_over_file: str) -> dict:
        out = dict()
        
        with open(lift_over_file, 'r') as f:
            for line in f:
                items = line.rstrip('\n').split('\t')
                key = items[0] + ":" + items[1]
                sv_id = items[-1]
                category = items[-2]
                if category == "Match":
                    value = items[2] + ":" + items[3]
                elif category == "Gap":
                    value = items[2] + ":" + items[3] + "-" + items[4] 
                else:
                    value = "-"
                if sv_id in out:
                    if key in out[sv_id]:
                        if category in out[sv_id][key]:
                            out[sv_id][key][category].append(value)
                        else:
                            out[sv_id][key][category] = [value]
                    else:
                        out[sv_id][key] = {category: [value]}
                else:
                    out[sv_id] = {key: {category: [value]}}
        
        return out
    
    liftover_dict = open_lift_over_file(args.liftover)
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        if args.feature == "liftover_GRCh38":
            print('\t'.join(header), "bp1_GRCh38_lf", "bp2_GRCh38_lf", "bp1_GRCh38_lf_cood",  "bp2_GRCh38_lf_cood", sep="\t", file = hout)
        elif args.feature == "liftover_chm13":
            print('\t'.join(header), "bp1_chm13_lf", "bp2_chm13_lf", "bp1_chm13_lf_cood",  "bp2_chm13_lf_cood", sep="\t", file = hout)
        else: 
            print('\t'.join(header), "bp1_prcgap_lf", "bp2_prcgap_lf", "bp1_prcgap_lf_cood",  "bp2_prcgap_lf_cood", sep="\t", file = hout)
        
        for F in dreader:
            tchr1, tpos1, tchr2, tpos2, tsvid = F["Chr_1"], int(F["Pos_1"]),  F["Chr_2"], int(F["Pos_2"]), F["SV_ID"]
            
            if "Match" in liftover_dict[tsvid][tchr1 + ":" + str(tpos1)]:
                record = "True"
            else:
                record = "False"
            
            if "Match" in liftover_dict[tsvid][tchr2 + ":" + str(tpos2)]:
                record =  record + "\tTrue"
            else:
                record = record + "\tFalse"
            
            if "Match" in liftover_dict[tsvid][tchr1 + ":" + str(tpos1)]:
                record = record + "\t" +  ",".join([str(i) for i in liftover_dict[tsvid][tchr1 + ":" + str(tpos1)]["Match"]])
            elif "Gap" in liftover_dict[tsvid][tchr1 + ":" + str(tpos1)]:
                record = record + "\t" + ",".join([str(i) for i in liftover_dict[tsvid][tchr1 + ":" + str(tpos1)]["Gap"]])
            else:
                record = record + "\t-"
                
            if "Match" in liftover_dict[tsvid][tchr2 + ":" + str(tpos2)]:
                record = record + "\t" + ",".join([str(i) for i in liftover_dict[tsvid][tchr2 + ":" + str(tpos2)]["Match"]])
            elif "Gap" in liftover_dict[tsvid][tchr2 + ":" + str(tpos2)]:
                record = record + "\t" + ",".join([str(i) for i in liftover_dict[tsvid][tchr2 + ":" + str(tpos2)]["Gap"]])
            else:
                record = record + "\t-"
            
            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_kmer(args):
    def open_support_reads(input_file: str) -> dict:
        out = dict()
        with open(input_file, 'r') as f:
            for line in f:
                items = line.rstrip('\n').split('\t')
                if items[7] in out:
                    out[items[7]].add(items[8])
                else:
                    out[items[7]] = {items[8]}

        return out

    def open_kmer_list(input_file: str) -> dict:
        out = dict()
        with open(input_file, 'r') as f:
            for line in f:
                items = line.rstrip('\n').split('\t')
                out[items[0]] = float(items[3])

        return out
    
    #support_read_hifi_dict = open_support_reads(args.support_reads_hifi)
    support_read_dict = open_support_reads(args.support_reads)
    #kmer_list_hifi = open_kmer_list(args.list_hifi)
    kmer_list = open_kmer_list(args.list)
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        #print('\t'.join(header), "kmer_ratio_HiFi", "kmer_ratio_ONT", sep="\t", file=hout)
        print('\t'.join(header), "kmer_ratio", sep="\t", file=hout)
        for F in dreader:
            sv_ids = F["Identical_SVs"].split(",")
            #t_kmer_ratio_hifi = list()
            t_kmer_ratio = list()
            for sv_id in sv_ids:
                #if sv_id[:2] == "1_":
                #    for read_id in support_read_hifi_dict[sv_id[2:]]:
                #        t_kmer_ratio_hifi.append(kmer_list_hifi[read_id])
                #else:
                for read_id in support_read_dict[sv_id]:
                    t_kmer_ratio.append(kmer_list[read_id])                    

            #record_hifi = f'{sum(t_kmer_ratio_hifi) / len(t_kmer_ratio_hifi):.3f}' if  len(t_kmer_ratio_hifi) != 0 else "-"
            record = f'{sum(t_kmer_ratio) / len(t_kmer_ratio):.3f}' if  len(t_kmer_ratio) != 0 else "-"
            #print('\t'.join(F.values()), record_hifi, record_ont, sep="\t", file = hout)
            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_misassembly(args) -> None:
    def tbx_annotation_misassembly(tchr1: str, tpos1: int, annotation_tbx, annotation_dist_margin = 1) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr1, max(tpos1 - annotation_dist_margin + 1, 0),
                tpos1 + annotation_dist_margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos1 >= int(record[1]) - annotation_dist_margin and \
                    int(tpos1) <= int(record[2]) + annotation_dist_margin:
                    return True

        return False
    
    annotation_tbx_hap1 = pysam.TabixFile(args.hap1_bed) if args.hap1_bed is not None else None
    annotation_tbx_hap2 = pysam.TabixFile(args.hap2_bed) if args.hap2_bed is not None else None

    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "bp1_misassembly", "bp2_misassembly", sep="\t", file = hout)

        for F in dreader:
            tchr1, tpos1, tchr2, tpos2 = F["Chr_1"], int(F["Pos_1"]), F["Chr_2"], int(F["Pos_2"])

            annotation_flag = tbx_annotation_misassembly(tchr1, tpos1, annotation_tbx_hap1) or tbx_annotation_misassembly(tchr1, tpos1, annotation_tbx_hap2)
            if annotation_flag:
                record = "True"
            else:
                record = "False"

            annotation_flag = tbx_annotation_misassembly(tchr2, tpos2, annotation_tbx_hap1) or tbx_annotation_misassembly(tchr2, tpos2, annotation_tbx_hap2)
            if annotation_flag:
                record = record + "\tTrue"
            else:
                record = record + "\tFalse"

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_centromere(args) -> None:
    def tbx_annotation_centromere(tchr: str, tpos: int, annotation_tbx, annotation_dist_margin = 1) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, max(tpos - annotation_dist_margin + 1, 0),
                tpos + annotation_dist_margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos >= int(record[1]) - annotation_dist_margin and \
                    int(tpos) <= int(record[2]) + annotation_dist_margin:
                    return record[3]

        return "-"
    
    annotation_tbx = pysam.TabixFile(args.centromere) if args.centromere is not None else None
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "bp1_centromere", "bp2_centromere", sep="\t", file = hout)

        for F in dreader:
            tchr1, tpos1, tchr2, tpos2 = F["Chr_1"], int(F["Pos_1"]), F["Chr_2"], int(F["Pos_2"])

            annotation_record_bp1 = tbx_annotation_centromere(tchr1, tpos1, annotation_tbx) 
            record = annotation_record_bp1 

            annotation_record_bp2 = tbx_annotation_centromere(tchr2, tpos2, annotation_tbx)
            record = record + "\t" + annotation_record_bp2

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_segdup(args) -> None:
    def tbx_annotation_segdup(tchr: str, tpos: int, annotation_tbx, annotation_dist_margin = 1) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, max(tpos - annotation_dist_margin + 1, 0),
                tpos + annotation_dist_margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos >= int(record[1]) - annotation_dist_margin and \
                    int(tpos) <= int(record[2]) + annotation_dist_margin:
                    return True

        return False
    
    annotation_tbx = pysam.TabixFile(args.segdup) if args.segdup is not None else None
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "bp1_segdup", "bp2_segdup", sep="\t", file = hout)

        for F in dreader:
            tchr1, tpos1, tchr2, tpos2 = F["Chr_1"], int(F["Pos_1"]), F["Chr_2"], int(F["Pos_2"])

            annotation_flag = tbx_annotation_segdup(tchr1, tpos1, annotation_tbx) 
            if annotation_flag:
                record = "True"
            else:
                record = "False"

            annotation_flag = tbx_annotation_segdup(tchr2, tpos2, annotation_tbx)
            if annotation_flag:
                record = record + "\tTrue"
            else:
                record = record + "\tFalse"

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_nanomonsv_other(args) -> None:
    class SVinfo:
        def __init__(self, chrom1, pos1, dir1, chrom2, pos2, dir2, sv_id):
            self.chrom1 = chrom1
            self.pos1 = pos1
            self.dir1 = dir1
            self.chrom2 = chrom2
            self.pos2 = pos2
            self.dir2 = dir2
            self.sv_id = sv_id

    def parse_nanomonsv_other(input_file) -> list:
        out = list()
        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                info = SVinfo(row["Chr_1"], int(row["Pos_1"]), row["Dir_1"], row["Chr_2"], int(row["Pos_2"]), row["Dir_2"], row["SV_ID"])
                out.append(info)

        return out
    
    info_other = parse_nanomonsv_other(args.other)
    margin = 100
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "nanomonsv_other", sep="\t", file = hout)

        for F in dreader:
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2 = F["Chr_1"], int(F["Pos_1"]), F["Dir_1"], F["Chr_2"], int(F["Pos_2"]), F["Dir_2"]
            record = "-"
            for h_i in info_other:
                if tchr1 != h_i.chrom1: continue
                if tchr2 != h_i.chrom2: continue
                if abs(tpos1 - h_i.pos1) > margin: continue
                if abs(tpos2 - h_i.pos2) > margin: continue
                if tdir1 != h_i.dir1: continue
                if tdir2 != h_i.dir2: continue
                record = h_i.sv_id
                break

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_gnomad(args) -> None:    
    def tbx_annotation_gnomad(tchr1: str, tpos1: int, tchr2: str, tpos2: int, annotation_tbx, margin=50) -> str:
        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr1, min(0, tpos1 - margin), tpos1 + margin)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if record[13] == "NA": continue
                if tchr1 == tchr2:
                    if tpos1 >= int(record[2]) - margin and tpos1 <= int(record[2]) + margin and \
                        tpos2 >= int(record[13]) - margin and tpos2 <= int(record[13]) + margin:
                        return record[3]
                else:
                    if tpos1 >= int(record[1]) - margin and tpos1 <= int(record[2]) + margin and \
                        tchr2 == record[9] and tpos2 >= int(record[13]) - margin and tpos2 <= int(record[13]) + margin:
                        return record[3]
        return "-"
    
    annotation_tbx = pysam.TabixFile(args.gnomad) if args.gnomad is not None else None
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "gnomAD", sep="\t", file = hout)

        for F in dreader:
            if F["bp1_GRCh38_lf"] == "False" or F["bp2_GRCh38_lf"] == "False":
                record = "-"
                print('\t'.join(F.values()), record, sep="\t", file = hout)
                continue
            p1 = F["bp1_GRCh38_lf_cood"].split(",")[0]
            p2 = F["bp2_GRCh38_lf_cood"].split(",")[0]
            tc1, tp1,  = p1.split(":")[0], int(p1.split(":")[1])
            tc2, tp2 = p2.split(":")[0], int(p2.split(":")[1])
            if tc1 not in ["chrX", "chrY"] and tc2 not in ["chrX", "chrY"] and int(tc1[3:]) > int(tc2[3:]):
                tchr1, tpos1 = tc2, tp2
                tchr2, tpos2 = tc1, tp1
            elif tc1 in ["chrX", "chrY"] or tc2 in ["chrX", "chrY"] and tc1 > tc2:
                tchr1, tpos1 = tc2, tp2
                tchr2, tpos2 = tc1, tp1
            else:
                tchr1, tpos1 = tc1, tp1
                tchr2, tpos2 = tc2, tp2
                
            record = tbx_annotation_gnomad(tchr1, tpos1, tchr2, tpos2, annotation_tbx) 

            print('\t'.join(F.values()), record, sep="\t", file = hout)
    
    
def arg_parser():
    # 1. GRch38/CHM13-based nanomonsv results annotation 
    # 2. Gene 
    # 3. RepeatMasker
    # 4. Contig size 
    # 5. Coordinate convertibility 
    # 6. k-mer ratio 
    # 7. misassembly 
    # 8. centromere
    # 9. segmental duplication
    # 10. nanomonsv other
    # 11. gnomAD
    
    parser = argparse.ArgumentParser(prog = "nanomonsv_annotation",
        description = "Add annotations to the result of nanomonsv")
    
    subparsers = parser.add_subparsers()
    # GRCh38/CHM13-based nanomonsv results annotation
    ref = subparsers.add_parser('ref', help="add annotation of GRCh38/CHM13-based results")
    ref.add_argument("--input", "-i", type = str,
                    help = "Path to the annotated nanomonsv result file")
    
    ref.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    ref.add_argument("--list", "-l", type = str,
                        help = "Path to the matched sv_id list file")

    ref.add_argument("--nanomonsv", "-n", type = str,
                        help = "Path to GRCh38/CHM13-based nanomonsv results")

    ref.add_argument("--feature", "-f", type = str, default = None,
                        help = "GRCh38_HiFi, GRCh38_ONT, CHM13_HiFi, or CHM13_ONT")
    
    ref.set_defaults(func=annotation_ref)
    
    # Gene annotation
    gene = subparsers.add_parser('gene', help="add annotation of genes")
    gene.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated nanomonsv result file")

    gene.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    gene.add_argument("--bed", "-b", type = str, default = None,
                        help = "Path to the tabix indexed liftoff bed file")
    
    gene.add_argument("--cgc", "-c", type = str, default = None,
                        help = "Path to the cancer gene census file")
    
    gene.set_defaults(func=annotation_gene)
    
    # RepeatMasker annotation
    rmsk = subparsers.add_parser('rmsk', help="add annotation of RepeatMasker result")
    rmsk.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated nanomonsv result file")

    rmsk.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    rmsk.add_argument("--bed", "-b", type = str, default = None,
                        help = "Path to the tabix indexed RepetMasker bed file")
    
    rmsk.set_defaults(func=annotation_rmsk)
    
    # Contig size annotation
    size = subparsers.add_parser('size', help="add annotation of contig size")
    size.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated nanomonsv result file")

    size.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    size.add_argument("--hap1_fasta", "-f", type = str, default = None,
                        help = "Path to the assembly hap1 fasta file")
    
    size.add_argument("--hap2_fasta", "-g", type = str, default = None,
                        help = "Path to the assembly hap2 fasta file")
    
    size.set_defaults(func=annotation_size)
    
    # Coordinate convertibility to GRCh38/CHM13 annotation
    conv = subparsers.add_parser('conv', help="add annotation of cordinate convertibility to GRCh38/CHM13")
    conv.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated nanomonsv result file")

    conv.add_argument("--output", "-o", type = str,
                        help = "Path to the output file") 
    
    conv.add_argument("--liftover", "-l", type = str, default = None,
                        help = "Path to the liftovered results")
    
    conv.add_argument("--feature", "-f", type = str, default = None,
                        help = "liftover_GRCh38, liftover_chm13 or liftover_prcgap")
    conv.set_defaults(func=annotation_conv)
    
    kmer = subparsers.add_parser('kmer', help="add annotation of mean unique kmer ratio which SV supporting reads have")
    kmer.add_argument("--input", "-i", type = str,
                        help = "Path to the nanomonsv result file")

    kmer.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    #kmer.add_argument("--support_reads_hifi", "-s", type=str,
    #                    help = "Path to the nanomonsv support reads (HiFi)")
    
    kmer.add_argument("--support_reads", "-t", type=str,
                        help = "Path to the nanomonsv support reads (ONT)")

    #kmer.add_argument("--list_hifi", "-l", type = str, default = None,
    #                    help = "Path to the k-mer ratio list file (HiFi)")

    kmer.add_argument("--list", "-m", type = str, default = None,
                        help = "Path to the k-mer ratio list file")

    kmer.set_defaults(func=annotation_kmer)
    
    misassembly = subparsers.add_parser('misassembly', help="add annotation whether SVs are in misassembly regions or not")
    misassembly.add_argument("--input", "-i", type = str,
                        help = "Path to the nanomonsv result file")

    misassembly.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    misassembly.add_argument("--hap1_bed", "-b", type = str, default = None,
                        help = "Path to the tabix indexed hap1 miassembly bed file")

    misassembly.add_argument("--hap2_bed", "-c", type = str, default = None,
                        help = "Path to the tabix indexed hap2 miassembly bed file")
    
    misassembly.set_defaults(func=annotation_misassembly)
    
    # Centromere annotation
    cen = subparsers.add_parser('cen', help="add annotation of centromere regions")
    cen.add_argument("--input", "-i", type = str,
                        help = "Path to the nanomonsv result file")
    
    cen.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    cen.add_argument("--centromere", "-s", type = str,
                        help = "Path to the centromere tabixed bed file")
    
    cen.set_defaults(func=annotation_centromere)
    
    # Segmental duplication annotation
    segdup = subparsers.add_parser('segdup', help="add annotation of segmental duplications")
    segdup.add_argument("--input", "-i", type = str,
                        help = "Path to the nanomonsv result file")
    
    segdup.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    segdup.add_argument("--segdup", "-s", type = str,
                        help = "Path to the segmental duplication tabixed bed file")
    
    segdup.set_defaults(func=annotation_segdup)

    # nanomonsv other annotation
    other = subparsers.add_parser('other', help="add annotation of nanomonsv other")
    other.add_argument("--input", "-i", type = str,
                        help = "Path to the nanomonsv result file")
    
    other.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    other.add_argument("--other", "-j", type = str,
                        help = "Path to the nanomonsv other result file")
    
    other.set_defaults(func=annotation_nanomonsv_other)
    
    # gnomAD annotation
    gnomad = subparsers.add_parser('gnomad', help="add annotation of gnomAD SVs")
    gnomad.add_argument("--input", "-i", type = str,
                        help = "Path to the nanomonsv result file")
    gnomad.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    gnomad.add_argument("--gnomad", "-g", type = str,
                        help = "Path to the gnomAD SVs tabixed bed file")
    gnomad.set_defaults(func=annotation_gnomad)

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
