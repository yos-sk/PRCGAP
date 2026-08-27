#! /usr/bin/env python3

import sys, csv
import pysam
import argparse
import gzip

def annotation_gene(args) -> None:
    def open_cgc_file(cgc_file: str) -> set:
        out = set()
        
        with open(cgc_file, 'r') as f:
            for i, line in enumerate(f):
                if i == 0: continue
                items = line.rstrip('\n').split('\t')
                out.add(items[0])
        
        return out

    def open_cmrg_file(cmrg_file: str) -> set:
        out = set()
        
        with open(cmrg_file, 'r') as f:
            for i, line in enumerate(f):
                if i == 0: continue
                items = line.rstrip('\n').split('\t')
                out.add(items[0])
        
        return out

    def get_transcript_ids(mane_summary_file: str) -> dict:
        """symbol -> {gene_id, [transcript_id]} from the MANE summary.

        Columns are located by name off the '#'-prefixed header rather than by
        position. Ensembl_Gene / Ensembl_nuc carry their version suffix and are
        kept as-is: a MANE release is tied to one Ensembl release, so a version
        that no longer matches the liftoff GTF means the two are out of step and
        should be updated together, not silently reconciled here.
        """
        out = dict()
        opener = gzip.open if mane_summary_file.endswith(".gz") else open
        with opener(mane_summary_file, "rt") as f:
            header = None
            for line in f:
                items = line.rstrip("\n").split("\t")
                if header is None:
                    header = {name.lstrip("#"): i for i, name in enumerate(items)}
                    for required in ("symbol", "Ensembl_Gene", "Ensembl_nuc"):
                        if required not in header:
                            raise KeyError(
                                f"{mane_summary_file}: MANE summary is missing "
                                f"column '{required}'")
                    continue
                gene_name = items[header["symbol"]]
                gene_id = items[header["Ensembl_Gene"]]
                transcript_id = items[header["Ensembl_nuc"]]
                if gene_name in out:
                    out[gene_name]["transcript_id"].append(transcript_id)
                else:
                    out[gene_name] = {"gene_id": gene_id, "transcript_id": [transcript_id]}
        return out
    
    def parse_gff_format(gff_column: str) -> dict:
        out = dict()
        items = gff_column.rstrip().split(";")
        for item in items:
            if len(item) == 0: continue
            key = item.split()[0]
            value = item.split()[1].strip("\"")
            out[key] = value

        return out
    
    version_mismatch = set()

    def tbx_annotation(tchr: str, tstart: int, tend: int, annotation_tbx, selected_transcript_id: dict) -> str:
        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend) 
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        output = dict()
        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart >= int(record[3]) and tend <= int(record[4]):
                    parsed_gff_info = parse_gff_format(record[8])
                    if record[2] == "gene": continue
                    if parsed_gff_info["gene_biotype"] != "protein_coding": continue
                    gene_name = parsed_gff_info["gene_name"] if "gene_name" in parsed_gff_info else parsed_gff_info["gene_id"]
                    gene_id = parsed_gff_info["gene_id"] + "." + parsed_gff_info["gene_version"]
                    transcript_id = parsed_gff_info["transcript_id"] + "." + parsed_gff_info["transcript_version"]
                    # Restrict to the MANE transcripts when the summary was
                    # given; without it every protein-coding transcript over the
                    # site is reported. A gene absent from MANE has no selected
                    # transcript and is dropped: reporting all of its transcripts
                    # instead would make the column hard to read.
                    if selected_transcript_id is not None:
                        selected = selected_transcript_id.get(gene_name)
                        if selected is None: continue
                        if (transcript_id not in selected["transcript_id"]
                                or gene_id != selected["gene_id"]):
                            # Same transcript at a different version means the
                            # MANE release and the liftoff GTF's Ensembl release
                            # are out of step; counted and reported once at the
                            # end so it cannot pass as "no gene here".
                            bare = transcript_id.split(".")[0]
                            if any(t.split(".")[0] == bare
                                   for t in selected["transcript_id"]):
                                version_mismatch.add(transcript_id)
                            continue
                    if gene_name in output:
                        if transcript_id in output[gene_name]:
                            if record[2] in ["five_prime_utr", "start_codon", "CDS", "stop_codon", "three_prime_utr"]:
                                output[gene_name][transcript_id] = record[2]
                            else:
                                output[gene_name][transcript_id] = record[2]
                        else:
                            output[gene_name][transcript_id] = record[2]
                    else:
                        output[gene_name] = {transcript_id: record[2]}
                                
        return output
    
    annotation_tbx = pysam.TabixFile(args.gff) if args.gff is not None else None
    cgc_gene_set = open_cgc_file(args.cgc) if args.cgc is not None else None
    cmrg_gene_set = open_cmrg_file(args.cmrg) if args.cmrg is not None else None

    selected_transcript_id = get_transcript_ids(args.mane) if args.mane is not None else None

    
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        print('\t'.join(header), "liftoff_gene", "liftoff_gene_info", "cgc", "cmrg",  sep="\t", file = hout)
        for F in dreader:
            tchr, tstart, tend = F["Contig"], int(F["Start"]), int(F["End"])
            annot_gene = tbx_annotation(tchr, tstart, tend, annotation_tbx, selected_transcript_id)
            if len(annot_gene) != 0:
                gene_names = ";".join([key for key in annot_gene])
                info = ""
                for i, gene_name in enumerate(annot_gene):
                    if i != 0:
                        info += ";"
                    for j, transcript_id in enumerate(annot_gene[gene_name]):
                        if j != 0:
                            info += "," + transcript_id + ":" + annot_gene[gene_name][transcript_id]
                        else: 
                            info += transcript_id + ":" + annot_gene[gene_name][transcript_id]
                record = gene_names + "\t" + info 
            else:
                record = "-\t-"
            
            if cgc_gene_set is None:
                # Optional database: not evaluated, not a negative.
                record = record + "\t-"
            else:
                output = set()
                for i, key in enumerate(annot_gene):
                    if key in cgc_gene_set:
                        output.add(key)
                if len(output) != 0:
                    record = record + "\t" + ",".join([k for k in output])
                else:
                    record = record + "\tFalse" 

            if cmrg_gene_set is None:
                # Optional database: not evaluated, not a negative.
                record = record + "\t-"
            else:
                output = set()
                for i, key in enumerate(annot_gene):
                    if key in cmrg_gene_set:
                        output.add(key)
                if len(output) != 0:
                    record = record + "\t" + ",".join([k for k in output])
                else:
                    record = record + "\tFalse" 

            print('\t'.join(F.values()), record, sep="\t", file = hout)

    if version_mismatch:
        print(f"WARNING: {len(version_mismatch)} transcripts matched a MANE "
              f"entry by ID but not by version, e.g. "
              f"{', '.join(sorted(version_mismatch)[:3])}. The MANE summary and "
              f"the Ensembl GTF used for liftoff are different releases; update "
              f"both together. Those genes carry no gene annotation.",
              file = sys.stderr)


def annotation_rmsk(args) -> None:
    def tbx_annotation_feature(tchr: str, tstart: int, tend: int, annotation_tbx) -> list:
        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart >= int(record[1]) and tend <= int(record[2]):
                    return [record[-2], record[-1]]

        return []
     
    annotation_tbx = pysam.TabixFile(args.bed) if args.bed is not None else None

    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "rmsk_class", "rmsk_type", sep="\t", file = hout)

        for F in dreader:
            tchr, tstart, tend = F["Contig"], int(F["Start"]), int(F["End"])

            annot_rmsk = tbx_annotation_feature(tchr, tstart, tend, annotation_tbx)

            if annot_rmsk != []:
                record = annot_rmsk[1] + "\t" + annot_rmsk[0]
            else:
                record = "-\t-"

            print('\t'.join(F.values()), record, sep="\t", file = hout)

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
        print('\t'.join(header), "contig_size", sep="\t", file = hout)
        for F in dreader:
            tchr = F["Contig"]
            record = str(contig_size[tchr]) 
            
            print('\t'.join(F.values()), record, sep="\t", file = hout)


def annotation_misassembly(args) -> None:
    def tbx_annotation_misassembly(tchr: str, tstart: int, tend: int, annotation_tbx) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart >= int(record[1]) and tend <= int(record[2]):
                    return True

        return False
    
    annotation_tbx_hap1 = pysam.TabixFile(args.hap1_bed) if args.hap1_bed is not None else None
    annotation_tbx_hap2 = pysam.TabixFile(args.hap2_bed) if args.hap2_bed is not None else None

    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "misassembly", sep="\t", file = hout)

        for F in dreader:
            tchr, tstart, tend = F["Contig"], int(F["Start"]), int(F["End"])

            annotation_flag = tbx_annotation_misassembly(tchr, tstart, tend, annotation_tbx_hap1) or tbx_annotation_misassembly(tchr, tstart, tend, annotation_tbx_hap2)
            if annotation_flag:
                record = "True"
            else:
                record = "False"

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_centromere(args) -> None:
     def tbx_annotation_centromere(tchr: str, tstart: int, tend: int, annotation_tbx) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart >= int(record[1]) and tend <= int(record[2]):
                    return record[3]

        return "-"
    
     annotation_tbx = pysam.TabixFile(args.centromere) if args.centromere is not None else None
     with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "centromere", sep="\t", file = hout)

        for F in dreader:
            tchr, tstart, tend = F["Contig"], int(F["Start"]), int(F["End"])

            record = tbx_annotation_centromere(tchr, tstart, tend, annotation_tbx) 

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_segdup(args) -> None:
     def tbx_annotation_segdup(tchr: str, tstart: int, tend: int, annotation_tbx) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart >= int(record[1]) and tend <= int(record[2]):
                    return record[7]

        return "-"
    
     annotation_tbx = pysam.TabixFile(args.segdup) if args.segdup is not None else None
     with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "segdup", "segdup_similarity", sep="\t", file = hout)

        for F in dreader:
            tchr, tstart, tend = F["Contig"], int(F["Start"]), int(F["End"])

            annot_record = tbx_annotation_segdup(tchr, tstart, tend, annotation_tbx) 
            if annot_record != "-":
                record = "True\t" + annot_record
            else:
                record = "False\t" + annot_record 

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_point_mutation_other(args) -> None:
    def tbx_annotation_other(tchr: str, tstart: int, tend: int, tref: str, talt: str, annotation_tbx) -> bool:

        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart == int(record[1]) - 1 and tref == record[3] and talt == record[4]:
                    return True

        return False

    annotation_tbx = pysam.TabixFile(args.other) if args.other is not None else None
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "point_mutation_other", sep="\t", file = hout)

        for F in dreader:
            tchr, tstart, tend, tref, talt = F["Contig"], int(F["Start"]), int(F["End"]), F["Ref"], F["Alt"]
            annotation_flag = tbx_annotation_other(tchr, tstart, tend, tref, talt, annotation_tbx) 
            if annotation_flag:
                record = "True"
            else:
                record = "False"

            print('\t'.join(F.values()), record, sep="\t", file = hout)

def annotation_gnomad(args) -> None:
    def reverse_complement(seq: str) -> str:
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def tbx_annotation_gnomad(tchr: str, tstart: int, tend: int, tref: str, talt: str, annotation_tbx) -> str:
        tabix_error_flag = False
        try:
            records = annotation_tbx.fetch(tchr, tstart, tend)
        except Exception as inst:
            print(f'{type(inst)}: {inst.args}', file = sys.stderr)
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tstart == int(record[1]) - 1 and tref == record[3] and talt == record[4]:
                    return ",".join(record[0:5])
                if tstart == int(record[1]) - 1 and tref == reverse_complement(record[3]) and talt == reverse_complement(record[4]):
                    return ",".join(record[0:5])

        return "-"
    
    annotation_tbx = pysam.TabixFile(args.gnomad) if args.gnomad is not None else None
    with open(args.input, 'r') as hin, open(args.output, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        print('\t'.join(header), "gnomAD", sep="\t", file = hout)

        for F in dreader:
            if "-" in F["GRCh38_pos"]:
                record = "-"
                print('\t'.join(F.values()), record, sep="\t", file = hout)
                continue
            tchr, tstart, tend, tref, talt = F["GRCh38_contig"], int(F["GRCh38_pos"]) - 1, int(F["GRCh38_pos"]), F["Ref"], F["Alt"]
            record = tbx_annotation_gnomad(tchr, tstart, tend, tref, talt, annotation_tbx) 

            print('\t'.join(F.values()), record, sep="\t", file = hout)
    
def arg_parser():
    # 1. Gene 
    # 2. RepeatMasker
    # 3. Contig size 
    # 4. misassembly 
    # 5. centromere
    # 6. segmental duplication
    # 7. point_mutation other
    # 8. gnomAD
    
    parser = argparse.ArgumentParser(prog = "point_mutation_annotation",
        description = "Add annotations to the result of point_mutation")
    
    subparsers = parser.add_subparsers()
    
    # Gene annotation
    gene = subparsers.add_parser('gene', help="add annotation of genes")
    gene.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated point_mutation result file")

    gene.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    gene.add_argument("--gff", "-g", type = str, default = None,
                        help = "Path to the tab-indexed liftoff gff file")
    
    gene.add_argument("--cgc", "-c", type = str, default = None,
                        help = "Path to the cancer gene census file")
    
    gene.add_argument("--cmrg", "-m", type = str, default = None,
                        help = "Path to the challenging medically relevant gene file")

    gene.add_argument("--mane", "-t", type = str, default = None,
                        help = "Path to the MANE summary file "
                               "(MANE.GRCh38.vX.Y.summary.txt.gz)")
    
    gene.set_defaults(func=annotation_gene)
    
    # RepeatMasker annotation
    rmsk = subparsers.add_parser('rmsk', help="add annotation of RepeatMasker result")
    rmsk.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated point_mutation result file")

    rmsk.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    rmsk.add_argument("--bed", "-b", type = str, default = None,
                        help = "Path to the tabix indexed RepetMasker bed file")
    
    rmsk.set_defaults(func=annotation_rmsk)
    
    # Contig size annotation
    size = subparsers.add_parser('size', help="add annotation of contig size")
    size.add_argument("--input", "-i", type = str,
                        help = "Path to the annotated point_mutation result file")

    size.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")

    size.add_argument("--hap1_fasta", "-f", type = str, default = None,
                        help = "Path to the assembly hap1 fasta file")
    
    size.add_argument("--hap2_fasta", "-g", type = str, default = None,
                        help = "Path to the assembly hap2 fasta file")
    
    size.set_defaults(func=annotation_size)
    
    
    # Misassembly
    misassembly = subparsers.add_parser('misassembly', help="add annotation whether SVs are in misassembly regions or not")
    misassembly.add_argument("--input", "-i", type = str,
                        help = "Path to the point_mutation result file")

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
                        help = "Path to the point_mutation result file")
    
    cen.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    cen.add_argument("--centromere", "-s", type = str,
                        help = "Path to the centromere tabixed bed file")
    
    cen.set_defaults(func=annotation_centromere)
    
    # Segmental duplication annotation
    segdup = subparsers.add_parser('segdup', help="add annotation of segmental duplications")
    segdup.add_argument("--input", "-i", type = str,
                        help = "Path to the point_mutation result file")
    
    segdup.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    segdup.add_argument("--segdup", "-s", type = str,
                        help = "Path to the segmental duplication tabixed bed file")
    
    segdup.set_defaults(func=annotation_segdup)

    # point mutation other annotation
    other = subparsers.add_parser('other', help="add annotation of point_mutation other")
    other.add_argument("--input", "-i", type = str,
                        help = "Path to the point_mutation result file")
    
    other.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    
    other.add_argument("--other", "-j", type = str,
                        help = "Path to the point_mutation other result file")
    
    other.set_defaults(func=annotation_point_mutation_other)
    
    # gnomAD annotation
    gnomad = subparsers.add_parser('gnomad', help="add annotation of gnomAD")
    gnomad.add_argument("--input", "-i", type = str,
                        help = "Path to the point_mutation result file")
    gnomad.add_argument("--output", "-o", type = str,
                        help = "Path to the output file")
    gnomad.add_argument("--gnomad", "-k", type = str,
                        help = "Path to the gnomAD tabixed VCF file")
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
