#!/usr/bin/env python3

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="""This script filters out
                                 nucleotide or protein sequences from one or
                                 more fasta files based on the presence of
                                 their IDs in the given GFF file.

                                 If the gene's ID is present in the GFF file
                                 the gene ID and its actual sequence will get
                                 printed out to stdout.""")
parser.add_argument("gff", metavar="gff_file",
                    help="file path to the gff file")
parser.add_argument("fasta_files", metavar="fasta_file", nargs="+",
                    help="file path(s) to the gene/protein fasta file(s)")
args = parser.parse_args()


"""Figure out if we are dealing with nucleotide or protein sequences."""
suffix = args.fasta_files[0][-3:]
allowed_start_codons = frozenset(["ATG", "GTG", "TTG"])

"""Parse the CDS gene IDs out of the given gff file."""
gene_to_start_type = {}
shortened_genes = {}
fr = open(args.gff, "r")
for line in fr:
    if not line or line[0] == "#" or line == "\n":
        continue
    start_type = ""
    fields = line.split("\t")
    if fields[2] != "CDS":
        continue
    gene_id = fields[-1].split(";", 1)[0].split("=", 1)[1]
    if "shortened" in fields[-1]:
        gene_id_parts = gene_id.rsplit("_", 2)
        shortened_info = fields[-1].split("shortened=", 1)[1].split(";", 1)[0]
        pos, coord = shortened_info.split()[1:]
        if pos == "start":
            gene_id_parts[1] = coord
        else:
            gene_id_parts[2] = coord
        original_gene_id = "_".join(gene_id_parts)
        shortened_genes[original_gene_id] = gene_id + "\t" + shortened_info
#        print("Found shortened info '" + shortened_info + "' for gene " +
#              gene_id + " (original: " + original_gene_id + ")",
#              file = sys.stderr)
        gene_id = original_gene_id
    if suffix == "faa":
        if "start_type" in fields[-1]:
            start_type = fields[-1].split("start_type=")[1].rstrip().split(";")[0]
    gene_to_start_type[gene_id] = start_type


def process_and_print_gene(fields, gene_id, seq):
#    print("papg got gene: " +  gene_id, file=sys.stderr)
#    print("seq:\n" + "\n".join([seq[0+i:i+70]
#                                for i in range(0, len(seq), 70)]),
#          file=sys.stderr)
    if gene_id in shortened_genes:
#        print("It needs to get shortened!", file=sys.stderr)
        new_gene_id, shortened_info = shortened_genes[gene_id].split("\t")
        old_gene_len = len(seq)
        gene_id_parts = new_gene_id.rsplit("_", 2)
        new_gene_len = (int(gene_id_parts[2]) -
                        int(gene_id_parts[1]) + 1)
        if suffix == "faa":
            new_gene_len = int(new_gene_len / 3)

        offset = old_gene_len - new_gene_len
        if suffix == "faa":
            offset += 1 # Since we removed the * from the end of the aa seq!
#        print("LEN: " +str(old_gene_len) + ", " +
#              str(new_gene_len) +
#              ", Offset: " + str(offset) + ", New gene: " +
#              new_gene_id, file = sys.stderr)
        seq = seq[offset:]
#        print("Shortened seq:\n" + "\n".join([seq[0+i:i+70]
#                                              for i in range(0, len(seq),
#                                                             70)]),
#              file=sys.stderr)
        if suffix == "faa":
            seq = "M" + seq[1:]
        fields.append("shortened=" + shortened_info)
        fields[0] = ">" + new_gene_id
        if "start" in shortened_info:
            fields[2] = gene_id_parts[1]
        else:
            fields[3] = gene_id_parts[2]
        if suffix == "fna":
            gene_id = new_gene_id
    if suffix == "fna":
        start_type = seq[0:3]
        if start_type not in allowed_start_codons:
            """This is a lazy implementation due to current time constraints.
            It obviously assumes that there's no bug in the CDS shortening
            parts of the pipeline.
            The preferred way would to get this information from the gene
            callers. Prodigal lists this information in it's GFF output, but
            the currently GeneMark unfortunately does not list it in its GFF
            output. It could get derived from its lst output (but that doesn't
            contain the score information). Maybe the new GeneMark version has
            this information provided with its GFF output..."""
            start_type = "Edge"
        gene_to_start_type[gene_id] = start_type
#        print("STORING: " + gene_id + " - " + start_type + "\n\n",
#              file = sys.stderr)
    fields.append("start_type=" + gene_to_start_type[gene_id])
    print(" # ".join(fields))
    print("\n".join(seq[0+i:i+70] for i in range(0, len(seq), 70)))



"""Now loop over the gene/protein fasta files."""
already_processed = set()
for fasta_file in args.fasta_files:
#    print("\n\nPROCESSING: " + fasta_file + "\n", file=sys.stderr)
    print_gene = False
    fields = []
    gene_id = ""
    seq = ""
    fr = open(fasta_file, "r")
    for line in fr:
        if not line or line[0] == "#" or line == "\n":
            continue
        if line[0] == ">":
            if print_gene:
                process_and_print_gene(fields, gene_id, seq)
            print_gene = False
            fields = line.rstrip().split(" # ")
            gene_id = fields[0][1:]
            if (gene_id in gene_to_start_type
                    and gene_id not in already_processed):
                print_gene = True
                seq = ""
                already_processed.add(gene_id)
                if gene_id in shortened_genes:
                    already_processed.add(shortened_genes[gene_id].split("\t")[0])
        else:
            if print_gene:
                seq += line.rstrip()
    if print_gene:
        process_and_print_gene(fields, gene_id, seq)
    fr.close()


"""Add start_type information to GFF."""
if suffix == "fna":
    tmp_file = os.path.dirname(os.path.abspath(args.gff)) + "/tmp.gff"
    fw = open(tmp_file, "w")
    fr = open(args.gff, "r")
    for line in fr:
        if not line or line[0] == "#" or line == "\n":
            continue
        fields = line.split("\t")
        if fields[2] == "CDS":
            attributes = fields[-1].rstrip().split(";")
            gene_id = attributes[0][3:].strip()
            attributes.append("start_type=" +
                                    gene_to_start_type[gene_id])
            fields[-1] = ";".join(attributes)
            line = "\t".join(fields) + "\n"
        fw.write(line)
    fr.close()
    fw.close()
    os.rename(tmp_file, args.gff)


