#!/usr/bin/env python3

import argparse
import os
import re
import sys
from datetime import datetime
from statistics import median, mean, stdev

parser = argparse.ArgumentParser()
parser.add_argument("fna_file", help="the final fna file")
parser.add_argument("gff_file", help="the gff file")

args = parser.parse_args()

""" Global variables. """
seq_data = {}
tool_data = {}
total_predicted_genes = 0
gaps_data = []

N_STRETCH = "N" * 100



"""
    Runs over the given fasta file and fills the global seq_data dictionary
    with the sequence name as key and the sequence string and its length as
    values.
"""
def get_sequences_and_their_lengths(fna_file):
    print(str(datetime.now()) + " - Parsing fna file...")
    seq_lines = []
    with open(fna_file, "r") as fr:
        for line in fr:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if seq_lines:
                    seq_data[seq_name]["seq"] = "".join(seq_lines)
                    seq_lines.clear()
                seq_name = line[1:]
                seq_data[seq_name] = {}
                seq_data[seq_name]["length"] = 0

                continue

            seq_lines.append(line)
            seq_data[seq_name]["length"] += len(line)
        """ Add last sequence. """
        if seq_lines:
            seq_data[seq_name]["seq"] = "".join(seq_lines)
            seq_lines = None

    print(str(datetime.now()) + " - \t...done.")




"""
    Calculates basic statistics over the given list of lengths and returns them
    in a dictionary.
"""
def get_stats_dict(lengths_list):
    len_stats = {}
    len_stats["Number of seqs"] = len(lengths_list)
    len_stats["Number of bps"] = sum(lengths_list)
    len_stats["Median length"] = round(median(lengths_list), 3)
    len_stats["Average length"] = round(mean(lengths_list), 3)
    len_stats["Length shortest seq"] = min(lengths_list)
    len_stats["Length longest seq"] = max(lengths_list)
    if len(lengths_list) == 1:
        len_stats["Standard deviation"] = 0
    else:
        len_stats["Standard deviation"] = round(stdev(lengths_list), 3)

    return len_stats




"""
    Checks if a gap (defined by the given start and stop coordinate on the
    given sequence) contains a long stretch of Ns (defined via the N_STRETCH
    variable).
"""
def gap_has_not_too_many_ns(seq_name, gap_start, gap_stop):
    """ The previous_feature_end is the gap start here due to the 0 offset. """
    intergenic_region = seq_data[seq_name]["seq"][gap_start:gap_stop]
    if intergenic_region.find(N_STRETCH) >= 0:
        return False
    return True




"""
    Check if the given gap length crosses the reporting threshhold (of 3000 bp
    withough containing a stretch of a 100 consecutive Ns or more).
"""
def check_for_significant_gap(seq_name, previous_feature_end, feature_start):
    gap_start = previous_feature_end + 1
    gap_end = feature_start - 1
    gap_length = gap_end - gap_start + 1
    if (gap_length > 3000 and gap_has_not_too_many_ns(seq_name,
                                                  previous_feature_end,
                                                  feature_start)):
        gap_info = {}
        gap_info["Seq name"] = seq_name
        gap_info["Gap start"] = str(gap_start)
        gap_info["Gap end"] = str(gap_end)
        gap_info["Gap length"] = str(gap_length)
        """ Add gap to global gaps list. """
        gaps_data.append(gap_info)




"""
    Runs over the given GFF file and stores tool and feature data in global
    tool_data dictionary. Additionally updates the total genes counter and
    memorizes too long intergenig regions (gaps).
"""
def parse_gff_file(gff_file):
    """ Getting structural annotation stats. """
    print(str(datetime.now()) + " - Parsing gff file...")
    global total_predicted_genes
    seq_name = ""
    previous_feature_end = 0
    coding_bps = []
    with open(gff_file, "r") as fr:
        for line in fr:
            line = line.rstrip()

            if not line:
                continue

            if re.search("Parent", line):
                continue

            fields = line.split("\t")

            """ New contig. """
            if fields[0] != seq_name:
                if seq_name != "":
                    """ Count total coding bps. """
                    seq_data[seq_name]["coding_bps"] = sum(coding_bps)
                seq_name = fields[0]
                previous_feature_end = 0
                coding_bps = [0] * seq_data[seq_name]["length"]

            """ Update global gene counter. """
            total_predicted_genes += 1

            """ Store Tool and Feature Type data. """
            tool = fields[1]
            feature_type = fields[2]
            if tool not in tool_data:
                tool_data[tool] = {}
            if feature_type not in tool_data[tool]:
                tool_data[tool][feature_type] = {}
            """ Remember the sequences on which the tool and feature type were
            found. """
            if "found_on_seqs" not in tool_data[tool][feature_type]:
                tool_data[tool][feature_type]["found_on_seqs"] = {}
            tool_data[tool][feature_type]["found_on_seqs"][seq_name] = 1

            """ Remember the lengths of each feature type predicted by each
            tool. """
            feature_start = int(fields[3])
            feature_end = int(fields[4])
            for i in range((feature_start - 1), feature_end):
                try:
                    coding_bps[i] = 1
                except IndexError:
                    print(f'Off-contig coordinate ({i+1}) reported by {tool}. '
                          f'Contig: {seq_name} , Length: {seq_data[seq_name]["length"]}', file=sys.stderr)
                    print('Aborting!')
                    sys.exit(1)
            feature_length = feature_end - feature_start + 1
            if "feature_lengths" in tool_data[tool][feature_type]:
                tool_data[tool][feature_type]["feature_lengths"].append(feature_length)
            else:
                tool_data[tool][feature_type]["feature_lengths"] = [feature_length]

            """ Check if there's a significantly long intergenig region between
            this and the previous feature/gene. """
            if feature_end < previous_feature_end:
                continue
            check_for_significant_gap(seq_name, previous_feature_end,
                                      feature_start)
            previous_feature_end = feature_end

        """ Add coding bps of last sequence to total coding bps. """
        seq_data[seq_name]["coding_bps"] = sum(coding_bps)

    print(str(datetime.now()) + " - \t...done.")




"""
    Writes out the given data to the given file writer in a formatted table.
"""
def write_table(file_writer, title, column_names, data):
    file_writer.write(title + "\n")
    """ Write top frame line line. """
    fw.write("="*sum([len(column_name)+3 for column_name in column_names])
             + "\n")
    fw.write("\t".join(column_names) + "\n")
    """ Write header row separator line. """
    fw.write("-"*sum([len(column_name)+3 for column_name in column_names])
             + "\n")
    for data_dict in data:
        file_writer.write(str(data_dict[column_names[0]]))
        for i in range(1, len(column_names)):
            file_writer.write("\t" + str(data_dict[column_names[i]]))
        file_writer.write("\n")
    """ Write bottom frame line of table. """
    fw.write("="*sum([len(column_name)+3 for column_name in column_names])
             + "\n")
    fw.write("\n\n")




#####################################################
### Beginn of actual script work (parsing, etc.). ###
#####################################################
""" Length stats of nucleotide fasta file. """
get_sequences_and_their_lengths(args.fna_file)

""" Parse the given gff file. """
parse_gff_file(args.gff_file)


stats_out_file = args.fna_file[:args.fna_file.rfind("_")]
stats_out_file += "_structural_annotation_stats.tsv"
fw = open(stats_out_file, "w")


print(str(datetime.now()) + " - Calculating stats for fna file now...")
lengths_list = []
for seq_name in seq_data.keys():
    lengths_list.append(seq_data[seq_name]["length"])
stats_dict = get_stats_dict(lengths_list)
lengths_list.clear()
stats_dict["Data type"] = "'final_fasta'"

""" Remember how many bps got processed in total. """
processed_bps = stats_dict["Number of bps"]

""" Start putting together 1st table. """
table_title = "Processed Sequences Statistics"
table_columns = ["Data type", "Number of seqs", "Number of bps", \
                 "Length shortest seq", "Length longest seq", \
                 "Average length", "Median length", \
                 "Standard deviation"]
table_data = []
table_data.append(stats_dict)

tmp_lengths_list = []
print(str(datetime.now()) +
      " - Calculating stats for sequences with genes now...")
total_coding_bps = 0
for seq_name in seq_data.keys():
    if "coding_bps" in seq_data[seq_name]:
        lengths_list.append(seq_data[seq_name]["length"])
        total_coding_bps += seq_data[seq_name]["coding_bps"]
    else:
        tmp_lengths_list.append(seq_data[seq_name]["length"])
stats_dict = get_stats_dict(lengths_list)
lengths_list.clear()
stats_dict["Data type"] = "'sequences_with_genes'"
table_data.append(stats_dict)


lengths_list = tmp_lengths_list
tmp_lengths_list = None
bps_on_sequences_without_genes = 0
if lengths_list:
    print(str(datetime.now()) +
          " - Calculating stats for sequences without genes now...")
    stats_dict = get_stats_dict(lengths_list)
    lengths_list.clear()
    stats_dict["Data type"] = "'sequences_without_genes'"
    table_data.append(stats_dict)
    bps_on_sequences_without_genes = stats_dict["Number of bps"]

""" Write first table. """
write_table(fw, table_title, table_columns, table_data)
table_columns.clear()


""" Start putting together 2nd table. """
table_title = "Predicted Genes Statistics"
table_columns.extend(["Feature type", "Prediction method", \
                      "Number of predicted features", \
                      "Number of seqs", "Number of bps", \
                      "Length shortest seq", "Length longest seq", \
                      "Average length", "Median length", \
                      "Standard deviation"])
table_data.clear()

for tool in tool_data.keys():
    for feature_type in tool_data[tool]:
        print(str(datetime.now()) + " - Calculating stats for " +
              feature_type + " predicted by " + tool + " now...")
        stats_dict = get_stats_dict(
                       tool_data[tool][feature_type]["feature_lengths"])
        stats_dict["Number of predicted features"] = stats_dict["Number of seqs"]
        stats_dict["Number of seqs"] = (
                len(tool_data[tool][feature_type]["found_on_seqs"].keys()))
        stats_dict["Prediction method"] = "'" + tool + "'"
        stats_dict["Feature type"] = "'" + feature_type + "'"
        table_data.append(stats_dict)
""" Write 2nd table to file. """
write_table(fw, table_title, table_columns, table_data)
table_columns.clear()


""" Work on third table. """
table_title = "General Quality Info"
table_columns.extend(["Coding density", "Genes per 1M bp", "Seqs per 1M bp"])
table_data.clear()

coding_density = round((total_coding_bps / processed_bps * 100), 2)
stats_dict.clear()
stats_dict["Coding density"] = str(coding_density) + "%"
genes_per_million_bp = round((total_predicted_genes * 1000000 / processed_bps),
                             2)
stats_dict["Genes per 1M bp"] = str(genes_per_million_bp)
seqs_per_million_bp = round((len(seq_data) * 1000000 / processed_bps), 2)
stats_dict["Seqs per 1M bp"] = str(seqs_per_million_bp)
table_data.append(stats_dict)
write_table(fw, table_title, table_columns, table_data)

if gaps_data:
    """ Work on fourth table. """
    table_title = "Long Intergenic Regions"
    table_columns.clear()
    table_columns.extend(["Seq name", "Gap start", "Gap end", "Gap length"])
    write_table(fw, table_title, table_columns, gaps_data)

fw.close()

print(str(datetime.now()) + " - All stats have been calculated and successfully stored.")


