# script to make count tables from the consensus reads of mtDNA for single files

import os
import pysam
import pandas as pd
import re
import numpy as np
import pickle
import argparse

# read command line inputs
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input BAM file", type=str)
  parser.add_argument("--output", "-o", help="Output pickle object", type=str)
  parser.add_argument("--missmatch_ratio", "-r", help="Proportion of missmatch reads allowed at the mutated site per UMI", type=float)
  parser.add_argument("--start", "-s", help="Genomic coordintate of the 1st nucleotide of mtDNA", type=int, default = 1)
  parser.add_argument("--end", "-e", help="Length of the mitochondrial genome", type=int, default = 16569)
  parser.add_argument("--contig", "-c", help="Contig symbol of the mitochondrial genome", type=str, default = 'MT')
  args = parser.parse_args()

# get arguments
bamfile = args.input
output = args.output
contig = args.contig
start = args.start
end = args.end
threshold = args.missmatch_ratio

# function to generate count tables
def count_reads(bamfile, contig, start, end, miss_threshold):

    # read BAM file
    bam = pysam.AlignmentFile(bamfile, "rb")

    # get cell barcode
    cell_barcode = re.search('.+/(.+).bam', bamfile).group(1)

    # create final dictionary
    count_table = {}

    # add cell_barcode entry to the count dictionary
    count_table[cell_barcode] = {"A_counts": np.array([0.0]*(end-start+1)), "T_counts": np.array([0.0]*(end-start+1)),
                                 "C_counts": np.array([0.0]*(end-start+1)), "G_counts": np.array([0.0]*(end-start+1)),
                                "A_qual": np.array([0.0]*(end-start+1)), "T_qual": np.array([0.0]*(end-start+1)),
                                 "C_qual": np.array([0.0]*(end-start+1)), "G_qual": np.array([0.0]*(end-start+1)),
                               "reads_umi": np.array([0.0]*(end-start+1)), "depth": np.array([0.0000001]*(end-start+1))}

    # Iterate through the positions of the mitochondrial genome
    for pos in bam.pileup(contig, start, end, max_depth = 100000000, truncate = True):

        # dicitonary to check for discrepancies in UMIs
        umi_dict = {}

        # get current position
        position = pos.reference_pos

        # omit reads which are outside the positions of interest
        if position < start or position > end:
            continue

        # Iterate through the aligned reads to each position
        for read in pos.pileups:

            # filter out splicing and deletion reads as they do not have query.position
            if read.is_del or read.is_refskip:
                continue

            # get nucleotide at the base of interest
            nt = read.alignment.query_sequence[read.query_position]

            # filter out reads which have N as nucleotide
            if nt == "N":
                continue

            # get total number of raw reads at the position of interest
            raw_reads = read.alignment.get_tag("cd")[read.query_position]

            if raw_reads == 0:
                continue

            # get number of reads which do not agree with consensus and compute mismatch ratio
            missmatch_reads = read.alignment.get_tag("ce")[read.query_position]

            # missmatch ratio
            missmatch_ratio = missmatch_reads/raw_reads

            # if missmatch > threshold rule out read of interest
            if missmatch_ratio > miss_threshold:
                continue

            # get read UMI
            umi = re.sub('^.+:' ,"", read.alignment.query_name)

            # get cell barcode (it is stored in the read name of the BAM file)
            cell_barcode = re.search(".+_(.+)_.+", read.alignment.query_name).group(1)

            # get base quality
            base_qual = read.alignment.query_qualities[read.query_position]

            # add UMI to rule out discrepancies between different amplicons
            if not cell_barcode + "_" + umi in umi_dict:

                # add entry
                umi_dict[cell_barcode + "_" + umi] = [nt, raw_reads, base_qual]

                # add one count to the corresponding nucleotide
                count_table[cell_barcode][nt + "_counts"][position] += 1

                # add reads supporting the UMI
                count_table[cell_barcode]["reads_umi"][position] += raw_reads

                # add base quality
                count_table[cell_barcode][nt + "_qual"][position] += base_qual

                # sum up read depth
                count_table[cell_barcode]["depth"][position] += 1

            else:
                # if the 2 amplicons agree we simply add more raw reads and keep the highest base quality
                if nt == umi_dict[cell_barcode + "_" + umi][0]:

                    # add reads supporting the UMI
                    count_table[cell_barcode]["reads_umi"][position] += raw_reads

                    # update base quality (I keep the highest one among all amplicons of the same UMI)
                    diff_quality = base_qual-umi_dict[cell_barcode + "_" + umi][2]
                    if diff_quality > 0:
                        count_table[cell_barcode][nt + "_qual"][position] += diff_quality
                        umi_dict[cell_barcode + "_" + umi][2] = base_qual

                    # add raw reads to the UMI dict
                    umi_dict[cell_barcode + "_" + umi][1] += raw_reads

                # if there is a disagreement
                else:

                    wrong_nt = umi_dict[cell_barcode + "_" + umi][0]

                    if wrong_nt == "N":
                        continue
                    else:
                        # avoid counting more contigs from this UMI
                        umi_dict[cell_barcode + "_" + umi][0] = "N"

                        # substract the count
                        count_table[cell_barcode][wrong_nt + "_counts"][position] -= 1

                        #substract raw reads from that UMI
                        wrong_reads = umi_dict[cell_barcode + "_" + umi][1]
                        count_table[cell_barcode]["reads_umi"][position] -= wrong_reads

                        #substract base quality
                        wrong_qual = umi_dict[cell_barcode + "_" + umi][2]
                        count_table[cell_barcode][wrong_nt + "_qual"][position] -= wrong_qual


    # compute the average number of reads/UMI by dividing by the coverage (nUMIs)
    count_table[cell_barcode]["reads_umi"] = np.divide(count_table[cell_barcode]["reads_umi"],
                                                       count_table[cell_barcode]["depth"])

    # compute the average number of reads/UMI by dividing by the coverage (nUMIs)
    count_table[cell_barcode]["A_qual"] = np.divide(count_table[cell_barcode]["A_qual"],
                                                    count_table[cell_barcode]["A_counts"] + 0.00001)
    count_table[cell_barcode]["T_qual"] = np.divide(count_table[cell_barcode]["T_qual"],
                                                    count_table[cell_barcode]["T_counts"]+ 0.00001)
    count_table[cell_barcode]["C_qual"] = np.divide(count_table[cell_barcode]["C_qual"],
                                                    count_table[cell_barcode]["C_counts"]+ 0.00001)
    count_table[cell_barcode]["G_qual"] = np.divide(count_table[cell_barcode]["G_qual"],
                                                    count_table[cell_barcode]["G_counts"]+ 0.00001)

    return(count_table)


# inde BAM file
pysam.index(bamfile)

# generate final table
count_table = count_reads(bamfile, contig, start, end, threshold)

# save count_table as pickle for further processing
with open(output, "wb") as file:
    pickle.dump(count_table, file, protocol=pickle.HIGHEST_PROTOCOL)

# close file
file.close()
