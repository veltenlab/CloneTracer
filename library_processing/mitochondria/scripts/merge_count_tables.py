# script to make count tables for the consensus reads of mtDNA

import os
import pysam
import pandas as pd
import re
import glob
import numpy as np
import pickle
import argparse

# read command line inputs
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--tables_directory", "-t", help="Input BAM file", type=str)
  parser.add_argument("--bams_directory", "-b", help="Input BAM file", type=str)
  parser.add_argument("--output", "-o", help="Output pickle object", type=str)
  parser.add_argument("--keep_temp", "-k", help="Boolean indicating to keep temporal tables and BAM files (debugging)", type=bool, default = False)
  args = parser.parse_args()

# get arguments
table_directory = args.tables_directory
bam_directory = args.bams_directory
output = args.output
keep_temp = args.keep_temp

# get single cell count tables
table_files = glob.glob(table_directory + "*pickle")

# put all dictionaries into a list
table_list = [pickle.load(open(file, "rb")) for file in table_files]

# save count_table as pickle for further processing
with open(output, "wb") as file:
    pickle.dump(table_list, file, protocol=pickle.HIGHEST_PROTOCOL)


# remove single-cell BAMs and tables unless otherwise specified

if not keep_temp:

    # # remove single-cell count tables
    # rm_tables = "rm " + table_directory
    #
    # os.system(rm_tables)

    # remove single-cell BAM files
    rm_bams = "rm " + bam_directory

    os.system(rm_bams)
