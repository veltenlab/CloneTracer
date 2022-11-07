# script to get the number of reads/UMI

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
  parser.add_argument("--output", "-o", help="Output text file", type=str)
  args = parser.parse_args()

# get arguments
bamfile = args.input
output = args.output

# open BAM file
bam = pysam.AlignmentFile(bamfile, "rb")

# write number of reads/UMI as txt file
with open(output, "w") as file:

    for read in bam.fetch(multiple_iterators = False):

        file.write(str(read.get_tag("cD")) + "\n")

# close file
file.close()
