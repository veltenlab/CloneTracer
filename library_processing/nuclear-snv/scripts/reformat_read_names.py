import argparse

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input bam file", type=str)
  parser.add_argument("--barcode", "-b", help= "Cell barcode from the single-cell BAM file", type=str)
  parser.add_argument("--output", "-o",help="Output BAM file", type=str)
  args = parser.parse_args()


import pysam
import os

# index input BAM file in order to be able to iterate through its reads
pysam.index(args.input)

# read BAM file
old_bam = pysam.AlignmentFile(args.input, 'rb', check_sq=False)

# open new BAM file where reads will be stored
new_bam = pysam.AlignmentFile(args.output, 'wb', template = old_bam)

# cell barcode
cell_barcode = args.barcode

# iterate through the reads and add an integer and the cell_barcode at the beginning of read name
count = 0

for read in old_bam.fetch(multiple_iterators=False, until_eof=True):

    count += 1

    read.query_name = str(count) + "_" + cell_barcode + "_" + read.query_name

    new_bam.write(read)


new_bam.close()
old_bam.close()

# remove the index of the old BAM file
old_index = args.input + ".bai"

rm_index = "rm " + old_index

os.system(rm_index)

print("Done! {} reads written".format(count))
