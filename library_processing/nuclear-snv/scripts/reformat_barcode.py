import argparse

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input bam file", type=str)
  parser.add_argument("--barcodes", "-b", help="List of valid barcodes output by cell_ranger", type=str)
  parser.add_argument("--output", "-o",help="Output BAM file", type=str)
  args = parser.parse_args()


import pysam

# Read in barcodes filtered by cellranger
with open(args.barcodes) as barcode_file_handle:
    content = barcode_file_handle.readlines()
bc = [x.strip() for x in content]

# read BAM file
old_bam = pysam.AlignmentFile(args.input, 'rb')

# open new BAM file where reads will be stored
new_bam = pysam.AlignmentFile(args.output, 'wb', template = old_bam)

# iterate through the reads and add '-1' to the cell barcode
# it only writes reads coming from cell barcodes present in the cellranger output
for read in old_bam.fetch():

    if read.has_tag('CB') and read.get_tag('CB') + '-1' in bc:

        # get current tag
        tag = read.get_tag('CB')

        # set the new tag
        read.set_tag('CB', tag + '-1')

        # write the read to new BAM file
        new_bam.write(read)
