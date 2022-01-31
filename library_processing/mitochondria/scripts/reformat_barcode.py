import argparse
import pysam
import multiprocessing

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input bam file", type=str)
  parser.add_argument("--barcodes", "-b", help="List of valid barcodes output by cell_ranger", type=str)
  parser.add_argument("--output", "-o",help="Output BAM file", type=str)
  args = parser.parse_args()

# function to reformat barcodes and filter for reads presnent in cellranger barcodes
def reformat_barcodes(old_bam, new_bam, barcodes):

     old_bam = pysam.AlignmentFile(old_bam, 'rb')

     new_bam = pysam.AlignmentFile(new_bam, "wb", template = old_bam)

     with open(barcodes) as barcode_file_handle:
         content = barcode_file_handle.readlines()
     bc = [x.strip() for x in content]

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


# read BAM file
old_bam = args.input
new_bam = args.output
barcodes = args.barcodes
#cores = args.cores

# reformat barcodes and filter reads
reformat_barcodes(old_bam, new_bam, barcodes)
