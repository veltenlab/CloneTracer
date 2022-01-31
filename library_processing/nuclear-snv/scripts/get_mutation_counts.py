import argparse
import pysam
import gzip
import pandas as pd
import os


#Get the mitochondrial sequence at each position in all the reads.

def get_count_table(input_bam, output, umi_tag, cell_tag, genes, barcodes):

  # read BAM file
  bam = pysam.AlignmentFile(input_bam, "rb")


  # read file with genes of interest. It contains gene name, chromosome and position
  hit_genes = pd.read_csv(genes)


  # read barcodes which passed the filters of cellranger and store them in a list
  barcodes = list(pd.read_csv(barcodes).iloc[:,0])


  # create an empty dictionary to store the arrays
  counts = {}


  # loop over the genomic sites targeted
  for i in hit_genes.index.values:

      # get contig name
      contig = str(hit_genes.iloc[i,1])


      # get position. -1 need to be added because coordinates in pysam are 0-indexed
      start = (hit_genes.iloc[i,2]-1)


      # collect the mutational site contig and position
      site = contig + ":" + str(hit_genes.iloc[i,2])


      # create an entry in the counts dictionary for the position of interest
      counts[site] = {}


      # loop over reads aligning to the position of interest.
      for position in bam.pileup(contig, start, start+1, max_depth = 1000000, truncate = True):

        # we only want to count nucleotides in the mutated site
        if not position.reference_pos == start:
              continue

        # iterate through the reads aligned in the position and count nucleotides
        for read in position.pileups:


          # filter out reads which do not have a cell barcode
          if not read.alignment.has_tag(cell_tag) or not  read.alignment.has_tag(umi_tag):
            continue


          # get cell barcode associated to the read
          cell = read.alignment.get_tag(cell_tag)


          # exclude reads from "cells" not present in the list of valid cell barcodes
          if not cell in barcodes:
              continue


          # Create and empty site for each cell where nucleotide counts are stored
          if not cell in counts[site]:
              counts[site][cell] =  {"A": 0,
                                     "C": 0,
                                     "T": 0,
                                     "G": 0,
                                     "N": 0,
                                     "indel": 0,
                                     "coverage": 0
                                    }


          # count indels
          if read.indel:
              counts[site][cell]["indel"] += 1


          # add count depending on the nucleotide present
          if not read.is_del and not read.is_refskip:
              counts[site][cell][read.alignment.query_sequence[read.query_position]] +=1


          # add one read count to the site to determine the final coverage
          counts[site][cell]["coverage"] += 1


  # open file to write the output
  f = open(output, "w")

  # write the file as csv.
  for site, data in counts.items():
    for cell, nucleotides in data.items():
        for base, count in nucleotides.items():
            f.write(site + "," + cell + "," + base + "," + str(count) +  "\n")

  f.close()




if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input_bam", help="Input bam file", type=str)
  parser.add_argument("--output", help="Output csv file", type=str)
  parser.add_argument("--cell_tag", help="Tag containing cell barcode (default: CB)", type=str, default ="CB")
  parser.add_argument("--umi_tag", help="Tag containing UMI barcode (default: UB)", type=str, default ="UB")
  parser.add_argument("--barcodes", help="File with valid cell barcodes one per line", type=str)
  parser.add_argument("--genes", help="csv file with the genomic sites targeted. The 2nd column should indicate the contig and the 3rd the position. It is key that the contig format is the equal to the BAM file", type=str)
  args = parser.parse_args()
  get_count_table(args.input_bam, args.output, args.umi_tag, args.cell_tag, args.genes, args.barcodes)
