# script to make count tables for mutations of interest

import os
import pysam
import pandas as pd
import re
import numpy
import pickle
import argparse

# read command line inputs
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input BAM file", type=str)
  parser.add_argument("--output", "-o", help="Output pickle object", type=str)
  parser.add_argument("--mutations", "-m", help="CSV file with mutations of interest. It should contain the following columns: symbol, contig, position, ref, alt", type=str)
  parser.add_argument("--missmatch_ratio", "-r", help="Proportion of missmatch reads allowed at the mutated site per UMI", type=float)
  args = parser.parse_args()

# get arguments
bamfile = args.input
output = args.output
variants = args.mutations
max_ratio = args.missmatch_ratio

# open BAM file
bam = pysam.AlignmentFile(bamfile, "rb")

# read file with mutation position
mutations = pd.read_csv(variants)


# function to process substitutions
def process_snv(bam, contig, start, miss_thres, ref, alt):

      # create dictionaries to store the mutational status, total raw reads and missmatch ratio of amplicons
      status, coverage, miss_reads = {}, {}, {}

      # loop over reads aligning to the position of interest.
      for position in bam.pileup(contig, start, start+1, max_depth = 1000000, truncate = True):

        # we only want to count nucleotides in the mutated site
        if not position.reference_pos == start:
              continue

        # iterate through the reads aligned in the position and count nucleotides
        for read in position.pileups:


            # filter out indels and splicing sites
            if read.is_del or read.is_refskip:
                continue

            # get total number of raw reads at the position of interest
            raw_reads = read.alignment.get_tag("cd")[read.query_position]

            # skip if no raw reads are found
            if raw_reads == 0:
              continue

            # get number of reads which do not agree with consensus and compute mismatch ratio
            missmatch_reads = read.alignment.get_tag("ce")[read.query_position]


            # missmatch ratio
            missmatch_ratio = missmatch_reads/raw_reads


            # if missmatch > 0.2 rule out read of interest
            if missmatch_ratio > miss_thres:
                continue

            # get nucleotide at the base of interest
            nt = read.alignment.query_sequence[read.query_position]


            # if the nucleotide is different from reference or alternative rule out the read
            if not nt in [ref, alt]:
                continue

            # get read UMI (it is encoded in the read name. There might be up to 4 entries per UMI, one per inner primer)
            umi = re.sub('^.+:' ,"", read.alignment.query_name)


            # get cell barcode (it is stored in the read name of the BAM file)
            cell_barcode = re.search(".+_(.+)_.+", read.alignment.query_name).group(1)


            if not cell_barcode in status:
                status[cell_barcode] = {}
                coverage[cell_barcode] = {}
                miss_reads[cell_barcode] = {}


            # if is the first amplicon from the UMI add it to the dictionaries
            if not umi in status[cell_barcode]:
                status[cell_barcode][umi] = []
                coverage[cell_barcode][umi] = []
                miss_reads[cell_barcode][umi] = []

            # add mutational status of UMI
            if nt == ref:
                status[cell_barcode][umi].append(0)
            elif nt == alt:
                status[cell_barcode][umi].append(1)

            # add number of raw reads from this UMI and amplicon
            coverage[cell_barcode][umi].append(raw_reads)

            # add missmatch ratio
            miss_reads[cell_barcode][umi].append(missmatch_reads)

        return status, coverage, miss_reads


# function to process indels
def process_indel(bam, contig, start, threshold, ref):

    # create dictionaries to store the mutational status, total raw reads and missmatch ratio of amplicons
      status, coverage, miss_reads = {}, {}, {}

    # loop over reads aligning to the position of interest.
      for position in bam.pileup(contig, start, start+1, max_depth = 100000000, truncate = True):


        # we only want to count nucleotides in the mutated site
        if not position.reference_pos == start:
              continue

        # iterate through the reads aligned in the position and count nucleotides
        for read in position.pileups:


            # filter out splicing sites
            if read.is_refskip:
                continue

            # get total number of raw reads at the position of interest
            raw_reads = read.alignment.get_tag("cd")[read.query_position]

            # skip if no raw reads are found
            if raw_reads == 0:
              continue


            # get number of reads which do not agree with consensus and compute mismatch ratio
            missmatch_reads = read.alignment.get_tag("ce")[read.query_position]


            # missmatch ratio
            missmatch_ratio = missmatch_reads/raw_reads


            # if missmatch > 0.2 rule out read of interest
            if missmatch_ratio > threshold:
                continue

            # get read UMI (it is encoded in the read name. There might be up to 4 entries per UMI, one per inner primer)
            umi = re.sub('^.+:' ,"", read.alignment.query_name)


            # get cell barcode (it is stored in the read name of the BAM file)
            cell_barcode = re.search(".+_(.+)_.+", read.alignment.query_name).group(1)


            if not cell_barcode in status:
                status[cell_barcode] = {}
                coverage[cell_barcode] = {}
                miss_reads[cell_barcode] = {}


            # if is the first amplicon from the UMI add it to the dictionaries
            if not umi in status[cell_barcode]:
                status[cell_barcode][umi] = []
                coverage[cell_barcode][umi] = []
                miss_reads[cell_barcode][umi] = []

            # if read is an indel then it's a mutation
            if read.indel:
                status[cell_barcode][umi].append(1)

                # add number of raw reads from this UMI and amplicon
                coverage[cell_barcode][umi].append(raw_reads)

                # add missmatch ratio
                miss_reads[cell_barcode][umi].append(missmatch_reads)

            # add reference reads
            elif read.alignment.query_sequence[read.query_position] == ref:
                status[cell_barcode][umi].append(0)

                # add number of raw reads from this UMI and amplicon
                coverage[cell_barcode][umi].append(raw_reads)

                # add missmatch ratio
                miss_reads[cell_barcode][umi].append(missmatch_reads)

        return status, coverage, miss_reads

# create dicitonary where counts will be stored
counts = {}

# get nucleotide count tables for each mutated site of interest
for i in mutations.index.values:

  # get contig name
  contig = mutations.iloc[i,1]

  # get position. -1 need to be added because coordinates in pysam are 0-indexed
  start = mutations.iloc[i,2]-1
  # collect the mutational site contig and position
  site = contig + ":" + str(mutations.iloc[i,2])

  # create dictionaries to store the mutational status, total raw reads and missmatch ratio of amplicons
  status, coverage, miss_reads = {}, {}, {}

  # get reference and alt alleles
  ref = mutations.iloc[i,3]
  alt = mutations.iloc[i,4]

  # in case of SNVs
  if len(alt) == 1:

        status, coverage, miss_reads = process_snv(bam, contig, start, max_ratio, ref, alt)

  # in case of deletions I have to select the reference nucleotide
  elif len(ref) > 1:

        # the reference has more than nucleotide, I select the last on which is where the deletion is positioned
        ref = list(ref)[1]

        status, coverage, miss_reads = process_indel(bam, contig, start, max_ratio, ref)


  # in case of indels is also a bit different
  elif len(alt) > 1:

        status, coverage, miss_reads = process_indel(bam, contig, start, max_ratio, ref)


   # aggregate cell counts and other metrics per cell
  for cell in status:

        # make sure that UMIs with more than one amplicon have the same allele
        # make counts of total reference and mutant counts

        if cell not in counts:
            counts[cell] = {site: {'ref': 0, 'alt': 0} }
        else:
            counts[cell][site] = {'ref': 0, 'alt': 0}

        # list to store inconsistent barcodes
        incons_barcodes = []

        for key, value in status[cell].items():

            state = list(set(value))

            # only take UMIs in which all amplicons have the same allele
            if(len(state) == 1):

                # count reference and alternative UMI counts
                if state[0] == 0:
                    counts[cell][site]['ref'] += 1
                else:
                    counts[cell][site]['alt'] += 1

            else:
                incons_barcodes.append(key)

        read_depth_mut = []
        read_depth_ref = []

        # add mean read depth at the position of interest for each allele
        for key, value in coverage[cell].items():

            # filter out barcodes which are inconsistent
            if key in incons_barcodes:
                continue

            if status[cell][key][0] == 0:
                read_depth_ref.append(numpy.sum(value))
            else:
                read_depth_mut.append(numpy.sum(value))

        counts[cell][site]["coverage_reference"] = numpy.mean(read_depth_ref)
        counts[cell][site]["coverage_mutant"] = numpy.mean(read_depth_mut)


        miss_mut = []
        miss_ref = []

        #add missmatch ratio of reference and mutant at the position of interest
        for key, value in miss_reads[cell].items():

            # filter out barcodes which are inconsistent
            if key in incons_barcodes:
                continue

            if status[cell][key][0] == 0:
                miss_ref.append(numpy.mean([int(a)/int(b) for a,b in zip(value, coverage[cell][key])]))
            else:
                miss_mut.append(numpy.mean([int(a)/int(b) for a,b in zip(value, coverage[cell][key])]))

        counts[cell][site]["miss_reference"] = numpy.mean(miss_ref)
        counts[cell][site]["miss_mutant"] = numpy.mean(miss_mut)

# save dictionary as pickle for further processing
with open(output, "wb") as file:
    pickle.dump(counts, file, protocol=pickle.HIGHEST_PROTOCOL)
