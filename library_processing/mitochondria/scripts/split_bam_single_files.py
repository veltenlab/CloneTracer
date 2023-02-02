# script to split aligned BAM file to single-cell BAM files for read deduplication

import argparse
import pysam
import gzip
import pandas as pd
import os
import multiprocessing
import math

# read command line inputs
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input_bam", "-i", help="Input bam file", type=str)
  parser.add_argument("--output_directory", "-o", help="Directory where single cell BAM files will be written", type=str)
  parser.add_argument("--cell_tag", "-t", help="Tag containing cell barcode (default: CB)", type=str, default ="CB")
  parser.add_argument("--barcodes", "-b", help="File with valid cell barcodes one per line", type=str)
  parser.add_argument("--present_barcodes", "-p", help="File with cell barcodes present in the BAM file (output from rule get_barcodes)", type=str)
  parser.add_argument("--max_files", "-m", help="Maximum number of opened files allowed by the operative system", type=int)
  parser.add_argument("--cores", "-c", help="Number of cores. More than one is recommended particularly if the number of cells is much larger than the number of maximum opened files the operative system allows", type=int)
  args = parser.parse_args()

# assign them to variables
bam_file = args.input_bam
output_directory = args.output_directory
cell_tag = args.cell_tag
bcfile = args.barcodes
barcodes = args.present_barcodes
cores = args.cores

if args.max_files-cores-2 <= 1:
    max_files = args.max_files

else:
    max_files = args.max_files-cores-2

# function to write reads
def write_bam(barcodes, outfiles, bamfile):

    # # read barcoded bam file
    bam = pysam.AlignmentFile(bamfile, "rb")

    # create single cell bam files and keep them open (in pysam as soon as you
    # close a BAM file you cannot add reads).
    bam_cell_files = [pysam.AlignmentFile(file, "wb", template = bam) for file in outfiles]

    # iterate through the reads and assign them to the single cell BAM files
    for read in bam.fetch(multiple_iterators=False):

        # remove reads which do not have a cell barcode tag
        if not read.has_tag("CB") or not read.has_tag("UB"):
            continue

        # get read cell tag
        cell_tag = read.get_tag("CB")

        # rule out reads whose barcode is not in the list of filtered barcodes
        if not cell_tag in barcodes:
            continue

        # get the position of the cell barcode
        index = barcodes.index(cell_tag)

        # get the corresponding single-cell BAM file
        sc_bam = bam_cell_files[index]

        # write read to the file
        sc_bam.write(read)

    # close single cell bam files
    for file in bam_cell_files:
        file.close()

    # close barcoded bam file
    bam.close()

# create directory if it doesn't exist
if not os.path.isdir(output_directory):
    os.makedirs(output_directory)

# Read in the barcodes filtered by cellranger
with open(bcfile) as barcode_file_handle:
    content = barcode_file_handle.readlines()
cellranger_barcodes = [x.strip() for x in content]

# read barcodes present in the file
with open(barcodes) as barcode_file:
    content = barcode_file.readlines()
present_barcodes = [x.strip() for x in content]

# filter for cellranger barcodes present in the BAM file
final_barcodes = [barcode for barcode in present_barcodes if barcode in cellranger_barcodes]

# create file paths for each single cell bam file
bam_file_paths = [output_directory + "/" + bc1 + ".bam" for bc1 in final_barcodes]

# make dicitonaries with barcodes and output files per process
barcodes_list, outfiles_list = [], []


if cores > 1:
    # if the total barcodes < max open files i just distribute the barcodes in the cores
    if len(final_barcodes) < max_files:

        files_core = round(len(final_barcodes)/cores)

        # add files to dictionary
        for i in range(cores):

            # the first core is different
            if i == 0:

                barcodes_list.append(final_barcodes[0:files_core*(i+1)])

                outfiles_list.append(bam_file_paths[0:files_core*(i+1)])

            elif not i == cores:

                barcodes_list.append(final_barcodes[files_core*i:files_core*(i+1)])

                outfiles_list.append(bam_file_paths[files_core*i:files_core*(i+1)])

            else:

                barcodes_list.append(final_barcodes[files_core*i:])

                outfiles_list.append(bam_file_paths[files_core*i:])
        # run parallel process
        with multiprocessing.Pool(cores) as pool:
            pool.starmap(write_bam, zip(barcodes_list, outfiles_list, [bam_file] * cores))

    # if there are more barcodes than max open files more iterations need to be done
    else:
        # compute how many single cell files/core should be processed.
        files_core = math.floor(max_files/cores)

        # get how many files/iteration will be open
        files_iteration = files_core*cores

        # get how many iterations should be computed
        iterations = math.ceil(len(final_barcodes)/files_core)

        # add files to dictionary
        for i in range(iterations):

            # the last core is different
            if i == 0:

                barcodes_list.append(final_barcodes[0:files_core*(i+1)])

                outfiles_list.append(bam_file_paths[0:files_core*(i+1)])

            elif not i == cores:

                barcodes_list.append(final_barcodes[files_core*i:files_core*(i+1)])

                outfiles_list.append(bam_file_paths[files_core*i:files_core*(i+1)])

            else:

                barcodes_list.append(final_barcodes[files_core*i:])

                outfiles_list.append(bam_file_paths[files_core*i:])

        # compute how many blocks of n = cores iterations have to be done
        blocks = math.ceil(iterations/cores)

        for i in range(blocks):

            if i == 0:

                # run parallel process
                with multiprocessing.Pool(cores) as pool:
                    pool.starmap(write_bam, zip(barcodes_list[0:cores],
                                                outfiles_list[0:cores],
                                                [bam_file] * cores))
            elif i != blocks-1:

                # run parallel process
                with multiprocessing.Pool(cores) as pool:
                    pool.starmap(write_bam, zip(barcodes_list[cores*i:cores*(i+1)],
                                                outfiles_list[cores*i:cores*(i+1)],
                                                [bam_file] * cores))

            else:

                # final iterations
                final_iter = iterations-((blocks-1)*cores)

                # run parallel process
                with multiprocessing.Pool(cores) as pool:
                    pool.starmap(write_bam, zip(barcodes_list[cores*i:],
                                                outfiles_list[cores*i:],
                                                [bam_file] * final_iter))

# in case only one core is provided
else:

    # compute number of iterations needed
    iterations = math.ceil(len(final_barcodes)/max_files)

    # add files to dictionary
    for i in range(iterations):

        # the last core is different
        if i == 0:

            write_bam(final_barcodes[0:max_files*(i+1)],
                      bam_file_paths[0:max_files*(i+1)],
                      bam_file)

        elif not i == iterations:

            write_bam(final_barcodes[max_files*i:max_files*(i+1)],
                      bam_file_paths[max_files*i:max_files*(i+1)],
                      bam_file)

        else:

            write_bam(final_barcodes[max_files*i:],
                      bam_file_paths[max_files*i:],
                      bam_file)


print("BAM split completed!")
