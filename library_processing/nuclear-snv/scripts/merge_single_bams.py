import argparse
import pysam
import gzip
import pandas as pd
import os
import glob
import re
import math
import multiprocessing

# read command line inputs
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input_directory", "-i", help="Directory where BAM files are stored", type=str)
  parser.add_argument("--output", "-o", help="Merged BAM file", type=str)
  parser.add_argument("--outdir", "-d", help="Directory where final BAM file will be stored", type=str)
  parser.add_argument("--max_files", "-m", help="Maximum number of files allowed to be open by OS", type=int)
  parser.add_argument("--cores", "-c", help="Number of cores. If > 10k cells are presented 2 are required (otherwise there would be too many open files)", type=int)
  parser.add_argument("--keep_temp", "-k", help="Boolean indicating whether to keep temporal folders", type=bool, default=False)
  args = parser.parse_args()

# list of initial BAM
bam_dir = args.input_directory
initial_bams = glob.glob(bam_dir + "*.bam")
merged_bam = args.output
out_dir = args.outdir
cores = args.cores
max_files = args.max_files-cores-2
max_files_samtools = 996
keep_temp = args.keep_temp

if args.max_files-cores-2 <= 1:
    max_files = args.max_files

else:
    max_files = args.max_files-cores-2

# function to merge single cell bam files
def merge_bams(bam_list, merged_bam):

    shell_merge = "samtools merge -n" + " " + merged_bam + " " + " ".join(bam_list)

    os.system(shell_merge)

# non-empty BAM list
final_bams = []

for file in initial_bams:

    # create index for the file
    pysam.index(file)

    # open the file to iterate through reads
    bam = pysam.AlignmentFile(file, "rb", check_sq =False)

    count = 0

    for read in bam.fetch(multiple_iterators=False, until_eof=True):

        count += 1

        if count > 0:

            final_bams.append(file)

            break

bams_list = []

# if the total barcodes < max open files i just distribute the barcodes in the cores
if len(final_bams) < max_files and len(final_bams)/cores < max_files_samtools:

    files_core = round(len(final_bams)/cores)

    temp_bams = [merged_bam + "_" + str(i) + ".bam" for i in range(cores)]

    # add files to dictionary
    for i in range(cores):

        # the first core is different
        if i == 0:

            bams_list.append(final_bams[0:files_core*(i+1)])

        elif not i == cores-1:

            bams_list.append(final_bams[files_core*i:files_core*(i+1)])

        else:

            bams_list.append(final_bams[files_core*i:])

    # run parallel process
    with multiprocessing.Pool(cores) as pool:
        pool.starmap(merge_bams, zip(bams_list, temp_bams))

    # create shell string to merge intermediate files into final BAM file
    shell_merge = "samtools merge -n" + " " + merged_bam + " " + " ".join(temp_bams)

    # execute merge command
    os.system(shell_merge)

    # remove intermediate BAM files
    rm = "rm " + " ".join(temp_bams)

    os.system(rm)

# if there are more barcodes than max open files more iterations need to be done
elif len(final_bams) > max_files and max_files_samtools*cores < max_files or len(final_bams) < max_files and max_files_samtools*cores < max_files:

    # each core will merge as many files as samtools merge allows (1000 in our cluster)
    files_core = max_files_samtools

    # iterations
    iterations = math.ceil(len(final_bams)/files_core)

    # get file names of temporal bam files
    temp_bams = [merged_bam + "_" + str(i) + ".bam" for i in range(iterations)]

    # add files to dictionary
    for i in range(iterations):

        # the first core is different
        if i == 0:

            bams_list.append(final_bams[0:files_core*(i+1)])

        elif not i == iterations-1:

            bams_list.append(final_bams[files_core*i:files_core*(i+1)])

        else:

            bams_list.append(final_bams[files_core*i:])

    # get how many blocks of n = cores iterations have to be computed.
    blocks = math.ceil(iterations/cores)

    # for each block of iterations = cores execute them in parallel
    for i in range(blocks):

        if i == 0:

            # run parallel process
            with multiprocessing.Pool(cores) as pool:
                pool.starmap(merge_bams, zip(bams_list[0:cores], temp_bams[0:cores]))

        elif not i == blocks-1:

            # run parallel process
            with multiprocessing.Pool(cores) as pool:
                pool.starmap(merge_bams, zip(bams_list[cores*i:cores*(i+1)], temp_bams[cores*i:cores*(i+1)]))

        else:

            # run parallel process
            with multiprocessing.Pool(cores) as pool:
                pool.starmap(merge_bams, zip(bams_list[cores*i:], temp_bams[cores*i:]))

    # create shell string to merge intermediate files into final BAM file
    shell_merge = "samtools merge -n" + " " + merged_bam + " " + " ".join(temp_bams)

    # execute merge command
    os.system(shell_merge)

    # remove intermediate BAM files
    rm = "rm " + " ".join(temp_bams)

    os.system(rm)

# in cases the max_files is very low or just 1 core is chosen
else:

    temp_bams = [re.sub(".bam", "",merged_bam) + "_" + str(i) + "_temp.bam" for i in range(len(final_bams)-1)]

    # do as many iterations as non-empty BAM files are present
    for i in range(len(final_bams)-1):

        # first iteration is a bit different
        if i == 0:

            merge_bams(final_bams[0:2], temp_bams[i])

        else:

            shell_merge = "samtools merge -n " + temp_bams[i] + " " + temp_bams[i-1] + " " + final_bams[i+1]

            os.system(shell_merge)

    # rename last bam file which contains all the reads
    mv_shell = "mv " + temp_bams[-1] + " " + merged_bam
    os.system(mv_shell)

    # remove intermediary bam files
    int_bams = glob.glob("*_temp.bam")
    rm_shell = "rm " + " ".join(int_bams)
    os.system(rm_shell)

# remove intermediate empty folders
if not keep_temp:

    temp_dir = out_dir + "/temp_bams"

    rm = "rm -r " + temp_dir

    os.system(rm)
