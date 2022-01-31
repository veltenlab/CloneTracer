import argparse

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input bam file", type=str)
  parser.add_argument("--errors_umi", "-e", help= "Missmatches allowed when sorting by UMI using GroupReadsByUmi", type=int, default=0)
  parser.add_argument("--min_reads", "-r", help= "Minimum reads for a UMI to be considered", type=int, default=1)
  parser.add_argument("--mapping_quality", "-m", help= "Minimum mapping quality for a read to be considered in GroupReadsByUmi", type=int, default=0)
  parser.add_argument("--output", "-o",help="Output BAM file", type=str)
  args = parser.parse_args()


import re
import os

# get arguments from parser
input_bam = args.input
output_bam = args.output
errors_umi = args.errors_umi
min_reads = args.min_reads
map_qual = args.mapping_quality

# get barcode of the cell of interest
barcode = re.search(".+split_bams/(.+).bam", input_bam).group(1)

# create directories to store intermediate BAM files and log files
# get parent directory
temp_directory = re.search("(.+)split.+", input_bam).group(1)

# get logs folder
log_dir = re.search("(.+)count.+", input_bam).group(1) + "logs/"

# new BAM and log directories
sorted_dir = temp_directory + "sorted/"
consensus_dir = temp_directory + "consensus/"
renamed_dir = temp_directory + "renamed/"
output_dir = temp_directory + "processed/"
sorted_log_dir = log_dir + "sort_bams/"
consensus_log_dir = log_dir + "consensus_bams/"
renamed_log_dir = log_dir + "rename_bams/"
sort_final_log_dir = log_dir + "sort_final_bams/"


# create directories if they do not exist
for dir in [sorted_dir, consensus_dir, renamed_dir, output_dir, sorted_log_dir, consensus_log_dir, renamed_log_dir, sort_final_log_dir]:

    if not os.path.isdir(dir):
        os.makedirs(dir)

# make log files for each of the steps
sorted_log = log_dir + "sort_bams/" + barcode + ".log"
consensus_log = log_dir + "consensus_bams/" + barcode + ".log"
renamed_log = log_dir + "rename_bams/" + barcode + ".log"
sort_final_log = log_dir + "sort_final_bams/" + barcode + ".log"

# create intermediate BAM files
sorted_bam = sorted_dir + barcode + ".bam"
consensus_bam = consensus_dir + barcode + ".bam"
renamed_bam = renamed_dir + barcode + ".bam"

# sort UMI
sort_umi = "fgbio GroupReadsByUmi -i " + input_bam + " -o " + sorted_bam + " -s Identity -e " + str(errors_umi) + " -m 0 -t UB 2> " + sorted_log
os.system(sort_umi)

# consensus read
consensus_seq = "fgbio CallMolecularConsensusReads -i " + sorted_bam + " -o " + consensus_bam + " -M " + str(min_reads) + " -t UB 2> " + consensus_log
os.system(consensus_seq)

# reformat read name
reformat_name = "python scripts/reformat_read_names.py -i " + consensus_bam + " -o " + renamed_bam + " -b " + barcode + " 2> " + renamed_log
os.system(reformat_name)

# sort BAM by query name
sort_bam = "picard SortSam INPUT=" + renamed_bam + " OUTPUT=" + output_bam + " SORT_ORDER=queryname TMP_DIR=./tmp 2> " + sort_final_log
os.system(sort_bam)

# remove temporary files
rm = " ".join(["rm", sorted_bam, consensus_bam, renamed_bam])
os.system(rm)
