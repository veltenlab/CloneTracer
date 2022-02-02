#!/usr/bin/env Rscript

# make BED file from list of amplified genes in mutation library

# check if packages are installed
message("Checking if required packages are installed")
for (i in c("optparse","BiocManager", "tidyverse", "rtracklayer", "Rsamtools")){
  
  if(!requireNamespace(i, quietly = TRUE)){
    
    
    if(i %in% c("optparse", "BiocManager", "tidyverse")){
      
      install.packages(i, repos="https://ftp.fau.de/cran/")
      
    }else{BiocManager::install(i)}
    
  }
  
}


# parse command line arguments ---------------------------------------------------------------------

library("optparse")

# create arguments list
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the table containing the genomic coordinates of the sites of interest (one site per row). 
              It should contain three columns: 
                CHROM: chromosome.
                POS: genomic position of the mutation 
                symbol: gene name",
                metavar = "character"),
  make_option(c("-b", "--bam"), type = "character", default = NULL,
              help = "unfiltered BAM file",
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output bed file", metavar = "character")
)


# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  
  
  if (is.null(opt[[arg]])) {
    
    
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
    
  }
  
}


# check that all required parameters are provided
required_args <- c("input", "bam", "output")


for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}



# functions to run the script------------------------------------------------------------------------

library(tidyverse)
library(rtracklayer)
library(Rsamtools)

# read csv file with mutations
table_genes <- read_csv(opt$input)

# make GRanges object
granges <- GRanges(seqnames = table_genes$CHROM,
                   ranges = table_genes$POS)


# load BAM file
bam <- BamFile(opt$bam)


# change seqlevelsStyle of GRanges to match BAM file
seqlevelsStyle(granges) <- seqlevelsStyle(bam)


# export GRanges as BED file
export(granges, con = opt$output, format = "bed")


message("BED file created !")


