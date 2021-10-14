#!/usr/bin/env Rscript

# Subset the annotated mutation list for primer design


# set the library folder
.libPaths("/nfs/users2/lvelten/sbeneyto/.conda/envs/exome/lib/R/library")


# parse command line arguments ---------------------------------------------------------------------

library(optparse)

# create arguments list
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "csv file with the mutations of interest. It must contain the following 3 columns: 
                      CHROM: chromosome,
                      POS: position,
                      symbol: gene name",
              metavar = "character"),
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "text file with gene names that want to be amplified and the position of the site
                      separated by an underscore. One mutation per line. Position is necessary since there
                      could be more than one mutation per gene", 
              metavar = "character"),
  make_option(c("-d", "--out_directory"), type = "character", default = NULL,
              help = "directory where table with selected genes will be stored", 
              metavar = "character"),
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
required_args <- c("input", "genes", "out_directory")


for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}



# code to run the script .............................................................................

# load packages
library(tidyverse)


# read table with all annotated mutations
mutations_unfilt <- read_csv(opt$input)


# load selected mutations
mutations <- read_delim(opt$gene, col_names = F, delim = "\n") %>%
              separate(X1, into = c("symbol", "position"))


# subset the selected mutations
mutations_filt <- mutations_unfilt %>% filter(symbol %in% mutations$symbol &
                                              POS %in% mutations$position)


# create a new directory to store info related to primers 

if(!dir.exists("primers")){
  
  dir.create(opt$out_directory)
  
}


# save table with selected variants in the primers folder
write_csv(mutations_filt, file = paste0(opt$out_directory, "selected_variants.csv"))







