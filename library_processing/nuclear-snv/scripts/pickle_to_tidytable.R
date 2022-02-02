#!/usr/bin/env Rscript

# convert count table from pickle python object to dataframe saved as RDS


message("Checking if required packages are installed")
for (i in c("optparse","tidyverse", "reticulate")){
  
  if(!requireNamespace(i, quietly = TRUE)){
    
      
    install.packages(i, repos="https://ftp.fau.de/cran/")
      
  }
  
}


# parse command line arguments ---------------------------------------------------------------------

library("optparse")

# create arguments list
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Count table as pickle object output by rule make_count_table",
              metavar = "character"),
  make_option(c("-v", "--variants"), type = "character", default = NULL,
              help = "csv file with amplified variants in the mutation library. It should contain the following columns:
              CHROM: chromosome,
              POS: position,
              symbol: gene name",
              metavar = "character"),
  make_option(c("-b", "--barcodes"), type = "character", default = NULL,
              help = "File with high quality barcodes output by cellranger (one per line with no column name) ", metavar = "character"),
  make_option(c("-r", "--reads_cell"), type = "character", default = NULL,
              help = "txt file with reads/cell output by rule get_read_cells of the mutation libarary pipeline ", metavar = "character"),  
  make_option(c("-t", "--min_mutant_umis"), type = "integer", default = 1,
              help = "Minimum number of UMIs to label a cell as mutant (Default: 1)", metavar = "integer"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output rds file", metavar = "character")
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
required_args <- c("input", "variants", "barcodes", "reads_cell", "min_mutant_umis","output")


for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}

# script to make count table as RDS -----------------------------------------------------------------------

# load packages
library(reticulate)
library(tidyverse)

# import pandas to load pickle object
pd <- import("pandas")


# import pickle object as nested R list
pickle <- pd$read_pickle(opt$input)


# load variants
mutations <- read_csv(opt$variants) 

# necessary in case we want to genotype mitochondrial sites from the transcriptome
if(mutations$CHROM[1] != "MT"){
  
  mutations <- mutations %>% mutate(CHROM = gsub("chr","",CHROM))
  
}
              

# load high-quality barcodes
barcodes <- read_csv(opt$barcodes,
                     col_names = F) %>% pull(X1)

# get reads/cell in the mutation library
reads_cell <- read_csv(opt$reads_cell,
                       col_names = F) %>% 
                separate(X1, into = c("reads_mut_lib", "cell_barcode"), sep = " ") %>% 
                mutate(cell_barcode = gsub(".+Z:", "", cell_barcode))

# get threshold UMIs to assign mutant cells
min_umis <- opt$min_mutant_umis


# get mutated sites in the same nomenclature as in the pickle object
sites <- paste(mutations$CHROM, mutations$POS, sep = ":")


# iterate through all cells
list_rows <- lapply(1:length(names(pickle)), function(i){
  
  # get cell barcode
  cell_barcode <- names(pickle)[i] 
  
  # create an entry for each of the mutations (if no info is available an NA will be filled in)
  mutation_list <- lapply(1:length(sites), function(x){

    # get mutation
    mut_site <- sites[x]

    # if the mutation is covered in the cell the info is stored
    if(mut_site %in% names(pickle[[cell_barcode]])){
      
      # get infor for the site of interest
      cell_data_list <- pickle[[cell_barcode]][[mut_site]]

      
      # transorm list to data frame
      cell_data_table <- data.frame(t(sapply(cell_data_list, c)))

      
    # if cell is not covered then a row of NAs is added
    }else{cell_data_table <- data.frame(ref = NA, alt = NA, coverage_reference = NA, coverage_mutant = NA, miss_reference = NA, miss_mutant = NA)}

    
    # add genomic coordinates
    cell_data_table["site"] <- mut_site

    # return final table
    cell_data_table

  })

  # bind rows for all genomic sites
  mut_table <- do.call("bind_rows", mutation_list)

  # add cell_barcode
  mut_table["cell_barcode"] <- cell_barcode
  
  mut_table
  
})


# merge single_cell tables into one dataframe
count_table_init <- do.call("bind_rows", list_rows)


# add cells which have no coverage for any site as dropout
dropout_barcodes <- barcodes[which(!barcodes %in% names(pickle))]


# make dropout dataframe to left_join
dropout_table <- data.frame(ref = NA, alt = NA, coverage_reference = NA, coverage_mutant = NA, miss_reference = NA, miss_mutant = NA,
                            site = rep(sites, times = length(dropout_barcodes)),
                            cell_barcode = rep(dropout_barcodes, each = length(sites)))

# make final count table
count_table <- count_table_init %>% bind_rows(dropout_table) %>% 
                separate(site, into = c("CHROM", "POS"), sep = ":") %>% 
                mutate(POS = as.integer(POS)) %>% 
                left_join(mutations %>% dplyr::select(-ref, -alt)) %>% 
                mutate(status = case_when(
                                          is.na(ref) & is.na(alt) ~ "dropout",
                                          alt >= min_umis ~ "mutant",
                                          ref > 0 ~ "reference")) %>% 
                left_join(reads_cell) %>% 
                mutate(umis = ref+alt,
                       coverage_reference = ifelse(is.nan(coverage_reference), 0, coverage_reference),
                       coverage_mutant = ifelse(is.nan(coverage_mutant), 0, coverage_mutant)) %>% 
                mutate(avg_coverage = (ref*coverage_reference+alt*coverage_mutant)/(ref+alt)) %>% 
                dplyr::select(cell_barcode, symbol, status,everything()) %>% 
                as_tibble()


# save final table as RDS
saveRDS(count_table, file = opt$output)





