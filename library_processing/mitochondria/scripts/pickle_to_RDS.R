#!/usr/bin/env Rscript

# script to conver the pickle count object to list of single-cell matrixes for variant calling as well as making a SummarisedExperiment object

# parse command line arguments ---------------------------------------------------------------------

library(optparse)

# create arguments list
option_list = list(
  make_option(c("-i", "--input_directory"), type = "character", default = NULL,
              help = "Directory where single-cell pickle objects are stored",
              metavar = "character"),
  make_option(c("-b", "--barcodes"), type = "character", default = NULL,
              help = "File with high quality barcodes output by cellranger (one per line with no column name) ", metavar = "character"),
  make_option(c("-s", "--summarised_experiment"), type = "character", default = NULL,
              help = "output rds file for summarised experiment object", metavar = "character"),
  make_option(c("-r", "--reads_cell"), type = "character", default = NULL,
              help = "txt file containing reads/cell as output by rule get_reads_cell", metavar = "character"),  
  make_option(c("-o", "--outRDS"), type = "character", default = NULL,
              help = "output rds file for list of single-cell count tables (compatible with mutaClone package functions)", metavar = "character"),
  make_option(c("-c", "--cores"), type = "integer", default = NULL,
              help = "number of cores", metavar = "character"),
  make_option(c("-k", "--keep_temp"), type = "logical", default = FALSE,
              help = "whether to keep temporal single-cell count_table.pickle and BAM files. Default: FALSE")
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
required_args <- c("input_directory", "barcodes", "summarised_experiment", "reads_cell", "outRDS", "cores", "keep_temp")


for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}

# functions for the script ----------------------------------------------------------------------

# make mito count tables from SummarisedExperiment object
make_mito_count_table <- function(object){
  
  
  
  # reshape matrices into an output compatible with mitoClone functions
  list_clean_count_tables <- mclapply(1:length(colnames(object)), function(x){
    
    
    table <- data.frame(A = assays(object)$A_counts[,colnames(object)[x]],
                        T = assays(object)$T_counts[,colnames(object)[x]],
                        C = assays(object)$C_counts[,colnames(object)[x]],
                        G = assays(object)$G_counts[,colnames(object)[x]]) %>% 
      as.matrix()
    
  }, mc.cores = 8)
  
  
  # add cell barcodes as column names
  names(list_clean_count_tables) <- colnames(object)
  
  
  # return list of count tables
  return(list_clean_count_tables)
  
  
}

# script to make count table as RDS -----------------------------------------------------------------------

# load packages
library(reticulate)
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(SummarizedExperiment)
library(Seurat)

# import pandas to load pickle object
pd <- import("pandas")

# get list of pickle objects
pickle_files <- list.files(opt$input_directory, pattern = "pickle", full.names = T)


# get cell barcodes from file names
cell_barcodes <- gsub(".+tables/(.+).pickle", "\\1", pickle_files)


# import one object to get the column names
pickle <- pd$read_pickle(pickle_files[1])


# load high-quality barcodes
barcodes <- read_csv(opt$barcodes,
                     col_names = F) %>% pull(X1)

# get reads/cell in the mitochondrial library
reads_cell <- read_csv(opt$reads_cell,
                       col_names = F) %>% 
                separate(X1, into = c("reads_mito_lib", "cell_barcode"), sep = " ") %>% 
                mutate(cell_barcode = gsub(".+Z:", "", cell_barcode)) %>% 
                filter(cell_barcode %in% barcodes)


count <- 0

# For each feature measured aggregate the data in a matrix to fit in a summarisedExperiment object 
list_matrix <- mclapply(1:lengths(pickle), function(i){
  
  # iterate through all cells
  list_cell <- lapply(1:length(pickle_files), function(x){
    
    count <- count +1
    
    # load pickle object
    pickle_cell <- pd$read_pickle(pickle_files[x])
    
    # transform into a sparse matrix
    as.sparse(as.data.frame(round(pickle_cell[[cell_barcodes[x]]][[i]])))
    

  })
  
  # create a matrix with cells as columns and mitochondrial positions as rows
  matrix <- do.call("cbind", list_cell)
  
  colnames(matrix) <- cell_barcodes
  
  rownames(matrix) <- as.character(1:nrow(matrix))
  
  matrix
  
}, mc.cores = opt$cores)

# add assay names as list names
names(list_matrix) <- names(pickle[[1]])


# make rownames as positions in the MT genome
row_data <-  GRanges(seqnames = "MT", 
                    ranges = IRanges(1:nrow(list_matrix[[1]]), width = 1))


# add cell barcodes as column names and library size as metadata
col_data <- data.frame(cell_barcode = cell_barcodes, row.names = cell_barcodes) %>% 
              left_join(reads_cell)


# crete Summarised experiment object
sum_object <- SummarizedExperiment(assays = list(
  
              "A_counts" = list_matrix[["A_counts"]], "T_counts" = list_matrix[["T_counts"]],
              "C_counts" = list_matrix[["C_counts"]], "G_counts" = list_matrix[["G_counts"]],
              "A_qual" = list_matrix[["A_qual"]], "T_qual" = list_matrix[["T_qual"]],
              "C_qual" = list_matrix[["C_qual"]], "G_qual" = list_matrix[["G_qual"]],
              "reads_umi" = list_matrix[["reads_umi"]], "nUMIs" = list_matrix[["depth"]]
  
), colData = col_data, rowData = row_data)



# save SummarisedExperiment object as RDS
saveRDS(sum_object, file = opt$summarised_experiment)


# make a list of count tables for further variant calling
list_count_tables <- make_mito_count_table(sum_object)


# save count table as RDS for variant calling
saveRDS(list_count_tables, file = opt$outRDS)


# remove temporal files unless otherwise specified
if(opt$keep_temp == F){
  
  # string to remove folder with temporal BAM files
  remove_bams <- paste0("rm -r ", gsub("/count_tables$", "", opt$input_directory))
  
  system(remove_bams)
  
}


