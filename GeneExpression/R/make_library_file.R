##############################   make_library_file.R   ############################## 
#' This script creates the libraries.csv file which specifies where are the RNA and ADT fastqs for cellranger

# parse command line arguments ---------------------------------------------------------------------
library(optparse)
library(tidyverse)
# create arguments list
option_list=list(
    make_option(c("-p", "--patient_name"), type="character", default=NULL,
        help="Name of the patient as specified in the Snakemake pipeline",
        metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="Path were output file will be stored",
        metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser){
    if (is.null(opt[[arg]])){
        print_help(opt_parser)
        stop(arg," argument is required!", call.=FALSE)
    }
}
## check that all required parameters are provided
required_args <- c("patient_name")
res <- lapply(required_args, check_required_args, opt=opt, opt_parser=opt_parser)

name <- opt$patient_name
# make paths to fastqs
path_wta <- file.path(getwd(), "raw_data", name, "wta")
path_adt <- file.path(getwd(), "raw_data", name, "ab_citeseq")
# make data frame 
table <- data.frame(fastqs=c(path_wta, path_adt),
                    sample=name,
                    library_type=c("Gene Expression", "Antibody Capture"))
# write csv
write_csv(table, file=opt$output)
