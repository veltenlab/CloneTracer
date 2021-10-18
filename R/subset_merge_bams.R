##############################   filter_variants.R   ############################## 
#' This script subsets BAM files for the filtered variants

## parse command line arguments ---------------------------------------------------------------------

library(optparse)
library(tidyverse)

# create arguments list
option_list=list(
    make_option(c("-b", "--bed"), type="character", default=NULL,
        help="BED file with the genomic coordinates of genes of interest", 
        metavar="character"),
    make_option(c("-p", "--patient_name"), type="character", default=NULL,
        help="patient name from which BAM files will be used in case they are already processed", 
        metavar="character"),
    make_option(c("-d", "--out_directory"), type="character", default=NULL,
        help="directory name where the BAM file will be stored", 
        metavar="character"),  
    make_option(c("-o", "--out_bam"), type="character", default=NULL,
        help="relative path to the subset and merged BAM file", 
        metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#' check_required_args
#' 
#' A function to check for required arguments
check_required_args <- function(arg, opt, opt_parser){
    if (is.null(opt[[arg]])){
        print_help(opt_parser)
        stop(arg," argument is required!", call.=FALSE)
    }
}
required_args <- c("patient_name", "out_bam", "bed", "out_directory")
res <- lapply(required_args, check_required_args, opt=opt, opt_parser=opt_parser)

## Script body .............................................................................
patient_name <- opt$patient_name

## check whether BAM file from the patient of interest is present
if(file.exists(paste0("../gene_expression/data/", patient_name, "/outs/possorted_genome_bam.bam"))){
    file <- paste0("../gene_expression/data/", patient_name, "/outs/possorted_genome_bam.bam")
    ## subset BAM file for the genes of interest
    system(paste("samtools view -h -b -L", opt$bed, file, ">", opt$out_bam))
    ## index subset BAM file
    system(paste("samtools index -b", opt$out_bam))
}else{
    sample_names <- c("AKLA1b", "AKLW1", "AKLK1", "AKLM1")
    ## get BAM files for the specified patients
    bams <- sapply(1:length(sample_names), function(i){
        paste0("../gene_expression/data/", sample_names[i], "/outs/possorted_genome_bam.bam")
    })
    ## subset bam files 
    for (file in bams){
        patient_name <- str_split(file, pattern="/")[[1]][4]
        system(paste("samtools view -h -b -L", opt$bed, file, ">", paste0(opt$out_directory, patient_name, ".bam")))
    }
    ## collect all subsetted bam files
    subset_bams <- list.files(opt$out_directory, pattern=".bam", full.names=T)
    ## merge BAM files
    system(paste("samtools merge", opt$out_bam, subset_bams[1], subset_bams[2], subset_bams[3], subset_bams[4]))
    ## create index file
    system(paste("samtools index -b", opt$out_bam))
    ## get individual bams to eliminate them
    small_bams <- list.files(path=opt$out_directory, pattern="AKL", full.names=T)
    ## remove individual bam files
    system(paste("rm", small_bams[1], small_bams[2], small_bams[3], small_bams[4]))
}