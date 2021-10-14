##############################   qc_exome.R   ############################## 
#' Generate cumulative coverage of target regions for exome sequencing data. 

## parse command line arguments ---------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

cat(date(),": Start creating QC plots.\n")
option_list <- list(
    make_option(c("-r", "--reference"), type="character", default=NULL,
        help= "Path to the QC table from the healthy reference sample generated the picard tool CollectHsMertics",
        metavar="character"),
    make_option(c("-t", "--tumor"), type="character", default=NULL,
        help="Path to the QC table from the tumor sample generated the picard tool CollectHsMertics", 
        metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="Path to the png file which will store the target coverage plot",
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

## check that all required parameters are provided
required_args <- c("reference", "tumor", "output")
res <- lapply(required_args, check_required_args, opt=opt, opt_parser=opt_parser)

####################################.............................................................................
# read QC summary file
qc_metrics <- str_split(read_lines(opt$tumor)[8], "\t")[[1]]  %>% 
    as_tibble() %>% 
    mutate(colnames=str_split(read_lines(opt$tumor)[7], "\t")[[1]]) %>% 
    filter(colnames %in% c("PCT_TARGET_BASES_1X", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", 
        "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_40X", 
        "PCT_TARGET_BASES_50X", "PCT_TARGET_BASES_100X")) %>% 
    mutate(sample="tumor_cells", value=as.double(value)) %>% 
    bind_rows(str_split(read_lines(opt$reference)[8], "\t")[[1]]  %>% 
    as_tibble() %>% 
    mutate(colnames=str_split(read_lines(opt$reference)[7], "\t")[[1]]) %>% 
    filter(colnames %in% c("PCT_TARGET_BASES_1X", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", 
        "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_40X", 
        "PCT_TARGET_BASES_50X", "PCT_TARGET_BASES_100X")) %>% 
    mutate(sample="Tcells", value=as.double(value))) %>% 
    mutate(colnames=factor(gsub(".*BASES\\_", "", colnames), 
        levels=c("1X", "2X", "10X","20X", "30X", "40X", "50X", "100X")),
        sample=as.factor(sample)) %>% 
    rename(coverage=colnames)

target_coverage <- ggplot(qc_metrics, aes(x=coverage, y=value, colour=sample, group=sample)) + 
    geom_point(aes(colour=sample)) + 
    geom_line() +
    theme_classic() +
    scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 1), expand=c(0,0)) +
    ylab("Fraction of targets covered")

ggsave(target_coverage, filename=opt$output)

cat(date(),": QC plot generated successfully!\n")
