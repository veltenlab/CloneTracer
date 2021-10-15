##############################   filter_variants.R   ############################## 
#' Filter output from an annotated vcf file. 
#' 
#' This file filters for variants which are present in exonic regions and are non-synonymous. Furthermore, we compute 
#' the difference between the allele frequency in the tumor and control cells and use it as aparameter to select 
#' interesting variants.
#' This approach works better than merely accepting the variants labelled as 'PASS' by Mutect2 particularly if 
#' T cells are used as controls.
#' Additionally if a variant does not fullfil this criteria, but affects an AML gene present in the myeloid 
#' Illumina panel, it is also added to the final table. This can be useful to identify pre-leukemic mutations, 
#' for example if T cells are used as controls.
#' Synonymous mutations, which occur in highly expressed genes (using Anne's gene expression data) and have a high 
#' difference in allele frequency between the tumor and the T cells, are also filtered for. These, if present, could
#' be used as markers (similar to mitochondrial mutations).
# parse command line arguments ---------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(ballgown))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Seurat))

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
        help = "vcf file annotated with annovar as outputed by the Mutect2 exome pipeline",
        metavar = "character"),
    make_option(c("-g", "--gtf_file"), type = "character", default = NULL,
        help = "gtf file from which the exon sequences will be extracted", 
        metavar = "character"),
    make_option(c("-o", "--outcsv"), type = "character", default = NULL,
        help = "relative path where the variant table will be saved (must be a .csv file)", 
        metavar = "character"),
    make_option(c("-b", "--outbed"), type = "character", default = NULL,
        help = "BED file with the genomic coordinates of the selected variants", 
        metavar = "character"),
    make_option(c("-s", "--seurat"), type = "character", default = NULL,
                help = "Path to a Seurat (RDS) file to check for gene expression", 
                metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
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
required_args <- c("input", "gtf_file","outcsv", "outbed")
res <- lapply(required_args, check_required_args, opt=opt, opt_parser=opt_parser)

## SCRIPT -------------------------------------------------------------------
## read file VCF to dataframe. Extract allele frequency and approximated read depth
vcf <- vcfR2tidy(read.vcfR(opt$input),
    format_fields = c("AF", "DP", "AD"))
variants <- vcf$fix  
          
## get allele frequencies in tcells and tumor cells for the sites of interest          
af <- vcf$gt %>% filter(ChromKey %in% variants$ChromKey & POS %in% variants$POS) %>% 
    pivot_wider(id_cols = c(ChromKey, POS), names_from = Indiv,
        values_from = gt_AF) %>% 
    mutate(prob_af_cd3_cells = as.double(cd3_cells),
        prob_af_tumor_cells = as.double(tumor_cells)) %>% 
    dplyr::select(-tumor_cells, -cd3_cells)

## get approximate read depth in the positions
coverage <- vcf$gt %>% filter(ChromKey %in% variants$ChromKey & POS %in% variants$POS) %>% 
    pivot_wider(id_cols = c(ChromKey, POS), names_from = Indiv,
        values_from = gt_DP) %>% 
    mutate(coverage = as.double(cd3_cells) + as.double(tumor_cells)) %>% 
    dplyr::rename(depth_tcells = cd3_cells, depth_tumor = tumor_cells) 

## get reference and mutant counts
read_af <- vcf$gt %>% 
    separate(gt_AD, into = c("ref_counts", "alt_counts"), sep = ",") %>% 
    mutate(af = as.double(alt_counts)/(as.double(alt_counts)+as.double(ref_counts))) %>% 
    filter(ChromKey %in% variants$ChromKey & POS %in% variants$POS) %>% 
    pivot_wider(id_cols = c(ChromKey, POS), names_from = Indiv,
        values_from = af) %>% 
    dplyr::rename(af_tcells = cd3_cells, af_tumor = tumor_cells)

## AML GENES ---------------------------------------------------------------------------------------
## get variants in AML genes present in the Illumina myeloid panel
aml_genes <- read_delim("../data/myeloid_panel_genes.txt", delim = "\n",
    col_names = FALSE) %>% pull(X1)

## get aml variants if present
if(any(aml_genes %in% variants$Gene.refGene)){
aml_variants <- variants %>%  filter(Func.refGene %in% c("exonic", "splicing"),
        ExonicFunc.refGene %in% c("nonsynonymous_SNV", "frameshift_insertion",
            "frameshift_deletion", "stopgain")) %>% 
    filter(FILTER %in% c("PASS", "germline", "normal_artifact", "germline;normal_artifact")) %>% 
    filter(Gene.refGene %in% aml_genes) %>%  
    left_join(af) %>%
    left_join(coverage) %>% 
    left_join(read_af) %>% 
    mutate(af_diff = (af_tumor - af_tcells)/af_tcells,
        cancer_gene = ifelse(cosmic == ".", "no", "yes"),
        cosmic_id = gsub(".+(COSV.+)\\\\x3bO.+", "\\1", cosmic)) %>% 
    dplyr::select(-TLOD, -cosmic) %>% 
    filter(af_diff > 0.1)
}

## SYNONYMOUS VARIANTS -----------------------------------------------------------------------
## filter all non-exonic or synonymous mutations
syn_variants <- variants %>% filter(Func.refGene %in% c("exonic", "splicing"),
        !ExonicFunc.refGene %in% c(".", "unknown", "startloss","nonsynonymous_SNV", 
            "frameshift_insertion", "frameshift_deletion", "stopgain")) %>% 
    filter(FILTER %in% c("PASS", "germline", "normal_artifact", "germline;normal_artifact")) 

if(!is.null(opt$seurat)){
    ## load Seurat object with merged count tables
    count_matrix <- readRDS(opt$seurat)
    ## compute total number of cells
    total_cells <- ncol(GetAssayData(count_matrix, slot = "counts"))
    ## get gene names for which variants are found
    gene_names <- syn_variants %>% pull(Gene.refGene) %>% unique() %>% as.character()
    ## get only gene names which are present in the Seurat object
    gene_names <- gene_names[which(gene_names %in% rownames(GetAssayData(count_matrix, slot = "counts")))]
    ## sum counts for the genes of interest and divide by the total number of cells
    high_genes <- GetAssayData(count_matrix, slot = "counts")[gene_names,] %>% 
        as.data.frame() %>% t() %>% as.data.frame() %>%  
        summarise_if(is.numeric, sum, na.rm = TRUE) %>% 
        pivot_longer(cols = everything(), 
            names_to = "gene", values_to = "raw_counts") %>% 
        mutate(counts_cell = round(raw_counts/total_cells, digits = 2)) %>% 
        dplyr::select(-raw_counts) %>% 
        dplyr::rename(symbol = gene) %>% 
        filter(counts_cell > 0.5) %>% 
        pull(symbol)
    ## filter synonymous variants that happen in low expressed genes
    syn_variants <- syn_variants %>% filter(Gene.refGene %in% high_genes)
}

## filter variants with low allele frequency in the tumor
final_syn_variants <- syn_variants %>% left_join(af) %>% left_join(coverage) %>% 
    left_join(read_af) %>% 
    mutate(af_diff = (af_tumor - af_tcells)/af_tcells) %>% 
    mutate(cancer_gene = ifelse(cosmic == ".", "no", "yes"),
        cosmic_id = gsub(".+(COSV.+)\\\\x3bO.+", "\\1", cosmic)) %>% 
    dplyr::select(-TLOD, -cosmic) %>% 
    filter(af_diff > 2 & af_tumor > 0.1 & depth_tcells > 10 | FILTER == "PASS" & depth_tcells > 10)

## FINAL TABLE ---------------------------------------------------------------------------------
#' select final variants based on the difference between AF in the tumor and Tcells or in the 'PASS' filter
#' with sufficient coverage
final_variants <- variants %>% 
    filter(Func.refGene %in% c("exonic", "splicing"),
        ExonicFunc.refGene %in% c("nonsynonymous_SNV", "frameshift_insertion",
            "frameshift_deletion", "stopgain")) %>% 
    filter(FILTER %in% c("PASS", "germline", "normal_artifact", "germline;normal_artifact")) %>%
    left_join(af) %>% 
    left_join(coverage) %>% 
    left_join(read_af) %>% 
    mutate(af_diff = (af_tumor - af_tcells)/af_tcells) %>% 
    filter((af_diff > 2 & af_tumor > 0.1 & depth_tcells > 10) | (FILTER == "PASS" & depth_tcells > 10)) %>% 
    mutate(cancer_gene = ifelse(cosmic == ".", "no", "yes"),
        cosmic_id = gsub(".+(COSV.+)\\\\x3bO.+", "\\1", cosmic)) %>% 
    dplyr::select(-TLOD, -cosmic) %>% 
    rbind(aml_variants) %>% 
    rbind(final_syn_variants) %>% 
    distinct() %>% 
    dplyr::select(Gene.refGene, CHROM, POS, REF, ALT, ExonicFunc.refGene, af_tcells, af_tumor, 
        depth_tcells, depth_tumor, prob_af_cd3_cells, prob_af_tumor_cells,
        cancer_gene, cosmic_id, FILTER, coverage, af_diff) %>% 
    arrange(desc(af_diff))

## set colnames
colnames(final_variants) <- c("symbol", "CHROM", "POS", "REF", "ALT", 
    "consequence", "af_tcells","af_tumor", "depth_tcells","depth_tumor_cells",
    "prob_af_tcells", "prob_af_tumor", "cancer_gene", "cosmic_id","FILTER", "coverage", "af_diff")

## write final table to csv file
write_csv(final_variants, file = opt$outcsv)

## read gtf file
transcripts <- gffReadGR(opt$gtf)

## get transcript coordinates for the genes in the variant table
filtered_transcripts <- transcripts[elementMetadata(transcripts)[,"gene_name"] %in% final_variants$symbol]

## change seqlevels to UCSC
seqlevelsStyle(filtered_transcripts) <- "UCSC"

## eliminate metadata column in order to save the coordinates into a BED file
mcols(filtered_transcripts) <- NULL

## export gene coordinates as BED file
export(filtered_transcripts, con = opt$outbed, format = "bed")