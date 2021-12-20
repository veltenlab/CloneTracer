#!/usr/bin/env Rscript

# Design inner and outer primers for specific amplification of target sites (GoT)


# set the library folder
.libPaths("/nfs/users2/lvelten/sbeneyto/.conda/envs/exome/lib/R/library")


# parse command line arguments ---------------------------------------------------------------------

library(optparse)

# create arguments list
option_list = list(
  make_option(c("-i", "--input_csv"), type = "character", default = NULL,
              help = "csv file with the mutations of interest. It must contain the following 3 columns: 
                      CHROM: chromosome,
                      POS: position,
                      symbol: gene name",
              metavar = "character"),
  make_option(c("-b", "--subsetted_bam"), type = "character", default = NULL,
              help = "BAM file subsetted for the genes of interest", 
              metavar = "character"),
  make_option(c("-g", "--gtf_file"), type = "character", default = NULL,
              help = "GTF file from ensembl", 
              metavar = "character"),  
  make_option(c("-p", "--patient_name"), type = "character", default = NULL,
              help = "Name of the patient sample as specified in the Snakemake pipeline", 
              metavar = "character"),  
  make_option(c("-t", "--cores"), type = "integer", default = 8,
              help = "Number of cores", metavar = "character"),
  make_option(c("-d", "--out_directory"), type = "character", default = NULL,
              help = "directory where bed files will be saved", 
              metavar = "character"),  
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "relative path to the final csv file with variant information and distance to the polyA tail", 
              metavar = "character"),
  make_option(c("-e", "--exome_format"), type = "logical", default = TRUE, 
              help = "boolean indicating whether the input csv file contains all columns from the variant calling
              pipeline. FALSE means that the file contains the columns CHROM, POS and symbol but may miss some of the remaining ones")
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
required_args <- c("input_csv", "subsetted_bam", "patient_name","out_directory", 
                   "cores","gtf_file", "outfile", "exome_format")


for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}



# code to run the script .............................................................................

# load packages
library(tidyverse)
library(Seurat)
library(BiocParallel)
library(rtracklayer)
library(BSgenome)
library(TAPseq)
library(ballgown)
library(purrr)
library(BSgenome.Hsapiens.UCSC.hg38)



# read table with variants
raw_variants <- read_csv(opt$input_csv)


# extract symbols from genes of interest
gene_names <- raw_variants$symbol %>% unique()


# read gene annotation from gtf file to data frame
gtf_genes <- gffReadGR(opt$gtf_file)


# GET GENE EXPRESSION ---------------------------------------------------------------------

# check whether BAM file from the sample of interest is exists and create a Seurat object
if(file.exists(paste0("../gene_expression/data/", opt$patient, "/outs/filtered_feature_bc_matrix/matrix.mtx.gz"))){
  
  # create seurat object
  count_data <- Read10X(data.dir = paste0("../gene_expression/data/", opt$patient, "/outs/filtered_feature_bc_matrix/"))
  count_matrix <- CreateSeuratObject(counts = count_data$`Gene Expression`)
  
}else{
  
  # load Seurat object with merged count tables from Anne's patients
  count_matrix <- readRDS("data/reference_data/seurat_object_anne_raw.rds")
  
}

# compute total number of cells
total_cells <- ncol(GetAssayData(count_matrix, slot = "counts"))


# filter genes which are present in the counts matrix
gene_names <- gene_names[which(gene_names %in% rownames(count_matrix))]


# get raw counts for genes of interest
counts <- count_matrix@assays$RNA@counts[gene_names,] %>% 
             rowSums()


# divide counts by the total number of cells
gene_expression <- data.frame(symbol = names(counts),
                              counts_cell = counts/total_cells)


# DISTANCE TO polyA ---------------------------------------------------------------------


# get gene names
genes_mut_exons <- raw_variants$symbol %>% unique()


# put mutations in GRanges object
mutations <- unlist(GRangesList(lapply(1:nrow(raw_variants), function(i){
  
  range <- GRanges(seqnames = raw_variants$CHROM[i], ranges = raw_variants$POS[i])
  
})))



# get exons from genes of interest
target_genes <- gtf_genes[mcols(gtf_genes)[,"type"] == "exon" &
                            mcols(gtf_genes)[,"gene_name"] %in% genes_mut_exons] %>% 
                  sort()


# make seqlevels style equal
seqlevelsStyle(mutations) <- seqlevelsStyle(target_genes) <- "UCSC"


# get exons which overlap with the mutations of interest
exons_overlap <- to(findOverlaps(mutations, target_genes))


# extract overlapping transcripts
transcripts_overlap <- elementMetadata(target_genes[exons_overlap][,"transcript_id"]) %>%
                          unique() %>% as.data.frame() %>% pull()


# filter out transcripts which do not overlap with the mutations of interest.
target_genes <- target_genes[elementMetadata(target_genes)[,"transcript_id"] %in% transcripts_overlap]



# remove exons longer than 1kb. The 1st and last segments are excluded as they
# corresponded to the 3'-UTR and 5'-UTR 

filtered_exons <- mclapply(1:length(transcripts_overlap), function(i){
  
  # get exons from a particular transcript
  exons <- target_genes[elementMetadata(target_genes)[,"transcript_id"] == transcripts_overlap[i]]
  
  
  # check if any of the exons is longer than 1kb
  if(length(exons) == length(unique(c(1, which(width(ranges(exons))<1000), length(exons))))){
    
    filtered_ranges <- exons
    
    exons_excluded <- FALSE
    
    # if large exons are present and do not overlap with the mutation they are eliminated
  }else{exons_excluded <- TRUE
  
  
  # determine if any of the mutations fall into the exons
  overlap <- to(findOverlaps(mutations, exons))
  
  
  # filter out long exons
  filtered_ranges <- exons[unique(sort(c(1, which(width(ranges(exons))<1000), overlap, length(exons))))]
  
  
  }
  
  
  return(list(exons_filtered = exons_excluded, ranges = filtered_ranges))
  
}, mc.cores = opt$cores)


# merge filtered exons in single GRanges object
target_exons <- unlist(GRangesList(map(filtered_exons,2)))


# split GRanges exons by gene name
target_exons <- split(target_exons, f = target_exons$gene_name)


# extract genes for which exons were excluded in order to later on flag the genes
exon_flags <- target_genes[elementMetadata(target_genes)[,"transcript_id"] %in%
                             transcripts_overlap[which(unlist(map(filtered_exons, 1)))]] %>% 
                            .$gene_name %>% unique()


# set the number of cores to infer polyA site. 
#I recommend adjusting the number of workers based on the number of CPUs of the computer
register(MulticoreParam(workers = opt$cores))


# infer polyA sites. With the subsetted BAM it takes a few seconds
polyA_sites <- inferPolyASites(target_exons,
                               bam = opt$subsetted_bam, 
                               polyA_downstream = 50, by = 1,
                               wdsize = 100, min_cvrg = 100, parallel = TRUE)


# add gene names as metadata column to the GRanges object
polyA_sites$gene_name <- names(polyA_sites)


# split polyAs by gene
polyA_sites <- split(polyA_sites, f = polyA_sites$gene_name)


# function to select the smallest positive value
min_pos <- function(x){min(x[x>0])}


# select a unique polyA per gene
polyAs <- lapply(1:length(polyA_sites), function(i){
  
  
  gene_ends <- polyA_sites[[i]]
  
  
  # set polyA score threshold as half of the maximum score by the top ranked polyA
  threshold <- max(gene_ends$score)/2
  
  
  # set a minimum threshold of 400 
  #if (threshold < 400){threshold <- 400}
  
  # when only one polyA has a score a above the threshold
  if(length(which(elementMetadata(gene_ends)[,"score"] > threshold)) == 1){
    
    
    tail_candidates <- gene_ends[which(elementMetadata(gene_ends)[,"score"] > threshold)]
    
    mutation_position <- raw_variants %>% filter(symbol %in% names(gene_ends)) %>% 
      pull(POS) %>% as.integer()
    
    polyA_flag <- FALSE
    
    # the strand is important when it comes to computing distances
    if(strand(tail_candidates)@values == "+"){
      
      
      distances <- ranges(tail_candidates)@start-mutation_position
      
      
      # if the inferred polyA  upstream of the mutation a flag is raised
      if(min_pos(distances) == Inf){
        
        
        selected_tail <- NULL
        
        # if the polyA is downstream then is selected
      }else{selected_tail <- tail_candidates}
      
    }else{
      
      distances <- (ranges(tail_candidates)@start-mutation_position)*(-1)
      
      # if all inferred polyA is upstream of the mutation a flag is raised
      if(min_pos(distances) == Inf){
        
        selected_tail <- NULL
        
        # if the polyA is downstream then is selected
      }else{selected_tail <- tail_candidates}      
    }
    
    
    # if 2 or more polyAs have a score higher than 50 I pick the one closer to the mutation site.
  }else if(length(which(elementMetadata(gene_ends)[,"score"] > threshold)) > 1){
    
    tail_candidates <- gene_ends[which(elementMetadata(gene_ends)[,"score"] > threshold)]
    
    
    mutation_position <- raw_variants %>% filter(symbol %in% names(gene_ends)) %>% 
      pull(POS) %>% as.integer()
    
    polyA_flag <- TRUE
    
    # if there are more than 1 mutation in a particular gene I take the most downstream mutation
    if(length(mutation_position) > 1){
      
      
      # the strand is important when it comes to computing distances
      if(strand(tail_candidates)@values == "+"){
        
        # take the most downstream mutation
        mutation_position <- max(mutation_position)
        
        
        distances <- ranges(tail_candidates)@start-mutation_position
        
        
        # if all the inferred polyAs are upstream of the mutations a flag is raised
        if(min_pos(distances) == Inf){
          
          polyA_flag <- FALSE
          
          
          selected_tail <- NULL
          
          # if the polyA is downstream the is selected
        }else{selected_tail <- tail_candidates[which(distances %in% min_pos(distances))]}
        
        
      }else{
        # the most downstream mutation is in this case the one with a lower position due to the - strand
        mutation_position <- min(mutation_position)
        
        
        distances <- (ranges(tail_candidates)@start-mutation_position)*(-1)
        
        
        # if all the inferred polyAs are upstream of the mutations a flag is raised
        if(min_pos(distances) == Inf){
          
          polyA_flag <- FALSE
          
          
          selected_tail <- NULL
          
          # if the polyA is downstream the is selected
        }else{selected_tail <- tail_candidates[which(distances %in% min_pos(distances))]}
        
      } 
      
      
      # when there is only one mutation/gene, it gets automatically selected.  
    }else{
      
      # the strand is important when it comes to computing distances
      if(strand(tail_candidates)@values == "+"){
        
        distances <- ranges(tail_candidates)@start-mutation_position
        
        # if all the inferred polyAs are upstream of the mutations a flag is raised
        if(min_pos(distances) == Inf){
          
          polyA_flag <- FALSE
          
          
          selected_tail <- NULL
          
          # if the polyA is downstream the is selected
        }else{selected_tail <- tail_candidates[which(distances %in% min_pos(distances))]}
        
      }else{
        
        distances <- (ranges(tail_candidates)@start-mutation_position)*(-1)
        
        # if all the inferred polyAs are upstream of the mutation a flag is raised
        if(min_pos(distances) == Inf){
          
          polyA_flag <- FALSE
          
          
          selected_tail <- NULL
          
          # if the polyA is downstream the is selected
        }else{selected_tail <- tail_candidates[which(distances %in% min_pos(distances))]}                
      } 
      
    }
    
  }else{return(NULL)}
  
  
  return(list(ranges = selected_tail, flags = polyA_flag))
  
})



# select genes with a polyA detected
filt_polyAs <- polyAs[which(lengths(map(polyAs, 1)) == 1)]


# put all polyAs in one GRanges object
final_polyAs <-  unlist(GRangesList(map(filt_polyAs,1)))



# genes for which no polyA was found
genes_no_polyA <- genes_mut_exons[!genes_mut_exons %in% final_polyAs$gene_name]


# truncate transcripts at the inferred polyA sites
target_transcripts <- GenomicRanges::sort(GenomicRanges::reduce(truncateTxsPolyA(target_exons, 
                                                                                 polyA_sites = final_polyAs, 
                                                                                 parallel = TRUE,
                                                                                 transcript_id = "transcript_id"))) 


# get genes with multiple_polyA flag
polyAs_flag <- names(polyA_sites)[which(unlist(map(polyAs,2)))]


print("For the following genes more than one polyA was detected:")
print(polyAs_flag)


print("For the following genes no polyA was detected:")
print(genes_no_polyA)


# EXPORT BED FILES ----------------------------------------------------------------------

# create directory where bed files will be stored

if(file.exists(opt$out_directory) == F){
  
  dir.create(opt$out_directory)
  
}

## It is recommended to check the PolyA tails in IGV manually. 


# add gene names to the polyA object (only name and score columns are exported to the final bed file)
final_polyAs$name <- final_polyAs$gene_name


# export polyAs to bed file to manually check in IGV
export(unlist(polyA_sites), con = paste0(opt$out_directory, "polyA_candidates.bed"), format = "bed")


# export target gene annotations to load in IGV
export(unlist(target_exons), con = paste0(opt$out_directory, "GoT_hits.gtf"), format = "gtf")


# export mutations to BED file in order to load in IGV
mut_bed <- GRanges(seqnames = raw_variants$CHROM,
                   ranges = IRanges(start = as.integer(raw_variants$POS), 
                                    end = as.integer(raw_variants$POS)))
export(mut_bed, con = paste0(opt$out_directory, "mutations.bed"), format = "bed")


# export final selection of polyA sites
export(final_polyAs, con = paste0(opt$out_directory, "polyA_final_selection.bed"), format = "bed")


# COMPUTE DISTANCE TO THE polyA tail ---------------------------------------------------------


# locate the mutation among the exons
list_regions <- lapply(1:length(raw_variants$symbol), function(i){
  
  
  # if the gene is not found in the gtf file then no polyA distance can be computed
  if(raw_variants$symbol[[i]] %in% names(target_transcripts)){
    
    # Get exons from the gene of interest and merge overlapping sequences.
    exons <- target_transcripts[[raw_variants$symbol[i]]]
    seqlevelsStyle(exons) <- "UCSC"
    
    
    # Find where the mutation of interest is in the exon ranges
    mutation_site <- GRanges(seqnames = raw_variants$CHROM[i], ranges = raw_variants$POS[i])
    seqlevelsStyle(mutation_site) <- "UCSC"
    overlap <- to(findOverlaps(mutation_site, exons))
    
    
    # create a new GRanges object with exomic information upstream of the mutation of interest
    # the strand information is important in order to generate the sequence template
    # if the mutation is not mapped into any exonic region it should be reported
    
    if(length(overlap) == 1){
      
      if(runValue(strand(exons[1])) == "-"){
        
        
        if(overlap == 1){
          
          
          filtered_region <- GRanges(seqnames = runValue(seqnames(exons)), 
                                     ranges = IRanges(start = start(exons[overlap]), 
                                                      end = start(mutation_site)),
                                     strand = as.character(strand(exons[1])@values))
          
        }else{
          
          filtered_region <- c(exons[1:(overlap-1)],
                               GRanges(seqnames = runValue(seqnames(exons)), 
                                       ranges = IRanges(start = start(exons[overlap]), 
                                                        end = start(mutation_site)),
                                       strand = as.character(strand(exons[1])@values)))
        }
        
      }else{
        if(overlap == length(exons)){
          
          filtered_region <- GRanges(seqnames = runValue(seqnames(exons)), 
                                     ranges = IRanges(start = start(mutation_site), 
                                                      end = end(exons[overlap])),
                                     strand = as.character(strand(exons[1])@values))
          
        }else{
          
          filtered_region <- c(GRanges(seqnames = runValue(seqnames(exons)), 
                                       ranges = IRanges(start = start(mutation_site), 
                                                        end = end(exons[overlap])),
                                       strand = as.character(strand(exons[1])@values)),
                               exons[(overlap+1):length(exons)])
        }
      }
      
    }else{message(paste0(raw_variants$symbol[i], " ", i,
                         ": Mutation does not overlap with exonic regions! Check polyA sites in IGV."))}
    
    
  }else{filtered_region <- GRanges()}
  
  
})


# add gene names 
names(list_regions) <- paste(raw_variants$symbol, raw_variants$CHROM, raw_variants$POS, 
                             sep = "_")


# remove genes for which no overlap between exons and mutation was found
list_regions <- list_regions[lengths(list_regions) != 0]


# get sequence from the generated regions
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
txs_seqs <- getTxsSeq(GRangesList(list_regions), genome = hg38)


# extract distance to the 3'-end for each gene
distance_table <- tibble(symbol = names(txs_seqs),
                         distance_3_end = width(txs_seqs)) %>% 
                    separate(symbol, into = c("symbol", "CHROM","POS"), sep = "_") %>% 
                    mutate(POS = as.integer(POS))



# ANNOTATED TABLE --------------------------------------------------------------

annotated_variants <- raw_variants %>% 
                            mutate(POS = as.integer(POS)) %>% 
                            left_join(distance_table) %>% 
                            left_join(gene_expression) 


# add flags column
flag_vector <- sapply(1:length(unique(raw_variants$symbol)), function(x){
  
  x <- unique(raw_variants$symbol)[x]
  
  flag_vector <- vector()
  
  if(x %in% exon_flags){
    
    
    flag_vector <- "long_exon_excluded"
    
  }
  
  if(!x %in% names(final_polyAs)){
    
    flag_vector <- c(flag_vector, "no_polyA_found")
  }
  
  if(x %in% polyAs_flag){
    
    flag_vector <- c(flag_vector, "multiple_polyAs")
    
  }
  
  return(paste(flag_vector, collapse = ";"))
  
})

flag_table <- data.frame(flag = flag_vector,
                         symbol = unique(raw_variants$symbol))

# order columns in the table depending on the format
if(opt$exome_format){
  
  
  # add flags to the table
  annotated_variants <- annotated_variants %>% left_join(flag_table) %>% 
                          mutate(flag = if_else(counts_cell < 0.1, 
                                                paste(flag, "low_expressed", sep = ";"), flag),
                                 flag = if_else(flag == "", "PASS", flag)) %>% 
                          select(symbol, CHROM, POS, REF, ALT, consequence,
                                 counts_cell, distance_3_end, af_tcells, af_tumor, depth_tcells,depth_tumor_cells,
                                 prob_af_tcells, prob_af_tumor, flag,
                                 cancer_gene, cosmic_id, FILTER, coverage) %>% 
                          arrange(distance_3_end)
  
}else{
  
  # add flags to the table
  annotated_variants <- annotated_variants %>% left_join(flag_table) %>% 
                          mutate(flag = if_else(counts_cell < 0.1, 
                                                paste(flag, "low_expressed", sep = ";"), flag),
                                 flag = if_else(flag == "", "PASS", flag)) %>% 
                          select(symbol, CHROM, POS,
                                 counts_cell, distance_3_end, flag,everything()) %>% 
                          arrange(distance_3_end)
  
}

# save the annotated table
write_csv(annotated_variants, file = opt$outfile)
