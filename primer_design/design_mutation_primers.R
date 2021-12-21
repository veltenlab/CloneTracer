#!/usr/bin/env Rscript

# Design outer, middle and inner(staggered) primers for specific amplification of mutated sites (GoT)


# set the library folder
.libPaths("/nfs/users2/lvelten/sbeneyto/.conda/envs/rmutprimers/lib/R/library")


# parse command line arguments ---------------------------------------------------------------------

library("optparse")

# create arguments list
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the table containing the genomic coordinates of the sites of interest (one site per row). 
              It should contain four columns: 
                CHROM: chromosome.
                POS: genomic position of the mutation 
                symbol: gene name,
                distance_3_end: distance from the mutation to the 3-end of transcript",
              metavar = "character"),
  make_option(c("-o", "--outprefix"), type = "character", default = NULL,
              help = "Prefix for the output files. If a different directory is desired it should be specified 
              A csv and a BED file with primer information will be generated", metavar = "character"),
  make_option(c("-r", "--read_length"), type = "integer", default = NULL,
              help = "Read length of Read 2. It is recommended to be at least 75 bp since lower numbers could 
              make it difficult to generate specific primers due to the short template sequence", 
                      metavar = "character"),
  make_option(c("-g", "--gtf_file"), type = "character", default = NULL,
              help = "gtf file from which the exon sequences will be extracted", 
              metavar = "character"),
  make_option(c("-b", "--bed_file"), type = "character", default = NULL,
              help = "BED file with the genomic coordinates of the polyAs selected", 
              metavar = "character"),
  make_option(c("-p", "--patient_name"), type = "character", default = NULL,
              help = "Patient sample name to include in the final primer table", 
              metavar = "character"),
  make_option(c("-l", "--gene_list"), type = "character", default = NULL,
              help = "txt file with symbols of the selected genes. One symbol per line", 
              metavar = "character")
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
required_args <- c("input", "outprefix", "read_length", "gtf_file", "bed_file", "patient_name", "gene_list")


for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}



# functions to run the script .............................................................................

library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(BiocParallel)
library(TAPseq)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome)
library(mygene)
library(ballgown)


#function to generate template sequences based on a particular genomic site
get_GoT_regions <- function(genomic_sites, contig_list, gene_names, gtf_file, polyAs){
  
  
  # get transcripts
  gtf_genes <- gffReadGR(gtf_file)
  
  
  # get exons from genes of interest
  target_exons <- gtf_genes[mcols(gtf_genes)[,"type"] == "exon" &
                              mcols(gtf_genes)[,"gene_name"] %in% gene_names] %>% 
                    sort()
  
  # put them in a list split by gene name
  exons <- split(target_exons, f = target_exons$gene_name)
  
  
  # make contig names equal
  seqlevelsStyle(exons) <- "UCSC"
  
  
  # add the 'chr' string to contig list so that seqnames overlap between mutation and exons
  #contig_list <- paste0("chr", contig_list)
  
  
  # read polyA sites
  polyAs <- import(polyAs)[elementMetadata(import(polyAs))[,"name"] %in% gene_names]
  

  list_regions <- lapply(1:length(gene_names), function(i){
    
    
    # get exons from the gene of interest and merge overlapping sequences.
    gene_exons <- exons[[gene_names[i]]]
    

    # convert mutations site to GRanges object
    mutation_site <- GRanges(seqnames = contig_list[i], ranges = genomic_sites[i])

    
    # make seqlevels homogeneous
    seqlevelsStyle(gene_exons) <- seqlevelsStyle(mutation_site) <- "UCSC"
    
    
    # find which transcripts overlap with the mutation of interest
    overlap <- to(findOverlaps(mutation_site, gene_exons))
    
    
    # extract overlapping transcripts
    transcripts_overlap <- elementMetadata(gene_exons[overlap][,"transcript_id"]) %>%
                              unique() %>% as.data.frame() %>% pull()
    
    
    # filter out transcripts which do not overlap with the mutations of interest.
    gene_exons <- gene_exons[elementMetadata(gene_exons)[,"transcript_id"] %in% transcripts_overlap]
    
    
    # I annotated some transcript isoforms with non-canonical exons that I have been observing in some
    # genes. These will be eliminated to make sure that the primers do not fall into this unique exons
    blacklisted_transcripts <- readRDS("data/reference_data/blacklisted_transcripts.rds")
    
    
    # filter blacklisted transcripts
    gene_exons <- gene_exons[!elementMetadata(gene_exons)[, "transcript_id"] %in% blacklisted_transcripts]
    
    
    # if a polyA was selected, then only transcripts which overlap with it are chosen
    if(length(which(polyAs$name == gene_names[i])) != 0){
      
      
      gene_polyA <- polyAs[elementMetadata(polyAs)[,"name"] == gene_names[[i]]]
      
      
      # find which transcripts overlap with the selected polyA
      overlap_polyA <- to(findOverlaps(gene_polyA, gene_exons))
      
      
      # extract overlapping transcripts
      transcripts_overlap <- elementMetadata(gene_exons[overlap_polyA][,"transcript_id"]) %>%
                                unique() %>% as.data.frame() %>% pull()
      
      
      # filter out transcripts which do not overlap with the mutations of interest.
      gene_exons <- gene_exons[elementMetadata(gene_exons)[,"transcript_id"] %in% transcripts_overlap]
      
      
    }
    
    
    # merge ovelapping regions 
    gene_exons <- GenomicRanges::sort(GenomicRanges::reduce(gene_exons))
    
    
    # get a unique overlap range between the mutation and the list of exons
    overlap <- to(findOverlaps(mutation_site, gene_exons))
    
    
    # Create a new GRanges object with exomic information upstream of the mutation of interest.
    # The strand information is important in order to generate the sequence template.
    if(as.character(strand(gene_exons[1])@values) == "-"){
      
      # If the mutation is present in the most downstream exon (which for genes in the negative strand is the 1st exon)
      if(overlap == length(gene_exons)){
      
        
        filtered_region <- GenomicRanges::sort(GRanges(seqnames = contig_list[i], 
                                   ranges = IRanges(start = start(mutation_site), end = end(gene_exons[overlap])),
                                   strand = as.character(strand(gene_exons[1])@values)), decreasing = T)
        
      }else{
        
        
        filtered_region <- GenomicRanges::sort(c(GRanges(seqnames = contig_list[i], 
                                     ranges = IRanges(start = start(mutation_site), end = end(gene_exons[overlap])),
                                     strand = as.character(strand(gene_exons[1])@values)),
                             gene_exons[(overlap+1):length(gene_exons)]), decreasing = T)
        
      }
      
    }else{
      
      
      if(overlap == 1){
        
        
        filtered_region <- GenomicRanges::sort(GRanges(seqnames = contig_list[i], 
                                   ranges = IRanges(start = start(gene_exons[overlap]), end = start(mutation_site)),
                                   strand = as.character(strand(gene_exons[1])@values)), decreasing = F)
        
      }else{
        
        
        filtered_region <- GenomicRanges::sort(c(gene_exons[1:(overlap-1)], 
                             GRanges(seqnames = contig_list[i], 
                                     ranges = IRanges(start = start(gene_exons[overlap]), end = start(mutation_site)),
                                     strand = as.character(strand(gene_exons[1])@values))), decreasing = F)
        }
    }
    
  })
  
  return(list_regions)
  
}


# Code to design outer, middle and staggered inner primers close to mutations of interest ----------------

gene_names <- read_delim(opt$gene_list, delim = "\n", col_names = F) %>% 
                  pull(X1) %>% 
                  sub(x = ., pattern = "_\\d+", replacement = "")


gene_positions <- read_delim(opt$gene_list, delim = "\n", col_names = F) %>% 
                    pull(X1) %>% 
                    sub(x = ., pattern = "\\w+_", replacement = "")


# Read input table
hits_table <- read_csv(opt$input) %>% 
                filter(symbol %in% gene_names & POS %in% gene_positions)


# Check that the input table has the correct column namesquite
if(length(which(c("CHROM", "POS", "symbol") %in% colnames(hits_table)))<3){

    
  stop("The input table does not have the correct columns. It must include CHROM, POS and symbol columns. 
         See the help page for more details",
       call. = F)
  
} 


# generate template sequences which include the exons upstream of the genomic site and the 5'UTR
message("Generating template sequences")


# compute regions
regions_list <- get_GoT_regions(genomic_sites = hits_table$POS, contig_list = hits_table$CHROM, 
                                gene_names = hits_table$symbol, gtf_file = opt$gtf_file,
                                polyAs = opt$bed_file)


# put sequences into a GRangesList
ranges_list <- GRangesList(regions_list)


# name each item with gene name
names(ranges_list) <- hits_table$symbol


# load genome
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")


# get 5'-UTR and exons upstream of the mutation of interest 
sequences <- getTxsSeq(ranges_list, hg38)


# create outer primers. the read1, cell barcode and UMI are automatically added upstream of the 3'-end and account for 81bp 
# the outer primer lies 250-350 bp away from the mutation
# it can be that for short genes and mutations which are very upstream there is not enough RNA body to create the outer primers
outer_primers <- TAPseqInput(sequences,
                             product_size_range = c((90+250), (90+350)), 
                             primer_num_return = 5,
                             target_annot = ranges_list)
outer_primers <- designPrimers(outer_primers)


# create intermediate primers. We decided to include and intermediate PCR to increase specificity 
# this is particularly important for lowly expressed genes
# The primers lie between 100 and 200 bp away from the 3'-end of the transcript
message("Designing middle primers")
middle_primers <- TAPseqInput(sequences,
                              product_size_range = c((90+100), (90+200)), 
                              primer_num_return = 5,
                              target_annot = ranges_list)
middle_primers <- designPrimers(middle_primers)


# create inner primers. 
# this align between the mutation and a maximum of bp = read length-10 upstream 
message("Designing inner primers")
inner_primers <- TAPseqInput(sequences,
                             product_size_range = c(90, 90+opt$read_length-15), 
                             primer_num_return = 5,
                             target_annot = ranges_list)
inner_primers <- designPrimers(inner_primers)


# pick puter primers based on the penalty score
best_outer_primers <- pickPrimers(outer_primers, by = "penalty")


# pick middle primers based on the penalty score
best_middle_primers <- pickPrimers(middle_primers, by = "penalty")


# pick inner primers based on penalty score
best_inner_primers <- pickPrimers(inner_primers, by = "penalty")


message("Writing primers to output files")


# final data frame with column stating distance to the mutation
# the fragment length after each PCR is also specified
primers_table <- rbind(primerDataFrame(best_outer_primers) %>% mutate(primer_id = paste0(primer_id, "_outer")),
                       primerDataFrame(best_middle_primers) %>% mutate(primer_id = paste0(primer_id, "_middle")),
                       primerDataFrame(best_inner_primers) %>% mutate(primer_id = paste0(primer_id, "_inner"))) %>%
                        arrange(primer_id) %>% 
                        mutate(distance_mutation = pcr_product_size-81) %>%
                        dplyr::rename(symbol = seq_id) %>%
                        left_join(hits_table %>% dplyr::select(symbol, distance_3_end)) %>%
                        mutate(fragment_length = distance_mutation+distance_3_end+81) %>%
                        dplyr::select(-distance_3_end)


# write final table as csv
write_csv(primers_table, file = paste0(opt$outprefix, "primer_details.csv"))


# export primers as BED ranges
exportPrimerTrack(createPrimerTrack(best_outer_primers, color = "steelblue3"),
                  createPrimerTrack(best_middle_primers, color = "green"),
                  createPrimerTrack(best_inner_primers, color = "red"),
                  con = paste0(opt$outprefix,"primers.bed"))


# add commands to include the adaptor and stagger sequences into the inner primers
staggered_primers <- lapply(1:length(hits_table$symbol), function(i){


  # get sequence of outer primer
  outer_sequence <- primers_table %>% filter(symbol == hits_table$symbol[i]) %>% dplyr::slice(3) %>%
                      pull(sequence)

  # get sequence of middle primer
  middle_sequence <- primers_table %>% filter(symbol == hits_table$symbol[i]) %>% dplyr::slice(2) %>%
                      pull(sequence)


  # get sequence of inner primer
  inner_sequence <- primers_table %>% filter(symbol == hits_table$symbol[i]) %>% dplyr::slice(1) %>%
                      pull(sequence)


  # create staggered primers by adding partial Read2 handle + stagger
  staggered_primers <- c(paste0("CACCCGAGAATTCCA", inner_sequence),
                         paste0("CACCCGAGAATTCCAA", inner_sequence),
                         paste0("CACCCGAGAATTCCATT", inner_sequence),
                         paste0("CACCCGAGAATTCCACAT", inner_sequence))


  # generate a table with the final sequences
  complete_table <- data.frame(primer_name = paste(opt$patient_name, hits_table$symbol[i], 
                                                   c("inner_1", "inner_2", "inner_3", "inner_4", "middle", "outer"),
                                                   sep = "_"),
                               sequence = c(staggered_primers, middle_sequence, outer_sequence))

  complete_table

})


# create a table with primer sequences for all target mutations
complete_primer_table <- do.call("bind_rows", staggered_primers)


# export table as csv
write_csv(complete_primer_table, file = paste0(opt$outprefix, "primers_staggered.csv"))


# write a table with the selected mutations in the primers folder
write_csv(hits_table %>% dplyr::rename(ref = REF, alt = ALT), 
          file = paste0(opt$outprefix, "selected_variants.csv"))


message("Done!")

