library(vcfR)
library(tidyverse)
library(dplyr)
library(reader)
library(deepSNV)
library(VariantAnnotation)

# Read tsv with AML specific variants and extract cosmic ID

aml_tsv_sits <- read_delim("../../../../external_data/cosmic/AML_mutations.tsv", delim = "\t",
                      col_names = F) %>% pull(X18)

saveRDS(aml_tsv_sits, "aml_cosmic_ids.rds")

# Filter AML variants in VCF file

#############################################################################################################################

# Read complete COSMIC VCF

cosmic_vcf <- readVcf("../../../../external_data/cosmic/CosmicCodingMuts.vcf")

# Function to filter variants present in AML.

filter_aml_mutations <- function(x){
  
  aml_ids <- readRDS("aml_cosmic_ids.rds")
  
  as.vector(info(x)$LEGACY_ID %in% aml_ids)
  
}

# Set filter function for VCF

vcf_filter <- FilterRules(list(aml = filter_aml_mutations))

# Filter VCF file for AML mutations.

tbx <- TabixFile("../../../../external_data/cosmic/CosmicCodingMuts.vcf.gz")

filterVcf(tbx, 
          "hg38",
          "../../../../external_data/cosmic/aml_mutations.vcf",
          filters = vcf_filter,
          verbose = TRUE)

# Change chromosome annotation style 

aml_variants <- readVcf("../../../../external_data/cosmic/aml_mutations.vcf")

writeVcf(aml_variants, "../../../../external_data/cosmic/aml_mutations.vcf")

########################################################################################################################

# Subset AML mutations which are highly expressed 

aml_tsv <- read_delim("../../../../external_data/cosmic/AML.tsv", delim = "\t", col_names = F)

# Extract genes to compute gene expression

aml_genes <- gsub("^(.*?)_.*", "\\1", aml_tsv %>% pull(X1)) %>% unique()

# Select AML genes present in the MUTAseq expression object

aml_genes <- aml_genes[which(aml_genes %in% rownames(GetAssayData(mutaseq_seurat, slot = "counts", assay = "RNA") %>% 
                                                       as.data.frame()))]

# Get genes with high expression in MUTAseq patients

aml_gene_expression <- GetAssayData(mutaseq_seurat, slot = "counts", assay = "RNA") %>% 
                          as.data.frame() %>% t() %>% as.data.frame() %>% 
                          dplyr::select(aml_genes) %>% 
                          rownames_to_column(var = "patient") %>% 
                          separate(patient, into = c("name", "cell_type"), sep = "_", extra = "drop") %>% 
                          mutate(name = if_else(startsWith(name, "plate"), "healthy", name)) %>% 
                          group_by(name) %>% dplyr::select(-cell_type) %>% 
                          summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
                          column_to_rownames(var = "name") %>% t() %>% 
                          as.data.frame() %>% 
                          rownames_to_column(var = "symbol") %>% 
                          mutate_if(is.numeric, floor)

# Compute average expression across 5 samples

aml_gene_expression <- aml_gene_expression %>% 
                        mutate(avg_expression = rowMedians(aml_gene_expression %>% dplyr::select_if(is.numeric) %>% as.matrix()))

high_aml_genes <- aml_gene_expression %>% filter(avg_expression > 30) %>% 
                    pull(symbol)

saveRDS(high_aml_genes, "../../../../external_data/cosmic/high_aml_genes.rds")

###########################################################################################################################

# Subset cosmic AML variants which are highly expressed. 

all_aml_mutations <- readVcf("../../../../external_data/cosmic/aml_mutations.vcf.gz")

# Rename genes in the VCF to eliminate ensembl transcript information. 

info(all_aml_mutations)$GENE <- gsub("^(.*)_.*", "\\1", info(all_aml_mutations)$GENE) 

# Create a function to filter genes which are highly expressed. 

filter_high <- function(x){
  
  aml_genes <- readRDS("../../../../external_data/cosmic/high_aml_genes.rds")
  
  as.vector(info(x)$GENE %in% aml_genes)
  
}

vcf_filter <- FilterRules(list(aml = filter_high))

tbx <- TabixFile("../../../../external_data/cosmic/aml_mutations.vcf.gz")

filterVcf(tbx, 
          "hg38",
          "../../../../external_data/cosmic/aml_high_expressed_mutations.vcf",
          filters = vcf_filter,
          verbose = TRUE)
