library(vcfR)
library(tidyverse)
library(dplyr)
library(reader)
library(deepSNV)
library(GenomicRanges)
library(BiocParallel)
library(TAPseq)

# 
# gene annotations for chromosome 11 genomic region
data("chr11_genes")

# convert to GRangesList containing annotated exons per gene. for the sake of time we only use
# a small subset of the chr11 genes. use the full set for a more realistic example
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)
target_genes <- target_genes[18:27]

# K562 Drop-seq read data (this is just a small example file within the R package)
dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TAPseq")

# register backend for parallelization
register(MulticoreParam(workers = 5))

# infer polyA sites from Drop-seq data
polyA_sites <- inferPolyASites(target_genes, bam = dropseq_bam, polyA_downstream = 50,
                               wdsize = 100, min_cvrg = 1, parallel = TRUE)


# truncate transcripts at inferred polyA sites
truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = polyA_sites, parallel = TRUE)

library(BSgenome)

# human genome BSgenome object (needs to be istalled from Bioconductor)
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# get sequence for all truncated transcripts
txs_seqs <- getTxsSeq(truncated_txs, genome = hg38)


outer_primers <- TAPseqInput(txs_seqs, target_annot = truncated_txs,
                             product_size_range = c(350, 500), primer_num_return = 5)

# design 5 outer primers for each target gene
outer_primers <- designPrimers(outer_primers)
