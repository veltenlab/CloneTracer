# GeneExpression snakemake

This folder comprises the snakemake pipeline to execute the GeneExpression (CITE-seq or WTA) pipeline.

## Input

The pipeline requires raw, paired-end fastq input files. Please note that the files have to end with *_read1.fastq.gz* for read 1 and with *_read2.fastq.gz* for read 2. The pipeline accepts either samples from a WTA (in the *config.yml* file under *patients_wta*) experiment or CITE-seq data (under *patients_cite*). Please note that CITE-seq data requires both the gene expression and protein data.

## Output

The pipeline will generate a number of output files, including cellranger outputs, Seurat objects, projections to the reference map of [Triana et al.](https://www.biorxiv.org/content/10.1101/2021.03.18.435922v2), and diagnostic plots.
