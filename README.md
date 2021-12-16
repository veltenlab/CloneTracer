# MutaSeq v2

MutaSeq v2 is a methdology to amplify nuclear SNPs and the mitochondrial genome from 3' 10x Genomics cDNA libraries. This GitHub repository comprises all necessary scripts to process data generated using this method. Processing of raw files from nuclear SNPs and mitochondrial libraries are written as Snakemake pipelines separately and can be found in:

* Add mito pipeline as separate folder
* Add nuclear SNP pipeline as separate folder

Furthermore we have developed a Bayesian model to infer clonal hierarchies from scRNAseq data. All required scripts and detailed explanation on how to implement the model can be found in:

* [ClonalInference](Clonal inference)

