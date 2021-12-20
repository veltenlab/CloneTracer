# MutaSeq v2

MutaSeq v2 is a methdology to amplify nuclear SNVs and the mitochondrial genome from 3' 10x Genomics cDNA libraries. This GitHub repository comprises all necessary scripts to process data generated using this method. Processing of raw files from nuclear SNVs and mitochondrial libraries are written as Snakemake pipelines separately and can be found in:

* [Processing mitochondrial library](library_processing/mitochondrial_library)
* [Processing nuclear SNVs library](library_processing/nuclearSNVs_library)

Furthermore we have developed a Bayesian model to infer clonal hierarchies from scRNAseq data using nuclear and mitochondrial SNVs as well as CNVs. All required scripts and detailed explanation on how to implement the model can be found in:

* [Clonal inference](clonal_inference)

Finally, scripts to design primers targeting nuclear SNVs of interest can be found in:

* [Primer design for nuclar SNVs library](primer_design)

