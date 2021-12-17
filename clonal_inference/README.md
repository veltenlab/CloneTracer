# Inference of clonal hierarchies


## Description

Here we developed and bayesian model to infer clonal hierarchies from scRNAseq data using nuclear and mitochondrial SNVs as well as CNVs. 
The model uses stochastic variational inference to select the mutation tree with the highest evidence (lowest Evidence Lower Bound, ELBO).
Then it computes the posterior probability for each cell to belong to each clone in the inferred tree.


The model is implemented in [pyro](http://pyro.ai/) a probabilistic programming language written in Python which uses PyTorch as backend.  


## Installation and dependencies


The following python libraries are required to run the model: pyro, pytorch, pickle, pandas, numpy and seaborn. 

We have included a yml which contains all required packages to run the model. A conda environment can be created from it to run the model:

```
git clone https://github.com/veltenlab/MutaSeq-v2
cd clonal_inference
conda env create -f envs/clonal_inference.yml
```
## Input file

The model requires a config file in JSON format (see below how to create a JSON file in R or python, add links to headers). Input files for all patients 
in the manuscript are present in [input data](data).

### Required entries

* **M**: matrix of number of cells x number of mutations (rowsxcolumns) with mutant read counts for nuclear and mitochondrial SNVs. For CNVs it should comprise the total counts in the region of interest (e.g. for trisomy in chr8 it will contain all reads from genes in chr8).
* **N**: matrix number of cells x number of mutations (rowsxcolumns) with reference counts for nuclear and mitochondrial SNVs. For CNVs it should contain the total read counts. 
*  **mut_type**: vector of integers indicating the type of mutations (0: CNV, 1: nuclear SNV or 2: mitochondrial SNV)
*  **mut_names**: vector of strings indicating the mutation names
*  **cell_barcode**: vector with cell barcodes. It should match the order of the M and N matrices.


### Optional entries

* **r_cnv**: required when CNVs are present, vector of indicating whether the CNV is a gain or a loss (1.5 for gain, 0.5 for loss, 0 for SNVs).
* **bulk_M**: vector of counts supporting the mutant allele or CNV from bulk data (e.g. exome counts on mutant allele, number of cells with CNV from clinical karyotyping or FISH). In case of more than one sample 
* **bulk_N**: vector counts supporting the reference allele or diploid data from bulk data (e.g. exome counts on reference allele, number of cells diploid cells from clinical karyotyping or FISH)

