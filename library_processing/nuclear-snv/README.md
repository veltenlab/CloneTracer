# Processing of nuclear SNV libraries from CloneTracer

## Description

The workflow is written as a Snakemake pipeline can be divided into three main parts:

* Reads are tagged and aligned following the standard [DROPseq](https://mccarrolllab.org/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf) pipeline.
* Consensus reads are obtained by by UMI collapsing using fgbio tools [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html) and [CallMolecularConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html). 
* Reference and mutant reads are counted using [pysam](https://pysam.readthedocs.io/en/latest/api.html).

## Dependencies 

We provide a conda environment with all the required package dependencies. It can be installed as follows:

```
git clone https://github.com/veltenlab/CloneTracer
cd library_processing/nuclear-snv
conda env create -f envs/CloneTracer.yml
conda activate clonetracer
```

## Input files

For each library a directory should be created inside the "raw_data" folder with the sample name (e.g. raw_data/P1/ for sample P1). The directory must contain the following files:

* Fastq files for read1 and read2. They must contain "\*R1*" and "\*R2*" in the file name, respectively (e.g. raw_data/P1/P1_R1.fastq.gz for read1 fastq).
* `barcodes.tsv`: list of barcodes present in gene expression library (one barcode per line).
* `selected_variants.csv`: .csv file with gene name, chromosome and position of the amplified mutations as follows:

| symbol      | CHROM  | POS
| ----------- | ------ |-----------
| KRAS        | chr12  | 25245350
| NRAS        | chr1   | 114713909
| IDH2        | chr15  | 90088702

## Fill config.yml

In order to execute the Snakemake pipeline the following entries should be filled in the [config.yml](config.yml) file:

* **patient_id**: sample names. It must be equal to the sample directory created in the "raw_data" folder
* **
[gtf file from ensembl](http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.chr.gtf.gz)
