# Processing of mitochondrial libraries from CloneTracer

## Description

The workflow is written as a Snakemake pipeline and can be divided into three main parts:

* Reads are tagged and aligned to the mitochondrial genome following the standard [DROPseq](https://mccarrolllab.org/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf) pipeline.
* Consensus reads are obtained by UMI collapsing using fgbio tools [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html) and [CallMolecularConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html). 
* Single-cell nucleotide counts are generated using [pysam](https://pysam.readthedocs.io/en/latest/api.html).

## Dependencies 

We provide a conda environment with all the required package dependencies. It can be installed as follows:

```
git clone https://github.com/veltenlab/CloneTracer
cd CloneTracer/library_processing/mitochondria
conda env create -f envs/CloneTracer.yml
conda activate clonetracer
```
This conda environment also contains all packages to process nuclear libraries as described in [processing nuclear libraries](../nuclear-snv) section.

## Input files

For each library a directory should be created inside the "raw_data" folder with the sample name (e.g. `raw_data/P1` for sample P1). The directory must contain the following files:

* Fastq files for read1 and read2. They must contain "\*R1*" and "\*R2*" in the file name, respectively (e.g. `raw_data/P1/P1_R1.fastq.gz` for read1 fastq).
* `barcodes.tsv`: list of barcodes present in gene expression library (one barcode per line). It must follow the cellranger format of barcode-1.

## Fill config.yml

In order to execute the Snakemake pipeline the following entries should be filled in the [config.yml](config.yml) file:

* **patient_ids**: sample names. It must be equal to the sample directory created in the "raw_data" folder.
* **genomic_references**:
  - **fasta**: path to fasta file of the reference genome. Mitochondrial human reference genome can be [downloaded here](http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz).
  - **genome_dir**: directory where STAR index files for the reference genome are stored. See this [tutorial](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html) on how to index a reference genome using STAR aligner. 
  - **gtf**: [gtf file from ensembl](http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.chr.gtf.gz).
  - **dict**: sequence dictionary for the reference genome generated using [CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-) tool from picard.
  - **refFlat**: refFlat file from gtf. It can be generated using [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html) from UCSC tools. 

## Execute pipeline

### Command-line execution

We strongly recommend to run the workflow in a cluster environment with at least 8 available cores. To execute the snakemake pipeline in a SGE/UGE cluster use the following command:

```
snakemake --cluster "qsub -l h_rt=10:00:00 -l virtual_free=60G -pe 8 -e log_files/stderr.txt -o log_files/stdout.txt" --use-conda --keep-going -j 1000 --latency-wait 120
```

For other cluster environments replace the commands in brackets for the corresponding cluster submission command and resource specification flags.

### Snakemake profiles

Snakemake [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html) enable the rule-specific allocation of resources. We provide a profile for the execution of the pipeline in a SGE/UGE cluster:

```
snakemake --profile profiles/mitochondria/
```

In this [github repository](https://github.com/Snakemake-Profiles) you can find snakemake profiles for other cluster environments. 

## Output

The pipeline produces 2 main output files:

```
data/<sample_name>/mito_counts/count_table/<sample_name>_count_table.rds
data/<sample_name>/mito_counts/count_table/<sample_name>_sum_experiment.rds

```

The `count_table.rds` object contains a list with one entry per cell filled with matrix of nuceleotide counts over the mitochondrial genome. This format is compatible with functions `mutationCallsFromBlacklist` and `mutationCallsFromCohort` from the [mitoClone](https://github.com/veltenlab/mitoClone) package for mitochondrial variant discovery.

The `sum_experiment.rds` file contains a SummarizedExperiment object which contains nucleotide counts for each single cell and position as well as additional information such as average base quality and average number of reads supporting each UMI. 
