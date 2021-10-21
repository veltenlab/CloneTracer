# ExomePipeline snakemake

This folder comprises the snakemake pipeline to execute the Exome pipeline.

## Input

The pipeline requires raw, paired-end fastq input files. Please note that the files have to end with *_R1.fastq.gz* for read 1 and with *_R3.fastq.gz* for read 2. Files can further be split into two different lanes/flow cells (to be specified in the *config.yml* under *patients_novaseq*) and otherwise under patients_nextseq.

## Output

The pipeline will generate a number of output files, such as primer pairs for TAP-seq, as well as further exploratory statistics and plots in the folders *data* and *results*.
