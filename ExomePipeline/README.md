# ExomePipeline

## A collection of scripts for extracting mutations from exome data, together with the construction of amplicons (primers) for analysis in TAP-seq.

## Dependencies:
* ANNOVAR: https://annovar.openbioinformatics.org
* samtools: https://sourceforge.net/projects/samtools/files/samtools/. Samtools also has to be added to your PATH variable and thus executable from within R.
* gatk Version-4.1.9.0: https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip

## Required Resources:
* (UCSC) Human reference genome: hg38.fa [https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
* dbSNP: dbSNP\_exome.vcf.gz [https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz)
* exomeBaits for PicardTools
* targets.interval_list for PicardTools
* gene annotation: Homo_sapiens.GRCh38.100.chr.gtf [ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.chr.gtf.gz](ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.chr.gtf.gz)
