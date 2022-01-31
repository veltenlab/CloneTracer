import numpy
import HTSeq
import itertools
import collections
import argparse
import csv
import os
import pandas as pd


# function to create and empty GenomicRanges object with exons mitochondria and the genes tagged
def tag_genome(gtf_whole_genome, gene_hits):

    # Create and empty GenomicInterval object which spans the whole genome and mtDNA
    whole_genome = HTSeq.GenomicArrayOfSets("auto", stranded=True)


    # load gene annotation from gtf file
    genome_annotation = HTSeq.GFF_Reader(gtf_whole_genome)


    #load csv file with targeted genes and position of the mutations of interest and pull the column of gene names as set
    mutation_genes = set((pd.read_csv(gene_hits)['symbol']))


    # annotate genes and exons and genes of interest the GenomicInterval object
    for feature in genome_annotation:
        if feature.type == "exon":
            whole_genome[feature.iv] += "exon"
        elif feature.type == "gene":
            whole_genome[feature.iv] += "gene"
        elif feature.attr["gene_name"] in mutation_genes:
            whole_genome[feature.iv] += feature.attr["gene_name"]


    # annotate mitochondrial genome in the GenomicInterval object
    whole_genome[HTSeq.GenomicInterval( "MT", 0, 16570, "+")] += "mito"
    whole_genome[HTSeq.GenomicInterval( "MT", 0,  16570, "-")] += "mito"


    return({'genome': whole_genome, 'genes': mutation_genes})


# count reads which fall in the predefined features
def count_reads(bam_file, annotated_genome, gene_names):

    # read BAM file
    bam = HTSeq.BAM_Reader(bam_file)


    # create a list to store counts. The difference with a normal list is that one can use '+=1'
    # even if the element is not present in the list
    counts = collections.Counter()

    # count reads
    for read in bam:

        # count unmapped reads
        if not read.aligned:

            counts["unmapped"] += 1
            continue

        # a set is an unordered list which has unique elements
        tags = set()

        for iv, val in annotated_genome[read.iv].steps():

            # unite all tags of the genomic interval covered by the read
            tags |= val

        # check if the read aligns to any of the amplified genes
        if bool(tags & gene_names):
            counts[list(tags & gene_names)[0]] += 1

        # count reads aligning to the mitochondrial genome
        elif "mito" in tags:
            counts["mito"] += 1


        elif "exon" in tags:
            counts["exon"] += 1


        elif "gene" in tags:
            counts["gene"] += 1


        elif len(tags) == 0:
            counts["other"] += 1

    return(counts)

def global_function(gtf_whole_genome, bam_file, gene_hits, output):

    annotated_array = tag_genome(gtf_whole_genome, gene_hits)

    counts_table = count_reads(bam_file, annotated_array['genome'], annotated_array['genes'])

    with open(output,'w') as csvfile:
        writer=csv.writer(csvfile)
        writer.writerows(counts_table.items())

if __name__ == "__main__":

      parser = argparse.ArgumentParser()
      parser.add_argument("--bam", help="Input BAM file", type=str)
      parser.add_argument("--output", help="Output csv file", type=str)
      parser.add_argument("--genes", help="csv file with the targeted genes. It should contain 1 column named 'symbol' with the gene names", type=str)
      parser.add_argument("--gtf_genome", help="GTF file covering the entire genome", type=str)
      args = parser.parse_args()
      global_function(args.gtf_genome, args.bam, args.genes, args.output)
