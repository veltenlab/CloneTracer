# Rules to analyse exome sequencing from T cells and tumor cells. It outputs a vcf with potential annotated variants

### input, output and shell paths are all relative to the project directory ###

configfile: "config.yml"

rule all:
    input:
        expand("data/{patient}/{population}/duplicates_report.txt",
                patient = config["patient_names"]["mutations"],
                population = config["sample_population"]),
        expand("data/{patient}/{population}/recalibration_files.table",
                patient = config["patient_names"]["mutations"],
                population = config["sample_population"]),
        expand("data/{patient}/{population}/recalibrated_{population}.bam",
                patient = config["patient_names"]["mutations"],
                population = config["sample_population"]),
        expand("results/{patient}/qc_exome/{population}_summary_recalibration.pdf",
                patient = config["patient_names"]["mutations"],
                population = config["sample_population"]),
        expand("data/{patient}/{population}/qc_metrics.txt",
                patient = config["patient_names"]["mutations"],
                population = config["sample_population"]),
        expand("results/{patient}/qc_exome/qc_coverage.png",
                patient = config["patient_names"]["mutations"],
                population = config["sample_population"]),
        expand("data/{patient}/variants/variants.hg38_multianno.vcf",
               patient = config["patient_names"]["mutations"]),
        expand("results/{patient}/annotated_variants.csv",
                patient = config["patient_names"]["mutations"]),
        expand("results/{patient_primers}/primers/primer_details.csv",
               patient_primers = config["patient_names"]["primers"])

import glob

# this is necessary to avoid ambiguous rules
wildcard_constraints:
    patients_nextseq='|'.join(config["patients_nextseq"]),
    patients_novaseq='|'.join(config["patients_novaseq"])

# rule to merge fastqs from NovaSeq sequencer
rule merge_fastqs:
    input:
        fastqs_read1 = lambda wildcards: sorted(glob.glob("raw_data/{patients_novaseq}/{population}/*R1*".format(patients_novaseq = wildcards.patients_novaseq,
                                                                                                               population = wildcards.population))),
        fastqs_read2 = lambda wildcards: sorted(glob.glob("raw_data/{patients_novaseq}/{population}/*R3*".format(patients_novaseq = wildcards.patients_novaseq,
                                                                                                               population = wildcards.population)))
    output:
        merged_read1 = temp("raw_data/{patients_novaseq}/{population}/read1.fastq.paired.fq.gz"),
        merged_read3 = temp("raw_data/{patients_novaseq}/{population}/read3.fastq.paired.fq.gz")
    log:
        "raw_data/{patients_novaseq}/{population}/logs/merge_fastqs.log"
    shell:
        "cat {input.fastqs_read1} > {output.merged_read1} ; "
        "cat {input.fastqs_read2} > {output.merged_read3} "
        "2> {log}"

# sort fastq (for some reason fastqs in P270 and P221 are not read sorted...)
rule sort_fastqs:
    input:
        fastq1 = lambda wildcards: glob.glob('raw_data/{patients_nextseq}/{population}/*R1*'.format(patients_nextseq=wildcards.patients_nextseq,
                                                                                                    population=wildcards.population)),
        fastq2 = lambda wildcards: glob.glob('raw_data/{patients_nextseq}/{population}/*R3*'.format(patients_nextseq=wildcards.patients_nextseq,
                                                                                                    population=wildcards.population))
    output:
        sorted_read1 = temp("raw_data/{patients_nextseq}/{population}/read1.fastq.paired.fq.gz"),
        sorted_read3 = temp("raw_data/{patients_nextseq}/{population}/read3.fastq.paired.fq.gz")
    params:
        gz_read1 = "raw_data/{patients_nextseq}/{population}/read1.fastq",
        gz_read3 = "raw_data/{patients_nextseq}/{population}/read3.fastq",
        uncomp_read1 = "raw_data/{patients_nextseq}/{population}/read1.fastq.paired.fq",
        uncomp_read3 = "raw_data/{patients_nextseq}/{population}/read3.fastq.paired.fq"
    priority: 10
    log:
        "raw_data/{patients_nextseq}/{population}/logs/sort_fasqs.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "gunzip -c {input.fastq1} > {params.gz_read1} ; "
        "gunzip -c {input.fastq2} > {params.gz_read3} ; "
        "fastq_pair -t 500000000 {params.gz_read1} {params.gz_read3} ; "
        "gzip {params.uncomp_read1} ; "
        "gzip {params.uncomp_read3} ; "
        "rm {params.gz_read1} {params.gz_read3} ; "
        "rm raw_data/{wildcards.patients_nextseq}/{wildcards.population}/*single* "
        "2> {log}"


# align reads to reference genome using bwa
rule align_reads:
    input:
        fastq1 = "raw_data/{patient}/{population}/read1.fastq.paired.fq.gz",
        fastq2 = "raw_data/{patient}/{population}/read3.fastq.paired.fq.gz"
    output:
        temp("data/{patient}/{population}/raw_aligned.sam")
    params:
        reference_genome = config["ref_files"]["genome"],
        threads = config["bwa_resources"]["threads"]
    log:
        "data/{patient}/{population}/logs/alignment.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "bwa mem -t {params.threads} {params.reference_genome} {input.fastq1} {input.fastq2} > {output} "
        "2> {log}"


# Convert sam to bam:
rule sam_to_bam:
    input:
        "data/{patient}/{population}/raw_aligned.sam"
    output:
        temp("data/{patient}/{population}/raw_aligned.bam")
    log:
        "data/{patient}/{population}/logs/sam_to_bam.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "samtools view -S -b {input} > {output}"


# sort BAM file
rule sort_bam:
    input:
        "data/{patient}/{population}/raw_aligned.bam"
    output:
        temp("data/{patient}/{population}/aligned_sorted.bam")
    log:
        "data/{patient}/{population}/logs/sort_bam.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "picard SortSam "
        "INPUT={input} "
        "OUTPUT={output} "
        "SORT_ORDER=coordinate "
        "TMP_DIR=./tmp "
        "2> {log}"


# index BAM file
rule index_bam:
    input:
        "data/{patient}/{population}/aligned_sorted.bam"
    output:
        temp("data/{patient}/{population}/aligned_sorted.bai")
    log:
        "data/{patient}/{population}/logs/index_bam.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "samtools index -b "
        "{input} "
        "{output} "
        "TMP_DIR=./tmp "
        "2> {log}"

# Add readgroup tags with the sample population (required for the downstream tools
rule tag_bam:
    input:
        bam = "data/{patient}/{population}/aligned_sorted.bam",
        bai = "data/{patient}/{population}/aligned_sorted.bai"
    output:
        temp("data/{patient}/{population}/aligned_tagged.bam")
    log:
        "data/{patient}/{population}/logs/tag_bam.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "picard AddOrReplaceReadGroups "
        "INPUT={input.bam} "
        "OUTPUT={output} "
        "RGLB={wildcards.population} "
        "RGPL={wildcards.population} "
        "RGPU={wildcards.population} "
        "RGSM={wildcards.population} "
        "TMP_DIR=./tmp "
        "2> {log}"


# Mark duplicated read pairs
rule mark_duplicates:
    input:
        bam = "data/{patient}/{population}/aligned_tagged.bam",
        bai = "data/{patient}/{population}/aligned_sorted.bai"
    output:
        bam = temp("data/{patient}/{population}/marked_duplicates.bam"),
        summary = "data/{patient}/{population}/duplicates_report.txt"
    params:
        gatk_jar = config["jar_files"]["gatk"],
        initial_ram = config["java_specifications"]["initial_ram"],
        max_ram = config["java_specifications"]["max_ram"]
    log:
        "data/{patient}/{population}/logs/tag_duplicates.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "picard MarkDuplicates "
        "I={input.bam} "
        "O={output.bam} "
        "M={output.summary} "
        "TMP_DIR=./tmp "
        "2> {log}"


# sort bam after marking duplicates
rule sort_marked_bam:
    input:
        "data/{patient}/{population}/marked_duplicates.bam"
    output:
        temp("data/{patient}/{population}/marked_duplicates_sorted.bam")
    log:
        "data/{patient}/{population}/logs/sort_marked_bam.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "picard SortSam "
        "I={input} "
        "O={output} "
        "SORT_ORDER=coordinate "
        "TMP_DIR=./tmp "
        "2> {log}"

# Generate base quality score recalibration table
rule recalibration_table:
    input:
        bam = "data/{patient}/{population}/marked_duplicates_sorted.bam"
    output:
        table = "data/{patient}/{population}/recalibration_files.table"
    params:
        reference = config["ref_files"]["genome"],
        vcf = config["ref_files"]["dbSNP"],
        gatk_jar = config["jar_files"]["gatk"],
        initial_ram = config["java_specifications"]["initial_ram"],
        max_ram = config["java_specifications"]["max_ram"]
    log:
        "data/{patient}/{population}/logs/recalibration_table.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "java -Djava.io.tmpdir=./tmp -Xms{params.initial_ram}G -Xmx{params.max_ram}G -jar {params.gatk_jar} BaseRecalibrator "
        "-I {input.bam} "
        "-R {params.reference} "
        "--known-sites {params.vcf} "
        "-O {output.table} "
        "--tmp-dir ./tmp "
        "2> {log}"

# Recalibrate base quality scores
rule recalibrate_scores:
    input:
        bam = "data/{patient}/{population}/marked_duplicates_sorted.bam",
        recalibration_table = "data/{patient}/{population}/recalibration_files.table"
    output:
        bam = "data/{patient}/{population}/recalibrated_{population}.bam"
    params:
        reference = config["ref_files"]["genome"],
        gatk_jar = config["jar_files"]["gatk"],
        initial_ram = config["java_specifications"]["initial_ram"],
        max_ram = config["java_specifications"]["max_ram"]
    log:
        "data/{patient}/{population}/logs/recalibrate_scores.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "java -Djava.io.tmpdir=./tmp -Xms{params.initial_ram}G -Xmx{params.max_ram}G -jar {params.gatk_jar} ApplyBQSR "
        "-I {input.bam} "
        "-R {params.reference} "
        "-bqsr-recal-file {input.recalibration_table} "
        "-O {output.bam} "
        "--tmp-dir ./tmp "
        "2> {log}"

# Summary of the base score recalibration step
rule summary_recalibration:
    input:
        table = "data/{patient}/{population}/recalibration_files.table"
    output:
        plots = "results/{patient}/qc_exome/{population}_summary_recalibration.pdf"
    log:
        "data/{patient}/{population}/logs/recalibration_summary.log"
    params:
        gatk_jar = config["jar_files"]["gatk"],
        initial_ram = config["java_specifications"]["initial_ram"],
        max_ram = config["java_specifications"]["max_ram"]
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "module use /software/as/el7.2/EasyBuild/CRG/modules/all ; "
        "module load R ; "
        "java -Djava.io.tmpdir=./tmp -Xms{params.initial_ram}G -Xmx{params.max_ram}G -jar {params.gatk_jar} AnalyzeCovariates "
        "-bqsr {input.table} "
        "-plots {output.plots} "
        "--tmp-dir ./tmp "
        "2> {log}"

# Get data to make cumulative coverage plot
rule qc_target_coverage:
    input:
        bam = "data/{patient}/{population}/recalibrated_{population}.bam"
    output:
        table = "data/{patient}/{population}/qc_metrics.txt"
    params:
        reference = config["ref_files"]["genome"],
        baits = config["ref_files"]["exome_baits"],
        targets = config["ref_files"]["exome_targets"]
    log:
        "data/{patient}/{population}/logs/qc_target_coverage.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "picard CollectHsMetrics "
        "I={input.bam} "
        "O={output.table} "
        "R={params.reference} "
        "BAIT_INTERVALS={params.baits} "
        "TARGET_INTERVALS={params.targets} "
        "2>{log}"

# Plot QC target coverage
rule plot_target_coverage:
    input:
        table_cd3 = "data/{patient}/cd3_cells/qc_metrics.txt",
        table_tumor = "data/{patient}/tumor_cells/qc_metrics.txt"
    output:
        "results/{patient}/qc_exome/qc_coverage.png"
    log:
        "results/{patient}/logs/qc_coverage_plot.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "Rscript scripts/qc_exome.R "
        "-r {input.table_cd3} "
        "-t {input.table_tumor} "
        "-o {output} "
        "2> {log}"


# Detect somatic variant sites with mutect2 including germline variants
rule variants_mutect2:
    input:
        tcell_bam = "data/{patient}/cd3_cells/recalibrated_cd3_cells.bam",
        tumor_bam = "data/{patient}/tumor_cells/recalibrated_tumor_cells.bam"
    output:
        vcf = temp("data/{patient}/variants/variant_sites_m2.vcf.gz"),
        stats = "data/{patient}/variants/variant_sites_m2.vcf.gz.stats"
    params:
        name_tumor = config["sample_population"][1],
        name_tcells = config["sample_population"][0],
        reference = config["ref_files"]["genome"],
        gatk_jar = config["jar_files"]["gatk"],
        initial_ram = config["java_specifications"]["initial_ram"],
        max_ram = config["java_specifications"]["max_ram"]
    log:
        "data/{patient}/variants/logs/variants_mutect2.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "java -Djava.io.tmpdir=./tmp -Xms{params.initial_ram}G -Xmx{params.max_ram}G -jar {params.gatk_jar} Mutect2 "
        "-R {params.reference} "
        "-I {input.tcell_bam} "
        "-I {input.tumor_bam} "
        "-tumor {params.name_tumor} "
        "-normal {params.name_tcells} "
        "--genotype-germline-sites "
        "-O {output.vcf} "
        "--tmp-dir ./tmp "
        "2> {log}"

# Filter variants from default mutect variant calling
rule filter_variants:
    input:
        vcf = "data/{patient}/variants/variant_sites_m2.vcf.gz"
    output:
        vcf = temp("data/{patient}/variants/filtered_variants.vcf.gz")
    params:
        reference = config["ref_files"]["genome"],
        gatk_jar = config["jar_files"]["gatk"],
        initial_ram = config["java_specifications"]["initial_ram"],
        max_ram = config["java_specifications"]["max_ram"]
    log:
        "data/{patient}/variants/logs/filter_variants.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "java -Djava.io.tmpdir=./tmp -Xms{params.initial_ram}G -Xmx{params.max_ram}G -jar {params.gatk_jar} FilterMutectCalls "
        "-R {params.reference} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--tmp-dir ./tmp "
        "2> {log}"

# Functionally annotate filtered somatic variants using ANNOVAR. This adds gene annotation, mutation consequence and cosmic annotation.
rule annotate_variants:
    input:
        vcf = "data/{patient}/variants/filtered_variants.vcf.gz"
    output:
        table = "data/{patient}/variants/variants.hg38_multianno.vcf"
    params:
        script = config["annovar_files"]["annovar_annotate"],
        ref_genome = config["ref_genome"],
        annovar_database = config["annovar_files"]["annovar_database"]
    log:
        "data/{patient}/variants/logs/annotate_variants.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "perl {params.script} "
        "{input.vcf} "
        "{params.annovar_database} "
        "-out data/{wildcards.patient}/variants/variants "
        "-buildver {params.ref_genome} "
        "-remove "
        "-protocol refGene,cosmic "
        "-operation g,f "
        "-polish "
        "-vcfinput "
        "2>{log}"

# make final variant table and create BED file with coordinates of the genes of interest
rule variant_table:
    input:
        "data/{patient}/variants/variants.hg38_multianno.vcf"
    output:
        csv = temp("results/{patient}/variants_table.csv"),
        bed = "results/{patient}/genomic_files/gene_coordinates.bed"
    params:
        gtf = config["ref_files"]["gtf"]
    log:
        "results/{patient}/logs/final_variant_table.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "Rscript scripts/filter_variants.R "
        "-i {input} "
        "-g {params.gtf} "
        "-o {output.csv} "
        "-b {output.bed} "
        "2> {log}"


#rule subset gene expression BAM file to compute counts/cell for each mutated gene
rule subset_bams:
    input:
        bed = "results/{patient}/genomic_files/gene_coordinates.bed"
    output:
        bam = "results/{patient}/genomic_files/subset.bam"
    params:
        patients_bams = config["ref_files"]["bam_polyA"],
        directory = "results/{patient}/genomic_files/"
    log:
        "results/{patient}/logs/subset_merge_bams.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "Rscript scripts/subset_merge_bams.R "
        "-b {input.bed} "
        "-d {params.directory} "
        "-p {wildcards.patient} "
        "-o {output.bam} "
        "2> {log}"

# rule to compute distance from the filtered mutations to the polyA tail
rule compute_distance:
    input:
        csv = "results/{patient}/variants_table.csv",
        bam = "results/{patient}/genomic_files/subset.bam"
    output:
        "results/{patient}/annotated_variants.csv"
    params:
        gtf = config["ref_files"]["gtf"],
        directory = "results/{patient}/genomic_files/"
    log:
        "results/{patient}/logs/distance_polyA.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "Rscript scripts/distance_polyA.R "
        "-i {input.csv} "
        "-b {input.bam} "
        "-g {params.gtf} "
        "-p {wildcards.patient} "
        "-d {params.directory} "
        "-o {output} "
        "2> {log}"


# rule to generate outer, middle and inner (staggered) primers
# it needs selected_variants.txt file with the selected mutations one per line.
# it has to follow the format genesymbol_position \n.
rule design_primers:
    input:
        mutations = "results/{patient_primers}/annotated_variants.csv",
        selected_variants = "results/{patient_primers}/primers/selected_variants.txt"
    output:
        primers = "results/{patient_primers}/primers/primer_details.csv"
    params:
        gtf = config["ref_files"]["gtf"],
        bed = "results/{patient_primers}/genomic_files/polyA_final_selection.bed",
        read_length = config["read_length"]["default"]
    log:
        "data/{patient_primers}/logs/design_primers.log"
    conda:
        "../snakemake/envs/exome.yml"
    shell:
        "Rscript scripts/design_mutation_primers.R "
        "-i {input.mutations} "
        "-l {input.selected_variants} "
        "-o results/{wildcards.patient_primers}/primers/ "
        "-r {params.read_length} "
        "-g {params.gtf} "
        "-b {params.bed} "
        "-p {wildcards.patient_primers} "
        "2> {log}"
