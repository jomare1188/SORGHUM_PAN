#!/usr/bin/env python

# An Automated Pipeline to Transcriptome Assembly and Quality Assessment
#
# Author: Felipe Vaz Peres
#
# Preparation:
# 1) Setup 'config.yaml' file  with softwares path
# 2) Create 'samples_{genotype}.csv' file and add SRA identifiers (e.g SRR5258954,SRR5258955,SRR5258994,SRR5258995)
# 3) Create 'parts.csv' file and add the value of parts you want the kraken file to be split (e.g 00,01,02 for 3 equal parts)
# 4) Setup 'GENOTYPE' variable with the genotype name (e.g GENOTYPE=QN05-1509)
#
# Usage (remove the -n to dont show dry-run and start the jobs):
# 1) Load modules: BUSCO/3.0; transrate/1.0.3
# 2) Run the following command:
# snakemake -np -s Snakefile \
# --cluster "qsub -q all.q -V -cwd -l h={params.server} -pe smp {threads} -l mem_free={resources.mem_free}G" \
# --jobs 10
#
# Build DAG:
#
# snakemake -s Snakefile --dag | dot -Tsvg > dag.svg
import pandas as pd
import yaml

sample_lists = pd.read_csv("samples_RTx430.txt", sep="\t")
genotypes = sample_lists["genotype"]
samples = sample_lists["sample"]

rule all:
	input:
		[f"../results/{genotypes[i]}/1_raw_reads_in_fastq_format/{samples[i]}_1.fastq" for i in range(len(genotypes))],
		[f"../results/{genotypes[i]}/1_raw_reads_in_fastq_format/{samples[i]}_2.fastq" for i in range(len(genotypes))],
#		"../results/RTx430/salmon_index",
		[f"../results/{genotypes[i]}/2_trimmed/paired_{samples[i]}_1.fastq" for i in range(len(genotypes))]

rule download_fastq:
	output:
		R1 = "../results/{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq",
		R2 = "../results/{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq"
	threads: 1
	resources:
		load=100
	params:
		genotype="{genotype}",
		server="neotera",
		jobname="download"
	conda:
		"conda_envs/download.yml"
	log:
		"../results/{genotype}/1_raw_reads_in_fastq_format/{sample}.log"
	shell:
		"""
		cd ../results/{wildcards.genotype}/1_raw_reads_in_fastq_format && \
		ffq --ftp {wildcards.sample} | grep -Eo '\"url\": \"[^\"]*\"' | grep -o '\"[^\"]*\"$' | xargs wget && \
		gzip -dc < {wildcards.sample}_1.fastq.gz > {wildcards.sample}_1.fastq && \
		gzip -dc < {wildcards.sample}_2.fastq.gz > {wildcards.sample}_2.fastq && \
		python2.7 ./../../good_headers.py -1 {wildcards.sample}_2.fastq -2 {wildcards.sample}_1.fastq && \
		rm {wildcards.sample}_1.fastq.gz && \
		rm {wildcards.sample}_2.fastq.gz && \
		mv sraheaderfixed_{wildcards.sample}_2.fastq {wildcards.sample}_2.fastq && \
		mv sraheaderfixed_{wildcards.sample}_1.fastq {wildcards.sample}_1.fastq && \
		cd -
		"""

rule trimmomatic:
	input:
		R1 = "../results/{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq",
		R2 = "../results/{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq"
	output:
		P1 = "../results/{genotype}/2_trimmed/paired_{sample}_1.fastq",
		P2 = "../results/{genotype}/2_trimmed/paired_{sample}_2.fastq",
		U1 = "../results/{genotype}/2_trimmed/unpaired_{sample}_1.fastq",
		U2 = "../results/{genotype}/2_trimmed/unpaired_{sample}_2.fastq"
	params:
		jobname="trimmomatic",
		server="neotera"
	resources:
		load=1
	conda:
		"conda_envs/TRIMMOMATIC.yml"
	threads: 20
	log:
		"../results/{genotype}/trimmomatic/{sample}_log.txt"
	shell:
		"trimmomatic PE -trimlog trimmomatic.log -threads {threads} {input.R1} {input.R2} {output.P1} {output.U1} {output.P2} {output.U2} MINLEN:80"

rule salmon_index:
        input:
                transcriptome="../data/GCA_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
        output:
                salmon_index=directory("../results/BTx623/salmon_index")
        params:
                jobname="index",
		server="neotera"
        resources:
                load=1
        conda:
                "conda_envs/salmon.yml"
        threads: 20
        log:
                "../results/BTx623/salmon_index/log.txt"
        shell:
                "salmon index -t {input.transcriptome} -p {threads} -i {output.salmon_index} 2> {log}"
