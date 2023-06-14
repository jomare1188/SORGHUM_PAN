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

sample_lists = pd.read_csv("samples_Rio.tsv", sep="\t")

genotypes = sample_lists["genotype"]
samples = sample_lists["sample"]

rule all:
	input:
		[f"MyAssembly_{genotypes[i]}/1_raw_reads_in_fastq_format/{samples[i]}_1.fastq" for i in range(len(genotypes))],
		[f"MyAssembly_{genotypes[i]}/1_raw_reads_in_fastq_format/{samples[i]}_2.fastq" for i in range(len(genotypes))],
		[f"MyAssembly_{genotypes[i]}/stranded/{samples[i]}/aux_info/meta_info.json" for i in range(len(genotypes))],
		"stranded_samples_Rio.tsv"
#		[f"samples_{genotypes[i]}x2.csv" for i in range(len(genotypes))]

rule download_fastq:
	output:
		R1 = "MyAssembly_{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq",
		R2 = "MyAssembly_{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq"
	threads: 1
	resources:
		load=50
	params:
		genotype="{genotype}",
		server="figsrv",
		jobname="download"
	conda:
		"conda_envs/download.yml"
	log:
		"MyAssembly_{genotype}/logs/download_fastq/{sample}.log"
	shell:
		"""
		cd MyAssembly_{params.genotype}/1_raw_reads_in_fastq_format && \
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
rule get_stranded:
	input:
                salmon_index = "salmon_index/",
                R1 = "MyAssembly_{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq",
                R2 = "MyAssembly_{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq"
	output:
		stranded = "MyAssembly_{genotype}/stranded/{sample}/aux_info/meta_info.json",
		outdir=directory("MyAssembly_{genotype}/stranded/{sample}")
	conda:
		"conda_envs/salmon.yml"
	params:
                server="figsrv",
		jobname="get_stranded"
	resources:
                load=1
	threads: 20
	log:
                "MyAssembly_{genotype}/logs/stranded/{sample}.log"
	shell:
		"salmon quant -i {input.salmon_index} -p {threads} -l A -1 {input.R1} -2 {input.R2} -o {output.outdir} 2> {log}"
rule filter_stranded:
	input:
		stranded = [f"MyAssembly_{genotypes[i]}/stranded/{samples[i]}/aux_info/meta_info.json"  for i in range(len(genotypes))]
	output:
		"stranded_samples_Rio.tsv"
	threads: 1
	resources:
		load=1
	conda:
		"conda_envs/jq.yml"
	params:
		server="figsrv",
		jobname="filter_stranded"
	shell:
		"""
		jq -r '.library_types[]' {input.stranded} > lib.txt && \
		sed  '1i lib' lib.txt > lib2.txt && \
		paste samples_Rio.tsv lib2.txt > stranded_status_Rio.tsv && \
		awk '$3 ~ /S/' stranded_status_Rio.tsv > stranded_samples_Rio.tsv
		"""

rule salmon_index:
        input:
                transcriptome="rna.fna"
        output:
                salmon_index=directory("salmon_index")
        params:
                jobname="index"
        resources:
                load=1
        conda:
                "conda_envs/salmon.yml"
        threads: 30
        log:
                "salmon_index/log.txt"
        shell:
                "salmon index -t {input.transcriptome} -p {threads} -i {output.salmon_index} 2> {log}"
