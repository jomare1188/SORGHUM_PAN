import pandas as pd

GENOTYPE=["BAZ9504", "BTx623", "cv.BRS330", "Della", "DKS-3707", "DKS-4420", "keller", "M35-1", "Mota", "R9188", "Rioref", "RTx430", "S085", "SM100", "TAM428", "SC187", "Tx2737", "Tx3362", "Tx378", "Tx7000", "Rio", "NSL365694", "PI651187"]

rule all:
        input:
                expand("../results/{genotype}_busco", genotype=GENOTYPE)
rule busco:
        input:
                transcriptome="../../assemblies/contig_300/{genotype}.fasta"
        output:
                busco=directory("../results/{genotype}_busco")
        conda:
                "/Storage/data1/jorge.munoz/SORGHUM_PAN/snakemake/conda_envs/busco.yml"
        threads: 5
	params:
                jobname="{genotype}_busco"
        shell:
            "busco -i {input.transcriptome} -o {output.busco} -c {threads} -m transcriptome -l ../../../busco_db/embryophyta_odb10"

