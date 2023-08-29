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

#rule transrate:
#        input:
#                transcriptome="MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25_and_k31.Trinity.merged.final.fasta",
#                ref="rna.fna"
#        output:
#                transrate=directory("MyAssembly_{genotype}/9_transrate/assemblies.csv")
#        resources:
#                load=1
#        threads: 40
#        log:
#                "MyAssembly_{genotype}/logs/transrate/{genotype}.log"
#        conda:
#                "/data/j/TRANSCRIPTOMATOR/transrate.yml"
#        shell:
#                "transrate --assembly {input.transcriptome} --reference {input.ref} --threads {threads} --output {output.transrate} 2> {log}"


