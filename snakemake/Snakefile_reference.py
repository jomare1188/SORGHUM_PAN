import pandas as pd

samples = pd.read_csv("samples_Rio.csv")
parts = pd.read_csv("parts.csv")
GENOTYPE='Rio'

rule all:
        input:
                #expand("MyAssembly_{genotype}/6_aligned_reads_bam/{sample}.bam", genotype=GENOTYPE, sample=samples),
                #expand("MyAssembly_{genotype}/7_aligned_reads_bam/{genotype}_merged_sorted.bam", genotype=GENOTYPE, sample=samples),
                #expand("MyAssembly_{genotype}/7_transcriptome/trinity_{genotype}_k25.fasta", genotype=GENOTYPE)
                #expand("MyAssembly_{genotype}/5_trimmed_reads_kraken_reports/{sample}.trimmed.kraken", genotype=GENOTYPE, sample=samples),
                #expand("MyAssembly_{genotype}/6_contamination_removal/{sample}.trimmed.filtered.total.R1.fastq", genotype=GENOTYPE, sample=samples),
                #expand("MyAssembly_{genotype}/6_contamination_removal/{sample}.trimmed.filtered.total.R2.fastq", genotype=GENOTYPE, sample=samples),
                #expand("MyAssembly_{genotype}/6_contamination_removal/{sample}.trimmed.unclassified.total.R1.fastq", genotype=GENOTYPE, sample=samples),
                #expand("MyAssembly_{genotype}/6_contamination_removal/{sample}.trimmed.unclassified.total.R2.fastq", genotype=GENOTYPE, sample=samples),
                 expand("MyAssembly_{genotype}/8_busco/My_Assembly_{genotype}_k25_and_k31_busco", genotype=GENOTYPE),
                 expand("MyAssembly_{genotype}/9_transrate/assemblies.csv", genotype=GENOTYPE)
                #expand("MyAssembly_{genotype}/10_salmon/quant/quant.sf", genotype=GENOTYPE)
rule fastqc:
        output:
                html_1= "MyAssembly_{genotype}/2_raw_reads_fastqc_reports/{sample}_1_fastqc.html",
                zip_1 = "MyAssembly_{genotype}/2_raw_reads_fastqc_reports/{sample}_1_fastqc.zip",
                html_2= "MyAssembly_{genotype}/2_raw_reads_fastqc_reports/{sample}_2_fastqc.html",
                zip_2 = "MyAssembly_{genotype}/2_raw_reads_fastqc_reports/{sample}_2_fastqc.zip"
        threads: 1
        resources:
                load=1
        params:
                genotype="{genotype}"
	conda:
		"/data/j/TRANSCRIPTOMATOR/fastqc.yml"
	log:
		"MyAssembly_{genotype}/logs/fastqc/{sample}.log"
	shell:
		"fastqc -f fastq MyAssembly_{params.genotype}/1_raw_reads_in_fastq_format/{wildcards.sample}_1.fastq -t {threads} -o MyAssembly_{params.genotype}/2_raw_reads_fastqc_reports 2> {log};"
		"fastqc -f fastq MyAssembly_{params.genotype}/1_raw_reads_in_fastq_format/{wildcards.sample}_2.fastq -t {threads} -o MyAssembly_{params.genotype}/2_raw_reads_fastqc_reports 2> {log}"

rule bbduk:
        input:
                "MyAssembly_{genotype}/2_raw_reads_fastqc_reports/{sample}_1_fastqc.html",
        output:
                R1 = "MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.R1.fastq",
                R2 = "MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.R2.fastq",
                refstats = "MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.refstats",
                stats = "MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.stats"
        log:
                "MyAssembly_{genotype}/logs/bbduk/{sample}.log"
        threads: 20
        resources:
                load=20
        params:
                R1 = "MyAssembly_{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq",
                R2 = "MyAssembly_{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq",
	conda:
		"/data/j/TRANSCRIPTOMATOR/bbduk.yml"
	shell:
                "bbduk.sh -Xmx40g threads={threads} in1={params.R1} in2={params.R2} "
		"refstats={output.refstats} stats={output.stats} "
                "out1={output.R1} out2={output.R2} "
                "ref=/data/j/TRANSCRIPTOMATOR/dbs/adapters.fa,"
                "/data/j/TRANSCRIPTOMATOR/dbs/rfam-5.8s-database-id98.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/silva-bac-16s-id90.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/rfam-5s-database-id98.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/silva-bac-23s-id98.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/silva-arc-16s-id95.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/silva-euk-18s-id95.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/silva-arc-23s-id98.fasta,"
                "/data/j/TRANSCRIPTOMATOR/dbs/silva-euk-28s-id98.fasta "
                "minlength=75 qtrim=w trimq=20 tpe tbo 2> {log}"

rule fastqc_after_bbduk:
        input:
                R1 = "MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.R1.fastq",
                R2 = "MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.R2.fastq"
        output:
                html_1= "MyAssembly_{genotype}/4_trimmed_reads_fastqc_reports/{sample}.trimmed.R1_fastqc.html",
                zip_1 = "MyAssembly_{genotype}/4_trimmed_reads_fastqc_reports/{sample}.trimmed.R1_fastqc.zip",
                html_2= "MyAssembly_{genotype}/4_trimmed_reads_fastqc_reports/{sample}.trimmed.R2_fastqc.html",
                zip_2 = "MyAssembly_{genotype}/4_trimmed_reads_fastqc_reports/{sample}.trimmed.R2_fastqc.zip"
        threads: 1
        resources:
                load=1
        params:
                genotype="{genotype}",
	conda:
		"/data/j/TRANSCRIPTOMATOR/fastqc.yml"
	log:
                "MyAssembly_{genotype}/logs/fastqc_after_bbduk/{sample}.log"
	shell:
		"fastqc -t {threads} -f fastq {input.R1} -o MyAssembly_{params.genotype}/4_trimmed_reads_fastqc_reports 2> {log};"
		"fastqc -t {threads} -f fastq {input.R2} -o MyAssembly_{params.genotype}/4_trimmed_reads_fastqc_reports 2> {log}"

rule hisat2_index:
	input:
		reference_genome="/data/j/TRANSCRIPTOMATOR/REFERENCES/{genotype}/GCA_015952705.1_SbicolorRio_v2_genomic.fna"
	output:
		"/data/j/TRANSCRIPTOMATOR/REFERENCES/{genotype}/index.1.ht2"
	params:
		out="/data/j/TRANSCRIPTOMATOR/REFERENCES/{genotype}/index"
	threads: 40
	conda:
		"/data/j/TRANSCRIPTOMATOR/hisat2.yml"
	resources:
		load=1
	log:
		"MyAssembly_{genotype}/logs/hisat_index.log"
	shell:
		"hisat2-build -p {threads} {input.reference_genome} {params.out}"

rule align_reads_hisat2:
    input:
        forward_reads="MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.R1.fastq",
        reverse_reads="MyAssembly_{genotype}/3_trimmed_reads/{sample}.trimmed.R2.fastq",
        reference_genome_index="/data/j/TRANSCRIPTOMATOR/REFERENCES/{genotype}/index.1.ht2"
    output:
        "MyAssembly_{genotype}/5_aligned_reads/{sample}.sam"
    threads: 3
    resources:
        load=1
    params:
        genotype="{genotype}",
	index="/data/j/TRANSCRIPTOMATOR/REFERENCES/{genotype}/index"
    conda:
        "/data/j/TRANSCRIPTOMATOR/hisat2.yml"
    log:
        "MyAssembly_{genotype}/logs/hisat2_align/{sample}.log"
    shell:
        "mkdir -p MyAssembly_{wildcards.genotype}/5_aligned_reads && hisat2 -x {params.index} -1 {input.forward_reads} -2 {input.reverse_reads} -S {output} --threads {threads} 2> {log};"

rule sam_to_bam:
    input:
        sam_file="MyAssembly_{genotype}/5_aligned_reads/{sample}.sam"
    output:
        bam_files="MyAssembly_{genotype}/6_aligned_reads_bam/{sample}.bam"
    params:
        genotype="{genotype}"
    threads: 3
    conda:
        "/data/j/TRANSCRIPTOMATOR/hisat2.yml"
    log:
        "MyAssembly_{genotype}/logs/sam_to_bam/{sample}.log"
    shell:
        "samtools view -@ {threads} -bS {input.sam_file} > {output.bam_files} 2> {log}"

rule merge_and_sort_bam:
    input:
        bam_files=expand("MyAssembly_{genotype}/6_aligned_reads_bam/{sample}.bam", genotype=GENOTYPE, sample=samples)
    output:
        merged_bam="MyAssembly_{genotype}/7_trinity_assembly/{genotype}_merged_sorted.bam"
    params:
        genotype="{genotype}"
    threads: 48
    conda:
        "/data/j/TRANSCRIPTOMATOR/hisat2.yml"
    log:
        "MyAssembly_{genotype}/logs/merge_and_sort_bam.log"
    shell:
        "samtools merge -@ {threads} -f merged.bam {input.bam_files} && samtools sort -@ {threads} -o {output.merged_bam} merged.bam && rm merged.bam 2> {log}"

rule trinity_k25:
    input:
        bam_file="MyAssembly_{genotype}/7_trinity_assembly/{genotype}_merged_sorted.bam"
    output:
        "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25.Trinity-GG.fasta"
    threads: 30
    resources:
        load=1
    params:
        genotype="{genotype}"
    log:
        "MyAssembly_{genotype}/logs/trinity_assembly_k25.log"
    conda:
        "/data/j/TRANSCRIPTOMATOR/singularity.yml"
    shell:
        """
        singularity exec -B /data/j/TRANSCRIPTOMATOR/ -e /data/j/TRANSCRIPTOMATOR/trinityrnaseq.v2.15.0.simg Trinity --genome_guided_bam {input.bam_file} --genome_guided_max_intron 10000 --CPU {threads} --output /data/j/TRANSCRIPTOMATOR/MyAssembly_{params.genotype}/7_trinity_assembly/MyAssembly_{params.genotype}_trinity_k25  --max_memory 150G --full_cleanup --min_contig_length 200 --SS_lib_type RF --__KMER_SIZE 25 2> {log}
        """

rule trinity_k31:
    input:
        bam_file="MyAssembly_{genotype}/7_trinity_assembly/{genotype}_merged_sorted.bam"
    output:
        "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k31.Trinity-GG.fasta"
    threads: 30
    resources:
        load=1
    params:
        genotype="{genotype}"
    log:
        "MyAssembly_{genotype}/logs/trinity_assembly_k31.log"
    conda:
        "/data/j/TRANSCRIPTOMATOR/singularity.yml"
    shell:
        """
        singularity exec -B /data/j/TRANSCRIPTOMATOR/ -e /data/j/TRANSCRIPTOMATOR/trinityrnaseq.v2.15.0.simg Trinity --genome_guided_bam {input.bam_file} --genome_guided_max_intron 10000 --CPU {threads} --output /data/j/TRANSCRIPTOMATOR/MyAssembly_{params.genotype}/7_trinity_assembly/MyAssembly_{params.genotype}_trinity_k31 --max_memory 150G --full_cleanup --min_contig_length 200 --SS_lib_type RF --__KMER_SIZE 31 2> {log}
        """

rule cd_hit_est:
        input:
                k25 = "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25.Trinity-GG.fasta",
                k31 = "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k31.Trinity-GG.fasta"
        output:
                mod_k25 = "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25.Trinity.mod.fasta",
                mod_k31 = "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k31.Trinity.mod.fasta",
                merged_mod = "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25_and_k31.Trinity.merged.mod.fasta",
                final_cd_hit_est = "MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25_and_k31.Trinity.merged.final.fasta"
        params:
                genotype="{genotype}"
        threads: 48
        resources:
                load=1
        log:
                "MyAssembly_{genotype}/logs/cd_hit_est_transcriptomes/{genotype}.log"
	conda:
		"/data/j/TRANSCRIPTOMATOR/cd-hit.yml"
	shell:
		"sed 's/>/>k25_{params.genotype}_/' {input.k25} > {output.mod_k25};"
                "sed 's/>/>k31_{params.genotype}_/' {input.k31} > {output.mod_k31};"
                "cat {output.mod_k25} {output.mod_k31} > {output.merged_mod};"
                "cd-hit-est -i {output.merged_mod} -o {output.final_cd_hit_est} -c 1 -n 11 -T {threads} -M 0 -d 0 -r 0 -g 1"

rule busco:
        input:
                transcriptome="MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25_and_k31.Trinity.merged.final.fasta"
        output:
                busco="MyAssembly_{genotype}/8_busco/My_Assembly_{genotype}_k25_and_k31_busco"
        resources:
                load=1
        threads: 48
        log:
                "MyAssembly_{genotype}/logs/busco/{genotype}.log"
	conda:
		"/data/j/TRANSCRIPTOMATOR/busco.yml"
	shell:
		"busco -i {input.transcriptome} -o {output.busco} -c {threads} -m transcriptome -l /data/j/TRANSCRIPTOMATOR/dbs/busco_db/embryophyta_odb10 > {log}"

rule transrate:
        input:
                transcriptome="MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25_and_k31.Trinity.merged.final.fasta",
                ref="rna.fna"
        output:
                transrate=directory("MyAssembly_{genotype}/9_transrate/assemblies.csv")
        resources:
                load=1
        threads: 48
        log:
                "MyAssembly_{genotype}/logs/transrate/{genotype}.log"
	conda:
		"/data/j/TRANSCRIPTOMATOR/transrate.yml"
	shell:
		"transrate --assembly {input.transcriptome} --reference {input.ref} --threads {threads} --output {output.transrate} 2> {log}"

#rule salmon_index:
#        input:
#                transcriptome="MyAssembly_{genotype}/7_trinity_assembly/MyAssembly_{genotype}_trinity_k25_and_k31.Trinity.merged.final.fasta"
#        output:
#                salmon_index=directory("MyAssembly_{genotype}/quant_salmon_index/")
#        resources:
#                load=1
#        threads: 40
#        log:
#                "logs/{genotype}_salmon_index.log"
#	conda:
#		"/data/j/TRANSCRIPTOMATOR/salmon.yml"
#	shell:
#		"salmon index -t {input.transcriptome} -p {threads} -i {output.salmon_index} 2> {log}"
#
#rule salmon_quant:
#        input:
#                salmon_index = directory("MyAssembly_{genotype}/quant_salmon_index/"),
#                filtered_R1 = expand("MyAssembly_{{genotype}}/6_contamination_removal/{sample}.trimmed.filtered.total.R1.fastq", sample=samples),
#                filtered_R2 = expand("MyAssembly_{{genotype}}/6_contamination_removal/{sample}.trimmed.filtered.total.R2.fastq", sample=samples),
#                unclassified_R1 = expand("MyAssembly_{{genotype}}/6_contamination_removal/{sample}.trimmed.unclassified.total.R1.fastq", sample=samples),
#                unclassified_R2 = expand("MyAssembly_{{genotype}}/6_contamination_removal/{sample}.trimmed.unclassified.total.R2.fastq", sample=samples)
#        output:
#                salmon_quant=directory("MyAssembly_{genotype}/10_salmon/quant/quant.sf")
#        params:
#                filtered_total_R1=' '.join(filtered_total_R1),
#                filtered_total_R2=' '.join(filtered_total_R2),
#                unclassified_total_R1=' '.join(unclassified_total_R1),
#                unclassified_total_R2=' '.join(unclassified_total_R2),
#        resources:
#                load=1
#        threads: 40
#        log:
#                "MyAssembly_{genotype}/logs/salmon_quant/{genotype}.log"
#	conda:
#		"/data/j/TRANSCRIPTOMATOR/salmon.yml"
#	shell:
#		"salmon quant -i {input.salmon_index} -l A -1 {params.filtered_total_R1} {params.unclassified_total_R1} -2 {params.filtered_total_R2} {params.unclassified_total_R2} -p {threads} --validateMappings -o {output.salmon_quant} 2> {log}"
#
