samples=["BAZ9504", "BTx623", "cv.BRS330", "Della", "DKS-3707", "DKS-4420", "keller", "M35-1", "Mota", "R9188", "Rioref", "RTx430", "S085", "SM100", "TAM428", "SC187", "Tx2737", "Tx3362", "Tx378", "Tx7000"]
parts=expand("{sample}.part-{n}",sample=samples, n=range(1,11))

rule all:
    input:
        expand("/home/dmpachon/pfam_inflation_inflation/results/parts/{part}.pfam.tblout",
               part=parts)

rule hmmsearch_search:
    input:
        file="/home/dmpachon/pfam_inflation_inflation/data/parts/{part}.pep"
    output:
        tblout="/home/dmpachon/pfam_inflation_inflation/results/parts/{part}.pfam.tblout",
        pfamtblout="/home/dmpachon/pfam_inflation_inflation/results/parts/{part}.pfam.pfamtblout",
        out="/home/dmpachon/pfam_inflation_inflation/results/parts/{part}.pfam.out",
        domtblout="/home/dmpachon/pfam_inflation_inflation/results/parts/{part}.pfam.domtblout"
    threads: 6
    shell:
        """
        module load Bio/HMMER3/3.3.2
        hmmscan -o {output.out} --tblout {output.tblout} --pfamtblout {output.pfamtblout} --domtblout {output.domtblout} --cut_ga --cpu {threads} /home/dmpachon/pfam_inflation_inflation/data/Pfam-A.hmm {input}
        module unload Bio/HMMER3/3.3.2 
        """

