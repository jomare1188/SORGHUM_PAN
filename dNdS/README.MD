# Here you will find a Snakemake pipeline to get dN/dS from Sorghum pan-transcriptome data

## Requirements 
### Software
* Snakemake version xx
* MAFFT version xx
* pal2nal version xx
* CODEML version xx
* R version xx

### Files
* Sample file with a list of orthgroups to be analyzed (the orthogroups follow this pattern OGXXXXXXX)
* CDS fasta for pan-transcriptome ( ./../data/cds_pergenotype/all.cds.fasta)
* Proteins orthogroups fasta (from orthofinder: OGXXXXXXX.fa)
* Pan-transcriptome classification table (in this case ../../sorghumxxxxx/xxxx/panTranscriptomeClassificationTable_I2.7.tsv )
* Template configuration file for CODEML: template_CODEML.ctl

## Overview

This repository aims to calculate the Non-synonymous to Synonymous (dN/dS) mutation ratio across various elements of the sorghum pan-transcriptome. The process involves several key steps:

1) Protein Orthogroups Identification: We begin by utilizing the protein orthogroups identified in the sorghum pan-transcriptome.

2) Protein Alignment (MAFFT): For each orthogroup, we perform protein alignment using MAFFT.

3) Nucleotide Alignment (pal2nal): Using the CDS sequences and protein ids in each orthogroup we build a cds fasta for each orthogrouo (ensuring identical IDs between CDS and orthogroups), with the protein alignment and the correspondant cds fasta we generate a protein-informed nucleotide alignment through the pal2nal tool.

4) dN/dS Calculation (CODEML): Employing CODEML with runmode = -2, we calculate the dN/dS ratio from the informed nucleotide alignment.

5) Data Aggregation: We aggregate the results into a table, providing dN/dS values for each orthogroup and classifying them into pan-transcriptome categories (hard-core, soft-core, exclusive, and accessory).

6) Result Visualization: Finally, an R script is employed to generate visualizations, enhancing the interpretability of the calculated dN/dS ratios.

Feel free to explore the code and adapt the workflow based on your specific needs. Contributions and suggestions are welcome!
