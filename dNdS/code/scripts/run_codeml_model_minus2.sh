#!/usr/bin/env bash

orthogroups_dir=/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/WorkingDirectory/OrthoFinder/Results_I2.7/Orthogroup_Sequences/
tree=/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/WorkingDirectory/OrthoFinder/Results_I2.7/Gene_Trees/${snakemake_params[orthogroup]}_tree.txt

# prepare configuration file for CODEML
yes | cp template_CODEML.ctl ../data/${snakemake_params[orthogroup]}_CODEML.ctl
echo seqfile = ${snakemake_input[pal2nal]} >> ../data/${snakemake_params[orthogroup]}_CODEML.ctl
echo outfile = /home/dmpachon/dNdS/results/${snakemake_params[orthogroup]}_codeml.txt >> ../data/${snakemake_params[orthogroup]}_CODEML.ctl

# Actually run CODEML
mkdir -p /home/dmpachon/dNdS/results/${snakemake_params[orthogroup]}_codeml && cd /home/dmpachon/dNdS/results/${snakemake_params[orthogroup]}_codeml && /home/dmpachon/dNdS/data/paml4.9j/bin/codeml /home/dmpachon/dNdS/data/${snakemake_params[orthogroup]}_CODEML.ctl && cd -
# get omega  code
# get average of pairwise omega
omega=$(awk -f scripts/filter.awk /Storage/data1/hellen.silva/db-extraction/dNdS/results/${snakemake_params[orthogroup]}_codeml/rst) > ${snakemake_output[omega]}
echo $omega,${snakemake_params[orthogroup]} > ${snakemake_output[omega]}

