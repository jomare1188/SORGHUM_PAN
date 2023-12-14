#!/usr/bin/env bash

software=/home/dmpachon/dNdS/data/pal2nal.v14/pal2nal.pl
orthogroups_dir=/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/WorkingDirectory/OrthoFinder/Results_I2.7/Orthogroup_Sequences/

## get cds files
grep ">" ${orthogroups_dir}/${snakemake_params[orthogroup]}.fa | sed 's/>//g' > ../data/${snakemake_params[orthogroup]}.ids
grep --no-group-separator -A1 -w -f ../data/${snakemake_params[orthogroup]}.ids ../data/cds_pergenotype/all.cds.fasta > ${snakemake_output[cds]}
echo total cds $(wc -l ../data/${snakemake_params[orthogroup]}.ids | cut -f1 -d" ")    found cds $(grep -c ">" ../data/${snakemake_params[orthogroup]}.cds.fa) for ${snakemake_params[orthogroup]}.fa

## run pal2nal
${software}  "${snakemake_input[mafft_alignmet]}" "${snakemake_output[cds]}" ${snakemake_params[opts]}  > "${snakemake_output[pal2nal]}" 2>> "${snakemake_log[0]}"
