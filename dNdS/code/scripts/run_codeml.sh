#!/usr/bin/env bash

orthogroups_dir=/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/WorkingDirectory/OrthoFinder/Results_I2.7/Orthogroup_Sequences/
tree=/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/WorkingDirectory/OrthoFinder/Results_I2.7/Gene_Trees/${snakemake_params[orthogroup]}_tree.txt

# prepare configuration file for CODEML
yes | cp template_CODEML.ctl ../data/${snakemake_params[orthogroup]}_CODEML.ctl
echo seqfile = ${snakemake_input[pal2nal]} >> ../data/${snakemake_params[orthogroup]}_CODEML.ctl
echo treefile = $tree >> ../data/${snakemake_params[orthogroup]}_CODEML.ctl
echo outfile = /home/dmpachon/dNdS/results/${snakemake_params[orthogroup]}_codeml.txt >> ../data/${snakemake_params[orthogroup]}_CODEML.ctl

# paml 4.10 needs in the first line in the tree the number of species and the number of trees (? always one in this caes)
## for some mysterious reason this does not work in snakemake, maybe bash stric mode
#if  ![[ $(head -n 1 $tree) =~ ^[0-9] ]]; then
#	head=$(echo $(grep -c ">" ../data/${snakemake_params[orthogroup]}.cds.fa))
#	sed -i "1s/^/$(echo $head 1)\n/g" $tree
#fi
##

# now we installed the last version of palm ;)
# actually run CODEML
mkdir -p /home/dmpachon/dNdS/results/${snakemake_params[orthogroup]}_codeml && cd /home/dmpachon/dNdS/results/${snakemake_params[orthogroup]}_codeml && /home/dmpachon/dNdS/data/paml4.9j/bin/codeml /home/dmpachon/dNdS/data/${snakemake_params[orthogroup]}_CODEML.ctl && cd -
# get omega
omega=$(grep omega ../results/${snakemake_params[orthogroup]}_codeml.txt | rev | cut -f1 -d" " | rev) > ${snakemake_output[omega]}
echo $omega,${snakemake_params[orthogroup]} > ${snakemake_output[omega]}

