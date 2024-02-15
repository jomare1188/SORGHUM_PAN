#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=60
#SBATCH --mem=150gb
#SBATCH --job-name=RTx430_longest

source ~/.bashrc
conda activate salmon
## Make index
threads=60
#transcriptome=../data/refs/RTx430/GCA_003482435.1_Corteva_Sorghum_ONT_TX430_1.0_genomic.fna
index=../data/indexes/Rio/longest_cds_per_orthogroups
res_dir=../results/RTx430/longest
mkdir -p $res_dir

## Make transcriptome index
#salmon index -t ${transcriptome} -p ${threads} -i ${index}

## Quantification with salmon
rm -f salmon_commands_longest_RTx430
for id in $(cat ../data/RTx430/ids_samples_RTx430.txt)
do
	echo salmon quant -i ${index} -l A -1 ../data/RTx430/paired_${id}_1.fastq -2 ../data/RTx430/paired_${id}_2.fastq -p 10 -o ${res_dir}/${id} >> salmon_commands_longest_RTx430
done
cat salmon_commands_longest_RTx430 | parallel -j ${threads}
