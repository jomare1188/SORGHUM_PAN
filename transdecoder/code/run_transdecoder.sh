#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 23


module load transdecoder/5.5.0

#ls ./../../assemblies/ | grep .fasta > ./../../assemblies/transcriptomes_list
mkdir -p ./../results

rm transdecoder_commands
for i in $(cat ./../../assemblies/contig_300/genotypes)
do
	echo TransDecoder.LongOrfs -S -t ./../../assemblies/contig_300/${i}.fasta --output_dir ./../results/${i}.pep >> transdecoder_commands
done

cat transdecoder_commands | parallel -j23
