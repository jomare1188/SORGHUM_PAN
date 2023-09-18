#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 45
#$ -l hostname=neotera

mkdir -p ./../data

# copy and rename .pep of transcoder to a new folder with $GENOTYPE.pep name to run orthofinder
#for i in $(cat ./../data/contig_300/genotypes)
#do
#	cp -v ./../../transdecoder/results/${i}.pep/longest_orfs.pep ./../data/contig_300/$i.pep
#done

source /home/jorge.munoz/.bashrc
conda activate ORTHOFINDER

mkdir -p ./../results

# Run orthofinder to create diamond db with ultrasensitive options (need to configure scripts_of/config.json) and for i=2 

orthofinder -S diamond_ultra_sens3 -I 2 -t 45 -a 30 -f ./../data/contig_300 -o ./../results/I2

# Run orthofinder for many inflation values using the previuos diamond results
#for i in 1.1 1.3 1.5 1.8 2 2.5 3 4 5 6 7 10 15 20 30 40 50 100
#do
#	orthofinder -t 128 -a 128 -S diamond_ultra_sens3 -I ${i} -X -b ./../results/results_2/Results_*/WorkingDirectory/ 

#done

