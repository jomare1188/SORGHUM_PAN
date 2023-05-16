#!/bin/bash
#PBS -q par16
#PBS -l nodes=1:ppn=16
#PBS -N orthofinder

cd $PBS_O_WORKDIR
# get genotype list in transcoder folder
ls ./../../transdecoder/results/ | grep -v -w genotype_list > ./../../transdecoder/results/genotype_list
mkdir -p ./../data
# copy and rename .pep of transcoder to a new folder with $GENOTYPE.pep name to run orthofinder
for i in $(cat ./../../transdecoder/results/genotype_list)
do
	cp ./../../transdecoder/results/$i/longest_orfs.pep ./../data/$i.pep
done
conda activate ORTHOFINDER

# Run orthofinder for only one inflation value
#orthofinder -I 1.8 -t 16 -a 8 -f ./../data/ -o ./../results

# Run orthofinder for many inflation values
for i in {1 1.5 1.8 2 2.5 3 4 5 6 7 10 15 20 30 40 50 100}
do
	orthofinder -I ${i} -t 16 -a 8 -f ./../data/ -o ./../results_${i}
done


