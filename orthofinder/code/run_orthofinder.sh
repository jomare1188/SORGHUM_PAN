#!/bin/bash
#PBS -q par128
#PBS -l nodes=1:ppn=128
#PBS -N orthofinder

cd $PBS_O_WORKDIR
mkdir -p ./../data

# copy and rename .pep of transcoder to a new folder with $GENOTYPE.pep name to run orthofinder
for i in $(cat ./../data/genotype_list)
do
	cp ./../../transdecoder/results/$i/longest_orfs.pep ./../data/$i.pep
done

conda activate ORTHOFINDER

mkdir -p ./../results
# Run orthofinder for many inflation values
for i in 1 1.3 1.5 1.8 2 2.5 3 4 5 6 7 10 15 20 30 40 50 100
do
	orthofinder -I ${i} -t 128 -a 64 -f ./../data/ -o ./../results/results_${i}
done


