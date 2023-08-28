#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 8
rm run_me
for transcriptome in $(cat not_conting_300.list)
do	
	/Storage/data1/jorge.munoz/SCPT/snakemakePipeline/extractSequencesWithMinLength.pl -f $transcriptome -m 301 >> run_me
done

parallel -j8 -a run_me
