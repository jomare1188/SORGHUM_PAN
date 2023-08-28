#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --array=1-20
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6gb
#SBATCH --job-name=PANNZER2

GENOTYPE=`cat ../data/genotype_list.txt | awk -F "." '{print $1}' | head -n ${SLURM_ARRAY_TASK_ID} | tail -n1`
GENOTYPE_FILE=/home/dmpachon/pfam_inflation_inflation/data/${GENOTYPE}.pep
PANZZER=../SANSPANZ.3/runsanspanz.py
SPECIE="Saccharum officinarum x Saccharum spontaneum"

#/home/dmpachon/miniconda3/condabin/conda activate panzzer

source ~/.bashrc
conda activate panzzer
python ${PANZZER} -R -o ",../results/${GENOTYPE}_DE.out,../results/${GENOTYPE}_GO.out,../results/${GENOTYPE}_anno.out" -s "${SPECIE}" < ${GENOTYPE_FILE}



