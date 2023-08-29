#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

source /home/jorge.munoz/.bashrc
conda activate snake
snakemake -p -s snake_transrate.py --keep-going --use-conda --rerun-incomplete --jobs 12  --cluster "qsub -q all.q -V -cwd -pe smp {threads} -N {params.jobname}"

