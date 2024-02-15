#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -N snake_master
#$ -pe smp 1

module load miniconda3
conda activate snake
snakemake -p -s Snakefile_download.py --resources load=100 --cluster "qsub -q all.q -V -cwd -l h={params.server} -pe smp {threads} -N {params.jobname}" --jobs 10 --keep-going --conda-frontend mamba --use-conda

