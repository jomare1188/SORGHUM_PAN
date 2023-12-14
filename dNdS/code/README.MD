# To run the pipeline in slurm you need to load R 4.3 and be sure to have the MAFFT.yaml

# CAUTION this pipe line will create a swarn of jobs, ~500k using 4orMoreOG_samples.txt, you should configure load and --jobs parameters seriously
# The 
module load HPC/R/4.3.1
snakemake -p -s Snakefile.py --resources load=50 --cluster "sbatch --partition=long --cpus-per-task=1 --mem={params.mem} --ntasks-per-node={threads} --job-name={name}" --use-conda --jobs 100 --rerun-incomplete --keep-going --scheduler greedy 
