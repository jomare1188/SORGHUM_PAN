# Pipeline to used to mount transcriptomes in Sorghum pan-transcriptome
 There are three snakefiles each one is used for a different process: Snakefile_download.py recievies the the genotype and accessions information and download the data and check the strandness of each accession, Snakefile_denovo.py and Snakefile_reference.py recieves the genotype and accession information and mount the transcriptome and give some quality metrics. Snakefile_reference.py uses trinty reference guided strategy to mount transcriptome, Snakefile_denovo.py uses de novo trinity mount strategy.
# Snakefile_download.py
## Sample file
 you need to specify a sample file: in the snakefile line 26 containing in the first column the genotype and separated by a tab in the second column the accession:
 |genotype |       sample|
 | ------ | ------ |
 |Rio   |  SRR6188852|
 |Rio|     SRR6257132|
 |Rio |   SRR6257135|
 |Rio |   SRR6257137|
 |Rio |   SRR6257138|

you need to configure the rule filter_stranded lines 104 and 105 to coincide with the name the sample file and genotype.
 ## conda envs
 you can find the yml configuration files in conda_envs directory
 - download.yml: python2.7, ffq 
 - salmon.yml: salmon
 - jq.yml: jq
 ## output 
 the pipeline is going to give you two files
 In this case stranded_status_Rio.tsv which contains the library type given by salmon.
 |genotype|        sample | lib|
 | --- | ----------------|---- |
 |Rio|     SRR6188852    |  ISF|
 |Rio|     SRR6257132    | ISR|
 |Rio|     SRR6257135    |  IU|
 |Rio|    SRR6257137     | ISR|
 |Rio|     SRR6257138    |  IU|

and stranded_samples_Rio.tsv with the information of only stranded samples (the file given by the pipeline doesn't have heaeder).

 |genotype|        sample | lib|
 | --- | ----------------|---- |
 |Rio|     SRR6188852    |  ISF|
 |Rio|     SRR6257132    | ISR|
 |Rio|    SRR6257137     | ISR|

 we need the information of the stranded accession for the nexts steps in csv format, we can extract this information from the stranded_samples_Rio.tsv using 
 ```sh
 cut -f1 stranded_samples_Rio.tsv | sed ':a;N;$!ba;s/\n/,/g' > samples_Rio2.csv
 ```
 ## To run the pipelile
 you shuld personalize configuration directly in Snakemake_download.py, it is important to set properly jobs and cores limit acordingly to your resources
 ```sh
 snakemake -p -s Snakemake_download.py --rerun-incomplete --jobs 48 --keep-going --use-conda --resources 100 --cores 48 --latency-wait 60
 ```
 since this pipeline depends on data hosted by third parties is usuall that downloads breaks or simply don't work. Each time the pipelines breaks you can restart from  the point it broke with the same command used to run the pipeline.


