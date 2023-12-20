import pandas as pd
import os

# import ORTHOGROUPS
ORTHOGROUPS = pd.read_csv("4orMoreOG_samples.txt", names = ["OG"])

rule all:
    input:
            #expand("/home/dmpachon/dNdS/results/{orthogroup}_mafft_alignment.fa",orthogroup=ORTHOGROUPS.OG),
            #expand("/home/dmpachon/dNdS/results/{orthogroup}_pal2nal.paml", orthogroup=ORTHOGROUPS.OG),
            #expand("/home/dmpachon/dNdS/results/{orthogroup}_omega.txt", orthogroup=ORTHOGROUPS.OG),
            #"/home/dmpachon/dNdS/results/dnds_sorghum_pan.csv",
            "/home/dmpachon/dNdS/results/plot_dnds.png"

rule mafft:
        input:
           "/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/WorkingDirectory/OrthoFinder/Results_I2.7/Orthogroup_Sequences/{orthogroup}.fa"
        output:
            "/home/dmpachon/dNdS/results/{orthogroup}_mafft_alignment.fa"
        threads: 2
        conda:
            "envs/mafft.yaml"
        resources:
                load = 2
        name: "mafft_alignment"
        params:
                opts="--maxiterate 100",
                mem="4gb"
        log:
                "logs/mafft_{orthogroup}.log"
        script:
            "scripts/run_mafft.sh"
rule pal2nal:
        input:
            mafft_alignmet="/home/dmpachon/dNdS/results/{orthogroup}_mafft_alignment.fa"
        output:
            pal2nal="/home/dmpachon/dNdS/results/{orthogroup}_pal2nal.paml",
            cds="/home/dmpachon/dNdS/data/{orthogroup}.cds.fa"
        threads: 1
        conda:
            "envs/mafft.yaml"
        resources:
            load = 1
        name: "pal2nal"
        params:
            opts="-output paml",
            orthogroup="{orthogroup}",
            mem="4gb"
        log:
            "logs/pal2nal_{orthogroup}.log"
        script:
            "scripts/run_pal2nal.sh"

rule codeml:
        input:
            pal2nal="/home/dmpachon/dNdS/results/{orthogroup}_pal2nal.paml"
        output:
            omega="/home/dmpachon/dNdS/results/{orthogroup}_omega.txt"
        threads: 1
        conda:
            "envs/mafft.yaml"
        resources:
            load = 1
        name: "codeml"
        params:
            orthogroup="{orthogroup}",
            mem="4gb"
        log:
            "logs/codeml_{orthogroup}.log"
        script:
            "scripts/run_codeml_model_minus2.sh"

rule merge_omegas:
        input:
            expand("/home/dmpachon/dNdS/results/{orthogroup}_omega.txt", orthogroup=ORTHOGROUPS.OG)
        output:
            "/home/dmpachon/dNdS/results/dnds_sorghum_pan.csv"
        threads: 1
        conda:
            "envs/mafft.yaml"
        resources:
            load = 1
        name: "dn.ds_merge_tables"
        params:
            mem="20gb"
        log:
            "logs/merge_tables.log"
        script:
            "scripts/run_merge_table.sh"
rule make_plots:
        input:
            "/home/dmpachon/dNdS/results/dnds_sorghum_pan.csv",
            "/home/dmpachon/dNdS/data/panTranscriptomeClassificationTable_I2.7.tsv"
        output:
            "/home/dmpachon/dNdS/results/plot_dnds.png",
            "/home/dmpachon/dNdS/results/dnds_panclass_table.csv"
        threads: 1
        conda:
            "envs/mafft.yaml"
        resources:
            load = 1
        name: "make_dnds_plots"
        params:
            mem="8gb"
        log:
            "logs/make_plot.log"
        script:
            "scripts/make_plots.R"
