#!/usr/bin/env bash

sed -i 's/*//g' "${snakemake_input[0]}"
mafft ${snakemake_params[opts]} --thread ${snakemake[threads]}  "${snakemake_input[0]}"  > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"
