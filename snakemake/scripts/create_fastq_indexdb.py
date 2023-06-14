#!/usr/bin/env python

'''				Description
		This script creates a indexdb for FASTQ files!
		This step is necessary to run remove_contaminants.py
'''

import argparse
import sys
import os
from Bio import SeqIO

# Creating arguments

parser = argparse.ArgumentParser(prog='create_indexdb.py', description='create indexdb for FASTQ files', add_help=True)
parser.add_argument('-R1', dest='R1_file', metavar='<R1 file>', help='R1 FASTQ file', required=True)
parser.add_argument('-R2', dest='R2_file', metavar='<R1 file>', help='R2 FASTQ file', required=True)
parser.add_argument('-o', dest='output_path', metavar='<output path>', help='path to output files')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')

# Getting arguments

args = parser.parse_args()
R1 = args.R1_file
R2 = args.R2_file
output_path = args.output_path
cwd = os.getcwd()

#print(R1)
#print(R2)
#teste/3_trimmed_reads/SRR4014643.trimmed.R1.fastq
#teste/3_trimmed_reads/SRR4014643.trimmed.R2.fastq

if not os.path.exists(output_path):
	os.makedirs(output_path)

# Generate index names 

index_R1 = os.path.basename(R1[:-5] + "index")
index_R2 = os.path.basename(R2[:-5] + "index")
#print(index_R1 + "basename + index")
#print(index_R2 + "basename + index")
#SRR4014643.trimmed.R1.index
#SRR4014643.trimmed.R2.index

index_R1_output = os.path.join(output_path, index_R1)
index_R2_output = os.path.join(output_path, index_R2)
#print(index_R1_output + "join output + basename + index")
#print(index_R2_output + "join output + basename + index")
#teste/6_contamination_removal/index/SRR4014643.trimmed.R1.index
#teste/6_contamination_removal/index/SRR4014643.trimmed.R2.index

#SeqIO.index_db(idx_name, file, format)
record_R1_index_db = SeqIO.index_db(index_R1_output, R1, "fastq")
record_R2_index_db = SeqIO.index_db(index_R2_output, R2, "fastq")

#record_R1_index_db = SeqIO.index_db(index_R1_output, R1, "fastq")
#record_R2_index_db = SeqIO.index_db(index_R2_output, R2, "fastq")
