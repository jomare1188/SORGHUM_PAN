# Load required library
library(seqinr)

# Specify file paths
fasta_file <- "/Storage/data1/jorge.munoz/SORGHUM_PAN/assemblies/contig_300/longest_transcript_per_orthogroup/data/all_transcripts.fa"
id_list_file <- "/Storage/data1/jorge.munoz/SORGHUM_PAN/assemblies/contig_300/longest_transcript_per_orthogroup/data/longest_cds_no_p.ids"
output_file <- "/Storage/data1/jorge.munoz/SORGHUM_PAN/assemblies/contig_300/longest_transcript_per_orthogroup/results/transcript_of_longest_cds_per_orthogroup.fasta"

# Read the list of IDs
id_list <- scan(id_list_file, what = character())

# Read the FASTA file
fasta_sequences <- read.fasta(fasta_file)

# Extract sequences based on the ID list
extracted_sequences <- fasta_sequences[names(fasta_sequences) %in% id_list]

# Write the extracted sequences to a new FASTA file
write.fasta(sequences = extracted_sequences, names = names(extracted_sequences), file.out = output_file)

# Print a message indicating successful completion
cat("Sequences extracted successfully and written to", output_file, "\n")

