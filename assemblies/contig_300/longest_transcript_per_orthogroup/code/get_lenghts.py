from Bio import SeqIO

# Define the paths to the input files
fasta_file = "../data/all.cds_in_orthogroups.fasta"
ortho_file = "../data/Orthogroups_for_longest_cds_cogeqc2.tsv"

# Create a dictionary to store the sequences in the FASTA file
seq_dict = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    seq_dict[record.id] = len(record.seq)

# Open the orthogroup file for reading and a new file for writing
with open(ortho_file, "r") as input_file, open("../data/orthogroups_with_lengths.tsv", "w") as output_file:
    # Iterate over each line in the input file
    for line in input_file:
        # Remove any whitespace at the beginning and end of the line
        line = line.strip()
        # Split the line into columns
        columns = line.split("\t")
        # Get the ID of the sequence corresponding to this line
        seq_id = columns[1]
        # Get the length of the sequence from the dictionary
        seq_len = seq_dict.get(seq_id, "NA")
        # Add the sequence length as a new column to the line
        line_with_length = line + "\t" + str(seq_len) + "\n"
        # Write the line with the added column to the output file
        output_file.write(line_with_length)

