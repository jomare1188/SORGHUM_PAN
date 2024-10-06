from Bio import SeqIO
import pandas as pd

# Ler a tabela de classificação original
classification_df = pd.read_csv('/Storage/data1/jorge.munoz/SORGHUM_PAN/pantranscriptome/stats/results/panTranscriptomeClassificationTable_I2.7.tsv', sep='\s+', header=None)
classification_df.columns = ['Group', 'Orthogroup', 'cds_id']
classification_df['transcript_id'] = classification_df['cds_id'].str.replace(r'\.p\d+', '', regex=True)

#print(classification_df.head())
#print(f"Número de colunas : {classification_df.shape[1]}")


# 2. dicionários para armazenar o comprimento dos transcritos e CDS
transcript_lengths = {}
cds_lengths = {}

# 3. Calcular o comprimento dos transcritos
with open('/Storage/data1/jorge.munoz/SORGHUM_PAN/assemblies/contig_300/longest_transcript_per_orthogroup/data/all_transcripts.fa') as transcripts_fasta:
    for record in SeqIO.parse(transcripts_fasta, "fasta"):
        transcript_lengths[record.id] = len(record.seq)

# 4. Calcular o comprimento das CDS
with open('/Storage/data1/jorge.munoz/SORGHUM_PAN/assemblies/contig_300/longest_transcript_per_orthogroup/data/all.cds.fasta') as cds_fasta:
    for record in SeqIO.parse(cds_fasta, "fasta"):
        cds_lengths[record.id] = len(record.seq)

# 5. DataFrame com Group, Orthogroup, ID, Transcript_Length e CDS_Length
new_df = classification_df[['Group', 'Orthogroup', 'cds_id', 'transcript_id']].copy()  # Copiar as 3 primeiras colunas
new_df['Transcript_Length'] = new_df['transcript_id'].map(transcript_lengths)  # Adicionar a coluna de comprimentos dos transcritos
new_df['CDS_Length'] = new_df['cds_id'].map(cds_lengths)  # Adicionar a coluna de comprimentos dos CDS

# 6. tabela com as colunas criadas
new_df.to_csv('Classification_lengths_test.tsv', sep='\t', index=False)
