from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
import csv



yeast_genome_file = "yeast/yeast_genome.fasta.fasta"
yeast_ptt_file = "yeast/yeast_protein.ptt.ptt"





def extract_protein_sequences(yeast_genome):
    protein_sequences = {}
    for record in list(SeqIO.parse(yeast_genome, "fasta")):
        protein_id = record.id
        protein_sequence = str(record.seq.translate())
        protein_sequences[protein_id] = protein_sequence
    return protein_sequences


print(extract_protein_sequences(yeast_genome_file))

def parse_ptt_file(ptt_file):
    locus_tags = {}
    with open(ptt_file) as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            fields = line.strip().split("\t")
            locus_tag = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            locus_tags[locus_tag] = (start, end)
    return locus_tags






def calculate_molecular_weight(protein_sequence):
    return molecular_weight(protein_sequence)


# In[35]:
#
#
# with open("protein_molecular_weights.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(["Locus Tag", "Molecular Weight"])
#
#
#     for locus_tag, protein_sequence in zip(yeast_locus_tags, yeast_protein_sequences):
#         mw = calculate_molecular_weight(protein_sequence)
#         writer.writerow([locus_tag, mw])
#
#
#     for locus_tag, protein_sequence in zip(ecoli_locus_tags, ecoli_protein_sequences):
#         mw = calculate_molecular_weight(protein_sequence)
#         writer.writerow([locus_tag, mw])
#
#
# # In[ ]:
#



