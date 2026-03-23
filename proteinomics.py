from Bio.PDB import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




# Load the structure of the protein
parser = PDBParser()
structure = parser.get_structure("PROT", "6lu7.pdb")

# Get the active site residues of the protein
chain = structure[0]['A']
active_site = [res for res in chain if res.id[1] in active_site_residue_numbers]

# Get the sequence of the protein
seq = ""
for res in chain:
    seq += res.get_resname()

# Create a mutable sequence
seq = Seq(seq, generic_protein)
seq = seq.tomutable()

# Mutate the active site residues to all possible amino acids
amino_acids = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"]
for i in range(10):
    for aa in amino_acids:
        seq[active_site[i]] = aa
        record = SeqRecord(seq)
        print(f"Active Site Residue {i+1} mutated to {aa}  => {record.seq}")
