from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import pairwise2



# python <dir>\assBioPython.py filter_nbases AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGNNNCA
# AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGCA
def Filter(dna):
    d = Seq(dna)
    FN = d.replace('N', '')
    return FN


# python <dir>\assBioPython.py calc_nbases AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGNNNCA
# 3
def nbases(dna):
    d = Seq(dna)
    nbases = d.count('N')
    return nbases


# python <dir>\assBioPython.py reverse_complement AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGNNNCA
def Revcomp(dna):
    d = Seq(dna)
    RC = d.reverse_complement()
    return RC


# python <dir>\assBioPython.py gc AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGNNNCA
# G count C count / len(seq) *
def gcpercentage(dna):
    z = GC(dna)
    return z


# python <dir>\assBioPython.py transcribe AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGNNNCA
def Transcribe(dna):
    d = Seq(dna)
    transcribe = d.transcribe()
    return transcribe


# python <dir>\assBioPython.py is_valid AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGCA dna
def valid(seq, typee):
    d = 'ACGT'
    rna = 'ACGU'
    protein = 'ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz'
    valid = True
    if typee == 'DNA' or typee == 'dna':
        for base in seq:
            if base not in d:
                valid = False
        return valid
    elif typee == 'RNA' or typee == 'rna':
        for base in seq:
            if base not in rna:
                valid = False
        return valid
    elif typee == 'PROTEIN' or typee == 'protein':
        for base in seq:
            if base not in protein:
                valid = False
        return valid
    print('The type incorrect')

# AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGCA
# AGCTGACTGACTACGTCGAGTCTA
# there is output file have the alll alignments
# python <dir>assBioPython.py seq_alignment AGCTGACTGACTACGTCGAGTCTACGTCGAGTCGTACGCA AGCTGACTGACTACGTCGAGTCTA -o output.txt
def seq_alignment(seq1, seq2, file_name=None):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    if file_name == None:
        for alignment in alignments:
            print(pairwise2.format_alignment(*alignment))
    else:
        f = open(file_name, 'a')
        for alignment in alignments:
            f.write(pairwise2.format_alignment(*alignment))
        f.close()
