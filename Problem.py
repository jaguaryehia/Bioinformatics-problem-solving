from collections import Counter
from Bio.Seq import Seq

"""
Problem 1:

    A string is simply an ordered collection of symbols selected from some alphabet and formed into a word;
    the length of a string is the number of symbols that it contains.
    An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is
    "ATGCTTCAGAAAGGTCTTACG."

    Given: A DNA string s of length at most 1000 nt.

    Return: Four integers (separated by spaces) counting the respective number of times
    that the symbols 'A', 'C', 'G', and 'T' occur in s

    Sample Dataset

            AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC

    Sample Output

            20 12 17 21
"""


def count_dna_symbols(dna_string):
    return Counter(Seq(dna_string))['A'], Counter(Seq(dna_string))['C'], Counter(Seq(dna_string))['G'], \
           Counter(Seq(dna_string))['T']


print(count_dna_symbols("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))

"""

    Problem 2:

        An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
        Given a DNA string t
        corresponding to a coding strand, its transcribed RNA string u is formed by replacing 
        all occurrences of 'T' in t with 'U' in u

        Given: A DNA string t having length at most 1000 nt.

        Return: The transcribed RNA string of t

        Sample Dataset

                GATGGAACTTGACTACGTAAATT

        Sample Output

                GAUGGAACUUGACUACGUAAAUU
"""


def transcribe_dna_to_rna(dna_string):
    # Convert the DNA string to a Biopython Seq object for easier manipulation
    seq = Seq(dna_string)

    # Transcribe the DNA sequence to RNA by replacing 'T' with 'U'
    rna_string = seq.transcribe()

    # Return the transcribed RNA string as a regular Python string
    return str(rna_string)


# Sample DNA string
sample_dna_string = "GATGGAACTTGACTACGTAAATT"

# Call the function with the sample DNA string
result = transcribe_dna_to_rna(sample_dna_string)

# Print the result
print(result)

"""

    Problem 3:

    In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.
    The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s
    , then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

        Given: A :DNA string s of length at most 1000 b
    p.

    Re
    turn: The reverse complement sc of s

    Sample Dataset

 
            AAAACCCGGT

    Sample Output

            ACCGGGTTTT
"""


def reverse_complement(dna_string):
    # Convert the DNA string to a Biopython Seq object for easier manipulation
    seq = Seq(dna_string)

    # Get the reverse complement of the DNA sequence
    reverse_complement_seq = seq.reverse_complement()

    # Return the reverse complement as a regular Python string
    return str(reverse_complement_seq)


# Sample DNA string
sample_dna_string = "AAAACCCGGT"

# Call the function with the sample DNA string
result = reverse_complement(sample_dna_string)

# Print the result
print(result)

"""
Problem 4:
    
    A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. 
    Sequences can be finite or infinite. Two examples are the 
    finite sequence (π,−2–√,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,…). 
    We use the notation an to represent the n-th term of a sequence.
    A recurrence relation is a way of defining the terms of a sequence with respect to the 
    values of previous terms. In the case of Fibonacci's rabbits from the introduction, 
    any given month will contain the rabbits that were alive the previous month, plus any new offspring.
    A key observation is that the number of offspring in any month is equal to the number of rabbits that were 
    alive two months prior. As a result, if Fn
    represents the number of rabbit pairs alive after the n-th month, then we obtain 
    the Fibonacci sequence having terms Fn that are defined by the recurrence relation 
    Fn=Fn−1+Fn−2 (with F1=F2=1 to initiate the sequence). 
    Although the sequence bears Fibonacci's name, 
    it was known to Indian mathematicians over two millennia ago.
    When finding the n-th term of a sequence defined by a recurrence relation, 
    we can simply use the recurrence relation to generate terms for progressively larger values of n
    This problem introduces us to the computational technique of dynamic programming, 
    which successively builds up solutions by using the answers to smaller cases.

    Given: Positive integers n≤40 and k≤5

    Return: The total number of rabbit pairs that will be present after n
    months, if we begin with 1 pair and in each generation, every pair of 
    reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
    
    Sample Dataset
    
            5 3
    
    Sample Output
    
            19
"""

"""
count RNA Nucleotides

"""


def count_rna_nucleotides(rna_string):
    # Convert the RNA string to a Biopython Seq object for easier manipulation
    seq = Seq(rna_string)

    # Use Counter to count the occurrences of each RNA nucleotide in the sequence
    nucleotide_counts = Counter(seq)

    # Return the counts of 'A', 'C', 'G', and 'U' as integers
    return nucleotide_counts['A'], nucleotide_counts['C'], nucleotide_counts['G'], nucleotide_counts['U']


# Sample RNA string
sample_rna_string = "GAUGGAACUUGACUACGUAAAUU"

# Call the function with the sample RNA string
result = count_rna_nucleotides(sample_rna_string)

# Print the result in the required format
print(result)

"""
Problem 6:
    GC Content
"""

from Bio.SeqUtils import GC


def calculate_gc_content(dna_sequence):
    return GC(dna_sequence)


"""
    Problem 6

    Given two strings s and t of equal length, the Hamming distance between s and t,
    denoted dH(s,t), is the number 
    of corresponding symbols that differ in s and t
    See Figure 2.

    Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

    Return: The Hamming distance dH(s,t).
    Sample Dataset
        GAGCCTACTAACGGGAT
        CATCGTAATGACGGCCT
    Sample Output
        7
"""


def hamming_distance(s, t):
    return sum(a != b for a, b in zip(s, t))


# Sample DNA strings
s = "GAGCCTACTAACGGGAT"
t = "CATCGTAATGACGGCCT"

# Call the function with the sample DNA strings
result = hamming_distance(s, t)

# Print the result
print(result)

"""
The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

Given: An RNA string s

corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The protein string encoded by s

.
Sample Dataset

AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

Sample Output

MAMAPRTEINSTRING

"""

from Bio.Seq import Seq


def translate_rna_to_protein(rna_string):
    # Convert the RNA string to a Biopython Seq object for easier manipulation
    seq = Seq(rna_string)

    # Translate the RNA sequence into a protein sequence using the RNA codon table
    protein_string = seq.translate()

    # Return the protein string as a regular Python string
    return str(protein_string)


# Sample RNA string
sample_rna_string = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"

# Call the function with the sample RNA string
result = translate_rna_to_protein(sample_rna_string)

# Print the result
print(result)


def total_subsets(n):
    # Calculate the total number of subsets modulo 1,000,000
    total_subsets = 2 ** n % 1_000_000
    return total_subsets

# Sample input value
n = 3

# Call the function with the sample input
result = total_subsets(n)

# Print the result
print(result)


"""
Edit Distance Alignment
url: http://rosalind.info/problems/edta/

Given: Two protein strings s and t in FASTA format (with each string having length at most 1000 aa).
Return: The edit distance dE(s,t) followed by two augmented strings s′ and t′ representing an optimal alignment of s and t.
"""

import pprint
from Bio import SeqIO

def EditDistanceAlignment(s, t):
    m, n = len(s), len(t)
    if m*n==0:
        return m+n
    DP = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1):
        DP[i][0] = i
    for j in range(n+1):
        DP[0][j] = j
    for i in range(1, m+1):
        for j in range(1, n+1):
            left = DP[i-1][j] + 1
            down = DP[i][j-1] + 1
            left_down = DP[i-1][j-1]
            if s[i-1] != t[j-1]:
                left_down += 1
            DP[i][j] = min(left, down, left_down)
    # pprint.pprint(DP)
    edit_distance = DP[m][n]
    print(edit_distance)

    s_, t_ = "", ""
    i, j = 0, 0
    i, j = len(s), len(t)
    while (i>0 and j>0):
        left = DP[i][j-1]
        top = DP[i-1][j]
        left_top = DP[i-1][j-1]
        min_ = min(left, top, left_top)
        if DP[i][j]==min_:
            s_ = s[i-1]+s_
            t_ = t[j-1]+t_
            i -= 1
            j -= 1
        else:
            if (min_==left and min_==top) or (min_!=left and min_!=top):
                s_ = s[i-1]+s_
                t_ = t[j-1]+t_
                i -= 1
                j -= 1
            elif min_!=left and min_==top:
                s_ = s[i-1]+s_
                t_ = "-"+t_
                i -= 1
            elif min_==left and min_!=top:
                s_ = "-"+s_
                t_ = t[j-1]+t_
                j -= 1
    print(s_)
    print(t_)
    return DP





seq_name, seq_string = [], []
with open("../data/rosalind_edta.txt", "r") as fa:
        for seq_record in SeqIO.parse(fa, "fasta"):
            seq_name.append(seq_record.name)
            seq_string.append(str(seq_record.seq))
        s, t = seq_string
res = EditDistanceAlignment(s, t)

from Bio import SeqIO
from Bio.SeqUtils import GC

def read_fasta(file_path):
    """
    Read a FASTA file and return the sequences as a list of tuples (sequence_id, sequence).

    Parameters:
        file_path (str): The path to the FASTA file.

    Returns:
        list: List of tuples (sequence_id, sequence).
    """
    sequences = []
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append((record.id, str(record.seq)))
    return sequences

def fastaqc(sequences):
    """
    Perform FASTAqc on a list of sequences.

    Parameters:
        sequences (list): List of tuples (sequence_id, sequence).

    Returns:
        dict: Dictionary containing sequence statistics.
    """
    stats = {}
    for seq_id, sequence in sequences:
        stats[seq_id] = {
            "Length": len(sequence),
            "GC Content": GC(sequence)
        }
    return stats

# Example usage
fasta_file_path = "example.fasta"
sequences = read_fasta(fasta_file_path)
qc_results = fastaqc(sequences)

# Print FASTAqc results
for seq_id, stats in qc_results.items():
    print(f"Sequence ID: {seq_id}")
    print(f"Length: {stats['Length']}")
    print(f"GC Content: {stats['GC Content']:.2f}%")
    print("--------------")


