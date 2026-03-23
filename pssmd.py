import numpy as np
from math import log2

def read_data(file_name):
    with open(file_name, 'r') as file:
        num_seq, seq_len = map(int, file.readline().split())
        seqs = [file.readline().strip() for _ in range(num_seq)]
    return seqs


def print_seqs(seqs):
    print("Aligned Sequences:")
    for seq in seqs:
        print(seq)


def calc_pssm(sequences):

    t = len(sequences)
    n = len(sequences[0])
    total_no_pos = t*n

    pssm_matrix = [[0] * 4 for _ in range(n)] 

    for i in range(n):
        column = [sequence[i] for sequence in sequences]
        for j in range(4):
            count = column.count('ACGT'[j])
            pssm_matrix[i][j] = count / t

    pssm_matrix = list(zip(*pssm_matrix))
    for index , row in enumerate (pssm_matrix):
        pssm_matrix[index] = [log2(val / ((sum(row)*t) / total_no_pos)) if val != 0 else 0 for val in row ]

    return pssm_matrix



def print_pssm(pssm):
    print("\nPSSM Matrix:")
    for row in pssm:
        print(' '.join(f'{value:.2f}' for value in row))


def calc_prob(seq, pssm):
    prob = 0
    seq = seq.upper()
    for i, nuc in enumerate(seq):
        if nuc in 'ACGT':
            index = 'ACGT'.index(nuc)
            prob += pssm[index][i]
    return round(prob, 2)


file_name = 'PSSMData.txt'
seqs = read_data(file_name)

print_seqs(seqs)
pssm = calc_pssm(seqs)
print_pssm(pssm)

new_seq = input("Enter a new sequence of length {}: ".format(len(seqs[0])))
prob = calc_prob(new_seq, pssm)
print("Probability of the new sequence joining the rest of the aligned DNA sequences:", prob)
