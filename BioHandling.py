"""
    pipline for making Big Biopython library for bioinformatics:
        1. Biopython library:
            1. handling errors and parsing files
            2. functions that can be edited and make it more optimized with C++
            3. see SWs that has more functionalities for it
            4. make ML for analysis and predicting the blank datasets or predicting the results

        2. transcriptomics for bioinformatics:
            1. see how is R working with transcriptomics and the researches for it
            2. make class for it and make all functions for transcriptomics analysis
            3. see SWs that has more functionalities for it
            4. make sure of handling and more optimized with C++
            5. make ML for analysis and predicting the blank datasets or predicting the results

        3. proteinomics for bioinformatics:
            1. see how is R working with proteinomics and the structure of the protein analysis
            2. make class for it and make all functions for proteinomics analysis
            3. see SWs that has more functionalities for it
            4. make sure of handling and more optimized with C++
            5. make ML for analysis and predicting the blank datasets or predicting the results

        4. metabolomics for bioinformatics:
            1. see how is R working with proteinomics and the structure of the metadata analysis
            2. make class for it and make all functions for metabolomics analysis
            3. see SWs that has more functionalities for it
            4. make sure of handling and more optimized with C++
            5. make ML for analysis and predicting the blank datasets or predicting the results

        5. genomics for bioinformatics:
            1. see how is R working with genomics and the researches for it
            2. make class for it and make all functions for genomics analysis
            3. see SWs that has more functionalities for it
            4. make sure of handling and more optimized with C++
            5. make ML for analysis and predicting the blank datasets or predicting the results
"""
import gzip
from collections import Counter
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt


def getGCContent(seq):
    return GC(seq)


def dna_Reverse_Complement(dna):
    return Seq(dna).upper().reverse_complement_rna()


def dna_Complement(dna):
    return Seq(dna).upper().complement()


def dna_Translation(dna):
    return Seq(dna).upper().translate()


def dna_transcription(seq):
    return Seq(seq).upper().transcribe()


def dna_back_Transcribe(dna):
    return Seq(dna).upper().back_transcribe()


def fastq_N_uncalled_base(recs):
    n_cnt = Counter()
    for rec in recs:
        for i, letter in enumerate(rec.seq):
            pos = i + 1
            if letter == 'N':
                n_cnt[pos] += 1
    seq_len = max(n_cnt.keys())
    positions = range(1, seq_len + 1)
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.plot(positions, [n_cnt[x] for x in positions])
    ax.set_xlim(1, seq_len)


def get_FASTQ_info(recs):
    rec = next(recs)
    print(rec)
    print(rec.id, rec.description, rec.seq)
    print(rec.letter_annotations)
    count = None
    for rec in recs:
        count = Counter(rec.seq)
    total = sum(count.values())
    for letter, count in count.items():
        print('%s: %.2f %d' % (letter, 100. * count / total, count))


def read_FASTQ_File(filename):
    print(f'Reading FASTQ...{filename}')
    if filename.endswith('.fastq.gz'):
        try:
            f = SeqIO.parse(gzip.open(filename, 'rt', encoding='utf-8'), 'fastq')
        except OSError:
            f = SeqIO.parse(filename, 'fastq')
    elif filename.endswith('.fastq'):
        f = SeqIO.parse(filename, 'fastq')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))
    return f


# def get_FASTA_info(recs):


def read_FASTA_File(filename):
    print(f'Reading FASTA...{filename}')
    if filename.endswith('.fasta.gz'):
        try:
            f = SeqIO.parse(gzip.open(filename, 'rt', encoding='utf-8'), 'fasta')
        except OSError:
            f = SeqIO.parse(filename, 'fasta')
    elif filename.endswith('.fasta'):
        f = SeqIO.parse(filename, 'fasta')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))
    return f

# def alignmentFromFasta(fastaFileName, seq):


