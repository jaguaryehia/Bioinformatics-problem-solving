from Bio import SeqIO
from Bio.Seq import *
from Bio.SeqUtils import *
from Bio.PDB import *




# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-i', '--input', required=True)
#     parser.add_argument('-o', '--output', required=True)
#     args = parser.parse_args()
#     input_file = args.input
#     output_file = args.output
#     with open(output_file, 'w') as output_handle:
#         with open(input_file, 'r') as input_handle:
#             for record in SeqIO.parse(input_handle, 'fasta'):
#                 output_handle.write('>{}\n{}\n'.format(record.id, record.seq))
#                 output_handle.flush()

def GC_content(sequence):
    sq = Seq(sequence)
    content = GC(sq)
    return content


def transcription(sequence):
    coding_dna = Seq(sequence)
    messenger_rna = coding_dna.transcribe()
    return messenger_rna


def translation_from_dna(sequence):
    seque = Seq(sequence)
    protein = seque.translate()
    return protein


def translation_from_mrna(messenger_rna):
    seque = Seq(messenger_rna)
    protein = seque.translate()
    return protein


def reverse_complement(sequence):
    seq = Seq(sequence)
    seq = seq.reverse_complement()
    return seq


def back_transcription(messenger_rna):
    seque = Seq(messenger_rna)
    coding_dna = seque.back_transcribe()
    return coding_dna


