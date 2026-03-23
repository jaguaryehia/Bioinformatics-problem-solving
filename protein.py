from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
import pandas as pd

yeast_genome_file = "protein.faa"
yeast_ptt_file = "scc.csv"

"""
    this function read filename which is the fasta file 
        input:
            :param filename - fasta file : as a string
        output:
            :return recs : dictionary from fasta records
"""


def read_FASTA_File(filename):
    print("....Reading FASTA file...." + filename)
    recs = dict()
    records = SeqIO.parse(filename, "fasta")
    for i in records:
        recs[i.id] = i.seq
    return recs


"""
    this function is calculating the molecular weight 
        input:
            :param ID and seq - dictionary : as dictionary
        output:
            :return dictionary : replace the seq with the MW
"""


def calculate_molecular_weight(dictionary):
    for k, v in dictionary.items():  # Use .items() to iterate over key-value pairs
        dictionary[k] = molecular_weight(v, seq_type='protein')
    return dictionary


"""
    this function read the file csv with pandas to make it dataframe 
        input:
            :param filename - csv file : as a csv file
        output:
            :return dataframe : dataframe which is containing Protein product and Locus tag to compare them later 
"""


def parse_ptt_file(ptt_file):
    df = pd.read_csv(ptt_file)
    df2 = df[['Protein product', 'Locus tag']]
    return df2


"""
    this function comparing the dataframe and the dictionary that have th MW 
        input:
            :param  dict() - dictionary  : as dictionary which have the id for comparing
            :param dataframe - df : the dataframe which has the Protein product to compare
        output:
            :return dataframe : which has the protein product Locus tag MW
"""


def comparing(dictionary, df):
    for key in dictionary.keys():
        if key in df['Protein product'].values:
            df.loc[df['Protein product'] == key, 'MW'] = dictionary[key]
    df.to_csv('yeast.csv')
    return df


print(comparing(calculate_molecular_weight(read_FASTA_File(yeast_genome_file)), parse_ptt_file(yeast_ptt_file)))

