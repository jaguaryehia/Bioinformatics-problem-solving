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