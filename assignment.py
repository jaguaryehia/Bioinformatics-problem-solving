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