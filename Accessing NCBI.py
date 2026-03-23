from Bio import Entrez, SeqIO

Entrez.email = ""

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
