#SARS-coV-2 complete genome data analysis

from Bio import SeqIO


#Reading an loading FASTA format
for record in SeqIO.parse("sequence.fasta","fasta"):
    print(record)
#print(record.id)
#print(record.description)

#Reading sequence in a FASTA file
dna_record = SeqIO.read("sequence.fasta","fasta")
#print(dna_record)
dna_seq = dna_record.seq
#print(dna_seq)

#Reading GenBank file same as FASTA file
#print("From Genbank file")
for record in SeqIO.parse("sequence.gb","genbank"):
    print(record)
#print(record.id)
#print(record.description)

#Reading sequence in a genbank file
dna_record = SeqIO.read("sequence.gb","gb")
#print(dna_record)
dna_seq = dna_record.seq
#print(dna_seq)

