from Bio.Seq import Seq
from Bio.SeqUtils import seq3, seq1 #Converts aa to 3 letter version
from Bio.Data import CodonTable

dna_seq = Seq("ATGATCTCGTAA")
print(dna_seq, len(dna_seq))

#Complementary DNA sequence
print(dna_seq.complement())

#A T has 2 hydrogen bonds
#G C has 3 hydrogen bonds

#Reverse Complement DNA sequence, same as complementary sequence but in reverse
print(dna_seq.reverse_complement())

#mRNA
mrna = dna_seq.transcribe()
print(mrna)

#Method 1: Protein/Amino Acids
aa = mrna.translate()
print(aa)

#Method 2: Direct Translation
aa = dna_seq.translate()
print(aa)

#Custom stop codon
print(mrna.translate(stop_symbol="@"))

#Back transcribe to DNA, same as dna_seq
print(mrna.back_transcribe())

#Join steps
print(dna_seq.transcribe().translate())

#3 letter version of protein
three_letter_aa = seq3(str(aa))
print(three_letter_aa)

#Back to 1 letter version of protein same as original aa
one_letter_aa = seq1(str(three_letter_aa))
print(one_letter_aa)

#Methods on Bio.Data
print(dir(CodonTable))

#CodonTable for DNA by name
print(CodonTable.unambiguous_dna_by_name["Standard"])

#CodonTable for RNA by name
print(CodonTable.unambiguous_rna_by_name["Standard"])
