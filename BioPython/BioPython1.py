import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC 

#https://www.youtube.com/watch?v=t0Z0FP8t0_4

'''1. Creates coding strand object for coding DNA strand
   2. Transcription of coding DNA strand to mRNA
   3. Translation of mRNA to amino acid sequence

   Also don't forget DNA is double stranded where the other strand is
   the reverse complement of the other. Whereas RNA is single stranded.'''

#DNA strand: ATGTTACACTCCCGATGA

#Sequence object stored in cdna variable
dna = Seq("ATGTTACACTCCCGATGA",IUPAC.unambiguous_dna)
print(dna)

mrna = dna.transcribe() #This function only works on coding DNA strand,
print(mrna)

protein = mrna.translate()
print(protein) # *-indicates a stop codon

#but in reality transcribe method only works with a template strand, so can convert it into a coding strand

#template_dna variable replaced with the template sequence and use reverse_complement method for transcription
#coding_dna = template_dna.reverse_complement()

cdna2 = Seq("ATGCTGTTCGCCCTATAGTGTCTAAGCTAG", IUPAC.unambiguous_dna)
mrna2 = cdna2.transcribe()
print(mrna2)

protein2 = mrna2.translate()
print(protein2)

#first in frame stop codon indicated by to_stop parameter, which does not translate stop codon itself
protein2 = mrna2.translate(to_stop=True)
print(protein2)
