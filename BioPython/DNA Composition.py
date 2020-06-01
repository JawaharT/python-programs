from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC123, GC_skew, xGC_skew #xGC_skew needs tkinter
from Bio.SeqUtils import nt_search

dna = Seq("ATGATCTCGTAA")
print(dna)

#GC Content
print(GC(dna))

#Method 2
print(dna.count("G"))
def count_GC_2(dna):
    return float(dna.upper().count("G") + dna.upper().count("C"))/len(dna) * 100
print(count_GC_2(dna))

#Method 3
def count_GC_3(dna):
    return float(len([base for base in dna.upper() if base in "GC"])/len(dna)) * 100
print(count_GC_3(dna))

#No biopython function for AT content

def at_content(dna):
    return float(dna.upper().count("A") + dna.upper().count("T"))/len(dna) * 100
print(at_content(dna))

#Merge both functions into 1 function
def content(dna, at_or_gc):
    return float(dna.upper().count(at_or_gc[0]) + dna.upper().count(at_or_gc[1]))/len(dna) * 100
print("AT Content: {0}%, GC Content {1}% in example sequence".format(content(dna, "AT"), content(dna, "GC")))

#Can use Tm_Wallace: 'Rule of thumb', Tm_GC or Tm_NN parameters to find melting point
#Use GC Content to measure melting point

#Get melting point using Wallace
print(mt.Tm_Wallace(dna))

#Get melting point using GC content
print(mt.Tm_GC(dna))

#Get melting using GC content
def get_metrics(dna):
    return "GC: {0}, AT: {1}, Melting Point {2}".format(GC(dna), at_content(dna), mt.Tm_GC(dna))
print(get_metrics(dna))

#GC Skew
'''Check when G and C are over or under abundant in a particular region of DNA or RNA.
Helps identify DNA leading strand
GC skew pos = leading, GC skew neg = nagging'''

#GC content overall, 1st, 2nd and 3rd positions returned
print(GC123(dna))

#GC_skew, default: checks 1st 100 bases, 0 means no over or under abundant GC bases
print(GC_skew("ATGGGGTCCCGCTC"))

#Subsequences: search for smaller DNA within larger DNA returns substring and index it was found in
main = Seq("ACTATT")
sub = Seq("ACT")
print(nt_search(str(main), str(sub)))
print(nt_search(str(main), "ATT"))

