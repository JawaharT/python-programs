'''Similarity refers to resemblance of 2 sequences in comparison,
minimal number of edit operations (inserts, deletes, substitutions)
in order to transform seq into exact copy of other seq being aligned distance.'''

'''Identity refers to number of characters that exactly match in both seq,
gaps are not counted, measurement is relational to shorter of 2 seq,
effect that sequence identity is not transitive,
i.e if sequence A=B and B=C then A is not neccessarily = C
(in terms of identity distance measure):
A: AAGGCTT, B: AAGGC, C: AAGGCAT'''

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2
from Levenshtein import distance
import matplotlib.pyplot as plt
import numpy as np

seqA = Seq("AAGGCTT")
seqB = Seq("AAGGC")
seqC = Seq("AAGGCAT")

AvB = pairwise2.align.localxx(seqA,seqB,one_alignment_only=True,score_only=True)
AvC = pairwise2.align.localxx(seqA,seqC,one_alignment_only=True,score_only=True)
BvC = pairwise2.align.localxx(seqB,seqC,one_alignment_only=True,score_only=True)

#Seq A and B are 100% identical as a result
#Hence identity not equivalent to being the same
print("AVB", AvB/len(seqB)*100) #100% identical to A
print(seqA == seqB) #False so identical is not the same as similarity
print("AVC", AvC/len(seqC)*100) #85% identical to A
print("BVC", BvC/len(seqC)*100) #71% identitcal to B

#Checking similarity using Hamming distance, Dotplots and
'''Hamming distance: between 2 strings of equal length is the number of
positions at which the corresponding symbols are different, measures
min mumber of substitutions for 1 strings to change into the other,
i.e. edit or Levenshtein distance - quantifying how 2 strings are dissimilar
used for error correction/detection, quantify similarity of DNA seq
Levenshtein distance works for unequal length strings unlike Hamming'''

#Hamming distance custom function, not good for unequal length sequences
def hamming_dist(lhs,rhs):
    return len([(x,y) for x,y in zip(lhs,rhs) if x != y]) #returns the number of disimilar bases in the seqs
print("Hamming distance:", hamming_dist(Seq("ACTAT"), Seq("ACTTA")))

#Levenshtein distance
print("Levenshtein distance:", distance(str(seqA),str(seqB)))

'''Dotplot to visualise comparision between 2 biological sequences
and identify regions of close similarity, dots placed in places of similarity
1 seq on x and other on y axis, when residues of both sequences match at the same location
of plot a dot is drawn on corresponding position

Can be used to visually inspect sequences for repeated/inverted sequences,
regions of low sequence complexity, similar regions, repeated sequences,
sequence rearrangments, RNA structure, gene order'''

def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n,m = len(seq1),len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t,seq1,seq2,nonblank=chr(0x25A0),blank=" "):
    print(" |" + seq2)
    print("-"*(2*len(seq2)))
    for label,row in zip(seq1,M):
        line = "".join(nonblank if s < t else blank for s in row)
        print(label + "|" + line)

def dotplot(seq1,seq2,k=1,t=1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M,t,seq1,seq2)

seq1 = Seq("ACTTAG")
seq2 = Seq("AC")
dotplot(seq1,seq2)

#Identical sequences show us a diagonal line like example below
dotplot("TREE", "TREE")

#Making dotplot better looking
def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    #x-axis has list of sequence 2
    xt = plt.xticks(np.arange(len(list(seq2))),list(seq2))
    #y-axis has list of sequence 1
    yt = plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()
#dotplotx(seq1,seq2)

dna1 = Seq("ATGATCTCGTAA", generic_dna)
dna2 = Seq("ATTATCGTCGTAA", generic_dna)
dotplot(dna1,dna2)
dotplotx(dna1,dna2)
#Straight diagonal line shown for any sequence that have similar sub sequences
dotplotx(dna1,dna1)
