#Sequence alignment is a method of arranging DNA, RNA, aa or protein sequences
#to identify regions of similarity

#Similarities identified include functional, structural or evolutionary (homology)
#Homology: common ancestor

#Terms: Matches, mismatches or gap, e.g below
#seq 1 = A C T C G
#seq 2 = A T T C _
#Match: A and A, T and T, C and C, mismatch: C and T, Gap: G and _

#Types of alignment: Global and Local
#Global - finds best condordnace/agreement between all characters in both seq
#         Mostly from end to end (5' to 3'), by needle
#Local - finds the subsequences that align best by water

'''Local is used: 2 sequences have a small matched region and different lengths
                  overlapping sequences, 1 is subsequence over the other
                  Blast and EMBOSS are tools for alignment'''

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq

seq1 = Seq("ACTCG")
seq2 = Seq("ATTCG")

#Global alignment
alignments = pairwise2.align.globalxx(seq1,seq2)
print(alignments)

#To visualise alignments and the score, only the 1st one is shown
#. for mismatch, - for gap extension and | for match
print(format_alignment(*alignments[0]))

#View all global alignments
print("Global alignments:")
for a in alignments:
    print(format_alignment(*a))

#View all local alignments
local_alignments = pairwise2.align.localxx(seq1,seq2)
print("Local alignments:")
for l in local_alignments:
    print(format_alignment(*l))

#Only scores printed
print("Only score of 1st global alignment printed")
alignments2 = pairwise2.align.globalxx(seq1,seq2,one_alignment_only=True, score_only=True)
print(alignments2)

#Similarity (percentage) between sequences using global alignment
#if > 50% then characterised as similar
print("Using global alignment: {0}".format(alignments2/len(seq1)*100))

#Similarity (percentage) between sequences local alignment
#if > 50% then characterised as similar
local_alignments2 = pairwise2.align.localxx(seq1,seq2,one_alignment_only=True, score_only=True)
print("Using local alignment: {0}".format(local_alignments2/len(seq1)*100))

#Find all possible global alignments with max similarity
#Match: 2 points, mismatch: -1 point, gap: -0.5 point, for each gap extension: -0.1 point
#Global alignment with max size
glb_alignment = pairwise2.align.globalms(seq1,seq2,2,-1,-0.5,-0.1)
print("Global alignment with score based on points:")
for a in glb_alignment:
    print(format_alignment(*a))

