from .utils import Seq
from .io import read_sequence, read_all_sequences, read_sub_matrix, read_pairs
from .algs import smithwaterman, score, get_false_pos_rate, roc
import sys
import pandas as pd
import numpy as np

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 3:
    print("Usage: python -m smith_waterman <seq directory> <scoring matrix filename>")
    sys.exit(0)
# read in seqences, substitution matrix
seqs = read_all_sequences(sys.argv[1])
sub_matrix = read_sub_matrix(sys.argv[2])
# run smith waterman on two sample sequences
results = smithwaterman(seqs[3],seqs[4], sys.argv[2],10,2)
print(results.get('aligned_seq1'))
print(results.get('aligned_seq2'))
print(results.get('alignment_score'))
# find the alignment score of the two sequences
print(score(results.get('aligned_seq1'),results.get('aligned_seq2'),sys.argv[2],10,2))
print('False Positive Rate = ')
print(get_false_pos_rate(sys.argv[2],10,2,'Pospairs.txt','Negpairs.txt',0.7))

"""
# find optimal gap opening and gap extension penalties
false_positive_rates = list()
for open in (1,5,10,15,20):
    for extend in (1,3,5):
        rate = get_false_pos_rate(sub_matrix,open,extend,'Pospairs.txt','Negpairs.txt',0.7)
        false_positive_rates.append((open,extend,rate))
print(false_positive_rates)

# find the best substituion matrix
false_positive_rates = list()
for file in ():
    sub = read_sub_matrix(file)
    rate = get_false_pos_rate(sub,10,2,'Pospairs.txt','Negpairs.txt')
    false_positive_rates.append((sub,rate))
print(false_positive_rates)
"""        
print(roc(sys.argv[2],10,2,'Pospairs.txt','Negpairs.txt',[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], normalize=False))
print(roc(sys.argv[2],10,2,'Pospairs.txt','Negpairs.txt',[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], normalize=True))
        
        