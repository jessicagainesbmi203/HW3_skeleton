import numpy as np
from smith_waterman.algs import smithwaterman, score
from smith_waterman.io import read_sequence
from smith_waterman.utils import Seq
import re

def test_roc():
    return None

def test_smithwaterman():
    seq1 = read_sequence('sequences/prot-0018.fa')
    seq2 = read_sequence('sequences/prot-0018.fa')
    alignment = smithwaterman(seq1,seq2, 'BLOSUM50',10,2)
    assert alignment.get('aligned_seq1') == list(seq1.sequence)
    assert alignment.get('aligned_seq2') == list(seq2.sequence)
    alignment = smithwaterman(Seq('blank1'),Seq('blank2'),'BLOSUM50',10,2)
    assert alignment.get('aligned_seq1') == ''
    assert alignment.get('aligned_seq2') == ''
    assert alignment.get('alignment_score') == 0

def test_scoring():
    seq1 = read_sequence('sequences/prot-0018.fa')
    seq2 = read_sequence('sequences/prot-0022.fa')
    alignment = smithwaterman(seq1,seq2, 'BLOSUM50',10,2)
    assert(score(alignment.get('aligned_seq1'),alignment.get('aligned_seq2'),'BLOSUM50',10,2) == -33.0)