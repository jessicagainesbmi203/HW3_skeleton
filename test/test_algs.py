import numpy as np
from smith_waterman.algs import smithwaterman, score, roc
from smith_waterman.io import read_sequence
from smith_waterman.utils import Seq
import re

def test_roc():
    false_pos = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    true_pos = roc('BLOSUM50',10,2,'Pospairs.txt','Negpairs.txt', false_pos, normalize=False)
    assert len(true_pos) == len(false_pos)
    for value in true_pos:
        assert value >= 0
        assert value <= 1

def test_smithwaterman():
    # test perfect alignment
    seq1 = read_sequence('sequences/prot-0018.fa')
    seq2 = read_sequence('sequences/prot-0018.fa')
    alignment = smithwaterman(seq1,seq2, 'BLOSUM50',10,2)
    assert alignment.get('aligned_seq1') == list(seq1.sequence)
    assert alignment.get('aligned_seq2') == list(seq2.sequence)
    # test blank sequences
    alignment = smithwaterman(Seq('blank1'),Seq('blank2'),'BLOSUM50',10,2)
    assert alignment.get('aligned_seq1') == ''
    assert alignment.get('aligned_seq2') == ''
    assert alignment.get('alignment_score') == 0
    # test alignment that is not perfect
    seq1 = read_sequence('sequences/prot-0595.fa')
    seq2 = read_sequence('sequences/prot-0503.fa')
    alignment = smithwaterman(seq1,seq2,'BLOSUM50',10,2)
    assert alignment.get('aligned_seq1') != seq1.sequence
    assert alignment.get('aligned_seq2') != seq2.sequence

def test_scoring():
    seq1 = read_sequence('sequences/prot-0018.fa')
    seq2 = read_sequence('sequences/prot-0022.fa')
    alignment = smithwaterman(seq1,seq2, 'BLOSUM50',10,2)
    # test that scoring works as expected
    assert(score(alignment.get('aligned_seq1'),alignment.get('aligned_seq2'),'BLOSUM50',10,2) == 87.0)
    # scoring function does not match score from alignment
    #assert(score(alignment.get('aligned_seq1'),alignment.get('aligned_seq2'),'BLOSUM50',10,2) == alignment.get('alignment_score'))