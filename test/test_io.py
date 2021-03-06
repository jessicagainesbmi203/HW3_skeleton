from smith_waterman.io import read_sequence, read_all_sequences, read_sub_matrix
from smith_waterman.utils import Seq
import numpy as np

def test_read_sequence():
    seq = read_sequence('sequences/prot-0004.fa')
    correct_seq = 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
    assert seq.__repr__() == 'prot-0004'
    assert seq.sequence == correct_seq
    
def test_read_all_sequences():
    correct_seq = 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
    seqs = read_all_sequences('sequences')
    assert len(seqs) == 182
    
def test_read_sub_matrix():
    sub_matrix = read_sub_matrix('BLOSUM50')
    assert sub_matrix.shape == (24,24)
    assert not sub_matrix.isnull().values.all()