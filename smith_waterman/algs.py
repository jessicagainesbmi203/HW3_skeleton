import numpy as np
import pandas as pd
from .io import read_sub_matrix

def smithwaterman(seq1,seq2,filepath):
    """
    Perform sequence alignment according to the Smith-Waterman Algorithm
    Inputs: Two sequences to compare, filepath to desired substitution matrix
    Outputs:
    """
    scoring_matrix = pd.DataFrame(columns=list(seq1.sequence),index=list(seq2.sequence))
    print(scoring_matrix)
    sub_matrix = read_sub_matrix(filepath)
    
    return None


def score():
    return None

def roc():
    """
    Create a Receiver Operator Curve (ROC)
    X-axis : fraction of false positives
    Y-axis : fraction of true positives
    """
    return None
