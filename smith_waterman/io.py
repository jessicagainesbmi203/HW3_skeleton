import os
import glob
from .utils import Seq
import pandas as pd
import re

def read_sequence(filepath):
    """
    Read a single sequence from a .fa file
    Input: Filepath of the .fa file describing the sequence
    Output: a Seq object with the name and sequence from the file
    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".fa":
        raise IOError("%s is not a FA file"%filepath)
    seq = Seq(name[0])
    build_seq = ''
    # open fa file
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in f:
            if line[0:1] != '>':
                build_seq = build_seq + line.strip('\n')
    build_seq = build_seq.replace('x','')
    seq.sequence = build_seq
    return seq
    
def read_all_sequences(dir):
    """
    Read all sequences from a directory into a list of Seq objects
    Input: the directory containing .fa files describing the sequences
    Output: a list of Seq objects
    """
    files = glob.glob(dir + '/*.fa')

    sequences = []
    # iterate over each .fa file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):

        sequences.append(read_sequence(filepath))

    print("Read in %d sequences"%len(sequences))

    return sequences
    
def read_sub_matrix(filepath):
    """
    Read in a substitution matrix
    Inputs: filepath: filepath to document with the substitution matrix
    Outputs: a dataframe containing the substitution matrix
    """
    i = 0
    df = pd.DataFrame()
    index = []
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in f:
            if line[0:1] != "#":
                if i == 0:
                    df = pd.DataFrame(columns=line.split())
                    index = line.split()
                else:
                    df.loc[i,:] = list(map(int,line.split()))
                i += 1
        df.index = index
    return df

def read_pairs(filepath):
    """
    Read in a file containing pairs of sequences
    Inputs: filepath: the file path to document containing pairs of sequences
    Output: a list of pairs of sequences
    """
    pairs = list()
    with open(filepath, "r") as f:
        for line in f:
            seqs = line.split()
            seq1 = read_sequence(seqs[0])
            seq2 = read_sequence(seqs[1])
            pairs.append([seq1,seq2])
    return pairs




