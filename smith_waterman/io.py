import os
import glob
from .utils import Seq
import pandas as pd

def read_sequence(filepath):

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
    seq.sequence = build_seq
    return seq
    
def read_all_sequences(dir):
    files = glob.glob(dir + '/*.fa')

    sequences = []
    # iterate over each .fa file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):

        sequences.append(read_sequence(filepath))

    print("Read in %d sequences"%len(sequences))

    return sequences
    
def read_scoring_matrix(filepath):
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
                    df.loc[i,:] = line.split()
                i += 1
        df.index = index
    return df
