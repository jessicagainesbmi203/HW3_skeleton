import os
import glob

def read_sequence(filepath):

    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".fa":
        raise IOError("%s is not a FA file"%filepath)
    sequence = ''
    # open fa file
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in f:
            if line[0:1] != '>':
                sequence = sequence + line.strip('\n')
    return sequence