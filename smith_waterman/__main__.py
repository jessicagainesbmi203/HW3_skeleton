from .utils import Seq
from .io import read_sequence, read_all_sequences, read_sub_matrix
from .algs import smithwaterman
import sys

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 3:
    print("Usage: python -m smith_waterman <seq directory> <scoring matrix filename>")
    sys.exit(0)
    
seqs = read_all_sequences(sys.argv[1])
#print(seqs)

sub_matrix = read_sub_matrix(sys.argv[2])
#print(sub_matrix)

smithwaterman(seqs[0],seqs[1], sys.argv[2])