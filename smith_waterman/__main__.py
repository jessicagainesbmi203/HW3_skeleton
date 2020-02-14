from .utils import Seq
from .io import read_sequence, read_all_sequences
from .algs import smithwaterman
import sys

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 2:
    print("Usage: python -m smith_waterman <seq directory>")
    sys.exit(0)
    
seqs = read_all_sequences(sys.argv[1])
print(seqs)