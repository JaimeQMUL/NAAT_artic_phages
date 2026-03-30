from Bio import AlignIO
from collections import Counter
import math
import matplotlib.pyplot as plt
import numpy as np

alignment_file = "../../uvsx/data/references/aligned.fasta"  # output from MAFFT
alignment = AlignIO.read(alignment_file, "fasta")

print(f"Number of sequences: {len(alignment)}")
print(f"Alignment length: {alignment.get_alignment_length()}")