#!/usr/bin/env python3
from Bio import AlignIO
import numpy as np
import sys

VALID_NT = {"A", "C", "G", "T"}

def calculate_entropy(alignment):
    entropy = []
    for col in range(alignment.get_alignment_length()):
        counts = {}
        for nt in alignment[:, col]:
            base = nt.upper()
            if base in VALID_NT:
                counts[base] = counts.get(base, 0) + 1
        total = sum(counts.values())
        if total == 0:
            entropy.append(0.0)
            continue
        s = -sum((c/total) * np.log2(c/total) for c in counts.values())
        entropy.append(s)
    return np.mean(entropy)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <alignment.fasta>")
        sys.exit(1)
    aln = AlignIO.read(sys.argv[1], "fasta")
    mean_entropy = calculate_entropy(aln)
    print(f"{sys.argv[1]}\t{mean_entropy:.4f}")