"""
Shared utilities: substitution matrix loading, FASTA I/O.
"""

import os
import numpy as np


def load_blosum62():
    """
    Load BLOSUM62 substitution matrix from Biopython.
    Returns a dict-of-dicts: blosum[aa1][aa2] -> score.
    """
    from Bio.Align import substitution_matrices
    mat = substitution_matrices.load("BLOSUM62")

    blosum = {}
    for aa1 in mat.alphabet:
        blosum[aa1] = {}
        for aa2 in mat.alphabet:
            blosum[aa1][aa2] = mat[aa1, aa2]
    return blosum


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def write_fasta(sequences, filepath):
    """
    Write sequences to a FASTA file.
    sequences: dict of name -> sequence string (may contain gaps).
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def read_fasta(filepath):
    """
    Read a FASTA file. Returns dict of name -> sequence string.
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        sequences[current_name] = "".join(current_seq)

    return sequences


def percent_identity(seq1, seq2):
    """
    Compute percent identity between two aligned sequences (same length, with gaps).
    """
    assert len(seq1) == len(seq2), "Aligned sequences must have equal length"
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-")
    aligned_positions = sum(1 for a, b in zip(seq1, seq2) if a != "-" or b != "-")
    if aligned_positions == 0:
        return 0.0
    return matches / aligned_positions


def jukes_cantor_distance(p):
    """
    Jukes-Cantor correction for protein sequences.
    p = fraction of differing positions (1 - fraction_identical).
    Returns evolutionary distance, or 3.0 as cap for saturated divergence.
    """
    if p >= 19 / 20:
        return 3.0  # saturated
    import math
    d = -19 / 20 * math.log(1 - 20 / 19 * p)
    return d
