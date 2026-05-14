"""
Conservation scoring from multiple sequence alignments.
"""

import math
import numpy as np
from collections import Counter


def shannon_entropy(column):
    """
    Compute Shannon entropy (in bits) for a single MSA column.

    Parameters
    ----------
    column : str or list
        Amino acid characters at one alignment position (gaps excluded).

    Returns
    -------
    float
        Shannon entropy H. Lower = more conserved.
    """
    # Exclude gaps
    residues = [c for c in column if c != "-"]
    if len(residues) == 0:
        return 0.0

    counts = Counter(residues)
    total = len(residues)
    H = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            H -= p * math.log2(p)
    return H


def conservation_scores(alignment, reference_name=None):
    """
    Compute per-residue conservation scores from an MSA.

    Score = 1 - H / H_max, where H_max = log2(20).
    Score of 1.0 = perfectly conserved; 0.0 = maximum entropy.

    Parameters
    ----------
    alignment : dict
        name -> aligned sequence (all same length).
    reference_name : str, optional
        If provided, map scores to ungapped positions of this sequence.
        Columns where the reference has a gap are excluded.

    Returns
    -------
    scores : np.ndarray
        Conservation scores. If reference_name is given, length = length
        of reference sequence (ungapped). Otherwise, length = alignment width.
    positions : list of int or None
        If reference_name given, list of 1-based residue positions.
        Otherwise None.
    """
    seqs = list(alignment.values())
    aln_length = len(seqs[0])
    H_max = math.log2(20)

    # Compute entropy at every column
    col_entropies = []
    for col_idx in range(aln_length):
        column = [s[col_idx] for s in seqs]
        H = shannon_entropy(column)
        col_entropies.append(H)

    if reference_name is not None:
        ref_seq = alignment[reference_name]
        scores = []
        positions = []
        res_num = 0
        for col_idx in range(aln_length):
            if ref_seq[col_idx] != "-":
                res_num += 1
                score = 1.0 - col_entropies[col_idx] / H_max
                scores.append(score)
                positions.append(res_num)
        return np.array(scores), positions
    else:
        scores = np.array([1.0 - H / H_max for H in col_entropies])
        return scores, None
