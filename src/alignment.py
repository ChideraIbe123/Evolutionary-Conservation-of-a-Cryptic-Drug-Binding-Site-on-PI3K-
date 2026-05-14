"""
Sequence alignment algorithms implemented from scratch.
- Needleman-Wunsch with affine gap penalties
- Profile-profile alignment
- Progressive multiple sequence alignment
"""

import numpy as np
from collections import Counter


def needleman_wunsch(seq1, seq2, sub_matrix, gap_open=-10, gap_extend=-0.5):
    """
    Global pairwise alignment with affine gap penalties.

    Uses three DP matrices:
      M[i,j]  = best score ending with seq1[i] aligned to seq2[j]
      Ix[i,j] = best score ending with a gap in seq2 (insertion in seq1)
      Iy[i,j] = best score ending with a gap in seq1 (insertion in seq2)

    Parameters
    ----------
    seq1, seq2 : str
        Protein sequences (ungapped).
    sub_matrix : dict of dict
        Substitution scores, e.g. BLOSUM62[aa1][aa2].
    gap_open : float
        Penalty for opening a new gap (negative).
    gap_extend : float
        Penalty for extending an existing gap (negative).

    Returns
    -------
    aligned1, aligned2 : str
        Aligned sequences with '-' for gaps.
    score : float
        Optimal alignment score.
    """
    n, m = len(seq1), len(seq2)
    NEG_INF = float("-inf")

    # Initialize DP matrices
    M = np.full((n + 1, m + 1), NEG_INF)
    Ix = np.full((n + 1, m + 1), NEG_INF)
    Iy = np.full((n + 1, m + 1), NEG_INF)

    # Traceback matrices: 0=M, 1=Ix, 2=Iy
    trM = np.zeros((n + 1, m + 1), dtype=int)
    trIx = np.zeros((n + 1, m + 1), dtype=int)
    trIy = np.zeros((n + 1, m + 1), dtype=int)

    # Base cases
    M[0, 0] = 0
    for i in range(1, n + 1):
        Ix[i, 0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        Iy[0, j] = gap_open + (j - 1) * gap_extend

    # Fill
    for i in range(1, n + 1):
        aa1 = seq1[i - 1]
        for j in range(1, m + 1):
            aa2 = seq2[j - 1]

            # Match/mismatch score
            s = sub_matrix.get(aa1, {}).get(aa2, -1)

            # M[i,j]: align aa1 with aa2
            candidates_m = [
                (M[i-1, j-1] + s, 0),
                (Ix[i-1, j-1] + s, 1),
                (Iy[i-1, j-1] + s, 2),
            ]
            best = max(candidates_m, key=lambda x: x[0])
            M[i, j] = best[0]
            trM[i, j] = best[1]

            # Ix[i,j]: gap in seq2 (consume seq1[i])
            candidates_ix = [
                (M[i-1, j] + gap_open, 0),
                (Ix[i-1, j] + gap_extend, 1),
            ]
            best = max(candidates_ix, key=lambda x: x[0])
            Ix[i, j] = best[0]
            trIx[i, j] = best[1]

            # Iy[i,j]: gap in seq1 (consume seq2[j])
            candidates_iy = [
                (M[i, j-1] + gap_open, 0),
                (Iy[i, j-1] + gap_extend, 2),
            ]
            best = max(candidates_iy, key=lambda x: x[0])
            Iy[i, j] = best[0]
            trIy[i, j] = best[1]

    # Find best ending score
    end_scores = [(M[n, m], 0), (Ix[n, m], 1), (Iy[n, m], 2)]
    score, current_mat = max(end_scores, key=lambda x: x[0])

    # Traceback
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        # Boundary guards: if one index is 0, force the appropriate gap state
        if i == 0:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
            current_mat = 2  # Iy
            continue
        if j == 0:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
            current_mat = 1  # Ix
            continue

        if current_mat == 0:  # M
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            current_mat = trM[i, j]
            i -= 1
            j -= 1
        elif current_mat == 1:  # Ix
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            current_mat = trIx[i, j]
            i -= 1
        else:  # Iy
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            current_mat = trIy[i, j]
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), score


def compute_distance_matrix(sequences, sub_matrix, gap_open=-10, gap_extend=-0.5):
    """
    Compute all-vs-all pairwise distance matrix from NW alignments.

    Parameters
    ----------
    sequences : dict
        name -> sequence string (ungapped).

    Returns
    -------
    dist_matrix : np.ndarray
        Symmetric distance matrix.
    labels : list
        Sequence names in matrix order.
    """
    from src.utils import jukes_cantor_distance

    labels = list(sequences.keys())
    n = len(labels)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            aln1, aln2, _ = needleman_wunsch(
                sequences[labels[i]], sequences[labels[j]],
                sub_matrix, gap_open, gap_extend
            )
            # Compute fraction of differing aligned positions
            aligned_pos = sum(1 for a, b in zip(aln1, aln2) if a != "-" or b != "-")
            identical = sum(1 for a, b in zip(aln1, aln2) if a == b and a != "-")
            if aligned_pos > 0:
                p = 1.0 - identical / aligned_pos
            else:
                p = 1.0
            d = jukes_cantor_distance(p)
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    return dist_matrix, labels


def _profile_from_sequences(aligned_seqs):
    """
    Convert a list of aligned sequences (same length) into a profile.
    Profile is a list of Counter dicts, one per column.
    """
    if not aligned_seqs:
        return []
    length = len(aligned_seqs[0])
    profile = []
    for col_idx in range(length):
        counts = Counter()
        for seq in aligned_seqs:
            aa = seq[col_idx]
            if aa != "-":
                counts[aa] += 1
        profile.append(counts)
    return profile


def _profile_score(prof1_col, prof2_col, sub_matrix):
    """Average substitution score between two profile columns."""
    total_score = 0.0
    total_pairs = 0
    for aa1, c1 in prof1_col.items():
        for aa2, c2 in prof2_col.items():
            s = sub_matrix.get(aa1, {}).get(aa2, -1)
            total_score += s * c1 * c2
            total_pairs += c1 * c2
    if total_pairs == 0:
        return 0.0
    return total_score / total_pairs


def profile_align(seqs1, seqs2, sub_matrix, gap_open=-10, gap_extend=-0.5):
    """
    Align two groups of aligned sequences (profiles) using NW with affine gaps.

    Parameters
    ----------
    seqs1 : list of str
        First group of aligned sequences (same length within group).
    seqs2 : list of str
        Second group of aligned sequences.

    Returns
    -------
    merged : list of str
        All sequences aligned together (seqs1 first, then seqs2).
    """
    prof1 = _profile_from_sequences(seqs1)
    prof2 = _profile_from_sequences(seqs2)
    len1, len2 = len(prof1), len(prof2)

    NEG_INF = float("-inf")
    M = np.full((len1 + 1, len2 + 1), NEG_INF)
    Ix = np.full((len1 + 1, len2 + 1), NEG_INF)
    Iy = np.full((len1 + 1, len2 + 1), NEG_INF)

    trM = np.zeros((len1 + 1, len2 + 1), dtype=int)
    trIx = np.zeros((len1 + 1, len2 + 1), dtype=int)
    trIy = np.zeros((len1 + 1, len2 + 1), dtype=int)

    M[0, 0] = 0
    for i in range(1, len1 + 1):
        Ix[i, 0] = gap_open + (i - 1) * gap_extend
    for j in range(1, len2 + 1):
        Iy[0, j] = gap_open + (j - 1) * gap_extend

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            s = _profile_score(prof1[i-1], prof2[j-1], sub_matrix)

            candidates_m = [
                (M[i-1, j-1] + s, 0),
                (Ix[i-1, j-1] + s, 1),
                (Iy[i-1, j-1] + s, 2),
            ]
            best = max(candidates_m, key=lambda x: x[0])
            M[i, j] = best[0]
            trM[i, j] = best[1]

            candidates_ix = [
                (M[i-1, j] + gap_open, 0),
                (Ix[i-1, j] + gap_extend, 1),
            ]
            best = max(candidates_ix, key=lambda x: x[0])
            Ix[i, j] = best[0]
            trIx[i, j] = best[1]

            candidates_iy = [
                (M[i, j-1] + gap_open, 0),
                (Iy[i, j-1] + gap_extend, 2),
            ]
            best = max(candidates_iy, key=lambda x: x[0])
            Iy[i, j] = best[0]
            trIy[i, j] = best[1]

    end_scores = [(M[len1, len2], 0), (Ix[len1, len2], 1), (Iy[len1, len2], 2)]
    _, current_mat = max(end_scores, key=lambda x: x[0])

    # Traceback — build column operations
    ops = []  # list of ("match", i, j) or ("gap1", j) or ("gap2", i)
    i, j = len1, len2
    while i > 0 or j > 0:
        # Boundary guards
        if i == 0:
            ops.append(("gap1", j - 1))
            j -= 1
            current_mat = 2
            continue
        if j == 0:
            ops.append(("gap2", i - 1))
            i -= 1
            current_mat = 1
            continue

        if current_mat == 0:
            ops.append(("match", i - 1, j - 1))
            current_mat = trM[i, j]
            i -= 1
            j -= 1
        elif current_mat == 1:
            ops.append(("gap2", i - 1))
            current_mat = trIx[i, j]
            i -= 1
        else:
            ops.append(("gap1", j - 1))
            current_mat = trIy[i, j]
            j -= 1

    ops.reverse()

    # Build merged alignment
    result1 = [[] for _ in seqs1]
    result2 = [[] for _ in seqs2]

    for op in ops:
        if op[0] == "match":
            ci, cj = op[1], op[2]
            for k, seq in enumerate(seqs1):
                result1[k].append(seq[ci])
            for k, seq in enumerate(seqs2):
                result2[k].append(seq[cj])
        elif op[0] == "gap2":
            ci = op[1]
            for k, seq in enumerate(seqs1):
                result1[k].append(seq[ci])
            for k in range(len(seqs2)):
                result2[k].append("-")
        else:  # gap1
            cj = op[1]
            for k in range(len(seqs1)):
                result1[k].append("-")
            for k, seq in enumerate(seqs2):
                result2[k].append(seq[cj])

    merged = ["".join(r) for r in result1] + ["".join(r) for r in result2]
    return merged


def progressive_msa(sequences, guide_tree_newick, sub_matrix,
                    gap_open=-10, gap_extend=-0.5):
    """
    Progressive multiple sequence alignment using a guide tree.

    Parameters
    ----------
    sequences : dict
        name -> sequence string (ungapped).
    guide_tree_newick : str
        Newick-format guide tree with leaf names matching sequence keys.
    sub_matrix : dict of dict
        Substitution matrix.

    Returns
    -------
    alignment : dict
        name -> aligned sequence string (with gaps).
    """
    from Bio import Phylo
    from io import StringIO

    tree = Phylo.read(StringIO(guide_tree_newick), "newick")

    # Map each clade to a list of (name, aligned_seq)
    alignment_cache = {}

    def _get_alignment(clade):
        clade_id = id(clade)
        if clade_id in alignment_cache:
            return alignment_cache[clade_id]

        if clade.is_terminal():
            name = clade.name
            result = [(name, sequences[name])]
            alignment_cache[clade_id] = result
            return result

        children = clade.clades
        left_aln = _get_alignment(children[0])
        right_aln = _get_alignment(children[1])

        left_seqs = [s for _, s in left_aln]
        right_seqs = [s for _, s in right_aln]

        merged = profile_align(left_seqs, right_seqs, sub_matrix,
                               gap_open, gap_extend)

        left_names = [n for n, _ in left_aln]
        right_names = [n for n, _ in right_aln]
        all_names = left_names + right_names

        result = list(zip(all_names, merged))
        alignment_cache[clade_id] = result
        return result

    final = _get_alignment(tree.root)
    return {name: seq for name, seq in final}


def upgma_tree(dist_matrix, labels):
    """
    Simple UPGMA tree construction for use as a fallback guide tree.
    Returns a Newick string.
    """
    n = len(labels)
    # Work with copies
    dm = dist_matrix.copy()
    nodes = list(labels)
    sizes = [1] * n

    while len(nodes) > 1:
        # Find minimum distance pair
        min_d = float("inf")
        mi, mj = 0, 1
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if dm[i, j] < min_d:
                    min_d = dm[i, j]
                    mi, mj = i, j

        # Branch length
        bl = min_d / 2

        # New node label
        new_node = f"({nodes[mi]}:{bl:.6f},{nodes[mj]}:{bl:.6f})"
        new_size = sizes[mi] + sizes[mj]

        # Compute distances from new node to all others
        new_dm = np.zeros((len(nodes) - 1, len(nodes) - 1))
        new_nodes = []
        new_sizes = []

        keep = [k for k in range(len(nodes)) if k != mi and k != mj]
        for idx, k in enumerate(keep):
            new_nodes.append(nodes[k])
            new_sizes.append(sizes[k])
        new_nodes.append(new_node)
        new_sizes.append(new_size)

        # Fill distances for kept nodes
        for a_idx, a in enumerate(keep):
            for b_idx, b in enumerate(keep):
                new_dm[a_idx, b_idx] = dm[a, b]

        # Distance from new node to each kept node
        new_idx = len(keep)
        for a_idx, a in enumerate(keep):
            d = (sizes[mi] * dm[mi, a] + sizes[mj] * dm[mj, a]) / new_size
            new_dm[a_idx, new_idx] = d
            new_dm[new_idx, a_idx] = d

        dm = new_dm
        nodes = new_nodes
        sizes = new_sizes

    return nodes[0] + ";"
