"""
Neighbor-Joining algorithm implemented from scratch (Saitou & Nei, 1987).
"""

import numpy as np


def neighbor_joining(distance_matrix, labels):
    """
    Construct a phylogenetic tree using the Neighbor-Joining algorithm.

    Parameters
    ----------
    distance_matrix : np.ndarray
        Symmetric distance matrix (n x n).
    labels : list of str
        Taxon/sequence names corresponding to matrix rows/columns.

    Returns
    -------
    newick : str
        Tree in Newick format.
    """
    n = len(labels)
    if n < 2:
        raise ValueError("Need at least 2 taxa")
    if n == 2:
        d = distance_matrix[0, 1]
        return f"({labels[0]}:{d/2:.6f},{labels[1]}:{d/2:.6f});"

    # Work with mutable copies
    dm = distance_matrix.astype(float).copy()
    nodes = list(labels)

    while len(nodes) > 2:
        k = len(nodes)

        # Step 1: Compute net divergence r_i
        r = np.sum(dm, axis=1)

        # Step 2: Compute Q-matrix
        Q = np.full((k, k), float("inf"))
        for i in range(k):
            for j in range(i + 1, k):
                Q[i, j] = (k - 2) * dm[i, j] - r[i] - r[j]
                Q[j, i] = Q[i, j]

        # Step 3: Find the pair with minimum Q
        # Get indices of minimum Q value
        min_val = float("inf")
        mi, mj = 0, 1
        for i in range(k):
            for j in range(i + 1, k):
                if Q[i, j] < min_val:
                    min_val = Q[i, j]
                    mi, mj = i, j

        # Step 4: Compute branch lengths to new node u
        d_ij = dm[mi, mj]
        if k > 2:
            branch_i = d_ij / 2 + (r[mi] - r[mj]) / (2 * (k - 2))
            branch_j = d_ij - branch_i
        else:
            branch_i = d_ij / 2
            branch_j = d_ij / 2

        # Ensure non-negative branch lengths
        branch_i = max(branch_i, 0.0)
        branch_j = max(branch_j, 0.0)

        # Create new node
        new_node = f"({nodes[mi]}:{branch_i:.6f},{nodes[mj]}:{branch_j:.6f})"

        # Step 5: Compute distances from new node u to all remaining taxa
        keep = [idx for idx in range(k) if idx != mi and idx != mj]
        new_dists = []
        for idx in keep:
            d_uk = (dm[mi, idx] + dm[mj, idx] - d_ij) / 2
            d_uk = max(d_uk, 0.0)
            new_dists.append(d_uk)

        # Step 6: Build new distance matrix
        new_k = len(keep) + 1
        new_dm = np.zeros((new_k, new_k))

        # Fill distances among kept nodes
        for a_new, a_old in enumerate(keep):
            for b_new, b_old in enumerate(keep):
                new_dm[a_new, b_new] = dm[a_old, b_old]

        # Fill distances to new node (last row/column)
        u_idx = new_k - 1
        for a_new, d_uk in enumerate(new_dists):
            new_dm[a_new, u_idx] = d_uk
            new_dm[u_idx, a_new] = d_uk

        # Update working copies
        new_nodes = [nodes[idx] for idx in keep] + [new_node]
        dm = new_dm
        nodes = new_nodes

    # Final two nodes
    d = dm[0, 1]
    newick = f"({nodes[0]}:{d/2:.6f},{nodes[1]}:{d/2:.6f});"
    return newick
