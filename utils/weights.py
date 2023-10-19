import numpy as np
# import numba
# from numba import prange

# @numba.jit(nopython=True, fastmath=True, parallel=True)
# def calc_num_cluster_members_nogaps_parallel(matrix, identity_threshold, invalid_value):
#     """
#     From EVCouplings: https://github.com/debbiemarkslab/EVcouplings
#     Calculate number of sequences in alignment
#     within given identity_threshold of each other
#     Parameters
#     ----------
#     matrix : np.array
#         N x L matrix containing N sequences of length L.
#         Matrix must be mapped to range(0, num_symbols) using
#         map_matrix function
#     identity_threshold : float
#         Sequences with at least this pairwise identity will be
#         grouped in the same cluster.
#     invalid_value : int
#         Value in matrix that is considered invalid, e.g. gap or lowercase character.
#     Returns
#     -------
#     np.array
#         Vector of length N containing number of cluster
#         members for each sequence (inverse of sequence
#         weight)
#     """
#     N, L = matrix.shape
#     L = 1.0 * L
#
#     # Empty sequences are filtered out before this function and are ignored
#     # minimal cluster size is 1 (self)
#     num_neighbors = np.ones((N))
#     L_non_gaps = L - np.sum(matrix == invalid_value, axis=1)  # Edit: From EVE, use the non-gapped length
#     # compare all pairs of sequences
#     # Edit: Rewrote loop without any dependencies between inner and outer loops, so that it can be parallelized
#     for i in prange(N):
#         num_neighbors_i = 1  # num_neighbors_i = 0  # TODO why did I make this 0 again? Probably because I thought I'd have to count i == j
#         for j in range(N):
#             if i == j:
#                 continue
#             pair_matches = 0
#             for k in range(L):  # This should hopefully be vectorised by numba
#                 if matrix[i, k] == matrix[j, k] and matrix[
#                     i, k] != invalid_value:  # Edit(Lood): Don't count gaps as matches
#                     pair_matches += 1
#             # Edit(Lood): Calculate identity as fraction of non-gapped positions (so this similarity is asymmetric)
#             # Note: Changed >= to > to match EVE / DeepSequence code
#             if pair_matches / L_non_gaps[i] > identity_threshold:
#                 num_neighbors_i += 1
#
#         num_neighbors[i] = num_neighbors_i
#
#     return num_neighbors

def calc_weights_evcouplings(matrix_mapped, identity_threshold, empty_value, num_cpus=1):
    """
        From EVCouplings: https://github.com/debbiemarkslab/EVcouplings
        Calculate weights for sequences in alignment by
        clustering all sequences with sequence identity
        greater or equal to the given threshold.
        Parameters
        ----------
        identity_threshold : float
            Sequence identity threshold
        """
    empty_idx = is_empty_sequence_matrix(matrix_mapped,
                                         empty_value=empty_value)  # e.g. sequences with just gaps or lowercase, no valid AAs
    N = matrix_mapped.shape[0]

    # Original EVCouplings code structure, plus gap handling
    if num_cpus != 1:
        num_cluster_members = calc_num_cluster_members_nogaps_parallel(matrix_mapped[~empty_idx], identity_threshold,
                                                                       invalid_value=empty_value)
    # Empty sequences: weight 0
    weights = np.zeros((N))
    weights[~empty_idx] = 1.0 / num_cluster_members
    return weights

def is_empty_sequence_matrix(matrix, empty_value):
    assert len(matrix.shape) == 2, f"Matrix must be 2D; shape={matrix.shape}"
    assert isinstance(empty_value, (int, float)), f"empty_value must be a number; type={type(empty_value)}"
    # Check for each sequence if all positions are equal to empty_value
    empty_idx = np.all((matrix == empty_value), axis=1)
    return empty_idx
