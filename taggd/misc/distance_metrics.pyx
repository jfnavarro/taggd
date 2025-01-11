# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
"""
Some functions to compute distance
between sequences
"""

cimport numpy as np
import numpy as np

cdef int hamming_distance(str seq1, str seq2, int limit=0):
    """
    Calculates the Hamming distance between two equal-length sequences.

    Args:
        seq1: The first sequence.
        seq2: The second sequence.
        limit: The maximum distance limit. If the distance exceeds this limit,
            the calculation aborts and returns `limit + 1`.

    Returns:
        The Hamming distance between the sequences, or `limit + 1` if the limit is exceeded.
    """
    cdef int len1 = len(seq1)
    cdef int len2 = len(seq2)
    if len1 != len2:
        raise ValueError("Sequences must have equal lengths.")
    cdef int i = 0
    cdef int sum = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            sum += 1
        if limit > 0 and sum > limit:
            return limit + 1
    return sum

cdef int levenshtein_distance(str seq1, str seq2, int limit=0):
    """
    Calculates the Levenshtein distance between two sequences.

    The Levenshtein distance allows for substitutions, insertions, and deletions.
    The sequences do not need to have equal lengths.

    Args:
        seq1: The first sequence.
        seq2: The second sequence.
        limit: The maximum distance limit. If the calculated distance exceeds
            this limit, the calculation aborts and returns `limit + 1`.

    Returns:
        The Levenshtein distance between the sequences, or `limit + 1` if the limit is exceeded.
    """
    cdef int len1 = len(seq1)
    cdef int len2 = len(seq2)
    if len1 == 0:
        return len2
    if len2 == 0:
        return len1

    # Memoryview arrays for the rows
    cdef int[:] one_ago = np.zeros(len2 + 1, dtype=np.int32)
    cdef int[:] this_row = np.arange(len2 + 1, dtype=np.int32)
    cdef int x, y
    cdef int del_cost, add_cost, sub_cost

    for x in range(len1):
        # Swap rows
        one_ago, this_row = this_row, one_ago

        # Initialize the first column of the current row
        this_row[0] = x + 1

        for y in range(len2):
            del_cost = one_ago[y + 1] + 1  # Deletion cost
            add_cost = this_row[y] + 1     # Insertion cost
            sub_cost = one_ago[y] + (seq1[x] != seq2[y])  # Substitution cost

            # Take the minimum of the three operations
            this_row[y + 1] = min(del_cost, add_cost, sub_cost)

        # Early stopping if the minimum distance exceeds the limit
        if limit > 0 and min(this_row) > limit:
            return limit + 1

    # Return the distance at the bottom-right corner
    return this_row[len2]

cdef int subglobal_distance(str s1, str s2):
    """
    Calculates the edit distance for a sub-global alignment of two sequences.

    Args:
        s1: The longer sequence (probe).
        s2: The shorter sequence (query).

    Returns:
        The minimum edit distance for aligning `s2` to `s1`.
    """
    cdef int xLen = len(s1)
    cdef int yLen = len(s2)
    if xLen < yLen:
        raise ValueError("Sub-global edit distance is undefined for sequences "
                         "where the probe is shorter than the aligned sequence.")

    cdef int x, y

    # Declare and initialize the DP array
    cdef np.ndarray[np.uint32_t, ndim=2] d = np.empty((xLen + 1, yLen + 1), dtype=np.uint32)

    # Initialize DP matrix
    for x in range(0, xLen + 1):
        d[x, 0] = 0  # No penalty for aligning with an empty query
    for y in range(1, yLen + 1):
        d[0, y] = y  # Full penalty for unaligned query

    # Fill the DP matrix
    for x in range(1, xLen + 1):
        for y in range(1, yLen + 1):
            # Update using insertion, deletion, or substitution costs
            d[x, y] = min(
                d[x - 1, y] + 1,  # Deletion
                d[x, y - 1] + 1,  # Insertion
                d[x - 1, y - 1] + (s1[x - 1] != s2[y - 1])  # Substitution
            )

    # Find the minimum distance for aligning s2 anywhere in s1
    cdef int mini = d[1, yLen]
    for x in range(2, xLen + 1):
        mini = min(mini, d[x, yLen])

    return mini

# Wrappers to be able to test
cpdef int py_hamming_distance(str seq1, str seq2, int limit=0):
    return hamming_distance(seq1, seq2, limit)

cpdef int py_levenshtein_distance(str seq1, str seq2, int limit=0):
    return levenshtein_distance(seq1, seq2, limit)

cpdef int py_subglobal_distance(str s1, str s2):
    return subglobal_distance(s1, s2)
