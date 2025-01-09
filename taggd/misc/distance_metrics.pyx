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
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        limit (int): The maximum distance limit. If the distance exceeds this limit,
            the calculation aborts and returns `limit + 1`.

    Returns:
        int: The Hamming distance between the sequences, or `limit + 1` if the limit is exceeded.
    """
    cdef int i = 0
    cdef int sum = 0
    for i in xrange(len(seq1)):
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
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        limit (int): The maximum distance limit. If the calculated distance exceeds
            this limit, the calculation aborts and returns `limit + 1`.

    Returns:
        int: The Levenshtein distance between the sequences, or `limit + 1` if the limit is exceeded.
    """
    cdef list one_ago = None
    cdef list this_row = range(1, len(seq2) + 1) + [0]
    cdef int x
    cdef int y
    cdef int del_cost
    cdef int add_cost
    cdef int sub_cost
    for x in xrange(len(seq1)):
        one_ago = this_row
        this_row = [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            del_cost = one_ago[y] + 1
            add_cost = this_row[y - 1] + 1
            sub_cost = one_ago[y - 1] + (seq1[x] != seq2[y])
            this_row[y] = min(del_cost, add_cost, sub_cost)
        if limit > 0 and x > limit and min(this_row) > limit:
            return limit + 1
    return this_row[len(seq2) - 1]

cdef int subglobal_distance(str s1, str s2):
    """
    Calculates the edit distance for a sub-global alignment of two sequences.

    The sub-global alignment computes the edit distance of `s2` aligned against `s1`,
    where mismatches and indels both score as 1. Overhanging parts of `s1` are not included in the distance calculation.

    Args:
        s1 (str): The longer sequence (probe).
        s2 (str): The shorter sequence (query).

    Returns:
        int: The minimum edit distance for aligning `s2` to `s1`.
    """
    cdef int xLen = len(s1)
    cdef int yLen = len(s2)
    if xLen < yLen:
        raise ValueError("Sub-global edit distance is undefined for sequences " \
                         "where the probe is shorter than the aligned sequence.")
    cdef int x
    cdef int y
    cdef np.ndarray[np.uint32_t, ndim=2] d = np.empty([xLen+1, yLen+1], dtype=np.uint32)

    # Initialize array
    for x in xrange(0, xLen+1):
        d[x,0] = 0
    for y in xrange(1, yLen+1):
        d[0,y] = y # To ensure all of s2 is spanned.

    # Perform DP.
    for x in xrange(1, xLen+1):
        # Fill matrix.
        for y in xrange(1, yLen+1):
            d[x,y] = min( min(d[x-1,y]+1, d[x,y-1]+1), d[x-1,y-1] + int(s1[x-1] != s2[y-1]) )

    # Find min for sub-global alignment so that all of s2 is covered,
    # but not necessarily all of s1 sequence.
    cdef int mini = 1000000
    cdef int iPos = 0
    cdef int i = xLen
    cdef int j
    while i > 0:
        if d[i,yLen] < mini:
            mini = d[i,yLen]
            iPos = i
        i -= 1

    # Return min distance
    return mini
