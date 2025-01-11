# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
"""
Contains utilities for working with k-mer chunks of the barcodes or other sequences.
"""
from cpython cimport bool
from collections import defaultdict


def default_dict_factory():
    return defaultdict(list)


cpdef object get_kmers_dicts(list seqs, int k, bool round_robin=False, int slider_increment=1):
    """
    Generates dictionaries for k-mers of a list of sequences.

    The last k-mer is always included, irrespective of the slider increment.

    Args:
        seqs: The input sequences.
        k: The k-mer length.
        round_robin: Whether to treat the sequences as circular.
        slider_increment: The step size for sliding the k-mer window.

    Returns:
        A dictionary where:
            - The keys are k-mers.
            - The values are dictionaries mapping sequences to lists of k-mer offsets.
    """
    cdef object kmer2seq = defaultdict(default_dict_factory)
    cdef list kmers_offsets
    cdef str seq
    cdef str kmer
    cdef int offset

    # Iterate through each sequence and generate k-mers using get_kmers
    for seq in seqs:
        kmers_offsets = get_kmers(seq, k, round_robin, slider_increment)
        for kmer, offset in kmers_offsets:
            kmer2seq[kmer][seq].append(offset)

    # Important to be able to generate KeyError
    kmer2seq.default_factory = None
    return kmer2seq


cdef list get_kmers(str seq, int k, bool round_robin=False, int slider_increment=1):
    """
    Generates the k-mers of a sequence as a list of k-mer and offset tuples.

    The last k-mer is always included, irrespective of the slider increment.

    Args:
        seq: The input sequence.
        k: The k-mer length.
        round_robin: Whether to treat the sequence as circular.
        slider_increment: The step size for sliding the k-mer window.

    Returns:
        A list of tuples where each tuple contains a k-mer
        and its offset in the sequence.
    """
    cdef list kmer_list = list()
    cdef str seqq = seq + seq[0:k-1] if round_robin else seq
    cdef str kmer
    cdef int i

    # Simply compute kmers for the sequence
    for i in range(0, len(seqq)-k+1, slider_increment):
        kmer = seqq[i:i+k]
        kmer_list.append((kmer, i))

    # Special treatment of last
    if len(seqq) % slider_increment != 0:
        i = len(seqq)-k
        kmer = seqq[i:len(seqq)]
        kmer_list.append((kmer, i))

    return kmer_list


# Wrapper for tests
cpdef list py_get_kmers(str seq, int k, bool round_robin=False, int slider_increment=1):
    return get_kmers(seq, k, round_robin, slider_increment)
