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
        seqs (List[str]): The input sequences.
        k (int): The k-mer length.
        round_robin (bool): Whether to treat the sequences as circular.
        slider_increment (int): The step size for sliding the k-mer window.

    Returns:
        Dict[str, Dict[str, List[int]]]: A dictionary where:
            - The keys are k-mers.
            - The values are dictionaries mapping sequences to lists of k-mer offsets.
    """
    cdef object kmer2seq = defaultdict(default_dict_factory)
    cdef str seq
    cdef str seqq
    cdef str kmer
    cdef int i
    for seq in seqs:
        # Adjust barcode if round robin
        seqq = seq + seq[0:(k-1)] if round_robin else seq
        # Create the Kmers of length k with slider increment
        for i in range(0, len(seqq)-k+1, slider_increment):
            kmer = seqq[i:i+k]
            kmer2seq[kmer][seq].append(i)
        # Special treatment of last in case we would skip it during incrementation
        if len(seqq) % slider_increment != 0:
            i = len(seqq)-k
            kmer = seqq[i:len(seqq)]
            kmer2seq[kmer][seq].add(i)
    # Important to be able to generate KeyError
    kmer2seq.default_factory = None
    return kmer2seq

cdef list get_kmers(str seq, int k, bool round_robin=False, int slider_increment=0):
    """
    Generates the k-mers of a sequence as a list of k-mer and offset tuples.

    The last k-mer is always included, irrespective of the slider increment.

    Args:
        seq (str): The input sequence.
        k (int): The k-mer length.
        round_robin (bool): Whether to treat the sequence as circular.
        slider_increment (int): The step size for sliding the k-mer window.

    Returns:
        List[Tuple[str, int]]: A list of tuples where each tuple contains a k-mer
        and its offset in the sequence.
    """
    cdef list kmer_list = list()
    cdef str seqq = seq + seq[0:(k-1)] if round_robin else seq
    cdef str kmer
    cdef int i
    # Simply compute kmers for the sequence
    # TODO this function could be used in get_kmers_dictst to avoid code duplication
    for i in range(0, len(seqq)-k+1, slider_increment):
        kmer = seqq[i:i+k]
        kmer_list.append((kmer, i))
    # Special treatment of last
    if len(seqq) % slider_increment != 0:
        i = len(seqq)-k
        kmer = seqq[i:len(seqq)]
        kmer_list.append((kmer, i))
    return kmer_list
