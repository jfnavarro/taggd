"""
Contains functions for demultiplexing records.
It uses a multi-thread approach where several file descriptors
are open and every i line of the input file is processed.
The resulting file are then merged
"""

import random
from cpython cimport bool
import taggd.core.match as match
from taggd.core.demultiplex_search_functions cimport get_candidates, get_distances, get_top_hits

cdef str trim_helpers(str seq, list trim_sequences):
    """Simply helper function to remove
    helper sequences from a barcode"""
    cdef int prev_start = 0
    cdef int prev_end = 0
    for start, end in trim_sequences:
        offset = prev_end - prev_start
        seq = seq[:(start-offset)] + seq[(end-offset):]
        prev_start = start
        prev_end = end
    return seq

cdef list demultiplex_record(object rec): # TODO pass needed arguments
    """
    Demultiplexes a record and returns a list of match objects
    (only more than one if ambiguous).
    """
    # Define local variables to speed up
    cdef str bcseq = None
    cdef int dist = 0
    cdef str read_barcode

    # Try perfect hit first.
    if barcode_tag is None:
        sequence = rec.sequence
    else:
        try:
            sequence = {tag:value for tag,value in rec.attributes["tags"]}[barcode_tag]
        except KeyError:
            raise ValueError(f"The specified SAM/BAM tag {barcode_tag} is not present for record {rec.annotation}")
    read_barcode = sequence[start_position:(start_position+barcode_length)]
    if trim_sequences is not None:
        read_barcode = trim_helpers(read_barcode)

    if read_barcode in true_barcodes:
        return [match.Match(rec, match.MATCHED_PERFECTLY, read_barcode, 0)]

    # Homopolymer filter.
    for filter in homopolymers:
        if filter in read_barcode:
            return [match.Match(rec, match.UNMATCHED, "-", -1)]

    # Include overhang.
    if pre_overhang != 0 or post_overhang != 0:
        read_barcode = sequence[(start_position - pre_overhang):min(len(sequence), \
                                    (start_position + barcode_length + post_overhang))]

    # Narrow down hits.
    cdef list candidates = get_candidates(read_barcode)
    cdef qual_hits = get_distances(read_barcode, candidates)
    cdef top_hits = get_top_hits(qual_hits)

    if not top_hits:
        # UNMATCHED
        return [match.Match(rec, match.UNMATCHED, "-", -1)]

    random.shuffle(top_hits)
    if len(top_hits) == 1 or multiple_hits_keep_one:
        # Unambiguous match
        bcseq, dist = top_hits[0]
        return [match.Match(rec, match.MATCHED_UNAMBIGUOUSLY, bcseq, dist)]

    # Add the rest as ambiguous match
    return [match.Match(rec, match.MATCHED_AMBIGUOUSLY, bcseq, dist) for bcseq,dist in top_hits]
