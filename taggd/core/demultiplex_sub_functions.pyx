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
import taggd.constants as constants

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

cpdef list demultiplex_record(
    object rec,
    str barcode_tag,
    object true_barcodes,
    list homopolymers,
    int start_position,
    int barcode_length,
    list trim_sequences,
    int pre_overhang,
    int post_overhang,
    int k,
    object kmer2seq,
    int metric_choice,
    int max_edit_distance,
    float ambiguity_factor,
    int slider_increment,
    bool no_offset_speedup,
    bool multiple_hits_keep_one
):
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
            sequence = {tag: value for tag, value in rec.attributes["tags"]}[barcode_tag]
        except KeyError:
            raise ValueError(f"The specified SAM/BAM tag {barcode_tag} is not present for record {rec.annotation}")
    read_barcode = sequence[start_position:(start_position+barcode_length)]
    if trim_sequences is not None:
        read_barcode = trim_helpers(read_barcode, trim_sequences)

    if read_barcode in true_barcodes:
        return [match.Match(constants.MATCHED_PERFECTLY, read_barcode, 0)]

    # Homopolymer filter.
    for filter in homopolymers:
        if filter in read_barcode:
            return [match.Match(constants.UNMATCHED, "-", -1)]

    # Include overhang.
    if pre_overhang != 0 or post_overhang != 0:
        start = start_position - pre_overhang
        stop = min(len(sequence), start_position + barcode_length + post_overhang)
        read_barcode = sequence[start:stop]

    # Narrow down hits.
    cdef list candidates = get_candidates(
        read_barcode,
        k,
        slider_increment,
        kmer2seq,
        no_offset_speedup,
        pre_overhang,
        post_overhang,
        max_edit_distance
    )
    cdef qual_hits = get_distances(read_barcode, candidates, metric_choice, max_edit_distance)
    cdef top_hits = get_top_hits(qual_hits, ambiguity_factor)

    if not top_hits:
        # UNMATCHED
        return [match.Match(constants.UNMATCHED, "-", -1)]

    random.shuffle(top_hits)
    if len(top_hits) == 1 or multiple_hits_keep_one:
        # Unambiguous match
        bcseq, dist = top_hits[0]
        return [match.Match(constants.MATCHED_UNAMBIGUOUSLY, bcseq, dist)]

    # Add the rest as ambiguous match
    return [match.Match(constants.MATCHED_AMBIGUOUSLY, bcseq, dist) for bcseq, dist in top_hits]
