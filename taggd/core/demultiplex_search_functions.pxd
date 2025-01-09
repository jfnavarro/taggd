from cpython cimport bool

cdef list get_candidates(
    str read_barcode,
    int k,
    int slider_increment,
    object kmer2seq,
    bool no_offset_speedup,
    int pre_overhang,
    int post_overhang,
    int max_edit_distance
)

cdef list get_distances(str read_barcode, list candidates, int metric_choice, int max_edit_distance)

cdef list get_top_hits(list qual_hits, float ambiguity_factor)
