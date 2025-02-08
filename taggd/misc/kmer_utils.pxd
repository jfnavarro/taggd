from cpython cimport bool

cdef list get_kmers(str seq, int k, bool round_robin=?, int slider_increment=?)
