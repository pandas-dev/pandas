from libcpp.map cimport map
from numpy cimport ndarray, int32_t

ctypedef int32_t i4

def map_indices(ndarray[i4] values):
    cdef:
        i4 i, n
        map[i4, i4] mapping

    mapping = map[i4, i4]()

    n = len(values)
    for i in range(n):
        mapping[i] = values[i]
