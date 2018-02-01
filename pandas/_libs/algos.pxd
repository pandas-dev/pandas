from util cimport numeric
from numpy cimport float64_t, double_t

cpdef numeric kth_smallest(numeric[:] a, Py_ssize_t k) nogil

cdef inline Py_ssize_t swap(numeric *a, numeric *b) nogil:
    cdef numeric t

    # cython doesn't allow pointer dereference so use array syntax
    t = a[0]
    a[0] = b[0]
    b[0] = t
    return 0

cdef:
    int TIEBREAK_AVERAGE = 0
    int TIEBREAK_MIN = 1
    int TIEBREAK_MAX = 2
    int TIEBREAK_FIRST = 3
    int TIEBREAK_FIRST_DESCENDING = 4
    int TIEBREAK_DENSE = 5
