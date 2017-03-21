from util cimport numeric
from numpy cimport float64_t, double_t

cpdef numeric kth_smallest(numeric[:] a, Py_ssize_t k)

cdef inline Py_ssize_t swap(numeric *a, numeric *b) nogil except -1:
    cdef numeric t

    # cython doesn't allow pointer dereference so use array syntax
    t = a[0]
    a[0] = b[0]
    b[0] = t
    return 0


cdef inline kth_smallest_c(float64_t* a, Py_ssize_t k, Py_ssize_t n):
    cdef:
        Py_ssize_t i, j, l, m
        double_t x, t

    l = 0
    m = n -1
    while (l<m):
        x = a[k]
        i = l
        j = m

        while 1:
            while a[i] < x: i += 1
            while x < a[j]: j -= 1
            if i <= j:
                swap(&a[i], &a[j])
                i += 1; j -= 1

            if i > j: break

        if j < k: l = i
        if k < i: m = j
    return a[k]
