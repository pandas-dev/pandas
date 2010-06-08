
cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cdef inline _checknull(object val):
    if isinstance(val, float):
        return val != val or val == INF or val == NEGINF
    else:
        return val is None

cpdef checknull(object val):
    return _checknull(val)

def isnullobj(ndarray input):
    cdef int i, length
    cdef object val
    cdef ndarray[npy_int8, ndim=1] result
    cdef flatiter iter

    length = PyArray_SIZE(input)

    result = <ndarray> np.zeros(length, dtype=np.int8)

    iter= PyArray_IterNew(input)

    for i from 0 <= i < length:
        val = PyArray_GETITEM(input, PyArray_ITER_DATA(iter))

        if _checknull(val):
            result[i] = 1

        PyArray_ITER_NEXT(iter)

    return result
