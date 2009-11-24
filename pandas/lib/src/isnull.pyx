
cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cdef inline _checknull(object val):
    return val is None or val != val or val == INF or val == NEGINF

cdef ndarray _isnullobj(input):
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

def isnull(input):
    '''
    Replacement for numpy.isnan / -numpy.isfinite which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
    cdef ndarray[npy_int8, ndim=1] result

    if isinstance(input, np.ndarray):
        if input.dtype.kind in ('O', 'S'):
            result = _isnullobj(input)

            return result.astype(np.bool)
        else:
            return -np.isfinite(input)
    else:
        return _checknull(input)

def notnull(input):
    '''
    Replacement for numpy.isfinite / -numpy.isnan which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
    if isinstance(input, np.ndarray):
        return -isnull(input)
    else:
        return not bool(_checknull(input))

