cdef double __add(double a, double b):
    return a + b
cdef double __sub(double a, double b):
    return a - b
cdef double __div(double a, double b):
    return a / b
cdef double __mul(double a, double b):
    return a * b
cdef double __eq(double a, double b):
    return a == b
cdef double __ne(double a, double b):
    return a != b
cdef double __lt(double a, double b):
    return a < b
cdef double __gt(double a, double b):
    return a > b
cdef double __pow(double a, double b):
    return a ** b

ctypedef double (* double_func)(double a, double b)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef ndarray _applyFunc(double_func func, ndarray index, ndarray ao,
                        ndarray bo, dict aMap, dict bMap):
    '''
    C function taking a function pointer for quickly adding two Series objects.
    '''
    cdef ndarray result
    cdef double *result_data, *a_data, *b_data
    cdef int length
    cdef Py_ssize_t i, aidx, bidx
    cdef double nan
    cdef object idx
    cdef ndarray[object, ndim=1] ibuf
    cdef ndarray[double_t, ndim=1] A, B

    A = ao
    B = bo

    ibuf = index

    nan = <double> np.NaN
    length = len(index)
    result = <ndarray> np.empty(length, dtype=float)
    result_data = <double *> result.data

    for i from 0 <= i < length:
        idx = ibuf[i]

        if idx not in aMap or idx not in bMap:
            result_data[i] = nan
            continue

        aidx = aMap[idx]
        bidx = bMap[idx]
        result_data[i] = func(A[aidx], B[bidx])

    return result

def combineFunc(object name, ndarray index, ndarray ao,
                ndarray bo, dict aMap, dict bMap):
    '''
    Combine two series (values and index maps for each passed in) using the
    indicated function.
    '''
    if name == "__add__":
        return _applyFunc(__add, index, ao, bo, aMap, bMap)
    elif name == "__sub__":
        return _applyFunc(__sub, index, ao, bo, aMap, bMap)
    elif name == "__div__":
        return _applyFunc(__div, index, ao, bo, aMap, bMap)
    elif name == "__mul__":
        return _applyFunc(__mul, index, ao, bo, aMap, bMap)
    elif name == "__eq__":
        return _applyFunc(__eq, index, ao, bo, aMap, bMap)
    elif name == "__ne__":
        return _applyFunc(__ne, index, ao, bo, aMap, bMap)
    elif name == "__lt__":
        return _applyFunc(__lt, index, ao, bo, aMap, bMap)
    elif name == "__gt__":
        return _applyFunc(__gt, index, ao, bo, aMap, bMap)
    elif name == "__pow__":
        return _applyFunc(__pow, index, ao, bo, aMap, bMap)
    else:
        raise Exception('bad funcname requested of Cython code')
