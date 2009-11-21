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

cdef ndarray _applyFunc(double_func func, ndarray index, object ao,
                        object bo, dict aMap, dict bMap):
    '''
    C function taking a function pointer for quickly adding two Series objects.
    '''
    cdef ndarray A, B, result
    cdef double *result_data
    cdef int i, length
    cdef flatiter itera, iterb, iteridx
    cdef double nan
    cdef object idx

    # This is EXTREMELY important, otherwise you will get very
    # undesired results
    A = PyArray_ContiguousFromAny(ao, NPY_DOUBLE, 1, 1)
    B = PyArray_ContiguousFromAny(bo, NPY_DOUBLE, 1, 1)

    nan = <double> np.NaN
    length = PyArray_SIZE(index)

    result = <ndarray> np.empty(length, np.float64)
    result_data = <double *>result.data

    itera = <flatiter> PyArray_IterNew(A)
    iterb = <flatiter> PyArray_IterNew(B)
    iteridx = PyArray_IterNew(index)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, <void *> iteridx.dataptr)
        PyArray_ITER_NEXT(iteridx)

        if idx not in aMap or idx not in bMap:
            result_data[i] = nan
            continue

        result_data[i] = func((<double *>A.data)[aMap[idx]],
                            (<double *>B.data)[bMap[idx]])

    return result

def combineFunc(object name, ndarray index, object ao,
                object bo, dict aMap, dict bMap):
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
