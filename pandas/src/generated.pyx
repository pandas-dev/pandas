@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_float64(ndarray[float64_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[float64_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            outbuf[i] = NaN
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_object(ndarray[object] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[object] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            outbuf[i] = NaN
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_int32(ndarray[int32_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[int32_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            raise ValueError('No NA values allowed')
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_int64(ndarray[int64_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[int64_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            raise ValueError('No NA values allowed')
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_bool(ndarray[uint8_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[uint8_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            raise ValueError('No NA values allowed')
        else:
            outbuf[i] = values[idx]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]

        if idx == -1:
            for j from 0 <= j < k:
                outbuf[i, j] = NaN
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_object(ndarray[object, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]

        if idx == -1:
            for j from 0 <= j < k:
                outbuf[i, j] = NaN
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]

        if idx == -1:
            for j from 0 <= j < k:
                raise ValueError('No NA values allowed')
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]

        if idx == -1:
            for j from 0 <= j < k:
                raise ValueError('No NA values allowed')
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]

        if idx == -1:
            for j from 0 <= j < k:
                raise ValueError('No NA values allowed')
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for j from 0 <= j < k:
        idx = indexer[j]

        if idx == -1:
            for i from 0 <= i < n:
                outbuf[i, j] = NaN
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_object(ndarray[object, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for j from 0 <= j < k:
        idx = indexer[j]

        if idx == -1:
            for i from 0 <= i < n:
                outbuf[i, j] = NaN
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for j from 0 <= j < k:
        idx = indexer[j]

        if idx == -1:
            for i from 0 <= i < n:
                raise ValueError('No NA values allowed')
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for j from 0 <= j < k:
        idx = indexer[j]

        if idx == -1:
            for i from 0 <= i < n:
                raise ValueError('No NA values allowed')
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for j from 0 <= j < k:
        idx = indexer[j]

        if idx == -1:
            for i from 0 <= i < n:
                raise ValueError('No NA values allowed')
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]


