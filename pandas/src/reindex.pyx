def pad_inplace_float64(ndarray[float64_t] values,
                        ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be padded otherwise False

    e.g.
    pad_inplace_float64(values, isnull(values))
    '''
    cdef:
        Py_ssize_t i, n
        float64_t val

    n = len(values)
    val = NaN
    for i from 0 <= i < n:
        if mask[i]:
            values[i] = val
        else:
            val = values[i]

def get_pad_indexer(ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be padded otherwise False

    e.g.
    pad_inplace_float64(values, isnull(values))
    '''
    cdef:
        Py_ssize_t i, n
        int32_t idx
        ndarray[int32_t] indexer

    n = len(mask)
    indexer = np.empty(n, dtype=np.int32)

    idx = 0
    for i from 0 <= i < n:
        if not mask[i]:
            idx = i
        indexer[i] = idx

    return indexer

def get_backfill_indexer(ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be padded otherwise False

    e.g.
    pad_inplace_float64(values, isnull(values))
    '''
    cdef:
        Py_ssize_t i, n
        int32_t idx
        ndarray[int32_t] indexer

    n = len(mask)
    indexer = np.empty(n, dtype=np.int32)

    idx = n - 1
    i = n - 1
    while i >= 0:
        if not mask[i]:
            idx = i
        indexer[i] = idx
        i -= 1

    return indexer

def backfill_inplace_float64(ndarray[float64_t] values,
                             ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be backfilled otherwise False
    '''
    cdef:
        Py_ssize_t i, n
        float64_t val

    n = len(values)
    val = NaN
    i = n - 1
    while i >= 0:
        if mask[i]:
            values[i] = val
        else:
            val = values[i]
        i -= 1

def ordered_left_join(ndarray[object] left, ndarray[object] right):
    # cdef dict right_map = map_indices_buf(right)
    # return merge_indexer(left, right_map)
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[int32_t] indexer
        ndarray[uint8_t] mask
        object val

    i = 0
    j = 0
    n = len(left)
    k = len(right)

    indexer = np.zeros(n, dtype=np.int32)
    mask = np.ones(n, dtype=np.uint8)

    for i from 0 <= i < n:
        val = left[i]

        while j < k and right[j] < val:
            j += 1

        if j == k:
            break

        if val == right[j]:
            indexer[i] = j
            mask[i] = 0

    return indexer, mask.view(np.bool_)

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_2d(ndarray[int64_t] left, ndarray[int64_t] right,
                 ndarray[float64_t, ndim=2] lvalues,
                 ndarray[float64_t, ndim=2] rvalues,
                 ndarray[float64_t, ndim=2] out):
    cdef:
        Py_ssize_t i, j, k, nright, nleft, kright, kleft
        int64_t val

    nleft, kleft = (<object> lvalues).shape
    nright, kright = (<object> rvalues).shape

    j = 0
    for i from 0 <= i < nleft:
        for k from 0 <= k < kleft:
            out[i, k] = lvalues[i, k]

        val = left[i]

        while j < nright and right[j] < val:
            j += 1

        if j == nright:
            for k from kleft <= k < kleft + kright:
                out[i, k] = NaN
            continue

        if val == right[j]:
            for k from kleft <= k < kleft + kright:
                out[i, k] = rvalues[j, k - kleft]
        else:
            for k from kleft <= k < kleft + kright:
                out[i, k] = NaN

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_1d(ndarray[int64_t] left, ndarray[int64_t] right,
                 ndarray[float64_t] lvalues,
                 ndarray[float64_t] rvalues,
                 ndarray[float64_t, ndim=2] out):
    cdef:
        Py_ssize_t i, j, nright, nleft
        int64_t val

    nleft = len(lvalues)
    nright = len(rvalues)

    j = 0
    for i from 0 <= i < nleft:
        out[i, 0] = lvalues[i]

        val = left[i]

        while j < nright and right[j] < val:
            j += 1

        if j == nright:
            out[i, 1] = NaN
            continue

        if val == right[j]:
            out[i, 1] = rvalues[j]
        else:
            out[i, 1] = NaN


@cython.wraparound(False)
@cython.boundscheck(False)
def take_join_contiguous(ndarray[float64_t, ndim=2] lvalues,
                         ndarray[float64_t, ndim=2] rvalues,
                         ndarray[int32_t] lindexer,
                         ndarray[int32_t] rindexer,
                         ndarray out):
    cdef:
        Py_ssize_t i, j, rk, lk, n, lidx, ridx
        float64_t *outbuf

    assert(out.flags.contiguous)

    outbuf = <float64_t*> out.data

    n = len(lindexer)
    lk = lvalues.shape[1]
    rk = rvalues.shape[1]

    for i from 0 <= i < n:
        lidx = lindexer[i]
        ridx = rindexer[i]

        if lidx == -1:
            for j from 0 <= j < lk:
                outbuf[0] = NaN
                outbuf = outbuf + 1
        else:
            for j from 0 <= j < lk:
                outbuf[0] = lvalues[lidx, j]
                outbuf = outbuf + 1

        if lidx == -1:
            for j from 0 <= j < rk:
                outbuf[0] = NaN
                outbuf = outbuf + 1
        else:
            for j from 0 <= j < rk:
                outbuf[0] = rvalues[ridx, j]
                outbuf = outbuf + 1

@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_list(list values, dict oldMap):
    cdef int i, j, length, newLength
    cdef object idx
    cdef ndarray[int32_t] fill_vec

    newLength = len(values)
    fill_vec = np.empty(newLength, dtype=np.int32)
    for i from 0 <= i < newLength:
        idx = values[i]
        if idx in oldMap:
            fill_vec[i] = oldMap[idx]
        else:
            fill_vec[i] = -1

    return fill_vec
