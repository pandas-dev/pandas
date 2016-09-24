"""
Template for each `dtype` helper function using 1-d template

# 1-d template
- map_indices
- pad
- pad_1d
- pad_2d
- backfill
- backfill_1d
- backfill_2d
- is_monotonic
- arrmap

WARNING: DO NOT edit .pxi FILE directly, .pxi is generated from .pxi.in
"""

#----------------------------------------------------------------------
# 1-d template
#----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_float64(ndarray[float64_t] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_float64(ndarray[float64_t] old, ndarray[float64_t] new,
                 limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef float64_t cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_float64(ndarray[float64_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef float64_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_float64(ndarray[float64_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef float64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

"""
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
"""


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_float64(ndarray[float64_t] old, ndarray[float64_t] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef float64_t cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_float64(ndarray[float64_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef float64_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1, -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_float64(ndarray[float64_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef float64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1, -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_float64(ndarray[float64_t] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        float64_t prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and arr[0] == iNaT:
        return False, False, True

    with nogil:
        prev = arr[0]
        for i in range(1, n):
            cur = arr[i]
            if timelike and cur == iNaT:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if cur < prev:
                is_monotonic_inc = 0
            elif cur > prev:
                is_monotonic_dec = 0
            elif cur == prev:
                is_unique = 0
            else:
                # cur or prev is NaN
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if not is_monotonic_inc and not is_monotonic_dec:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            prev = cur
    return is_monotonic_inc, is_monotonic_dec, \
           is_unique and (is_monotonic_inc or is_monotonic_dec)


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_float64(ndarray[float64_t] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_float32(ndarray[float32_t] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_float32(ndarray[float32_t] old, ndarray[float32_t] new,
                 limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef float32_t cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_float32(ndarray[float32_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef float32_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_float32(ndarray[float32_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef float32_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

"""
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
"""


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_float32(ndarray[float32_t] old, ndarray[float32_t] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef float32_t cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_float32(ndarray[float32_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef float32_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1, -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_float32(ndarray[float32_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef float32_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1, -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_float32(ndarray[float32_t] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        float32_t prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and arr[0] == iNaT:
        return False, False, True

    with nogil:
        prev = arr[0]
        for i in range(1, n):
            cur = arr[i]
            if timelike and cur == iNaT:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if cur < prev:
                is_monotonic_inc = 0
            elif cur > prev:
                is_monotonic_dec = 0
            elif cur == prev:
                is_unique = 0
            else:
                # cur or prev is NaN
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if not is_monotonic_inc and not is_monotonic_dec:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            prev = cur
    return is_monotonic_inc, is_monotonic_dec, \
           is_unique and (is_monotonic_inc or is_monotonic_dec)


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_float32(ndarray[float32_t] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_object(ndarray[object] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_object(ndarray[object] old, ndarray[object] new,
                 limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef object cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_object(ndarray[object] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef object val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_object(ndarray[object, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef object val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

"""
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
"""


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_object(ndarray[object] old, ndarray[object] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef object cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_object(ndarray[object] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef object val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1, -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_object(ndarray[object, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef object val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1, -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_object(ndarray[object] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        object prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and arr[0] == iNaT:
        return False, False, True

    
    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if timelike and cur == iNaT:
            is_monotonic_inc = 0
            is_monotonic_dec = 0
            break
        if cur < prev:
            is_monotonic_inc = 0
        elif cur > prev:
            is_monotonic_dec = 0
        elif cur == prev:
            is_unique = 0
        else:
            # cur or prev is NaN
            is_monotonic_inc = 0
            is_monotonic_dec = 0
            break
        if not is_monotonic_inc and not is_monotonic_dec:
            is_monotonic_inc = 0
            is_monotonic_dec = 0
            break
        prev = cur
    return is_monotonic_inc, is_monotonic_dec, \
           is_unique and (is_monotonic_inc or is_monotonic_dec)


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_object(ndarray[object] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_int32(ndarray[int32_t] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_int32(ndarray[int32_t] old, ndarray[int32_t] new,
                 limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef int32_t cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_int32(ndarray[int32_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef int32_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_int32(ndarray[int32_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int32_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

"""
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
"""


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_int32(ndarray[int32_t] old, ndarray[int32_t] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef int32_t cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_int32(ndarray[int32_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef int32_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1, -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_int32(ndarray[int32_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int32_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1, -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_int32(ndarray[int32_t] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        int32_t prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and arr[0] == iNaT:
        return False, False, True

    with nogil:
        prev = arr[0]
        for i in range(1, n):
            cur = arr[i]
            if timelike and cur == iNaT:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if cur < prev:
                is_monotonic_inc = 0
            elif cur > prev:
                is_monotonic_dec = 0
            elif cur == prev:
                is_unique = 0
            else:
                # cur or prev is NaN
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if not is_monotonic_inc and not is_monotonic_dec:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            prev = cur
    return is_monotonic_inc, is_monotonic_dec, \
           is_unique and (is_monotonic_inc or is_monotonic_dec)


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_int32(ndarray[int32_t] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_int64(ndarray[int64_t] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_int64(ndarray[int64_t] old, ndarray[int64_t] new,
                 limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef int64_t cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_int64(ndarray[int64_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef int64_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_int64(ndarray[int64_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

"""
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
"""


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_int64(ndarray[int64_t] old, ndarray[int64_t] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef int64_t cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_int64(ndarray[int64_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef int64_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1, -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_int64(ndarray[int64_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1, -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_int64(ndarray[int64_t] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        int64_t prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and arr[0] == iNaT:
        return False, False, True

    with nogil:
        prev = arr[0]
        for i in range(1, n):
            cur = arr[i]
            if timelike and cur == iNaT:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if cur < prev:
                is_monotonic_inc = 0
            elif cur > prev:
                is_monotonic_dec = 0
            elif cur == prev:
                is_unique = 0
            else:
                # cur or prev is NaN
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if not is_monotonic_inc and not is_monotonic_dec:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            prev = cur
    return is_monotonic_inc, is_monotonic_dec, \
           is_unique and (is_monotonic_inc or is_monotonic_dec)


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_int64(ndarray[int64_t] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_bool(ndarray[uint8_t] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_bool(ndarray[uint8_t] old, ndarray[uint8_t] new,
                 limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef uint8_t cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_bool(ndarray[uint8_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef uint8_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_bool(ndarray[uint8_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef uint8_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

"""
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
"""


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_bool(ndarray[uint8_t] old, ndarray[uint8_t] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef uint8_t cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_bool(ndarray[uint8_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef uint8_t val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1, -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_bool(ndarray[uint8_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef uint8_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1, -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_bool(ndarray[uint8_t] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        uint8_t prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and arr[0] == iNaT:
        return False, False, True

    with nogil:
        prev = arr[0]
        for i in range(1, n):
            cur = arr[i]
            if timelike and cur == iNaT:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if cur < prev:
                is_monotonic_inc = 0
            elif cur > prev:
                is_monotonic_dec = 0
            elif cur == prev:
                is_unique = 0
            else:
                # cur or prev is NaN
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            if not is_monotonic_inc and not is_monotonic_dec:
                is_monotonic_inc = 0
                is_monotonic_dec = 0
                break
            prev = cur
    return is_monotonic_inc, is_monotonic_dec, \
           is_unique and (is_monotonic_inc or is_monotonic_dec)


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_bool(ndarray[uint8_t] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)

#----------------------------------------------------------------------
# put template
#----------------------------------------------------------------------


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_float64(ndarray[float64_t, ndim=2] arr,
                     ndarray[float64_t, ndim=2] out,
                     Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]


def put2d_float64_float64(ndarray[float64_t, ndim=2, cast=True] values,
                                 ndarray[int64_t] indexer, Py_ssize_t loc,
                                 ndarray[float64_t] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_float32(ndarray[float32_t, ndim=2] arr,
                     ndarray[float32_t, ndim=2] out,
                     Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]


def put2d_float32_float32(ndarray[float32_t, ndim=2, cast=True] values,
                                 ndarray[int64_t] indexer, Py_ssize_t loc,
                                 ndarray[float32_t] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_int8(ndarray[int8_t, ndim=2] arr,
                     ndarray[float32_t, ndim=2] out,
                     Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]


def put2d_int8_float32(ndarray[int8_t, ndim=2, cast=True] values,
                                 ndarray[int64_t] indexer, Py_ssize_t loc,
                                 ndarray[float32_t] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_int16(ndarray[int16_t, ndim=2] arr,
                     ndarray[float32_t, ndim=2] out,
                     Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]


def put2d_int16_float32(ndarray[int16_t, ndim=2, cast=True] values,
                                 ndarray[int64_t] indexer, Py_ssize_t loc,
                                 ndarray[float32_t] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_int32(ndarray[int32_t, ndim=2] arr,
                     ndarray[float64_t, ndim=2] out,
                     Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]


def put2d_int32_float64(ndarray[int32_t, ndim=2, cast=True] values,
                                 ndarray[int64_t] indexer, Py_ssize_t loc,
                                 ndarray[float64_t] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_int64(ndarray[int64_t, ndim=2] arr,
                     ndarray[float64_t, ndim=2] out,
                     Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]


def put2d_int64_float64(ndarray[int64_t, ndim=2, cast=True] values,
                                 ndarray[int64_t] indexer, Py_ssize_t loc,
                                 ndarray[float64_t] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]

#----------------------------------------------------------------------
# ensure_dtype
#----------------------------------------------------------------------

cdef int PLATFORM_INT = (<ndarray> np.arange(0, dtype=np.intp)).descr.type_num

cpdef ensure_platform_int(object arr):
    # GH3033, GH1392
    # platform int is the size of the int pointer, e.g. np.intp
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == PLATFORM_INT:
            return arr
        else:
            return arr.astype(np.intp)
    else:
        return np.array(arr, dtype=np.intp)

cpdef ensure_object(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_OBJECT:
            return arr
        else:
            return arr.astype(np.object_)
    elif hasattr(arr, 'asobject'):
        return arr.asobject
    else:
        return np.array(arr, dtype=np.object_)

cpdef ensure_float64(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_FLOAT64:
            return arr
        else:
            return arr.astype(np.float64)
    else:
        return np.array(arr, dtype=np.float64)

cpdef ensure_float32(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_FLOAT32:
            return arr
        else:
            return arr.astype(np.float32)
    else:
        return np.array(arr, dtype=np.float32)

cpdef ensure_int8(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_INT8:
            return arr
        else:
            return arr.astype(np.int8)
    else:
        return np.array(arr, dtype=np.int8)

cpdef ensure_int16(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_INT16:
            return arr
        else:
            return arr.astype(np.int16)
    else:
        return np.array(arr, dtype=np.int16)

cpdef ensure_int32(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_INT32:
            return arr
        else:
            return arr.astype(np.int32)
    else:
        return np.array(arr, dtype=np.int32)

cpdef ensure_int64(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_INT64:
            return arr
        else:
            return arr.astype(np.int64)
    else:
        return np.array(arr, dtype=np.int64)
