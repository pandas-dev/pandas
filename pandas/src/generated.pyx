
cimport numpy as np
cimport cython

from numpy cimport *

from cpython cimport (PyDict_New, PyDict_GetItem, PyDict_SetItem,
                      PyDict_Contains, PyDict_Keys,
                      Py_INCREF, PyTuple_SET_ITEM,
                      PyTuple_SetItem,
                      PyTuple_New)
from cpython cimport PyFloat_Check
cimport cpython

import numpy as np
isnan = np.isnan
cdef double NaN = <double> np.NaN
cdef double nan = NaN

from datetime import datetime as pydatetime

# this is our datetime.pxd
from datetime cimport *

from khash cimport *

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

ctypedef unsigned char UChar

cimport util
from util cimport is_array, _checknull, _checknan

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()
import_ufunc()

cdef int PLATFORM_INT = (<ndarray> np.arange(0, dtype=np.int_)).descr.type_num

cpdef ensure_platform_int(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == PLATFORM_INT:
            return arr
        else:
            return arr.astype(np.int_)
    else:
        return np.array(arr, dtype=np.int_)



cpdef ensure_float64(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_FLOAT64:
            return arr
        else:
            return arr.astype(np.float64)
    else:
        return np.array(arr, dtype=np.float64)


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


cpdef ensure_object(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_OBJECT:
            return arr
        else:
            return arr.astype(np.object_)
    else:
        return np.array(arr, dtype=np.object_)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_float64(ndarray[float64_t] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_object(ndarray[object] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_int32(ndarray[int32_t] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_int64(ndarray[int64_t] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_bool(ndarray[uint8_t] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
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
def pad_inplace_float64(ndarray[float64_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef float64_t val
    cdef int lim, fill_count = 0

    N = len(values)

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
def pad_inplace_object(ndarray[object] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef object val
    cdef int lim, fill_count = 0

    N = len(values)

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
def pad_inplace_int32(ndarray[int32_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef int32_t val
    cdef int lim, fill_count = 0

    N = len(values)

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
def pad_inplace_int64(ndarray[int64_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef int64_t val
    cdef int lim, fill_count = 0

    N = len(values)

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
def pad_inplace_bool(ndarray[uint8_t] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef uint8_t val
    cdef int lim, fill_count = 0

    N = len(values)

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
def backfill_inplace_float64(ndarray[float64_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef float64_t val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
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
def backfill_inplace_object(ndarray[object] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef object val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
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
def backfill_inplace_int32(ndarray[int32_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef int32_t val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
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
def backfill_inplace_int64(ndarray[int64_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef int64_t val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
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
def backfill_inplace_bool(ndarray[uint8_t] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef uint8_t val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
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
@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_object(ndarray[object, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef object val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

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
@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_int32(ndarray[int32_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int32_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

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
@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_int64(ndarray[int64_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

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
@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_bool(ndarray[uint8_t, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef uint8_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

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

@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_float64(ndarray[float64_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef float64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
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
def backfill_2d_inplace_object(ndarray[object, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef object val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
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
def backfill_2d_inplace_int32(ndarray[int32_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int32_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
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
def backfill_2d_inplace_int64(ndarray[int64_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef int64_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
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
def backfill_2d_inplace_bool(ndarray[uint8_t, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef uint8_t val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]

@cython.wraparound(False)
def take_1d_float64(ndarray[float64_t] values,
                     ndarray[int64_t] indexer,
                     out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[float64_t] outbuf
        float64_t fv

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    if False and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                raise ValueError('No NA values allowed')
            else:
                outbuf[i] = values[idx]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                outbuf[i] = fv
            else:
                outbuf[i] = values[idx]

@cython.wraparound(False)
def take_1d_object(ndarray[object] values,
                     ndarray[int64_t] indexer,
                     out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[object] outbuf
        object fv

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    if False and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                raise ValueError('No NA values allowed')
            else:
                outbuf[i] = values[idx]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                outbuf[i] = fv
            else:
                outbuf[i] = values[idx]

@cython.wraparound(False)
def take_1d_int32(ndarray[int32_t] values,
                     ndarray[int64_t] indexer,
                     out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[int32_t] outbuf
        int32_t fv

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                raise ValueError('No NA values allowed')
            else:
                outbuf[i] = values[idx]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                outbuf[i] = fv
            else:
                outbuf[i] = values[idx]

@cython.wraparound(False)
def take_1d_int64(ndarray[int64_t] values,
                     ndarray[int64_t] indexer,
                     out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[int64_t] outbuf
        int64_t fv

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                raise ValueError('No NA values allowed')
            else:
                outbuf[i] = values[idx]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                outbuf[i] = fv
            else:
                outbuf[i] = values[idx]

@cython.wraparound(False)
def take_1d_bool(ndarray[uint8_t] values,
                     ndarray[int64_t] indexer,
                     out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[uint8_t] outbuf
        uint8_t fv

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                raise ValueError('No NA values allowed')
            else:
                outbuf[i] = values[idx]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                outbuf[i] = fv
            else:
                outbuf[i] = values[idx]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_float64(ndarray[float64_t] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        float64_t prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique
@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_object(ndarray[object] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        object prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique
@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_int32(ndarray[int32_t] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        int32_t prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique
@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_int64(ndarray[int64_t] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        int64_t prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique
@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_bool(ndarray[uint8_t] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        uint8_t prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_float64(ndarray[float64_t] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_object(ndarray[object] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_int32(ndarray[int32_t] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_int64(ndarray[int64_t] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_bool(ndarray[uint8_t] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result


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
def arrmap_bool(ndarray[uint8_t] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if False and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    raise ValueError('No NA values allowed')
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_object(ndarray[object, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf
        object fv

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if False and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    raise ValueError('No NA values allowed')
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf
        int32_t fv

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    raise ValueError('No NA values allowed')
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf
        int64_t fv

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    raise ValueError('No NA values allowed')
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf
        uint8_t fv

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    raise ValueError('No NA values allowed')
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    outbuf[i, j] = values[idx, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf
        float64_t fv

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if False and _checknan(fill_value):
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    raise ValueError('No NA values allowed')
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]
    else:
        fv = fill_value
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    outbuf[i, j] = fv
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_object(ndarray[object, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf
        object fv

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if False and _checknan(fill_value):
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    raise ValueError('No NA values allowed')
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]
    else:
        fv = fill_value
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    outbuf[i, j] = fv
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf
        int32_t fv

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    raise ValueError('No NA values allowed')
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]
    else:
        fv = fill_value
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    outbuf[i, j] = fv
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf
        int64_t fv

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    raise ValueError('No NA values allowed')
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]
    else:
        fv = fill_value
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    outbuf[i, j] = fv
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int64_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf
        uint8_t fv

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if True and _checknan(fill_value):
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    raise ValueError('No NA values allowed')
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]
    else:
        fv = fill_value
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    outbuf[i, j] = fv
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int64_t] idx0,
                           ndarray[int64_t] idx1,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out


    if False and _checknan(fill_value):
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    raise ValueError('No NA values allowed')
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        raise ValueError('No NA values allowed')
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]
    else:
        fv = fill_value
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        outbuf[i, j] = fv
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_object(ndarray[object, ndim=2] values,
                           ndarray[int64_t] idx0,
                           ndarray[int64_t] idx1,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf
        object fv

    n = len(idx0)
    k = len(idx1)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out


    if False and _checknan(fill_value):
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    raise ValueError('No NA values allowed')
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        raise ValueError('No NA values allowed')
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]
    else:
        fv = fill_value
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        outbuf[i, j] = fv
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int64_t] idx0,
                           ndarray[int64_t] idx1,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf
        int32_t fv

    n = len(idx0)
    k = len(idx1)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out


    if True and _checknan(fill_value):
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    raise ValueError('No NA values allowed')
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        raise ValueError('No NA values allowed')
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]
    else:
        fv = fill_value
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        outbuf[i, j] = fv
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int64_t] idx0,
                           ndarray[int64_t] idx1,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf
        int64_t fv

    n = len(idx0)
    k = len(idx1)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out


    if True and _checknan(fill_value):
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    raise ValueError('No NA values allowed')
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        raise ValueError('No NA values allowed')
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]
    else:
        fv = fill_value
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        outbuf[i, j] = fv
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int64_t] idx0,
                           ndarray[int64_t] idx1,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf
        uint8_t fv

    n = len(idx0)
    k = len(idx1)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out


    if True and _checknan(fill_value):
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    raise ValueError('No NA values allowed')
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        raise ValueError('No NA values allowed')
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]
    else:
        fv = fill_value
        for i in range(n):
            idx = idx0[i]
            if idx == -1:
                for j in range(k):
                    outbuf[i, j] = fv
            else:
                for j in range(k):
                    if idx1[j] == -1:
                        outbuf[i, j] = fv
                    else:
                        outbuf[i, j] = values[idx, idx1[j]]


@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_unique_float64(ndarray[float64_t] left,
                                      ndarray[float64_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        float64_t lval, rval

    i = 0
    j = 0
    nleft = len(left)
    nright = len(right)

    indexer = np.empty(nleft, dtype=np.int64)
    while True:
        if i == nleft:
            break

        if j == nright:
            indexer[i] = -1
            i += 1
            continue

        rval = right[j]

        while i < nleft - 1 and left[i] == rval:
            indexer[i] = j
            i += 1

        if left[i] == right[j]:
            indexer[i] = j
            i += 1
            while i < nleft - 1 and left[i] == rval:
                indexer[i] = j
                i += 1
            j += 1
        elif left[i] > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_unique_object(ndarray[object] left,
                                      ndarray[object] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        object lval, rval

    i = 0
    j = 0
    nleft = len(left)
    nright = len(right)

    indexer = np.empty(nleft, dtype=np.int64)
    while True:
        if i == nleft:
            break

        if j == nright:
            indexer[i] = -1
            i += 1
            continue

        rval = right[j]

        while i < nleft - 1 and left[i] == rval:
            indexer[i] = j
            i += 1

        if left[i] == right[j]:
            indexer[i] = j
            i += 1
            while i < nleft - 1 and left[i] == rval:
                indexer[i] = j
                i += 1
            j += 1
        elif left[i] > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_unique_int32(ndarray[int32_t] left,
                                      ndarray[int32_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        int32_t lval, rval

    i = 0
    j = 0
    nleft = len(left)
    nright = len(right)

    indexer = np.empty(nleft, dtype=np.int64)
    while True:
        if i == nleft:
            break

        if j == nright:
            indexer[i] = -1
            i += 1
            continue

        rval = right[j]

        while i < nleft - 1 and left[i] == rval:
            indexer[i] = j
            i += 1

        if left[i] == right[j]:
            indexer[i] = j
            i += 1
            while i < nleft - 1 and left[i] == rval:
                indexer[i] = j
                i += 1
            j += 1
        elif left[i] > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_unique_int64(ndarray[int64_t] left,
                                      ndarray[int64_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        int64_t lval, rval

    i = 0
    j = 0
    nleft = len(left)
    nright = len(right)

    indexer = np.empty(nleft, dtype=np.int64)
    while True:
        if i == nleft:
            break

        if j == nright:
            indexer[i] = -1
            i += 1
            continue

        rval = right[j]

        while i < nleft - 1 and left[i] == rval:
            indexer[i] = j
            i += 1

        if left[i] == right[j]:
            indexer[i] = j
            i += 1
            while i < nleft - 1 and left[i] == rval:
                indexer[i] = j
                i += 1
            j += 1
        elif left[i] > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer



def left_join_indexer_float64(ndarray[float64_t] left,
                              ndarray[float64_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        float64_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[float64_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.float64)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    i += 1
                    count += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                count += 1
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer


def left_join_indexer_object(ndarray[object] left,
                              ndarray[object] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        object lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[object] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=object)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    i += 1
                    count += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                count += 1
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer


def left_join_indexer_int32(ndarray[int32_t] left,
                              ndarray[int32_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int32_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[int32_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.int32)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    i += 1
                    count += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                count += 1
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer


def left_join_indexer_int64(ndarray[int64_t] left,
                              ndarray[int64_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int64_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[int64_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.int64)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    i += 1
                    count += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                count += 1
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer


@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_float64(ndarray[float64_t] left,
                                ndarray[float64_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        float64_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[float64_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        count = nright
    elif nright == 0:
        count = nleft
    else:
        while True:
            if i == nleft:
                count += nright - j
                break
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                count += 1
                j += 1

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.float64)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        for j in range(nright):
            lindexer[j] = -1
            rindexer[j] = j
            result[j] = right[j]
    elif nright == 0:
        for i in range(nright):
            lindexer[i] = i
            rindexer[i] = -1
            result[i] = left[i]
    else:
        while True:
            if i == nleft:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    count += 1
                    j += 1
                break
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    count += 1
                    i += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                count += 1
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                count += 1
                j += 1

    return result, lindexer, rindexer

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_object(ndarray[object] left,
                                ndarray[object] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        object lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[object] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        count = nright
    elif nright == 0:
        count = nleft
    else:
        while True:
            if i == nleft:
                count += nright - j
                break
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                count += 1
                j += 1

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=object)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        for j in range(nright):
            lindexer[j] = -1
            rindexer[j] = j
            result[j] = right[j]
    elif nright == 0:
        for i in range(nright):
            lindexer[i] = i
            rindexer[i] = -1
            result[i] = left[i]
    else:
        while True:
            if i == nleft:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    count += 1
                    j += 1
                break
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    count += 1
                    i += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                count += 1
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                count += 1
                j += 1

    return result, lindexer, rindexer

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_int32(ndarray[int32_t] left,
                                ndarray[int32_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        int32_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[int32_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        count = nright
    elif nright == 0:
        count = nleft
    else:
        while True:
            if i == nleft:
                count += nright - j
                break
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                count += 1
                j += 1

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.int32)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        for j in range(nright):
            lindexer[j] = -1
            rindexer[j] = j
            result[j] = right[j]
    elif nright == 0:
        for i in range(nright):
            lindexer[i] = i
            rindexer[i] = -1
            result[i] = left[i]
    else:
        while True:
            if i == nleft:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    count += 1
                    j += 1
                break
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    count += 1
                    i += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                count += 1
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                count += 1
                j += 1

    return result, lindexer, rindexer

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_int64(ndarray[int64_t] left,
                                ndarray[int64_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        int64_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[int64_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        count = nright
    elif nright == 0:
        count = nleft
    else:
        while True:
            if i == nleft:
                count += nright - j
                break
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                count += 1
                j += 1

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.int64)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        for j in range(nright):
            lindexer[j] = -1
            rindexer[j] = j
            result[j] = right[j]
    elif nright == 0:
        for i in range(nright):
            lindexer[i] = i
            rindexer[i] = -1
            result[i] = left[i]
    else:
        while True:
            if i == nleft:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    count += 1
                    j += 1
                break
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    count += 1
                    i += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                count += 1
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                count += 1
                j += 1

    return result, lindexer, rindexer


@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_float64(ndarray[float64_t] left,
                              ndarray[float64_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        float64_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[float64_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.float64)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = rval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_object(ndarray[object] left,
                              ndarray[object] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        object lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[object] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=object)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = rval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_int32(ndarray[int32_t] left,
                              ndarray[int32_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int32_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[int32_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.int32)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = rval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_int64(ndarray[int64_t] left,
                              ndarray[int64_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int64_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[int64_t] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=np.int64)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = rval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer


