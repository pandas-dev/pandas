
cimport numpy as np
cimport cython

from libc.string cimport memmove

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

from datetime import datetime as pydatetime

# this is our datetime.pxd
from datetime cimport *

from khash cimport *

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
cpdef map_indices_float32(ndarray[float32_t] index):
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
def is_monotonic_float32(ndarray[float32_t] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        float32_t prev, cur
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

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

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
def groupby_float32(ndarray[float32_t] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

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

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

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

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

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

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

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

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

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
def take_1d_bool_bool(ndarray[uint8_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[uint8_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        uint8_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_bool_object(ndarray[uint8_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[object] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        object fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = True if values[idx] > 0 else False

@cython.wraparound(False)
def take_1d_int8_int8(ndarray[int8_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int8_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int8_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int8_int32(ndarray[int8_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int32_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int32_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int8_int64(ndarray[int8_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int8_float64(ndarray[int8_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int16_int16(ndarray[int16_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int16_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int16_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int16_int32(ndarray[int16_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int32_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int32_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int16_int64(ndarray[int16_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int16_float64(ndarray[int16_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int32_int32(ndarray[int32_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int32_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int32_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int32_int64(ndarray[int32_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int32_float64(ndarray[int32_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int64_int64(ndarray[int64_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[int64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        int64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_int64_float64(ndarray[int64_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_float32_float32(ndarray[float32_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float32_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float32_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_float32_float64(ndarray[float32_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_float64_float64(ndarray[float64_t] values,
                              ndarray[int64_t] indexer,
                              ndarray[float64_t] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        float64_t fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]

@cython.wraparound(False)
def take_1d_object_object(ndarray[object] values,
                              ndarray[int64_t] indexer,
                              ndarray[object] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        object fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = values[idx]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_bool_bool(ndarray[uint8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[uint8_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        uint8_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            uint8_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(uint8_t) and
            sizeof(uint8_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(uint8_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_bool_object(ndarray[uint8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[object, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        object fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            object *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(object) and
            sizeof(object) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(object) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = True if values[idx, j] > 0 else False

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int8_int8(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int8_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int8_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            int8_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int8_t) and
            sizeof(int8_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int8_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int8_int32(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int32_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            int32_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int32_t) and
            sizeof(int32_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int32_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int8_int64(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            int64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int8_float64(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int16_int16(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int16_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int16_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            int16_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int16_t) and
            sizeof(int16_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int16_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int16_int32(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int32_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            int32_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int32_t) and
            sizeof(int32_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int32_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int16_int64(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            int64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int16_float64(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int32_int32(ndarray[int32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int32_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            int32_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int32_t) and
            sizeof(int32_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int32_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int32_int64(ndarray[int32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            int64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int32_float64(ndarray[int32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int64_int64(ndarray[int64_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            int64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(int64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int64_float64(ndarray[int64_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_float32_float32(ndarray[float32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float32_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            float32_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float32_t) and
            sizeof(float32_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float32_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_float32_float64(ndarray[float32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_float64_float64(ndarray[float64_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF True:
        cdef:
            float64_t *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(float64_t) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_object_object(ndarray[object, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[object, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        object fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF False:
        cdef:
            object *v, *o

        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(object) and
            sizeof(object) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(object) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = values[idx, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_bool_bool(ndarray[uint8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[uint8_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        uint8_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            uint8_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(uint8_t) and
            sizeof(uint8_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(uint8_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_bool_object(ndarray[uint8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[object, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        object fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            object *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(object) and
            sizeof(object) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(object) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = True if values[i, idx] > 0 else False

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int8_int8(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int8_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int8_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            int8_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int8_t) and
            sizeof(int8_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int8_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int8_int32(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int32_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            int32_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int32_t) and
            sizeof(int32_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int32_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int8_int64(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            int64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int8_float64(ndarray[int8_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int16_int16(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int16_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int16_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            int16_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int16_t) and
            sizeof(int16_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int16_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int16_int32(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int32_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            int32_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int32_t) and
            sizeof(int32_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int32_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int16_int64(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            int64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int16_float64(ndarray[int16_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int32_int32(ndarray[int32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int32_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            int32_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int32_t) and
            sizeof(int32_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int32_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int32_int64(ndarray[int32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            int64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int32_float64(ndarray[int32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int64_int64(ndarray[int64_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        int64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            int64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(int64_t) and
            sizeof(int64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(int64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int64_float64(ndarray[int64_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_float32_float32(ndarray[float32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float32_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            float32_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float32_t) and
            sizeof(float32_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float32_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_float32_float64(ndarray[float32_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            float64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_float64_float64(ndarray[float64_t, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        float64_t fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF True:
        cdef:
            float64_t *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(float64_t) and
            sizeof(float64_t) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(float64_t) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_object_object(ndarray[object, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[object, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        object fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF False:
        cdef:
            object *v, *o

        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(object) and
            sizeof(object) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(object) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = values[i, idx]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_bool_bool(ndarray[uint8_t, ndim=2] values,
                                    indexer,
                                    ndarray[uint8_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        uint8_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_bool_object(ndarray[uint8_t, ndim=2] values,
                                    indexer,
                                    ndarray[object, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        object fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = True if values[idx, idx1[j]] > 0 else False

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int8_int8(ndarray[int8_t, ndim=2] values,
                                    indexer,
                                    ndarray[int8_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int8_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int8_int32(ndarray[int8_t, ndim=2] values,
                                    indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int32_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int8_int64(ndarray[int8_t, ndim=2] values,
                                    indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int8_float64(ndarray[int8_t, ndim=2] values,
                                    indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int16_int16(ndarray[int16_t, ndim=2] values,
                                    indexer,
                                    ndarray[int16_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int16_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int16_int32(ndarray[int16_t, ndim=2] values,
                                    indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int32_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int16_int64(ndarray[int16_t, ndim=2] values,
                                    indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int16_float64(ndarray[int16_t, ndim=2] values,
                                    indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int32_int32(ndarray[int32_t, ndim=2] values,
                                    indexer,
                                    ndarray[int32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int32_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int32_int64(ndarray[int32_t, ndim=2] values,
                                    indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int32_float64(ndarray[int32_t, ndim=2] values,
                                    indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int64_int64(ndarray[int64_t, ndim=2] values,
                                    indexer,
                                    ndarray[int64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        int64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_int64_float64(ndarray[int64_t, ndim=2] values,
                                    indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_float32_float32(ndarray[float32_t, ndim=2] values,
                                    indexer,
                                    ndarray[float32_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float32_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_float32_float64(ndarray[float32_t, ndim=2] values,
                                    indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_float64_float64(ndarray[float64_t, ndim=2] values,
                                    indexer,
                                    ndarray[float64_t, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        float64_t fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_object_object(ndarray[object, ndim=2] values,
                                    indexer,
                                    ndarray[object, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        object fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = values[idx, idx1[j]]


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

@cython.wraparound(False)
@cython.wraparound(False)
def group_last_float64(ndarray[float64_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    for i in range(N):
        lab = labels[i]
        if lab < 0:
            continue

        counts[lab] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[lab, j] += 1
                resx[lab, j] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
@cython.wraparound(False)
@cython.wraparound(False)
def group_last_float32(ndarray[float32_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float32_t, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    for i in range(N):
        lab = labels[i]
        if lab < 0:
            continue

        counts[lab] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[lab, j] += 1
                resx[lab, j] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]

@cython.wraparound(False)
@cython.wraparound(False)
def group_last_bin_float64(ndarray[float64_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] resx, nobs

    nobs = np.zeros_like(out)
    resx = np.empty_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    for i in range(N):
        while b < ngroups - 1 and i >= bins[b]:
            b += 1

        counts[b] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[b, j] += 1
                resx[b, j] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
@cython.wraparound(False)
@cython.wraparound(False)
def group_last_bin_float32(ndarray[float32_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float32_t, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        ndarray[float32_t, ndim=2] resx, nobs

    nobs = np.zeros_like(out)
    resx = np.empty_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    for i in range(N):
        while b < ngroups - 1 and i >= bins[b]:
            b += 1

        counts[b] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[b, j] += 1
                resx[b, j] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_float64(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    for i in range(N):
        lab = labels[i]
        if lab < 0:
            continue

        counts[lab] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[lab, j] += 1
                if nobs[lab, j] == rank:
                    resx[lab, j] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_float32(ndarray[float32_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float32_t, ndim=2] values,
              ndarray[int64_t] labels, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    for i in range(N):
        lab = labels[i]
        if lab < 0:
            continue

        counts[lab] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[lab, j] += 1
                if nobs[lab, j] == rank:
                    resx[lab, j] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_bin_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] resx, nobs

    nobs = np.zeros_like(out)
    resx = np.empty_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    for i in range(N):
        while b < ngroups - 1 and i >= bins[b]:
            b += 1

        counts[b] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[b, j] += 1
                if nobs[b, j] == rank:
                    resx[b, j] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_bin_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] bins, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        ndarray[float32_t, ndim=2] resx, nobs

    nobs = np.zeros_like(out)
    resx = np.empty_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    for i in range(N):
        while b < ngroups - 1 and i >= bins[b]:
            b += 1

        counts[b] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[b, j] += 1
                if nobs[b, j] == rank:
                    resx[b, j] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_float64(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j]
@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_float32(ndarray[float32_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float32_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_bin_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b, nbins
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j]
@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_bin_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b, nbins
        float32_t val, count
        ndarray[float32_t, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_float64(ndarray[float64_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] prodx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    prodx[lab, j] *= val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                prodx[lab, 0] *= val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = prodx[i, j]
@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_float32(ndarray[float32_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float32_t, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] prodx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    prodx[lab, j] *= val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                prodx[lab, 0] *= val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = prodx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_bin_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] prodx, nobs

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    prodx[b, j] *= val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                prodx[b, 0] *= val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = prodx[i, j]
@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_bin_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        ndarray[float32_t, ndim=2] prodx, nobs

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    prodx[b, j] *= val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                prodx[b, 0] *= val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = prodx[i, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def group_var_float64(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, ct
        ndarray[float64_t, ndim=2] nobs, sumx, sumxx

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):

            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
                    sumxx[lab, j] += val * val
    else:
        for i in range(N):

            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]
            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val
                sumxx[lab, 0] += val * val


    for i in range(len(counts)):
        for j in range(K):
            ct = nobs[i, j]
            if ct < 2:
                out[i, j] = nan
            else:
                out[i, j] = ((ct * sumxx[i, j] - sumx[i, j] * sumx[i, j]) /
                             (ct * ct - ct))
@cython.wraparound(False)
@cython.boundscheck(False)
def group_var_float32(ndarray[float32_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float32_t, ndim=2] values,
              ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, ct
        ndarray[float32_t, ndim=2] nobs, sumx, sumxx

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):

            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
                    sumxx[lab, j] += val * val
    else:
        for i in range(N):

            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]
            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val
                sumxx[lab, 0] += val * val


    for i in range(len(counts)):
        for j in range(K):
            ct = nobs[i, j]
            if ct < 2:
                out[i, j] = nan
            else:
                out[i, j] = ((ct * sumxx[i, j] - sumx[i, j] * sumx[i, j]) /
                             (ct * ct - ct))

@cython.wraparound(False)
@cython.boundscheck(False)
def group_var_bin_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):

    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, ct
        ndarray[float64_t, ndim=2] nobs, sumx, sumxx

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
                    sumxx[b, j] += val * val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val
                sumxx[b, 0] += val * val

    for i in range(ngroups):
        for j in range(K):
            ct = nobs[i, j]
            if ct < 2:
                out[i, j] = nan
            else:
                out[i, j] = ((ct * sumxx[i, j] - sumx[i, j] * sumx[i, j]) /
                             (ct * ct - ct))
@cython.wraparound(False)
@cython.boundscheck(False)
def group_var_bin_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] bins):

    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, ct
        ndarray[float32_t, ndim=2] nobs, sumx, sumxx

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
                    sumxx[b, j] += val * val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val
                sumxx[b, 0] += val * val

    for i in range(ngroups):
        for j in range(K):
            ct = nobs[i, j]
            if ct < 2:
                out[i, j] = nan
            else:
                out[i, j] = ((ct * sumxx[i, j] - sumx[i, j] * sumx[i, j]) /
                             (ct * ct - ct))

@cython.wraparound(False)
@cython.boundscheck(False)
def group_mean_float64(ndarray[float64_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]
                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]
            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val

    for i in range(len(counts)):
        for j in range(K):
            count = nobs[i, j]
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j] / count
@cython.wraparound(False)
@cython.boundscheck(False)
def group_mean_float32(ndarray[float32_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float32_t, ndim=2] values,
               ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]
                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]
            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val

    for i in range(len(counts)):
        for j in range(K):
            count = nobs[i, j]
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j] / count


def group_mean_bin_float64(ndarray[float64_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int64_t] bins):
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape
    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val

    for i in range(ngroups):
        for j in range(K):
            count = nobs[i, j]
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j] / count

def group_mean_bin_float32(ndarray[float32_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float32_t, ndim=2] values,
                   ndarray[int64_t] bins):
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        ndarray[float32_t, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape
    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val

    for i in range(ngroups):
        for j in range(K):
            count = nobs[i, j]
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j] / count

@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_float64(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] minx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    if val < minx[lab, j]:
                        minx[lab, j] = val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                if val < minx[lab, 0]:
                    minx[lab, 0] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = minx[i, j]
@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_float32(ndarray[float32_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float32_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] minx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    if val < minx[lab, j]:
                        minx[lab, j] = val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                if val < minx[lab, 0]:
                    minx[lab, 0] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = minx[i, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_bin_float64(ndarray[float64_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] minx, nobs

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    if val < minx[b, j]:
                        minx[b, j] = val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                if val < minx[b, 0]:
                    minx[b, 0] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = minx[i, j]
@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_bin_float32(ndarray[float32_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float32_t, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        ndarray[float32_t, ndim=2] minx, nobs

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    if val < minx[b, j]:
                        minx[b, j] = val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                if val < minx[b, 0]:
                    minx[b, 0] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = minx[i, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_float64(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] maxx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    if val > maxx[lab, j]:
                        maxx[lab, j] = val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                if val > maxx[lab, 0]:
                    maxx[lab, 0] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = maxx[i, j]
@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_float32(ndarray[float32_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float32_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        ndarray[float32_t, ndim=2] maxx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    if val > maxx[lab, j]:
                        maxx[lab, j] = val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                if val > maxx[lab, 0]:
                    maxx[lab, 0] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = maxx[i, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_bin_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] maxx, nobs

    nobs = np.zeros_like(out)
    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    if val > maxx[b, j]:
                        maxx[b, j] = val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                if val > maxx[b, 0]:
                    maxx[b, 0] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = maxx[i, j]
@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_bin_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        ndarray[float32_t, ndim=2] maxx, nobs

    nobs = np.zeros_like(out)
    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    if val > maxx[b, j]:
                        maxx[b, j] = val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                if val > maxx[b, 0]:
                    maxx[b, 0] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = maxx[i, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def group_ohlc_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        float64_t vopen, vhigh, vlow, vclose, NA
        bint got_first = 0

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    if out.shape[1] != 4:
        raise ValueError('Output array must have 4 columns')

    NA = np.nan

    b = 0
    if K > 1:
        raise NotImplementedError
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                if not got_first:
                    out[b, 0] = NA
                    out[b, 1] = NA
                    out[b, 2] = NA
                    out[b, 3] = NA
                else:
                    out[b, 0] = vopen
                    out[b, 1] = vhigh
                    out[b, 2] = vlow
                    out[b, 3] = vclose
                b += 1
                got_first = 0

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                if not got_first:
                    got_first = 1
                    vopen = val
                    vlow = val
                    vhigh = val
                else:
                    if val < vlow:
                        vlow = val
                    if val > vhigh:
                        vhigh = val
                vclose = val

        if not got_first:
            out[b, 0] = NA
            out[b, 1] = NA
            out[b, 2] = NA
            out[b, 3] = NA
        else:
            out[b, 0] = vopen
            out[b, 1] = vhigh
            out[b, 2] = vlow
            out[b, 3] = vclose
@cython.wraparound(False)
@cython.boundscheck(False)
def group_ohlc_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float32_t val, count
        float32_t vopen, vhigh, vlow, vclose, NA
        bint got_first = 0

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    if out.shape[1] != 4:
        raise ValueError('Output array must have 4 columns')

    NA = np.nan

    b = 0
    if K > 1:
        raise NotImplementedError
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                if not got_first:
                    out[b, 0] = NA
                    out[b, 1] = NA
                    out[b, 2] = NA
                    out[b, 3] = NA
                else:
                    out[b, 0] = vopen
                    out[b, 1] = vhigh
                    out[b, 2] = vlow
                    out[b, 3] = vclose
                b += 1
                got_first = 0

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                if not got_first:
                    got_first = 1
                    vopen = val
                    vlow = val
                    vhigh = val
                else:
                    if val < vlow:
                        vlow = val
                    if val > vhigh:
                        vhigh = val
                vclose = val

        if not got_first:
            out[b, 0] = NA
            out[b, 1] = NA
            out[b, 2] = NA
            out[b, 3] = NA
        else:
            out[b, 0] = vopen
            out[b, 1] = vhigh
            out[b, 2] = vlow
            out[b, 3] = vclose

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
def left_join_indexer_unique_float32(ndarray[float32_t] left,
                                      ndarray[float32_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        float32_t lval, rval

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


def left_join_indexer_float32(ndarray[float32_t] left,
                              ndarray[float32_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        float32_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[float32_t] result

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
    result = np.empty(count, dtype=np.float32)

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
def outer_join_indexer_float32(ndarray[float32_t] left,
                                ndarray[float32_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        float32_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[float32_t] result

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
    result = np.empty(count, dtype=np.float32)

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
def inner_join_indexer_float32(ndarray[float32_t] left,
                              ndarray[float32_t] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        float32_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[float32_t] result

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
    result = np.empty(count, dtype=np.float32)

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


