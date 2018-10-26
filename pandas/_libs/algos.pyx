# -*- coding: utf-8 -*-

import cython
from cython import Py_ssize_t

from libc.stdlib cimport malloc, free
from libc.string cimport memmove
from libc.math cimport fabs, sqrt

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray,
                    NPY_INT64, NPY_INT32, NPY_INT16, NPY_INT8,
                    NPY_UINT64, NPY_UINT32, NPY_UINT16, NPY_UINT8,
                    NPY_FLOAT32, NPY_FLOAT64,
                    NPY_OBJECT,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t,
                    double_t)
cnp.import_array()


cimport util
from util cimport numeric, get_nat

from khash cimport (khiter_t,
                    kh_destroy_int64, kh_put_int64,
                    kh_init_int64, kh_int64_t,
                    kh_resize_int64, kh_get_int64)

import missing

cdef float64_t FP_ERR = 1e-13

cdef double NaN = <double> np.NaN
cdef double nan = NaN

cdef int64_t iNaT = get_nat()

tiebreakers = {
    'average': TIEBREAK_AVERAGE,
    'min': TIEBREAK_MIN,
    'max': TIEBREAK_MAX,
    'first': TIEBREAK_FIRST,
    'dense': TIEBREAK_DENSE,
}


cdef inline bint are_diff(object left, object right):
    try:
        return fabs(left - right) > FP_ERR
    except TypeError:
        return left != right


class Infinity(object):
    """ provide a positive Infinity comparison method for ranking """

    __lt__ = lambda self, other: False
    __le__ = lambda self, other: isinstance(other, Infinity)
    __eq__ = lambda self, other: isinstance(other, Infinity)
    __ne__ = lambda self, other: not isinstance(other, Infinity)
    __gt__ = lambda self, other: (not isinstance(other, Infinity) and
                                  not missing.checknull(other))
    __ge__ = lambda self, other: not missing.checknull(other)


class NegInfinity(object):
    """ provide a negative Infinity comparison method for ranking """

    __lt__ = lambda self, other: (not isinstance(other, NegInfinity) and
                                  not missing.checknull(other))
    __le__ = lambda self, other: not missing.checknull(other)
    __eq__ = lambda self, other: isinstance(other, NegInfinity)
    __ne__ = lambda self, other: not isinstance(other, NegInfinity)
    __gt__ = lambda self, other: False
    __ge__ = lambda self, other: isinstance(other, NegInfinity)


cpdef ndarray[int64_t, ndim=1] unique_deltas(int64_t[:] arr):
    """
    Efficiently find the unique first-differences of the given array.

    Parameters
    ----------
    arr : ndarray[in64_t]

    Returns
    -------
    result : ndarray[int64_t]
        result is sorted
    """
    cdef:
        Py_ssize_t i, n = len(arr)
        int64_t val
        khiter_t k
        kh_int64_t *table
        int ret = 0
        list uniques = []

    table = kh_init_int64()
    kh_resize_int64(table, 10)
    for i in range(n - 1):
        val = arr[i + 1] - arr[i]
        k = kh_get_int64(table, val)
        if k == table.n_buckets:
            kh_put_int64(table, val, &ret)
            uniques.append(val)
    kh_destroy_int64(table)

    result = np.array(uniques, dtype=np.int64)
    result.sort()
    return result


@cython.wraparound(False)
@cython.boundscheck(False)
def is_lexsorted(list_of_arrays: list) -> bint:
    cdef:
        Py_ssize_t i
        Py_ssize_t n, nlevels
        int64_t k, cur, pre
        ndarray arr
        bint result = True

    nlevels = len(list_of_arrays)
    n = len(list_of_arrays[0])

    cdef int64_t **vecs = <int64_t**> malloc(nlevels * sizeof(int64_t*))
    for i in range(nlevels):
        arr = list_of_arrays[i]
        assert arr.dtype.name == 'int64'
        vecs[i] = <int64_t*> cnp.PyArray_DATA(arr)

    # Assume uniqueness??
    with nogil:
        for i in range(1, n):
            for k in range(nlevels):
                cur = vecs[k][i]
                pre = vecs[k][i -1]
                if cur == pre:
                    continue
                elif cur > pre:
                    break
                else:
                    result = False
                    break
    free(vecs)
    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def groupsort_indexer(int64_t[:] index, Py_ssize_t ngroups):
    """
    compute a 1-d indexer that is an ordering of the passed index,
    ordered by the groups. This is a reverse of the label
    factorization process.

    Parameters
    ----------
    index: int64 ndarray
        mappings from group -> position
    ngroups: int64
        number of groups

    return a tuple of (1-d indexer ordered by groups, group counts)
    """

    cdef:
        Py_ssize_t i, loc, label, n
        ndarray[int64_t] counts, where, result

    counts = np.zeros(ngroups + 1, dtype=np.int64)
    n = len(index)
    result = np.zeros(n, dtype=np.int64)
    where = np.zeros(ngroups + 1, dtype=np.int64)

    with nogil:

        # count group sizes, location 0 for NA
        for i in range(n):
            counts[index[i] + 1] += 1

        # mark the start of each contiguous group of like-indexed data
        for i in range(1, ngroups + 1):
            where[i] = where[i - 1] + counts[i - 1]

        # this is our indexer
        for i in range(n):
            label = index[i] + 1
            result[where[label]] = i
            where[label] += 1

    return result, counts


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef numeric kth_smallest(numeric[:] a, Py_ssize_t k) nogil:
    cdef:
        Py_ssize_t i, j, l, m, n = a.shape[0]
        numeric x

    with nogil:
        l = 0
        m = n - 1

        while l < m:
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


# ----------------------------------------------------------------------
# Pairwise correlation/covariance


@cython.boundscheck(False)
@cython.wraparound(False)
def nancorr(ndarray[float64_t, ndim=2] mat, bint cov=0, minp=None):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        bint minpv
        ndarray[float64_t, ndim=2] result
        uint8_t[:, :] mask
        int64_t nobs = 0
        float64_t vx, vy, sumx, sumy, sumxx, sumyy, meanx, meany, divisor

    N, K = (<object> mat).shape

    if minp is None:
        minpv = 1
    else:
        minpv = <int>minp

    result = np.empty((K, K), dtype=np.float64)
    mask = np.isfinite(mat).view(np.uint8)

    with nogil:
        for xi in range(K):
            for yi in range(xi + 1):
                nobs = sumxx = sumyy = sumx = sumy = 0
                for i in range(N):
                    if mask[i, xi] and mask[i, yi]:
                        vx = mat[i, xi]
                        vy = mat[i, yi]
                        nobs += 1
                        sumx += vx
                        sumy += vy

                if nobs < minpv:
                    result[xi, yi] = result[yi, xi] = NaN
                else:
                    meanx = sumx / nobs
                    meany = sumy / nobs

                    # now the cov numerator
                    sumx = 0

                    for i in range(N):
                        if mask[i, xi] and mask[i, yi]:
                            vx = mat[i, xi] - meanx
                            vy = mat[i, yi] - meany

                            sumx += vx * vy
                            sumxx += vx * vx
                            sumyy += vy * vy

                    divisor = (nobs - 1.0) if cov else sqrt(sumxx * sumyy)

                    if divisor != 0:
                        result[xi, yi] = result[yi, xi] = sumx / divisor
                    else:
                        result[xi, yi] = result[yi, xi] = NaN

    return result

# ----------------------------------------------------------------------
# Pairwise Spearman correlation


@cython.boundscheck(False)
@cython.wraparound(False)
def nancorr_spearman(ndarray[float64_t, ndim=2] mat, Py_ssize_t minp=1):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        ndarray[float64_t, ndim=2] result
        ndarray[float64_t, ndim=1] maskedx
        ndarray[float64_t, ndim=1] maskedy
        uint8_t[:, :] mask
        int64_t nobs = 0
        float64_t vx, vy, sumx, sumxx, sumyy, mean, divisor

    N, K = (<object> mat).shape

    result = np.empty((K, K), dtype=np.float64)
    mask = np.isfinite(mat).view(np.uint8)

    for xi in range(K):
        for yi in range(xi + 1):
            nobs = 0
            for i in range(N):
                if mask[i, xi] and mask[i, yi]:
                    nobs += 1

            if nobs < minp:
                result[xi, yi] = result[yi, xi] = NaN
            else:
                maskedx = np.empty(nobs, dtype=np.float64)
                maskedy = np.empty(nobs, dtype=np.float64)
                j = 0
                for i in range(N):
                    if mask[i, xi] and mask[i, yi]:
                        maskedx[j] = mat[i, xi]
                        maskedy[j] = mat[i, yi]
                        j += 1
                maskedx = rank_1d_float64(maskedx)
                maskedy = rank_1d_float64(maskedy)

                mean = (nobs + 1) / 2.

                # now the cov numerator
                sumx = sumxx = sumyy = 0

                for i in range(nobs):
                    vx = maskedx[i] - mean
                    vy = maskedy[i] - mean

                    sumx += vx * vy
                    sumxx += vx * vx
                    sumyy += vy * vy

                divisor = sqrt(sumxx * sumyy)

                if divisor != 0:
                    result[xi, yi] = result[yi, xi] = sumx / divisor
                else:
                    result[xi, yi] = result[yi, xi] = NaN

    return result


# ----------------------------------------------------------------------

ctypedef fused algos_t:
    float64_t
    float32_t
    object
    int64_t
    int32_t
    int16_t
    int8_t
    uint64_t
    uint32_t
    uint16_t
    uint8_t


# TODO: unused; needed?
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices(algos_t[:] index):
    """
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    """
    cdef:
        Py_ssize_t i, length
        dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def pad(algos_t[:] old, algos_t[:] new, limit=None):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t, ndim=1] indexer
        algos_t cur, next
        int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
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


pad_float64 = pad["float64_t"]
pad_float32 = pad["float32_t"]
pad_object = pad["object"]
pad_int64 = pad["int64_t"]
pad_int32 = pad["int32_t"]
pad_int16 = pad["int16_t"]
pad_int8 = pad["int8_t"]
pad_uint64 = pad["uint64_t"]
pad_uint32 = pad["uint32_t"]
pad_uint16 = pad["uint16_t"]
pad_uint8 = pad["uint8_t"]
pad_bool = pad["uint8_t"]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace(algos_t[:] values, uint8_t[:] mask, limit=None):
    cdef:
        Py_ssize_t i, N
        algos_t val
        int lim, fill_count = 0

    N = len(values)

    # GH#2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
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


pad_inplace_float64 = pad_inplace["float64_t"]
pad_inplace_float32 = pad_inplace["float32_t"]
pad_inplace_object = pad_inplace["object"]
pad_inplace_int64 = pad_inplace["int64_t"]
pad_inplace_int32 = pad_inplace["int32_t"]
pad_inplace_uint64 = pad_inplace["uint64_t"]
pad_inplace_bool = pad_inplace["uint8_t"]


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace(algos_t[:, :] values, uint8_t[:, :] mask, limit=None):
    cdef:
        Py_ssize_t i, j, N, K
        algos_t val
        int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH#2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
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


pad_2d_inplace_float64 = pad_2d_inplace["float64_t"]
pad_2d_inplace_float32 = pad_2d_inplace["float32_t"]
pad_2d_inplace_object = pad_2d_inplace["object"]
pad_2d_inplace_int64 = pad_2d_inplace["int64_t"]
pad_2d_inplace_int32 = pad_2d_inplace["int32_t"]
pad_2d_inplace_uint64 = pad_2d_inplace["uint64_t"]
pad_2d_inplace_bool = pad_2d_inplace["uint8_t"]


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
def backfill(algos_t[:] old, algos_t[:] new, limit=None):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t, ndim=1] indexer
        algos_t cur, prev
        int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
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


backfill_float64 = backfill["float64_t"]
backfill_float32 = backfill["float32_t"]
backfill_object = backfill["object"]
backfill_int64 = backfill["int64_t"]
backfill_int32 = backfill["int32_t"]
backfill_int16 = backfill["int16_t"]
backfill_int8 = backfill["int8_t"]
backfill_uint64 = backfill["uint64_t"]
backfill_uint32 = backfill["uint32_t"]
backfill_uint16 = backfill["uint16_t"]
backfill_uint8 = backfill["uint8_t"]
backfill_bool = backfill["uint8_t"]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace(algos_t[:] values, uint8_t[:] mask, limit=None):
    cdef:
        Py_ssize_t i, N
        algos_t val
        int lim, fill_count = 0

    N = len(values)

    # GH#2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
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


backfill_inplace_float64 = backfill_inplace["float64_t"]
backfill_inplace_float32 = backfill_inplace["float32_t"]
backfill_inplace_object = backfill_inplace["object"]
backfill_inplace_int64 = backfill_inplace["int64_t"]
backfill_inplace_int32 = backfill_inplace["int32_t"]
backfill_inplace_uint64 = backfill_inplace["uint64_t"]
backfill_inplace_bool = backfill_inplace["uint8_t"]


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace(algos_t[:, :] values, uint8_t[:, :] mask, limit=None):
    cdef:
        Py_ssize_t i, j, N, K
        algos_t val
        int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH#2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
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


backfill_2d_inplace_float64 = backfill_2d_inplace["float64_t"]
backfill_2d_inplace_float32 = backfill_2d_inplace["float32_t"]
backfill_2d_inplace_object = backfill_2d_inplace["object"]
backfill_2d_inplace_int64 = backfill_2d_inplace["int64_t"]
backfill_2d_inplace_int32 = backfill_2d_inplace["int32_t"]
backfill_2d_inplace_uint64 = backfill_2d_inplace["uint64_t"]
backfill_2d_inplace_bool = backfill_2d_inplace["uint8_t"]


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap(algos_t[:] index, object func):
    cdef:
        Py_ssize_t length = index.shape[0]
        Py_ssize_t i = 0
        object[:] result = np.empty(length, dtype=np.object_)

    from pandas._libs.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)


arrmap_float64 = arrmap["float64_t"]
arrmap_float32 = arrmap["float32_t"]
arrmap_object = arrmap["object"]
arrmap_int64 = arrmap["int64_t"]
arrmap_int32 = arrmap["int32_t"]
arrmap_uint64 = arrmap["uint64_t"]
arrmap_bool = arrmap["uint8_t"]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic(ndarray[algos_t, ndim=1] arr, bint timelike):
    """
    Returns
    -------
    is_monotonic_inc, is_monotonic_dec, is_unique
    """
    cdef:
        Py_ssize_t i, n
        algos_t prev, cur
        bint is_monotonic_inc = 1
        bint is_monotonic_dec = 1
        bint is_unique = 1
        bint is_strict_monotonic = 1

    n = len(arr)

    if n == 1:
        if arr[0] != arr[0] or (timelike and <int64_t>arr[0] == iNaT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and <int64_t>arr[0] == iNaT:
        return False, False, True

    if algos_t is not object:
        with nogil:
            prev = arr[0]
            for i in range(1, n):
                cur = arr[i]
                if timelike and <int64_t>cur == iNaT:
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
    else:
        # object-dtype, identical to above except we cannot use `with nogil`
        prev = arr[0]
        for i in range(1, n):
            cur = arr[i]
            if timelike and <int64_t>cur == iNaT:
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

    is_strict_monotonic = is_unique and (is_monotonic_inc or is_monotonic_dec)
    return is_monotonic_inc, is_monotonic_dec, is_strict_monotonic


is_monotonic_float64 = is_monotonic["float64_t"]
is_monotonic_float32 = is_monotonic["float32_t"]
is_monotonic_object = is_monotonic["object"]
is_monotonic_int64 = is_monotonic["int64_t"]
is_monotonic_int32 = is_monotonic["int32_t"]
is_monotonic_int16 = is_monotonic["int16_t"]
is_monotonic_int8 = is_monotonic["int8_t"]
is_monotonic_uint64 = is_monotonic["uint64_t"]
is_monotonic_uint32 = is_monotonic["uint32_t"]
is_monotonic_uint16 = is_monotonic["uint16_t"]
is_monotonic_uint8 = is_monotonic["uint8_t"]
is_monotonic_bool = is_monotonic["uint8_t"]


# generated from template
include "algos_common_helper.pxi"
include "algos_rank_helper.pxi"
include "algos_take_helper.pxi"
