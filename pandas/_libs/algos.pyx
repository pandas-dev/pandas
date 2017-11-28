# cython: profile=False

cimport numpy as np
import numpy as np

cimport cython
from cython cimport Py_ssize_t

np.import_array()

cdef float64_t FP_ERR = 1e-13

cimport util

from libc.stdlib cimport malloc, free
from libc.string cimport memmove

from numpy cimport (ndarray,
                    NPY_INT64, NPY_UINT64, NPY_INT32, NPY_INT16, NPY_INT8,
                    NPY_FLOAT32, NPY_FLOAT64,
                    NPY_OBJECT,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t,
                    double_t)


cdef double NaN = <double> np.NaN
cdef double nan = NaN

from libc.math cimport sqrt, fabs

# this is our util.pxd
from util cimport numeric, get_nat

import missing

cdef int64_t iNaT = get_nat()

cdef:
    int TIEBREAK_AVERAGE = 0
    int TIEBREAK_MIN = 1
    int TIEBREAK_MAX = 2
    int TIEBREAK_FIRST = 3
    int TIEBREAK_FIRST_DESCENDING = 4
    int TIEBREAK_DENSE = 5

tiebreakers = {
    'average': TIEBREAK_AVERAGE,
    'min': TIEBREAK_MIN,
    'max': TIEBREAK_MAX,
    'first': TIEBREAK_FIRST,
    'dense': TIEBREAK_DENSE,
}


cdef inline are_diff(object left, object right):
    try:
        return fabs(left - right) > FP_ERR
    except TypeError:
        return left != right


class Infinity(object):
    """ provide a positive Infinity comparision method for ranking """

    __lt__ = lambda self, other: False
    __le__ = lambda self, other: self is other
    __eq__ = lambda self, other: self is other
    __ne__ = lambda self, other: self is not other
    __gt__ = lambda self, other: self is not other
    __ge__ = lambda self, other: True


class NegInfinity(object):
    """ provide a negative Infinity comparision method for ranking """

    __lt__ = lambda self, other: self is not other
    __le__ = lambda self, other: True
    __eq__ = lambda self, other: self is other
    __ne__ = lambda self, other: self is not other
    __gt__ = lambda self, other: False
    __ge__ = lambda self, other: self is other


@cython.wraparound(False)
@cython.boundscheck(False)
def is_lexsorted(list list_of_arrays):
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
        vecs[i] = <int64_t*> arr.data

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
def groupsort_indexer(ndarray[int64_t] index, Py_ssize_t ngroups):
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


cpdef numeric median(numeric[:] arr):
    """
    A faster median
    """
    cdef Py_ssize_t n = arr.size

    if n == 0:
        return np.NaN

    arr = arr.copy()

    if n % 2:
        return kth_smallest(arr, n // 2)
    else:
        return (kth_smallest(arr, n // 2) +
                kth_smallest(arr, n // 2 - 1)) / 2


# ----------------------------------------------------------------------
# Pairwise correlation/covariance


@cython.boundscheck(False)
@cython.wraparound(False)
def nancorr(ndarray[float64_t, ndim=2] mat, bint cov=0, minp=None):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        bint minpv
        ndarray[float64_t, ndim=2] result
        ndarray[uint8_t, ndim=2] mask
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
        ndarray[uint8_t, ndim=2] mask
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


# generated from template
include "algos_common_helper.pxi"
include "algos_rank_helper.pxi"
include "algos_take_helper.pxi"
