# cython: profile=False

from numpy cimport *
cimport numpy as np
import numpy as np

cimport cython

import_array()

cdef float64_t FP_ERR = 1e-13

cimport util

from libc.stdlib cimport malloc, free
from libc.string cimport memmove

from numpy cimport NPY_INT8 as NPY_int8
from numpy cimport NPY_INT16 as NPY_int16
from numpy cimport NPY_INT32 as NPY_int32
from numpy cimport NPY_INT64 as NPY_int64
from numpy cimport NPY_FLOAT16 as NPY_float16
from numpy cimport NPY_FLOAT32 as NPY_float32
from numpy cimport NPY_FLOAT64 as NPY_float64

from numpy cimport (int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float16_t, float32_t, float64_t)

int8 = np.dtype(np.int8)
int16 = np.dtype(np.int16)
int32 = np.dtype(np.int32)
int64 = np.dtype(np.int64)
float16 = np.dtype(np.float16)
float32 = np.dtype(np.float32)
float64 = np.dtype(np.float64)

cdef double NaN = <double> np.NaN
cdef double nan = NaN

cdef extern from "src/headers/math.h":
    double sqrt(double x) nogil
    double fabs(double) nogil

# this is our util.pxd
from util cimport numeric, get_nat

cimport lib
from lib cimport is_null_datetimelike
from pandas import lib

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


cdef inline Py_ssize_t swap(numeric *a, numeric *b) nogil except -1:
    cdef numeric t

    # cython doesn't allow pointer dereference so use array syntax
    t = a[0]
    a[0] = b[0]
    b[0] = t
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef numeric kth_smallest(numeric[:] a, Py_ssize_t k):
    cdef:
        Py_ssize_t i, j, l, m, n = a.size
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


cdef inline kth_smallest_c(float64_t* a, Py_ssize_t k, Py_ssize_t n):
    cdef:
        Py_ssize_t i, j, l, m
        double_t x, t

    l = 0
    m = n -1
    while (l<m):
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


# -------------- Min, Max subsequence

def max_subseq(ndarray[double_t] arr):
    cdef:
        Py_ssize_t i=0, s=0, e=0, T, n
        double m, S

    n = len(arr)

    if len(arr) == 0:
        return (-1, -1, None)

    m = arr[0]
    S = m
    T = 0

    for i in range(1, n):
        # S = max { S + A[i], A[i] )
        if (S > 0):
            S = S + arr[i]
        else:
            S = arr[i]
            T = i
        if S > m:
            s = T
            e = i
            m = S

    return (s, e, m)


def min_subseq(ndarray[double_t] arr):
    cdef:
        Py_ssize_t s, e
        double m

    (s, e, m) = max_subseq(-arr)

    return (s, e, -m)

#----------------------------------------------------------------------
# Pairwise correlation/covariance


@cython.boundscheck(False)
@cython.wraparound(False)
def nancorr(ndarray[float64_t, ndim=2] mat, cov=False, minp=None):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        ndarray[float64_t, ndim=2] result
        ndarray[uint8_t, ndim=2] mask
        int64_t nobs = 0
        float64_t vx, vy, sumx, sumy, sumxx, sumyy, meanx, meany, divisor

    N, K = (<object> mat).shape

    if minp is None:
        minp = 1

    result = np.empty((K, K), dtype=np.float64)
    mask = np.isfinite(mat).view(np.uint8)

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

            if nobs < minp:
                result[xi, yi] = result[yi, xi] = np.NaN
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
                    result[xi, yi] = result[yi, xi] = np.NaN

    return result

#----------------------------------------------------------------------
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
                result[xi, yi] = result[yi, xi] = np.NaN
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
                    result[xi, yi] = result[yi, xi] = np.NaN

    return result

#----------------------------------------------------------------------
# group operations


@cython.wraparound(False)
@cython.boundscheck(False)
def is_lexsorted(list list_of_arrays):
    cdef:
        int i
        Py_ssize_t n, nlevels
        int64_t k, cur, pre
        ndarray arr

    nlevels = len(list_of_arrays)
    n = len(list_of_arrays[0])

    cdef int64_t **vecs = <int64_t**> malloc(nlevels * sizeof(int64_t*))
    for i from 0 <= i < nlevels:
        arr = list_of_arrays[i]
        vecs[i] = <int64_t*> arr.data

    # Assume uniqueness??
    for i from 1 <= i < n:
        for k from 0 <= k < nlevels:
            cur = vecs[k][i]
            pre = vecs[k][i -1]
            if cur == pre:
                continue
            elif cur > pre:
                break
            else:
                return False
    free(vecs)
    return True


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
        for i from 0 <= i < n:
            counts[index[i] + 1] += 1

        # mark the start of each contiguous group of like-indexed data
        for i from 1 <= i < ngroups + 1:
            where[i] = where[i - 1] + counts[i - 1]

        # this is our indexer
        for i from 0 <= i < n:
            label = index[i] + 1
            result[where[label]] = i
            where[label] += 1

    return result, counts

# TODO: aggregate multiple columns in single pass
#----------------------------------------------------------------------
# first, nth, last


@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_object(ndarray[object, ndim=2] out,
                     ndarray[int64_t] counts,
                     ndarray[object, ndim=2] values,
                     ndarray[int64_t] labels,
                     int64_t rank):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab
        object val
        float64_t count
        ndarray[int64_t, ndim=2] nobs
        ndarray[object, ndim=2] resx

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty((<object> out).shape, dtype=object)

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
                out[i, j] = <object> nan
            else:
                out[i, j] = resx[i, j]


@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_bin_object(ndarray[object, ndim=2] out,
                         ndarray[int64_t] counts,
                         ndarray[object, ndim=2] values,
                         ndarray[int64_t] bins, int64_t rank):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[float64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.float64)
    resx = np.empty((<object> out).shape, dtype=object)

    if len(bins) == 0:
        return
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
def group_last_object(ndarray[object, ndim=2] out,
                      ndarray[int64_t] counts,
                      ndarray[object, ndim=2] values,
                      ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty((<object> out).shape, dtype=object)

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_last_bin_object(ndarray[object, ndim=2] out,
                          ndarray[int64_t] counts,
                          ndarray[object, ndim=2] values,
                          ndarray[int64_t] bins):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[float64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.float64)
    resx = np.empty((<object> out).shape, dtype=object)

    if len(bins) == 0:
        return
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

cdef inline float64_t _median_linear(float64_t* a, int n):
    cdef int i, j, na_count = 0
    cdef float64_t result
    cdef float64_t* tmp

    if n == 0:
        return NaN

    # count NAs
    for i in range(n):
        if a[i] != a[i]:
            na_count += 1

    if na_count:
        if na_count == n:
            return NaN

        tmp = <float64_t*> malloc((n - na_count) * sizeof(float64_t))

        j = 0
        for i in range(n):
            if a[i] == a[i]:
                tmp[j] = a[i]
                j += 1

        a = tmp
        n -= na_count

    if n % 2:
        result = kth_smallest_c( a, n / 2, n)
    else:
        result = (kth_smallest_c(a, n / 2, n) +
                  kth_smallest_c(a, n / 2 - 1, n)) / 2

    if na_count:
        free(a)

    return result


# generated from template
include "algos_common_helper.pxi"
include "algos_groupby_helper.pxi"
include "algos_rank_helper.pxi"
include "algos_take_helper.pxi"
