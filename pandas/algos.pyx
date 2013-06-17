from numpy cimport *
cimport numpy as np
import numpy as np

cimport cython

import_array()

cdef float64_t FP_ERR = 1e-13

cimport util

from libc.stdlib cimport malloc, free

from numpy cimport NPY_INT8 as NPY_int8
from numpy cimport NPY_INT16 as NPY_int16
from numpy cimport NPY_INT32 as NPY_int32
from numpy cimport NPY_INT64 as NPY_int64
from numpy cimport NPY_FLOAT16 as NPY_float16
from numpy cimport NPY_FLOAT32 as NPY_float32
from numpy cimport NPY_FLOAT64 as NPY_float64

int8 = np.dtype(np.int8)
int16 = np.dtype(np.int16)
int32 = np.dtype(np.int32)
int64 = np.dtype(np.int64)
float16 = np.dtype(np.float16)
float32 = np.dtype(np.float32)
float64 = np.dtype(np.float64)

cdef np.int8_t MINint8 = np.iinfo(np.int8).min
cdef np.int16_t MINint16 = np.iinfo(np.int16).min
cdef np.int32_t MINint32 = np.iinfo(np.int32).min
cdef np.int64_t MINint64 = np.iinfo(np.int64).min
cdef np.float16_t MINfloat16 = np.NINF
cdef np.float32_t MINfloat32 = np.NINF
cdef np.float64_t MINfloat64 = np.NINF

cdef np.int8_t MAXint8 = np.iinfo(np.int8).max
cdef np.int16_t MAXint16 = np.iinfo(np.int16).max
cdef np.int32_t MAXint32 = np.iinfo(np.int32).max
cdef np.int64_t MAXint64 = np.iinfo(np.int64).max
cdef np.float16_t MAXfloat16 = np.inf
cdef np.float32_t MAXfloat32 = np.inf
cdef np.float64_t MAXfloat64 = np.inf

cdef double NaN = <double> np.NaN
cdef double nan = NaN


cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b


cdef extern from "src/headers/math.h":
    double sqrt(double x)
    double fabs(double)
    int signbit(double)

from pandas import lib

include "skiplist.pyx"


cdef:
    int TIEBREAK_AVERAGE = 0
    int TIEBREAK_MIN = 1
    int TIEBREAK_MAX = 2
    int TIEBREAK_FIRST = 3
    int TIEBREAK_FIRST_DESCENDING = 4

tiebreakers = {
    'average' : TIEBREAK_AVERAGE,
    'min' : TIEBREAK_MIN,
    'max' : TIEBREAK_MAX,
    'first' : TIEBREAK_FIRST
}


# ctypedef fused pvalue_t:
#     float64_t
#     int64_t
#     object

# from cython cimport floating, integral

cdef _take_2d_float64(ndarray[float64_t, ndim=2] values,
                      object idx):
    cdef:
        Py_ssize_t i, j, N, K
        ndarray[Py_ssize_t, ndim=2, cast=True] indexer = idx
        ndarray[float64_t, ndim=2] result
        object val

    N, K = (<object> values).shape
    result = np.empty_like(values)
    for i in range(N):
        for j in range(K):
            result[i, j] = values[i, indexer[i, j]]
    return result

cdef _take_2d_int64(ndarray[int64_t, ndim=2] values,
                      object idx):
    cdef:
        Py_ssize_t i, j, N, K
        ndarray[Py_ssize_t, ndim=2, cast=True] indexer = idx
        ndarray[int64_t, ndim=2] result
        object val

    N, K = (<object> values).shape
    result = np.empty_like(values)
    for i in range(N):
        for j in range(K):
            result[i, j] = values[i, indexer[i, j]]
    return result

cdef _take_2d_object(ndarray[object, ndim=2] values,
                     object idx):
    cdef:
        Py_ssize_t i, j, N, K
        ndarray[Py_ssize_t, ndim=2, cast=True] indexer = idx
        ndarray[object, ndim=2] result
        object val

    N, K = (<object> values).shape
    result = values.copy()
    for i in range(N):
        for j in range(K):
            result[i, j] = values[i, indexer[i, j]]
    return result


def rank_1d_float64(object in_arr, ties_method='average', ascending=True,
                    na_option='keep'):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[float64_t] sorted_data, ranks, values
        ndarray[int64_t] argsorted
        float64_t val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
        bint keep_na = 0
    tiebreak = tiebreakers[ties_method]

    values = np.asarray(in_arr).copy()

    keep_na = na_option == 'keep'

    if ascending ^ (na_option == 'top'):
        nan_value = np.inf
    else:
        nan_value = -np.inf
    mask = np.isnan(values)
    np.putmask(values, mask, nan_value)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort()

    if not ascending:
        _as = _as[::-1]

    sorted_data = values.take(_as)
    argsorted = _as.astype('i8')

    for i in range(n):
        sum_ranks += i + 1
        dups += 1
        val = sorted_data[i]
        if (val == nan_value) and keep_na:
            ranks[argsorted[i]] = nan
            continue
        if i == n - 1 or fabs(sorted_data[i + 1] - val) > FP_ERR:
            if tiebreak == TIEBREAK_AVERAGE:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = sum_ranks / dups
            elif tiebreak == TIEBREAK_MIN:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i - dups + 2
            elif tiebreak == TIEBREAK_MAX:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i + 1
            elif tiebreak == TIEBREAK_FIRST:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = j + 1
            elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = 2 * i - j - dups + 2
            sum_ranks = dups = 0
    return ranks


def rank_1d_int64(object in_arr, ties_method='average', ascending=True,
                  na_option='keep'):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[int64_t] sorted_data, values
        ndarray[float64_t] ranks
        ndarray[int64_t] argsorted
        int64_t val
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    values = np.asarray(in_arr)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort()

    if not ascending:
        _as = _as[::-1]

    sorted_data = values.take(_as)
    argsorted = _as.astype('i8')

    for i in range(n):
        sum_ranks += i + 1
        dups += 1
        val = sorted_data[i]
        if i == n - 1 or fabs(sorted_data[i + 1] - val) > 0:
            if tiebreak == TIEBREAK_AVERAGE:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = sum_ranks / dups
            elif tiebreak == TIEBREAK_MIN:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i - dups + 2
            elif tiebreak == TIEBREAK_MAX:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i + 1
            elif tiebreak == TIEBREAK_FIRST:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = j + 1
            elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = 2 * i - j - dups + 2
            sum_ranks = dups = 0
    return ranks


def rank_2d_float64(object in_arr, axis=0, ties_method='average',
                    ascending=True, na_option='keep'):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, dups = 0
        ndarray[float64_t, ndim=2] ranks, values
        ndarray[int64_t, ndim=2] argsorted
        float64_t val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
        bint keep_na = 0

    tiebreak = tiebreakers[ties_method]

    keep_na = na_option == 'keep'

    in_arr = np.asarray(in_arr)

    if axis == 0:
        values = in_arr.T.copy()
    else:
        values = in_arr.copy()

    if ascending ^ (na_option == 'top'):
        nan_value = np.inf
    else:
        nan_value = -np.inf

    np.putmask(values, np.isnan(values), nan_value)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(axis=1, kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort(1)

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d_float64(values, _as)
    argsorted = _as.astype('i8')

    for i in range(n):
        dups = sum_ranks = 0
        for j in range(k):
            sum_ranks += j + 1
            dups += 1
            val = values[i, j]
            if val == nan_value and keep_na:
                ranks[i, argsorted[i, j]] = nan
                continue
            if j == k - 1 or fabs(values[i, j + 1] - val) > FP_ERR:
                if tiebreak == TIEBREAK_AVERAGE:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = sum_ranks / dups
                elif tiebreak == TIEBREAK_MIN:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j + 1
                elif tiebreak == TIEBREAK_FIRST:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = z + 1
                elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = 2 * j - z - dups + 2
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks


def rank_2d_int64(object in_arr, axis=0, ties_method='average',
                    ascending=True, na_option='keep'):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, dups = 0
        ndarray[float64_t, ndim=2] ranks
        ndarray[int64_t, ndim=2] argsorted
        ndarray[int64_t, ndim=2, cast=True] values
        int64_t val
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    if axis == 0:
        values = np.asarray(in_arr).T
    else:
        values = np.asarray(in_arr)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(axis=1, kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort(1)

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d_int64(values, _as)
    argsorted = _as.astype('i8')

    for i in range(n):
        dups = sum_ranks = 0
        for j in range(k):
            sum_ranks += j + 1
            dups += 1
            val = values[i, j]
            if j == k - 1 or fabs(values[i, j + 1] - val) > FP_ERR:
                if tiebreak == TIEBREAK_AVERAGE:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = sum_ranks / dups
                elif tiebreak == TIEBREAK_MIN:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j + 1
                elif tiebreak == TIEBREAK_FIRST:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = z + 1
                elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = 2 * j - z - dups + 2
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks


def rank_1d_generic(object in_arr, bint retry=1, ties_method='average',
                    ascending=True, na_option='keep'):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[float64_t] ranks
        ndarray sorted_data, values
        ndarray[int64_t] argsorted
        object val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
        bint keep_na = 0

    tiebreak = tiebreakers[ties_method]

    keep_na = na_option == 'keep'

    values = np.array(in_arr, copy=True)

    if values.dtype != np.object_:
        values = values.astype('O')

    if ascending ^ (na_option == 'top'):
        # always greater than everything
        nan_value = Infinity()
    else:
        nan_value = NegInfinity()

    mask = lib.isnullobj(values)
    np.putmask(values, mask, nan_value)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    try:
        _as = values.argsort()
    except TypeError:
        if not retry:
            raise

        valid_locs = (-mask).nonzero()[0]
        ranks.put(valid_locs, rank_1d_generic(values.take(valid_locs), 0,
                                              ties_method=ties_method,
                                              ascending=ascending))
        np.putmask(ranks, mask, np.nan)
        return ranks

    if not ascending:
        _as = _as[::-1]

    sorted_data = values.take(_as)
    argsorted = _as.astype('i8')

    for i in range(n):
        sum_ranks += i + 1
        dups += 1
        val = util.get_value_at(sorted_data, i)
        if val is nan_value and keep_na:
            ranks[argsorted[i]] = nan
            continue
        if (i == n - 1 or
            are_diff(util.get_value_at(sorted_data, i + 1), val)):
            if tiebreak == TIEBREAK_AVERAGE:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = sum_ranks / dups
            elif tiebreak == TIEBREAK_MIN:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i - dups + 2
            elif tiebreak == TIEBREAK_MAX:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i + 1
            elif tiebreak == TIEBREAK_FIRST:
                raise ValueError('first not supported for non-numeric data')
            sum_ranks = dups = 0
    return ranks

cdef inline are_diff(object left, object right):
    try:
        return fabs(left - right) > FP_ERR
    except TypeError:
        return left != right

_return_false = lambda self, other: False
_return_true = lambda self, other: True

class Infinity(object):

    __lt__ = _return_false
    __le__ = _return_false
    __eq__ = _return_false
    __ne__ = _return_true
    __gt__ = _return_true
    __ge__ = _return_true
    __cmp__ = _return_false

class NegInfinity(object):

    __lt__ = _return_true
    __le__ = _return_true
    __eq__ = _return_false
    __ne__ = _return_true
    __gt__ = _return_false
    __ge__ = _return_false
    __cmp__ = _return_true

def rank_2d_generic(object in_arr, axis=0, ties_method='average',
                    ascending=True, na_option='keep'):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, infs, dups = 0
        ndarray[float64_t, ndim=2] ranks
        ndarray[object, ndim=2] values
        ndarray[int64_t, ndim=2] argsorted
        object val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
        bint keep_na = 0

    tiebreak = tiebreakers[ties_method]

    keep_na = na_option == 'keep'

    in_arr = np.asarray(in_arr)

    if axis == 0:
        values = in_arr.T.copy()
    else:
        values = in_arr.copy()

    if values.dtype != np.object_:
        values = values.astype('O')

    if ascending ^ (na_option == 'top'):
        # always greater than everything
        nan_value = Infinity()
    else:
        nan_value = NegInfinity()

    mask = lib.isnullobj2d(values)
    np.putmask(values, mask, nan_value)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    try:
        _as = values.argsort(1)
    except TypeError:
        values = in_arr
        for i in range(len(values)):
            ranks[i] = rank_1d_generic(in_arr[i],
                                       ties_method=ties_method,
                                       ascending=ascending)
        if axis == 0:
            return ranks.T
        else:
            return ranks

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d_object(values, _as)
    argsorted = _as.astype('i8')

    for i in range(n):
        dups = sum_ranks = infs = 0
        for j in range(k):
            val = values[i, j]
            if val is nan_value and keep_na:
                ranks[i, argsorted[i, j]] = nan
                infs += 1
                continue
            sum_ranks += (j - infs) + 1
            dups += 1
            if j == k - 1 or are_diff(values[i, j + 1], val):
                if tiebreak == TIEBREAK_AVERAGE:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = sum_ranks / dups
                elif tiebreak == TIEBREAK_MIN:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j + 1
                elif tiebreak == TIEBREAK_FIRST:
                    raise ValueError('first not supported for '
                                     'non-numeric data')
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks

# def _take_indexer_2d(ndarray[float64_t, ndim=2] values,
#                      ndarray[Py_ssize_t, ndim=2, cast=True] indexer):
#     cdef:
#         Py_ssize_t i, j, N, K
#         ndarray[float64_t, ndim=2] result

#     N, K = (<object> values).shape
#     result = np.empty_like(values)
#     for i in range(N):
#         for j in range(K):
#             result[i, j] = values[i, indexer[i, j]]
#     return result


# Cython implementations of rolling sum, mean, variance, skewness,
# other statistical moment functions
#
# Misc implementation notes
# -------------------------
#
# - In Cython x * x is faster than x ** 2 for C types, this should be
#   periodically revisited to see if it's still true.
#
# -

def _check_minp(win, minp, N):
    if minp > win:
        raise ValueError('min_periods (%d) must be <= window (%d)'
                        % (minp, win))
    elif minp > N:
        minp = N + 1
    elif minp == 0:
        minp = 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')
    return minp

# original C implementation by N. Devillard.
# This code in public domain.
# Function :   kth_smallest()
# In       :   array of elements, # of elements in the array, rank k
# Out      :   one element
# Job      :   find the kth smallest element in the array

#             Reference:

#               Author: Wirth, Niklaus
#                Title: Algorithms + data structures = programs
#            Publisher: Englewood Cliffs: Prentice-Hall, 1976
# Physical description: 366 p.
#               Series: Prentice-Hall Series in Automatic Computation

def kth_smallest(ndarray[double_t] a, Py_ssize_t k):
    cdef:
        Py_ssize_t i,j,l,m,n
        double_t x, t

    n = len(a)

    l = 0
    m = n-1
    while (l<m):
        x = a[k]
        i = l
        j = m

        while 1:
            while a[i] < x: i += 1
            while x < a[j]: j -= 1
            if i <= j:
                t = a[i]
                a[i] = a[j]
                a[j] = t
                i += 1; j -= 1

            if i > j: break

        if j < k: l = i
        if k < i: m = j
    return a[k]

cdef inline kth_smallest_c(float64_t* a, Py_ssize_t k, Py_ssize_t n):
    cdef:
        Py_ssize_t i,j,l,m
        double_t x, t

    l = 0
    m = n-1
    while (l<m):
        x = a[k]
        i = l
        j = m

        while 1:
            while a[i] < x: i += 1
            while x < a[j]: j -= 1
            if i <= j:
                t = a[i]
                a[i] = a[j]
                a[j] = t
                i += 1; j -= 1

            if i > j: break

        if j < k: l = i
        if k < i: m = j
    return a[k]


def median(ndarray arr):
    '''
    A faster median
    '''
    cdef int n = len(arr)

    if len(arr) == 0:
        return np.NaN

    arr = arr.copy()

    if n % 2:
        return kth_smallest(arr, n / 2)
    else:
        return (kth_smallest(arr, n / 2) +
                kth_smallest(arr, n / 2 - 1)) / 2


# -------------- Min, Max subsequence

def max_subseq(ndarray[double_t] arr):
    cdef:
        Py_ssize_t i=0,s=0,e=0,T,n
        double m, S

    n = len(arr)

    if len(arr) == 0:
        return (-1,-1,None)

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

#-------------------------------------------------------------------------------
# Rolling sum

def roll_sum(ndarray[double_t] input, int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            sum_x += val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if val == val:
            nobs += 1
            sum_x += val

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                sum_x -= prev
                nobs -= 1

        if nobs >= minp:
            output[i] = sum_x
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling mean

def roll_mean(ndarray[double_t] input,
               int win, int minp):
    cdef:
        double val, prev, result, sum_x = 0
        Py_ssize_t nobs = 0, i, neg_ct = 0
        Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)
    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            sum_x += val
            if signbit(val):
                neg_ct += 1

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if val == val:
            nobs += 1
            sum_x += val
            if signbit(val):
                neg_ct += 1

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                sum_x -= prev
                nobs -= 1
                if signbit(prev):
                    neg_ct -= 1

        if nobs >= minp:
            result = sum_x / nobs
            if neg_ct == 0 and result < 0:
                # all positive
                output[i] = 0
            elif neg_ct == nobs and result > 0:
                # all negative
                output[i] = 0
            else:
                output[i] = result
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Exponentially weighted moving average

def ewma(ndarray[double_t] input, double_t com, int adjust):
    '''
    Compute exponentially-weighted moving average using center-of-mass.

    Parameters
    ----------
    input : ndarray (float64 type)
    com : float64

    Returns
    -------
    y : ndarray
    '''

    cdef double cur, prev, neww, oldw, adj
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    if N == 0:
        return output

    neww = 1. / (1. + com)
    oldw = 1. - neww
    adj = oldw

    if adjust:
        output[0] = neww * input[0]
    else:
        output[0] = input[0]

    for i from 1 <= i < N:
        cur = input[i]
        prev = output[i - 1]

        if cur == cur:
            if prev == prev:
                output[i] = oldw * prev + neww * cur
            else:
                output[i] = neww * cur
        else:
            output[i] = prev

    if adjust:
        for i from 0 <= i < N:
            cur = input[i]

            if cur == cur:
                output[i] = output[i] / (1. - adj)
                adj *= oldw
            else:
                if i >= 1:
                    output[i] = output[i - 1]

    return output

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
# Rolling variance

def roll_var(ndarray[double_t] input, int win, int minp, int ddof=1):
    cdef double val, prev, sum_x = 0, sum_xx = 0, nobs = 0
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            sum_x += val
            sum_xx += val * val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if val == val:
            nobs += 1
            sum_x += val
            sum_xx += val * val

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                sum_x -= prev
                sum_xx -= prev * prev
                nobs -= 1

        if nobs >= minp:
            # pathological case
            if nobs == 1:
                output[i] = 0
                continue

            val = (nobs * sum_xx - sum_x * sum_x) / (nobs * (nobs - ddof))
            if val < 0:
                val = 0

            output[i] = val
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling skewness

def roll_skew(ndarray[double_t] input, int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    # 3 components of the skewness equation
    cdef double A, B, C, R

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            x += val
            xx += val * val
            xxx += val * val * val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if val == val:
            nobs += 1
            x += val
            xx += val * val
            xxx += val * val * val

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                x -= prev
                xx -= prev * prev
                xxx -= prev * prev * prev

                nobs -= 1

        if nobs >= minp:
            A = x / nobs
            B = xx / nobs - A * A
            C = xxx / nobs - A * A * A - 3 * A * B

            R = sqrt(B)

            output[i] = ((sqrt(nobs * (nobs - 1.)) * C) /
                         ((nobs-2) * R * R * R))
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling kurtosis


def roll_kurt(ndarray[double_t] input,
               int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0, xxxx = 0
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    # 5 components of the kurtosis equation
    cdef double A, B, C, D, R, K

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1

            # seriously don't ask me why this is faster
            x += val
            xx += val * val
            xxx += val * val * val
            xxxx += val * val * val * val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if val == val:
            nobs += 1
            x += val
            xx += val * val
            xxx += val * val * val
            xxxx += val * val * val * val

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                x -= prev
                xx -= prev * prev
                xxx -= prev * prev * prev
                xxxx -= prev * prev * prev * prev

                nobs -= 1

        if nobs >= minp:
            A = x / nobs
            R = A * A
            B = xx / nobs - R
            R = R * A
            C = xxx / nobs - R - 3 * A * B
            R = R * A
            D = xxxx / nobs - R - 6*B*A*A - 4*C*A

            K = (nobs * nobs - 1.)*D/(B*B) - 3*((nobs-1.)**2)
            K = K / ((nobs - 2.)*(nobs-3.))

            output[i] = K
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling median, min, max

ctypedef double_t (* skiplist_f)(object sl, int n, int p)

cdef _roll_skiplist_op(ndarray arg, int win, int minp, skiplist_f op):
    cdef ndarray[double_t] input = arg
    cdef double val, prev, midpoint
    cdef IndexableSkiplist skiplist
    cdef Py_ssize_t nobs = 0, i

    cdef Py_ssize_t N = len(input)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    skiplist = IndexableSkiplist(win)

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist.insert(val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]

            if prev == prev:
                skiplist.remove(prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist.insert(val)

        output[i] = op(skiplist, nobs, minp)

    return output

from skiplist cimport *

def roll_median_c(ndarray[float64_t] arg, int win, int minp):
    cdef double val, res, prev
    cdef:
        int ret=0
        skiplist_t *sl
        Py_ssize_t midpoint, nobs = 0, i


    cdef Py_ssize_t N = len(arg)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    sl = skiplist_init(win)

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = arg[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist_insert(sl, val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = arg[i]

        if i > win - 1:
            prev = arg[i - win]

            if prev == prev:
                skiplist_remove(sl, prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist_insert(sl, val)

        if nobs >= minp:
            midpoint = nobs / 2
            if nobs % 2:
                res = skiplist_get(sl, midpoint, &ret)
            else:
                res = (skiplist_get(sl, midpoint, &ret) +
                       skiplist_get(sl, (midpoint - 1), &ret)) / 2
        else:
            res = NaN

        output[i] = res

    skiplist_destroy(sl)

    return output

def roll_median_cython(ndarray input, int win, int minp):
    '''
    O(N log(window)) implementation using skip list
    '''
    return _roll_skiplist_op(input, win, minp, _get_median)

# Unfortunately had to resort to some hackery here, would like for
# Cython to be able to get this right.

cdef double_t _get_median(object sl, int nobs, int minp):
    cdef Py_ssize_t midpoint
    cdef IndexableSkiplist skiplist = <IndexableSkiplist> sl
    if nobs >= minp:
        midpoint = nobs / 2
        if nobs % 2:
            return skiplist.get(midpoint)
        else:
            return (skiplist.get(midpoint) +
                    skiplist.get(midpoint - 1)) / 2
    else:
        return NaN

#----------------------------------------------------------------------

# Moving maximum / minimum code taken from Bottleneck under the terms
# of its Simplified BSD license
# https://github.com/kwgoodman/bottleneck

cdef struct pairs:
    double value
    int death

from libc cimport stdlib

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_max2(ndarray[float64_t] a, int window, int minp):
    "Moving max of 1d array of dtype=float64 along axis=0 ignoring NaNs."
    cdef np.float64_t ai, aold
    cdef Py_ssize_t count
    cdef pairs* ring
    cdef pairs* minpair
    cdef pairs* end
    cdef pairs* last
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef np.npy_intp *dims = [n0]
    cdef np.ndarray[np.float64_t, ndim=1] y = PyArray_EMPTY(1, dims,
		NPY_float64, 0)

    if window < 1:
        raise ValueError('Invalid window size %d'
                         % (window))

    if minp > window:
        raise ValueError('Invalid min_periods size %d greater than window %d'
                        % (minp, window))

    minp = _check_minp(window, minp, n0)

    window = min(window, n0)

    ring = <pairs*>stdlib.malloc(window * sizeof(pairs))
    end = ring + window
    last = ring

    minpair = ring
    ai = a[0]
    if ai == ai:
        minpair.value = ai
    else:
        minpair.value = MINfloat64
    minpair.death = window

    count = 0
    for i0 in range(n0):
        ai = a[i0]
        if ai == ai:
            count += 1
        else:
            ai = MINfloat64
        if i0 >= window:
            aold = a[i0 - window]
            if aold == aold:
                count -= 1
        if minpair.death == i0:
            minpair += 1
            if minpair >= end:
                minpair = ring
        if ai >= minpair.value:
            minpair.value = ai
            minpair.death = i0 + window
            last = minpair
        else:
            while last.value <= ai:
                if last == ring:
                    last = end
                last -= 1
            last += 1
            if last == end:
                last = ring
            last.value = ai
            last.death = i0 + window
        if count >= minp:
            y[i0] = minpair.value
        else:
            y[i0] = NaN

    for i0 in range(minp - 1):
        y[i0] = NaN

    stdlib.free(ring)
    return y

def roll_max(ndarray input, int win, int minp):
    '''
    O(N log(window)) implementation using skip list
    '''
    return _roll_skiplist_op(input, win, minp, _get_max)


cdef double_t _get_max(object skiplist, int nobs, int minp):
    if nobs >= minp:
        return <IndexableSkiplist> skiplist.get(nobs - 1)
    else:
        return NaN

def roll_min(ndarray input, int win, int minp):
    '''
    O(N log(window)) implementation using skip list
    '''
    return _roll_skiplist_op(input, win, minp, _get_min)

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_min2(np.ndarray[np.float64_t, ndim=1] a, int window, int minp):
    "Moving min of 1d array of dtype=float64 along axis=0 ignoring NaNs."
    cdef np.float64_t ai, aold
    cdef Py_ssize_t count
    cdef pairs* ring
    cdef pairs* minpair
    cdef pairs* end
    cdef pairs* last
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef np.npy_intp *dims = [n0]
    cdef np.ndarray[np.float64_t, ndim=1] y = PyArray_EMPTY(1, dims,
		NPY_float64, 0)

    if window < 1:
        raise ValueError('Invalid window size %d'
                         % (window))

    if minp > window:
        raise ValueError('Invalid min_periods size %d greater than window %d'
                        % (minp, window))

    window = min(window, n0)

    minp = _check_minp(window, minp, n0)

    ring = <pairs*>stdlib.malloc(window * sizeof(pairs))
    end = ring + window
    last = ring

    minpair = ring
    ai = a[0]
    if ai == ai:
        minpair.value = ai
    else:
        minpair.value = MAXfloat64
    minpair.death = window

    count = 0
    for i0 in range(n0):
        ai = a[i0]
        if ai == ai:
            count += 1
        else:
            ai = MAXfloat64
        if i0 >= window:
            aold = a[i0 - window]
            if aold == aold:
                count -= 1
        if minpair.death == i0:
            minpair += 1
            if minpair >= end:
                minpair = ring
        if ai <= minpair.value:
            minpair.value = ai
            minpair.death = i0 + window
            last = minpair
        else:
            while last.value >= ai:
                if last == ring:
                    last = end
                last -= 1
            last += 1
            if last == end:
                last = ring
            last.value = ai
            last.death = i0 + window
        if count >= minp:
            y[i0] = minpair.value
        else:
            y[i0] = NaN

    for i0 in range(minp - 1):
        y[i0] = NaN

    stdlib.free(ring)
    return y

cdef double_t _get_min(object skiplist, int nobs, int minp):
    if nobs >= minp:
        return <IndexableSkiplist> skiplist.get(0)
    else:
        return NaN

def roll_quantile(ndarray[float64_t, cast=True] input, int win,
                  int minp, double quantile):
    '''
    O(N log(window)) implementation using skip list
    '''
    cdef double val, prev, midpoint
    cdef IndexableSkiplist skiplist
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    skiplist = IndexableSkiplist(win)

    minp = _check_minp(win, minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist.insert(val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]

            if prev == prev:
                skiplist.remove(prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist.insert(val)

        if nobs >= minp:
            idx = int((quantile / 1.) * (nobs - 1))
            output[i] = skiplist.get(idx)
        else:
            output[i] = NaN

    return output

def roll_generic(ndarray[float64_t, cast=True] input, int win,
                 int minp, object func):
    cdef ndarray[double_t] output, counts, bufarr
    cdef Py_ssize_t i, n
    cdef float64_t *buf, *oldbuf

    if not input.flags.c_contiguous:
        input = input.copy('C')

    buf = <float64_t*> input.data

    n = len(input)
    if n == 0:
        return input

    minp = _check_minp(win, minp, n)
    output = np.empty(n, dtype=float)
    counts = roll_sum(np.isfinite(input).astype(float), win, minp)

    bufarr = np.empty(win, dtype=float)
    oldbuf = <float64_t*> bufarr.data

    n = len(input)
    for i from 0 <= i < int_min(win, n):
        if counts[i] >= minp:
            output[i] = func(input[int_max(i - win + 1, 0) : i + 1])
        else:
            output[i] = NaN

    for i from win <= i < n:
        buf = buf + 1
        bufarr.data = <char*> buf
        if counts[i] >= minp:
            output[i] = func(bufarr)
        else:
            output[i] = NaN

    bufarr.data = <char*> oldbuf

    return output


def roll_window(ndarray[float64_t, ndim=1, cast=True] input,
                ndarray[float64_t, ndim=1, cast=True] weights,
                int minp, bint avg=True, bint avg_wgt=False):
    """
    Assume len(weights) << len(input)
    """
    cdef:
        ndarray[double_t] output, tot_wgt, counts
        Py_ssize_t in_i, win_i, win_n, win_k, in_n, in_k
        float64_t val_in, val_win, c, w

    in_n = len(input)
    win_n = len(weights)
    output = np.zeros(in_n, dtype=float)
    counts = np.zeros(in_n, dtype=float)
    if avg:
        tot_wgt = np.zeros(in_n, dtype=float)

    minp = _check_minp(len(weights), minp, in_n)

    if avg_wgt:
        for win_i from 0 <= win_i < win_n:
            val_win = weights[win_i]
            if val_win != val_win:
                continue

            for in_i from 0 <= in_i < in_n - (win_n - win_i) + 1:
                val_in = input[in_i]
                if val_in == val_in:
                    output[in_i + (win_n - win_i) - 1] += val_in * val_win
                    counts[in_i + (win_n - win_i) - 1] += 1
                    tot_wgt[in_i + (win_n - win_i) - 1] += val_win

        for in_i from 0 <= in_i < in_n:
            c = counts[in_i]
            if c < minp:
                output[in_i] = NaN
            else:
                w = tot_wgt[in_i]
                if w == 0:
                    output[in_i] = NaN
                else:
                    output[in_i] /= tot_wgt[in_i]

    else:
        for win_i from 0 <= win_i < win_n:
            val_win = weights[win_i]
            if val_win != val_win:
                continue

            for in_i from 0 <= in_i < in_n - (win_n - win_i) + 1:
                val_in = input[in_i]

                if val_in == val_in:
                    output[in_i + (win_n - win_i) - 1] += val_in * val_win
                    counts[in_i + (win_n - win_i) - 1] += 1

        for in_i from 0 <= in_i < in_n:
            c = counts[in_i]
            if c < minp:
                output[in_i] = NaN
            elif avg:
                output[in_i] /= c

    return output


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
        # vecs[i] = <int64_t *> (<ndarray> list_of_arrays[i]).data

        arr = list_of_arrays[i]
        vecs[i] = <int64_t *> arr.data
    # assume uniqueness??

    for i from 1 <= i < n:
        for k from 0 <= k < nlevels:
            cur = vecs[k][i]
            pre = vecs[k][i-1]
            if cur == pre:
                continue
            elif cur > pre:
                break
            else:
                return False
    free(vecs)
    return True


@cython.boundscheck(False)
def groupby_indices(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] labels, counts, arr, seen
        int64_t loc
        dict ids = {}
        object val
        int64_t k

    ids, labels, counts = group_labels(values)
    seen = np.zeros_like(counts)

    # try not to get in trouble here...
    cdef int64_t **vecs = <int64_t **> malloc(len(ids) * sizeof(int64_t*))
    result = {}
    for i from 0 <= i < len(counts):
        arr = np.empty(counts[i], dtype=np.int64)
        result[ids[i]] = arr
        vecs[i] = <int64_t *> arr.data

    for i from 0 <= i < n:
        k = labels[i]

        # was NaN
        if k == -1:
            continue

        loc = seen[k]
        vecs[k][loc] = i
        seen[k] = loc + 1

    free(vecs)

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def group_labels(ndarray[object] values):
    '''
    Compute label vector from input values and associated useful data

    Returns
    -------
    '''
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] labels = np.empty(n, dtype=np.int64)
        ndarray[int64_t] counts = np.empty(n, dtype=np.int64)
        dict ids = {}, reverse = {}
        int64_t idx
        object val
        int64_t count = 0

    for i from 0 <= i < n:
        val = values[i]

        # is NaN
        if val != val:
            labels[i] = -1
            continue

        # for large number of groups, not doing try: except: makes a big
        # difference
        if val in ids:
            idx = ids[val]
            labels[i] = idx
            counts[idx] = counts[idx] + 1
        else:
            ids[val] = count
            reverse[count] = val
            labels[i] = count
            counts[count] = 1
            count += 1

    return reverse, labels, counts[:count].copy()


@cython.boundscheck(False)
@cython.wraparound(False)
def groupsort_indexer(ndarray[int64_t] index, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, loc, label, n
        ndarray[int64_t] counts, where, result

    # count group sizes, location 0 for NA
    counts = np.zeros(ngroups + 1, dtype=np.int64)
    n = len(index)
    for i from 0 <= i < n:
        counts[index[i] + 1] += 1

    # mark the start of each contiguous group of like-indexed data
    where = np.zeros(ngroups + 1, dtype=np.int64)
    for i from 1 <= i < ngroups + 1:
        where[i] = where[i - 1] + counts[i - 1]

    # this is our indexer
    result = np.zeros(n, dtype=np.int64)
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
    '''
    Only aggregates on axis=0
    '''
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
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[float64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.float64)
    resx = np.empty((<object> out).shape, dtype=object)

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
    '''
    Only aggregates on axis=0
    '''
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
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[float64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.float64)
    resx = np.empty((<object> out).shape, dtype=object)

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



#----------------------------------------------------------------------
# median

def group_median(ndarray[float64_t, ndim=2] out,
                 ndarray[int64_t] counts,
                 ndarray[float64_t, ndim=2] values,
                 ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, size
        ndarray[int64_t] _counts
        ndarray data
        float64_t* ptr
    ngroups = len(counts)
    N, K = (<object> values).shape

    indexer, _counts = groupsort_indexer(labels, ngroups)
    counts[:] = _counts[1:]

    data = np.empty((K, N), dtype=np.float64)
    ptr = <float64_t*> data.data

    take_2d_axis1_float64_float64(values.T, indexer, out=data)

    for i in range(K):
        # exclude NA group
        ptr += _counts[0]
        for j in range(ngroups):
            size = _counts[j + 1]
            out[j, i] = _median_linear(ptr, size)
            ptr += size


cdef inline float64_t _median_linear(float64_t* a, int n):
    cdef int i, j, na_count = 0
    cdef float64_t result
    cdef float64_t* tmp

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
        result = kth_smallest_c(a, n / 2, n)
    else:
        result = (kth_smallest_c(a, n / 2, n) +
                  kth_smallest_c(a, n / 2 - 1, n)) / 2

    if na_count:
        free(a)

    return result

include "join.pyx"
include "generated.pyx"
