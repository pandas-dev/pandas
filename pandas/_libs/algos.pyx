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
                    uint32_t, uint64_t, float32_t, float64_t)
cnp.import_array()


cimport pandas._libs.util as util
from pandas._libs.util cimport numeric, get_nat

from pandas._libs.khash cimport (
    khiter_t, kh_destroy_int64, kh_put_int64, kh_init_int64, kh_int64_t,
    kh_resize_int64, kh_get_int64)

import pandas._libs.missing as missing

cdef float64_t FP_ERR = 1e-13

cdef float64_t NaN = <float64_t>np.NaN

cdef int64_t NPY_NAT = get_nat()

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


class Infinity:
    """
    Provide a positive Infinity comparison method for ranking.
    """
    __lt__ = lambda self, other: False
    __le__ = lambda self, other: isinstance(other, Infinity)
    __eq__ = lambda self, other: isinstance(other, Infinity)
    __ne__ = lambda self, other: not isinstance(other, Infinity)
    __gt__ = lambda self, other: (not isinstance(other, Infinity) and
                                  not missing.checknull(other))
    __ge__ = lambda self, other: not missing.checknull(other)


class NegInfinity:
    """
    Provide a negative Infinity comparison method for ranking.
    """
    __lt__ = lambda self, other: (not isinstance(other, NegInfinity) and
                                  not missing.checknull(other))
    __le__ = lambda self, other: not missing.checknull(other)
    __eq__ = lambda self, other: isinstance(other, NegInfinity)
    __ne__ = lambda self, other: not isinstance(other, NegInfinity)
    __gt__ = lambda self, other: False
    __ge__ = lambda self, other: isinstance(other, NegInfinity)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef ndarray[int64_t, ndim=1] unique_deltas(const int64_t[:] arr):
    """
    Efficiently find the unique first-differences of the given array.

    Parameters
    ----------
    arr : ndarray[in64_t]

    Returns
    -------
    ndarray[int64_t]
        An ordered ndarray[int64_t]
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

    cdef int64_t **vecs = <int64_t**>malloc(nlevels * sizeof(int64_t*))
    for i in range(nlevels):
        arr = list_of_arrays[i]
        assert arr.dtype.name == 'int64'
        vecs[i] = <int64_t*>cnp.PyArray_DATA(arr)

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
def groupsort_indexer(const int64_t[:] index, Py_ssize_t ngroups):
    """
    Compute a 1-d indexer.

    The indexer is an ordering of the passed index,
    ordered by the groups.

    Parameters
    ----------
    index: int64 ndarray
        Mappings from group -> position.
    ngroups: int64
        Number of groups.

    Returns
    -------
    tuple
        1-d indexer ordered by groups, group counts.

    Notes
    -----
    This is a reverse of the label factorization process.
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
def kth_smallest(numeric[:] a, Py_ssize_t k) -> numeric:
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
def nancorr(const float64_t[:, :] mat, bint cov=0, minp=None):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        bint minpv
        ndarray[float64_t, ndim=2] result
        ndarray[uint8_t, ndim=2] mask
        int64_t nobs = 0
        float64_t vx, vy, sumx, sumy, sumxx, sumyy, meanx, meany, divisor

    N, K = (<object>mat).shape

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
def nancorr_spearman(const float64_t[:, :] mat, Py_ssize_t minp=1):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        ndarray[float64_t, ndim=2] result
        ndarray[float64_t, ndim=2] ranked_mat
        ndarray[float64_t, ndim=1] maskedx
        ndarray[float64_t, ndim=1] maskedy
        ndarray[uint8_t, ndim=2] mask
        int64_t nobs = 0
        float64_t vx, vy, sumx, sumxx, sumyy, mean, divisor

    N, K = (<object>mat).shape

    result = np.empty((K, K), dtype=np.float64)
    mask = np.isfinite(mat).view(np.uint8)

    ranked_mat = np.empty((N, K), dtype=np.float64)

    for i in range(K):
        ranked_mat[:, i] = rank_1d(mat[:, i])

    for xi in range(K):
        for yi in range(xi + 1):
            nobs = 0
            # Keep track of whether we need to recompute ranks
            all_ranks = True
            for i in range(N):
                all_ranks &= not (mask[i, xi] ^ mask[i, yi])
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
                        maskedx[j] = ranked_mat[i, xi]
                        maskedy[j] = ranked_mat[i, yi]
                        j += 1

                if not all_ranks:
                    maskedx = rank_1d(maskedx)
                    maskedy = rank_1d(maskedy)

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


def _validate_limit(nobs: int, limit=None) -> int:
    """
    Check that the `limit` argument is a positive integer.

    Parameters
    ----------
    nobs : int
    limit : object

    Returns
    -------
    int
        The limit.
    """
    if limit is None:
        lim = nobs
    else:
        if not util.is_integer_object(limit):
            raise ValueError('Limit must be an integer')
        if limit < 1:
            raise ValueError('Limit must be greater than 0')
        lim = limit

    return lim


@cython.boundscheck(False)
@cython.wraparound(False)
def pad(ndarray[algos_t] old, ndarray[algos_t] new, limit=None):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t, ndim=1] indexer
        algos_t cur, next_val
        int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer[:] = -1

    lim = _validate_limit(nright, limit)

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

        next_val = old[i + 1]

        while j < nright and cur <= new[j] < next_val:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next_val

    return indexer


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace(algos_t[:] values, const uint8_t[:] mask, limit=None):
    cdef:
        Py_ssize_t i, N
        algos_t val
        int lim, fill_count = 0

    N = len(values)

    # GH#2778
    if N == 0:
        return

    lim = _validate_limit(N, limit)

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
def pad_2d_inplace(algos_t[:, :] values, const uint8_t[:, :] mask, limit=None):
    cdef:
        Py_ssize_t i, j, N, K
        algos_t val
        int lim, fill_count = 0

    K, N = (<object>values).shape

    # GH#2778
    if N == 0:
        return

    lim = _validate_limit(N, limit)

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
def backfill(ndarray[algos_t] old, ndarray[algos_t] new, limit=None):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t, ndim=1] indexer
        algos_t cur, prev
        int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer[:] = -1

    lim = _validate_limit(nright, limit)

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
def backfill_inplace(algos_t[:] values, const uint8_t[:] mask, limit=None):
    cdef:
        Py_ssize_t i, N
        algos_t val
        int lim, fill_count = 0

    N = len(values)

    # GH#2778
    if N == 0:
        return

    lim = _validate_limit(N, limit)

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
def backfill_2d_inplace(algos_t[:, :] values,
                        const uint8_t[:, :] mask,
                        limit=None):
    cdef:
        Py_ssize_t i, j, N, K
        algos_t val
        int lim, fill_count = 0

    K, N = (<object>values).shape

    # GH#2778
    if N == 0:
        return

    lim = _validate_limit(N, limit)

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
def is_monotonic(ndarray[algos_t, ndim=1] arr, bint timelike):
    """
    Returns
    -------
    tuple
        is_monotonic_inc : bool
        is_monotonic_dec : bool
        is_unique : bool
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
        if arr[0] != arr[0] or (timelike and <int64_t>arr[0] == NPY_NAT):
            # single value is NaN
            return False, False, True
        else:
            return True, True, True
    elif n < 2:
        return True, True, True

    if timelike and <int64_t>arr[0] == NPY_NAT:
        return False, False, True

    if algos_t is not object:
        with nogil:
            prev = arr[0]
            for i in range(1, n):
                cur = arr[i]
                if timelike and <int64_t>cur == NPY_NAT:
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
            if timelike and <int64_t>cur == NPY_NAT:
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


# ----------------------------------------------------------------------
# rank_1d, rank_2d
# ----------------------------------------------------------------------

ctypedef fused rank_t:
    object
    float64_t
    uint64_t
    int64_t


@cython.wraparound(False)
@cython.boundscheck(False)
def rank_1d(rank_t[:] in_arr, ties_method='average',
            ascending=True, na_option='keep', pct=False):
    """
    Fast NaN-friendly version of ``scipy.stats.rankdata``.
    """
    cdef:
        Py_ssize_t i, j, n, dups = 0, total_tie_count = 0, non_na_idx = 0

        ndarray[rank_t] sorted_data, values

        ndarray[float64_t] ranks
        ndarray[int64_t] argsorted
        ndarray[uint8_t, cast=True] sorted_mask

        rank_t val, nan_value

        float64_t sum_ranks = 0
        int tiebreak = 0
        bint keep_na = 0
        bint isnan, condition
        float64_t count = 0.0

    tiebreak = tiebreakers[ties_method]

    if rank_t is float64_t:
        values = np.asarray(in_arr).copy()
    elif rank_t is object:
        values = np.array(in_arr, copy=True)

        if values.dtype != np.object_:
            values = values.astype('O')
    else:
        values = np.asarray(in_arr)

    keep_na = na_option == 'keep'

    if rank_t is object:
        mask = missing.isnaobj(values)
    elif rank_t is float64_t:
        mask = np.isnan(values)
    elif rank_t is int64_t:
        mask = values == NPY_NAT

        # create copy in case of NPY_NAT
        # values are mutated inplace
        if mask.any():
            values = values.copy()

    # double sort first by mask and then by values to ensure nan values are
    # either at the beginning or the end. mask/(~mask) controls padding at
    # tail or the head
    if rank_t is not uint64_t:
        if ascending ^ (na_option == 'top'):
            if rank_t is object:
                nan_value = Infinity()
            elif rank_t is float64_t:
                nan_value = np.inf
            elif rank_t is int64_t:
                nan_value = np.iinfo(np.int64).max

            order = (values, mask)
        else:
            if rank_t is object:
                nan_value = NegInfinity()
            elif rank_t is float64_t:
                nan_value = -np.inf
            elif rank_t is int64_t:
                nan_value = np.iinfo(np.int64).min

            order = (values, ~mask)
        np.putmask(values, mask, nan_value)
    else:
        mask = np.zeros(shape=len(values), dtype=bool)
        order = (values, mask)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    if rank_t is object:
        _as = np.lexsort(keys=order)
    else:
        if tiebreak == TIEBREAK_FIRST:
            # need to use a stable sort here
            _as = np.lexsort(keys=order)
            if not ascending:
                tiebreak = TIEBREAK_FIRST_DESCENDING
        else:
            _as = np.lexsort(keys=order)

    if not ascending:
        _as = _as[::-1]

    sorted_data = values.take(_as)
    sorted_mask = mask.take(_as)
    _indices = np.diff(sorted_mask.astype(int)).nonzero()[0]
    non_na_idx = _indices[0] if len(_indices) > 0 else -1
    argsorted = _as.astype('i8')

    if rank_t is object:
        # TODO: de-duplicate once cython supports conditional nogil
        for i in range(n):
            sum_ranks += i + 1
            dups += 1

            val = sorted_data[i]

            if rank_t is not uint64_t:
                isnan = sorted_mask[i]
                if isnan and keep_na:
                    ranks[argsorted[i]] = NaN
                    continue

            count += 1.0

            if rank_t is object:
                condition = (
                    i == n - 1 or
                    are_diff(sorted_data[i + 1], val) or
                    i == non_na_idx
                )
            else:
                condition = (
                    i == n - 1 or
                    sorted_data[i + 1] != val or
                    i == non_na_idx
                )

            if condition:

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
                    if rank_t is object:
                        raise ValueError('first not supported for non-numeric data')
                    else:
                        for j in range(i - dups + 1, i + 1):
                            ranks[argsorted[j]] = j + 1
                elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                    for j in range(i - dups + 1, i + 1):
                        ranks[argsorted[j]] = 2 * i - j - dups + 2
                elif tiebreak == TIEBREAK_DENSE:
                    total_tie_count += 1
                    for j in range(i - dups + 1, i + 1):
                        ranks[argsorted[j]] = total_tie_count
                sum_ranks = dups = 0

    else:
        with nogil:
            # TODO: why does the 2d version not have a nogil block?
            for i in range(n):
                sum_ranks += i + 1
                dups += 1

                val = sorted_data[i]

                if rank_t is not uint64_t:
                    isnan = sorted_mask[i]
                    if isnan and keep_na:
                        ranks[argsorted[i]] = NaN
                        continue

                count += 1.0

                if rank_t is object:
                    condition = (
                        i == n - 1 or
                        are_diff(sorted_data[i + 1], val) or
                        i == non_na_idx
                    )
                else:
                    condition = (
                        i == n - 1 or
                        sorted_data[i + 1] != val or
                        i == non_na_idx
                    )

                if condition:

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
                        if rank_t is object:
                            raise ValueError('first not supported for non-numeric data')
                        else:
                            for j in range(i - dups + 1, i + 1):
                                ranks[argsorted[j]] = j + 1
                    elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                        for j in range(i - dups + 1, i + 1):
                            ranks[argsorted[j]] = 2 * i - j - dups + 2
                    elif tiebreak == TIEBREAK_DENSE:
                        total_tie_count += 1
                        for j in range(i - dups + 1, i + 1):
                            ranks[argsorted[j]] = total_tie_count
                    sum_ranks = dups = 0

    if pct:
        if tiebreak == TIEBREAK_DENSE:
            return ranks / total_tie_count
        else:
            return ranks / count
    else:
        return ranks


def rank_2d(rank_t[:, :] in_arr, axis=0, ties_method='average',
            ascending=True, na_option='keep', pct=False):
    """
    Fast NaN-friendly version of ``scipy.stats.rankdata``.
    """
    cdef:
        Py_ssize_t i, j, z, k, n, dups = 0, total_tie_count = 0

        Py_ssize_t infs

        ndarray[float64_t, ndim=2] ranks
        ndarray[rank_t, ndim=2] values

        ndarray[int64_t, ndim=2] argsorted

        rank_t val, nan_value

        float64_t sum_ranks = 0
        int tiebreak = 0
        bint keep_na = 0
        float64_t count = 0.0
        bint condition, skip_condition

    tiebreak = tiebreakers[ties_method]

    keep_na = na_option == 'keep'

    if axis == 0:
        values = np.asarray(in_arr).T.copy()
    else:
        values = np.asarray(in_arr).copy()

    if rank_t is object:
        if values.dtype != np.object_:
            values = values.astype('O')

    if rank_t is not uint64_t:
        if ascending ^ (na_option == 'top'):
            if rank_t is object:
                nan_value = Infinity()
            elif rank_t is float64_t:
                nan_value = np.inf
            elif rank_t is int64_t:
                nan_value = np.iinfo(np.int64).max

        else:
            if rank_t is object:
                nan_value = NegInfinity()
            elif rank_t is float64_t:
                nan_value = -np.inf
            elif rank_t is int64_t:
                nan_value = NPY_NAT

        if rank_t is object:
            mask = missing.isnaobj2d(values)
        elif rank_t is float64_t:
            mask = np.isnan(values)
        elif rank_t is int64_t:
            mask = values == NPY_NAT

        np.putmask(values, mask, nan_value)

    n, k = (<object>values).shape
    ranks = np.empty((n, k), dtype='f8')

    if rank_t is object:
        try:
            _as = values.argsort(1)
        except TypeError:
            values = in_arr
            for i in range(len(values)):
                ranks[i] = rank_1d(in_arr[i], ties_method=ties_method,
                                   ascending=ascending, pct=pct)
            if axis == 0:
                return ranks.T
            else:
                return ranks
    else:
        if tiebreak == TIEBREAK_FIRST:
            # need to use a stable sort here
            _as = values.argsort(axis=1, kind='mergesort')
            if not ascending:
                tiebreak = TIEBREAK_FIRST_DESCENDING
        else:
            _as = values.argsort(1)

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d(values, _as)
    argsorted = _as.astype('i8')

    for i in range(n):
        if rank_t is object:
            dups = sum_ranks = infs = 0
        else:
            dups = sum_ranks = 0

        total_tie_count = 0
        count = 0.0
        for j in range(k):
            if rank_t is not object:
                sum_ranks += j + 1
                dups += 1

            val = values[i, j]

            if rank_t is not uint64_t:
                if rank_t is object:
                    skip_condition = (val is nan_value) and keep_na
                else:
                    skip_condition = (val == nan_value) and keep_na
                if skip_condition:
                    ranks[i, argsorted[i, j]] = NaN

                    if rank_t is object:
                        infs += 1

                    continue

            count += 1.0

            if rank_t is object:
                sum_ranks += (j - infs) + 1
                dups += 1

            if rank_t is object:
                condition = j == k - 1 or are_diff(values[i, j + 1], val)
            else:
                condition = j == k - 1 or values[i, j + 1] != val

            if condition:
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
                    if rank_t is object:
                        raise ValueError('first not supported for non-numeric data')
                    else:
                        for z in range(j - dups + 1, j + 1):
                            ranks[i, argsorted[i, z]] = z + 1
                elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = 2 * j - z - dups + 2
                elif tiebreak == TIEBREAK_DENSE:
                    total_tie_count += 1
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = total_tie_count
                sum_ranks = dups = 0
        if pct:
            if tiebreak == TIEBREAK_DENSE:
                ranks[i, :] /= total_tie_count
            else:
                ranks[i, :] /= count
    if axis == 0:
        return ranks.T
    else:
        return ranks


ctypedef fused diff_t:
    float64_t
    float32_t
    int8_t
    int16_t
    int32_t
    int64_t

ctypedef fused out_t:
    float32_t
    float64_t


@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d(diff_t[:, :] arr,
            out_t[:, :] out,
            Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy, start, stop
        bint f_contig = arr.is_f_contig()

    # Disable for unsupported dtype combinations,
    #  see https://github.com/cython/cython/issues/2646
    if (out_t is float32_t
            and not (diff_t is float32_t or diff_t is int8_t or diff_t is int16_t)):
        raise NotImplementedError
    elif (out_t is float64_t
          and (diff_t is float32_t or diff_t is int8_t or diff_t is int16_t)):
        raise NotImplementedError
    else:
        # We put this inside an indented else block to avoid cython build
        #  warnings about unreachable code
        sx, sy = (<object>arr).shape
        with nogil:
            if f_contig:
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


# generated from template
include "algos_common_helper.pxi"
include "algos_take_helper.pxi"
