import cython
from cython import Py_ssize_t
from cython cimport floating

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t, complex64_t, complex128_t)
cnp.import_array()

cdef extern from "numpy/npy_math.h":
    float64_t NAN "NPY_NAN"

from pandas._libs.util cimport numeric, get_nat

from pandas._libs.algos cimport (swap, TiebreakEnumType, TIEBREAK_AVERAGE,
                                 TIEBREAK_MIN, TIEBREAK_MAX, TIEBREAK_FIRST,
                                 TIEBREAK_DENSE)
from pandas._libs.algos import (take_2d_axis1_float64_float64,
                                groupsort_indexer, tiebreakers)

from pandas._libs.missing cimport checknull

cdef int64_t NPY_NAT = get_nat()
_int64_max = np.iinfo(np.int64).max

cdef float64_t NaN = <float64_t>np.NaN

cdef enum InterpolationEnumType:
    INTERPOLATION_LINEAR,
    INTERPOLATION_LOWER,
    INTERPOLATION_HIGHER,
    INTERPOLATION_NEAREST,
    INTERPOLATION_MIDPOINT


cdef inline float64_t median_linear(float64_t* a, int n) nogil:
    cdef:
        int i, j, na_count = 0
        float64_t result
        float64_t* tmp

    if n == 0:
        return NaN

    # count NAs
    for i in range(n):
        if a[i] != a[i]:
            na_count += 1

    if na_count:
        if na_count == n:
            return NaN

        tmp = <float64_t*>malloc((n - na_count) * sizeof(float64_t))

        j = 0
        for i in range(n):
            if a[i] == a[i]:
                tmp[j] = a[i]
                j += 1

        a = tmp
        n -= na_count

    if n % 2:
        result = kth_smallest_c( a, n // 2, n)
    else:
        result = (kth_smallest_c(a, n // 2, n) +
                  kth_smallest_c(a, n // 2 - 1, n)) / 2

    if na_count:
        free(a)

    return result


# TODO: Is this redundant with algos.kth_smallest
cdef inline float64_t kth_smallest_c(float64_t* a,
                                     Py_ssize_t k,
                                     Py_ssize_t n) nogil:
    cdef:
        Py_ssize_t i, j, l, m
        float64_t x, t

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_median_float64(ndarray[float64_t, ndim=2] out,
                         ndarray[int64_t] counts,
                         ndarray[float64_t, ndim=2] values,
                         ndarray[int64_t] labels,
                         Py_ssize_t min_count=-1):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, ngroups, size
        ndarray[int64_t] _counts
        ndarray[float64_t, ndim=2] data
        float64_t* ptr

    assert min_count == -1, "'min_count' only used in add and prod"

    ngroups = len(counts)
    N, K = (<object>values).shape

    indexer, _counts = groupsort_indexer(labels, ngroups)
    counts[:] = _counts[1:]

    data = np.empty((K, N), dtype=np.float64)
    ptr = <float64_t*>cnp.PyArray_DATA(data)

    take_2d_axis1_float64_float64(values.T, indexer, out=data)

    with nogil:

        for i in range(K):
            # exclude NA group
            ptr += _counts[0]
            for j in range(ngroups):
                size = _counts[j + 1]
                out[j, i] = median_linear(ptr, size)
                ptr += size


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cumprod_float64(float64_t[:, :] out,
                          const float64_t[:, :] values,
                          const int64_t[:] labels,
                          int ngroups,
                          bint is_datetimelike,
                          bint skipna=True):
    """
    Cumulative product of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : float64 array
        Array to store cumprod in.
    values : float64 array
        Values to take cumprod of.
    labels : int64 array
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        Always false, `values` is never datetime-like.
    skipna : bool
        If true, ignore nans in `values`.

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        float64_t val
        float64_t[:, :] accum
        int64_t lab

    N, K = (<object>values).shape
    accum = np.ones((ngroups, K), dtype=np.float64)

    with nogil:
        for i in range(N):
            lab = labels[i]

            if lab < 0:
                continue
            for j in range(K):
                val = values[i, j]
                if val == val:
                    accum[lab, j] *= val
                    out[i, j] = accum[lab, j]
                else:
                    out[i, j] = NaN
                    if not skipna:
                        accum[lab, j] = NaN
                        break


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cumsum(numeric[:, :] out,
                 numeric[:, :] values,
                 const int64_t[:] labels,
                 int ngroups,
                 is_datetimelike,
                 bint skipna=True):
    """
    Cumulative sum of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : array
        Array to store cumsum in.
    values : array
        Values to take cumsum of.
    labels : int64 array
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    skipna : bool
        If true, ignore nans in `values`.

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        numeric val
        numeric[:, :] accum
        int64_t lab

    N, K = (<object>values).shape
    accum = np.zeros((ngroups, K), dtype=np.asarray(values).dtype)

    with nogil:
        for i in range(N):
            lab = labels[i]

            if lab < 0:
                continue
            for j in range(K):
                val = values[i, j]

                if numeric == float32_t or numeric == float64_t:
                    if val == val:
                        accum[lab, j] += val
                        out[i, j] = accum[lab, j]
                    else:
                        out[i, j] = NaN
                        if not skipna:
                            accum[lab, j] = NaN
                            break
                else:
                    accum[lab, j] += val
                    out[i, j] = accum[lab, j]


@cython.boundscheck(False)
@cython.wraparound(False)
def group_shift_indexer(int64_t[:] out, const int64_t[:] labels,
                        int ngroups, int periods):
    cdef:
        Py_ssize_t N, i, j, ii
        int offset = 0, sign
        int64_t lab, idxer, idxer_slot
        int64_t[:] label_seen = np.zeros(ngroups, dtype=np.int64)
        int64_t[:, :] label_indexer

    N, = (<object>labels).shape

    if periods < 0:
        periods = -periods
        offset = N - 1
        sign = -1
    elif periods > 0:
        offset = 0
        sign = 1

    if periods == 0:
        with nogil:
            for i in range(N):
                out[i] = i
    else:
        # array of each previous indexer seen
        label_indexer = np.zeros((ngroups, periods), dtype=np.int64)
        with nogil:
            for i in range(N):
                # reverse iterator if shifting backwards
                ii = offset + sign * i
                lab = labels[ii]

                # Skip null keys
                if lab == -1:
                    out[ii] = -1
                    continue

                label_seen[lab] += 1

                idxer_slot = label_seen[lab] % periods
                idxer = label_indexer[lab, idxer_slot]

                if label_seen[lab] > periods:
                    out[ii] = idxer
                else:
                    out[ii] = -1

                label_indexer[lab, idxer_slot] = ii


@cython.wraparound(False)
@cython.boundscheck(False)
def group_fillna_indexer(ndarray[int64_t] out, ndarray[int64_t] labels,
                         ndarray[uint8_t] mask, object direction,
                         int64_t limit):
    """
    Indexes how to fill values forwards or backwards within a group.

    Parameters
    ----------
    out : array of int64_t values which this method will write its results to
        Missing values will be written to with a value of -1
    labels : array containing unique label for each group, with its ordering
        matching up to the corresponding record in `values`
    mask : array of int64_t values where a 1 indicates a missing value
    direction : {'ffill', 'bfill'}
        Direction for fill to be applied (forwards or backwards, respectively)
    limit : Consecutive values to fill before stopping, or -1 for no limit

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object
    """
    cdef:
        Py_ssize_t i, N
        int64_t[:] sorted_labels
        int64_t idx, curr_fill_idx=-1, filled_vals=0

    N = len(out)

    # Make sure all arrays are the same size
    assert N == len(labels) == len(mask)

    sorted_labels = np.argsort(labels, kind='mergesort').astype(
        np.int64, copy=False)
    if direction == 'bfill':
        sorted_labels = sorted_labels[::-1]

    with nogil:
        for i in range(N):
            idx = sorted_labels[i]
            if mask[idx] == 1:  # is missing
                # Stop filling once we've hit the limit
                if filled_vals >= limit and limit != -1:
                    curr_fill_idx = -1
                filled_vals += 1
            else:  # reset items when not missing
                filled_vals = 0
                curr_fill_idx = idx

            out[idx] = curr_fill_idx

            # If we move to the next group, reset
            # the fill_idx and counter
            if i == N - 1 or labels[idx] != labels[sorted_labels[i + 1]]:
                curr_fill_idx = -1
                filled_vals = 0


@cython.boundscheck(False)
@cython.wraparound(False)
def group_any_all(uint8_t[:] out,
                  const int64_t[:] labels,
                  const uint8_t[:] values,
                  const uint8_t[:] mask,
                  object val_test,
                  bint skipna):
    """
    Aggregated boolean values to show truthfulness of group elements.

    Parameters
    ----------
    out : array of values which this method will write its results to
    labels : array containing unique label for each group, with its
        ordering matching up to the corresponding record in `values`
    values : array containing the truth value of each element
    mask : array indicating whether a value is na or not
    val_test : str {'any', 'all'}
        String object dictating whether to use any or all truth testing
    skipna : boolean
        Flag to ignore nan values during truth testing

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object.
    The returned values will either be 0 or 1 (False or True, respectively).
    """
    cdef:
        Py_ssize_t i, N = len(labels)
        int64_t lab
        uint8_t flag_val

    if val_test == 'all':
        # Because the 'all' value of an empty iterable in Python is True we can
        # start with an array full of ones and set to zero when a False value
        # is encountered
        flag_val = 0
    elif val_test == 'any':
        # Because the 'any' value of an empty iterable in Python is False we
        # can start with an array full of zeros and set to one only if any
        # value encountered is True
        flag_val = 1
    else:
        raise ValueError("'bool_func' must be either 'any' or 'all'!")

    out[:] = 1 - flag_val

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0 or (skipna and mask[i]):
                continue

            if values[i] == flag_val:
                out[lab] = flag_val


# ----------------------------------------------------------------------
# group_add, group_prod, group_var, group_mean, group_ohlc
# ----------------------------------------------------------------------

ctypedef fused complexfloating_t:
    float64_t
    float32_t
    complex64_t
    complex128_t


@cython.wraparound(False)
@cython.boundscheck(False)
def _group_add(complexfloating_t[:, :] out,
               int64_t[:] counts,
               complexfloating_t[:, :] values,
               const int64_t[:] labels,
               Py_ssize_t min_count=0):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        complexfloating_t val, count
        complexfloating_t[:, :] sumx
        int64_t[:, :] nobs

    if len(values) != len(labels):
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    sumx = np.zeros_like(out)

    N, K = (<object>values).shape

    with nogil:
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
                    if (complexfloating_t is complex64_t or
                            complexfloating_t is complex128_t):
                        # clang errors if we use += with these dtypes
                        sumx[lab, j] = sumx[lab, j] + val
                    else:
                        sumx[lab, j] += val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] < min_count:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j]


group_add_float32 = _group_add['float32_t']
group_add_float64 = _group_add['float64_t']
group_add_complex64 = _group_add['float complex']
group_add_complex128 = _group_add['double complex']


@cython.wraparound(False)
@cython.boundscheck(False)
def _group_prod(floating[:, :] out,
                int64_t[:] counts,
                floating[:, :] values,
                const int64_t[:] labels,
                Py_ssize_t min_count=0):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        floating val, count
        floating[:, :] prodx
        int64_t[:, :] nobs

    if not len(values) == len(labels):
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    prodx = np.ones_like(out)

    N, K = (<object>values).shape

    with nogil:
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

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] < min_count:
                    out[i, j] = NAN
                else:
                    out[i, j] = prodx[i, j]


group_prod_float32 = _group_prod['float']
group_prod_float64 = _group_prod['double']


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def _group_var(floating[:, :] out,
               int64_t[:] counts,
               floating[:, :] values,
               const int64_t[:] labels,
               Py_ssize_t min_count=-1):
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        floating val, ct, oldmean
        floating[:, :] mean
        int64_t[:, :] nobs

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    mean = np.zeros_like(out)

    N, K = (<object>values).shape

    out[:, :] = 0.0

    with nogil:
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
                    oldmean = mean[lab, j]
                    mean[lab, j] += (val - oldmean) / nobs[lab, j]
                    out[lab, j] += (val - mean[lab, j]) * (val - oldmean)

        for i in range(ncounts):
            for j in range(K):
                ct = nobs[i, j]
                if ct < 2:
                    out[i, j] = NAN
                else:
                    out[i, j] /= (ct - 1)


group_var_float32 = _group_var['float']
group_var_float64 = _group_var['double']


@cython.wraparound(False)
@cython.boundscheck(False)
def _group_mean(floating[:, :] out,
                int64_t[:] counts,
                floating[:, :] values,
                const int64_t[:] labels,
                Py_ssize_t min_count=-1):
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        floating val, count
        floating[:, :] sumx
        int64_t[:, :] nobs

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    sumx = np.zeros_like(out)

    N, K = (<object>values).shape

    with nogil:
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

        for i in range(ncounts):
            for j in range(K):
                count = nobs[i, j]
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j] / count


group_mean_float32 = _group_mean['float']
group_mean_float64 = _group_mean['double']


@cython.wraparound(False)
@cython.boundscheck(False)
def _group_ohlc(floating[:, :] out,
                int64_t[:] counts,
                floating[:, :] values,
                const int64_t[:] labels,
                Py_ssize_t min_count=-1):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab
        floating val, count
        Py_ssize_t ngroups = len(counts)

    assert min_count == -1, "'min_count' only used in add and prod"

    if len(labels) == 0:
        return

    N, K = (<object>values).shape

    if out.shape[1] != 4:
        raise ValueError('Output array must have 4 columns')

    if K > 1:
        raise NotImplementedError("Argument 'values' must have only "
                                  "one dimension")
    out[:] = np.nan

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab == -1:
                continue

            counts[lab] += 1
            val = values[i, 0]
            if val != val:
                continue

            if out[lab, 0] != out[lab, 0]:
                out[lab, 0] = out[lab, 1] = out[lab, 2] = out[lab, 3] = val
            else:
                out[lab, 1] = max(out[lab, 1], val)
                out[lab, 2] = min(out[lab, 2], val)
                out[lab, 3] = val


group_ohlc_float32 = _group_ohlc['float']
group_ohlc_float64 = _group_ohlc['double']


@cython.boundscheck(False)
@cython.wraparound(False)
def group_quantile(ndarray[float64_t] out,
                   ndarray[int64_t] labels,
                   numeric[:] values,
                   ndarray[uint8_t] mask,
                   float64_t q,
                   object interpolation):
    """
    Calculate the quantile per group.

    Parameters
    ----------
    out : ndarray
        Array of aggregated values that will be written to.
    labels : ndarray
        Array containing the unique group labels.
    values : ndarray
        Array containing the values to apply the function against.
    q : float
        The quantile value to search for.

    Notes
    -----
    Rather than explicitly returning a value, this function modifies the
    provided `out` parameter.
    """
    cdef:
        Py_ssize_t i, N=len(labels), ngroups, grp_sz, non_na_sz
        Py_ssize_t grp_start=0, idx=0
        int64_t lab
        uint8_t interp
        float64_t q_idx, frac, val, next_val
        ndarray[int64_t] counts, non_na_counts, sort_arr

    assert values.shape[0] == N

    if not (0 <= q <= 1):
        raise ValueError(f"'q' must be between 0 and 1. Got '{q}' instead")

    inter_methods = {
        'linear': INTERPOLATION_LINEAR,
        'lower': INTERPOLATION_LOWER,
        'higher': INTERPOLATION_HIGHER,
        'nearest': INTERPOLATION_NEAREST,
        'midpoint': INTERPOLATION_MIDPOINT,
    }
    interp = inter_methods[interpolation]

    counts = np.zeros_like(out, dtype=np.int64)
    non_na_counts = np.zeros_like(out, dtype=np.int64)
    ngroups = len(counts)

    # First figure out the size of every group
    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab == -1:  # NA group label
                continue

            counts[lab] += 1
            if not mask[i]:
                non_na_counts[lab] += 1

    # Get an index of values sorted by labels and then values
    order = (values, labels)
    sort_arr = np.lexsort(order).astype(np.int64, copy=False)

    with nogil:
        for i in range(ngroups):
            # Figure out how many group elements there are
            grp_sz = counts[i]
            non_na_sz = non_na_counts[i]

            if non_na_sz == 0:
                out[i] = NaN
            else:
                # Calculate where to retrieve the desired value
                # Casting to int will intentionally truncate result
                idx = grp_start + <int64_t>(q * <float64_t>(non_na_sz - 1))

                val = values[sort_arr[idx]]
                # If requested quantile falls evenly on a particular index
                # then write that index's value out. Otherwise interpolate
                q_idx = q * (non_na_sz - 1)
                frac = q_idx % 1

                if frac == 0.0 or interp == INTERPOLATION_LOWER:
                    out[i] = val
                else:
                    next_val = values[sort_arr[idx + 1]]
                    if interp == INTERPOLATION_LINEAR:
                        out[i] = val + (next_val - val) * frac
                    elif interp == INTERPOLATION_HIGHER:
                        out[i] = next_val
                    elif interp == INTERPOLATION_MIDPOINT:
                        out[i] = (val + next_val) / 2.0
                    elif interp == INTERPOLATION_NEAREST:
                        if frac > .5 or (frac == .5 and q > .5):  # Always OK?
                            out[i] = next_val
                        else:
                            out[i] = val

            # Increment the index reference in sorted_arr for the next group
            grp_start += grp_sz


# ----------------------------------------------------------------------
# group_nth, group_last, group_rank
# ----------------------------------------------------------------------

ctypedef fused rank_t:
    float64_t
    float32_t
    int64_t
    uint64_t
    object


cdef inline bint _treat_as_na(rank_t val, bint is_datetimelike) nogil:
    if rank_t is object:
        # Should never be used, but we need to avoid the `val != val` below
        #  or else cython will raise about gil acquisition.
        raise NotImplementedError

    elif rank_t is int64_t:
        return is_datetimelike and val == NPY_NAT
    elif rank_t is uint64_t:
        # There is no NA value for uint64
        return False
    else:
        return val != val


# GH#31710 use memorviews once cython 0.30 is released so we can
#  use `const rank_t[:, :] values`
@cython.wraparound(False)
@cython.boundscheck(False)
def group_last(rank_t[:, :] out,
               int64_t[:] counts,
               ndarray[rank_t, ndim=2] values,
               const int64_t[:] labels,
               Py_ssize_t min_count=-1):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        rank_t val
        ndarray[rank_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs
        bint runtime_error = False

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    if rank_t is object:
        resx = np.empty((<object>out).shape, dtype=object)
    else:
        resx = np.empty_like(out)

    N, K = (<object>values).shape

    if rank_t is object:
        # TODO: De-duplicate once conditional-nogil is available
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                if not checknull(val):
                    # NB: use _treat_as_na here once
                    #  conditional-nogil is available.
                    nobs[lab, j] += 1
                    resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = resx[i, j]
    else:
        with nogil:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    if not _treat_as_na(val, True):
                        # TODO: Sure we always want is_datetimelike=True?
                        nobs[lab, j] += 1
                        resx[lab, j] = val

            for i in range(ncounts):
                for j in range(K):
                    if nobs[i, j] == 0:
                        if rank_t is int64_t:
                            out[i, j] = NPY_NAT
                        elif rank_t is uint64_t:
                            runtime_error = True
                            break
                        else:
                            out[i, j] = NAN

                    else:
                        out[i, j] = resx[i, j]

    if runtime_error:
        # We cannot raise directly above because that is within a nogil
        #  block.
        raise RuntimeError("empty group with uint64_t")


# GH#31710 use memorviews once cython 0.30 is released so we can
#  use `const rank_t[:, :] values`
@cython.wraparound(False)
@cython.boundscheck(False)
def group_nth(rank_t[:, :] out,
              int64_t[:] counts,
              ndarray[rank_t, ndim=2] values,
              const int64_t[:] labels, int64_t rank=1,
              Py_ssize_t min_count=-1):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        rank_t val
        ndarray[rank_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs
        bint runtime_error = False

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    if rank_t is object:
        resx = np.empty((<object>out).shape, dtype=object)
    else:
        resx = np.empty_like(out)

    N, K = (<object>values).shape

    if rank_t is object:
        # TODO: De-duplicate once conditional-nogil is available
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                if not checknull(val):
                    # NB: use _treat_as_na here once
                    #  conditional-nogil is available.
                    nobs[lab, j] += 1
                    if nobs[lab, j] == rank:
                        resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = resx[i, j]

    else:
        with nogil:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    if not _treat_as_na(val, True):
                        # TODO: Sure we always want is_datetimelike=True?
                        nobs[lab, j] += 1
                        if nobs[lab, j] == rank:
                            resx[lab, j] = val

            for i in range(ncounts):
                for j in range(K):
                    if nobs[i, j] == 0:
                        if rank_t is int64_t:
                            out[i, j] = NPY_NAT
                        elif rank_t is uint64_t:
                            runtime_error = True
                            break
                        else:
                            out[i, j] = NAN
                    else:
                        out[i, j] = resx[i, j]

    if runtime_error:
        # We cannot raise directly above because that is within a nogil
        #  block.
        raise RuntimeError("empty group with uint64_t")


@cython.boundscheck(False)
@cython.wraparound(False)
def group_rank(float64_t[:, :] out,
               rank_t[:, :] values,
               const int64_t[:] labels,
               int ngroups,
               bint is_datetimelike, object ties_method="average",
               bint ascending=True, bint pct=False, object na_option="keep"):
    """
    Provides the rank of values within each group.

    Parameters
    ----------
    out : array of float64_t values which this method will write its results to
    values : array of rank_t values to be ranked
    labels : array containing unique label for each group, with its ordering
        matching up to the corresponding record in `values`
    ngroups : int
        This parameter is not used, is needed to match signatures of other
        groupby functions.
    is_datetimelike : bool, default False
        unused in this method but provided for call compatibility with other
        Cython transformations
    ties_method : {'average', 'min', 'max', 'first', 'dense'}, default
        'average'
        * average: average rank of group
        * min: lowest rank in group
        * max: highest rank in group
        * first: ranks assigned in order they appear in the array
        * dense: like 'min', but rank always increases by 1 between groups
    ascending : boolean, default True
        False for ranks by high (1) to low (N)
        na_option : {'keep', 'top', 'bottom'}, default 'keep'
    pct : boolean, default False
        Compute percentage rank of data within each group
    na_option : {'keep', 'top', 'bottom'}, default 'keep'
        * keep: leave NA values where they are
        * top: smallest rank if ascending
        * bottom: smallest rank if descending

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object
    """
    cdef:
        TiebreakEnumType tiebreak
        Py_ssize_t i, j, N, K, grp_start=0, dups=0, sum_ranks=0
        Py_ssize_t grp_vals_seen=1, grp_na_count=0, grp_tie_count=0
        ndarray[int64_t] _as
        ndarray[float64_t, ndim=2] grp_sizes
        ndarray[rank_t] masked_vals
        ndarray[uint8_t] mask
        bint keep_na
        rank_t nan_fill_val

    if rank_t is object:
        raise NotImplementedError("Cant do nogil")

    tiebreak = tiebreakers[ties_method]
    keep_na = na_option == 'keep'
    N, K = (<object>values).shape
    grp_sizes = np.ones_like(out)

    # Copy values into new array in order to fill missing data
    # with mask, without obfuscating location of missing data
    # in values array
    masked_vals = np.array(values[:, 0], copy=True)
    if rank_t is int64_t:
        mask = (masked_vals == NPY_NAT).astype(np.uint8)
    else:
        mask = np.isnan(masked_vals).astype(np.uint8)

    if ascending ^ (na_option == 'top'):
        if rank_t is int64_t:
            nan_fill_val = np.iinfo(np.int64).max
        elif rank_t is uint64_t:
            nan_fill_val = np.iinfo(np.uint64).max
        else:
            nan_fill_val = np.inf
        order = (masked_vals, mask, labels)
    else:
        if rank_t is int64_t:
            nan_fill_val = np.iinfo(np.int64).min
        elif rank_t is uint64_t:
            nan_fill_val = 0
        else:
            nan_fill_val = -np.inf

        order = (masked_vals, ~mask, labels)
    np.putmask(masked_vals, mask, nan_fill_val)

    # lexsort using labels, then mask, then actual values
    # each label corresponds to a different group value,
    # the mask helps you differentiate missing values before
    # performing sort on the actual values
    _as = np.lexsort(order).astype(np.int64, copy=False)

    if not ascending:
        _as = _as[::-1]

    with nogil:
        # Loop over the length of the value array
        # each incremental i value can be looked up in the _as array
        # that we sorted previously, which gives us the location of
        # that sorted value for retrieval back from the original
        # values / masked_vals arrays
        for i in range(N):
            # dups and sum_ranks will be incremented each loop where
            # the value / group remains the same, and should be reset
            # when either of those change
            # Used to calculate tiebreakers
            dups += 1
            sum_ranks += i - grp_start + 1

            # Update out only when there is a transition of values or labels.
            # When a new value or group is encountered, go back #dups steps(
            # the number of occurrence of current value) and assign the ranks
            # based on the the starting index of the current group (grp_start)
            # and the current index
            if (i == N - 1 or
                    (masked_vals[_as[i]] != masked_vals[_as[i+1]]) or
                    (mask[_as[i]] ^ mask[_as[i+1]]) or
                    (labels[_as[i]] != labels[_as[i+1]])):
                # if keep_na, check for missing values and assign back
                # to the result where appropriate
                if keep_na and mask[_as[i]]:
                    for j in range(i - dups + 1, i + 1):
                        out[_as[j], 0] = NaN
                        grp_na_count = dups
                elif tiebreak == TIEBREAK_AVERAGE:
                    for j in range(i - dups + 1, i + 1):
                        out[_as[j], 0] = sum_ranks / <float64_t>dups
                elif tiebreak == TIEBREAK_MIN:
                    for j in range(i - dups + 1, i + 1):
                        out[_as[j], 0] = i - grp_start - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for j in range(i - dups + 1, i + 1):
                        out[_as[j], 0] = i - grp_start + 1
                elif tiebreak == TIEBREAK_FIRST:
                    for j in range(i - dups + 1, i + 1):
                        if ascending:
                            out[_as[j], 0] = j + 1 - grp_start
                        else:
                            out[_as[j], 0] = 2 * i - j - dups + 2 - grp_start
                elif tiebreak == TIEBREAK_DENSE:
                    for j in range(i - dups + 1, i + 1):
                        out[_as[j], 0] = grp_vals_seen

                # look forward to the next value (using the sorting in _as)
                # if the value does not equal the current value then we need to
                # reset the dups and sum_ranks, knowing that a new value is
                # coming up. the conditional also needs to handle nan equality
                # and the end of iteration
                if (i == N - 1 or
                        (masked_vals[_as[i]] != masked_vals[_as[i+1]]) or
                        (mask[_as[i]] ^ mask[_as[i+1]])):
                    dups = sum_ranks = 0
                    grp_vals_seen += 1
                    grp_tie_count += 1

                # Similar to the previous conditional, check now if we are
                # moving to a new group. If so, keep track of the index where
                # the new group occurs, so the tiebreaker calculations can
                # decrement that from their position. fill in the size of each
                # group encountered (used by pct calculations later). also be
                # sure to reset any of the items helping to calculate dups
                if i == N - 1 or labels[_as[i]] != labels[_as[i+1]]:
                    if tiebreak != TIEBREAK_DENSE:
                        for j in range(grp_start, i + 1):
                            grp_sizes[_as[j], 0] = (i - grp_start + 1 -
                                                    grp_na_count)
                    else:
                        for j in range(grp_start, i + 1):
                            grp_sizes[_as[j], 0] = (grp_tie_count -
                                                    (grp_na_count > 0))
                    dups = sum_ranks = 0
                    grp_na_count = 0
                    grp_tie_count = 0
                    grp_start = i + 1
                    grp_vals_seen = 1

        if pct:
            for i in range(N):
                # We don't include NaN values in percentage
                # rankings, so we assign them percentages of NaN.
                if out[i, 0] != out[i, 0] or out[i, 0] == NAN:
                    out[i, 0] = NAN
                elif grp_sizes[i, 0] != 0:
                    out[i, 0] = out[i, 0] / grp_sizes[i, 0]


# ----------------------------------------------------------------------
# group_min, group_max
# ----------------------------------------------------------------------

# TODO: consider implementing for more dtypes
ctypedef fused groupby_t:
    float64_t
    float32_t
    int64_t
    uint64_t


@cython.wraparound(False)
@cython.boundscheck(False)
def group_max(groupby_t[:, :] out,
              int64_t[:] counts,
              ndarray[groupby_t, ndim=2] values,
              const int64_t[:] labels,
              Py_ssize_t min_count=-1):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        groupby_t val, count, nan_val
        ndarray[groupby_t, ndim=2] maxx
        bint runtime_error = False
        int64_t[:, :] nobs

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)

    maxx = np.empty_like(out)
    if groupby_t is int64_t:
        # Note: evaluated at compile-time
        maxx[:] = -_int64_max
        nan_val = NPY_NAT
    elif groupby_t is uint64_t:
        # NB: We do not define nan_val because there is no such thing
        #  for uint64_t.  We carefully avoid having to reference it in this
        #  case.
        maxx[:] = 0
    else:
        maxx[:] = -np.inf
        nan_val = NAN

    N, K = (<object>values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                if not _treat_as_na(val, True):
                    # TODO: Sure we always want is_datetimelike=True?
                    nobs[lab, j] += 1
                    if val > maxx[lab, j]:
                        maxx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    if groupby_t is uint64_t:
                        runtime_error = True
                        break
                    else:
                        out[i, j] = nan_val
                else:
                    out[i, j] = maxx[i, j]

    if runtime_error:
        # We cannot raise directly above because that is within a nogil
        #  block.
        raise RuntimeError("empty group with uint64_t")


@cython.wraparound(False)
@cython.boundscheck(False)
def group_min(groupby_t[:, :] out,
              int64_t[:] counts,
              ndarray[groupby_t, ndim=2] values,
              const int64_t[:] labels,
              Py_ssize_t min_count=-1):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        groupby_t val, count, nan_val
        ndarray[groupby_t, ndim=2] minx
        bint runtime_error = False
        int64_t[:, :] nobs

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)

    minx = np.empty_like(out)
    if groupby_t is int64_t:
        minx[:] = _int64_max
        nan_val = NPY_NAT
    elif groupby_t is uint64_t:
        # NB: We do not define nan_val because there is no such thing
        #  for uint64_t.  We carefully avoid having to reference it in this
        #  case.
        minx[:] = np.iinfo(np.uint64).max
    else:
        minx[:] = np.inf
        nan_val = NAN

    N, K = (<object>values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                if not _treat_as_na(val, True):
                    # TODO: Sure we always want is_datetimelike=True?
                    nobs[lab, j] += 1
                    if val < minx[lab, j]:
                        minx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    if groupby_t is uint64_t:
                        runtime_error = True
                        break
                    else:
                        out[i, j] = nan_val
                else:
                    out[i, j] = minx[i, j]

    if runtime_error:
        # We cannot raise directly above because that is within a nogil
        #  block.
        raise RuntimeError("empty group with uint64_t")


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cummin(groupby_t[:, :] out,
                 groupby_t[:, :] values,
                 const int64_t[:] labels,
                 int ngroups,
                 bint is_datetimelike):
    """
    Cumulative minimum of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : array
        Array to store cummin in.
    values : array
        Values to take cummin of.
    labels : int64 array
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        groupby_t val, mval
        ndarray[groupby_t, ndim=2] accum
        int64_t lab

    N, K = (<object>values).shape
    accum = np.empty((ngroups, K), dtype=np.asarray(values).dtype)
    if groupby_t is int64_t:
        accum[:] = _int64_max
    elif groupby_t is uint64_t:
        accum[:] = np.iinfo(np.uint64).max
    else:
        accum[:] = np.inf

    with nogil:
        for i in range(N):
            lab = labels[i]

            if lab < 0:
                continue
            for j in range(K):
                val = values[i, j]

                if _treat_as_na(val, is_datetimelike):
                    out[i, j] = val
                else:
                    mval = accum[lab, j]
                    if val < mval:
                        accum[lab, j] = mval = val
                    out[i, j] = mval


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cummax(groupby_t[:, :] out,
                 groupby_t[:, :] values,
                 const int64_t[:] labels,
                 int ngroups,
                 bint is_datetimelike):
    """
    Cumulative maximum of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : array
        Array to store cummax in.
    values : array
        Values to take cummax of.
    labels : int64 array
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        groupby_t val, mval
        ndarray[groupby_t, ndim=2] accum
        int64_t lab

    N, K = (<object>values).shape
    accum = np.empty((ngroups, K), dtype=np.asarray(values).dtype)
    if groupby_t is int64_t:
        accum[:] = -_int64_max
    elif groupby_t is uint64_t:
        accum[:] = 0
    else:
        accum[:] = -np.inf

    with nogil:
        for i in range(N):
            lab = labels[i]

            if lab < 0:
                continue
            for j in range(K):
                val = values[i, j]

                if _treat_as_na(val, is_datetimelike):
                    out[i, j] = val
                else:
                    mval = accum[lab, j]
                    if val > mval:
                        accum[lab, j] = mval = val
                    out[i, j] = mval
