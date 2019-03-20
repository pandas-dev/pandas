# -*- coding: utf-8 -*-

import cython
from cython import Py_ssize_t
from cython cimport floating

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t)
cnp.import_array()


from pandas._libs.util cimport numeric, get_nat

from pandas._libs.algos cimport (swap, TiebreakEnumType, TIEBREAK_AVERAGE,
                                 TIEBREAK_MIN, TIEBREAK_MAX, TIEBREAK_FIRST,
                                 TIEBREAK_DENSE)
from pandas._libs.algos import (take_2d_axis1_float64_float64,
                                groupsort_indexer, tiebreakers)

cdef int64_t NPY_NAT = get_nat()

cdef float64_t NaN = <float64_t>np.NaN


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
                          bint is_datetimelike,
                          bint skipna=True):
    """
    Only transforms on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        float64_t val
        float64_t[:, :] accum
        int64_t lab

    N, K = (<object>values).shape
    accum = np.ones_like(values)

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
                 is_datetimelike,
                 bint skipna=True):
    """
    Only transforms on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        numeric val
        numeric[:, :] accum
        int64_t lab

    N, K = (<object>values).shape
    accum = np.zeros_like(values)

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
        int offset, sign
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
    """Indexes how to fill values forwards or backwards within a group

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
    """Aggregated boolean values to show truthfulness of group elements

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


@cython.wraparound(False)
@cython.boundscheck(False)
def _group_add(floating[:, :] out,
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
        floating[:, :] sumx, nobs

    if len(values) != len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
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
                if nobs[i, j] < min_count:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j]


group_add_float32 = _group_add['float']
group_add_float64 = _group_add['double']


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
        floating[:, :] prodx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
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
        floating[:, :] nobs, mean

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
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
        floating[:, :] sumx, nobs

    assert min_count == -1, "'min_count' only used in add and prod"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
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
                # Casting to int will intentionaly truncate result
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


# generated from template
include "groupby_helper.pxi"
