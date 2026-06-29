cimport cython
from cython cimport (
    Py_ssize_t,
    floating,
)
from libc.math cimport (
    NAN,
    isfinite,
    sqrt,
)
from libc.stdlib cimport (
    free,
    malloc,
)

import numpy as np

cimport numpy as cnp
from numpy cimport (
    complex64_t,
    complex128_t,
    float32_t,
    float64_t,
    int8_t,
    int64_t,
    intp_t,
    ndarray,
    uint8_t,
    uint64_t,
)

cnp.import_array()

from pandas._libs cimport util
from pandas._libs.algos cimport (
    get_rank_nan_fill_val,
    kth_smallest_c,
)

from pandas._libs.algos import (
    groupsort_indexer,
    rank_1d,
    take_2d_axis1_bool_bool,
    take_2d_axis1_float64_float64,
)

from pandas._libs.dtypes cimport (
    numeric_object_t,
    numeric_t,
)
from pandas._libs.missing cimport checknull


cdef int64_t NPY_NAT = util.get_nat()

cdef float64_t NaN = <float64_t>np.nan

cdef enum InterpolationEnumType:
    INTERPOLATION_LINEAR,
    INTERPOLATION_LOWER,
    INTERPOLATION_HIGHER,
    INTERPOLATION_NEAREST,
    INTERPOLATION_MIDPOINT


cdef float64_t median_linear_mask(
    float64_t* a,
    int n,
    uint8_t* mask,
    bint skipna=True
) noexcept nogil:
    cdef:
        int i, j, na_count = 0
        float64_t* tmp
        float64_t result

    if n == 0:
        return NaN

    # count NAs
    for i in range(n):
        if mask[i]:
            na_count += 1

    if na_count:
        if na_count == n or not skipna:
            return NaN

        tmp = <float64_t*>malloc((n - na_count) * sizeof(float64_t))
        if tmp is NULL:
            raise MemoryError()

        j = 0
        for i in range(n):
            if not mask[i]:
                tmp[j] = a[i]
                j += 1

        a = tmp
        n -= na_count

    result = calc_median_linear(a, n)

    if na_count:
        free(a)

    return result


cdef float64_t median_linear(
    float64_t* a,
    int n,
    bint is_datetimelike=False,
    bint skipna=True,
) noexcept nogil:
    cdef:
        int i, j, na_count = 0
        float64_t* tmp
        float64_t result

    if n == 0:
        return NaN

    # count NAs
    if is_datetimelike:
        for i in range(n):
            if a[i] == NPY_NAT:
                na_count += 1
    else:
        for i in range(n):
            if a[i] != a[i]:
                na_count += 1

    if na_count:
        if na_count == n or not skipna:
            return NaN

        tmp = <float64_t*>malloc((n - na_count) * sizeof(float64_t))
        if tmp is NULL:
            raise MemoryError()

        j = 0
        if is_datetimelike:
            for i in range(n):
                if a[i] != NPY_NAT:
                    tmp[j] = a[i]
                    j += 1
        else:
            for i in range(n):
                if a[i] == a[i]:
                    tmp[j] = a[i]
                    j += 1

        a = tmp
        n -= na_count

    result = calc_median_linear(a, n)

    if na_count:
        free(a)

    return result


cdef float64_t calc_median_linear(float64_t* a, int n) noexcept nogil:
    cdef:
        float64_t result

    if n % 2:
        result = kth_smallest_c(a, n // 2, n)
    else:
        result = (kth_smallest_c(a, n // 2, n) +
                  kth_smallest_c(a, n // 2 - 1, n)) / 2

    return result


ctypedef fused int64float_t:
    int64_t
    uint64_t
    float32_t
    float64_t


@cython.boundscheck(False)
@cython.wraparound(False)
def group_median_float64(
    ndarray[float64_t, ndim=2] out,
    ndarray[int64_t] counts,
    ndarray[float64_t, ndim=2] values,
    ndarray[intp_t] labels,
    Py_ssize_t min_count=-1,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint is_datetimelike=False,
    bint skipna=True,
) -> None:
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, ngroups, size
        ndarray[intp_t] _counts
        ndarray[float64_t, ndim=2] data
        ndarray[uint8_t, ndim=2] data_mask
        ndarray[intp_t] indexer
        float64_t* ptr
        uint8_t* ptr_mask
        float64_t result
        bint uses_mask = mask is not None

    assert min_count == -1, "'min_count' only used in sum and prod"

    ngroups = len(counts)
    N, K = (<object>values).shape

    indexer, _counts = groupsort_indexer(labels, ngroups)
    counts[:] = _counts[1:]

    data = np.empty((K, N), dtype=np.float64)
    ptr = <float64_t*>cnp.PyArray_DATA(data)

    take_2d_axis1_float64_float64(values.T, indexer, out=data)

    if uses_mask:
        data_mask = np.empty((K, N), dtype=np.uint8)
        ptr_mask = <uint8_t *>cnp.PyArray_DATA(data_mask)

        take_2d_axis1_bool_bool(mask.T, indexer, out=data_mask, fill_value=1)

        with nogil:

            for i in range(K):
                # exclude NA group
                ptr += _counts[0]
                ptr_mask += _counts[0]

                for j in range(ngroups):
                    size = _counts[j + 1]
                    result = median_linear_mask(ptr, size, ptr_mask, skipna)
                    out[j, i] = result

                    if result != result:
                        result_mask[j, i] = 1
                    ptr += size
                    ptr_mask += size

    else:
        with nogil:
            for i in range(K):
                # exclude NA group
                ptr += _counts[0]
                for j in range(ngroups):
                    size = _counts[j + 1]
                    out[j, i] = median_linear(ptr, size, is_datetimelike, skipna)
                    ptr += size


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def group_cumprod(
    int64float_t[:, ::1] out,
    const int64float_t[:, :] values,
    const intp_t[::1] labels,
    int ngroups,
    bint is_datetimelike,
    bint skipna=True,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
) -> None:
    """
    Cumulative product of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : np.ndarray[np.float64, ndim=2]
        Array to store cumprod in.
    values : np.ndarray[np.float64, ndim=2]
        Values to take cumprod of.
    labels : np.ndarray[np.intp]
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        Always false, `values` is never datetime-like.
    skipna : bool
        If true, ignore nans in `values`.
    mask : np.ndarray[uint8], optional
        Mask of values
    result_mask : np.ndarray[int8], optional
        Mask of out array

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        Py_ssize_t i, j, N, K
        int64float_t val, na_val
        int64float_t[:, ::1] accum
        intp_t lab
        uint8_t[:, ::1] accum_mask
        bint uses_mask = mask is not None

        int64float_t *out_row
        uint8_t *res_mask_row
        int64float_t *accum_row
        uint8_t *acc_mask_row

    N = values.shape[0]
    K = values.shape[1]

    accum_arr = np.ones((ngroups, K), dtype=np.asarray(out).dtype)
    accum = accum_arr

    na_val = _get_na_val(<int64float_t>0, is_datetimelike)

    accum_mask_arr = np.zeros((ngroups, K), dtype=np.uint8)
    accum_mask = accum_mask_arr

    with nogil:
        if uses_mask:
            for i in range(N):

                lab = labels[i]
                if lab < 0:
                    continue

                out_row = &out[i, 0]
                res_mask_row = &result_mask[i, 0]
                accum_row = &accum[lab, 0]
                acc_mask_row = &accum_mask[lab, 0]

                for j in range(K):
                    if mask[i, j]:
                        res_mask_row[j] = 1
                        out_row[j] = 0
                        if not skipna:
                            accum_row[j] = na_val
                            acc_mask_row[j] = 1
                    else:
                        if acc_mask_row[j]:
                            out_row[j] = na_val
                            res_mask_row[j] = 1
                        else:
                            accum_row[j] *= values[i, j]
                            out_row[j] = accum_row[j]
        else:
            for i in range(N):

                lab = labels[i]
                if lab < 0:
                    continue

                out_row = &out[i, 0]
                accum_row = &accum[lab, 0]
                acc_mask_row = &accum_mask[lab, 0]

                for j in range(K):
                    val = values[i, j]

                    if _treat_as_na(val, False):
                        out_row[j] = na_val
                        if not skipna:
                            accum_row[j] = na_val
                            acc_mask_row[j] = 1
                    else:
                        if acc_mask_row[j]:
                            out_row[j] = na_val
                        else:
                            accum_row[j] *= val
                            out_row[j] = accum_row[j]


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cumsum(
    int64float_t[:, ::1] out,
    ndarray[int64float_t, ndim=2] values,
    const intp_t[::1] labels,
    int ngroups,
    bint is_datetimelike,
    bint skipna=True,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
) -> None:
    """
    Cumulative sum of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : np.ndarray[ndim=2]
        Array to store cumsum in.
    values : np.ndarray[ndim=2]
        Values to take cumsum of.
    labels : np.ndarray[np.intp]
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    skipna : bool
        If true, ignore nans in `values`.
    mask : np.ndarray[uint8], optional
        Mask of values
    result_mask : np.ndarray[int8], optional
        Mask of out array

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        Py_ssize_t i, j, N, K
        int64float_t val, y, t, na_val
        int64float_t[:, ::1] accum, compensation
        uint8_t[:, ::1] accum_mask
        intp_t lab
        bint isna_entry, isna_prev = False
        bint uses_mask = mask is not None

    N, K = (<object>values).shape

    if uses_mask:
        accum_mask = np.zeros((ngroups, K), dtype="uint8")

    accum = np.zeros((ngroups, K), dtype=np.asarray(values).dtype)
    compensation = np.zeros((ngroups, K), dtype=np.asarray(values).dtype)

    na_val = _get_na_val(<int64float_t>0, is_datetimelike)

    with nogil:
        for i in range(N):
            lab = labels[i]

            if uses_mask and lab < 0:
                # GH#58811
                result_mask[i, :] = True
                out[i, :] = 0
                continue
            elif lab < 0:
                continue

            for j in range(K):
                val = values[i, j]

                if uses_mask:
                    isna_entry = mask[i, j]
                else:
                    isna_entry = _treat_as_na(val, is_datetimelike)

                if not skipna:
                    if uses_mask:
                        isna_prev = accum_mask[lab, j]
                    else:
                        isna_prev = _treat_as_na(accum[lab, j], is_datetimelike)

                    if isna_prev:
                        if uses_mask:
                            result_mask[i, j] = True
                            # Be deterministic, out was initialized as empty
                            out[i, j] = 0
                        else:
                            out[i, j] = na_val
                        continue

                if isna_entry:

                    if uses_mask:
                        result_mask[i, j] = True
                        # Be deterministic, out was initialized as empty
                        out[i, j] = 0
                    else:
                        out[i, j] = na_val

                    if not skipna:
                        if uses_mask:
                            accum_mask[lab, j] = True
                        else:
                            accum[lab, j] = na_val

                else:
                    # For floats, use Kahan summation to reduce floating-point
                    # error (https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
                    if int64float_t == float32_t or int64float_t == float64_t:
                        y = val - compensation[lab, j]
                        t = accum[lab, j] + y
                        compensation[lab, j] = t - accum[lab, j] - y
                    else:
                        t = val + accum[lab, j]

                    accum[lab, j] = t
                    out[i, j] = t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def group_shift_indexer(
    int64_t[::1] out,
    const intp_t[::1] labels,
    int ngroups,
    int periods,
) -> None:
    cdef:
        Py_ssize_t N, i, ii, lab
        int64_t idxer, idxer_slot, seen_count

        int64_t[::1] label_seen
        int64_t[:, ::1] label_indexer

        const intp_t *labels_ptr
        int64_t *out_ptr
        int64_t *seen_ptr
        int64_t *indexer_ptr


    N = labels.shape[0]

    if N == 0:
        return

    if periods == 0:
        with nogil:
            for i in range(N):
                out[i] = i
        return

    if ngroups == 0:
        with nogil:
            for i in range(N):
                out[i] = -1
        return

    label_seen_arr = np.zeros(ngroups, dtype=np.int64)
    label_seen = label_seen_arr

    label_indexer_arr = np.empty((ngroups, abs(periods)), dtype=np.int64)
    label_indexer = label_indexer_arr

    labels_ptr = &labels[0]
    out_ptr = &out[0]
    seen_ptr = &label_seen[0]
    indexer_ptr = &label_indexer[0, 0]

    # Split forward and reverse traversal outside the hot loop.
    if periods > 0:
        with nogil:
            for i in range(N):
                ii = i
                lab = labels_ptr[ii]

                if lab == -1:
                    out_ptr[ii] = -1
                    continue

                seen_count = seen_ptr[lab] + 1
                seen_ptr[lab] = seen_count

                idxer_slot = seen_count % periods

                if seen_count > periods:
                    idxer = indexer_ptr[lab * periods + idxer_slot]
                    out_ptr[ii] = idxer
                else:
                    out_ptr[ii] = -1

                indexer_ptr[lab * periods + idxer_slot] = ii
    else:
        periods = -periods
        with nogil:
            for i in range(N):
                ii = N - 1 - i
                lab = labels_ptr[ii]

                if lab == -1:
                    out_ptr[ii] = -1
                    continue

                seen_count = seen_ptr[lab] + 1
                seen_ptr[lab] = seen_count

                idxer_slot = seen_count % periods

                if seen_count > periods:
                    idxer = indexer_ptr[lab * periods + idxer_slot]
                    out_ptr[ii] = idxer
                else:
                    out_ptr[ii] = -1

                indexer_ptr[lab * periods + idxer_slot] = ii


@cython.wraparound(False)
@cython.boundscheck(False)
def group_fillna_indexer(
    Py_ssize_t[::1] out,
    const intp_t[::1] labels,
    const uint8_t[:] mask,
    int64_t limit,
    bint compute_ffill,
    int ngroups,
) -> None:
    """
    Indexes how to fill values forwards or backwards within a group.

    Parameters
    ----------
    out : np.ndarray[np.intp]
        Values into which this method will write its results.
    labels : np.ndarray[np.intp]
        Array containing unique label for each group, with its ordering
        matching up to the corresponding record in `values`.
    mask : np.ndarray[np.uint8]
        Indicating whether a value is na or not.
    limit : int64_t
        Consecutive values to fill before stopping, or -1 for no limit.
    compute_ffill : bint
        Whether to compute ffill or bfill.
    ngroups : int
        Number of groups, larger than all entries of `labels`.

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object
    """
    cdef:
        Py_ssize_t idx, N = len(out)
        intp_t label
        intp_t[::1] last = -1 * np.ones(ngroups, dtype=np.intp)
        intp_t[::1] fill_count = np.zeros(ngroups, dtype=np.intp)

    # Make sure all arrays are the same size
    assert N == len(labels) == len(mask)

    with nogil:
        # Can't use for loop with +/- step
        # https://github.com/cython/cython/issues/1106
        idx = 0 if compute_ffill else N-1
        for _ in range(N):
            label = labels[idx]
            if label == -1:  # na-group gets na-values
                out[idx] = -1
            elif mask[idx] == 1:  # is missing
                # Stop filling once we've hit the limit
                if limit != -1 and fill_count[label] >= limit:
                    out[idx] = -1
                else:
                    out[idx] = last[label]
                    fill_count[label] += 1
            else:
                fill_count[label] = 0  # reset items when not missing
                last[label] = idx
                out[idx] = idx

            if compute_ffill:
                idx += 1
            else:
                idx -= 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def group_any_all(
    int8_t[:, ::1] out,
    const int8_t[:, :] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask,
    str val_test,
    bint skipna,
    uint8_t[:, ::1] result_mask,
) -> None:
    """
    Aggregated boolean values to show truthfulness of group elements. If the
    input is a nullable type (result_mask is not None), the result will be computed
    using Kleene logic.

    Parameters
    ----------
    out : np.ndarray[np.int8]
        Values into which this method will write its results.
    labels : np.ndarray[np.intp]
        Array containing unique label for each group, with its
        ordering matching up to the corresponding record in `values`
    values : np.ndarray[np.int8]
        Containing the truth value of each element.
    mask : np.ndarray[np.uint8]
        Indicating whether a value is na or not.
    val_test : {'any', 'all'}
        String object dictating whether to use any or all truth testing
    skipna : bool
        Flag to ignore nan values during truth testing
    result_mask : ndarray[bool, ndim=2], optional
        If not None, these specify locations in the output that are NA.
        Modified in-place.

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object.
    The returned values will either be 0, 1 (False or True, respectively), or
    -1 to signify a masked position in the case of a nullable input.
    """
    cdef:
        Py_ssize_t i, j, N = len(labels), K = out.shape[1]
        intp_t lab
        int8_t flag_val, val
        bint uses_mask = result_mask is not None

        int8_t *out_row = NULL
        uint8_t *res_mask_row = NULL

    if val_test == "all":
        flag_val = 0
    elif val_test == "any":
        flag_val = 1
    else:
        raise ValueError("'val_test' must be either 'any' or 'all'!")

    out[:] = 1 - flag_val

    # Split invariant mask and skipna branches outside the inner loop.
    with nogil:
        if uses_mask:
            if skipna:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    out_row = &out[lab, 0]
                    res_mask_row = &result_mask[lab, 0]

                    for j in range(K):
                        if mask[i, j]:
                            continue
                        if values[i, j] == flag_val:
                            out_row[j] = flag_val
                            res_mask_row[j] = 0
            else:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    out_row = &out[lab, 0]
                    res_mask_row = &result_mask[lab, 0]

                    for j in range(K):
                        if mask[i, j]:
                            if out_row[j] != flag_val:
                                res_mask_row[j] = 1
                            continue

                        if values[i, j] == flag_val:
                            out_row[j] = flag_val
                            res_mask_row[j] = 0
        else:
            if skipna:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    out_row = &out[lab, 0]

                    for j in range(K):
                        if mask[i, j]:
                            continue
                        if values[i, j] == flag_val:
                            out_row[j] = flag_val
            else:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    out_row = &out[lab, 0]

                    for j in range(K):
                        if values[i, j] == flag_val:
                            out_row[j] = flag_val


# ----------------------------------------------------------------------
# group_sum, group_prod, group_var, group_mean, group_ohlc
# ----------------------------------------------------------------------

ctypedef fused mean_t:
    float64_t
    float32_t
    complex64_t
    complex128_t

ctypedef fused sum_t:
    mean_t
    int64_t
    uint64_t
    object


@cython.wraparound(False)
@cython.boundscheck(False)
def group_sum(
    sum_t[:, ::1] out,
    int64_t[::1] counts,
    ndarray[sum_t, ndim=2] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask,
    uint8_t[:, ::1] result_mask=None,
    Py_ssize_t min_count=0,
    bint is_datetimelike=False,
    object initial=0,
    bint skipna=True,
) -> None:
    """
    Only aggregates on axis=0 using Kahan summation
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        sum_t val, t, y, nan_val
        sum_t[:, ::1] sumx, compensation
        int64_t[:, ::1] nobs
        Py_ssize_t len_values = len(values), len_labels = len(labels)
        bint uses_mask = mask is not None
        bint isna_entry, isna_result

    if len_values != len_labels:
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    if initial == 0:
        # the below is equivalent to `np.zeros_like(out)` but faster
        sumx = np.zeros((<object>out).shape, dtype=(<object>out).base.dtype)
        compensation = np.zeros((<object>out).shape, dtype=(<object>out).base.dtype)
    else:
        # in practice this path is only taken for strings to use empty string as initial
        assert sum_t is object
        sumx = np.full((<object>out).shape, initial, dtype=object)
        # object code path does not use `compensation`

    N, K = (<object>values).shape
    if uses_mask:
        nan_val = 0
    elif is_datetimelike:
        nan_val = NPY_NAT
    elif sum_t is int64_t or sum_t is uint64_t:
        # This has no effect as int64 can't be nan. Setting to 0 to avoid type error
        nan_val = 0
    else:
        nan_val = NAN

    with nogil(sum_t is not object):
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1

            for j in range(K):
                val = values[i, j]

                if uses_mask:
                    isna_entry = mask[i, j]
                else:
                    isna_entry = _treat_as_na(val, is_datetimelike)

                if not skipna:
                    if uses_mask:
                        isna_result = result_mask[lab, j]
                    else:
                        isna_result = _treat_as_na(sumx[lab, j], is_datetimelike)

                    if isna_result:
                        # If sum is already NA, don't add to it. This is important for
                        # datetimelikebecause adding a value to NPY_NAT may not result
                        # in a NPY_NAT
                        continue

                if not isna_entry:
                    nobs[lab, j] += 1

                    if sum_t is object or sum_t is int64_t or sum_t is uint64_t:
                        # NB: this does not use 'compensation' like the non-object
                        #  and non-integer track does.
                        if nobs[lab, j] == 1:
                            # i.e. we haven't added anything yet; avoid TypeError
                            #  if e.g. val is a str and sumx[lab, j] is 0
                            t = val
                        else:
                            t = sumx[lab, j] + val
                        sumx[lab, j] = t

                    else:
                        y = val - compensation[lab, j]
                        t = sumx[lab, j] + y
                        compensation[lab, j] = t - sumx[lab, j] - y

                        # Handle float overflow
                        if (
                            sum_t is float32_t or sum_t is float64_t
                        ) and not isfinite(compensation[lab, j]):
                            # GH#53606; GH#60303
                            # If val is +/- infinity compensation is NaN
                            # which would lead to results being NaN instead
                            # of +/- infinity. We cannot use util.is_nan
                            # because of no gil
                            compensation[lab, j] = 0

                        # Handle complex overflow
                        if (
                            sum_t is complex64_t or sum_t is complex128_t
                        ) and not isfinite(compensation[lab, j].real):
                            compensation[lab, j].real = 0

                        if (
                            sum_t is complex64_t or sum_t is complex128_t
                        ) and not isfinite(compensation[lab, j].imag):
                            compensation[lab, j].imag = 0

                        sumx[lab, j] = t
                elif not skipna:
                    if uses_mask:
                        result_mask[lab, j] = True
                    else:
                        sumx[lab, j] = nan_val

    _check_below_mincount(
        out, uses_mask, result_mask, ncounts, K, nobs, min_count, sumx
    )


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
def group_prod(
    int64float_t[:, ::1] out,
    int64_t[::1] counts,
    const int64float_t[:, :] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    Py_ssize_t min_count=0,
    bint skipna=True,
) -> None:
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        int64float_t val, nan_val
        int64float_t[:, ::1] prodx
        int64_t[:, ::1] nobs
        bint uses_mask = mask is not None

        int64float_t *prodx_row
        int64_t *nobs_row

    if values.shape[0] != len(labels):
        raise ValueError("len(index) != len(labels)")

    N = values.shape[0]
    K = values.shape[1]

    nobs_arr = np.zeros((out.shape[0], K), dtype=np.int64)
    nobs = nobs_arr

    prodx_arr = np.ones((out.shape[0], K), dtype=np.asarray(out).dtype)
    prodx = prodx_arr
    nan_val = _get_na_val(<int64float_t>0, False)

    with nogil:
        if uses_mask:
            for i in range(N):

                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1

                prodx_row = &prodx[lab, 0]
                nobs_row = &nobs[lab, 0]

                for j in range(K):
                    if not mask[i, j]:
                        nobs_row[j] += 1
                        prodx_row[j] *= values[i, j]
                    elif not skipna:
                        result_mask[lab, j] = True

        else:
            for i in range(N):

                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1

                prodx_row = &prodx[lab, 0]
                nobs_row = &nobs[lab, 0]

                for j in range(K):
                    val = values[i, j]
                    if not _treat_as_na(val, False):
                        nobs_row[j] += 1
                        prodx_row[j] *= val
                    elif not skipna:
                        prodx_row[j] = nan_val

    _check_below_mincount(
        out, uses_mask, result_mask, ncounts, K, nobs, min_count, prodx
    )


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def group_var(
    floating[:, ::1] out,
    int64_t[::1] counts,
    ndarray[floating, ndim=2] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    int64_t ddof=1,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint is_datetimelike=False,
    str name="var",
    bint skipna=True,
) -> None:
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        floating val, ct, oldmean
        floating[:, ::1] mean
        int64_t[:, ::1] nobs
        Py_ssize_t len_values = len(values), len_labels = len(labels)
        bint isna_entry, isna_result, uses_mask = mask is not None
        bint is_std = name == "std"
        bint is_sem = name == "sem"

    assert min_count == -1, "'min_count' only used in sum and prod"

    if len_values != len_labels:
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    mean = np.zeros((<object>out).shape, dtype=(<object>out).base.dtype)

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

                if uses_mask:
                    isna_entry = mask[i, j]
                elif is_datetimelike:
                    # With group_var, we cannot just use _treat_as_na bc
                    #  datetimelike dtypes get cast to float64 instead of
                    #  to int64.
                    isna_entry = val == NPY_NAT
                else:
                    isna_entry = _treat_as_na(val, is_datetimelike)

                if not skipna:
                    if uses_mask:
                        isna_result = result_mask[lab, j]
                    elif is_datetimelike:
                        # With group_var, we cannot just use _treat_as_na bc
                        #  datetimelike dtypes get cast to float64 instead of
                        #  to int64.
                        isna_result = out[lab, j] == NPY_NAT
                    else:
                        isna_result = _treat_as_na(out[lab, j], is_datetimelike)

                    if isna_result:
                        # If aggregate is already NA, don't add to it. This is
                        # important for datetimelike because adding a value to NPY_NAT
                        # may not result in a NPY_NAT
                        continue

                if not isna_entry:
                    nobs[lab, j] += 1
                    oldmean = mean[lab, j]
                    mean[lab, j] += (val - oldmean) / nobs[lab, j]
                    out[lab, j] += (val - mean[lab, j]) * (val - oldmean)
                elif not skipna:
                    nobs[lab, j] = 0
                    if uses_mask:
                        result_mask[lab, j] = True
                    else:
                        out[lab, j] = NAN

        for i in range(ncounts):
            for j in range(K):
                ct = nobs[i, j]
                if ct <= ddof:
                    if uses_mask:
                        result_mask[i, j] = True
                    else:
                        out[i, j] = NAN
                else:
                    if is_std:
                        out[i, j] = sqrt(out[i, j] / (ct - ddof))
                    elif is_sem:
                        out[i, j] = sqrt(out[i, j] / (ct - ddof) / ct)
                    else:
                        # just "var"
                        out[i, j] /= (ct - ddof)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.cpow(True)
def group_skew(
    float64_t[:, ::1] out,
    int64_t[::1] counts,
    ndarray[float64_t, ndim=2] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    cdef:
        Py_ssize_t i, j, N, K, lab, ngroups = len(counts)
        int64_t[:, ::1] nobs
        Py_ssize_t len_values = len(values), len_labels = len(labels)
        bint isna_entry, uses_mask = mask is not None
        float64_t[:, ::1] M1, M2, M3
        float64_t delta, delta_n, term1, val
        int64_t n1, n
        float64_t ct

    if len_values != len_labels:
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)

    # M1, M2, and M3 correspond to 1st, 2nd, and third Moments
    M1 = np.zeros((<object>out).shape, dtype=np.float64)
    M2 = np.zeros((<object>out).shape, dtype=np.float64)
    M3 = np.zeros((<object>out).shape, dtype=np.float64)

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

                if uses_mask:
                    isna_entry = mask[i, j]
                else:
                    isna_entry = _treat_as_na(val, False)

                if not isna_entry:
                    # Running stats update based on RunningStats::Push from
                    #  https://www.johndcook.com/blog/skewness_kurtosis/
                    n1 = nobs[lab, j]
                    n = n1 + 1

                    nobs[lab, j] = n
                    delta = val - M1[lab, j]
                    delta_n = delta / n
                    term1 = delta * delta_n * n1

                    M1[lab, j] += delta_n
                    M3[lab, j] += term1 * delta_n * (n - 2) - 3 * delta_n * M2[lab, j]
                    M2[lab, j] += term1
                elif not skipna:
                    M1[lab, j] = NaN
                    M2[lab, j] = NaN
                    M3[lab, j] = NaN

        for i in range(ngroups):
            for j in range(K):
                ct = <float64_t>nobs[i, j]
                if ct < 3:
                    if result_mask is not None:
                        result_mask[i, j] = 1
                    out[i, j] = NaN
                elif M2[i, j] == 0:
                    out[i, j] = 0
                else:
                    out[i, j] = (
                        (ct * (ct - 1) ** 0.5 / (ct - 2))
                        * (M3[i, j] / M2[i, j] ** 1.5)
                    )


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.cpow(True)
def group_kurt(
    float64_t[:, ::1] out,
    int64_t[::1] counts,
    ndarray[float64_t, ndim=2] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    cdef:
        Py_ssize_t i, j, N, K, lab, ngroups = len(counts)
        int64_t[:, ::1] nobs
        Py_ssize_t len_values = len(values), len_labels = len(labels)
        bint isna_entry, uses_mask = mask is not None
        float64_t[:, ::1] M1, M2, M3, M4
        float64_t delta, delta_n, delta_n2, term1, val
        int64_t n1, n
        float64_t ct, num, den, adj

    if len_values != len_labels:
        raise ValueError("len(index) != len(labels)")

    nobs = np.zeros((<object>out).shape, dtype=np.int64)

    # M1, M2, M3 and M4 correspond to 1st, 2nd, 3rd and 4th Moments
    M1 = np.zeros((<object>out).shape, dtype=np.float64)
    M2 = np.zeros((<object>out).shape, dtype=np.float64)
    M3 = np.zeros((<object>out).shape, dtype=np.float64)
    M4 = np.zeros((<object>out).shape, dtype=np.float64)

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

                if uses_mask:
                    isna_entry = mask[i, j]
                else:
                    isna_entry = _treat_as_na(val, False)

                if not isna_entry:
                    # Running stats update based on RunningStats::Push from
                    #  https://www.johndcook.com/blog/skewness_kurtosis/
                    n1 = nobs[lab, j]
                    n = n1 + 1

                    nobs[lab, j] = n
                    delta = val - M1[lab, j]
                    delta_n = delta / n
                    delta_n2 = delta_n * delta_n
                    term1 = delta * delta_n * n1

                    M1[lab, j] += delta_n
                    M4[lab, j] += (term1 * delta_n2 * (n*n - 3*n + 3)
                                   + 6 * delta_n2 * M2[lab, j]
                                   - 4 * delta_n * M3[lab, j])
                    M3[lab, j] += term1 * delta_n * (n - 2) - 3 * delta_n * M2[lab, j]
                    M2[lab, j] += term1
                elif not skipna:
                    M1[lab, j] = NaN
                    M2[lab, j] = NaN
                    M3[lab, j] = NaN
                    M4[lab, j] = NaN

        for i in range(ngroups):
            for j in range(K):
                ct = <float64_t>nobs[i, j]
                if ct < 4:
                    if result_mask is not None:
                        result_mask[i, j] = 1
                    out[i, j] = NaN
                elif M2[i, j] == 0:
                    out[i, j] = 0
                else:
                    num = ct * (ct + 1) * (ct - 1) * M4[i, j]
                    den = (ct - 2) * (ct - 3) * M2[i, j] ** 2
                    adj = 3.0 * (ct - 1) ** 2 / ((ct - 2) * (ct - 3))
                    out[i, j] = num / den - adj


@cython.wraparound(False)
@cython.boundscheck(False)
def group_mean(
    mean_t[:, ::1] out,
    int64_t[::1] counts,
    const mean_t[:, :] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    bint is_datetimelike=False,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    """
    Compute the mean per label given a label assignment for each value.
    NaN values are ignored.

    Parameters
    ----------
    out : np.ndarray[floating]
        Values into which this method will write its results.
    counts : np.ndarray[int64]
        A zeroed array of the same shape as labels,
        populated by group sizes during algorithm.
    values : np.ndarray[floating]
        2-d array of the values to find the mean of.
    labels : np.ndarray[np.intp]
        Array containing unique label for each group, with its
        ordering matching up to the corresponding record in `values`.
    min_count : Py_ssize_t
        Only used in sum and prod. Always -1.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    mask : ndarray[bool, ndim=2], optional
        Mask of the input values.
    result_mask : ndarray[bool, ndim=2], optional
        Mask of the out array
    skipna : bool, optional
        If True, ignore nans in `values`.

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object.
    `counts` is modified to hold group sizes
    """

    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        mean_t val, count, y, t, nan_val, comp_val
        mean_t[:, ::1] sumx, compensation
        int64_t[:, ::1] nobs
        Py_ssize_t len_values = values.shape[0]
        Py_ssize_t len_labels = labels.shape[0]
        bint isna_entry, isna_result, uses_mask = mask is not None

        mean_t *sum_row
        mean_t *comp_row
        int64_t *nobs_row
        mean_t *out_row
        uint8_t *res_mask_row

    assert min_count == -1, "'min_count' only used in sum and prod"

    if len_values != len_labels:
        raise ValueError("len(index) != len(labels)")

    N = values.shape[0]
    K = values.shape[1]

    dt = np.asarray(out).dtype
    nobs_arr = np.zeros((out.shape[0], K), dtype=np.int64)
    nobs = nobs_arr
    sumx_arr = np.zeros((out.shape[0], K), dtype=dt)
    sumx = sumx_arr
    comp_arr = np.zeros((out.shape[0], K), dtype=dt)
    compensation = comp_arr

    if uses_mask:
        nan_val = 0
    elif is_datetimelike:
        nan_val = NPY_NAT
    else:
        nan_val = NAN

    with nogil:
        if uses_mask:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                nobs_row = &nobs[lab, 0]
                sum_row = &sumx[lab, 0]
                comp_row = &compensation[lab, 0]

                for j in range(K):
                    if not mask[i, j]:
                        nobs_row[j] += 1
                        y = values[i, j] - comp_row[j]
                        t = sum_row[j] + y
                        comp_val = t - sum_row[j] - y
                        if comp_val != comp_val:
                            comp_val = 0.
                        comp_row[j] = comp_val
                        sum_row[j] = t
                    elif not skipna:
                        # Keep datetimelike division from changing an NA result.
                        nobs_row[j] = 0
                        result_mask[lab, j] = True

        elif is_datetimelike:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                nobs_row = &nobs[lab, 0]
                sum_row = &sumx[lab, 0]
                comp_row = &compensation[lab, 0]

                for j in range(K):
                    val = values[i, j]
                    isna_entry = val == NPY_NAT

                    if not skipna:
                        isna_result = sum_row[j] == NPY_NAT
                        if isna_result:
                            continue

                    if not isna_entry:
                        nobs_row[j] += 1
                        y = val - comp_row[j]
                        t = sum_row[j] + y
                        comp_val = t - sum_row[j] - y
                        if comp_val != comp_val:
                            comp_val = 0.
                        comp_row[j] = comp_val
                        sum_row[j] = t
                    elif not skipna:
                        nobs_row[j] = 0
                        sum_row[j] = nan_val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                nobs_row = &nobs[lab, 0]
                sum_row = &sumx[lab, 0]
                comp_row = &compensation[lab, 0]

                for j in range(K):
                    val = values[i, j]
                    isna_entry = _treat_as_na(val, False)

                    if not skipna:
                        isna_result = _treat_as_na(sum_row[j], False)
                        if isna_result:
                            continue

                    if not isna_entry:
                        nobs_row[j] += 1
                        y = val - comp_row[j]
                        t = sum_row[j] + y
                        comp_val = t - sum_row[j] - y
                        if comp_val != comp_val:
                            comp_val = 0.
                        comp_row[j] = comp_val
                        sum_row[j] = t
                    elif not skipna:
                        nobs_row[j] = 0
                        sum_row[j] = nan_val

        if uses_mask:
            for i in range(ncounts):
                nobs_row = &nobs[i, 0]
                sum_row = &sumx[i, 0]
                out_row = &out[i, 0]
                res_mask_row = &result_mask[i, 0]
                for j in range(K):
                    count = nobs_row[j]
                    if count == 0:
                        res_mask_row[j] = 1
                    else:
                        out_row[j] = sum_row[j] / count
        else:
            for i in range(ncounts):
                nobs_row = &nobs[i, 0]
                sum_row = &sumx[i, 0]
                out_row = &out[i, 0]
                for j in range(K):
                    count = nobs_row[j]
                    if count == 0:
                        out_row[j] = nan_val
                    else:
                        out_row[j] = sum_row[j] / count
@cython.wraparound(False)
@cython.boundscheck(False)
def group_ohlc(
    int64float_t[:, ::1] out,
    int64_t[::1] counts,
    ndarray[int64float_t, ndim=2] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
) -> None:
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, N, K, lab
        int64float_t val
        uint8_t[::1] first_element_set
        bint isna_entry, uses_mask = mask is not None

    assert min_count == -1, "'min_count' only used in sum and prod"

    if len(labels) == 0:
        return

    N, K = (<object>values).shape

    if out.shape[1] != 4:
        raise ValueError("Output array must have 4 columns")

    if K > 1:
        raise NotImplementedError("Argument 'values' must have only one dimension")

    if int64float_t is float32_t or int64float_t is float64_t:
        out[:] = NAN
    else:
        out[:] = 0

    first_element_set = np.zeros((<object>counts).shape, dtype=np.uint8)
    if uses_mask:
        result_mask[:] = True

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab == -1:
                continue

            counts[lab] += 1
            val = values[i, 0]

            if uses_mask:
                isna_entry = mask[i, 0]
            else:
                isna_entry = _treat_as_na(val, False)

            if isna_entry:
                continue

            if not first_element_set[lab]:
                out[lab, 0] = out[lab, 1] = out[lab, 2] = out[lab, 3] = val
                first_element_set[lab] = True
                if uses_mask:
                    result_mask[lab] = False
            else:
                out[lab, 1] = max(out[lab, 1], val)
                out[lab, 2] = min(out[lab, 2], val)
                out[lab, 3] = val


@cython.boundscheck(False)
@cython.wraparound(False)
def group_quantile(
    float64_t[:, ::1] out,
    ndarray[numeric_t, ndim=1] values,
    const intp_t[::1] labels,
    const uint8_t[:] mask,
    const float64_t[:] qs,
    const int64_t[::1] starts,
    const int64_t[::1] ends,
    str interpolation,
    uint8_t[:, ::1] result_mask,
    bint is_datetimelike,
) -> None:
    """
    Calculate the quantile per group.

    Parameters
    ----------
    out : np.ndarray[np.float64, ndim=2]
        Array of aggregated values that will be written to.
    values : np.ndarray
        Array containing the values to apply the function against.
    labels : ndarray[np.intp]
        Array containing the unique group labels.
    qs : ndarray[float64_t]
        The quantile values to search for.
    starts : ndarray[int64]
        Positions at which each group begins.
    ends : ndarray[int64]
        Positions at which each group ends.
    interpolation : {'linear', 'lower', 'highest', 'nearest', 'midpoint'}
    result_mask : ndarray[bool, ndim=2] or None
    is_datetimelike : bool
        Whether int64 values represent datetime64-like values.

    Notes
    -----
    Rather than explicitly returning a value, this function modifies the
    provided `out` parameter.

    Uses kth_smallest_c (an O(n) quickselect) rather than a full O(n log n)
    argsort. NAs are pre-filtered into a temporary buffer so datetimelike
    NaT-ordering and float NaN comparisons are never an issue.
    """
    cdef:
        Py_ssize_t i, N = len(labels), ngroups, non_na_sz, k, nqs
        Py_ssize_t idx = 0
        Py_ssize_t grp_size
        InterpolationEnumType interp
        float64_t q_val, q_idx, frac, val, next_val
        bint uses_result_mask = result_mask is not None
        Py_ssize_t start, end
        numeric_t* tmp
        Py_ssize_t j

    assert values.shape[0] == N
    assert starts is not None
    assert ends is not None
    assert len(starts) == len(ends)

    if any(not (0 <= q <= 1) for q in qs):
        wrong = [x for x in qs if not (0 <= x <= 1)][0]
        raise ValueError(
            f"Each 'q' must be between 0 and 1. Got '{wrong}' instead"
        )

    inter_methods = {
        "linear": INTERPOLATION_LINEAR,
        "lower": INTERPOLATION_LOWER,
        "higher": INTERPOLATION_HIGHER,
        "nearest": INTERPOLATION_NEAREST,
        "midpoint": INTERPOLATION_MIDPOINT,
    }
    interp = inter_methods[interpolation]

    nqs = len(qs)
    ngroups = len(out)

    with nogil:
        for i in range(ngroups):
            start = starts[i]
            end = ends[i]

            # Count non-NA elements in this group using direct indexing
            # (avoids memoryview slicing, which is not nogil-safe).
            grp_size = end - start
            non_na_sz = 0
            for j in range(grp_size):
                if mask[start + j] == 0:
                    non_na_sz += 1

            if non_na_sz == 0:
                for k in range(nqs):
                    if uses_result_mask:
                        result_mask[i, k] = 1
                    else:
                        out[i, k] = NaN
            else:
                # Copy non-NA values into a temporary mutable buffer.
                # Pre-filtering NAs means kth_smallest_c's comparisons are
                # always valid (no NaN/NaT values), so is_datetimelike needs
                # no special handling here.
                tmp = <numeric_t*>malloc(non_na_sz * sizeof(numeric_t))
                if tmp is NULL:
                    raise MemoryError()

                j = 0
                for k in range(grp_size):
                    if mask[start + k] == 0:
                        tmp[j] = values[start + k]
                        j += 1

                for k in range(nqs):
                    q_val = qs[k]

                    # Calculate where to retrieve the desired value.
                    # Casting to int will intentionally truncate result.
                    idx = <int64_t>(q_val * <float64_t>(non_na_sz - 1))

                    # kth_smallest_c is an in-place O(n) quickselect: it
                    # rearranges tmp so that tmp[idx] holds the idx-th
                    # smallest element. Calling it again for a different
                    # index on the same (now partially-sorted) buffer is
                    # always correct because quickselect is correct on any
                    # permutation of the values.
                    val = kth_smallest_c(tmp, idx, non_na_sz)

                    # If requested quantile falls evenly on a particular
                    # index then write that index's value out. Otherwise
                    # interpolate.
                    q_idx = q_val * (non_na_sz - 1)
                    frac = q_idx % 1

                    if frac == 0.0 or interp == INTERPOLATION_LOWER:
                        out[i, k] = val
                    else:
                        # After the previous partition,
                        # tmp[idx+1..non_na_sz-1] are all >= tmp[idx], so
                        # kth_smallest_c correctly finds their minimum (the
                        # (idx+1)-th order statistic).
                        next_val = kth_smallest_c(tmp, idx + 1, non_na_sz)
                        if interp == INTERPOLATION_LINEAR:
                            out[i, k] = val + (next_val - val) * frac
                        elif interp == INTERPOLATION_HIGHER:
                            out[i, k] = next_val
                        elif interp == INTERPOLATION_MIDPOINT:
                            out[i, k] = (val + next_val) / 2.0
                        elif interp == INTERPOLATION_NEAREST:
                            if frac > .5 or (frac == .5 and idx % 2 == 1):
                                # If quantile lies in the middle of two
                                # indexes, take the even index, as np.quantile
                                out[i, k] = next_val
                            else:
                                out[i, k] = val

                free(tmp)


# ----------------------------------------------------------------------
# group_nth, group_last, group_rank
# ----------------------------------------------------------------------

ctypedef fused numeric_object_complex_t:
    numeric_object_t
    complex64_t
    complex128_t


cdef bint _treat_as_na(numeric_object_complex_t val,
                       bint is_datetimelike) noexcept nogil:
    if numeric_object_complex_t is object:
        with gil:
            return checknull(val)

    elif numeric_object_complex_t is int64_t:
        return is_datetimelike and val == NPY_NAT
    elif (
        numeric_object_complex_t is float32_t
        or numeric_object_complex_t is float64_t
        or numeric_object_complex_t is complex64_t
        or numeric_object_complex_t is complex128_t
    ):
        return val != val
    else:
        # non-datetimelike integer
        return False


cdef numeric_object_t _get_min_or_max(
    numeric_object_t val,
    bint compute_max,
    bint is_datetimelike,
):
    """
    Find either the min or the max supported by numeric_object_t; 'val' is a
    placeholder to effectively make numeric_object_t an argument.
    """
    return get_rank_nan_fill_val(
        not compute_max,
        val=val,
        is_datetimelike=is_datetimelike,
    )


cdef numeric_t _get_na_val(numeric_t val, bint is_datetimelike):
    cdef:
        numeric_t na_val

    if numeric_t == float32_t or numeric_t == float64_t:
        na_val = NaN
    elif numeric_t is int64_t and is_datetimelike:
        na_val = NPY_NAT
    else:
        # Used in case of masks
        na_val = 0
    return na_val


ctypedef fused mincount_t:
    numeric_object_t
    complex64_t
    complex128_t


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline void _check_below_mincount(
    mincount_t[:, ::1] out,
    bint uses_mask,
    uint8_t[:, ::1] result_mask,
    Py_ssize_t ncounts,
    Py_ssize_t K,
    const int64_t[:, ::1] nobs,
    int64_t min_count,
    mincount_t[:, ::1] resx,
) noexcept:
    """
    Check if the number of observations for a group is below min_count,
    and if so set the result for that group to the appropriate NA-like value.
    """
    cdef:
        Py_ssize_t i, j

    with nogil(mincount_t is not object):
        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] >= min_count:
                    out[i, j] = resx[i, j]
                else:
                    # if we are integer dtype, not is_datetimelike, and
                    #  not uses_mask, then getting here implies that
                    #  counts[i] < min_count, which means we will
                    #  be cast to float64 and masked at the end
                    #  of WrappedCythonOp._call_cython_op. So we can safely
                    #  set a placeholder value in out[i, j].
                    if uses_mask:
                        result_mask[i, j] = True
                        # set out[i, j] to 0 to be deterministic, as
                        #  it was initialized with np.empty. Also ensures
                        #  we can downcast out if appropriate.
                        out[i, j] = 0
                    elif (
                        mincount_t is float32_t
                        or mincount_t is float64_t
                        or mincount_t is complex64_t
                        or mincount_t is complex128_t
                    ):
                        out[i, j] = NAN
                    elif mincount_t is int64_t:
                        # Per above, this is a placeholder in
                        #  non-is_datetimelike cases.
                        out[i, j] = NPY_NAT
                    elif mincount_t is object:
                        out[i, j] = None
                    else:
                        # placeholder, see above
                        out[i, j] = 0


@cython.wraparound(False)
@cython.boundscheck(False)
def group_last(
    numeric_object_t[:, ::1] out,
    int64_t[::1] counts,
    const numeric_object_t[:, :] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask,
    uint8_t[:, ::1] result_mask=None,
    Py_ssize_t min_count=-1,
    bint is_datetimelike=False,
    bint skipna=True,
) -> None:
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        numeric_object_t val
        numeric_object_t[:, ::1] resx
        int64_t[:, ::1] nobs
        bint uses_mask = mask is not None

    if values.shape[0] != len(labels):
        raise AssertionError("len(index) != len(labels)")

    min_count = max(min_count, 1)

    N = values.shape[0]
    K = values.shape[1]

    nobs_arr = np.zeros((out.shape[0], K), dtype=np.int64)
    nobs = nobs_arr

    if numeric_object_t is object:
        resx_arr = np.empty((out.shape[0], K), dtype=object)
    else:
        resx_arr = np.empty((out.shape[0], K), dtype=np.asarray(out).dtype)
    resx = resx_arr

    with nogil(numeric_object_t is not object):
        if uses_mask:
            if skipna:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    for j in range(K):
                        if not mask[i, j]:
                            nobs[lab, j] += 1
                            resx[lab, j] = values[i, j]
            else:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    for j in range(K):
                        nobs[lab, j] += 1
                        resx[lab, j] = values[i, j]
                        result_mask[lab, j] = mask[i, j]
        else:
            if skipna:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    for j in range(K):
                        val = values[i, j]
                        if not _treat_as_na(val, is_datetimelike):
                            nobs[lab, j] += 1
                            resx[lab, j] = val
            else:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    for j in range(K):
                        nobs[lab, j] += 1
                        resx[lab, j] = values[i, j]

    _check_below_mincount(
        out, uses_mask, result_mask, ncounts, K, nobs, min_count, resx
    )


@cython.wraparound(False)
@cython.boundscheck(False)
def group_nth(
    numeric_object_t[:, ::1] out,
    int64_t[::1] counts,
    const numeric_object_t[:, :] values,
    const intp_t[::1] labels,
    const uint8_t[:, :] mask,
    uint8_t[:, ::1] result_mask=None,
    int64_t min_count=-1,
    int64_t rank=1,
    bint is_datetimelike=False,
    bint skipna=True,
) -> None:
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        numeric_object_t val
        numeric_object_t[:, ::1] resx
        int64_t[:, ::1] nobs
        bint uses_mask = mask is not None
        bint isna_entry

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    min_count = max(min_count, 1)
    nobs = np.zeros((<object>out).shape, dtype=np.int64)
    if numeric_object_t is object:
        resx = np.empty((<object>out).shape, dtype=object)
    else:
        resx = np.empty_like(out)

    N, K = (<object>values).shape

    with nogil(numeric_object_t is not object):
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                if skipna:
                    if uses_mask:
                        isna_entry = mask[i, j]
                    else:
                        isna_entry = _treat_as_na(val, is_datetimelike)
                    if isna_entry:
                        continue

                nobs[lab, j] += 1
                if nobs[lab, j] == rank:
                    resx[lab, j] = val
                    if uses_mask and not skipna:
                        result_mask[lab, j] = mask[i, j]

    _check_below_mincount(
        out, uses_mask, result_mask, ncounts, K, nobs, min_count, resx
    )


@cython.boundscheck(False)
@cython.wraparound(False)
def group_rank(
    float64_t[:, ::1] out,
    ndarray[numeric_object_t, ndim=2] values,
    const intp_t[::1] labels,
    int ngroups,
    bint is_datetimelike,
    str ties_method="average",
    bint ascending=True,
    bint pct=False,
    str na_option="keep",
    const uint8_t[:, :] mask=None,
) -> None:
    """
    Provides the rank of values within each group.

    Parameters
    ----------
    out : np.ndarray[np.float64, ndim=2]
        Values to which this method will write its results.
    values : np.ndarray of numeric_object_t values to be ranked
    labels : np.ndarray[np.intp]
        Array containing unique label for each group, with its ordering
        matching up to the corresponding record in `values`
    ngroups : int
        This parameter is not used, is needed to match signatures of other
        groupby functions.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    ties_method : {'average', 'min', 'max', 'first', 'dense'}, default 'average'
        * average: average rank of group
        * min: lowest rank in group
        * max: highest rank in group
        * first: ranks assigned in order they appear in the array
        * dense: like 'min', but rank always increases by 1 between groups
    ascending : bool, default True
        False for ranks by high (1) to low (N)
        na_option : {'keep', 'top', 'bottom'}, default 'keep'
    pct : bool, default False
        Compute percentage rank of data within each group
    na_option : {'keep', 'top', 'bottom'}, default 'keep'
        * keep: leave NA values where they are
        * top: smallest rank if ascending
        * bottom: smallest rank if descending
    mask : np.ndarray[bool] or None, default None

    Notes
    -----
    This method modifies the `out` parameter rather than returning an object
    """
    cdef:
        Py_ssize_t i, k, N
        ndarray[float64_t, ndim=1] result
        const uint8_t[:] sub_mask

    N = values.shape[1]

    for k in range(N):
        if mask is None:
            sub_mask = None
        else:
            sub_mask = mask[:, k]

        result = rank_1d(
            values=values[:, k],
            labels=labels,
            is_datetimelike=is_datetimelike,
            ties_method=ties_method,
            ascending=ascending,
            pct=pct,
            na_option=na_option,
            mask=sub_mask,
        )
        for i in range(len(result)):
            if labels[i] >= 0:
                out[i, k] = result[i]


# ----------------------------------------------------------------------
# group_min, group_max
# ----------------------------------------------------------------------


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef group_min_max(
    numeric_t[:, ::1] out,
    int64_t[::1] counts,
    const numeric_t[:, :] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    bint is_datetimelike=False,
    bint compute_max=True,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
):
    """
    Compute minimum/maximum  of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : np.ndarray[numeric_t, ndim=2]
        Array to store result in.
    counts : np.ndarray[int64]
        Input as a zeroed array, populated by group sizes during algorithm
    values : array
        Values to find column-wise min/max of.
    labels : np.ndarray[np.intp]
        Labels to group by.
    min_count : Py_ssize_t, default -1
        The minimum number of non-NA group elements, NA result if threshold
        is not met
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    compute_max : bint, default True
        True to compute group-wise max, False to compute min
    mask : ndarray[bool, ndim=2], optional
        If not None, indices represent missing values,
        otherwise the mask will not be used
    result_mask : ndarray[bool, ndim=2], optional
        If not None, these specify locations in the output that are NA.
        Modified in-place.
    skipna : bool, default True
        If True, ignore nans in `values`.

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    `counts` is modified to hold group sizes
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ngroups = len(counts)
        numeric_t val, nan_val
        numeric_t[:, ::1] group_min_or_max
        int64_t[:, ::1] nobs
        bint uses_mask = mask is not None

        const uint8_t *mask_row
        numeric_t *gmm_row
        int64_t *nobs_row

    if values.shape[0] != len(labels):
        raise AssertionError("len(index) != len(labels)")

    min_count = max(min_count, 1)

    N = values.shape[0]
    K = values.shape[1]

    nobs_arr = np.zeros((ngroups, K), dtype=np.int64)
    nobs = nobs_arr
    nan_val = _get_na_val(<numeric_t>0, is_datetimelike)

    gmm_arr = np.empty((ngroups, K), dtype=np.asarray(out).dtype)
    group_min_or_max = gmm_arr
    group_min_or_max[:] = _get_min_or_max(<numeric_t>0, compute_max, is_datetimelike)

    with nogil:
        if compute_max:
            if uses_mask:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1

                    mask_row = &mask[i, 0]
                    nobs_row = &nobs[lab, 0]
                    gmm_row = &group_min_or_max[lab, 0]

                    for j in range(K):
                        if not mask_row[j]:
                            nobs_row[j] += 1
                            if values[i, j] > gmm_row[j]:
                                gmm_row[j] = values[i, j]
                        elif not skipna:
                            result_mask[lab, j] = True
            else:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    nobs_row = &nobs[lab, 0]
                    gmm_row = &group_min_or_max[lab, 0]

                    for j in range(K):
                        val = values[i, j]
                        if not _treat_as_na(val, is_datetimelike):
                            nobs_row[j] += 1
                            if val > gmm_row[j]:
                                gmm_row[j] = val
                        elif not skipna:
                            if not _treat_as_na(gmm_row[j], is_datetimelike):
                                gmm_row[j] = nan_val
        else:
            if uses_mask:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    mask_row = &mask[i, 0]
                    nobs_row = &nobs[lab, 0]
                    gmm_row = &group_min_or_max[lab, 0]

                    for j in range(K):
                        if not mask_row[j]:
                            nobs_row[j] += 1
                            if values[i, j] < gmm_row[j]:
                                gmm_row[j] = values[i, j]
                        elif not skipna:
                            result_mask[lab, j] = True
            else:
                for i in range(N):

                    lab = labels[i]
                    if lab < 0:
                        continue

                    counts[lab] += 1
                    nobs_row = &nobs[lab, 0]
                    gmm_row = &group_min_or_max[lab, 0]

                    for j in range(K):
                        val = values[i, j]
                        if not _treat_as_na(val, is_datetimelike):
                            nobs_row[j] += 1
                            if val < gmm_row[j]:
                                gmm_row[j] = val
                        elif not skipna:
                            if not _treat_as_na(gmm_row[j], is_datetimelike):
                                gmm_row[j] = nan_val

    _check_below_mincount(
        out, uses_mask, result_mask, ngroups, K, nobs, min_count, group_min_or_max
    )


@cython.wraparound(False)
@cython.boundscheck(False)
def group_idxmin_idxmax(
    intp_t[:, ::1] out,
    int64_t[::1] counts,
    ndarray[numeric_object_t, ndim=2] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    bint is_datetimelike=False,
    const uint8_t[:, :] mask=None,
    str name="idxmin",
    bint skipna=True,
    uint8_t[:, ::1] result_mask=None,
):
    """
    Compute index of minimum/maximum of columns of `values`, in row groups `labels`.

    This function only computes the row number where the minimum/maximum occurs, we'll
    take the corresponding index value after this function.

    Parameters
    ----------
    out : np.ndarray[intp, ndim=2]
        Array to store result in.
    counts : np.ndarray[int64]
        Input as a zeroed array, populated by group sizes during algorithm
    values : np.ndarray[numeric_object_t, ndim=2]
        Values to find column-wise min/max of.
    labels : np.ndarray[np.intp]
        Labels to group by.
    min_count : Py_ssize_t, default -1
        The minimum number of non-NA group elements, NA result if threshold
        is not met.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    name : {"idxmin", "idxmax"}, default "idxmin"
        Whether to compute idxmin or idxmax.
    mask : ndarray[bool, ndim=2], optional
        If not None, indices represent missing values,
        otherwise the mask will not be used
    skipna : bool, default True
        Flag to ignore nan values during truth testing
    result_mask : ndarray[bool, ndim=2], optional
        If not None, these specify locations in the output that are NA.
        Modified in-place.

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    `counts` is modified to hold group sizes
    """
    cdef:
        Py_ssize_t i, j, N, K, lab
        numeric_object_t val
        numeric_object_t[:, ::1] group_min_or_max
        uint8_t[:, ::1] seen
        bint uses_mask = mask is not None
        bint isna_entry
        bint compute_max = name == "idxmax"

    assert name == "idxmin" or name == "idxmax"

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    N, K = (<object>values).shape

    if numeric_object_t is object:
        group_min_or_max = np.empty((<object>out).shape, dtype=object)
        seen = np.zeros((<object>out).shape, dtype=np.uint8)
    else:
        group_min_or_max = np.empty_like(out, dtype=values.dtype)
        seen = np.zeros_like(out, dtype=np.uint8)

    # Sentinel for no valid values.
    out[:] = -1

    with nogil(numeric_object_t is not object):
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            for j in range(K):
                if not skipna and out[lab, j] == -1:
                    # Once we've hit NA there is no going back
                    continue

                val = values[i, j]

                if uses_mask:
                    isna_entry = mask[i, j]
                else:
                    isna_entry = _treat_as_na(val, is_datetimelike)

                if isna_entry:
                    if not skipna or not seen[lab, j]:
                        out[lab, j] = -1
                else:
                    if not seen[lab, j]:
                        seen[lab, j] = True
                        group_min_or_max[lab, j] = val
                        out[lab, j] = i
                    elif compute_max:
                        if val > group_min_or_max[lab, j]:
                            group_min_or_max[lab, j] = val
                            out[lab, j] = i
                    else:
                        if val < group_min_or_max[lab, j]:
                            group_min_or_max[lab, j] = val
                            out[lab, j] = i


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
def group_max(
    numeric_t[:, ::1] out,
    int64_t[::1] counts,
    const numeric_t[:, :] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    bint is_datetimelike=False,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    """See group_min_max.__doc__"""
    group_min_max(
        out,
        counts,
        values,
        labels,
        min_count=min_count,
        is_datetimelike=is_datetimelike,
        compute_max=True,
        mask=mask,
        result_mask=result_mask,
        skipna=skipna,
    )


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
def group_min(
    numeric_t[:, ::1] out,
    int64_t[::1] counts,
    const numeric_t[:, :] values,
    const intp_t[::1] labels,
    Py_ssize_t min_count=-1,
    bint is_datetimelike=False,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    """See group_min_max.__doc__"""
    group_min_max(
        out,
        counts,
        values,
        labels,
        min_count=min_count,
        is_datetimelike=is_datetimelike,
        compute_max=False,
        mask=mask,
        result_mask=result_mask,
        skipna=skipna,
    )


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef group_cummin_max(
    numeric_t[:, ::1] out,
    const numeric_t[:, :] values,
    const uint8_t[:, :] mask,
    uint8_t[:, ::1] result_mask,
    const intp_t[::1] labels,
    int ngroups,
    bint is_datetimelike,
    bint skipna,
    bint compute_max,
):
    """
    Cumulative minimum/maximum of columns of `values`, in row groups `labels`.

    Parameters
    ----------
    out : np.ndarray[numeric_t, ndim=2]
        Array to store cummin/max in.
    values : np.ndarray[numeric_t, ndim=2]
        Values to take cummin/max of.
    mask : np.ndarray[bool] or None
        If not None, indices represent missing values,
        otherwise the mask will not be used
    result_mask : ndarray[bool, ndim=2], optional
        If not None, these specify locations in the output that are NA.
        Modified in-place.
    labels : np.ndarray[np.intp]
        Labels to group by.
    ngroups : int
        Number of groups, larger than all entries of `labels`.
    is_datetimelike : bool
        True if `values` contains datetime-like entries.
    skipna : bool
        If True, ignore nans in `values`.
    compute_max : bool
        True if cumulative maximum should be computed, False
        if cumulative minimum should be computed

    Notes
    -----
    This method modifies the `out` parameter, rather than returning an object.
    """
    cdef:
        numeric_t[:, ::1] accum
        Py_ssize_t i, j, N, K
        numeric_t val, mval, na_val
        uint8_t[:, ::1] seen_na
        intp_t lab
        bint na_possible
        bint uses_mask = mask is not None
        bint isna_entry
        bint check_na

        const uint8_t *mask_row = NULL
        numeric_t *out_row = NULL
        numeric_t *accum_row = NULL
        uint8_t *res_mask_row = NULL
        uint8_t *seen_na_row = NULL

    N = values.shape[0]
    K = values.shape[1]

    accum = np.empty((ngroups, K), dtype=np.asarray(values).dtype)
    accum[:] = _get_min_or_max(<numeric_t>0, compute_max, is_datetimelike)

    na_val = _get_na_val(<numeric_t>0, is_datetimelike)

    if uses_mask:
        na_possible = True
        na_val = 0
    elif numeric_t is float64_t or numeric_t is float32_t:
        na_possible = True
    elif is_datetimelike:
        na_possible = True
    else:
        na_possible = False

    check_na = not skipna and na_possible

    if na_possible:
        seen_na = np.zeros((ngroups, K), dtype=np.uint8)

    with nogil:
        if compute_max:
            if uses_mask:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    mask_row = &mask[i, 0]
                    out_row = &out[i, 0]
                    res_mask_row = &result_mask[i, 0]
                    accum_row = &accum[lab, 0]
                    if check_na:
                        seen_na_row = &seen_na[lab, 0]

                    for j in range(K):
                        if check_na and seen_na_row[j]:
                            res_mask_row[j] = 1
                            out_row[j] = 0
                        else:
                            val = values[i, j]
                            isna_entry = mask_row[j]
                            if not isna_entry:
                                mval = accum_row[j]
                                if val > mval:
                                    accum_row[j] = val
                                    mval = val
                                out_row[j] = mval
                            else:
                                if check_na:
                                    seen_na_row[j] = 1
                                out_row[j] = val
            else:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    out_row = &out[i, 0]
                    accum_row = &accum[lab, 0]
                    if check_na:
                        seen_na_row = &seen_na[lab, 0]

                    for j in range(K):
                        if check_na and seen_na_row[j]:
                            out_row[j] = na_val
                        else:
                            val = values[i, j]
                            isna_entry = _treat_as_na(val, is_datetimelike)
                            if not isna_entry:
                                mval = accum_row[j]
                                if val > mval:
                                    accum_row[j] = val
                                    mval = val
                                out_row[j] = mval
                            else:
                                if check_na:
                                    seen_na_row[j] = 1
                                out_row[j] = val
        else:
            if uses_mask:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    mask_row = &mask[i, 0]
                    out_row = &out[i, 0]
                    res_mask_row = &result_mask[i, 0]
                    accum_row = &accum[lab, 0]
                    if check_na:
                        seen_na_row = &seen_na[lab, 0]

                    for j in range(K):
                        if check_na and seen_na_row[j]:
                            res_mask_row[j] = 1
                            out_row[j] = 0
                        else:
                            val = values[i, j]
                            isna_entry = mask_row[j]
                            if not isna_entry:
                                mval = accum_row[j]
                                if val < mval:
                                    accum_row[j] = val
                                    mval = val
                                out_row[j] = mval
                            else:
                                if check_na:
                                    seen_na_row[j] = 1
                                out_row[j] = val
            else:
                for i in range(N):
                    lab = labels[i]
                    if lab < 0:
                        continue

                    out_row = &out[i, 0]
                    accum_row = &accum[lab, 0]
                    if check_na:
                        seen_na_row = &seen_na[lab, 0]

                    for j in range(K):
                        if check_na and seen_na_row[j]:
                            out_row[j] = na_val
                        else:
                            val = values[i, j]
                            isna_entry = _treat_as_na(val, is_datetimelike)
                            if not isna_entry:
                                mval = accum_row[j]
                                if val < mval:
                                    accum_row[j] = val
                                    mval = val
                                out_row[j] = mval
                            else:
                                if check_na:
                                    seen_na_row[j] = 1
                                out_row[j] = val


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cummin(
    numeric_t[:, ::1] out,
    const numeric_t[:, :] values,
    const intp_t[::1] labels,
    int ngroups,
    bint is_datetimelike,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    """See group_cummin_max.__doc__"""

    group_cummin_max(
        out,
        values,
        mask,
        result_mask,
        labels,
        ngroups,
        is_datetimelike,
        skipna,
        False,  # compute_max = False
    )


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cummax(
    numeric_t[:, ::1] out,
    const numeric_t[:, :] values,
    const intp_t[::1] labels,
    int ngroups,
    bint is_datetimelike,
    const uint8_t[:, :] mask=None,
    uint8_t[:, ::1] result_mask=None,
    bint skipna=True,
) -> None:
    """See group_cummin_max.__doc__"""
    group_cummin_max(
        out,
        values,
        mask,
        result_mask,
        labels,
        ngroups,
        is_datetimelike,
        skipna,
        True,  # compute_max = True
    )
