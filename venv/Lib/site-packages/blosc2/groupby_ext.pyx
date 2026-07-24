#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################
# cython: boundscheck=False, wraparound=False, initializedcheck=False

"""Cython group-reduce kernels for CTable group_by()."""

import numpy as np
cimport numpy as np

from libc.stdint cimport int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy


# ----------------------------------------------------------------------
# Group-reduce kernels
# ----------------------------------------------------------------------

def groupby_dense_i32_f64_sum(
    np.ndarray keys,
    np.ndarray values,
    np.ndarray valid,
    np.ndarray sums,
    np.ndarray present,
    bint skip_key_null=False,
    int32_t key_null=0,
    bint skip_value_nan=False,
):
    """Accumulate ``sum(values)`` by dense int32 keys.

    This is a low-level CTable group-by helper.  *keys*, *values*, and *valid*
    are same-length 1-D chunk arrays.  *sums* and *present* are dense group
    state arrays indexed directly by key value.  Keys must be non-negative and
    already fit in the state arrays.
    """
    if keys.ndim != 1 or values.ndim != 1 or valid.ndim != 1:
        raise ValueError("keys, values and valid must be 1-D arrays")
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.ndim != 1 or present.ndim != 1:
        raise ValueError("sums and present must be 1-D arrays")
    if keys.dtype != np.dtype(np.int32):
        raise TypeError("keys must have dtype int32")
    if values.dtype != np.dtype(np.float64):
        raise TypeError("values must have dtype float64")
    if valid.dtype != np.dtype(np.bool_):
        raise TypeError("valid must have dtype bool")
    if sums.dtype != np.dtype(np.float64):
        raise TypeError("sums must have dtype float64")
    if present.dtype != np.dtype(np.bool_):
        raise TypeError("present must have dtype bool")

    cdef int32_t[:] keys_view = keys
    cdef double[:] values_view = values
    cdef np.npy_bool[:] valid_view = valid
    cdef double[:] sums_view = sums
    cdef np.npy_bool[:] present_view = present
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef int32_t key
    cdef double value

    if present.shape[0] != sums.shape[0]:
        raise ValueError("present and sums must have the same length")

    with nogil:
        for i in range(n):
            if not valid_view[i]:
                continue
            key = keys_view[i]
            if skip_key_null and key == key_null:
                continue
            if key < 0 or key >= nstates:
                continue
            value = values_view[i]
            if skip_value_nan and value != value:
                continue
            sums_view[key] += value
            present_view[key] = 1
    return None


def groupby_dense_i32_f64_sum_checked(
    np.ndarray keys,
    np.ndarray values,
    np.ndarray valid,
    np.ndarray sums,
    np.ndarray present,
    bint skip_key_null=False,
    int32_t key_null=0,
    bint skip_value_nan=False,
):
    """Checked dense int32/float64 sum kernel.

    Returns ``0`` on success, ``-1`` if a negative non-null key is found, or
    ``max_key + 1`` when the dense state arrays need to be grown.  The state is
    not mutated unless the function returns ``0``.
    """
    if keys.ndim != 1 or values.ndim != 1 or valid.ndim != 1:
        raise ValueError("keys, values and valid must be 1-D arrays")
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.ndim != 1 or present.ndim != 1:
        raise ValueError("sums and present must be 1-D arrays")
    if keys.dtype != np.dtype(np.int32):
        raise TypeError("keys must have dtype int32")
    if values.dtype != np.dtype(np.float64):
        raise TypeError("values must have dtype float64")
    if valid.dtype != np.dtype(np.bool_):
        raise TypeError("valid must have dtype bool")
    if sums.dtype != np.dtype(np.float64):
        raise TypeError("sums must have dtype float64")
    if present.dtype != np.dtype(np.bool_):
        raise TypeError("present must have dtype bool")
    if present.shape[0] != sums.shape[0]:
        raise ValueError("present and sums must have the same length")

    cdef int32_t[:] keys_view = keys
    cdef double[:] values_view = values
    cdef np.npy_bool[:] valid_view = valid
    cdef double[:] sums_view = sums
    cdef np.npy_bool[:] present_view = present
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef int32_t key
    cdef int32_t max_key = -1
    cdef int ret = 0
    cdef double value

    with nogil:
        for i in range(n):
            if not valid_view[i]:
                continue
            key = keys_view[i]
            if skip_key_null and key == key_null:
                continue
            if key < 0:
                ret = -1
                break
            if key > max_key:
                max_key = key
        if ret == 0:
            if max_key < 0:
                ret = 0
            elif max_key >= nstates:
                ret = <int>max_key + 1
            else:
                for i in range(n):
                    if not valid_view[i]:
                        continue
                    key = keys_view[i]
                    if skip_key_null and key == key_null:
                        continue
                    value = values_view[i]
                    if skip_value_nan and value != value:
                        continue
                    sums_view[key] += value
                    present_view[key] = 1
    return ret


def groupby_dense_f64_integral_key_f64_sum_checked(
    np.ndarray keys,
    np.ndarray values,
    np.ndarray valid,
    np.ndarray sums,
    np.ndarray present,
    bint skip_key_nan=True,
    bint skip_value_nan=False,
):
    """Checked dense float64-integral-key/float64 sum kernel.

    Fast path for float keys that are exactly integral, finite and
    non-negative.  Returns ``0`` on success, ``-1`` if a key cannot be handled,
    or ``max_key + 1`` when the dense state arrays need to be grown.  The state is
    not mutated unless the function returns ``0``.
    """
    if keys.ndim != 1 or values.ndim != 1 or valid.ndim != 1:
        raise ValueError("keys, values and valid must be 1-D arrays")
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.ndim != 1 or present.ndim != 1:
        raise ValueError("sums and present must be 1-D arrays")
    if keys.dtype != np.dtype(np.float64):
        raise TypeError("keys must have dtype float64")
    if values.dtype != np.dtype(np.float64):
        raise TypeError("values must have dtype float64")
    if valid.dtype != np.dtype(np.bool_):
        raise TypeError("valid must have dtype bool")
    if sums.dtype != np.dtype(np.float64):
        raise TypeError("sums must have dtype float64")
    if present.dtype != np.dtype(np.bool_):
        raise TypeError("present must have dtype bool")
    if present.shape[0] != sums.shape[0]:
        raise ValueError("present and sums must have the same length")

    cdef double[:] keys_view = keys
    cdef double[:] values_view = values
    cdef np.npy_bool[:] valid_view = valid
    cdef double[:] sums_view = sums
    cdef np.npy_bool[:] present_view = present
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef double key_f
    cdef int64_t key_i
    cdef int64_t max_key = -1
    cdef int ret = 0
    cdef double value

    with nogil:
        for i in range(n):
            if not valid_view[i]:
                continue
            key_f = keys_view[i]
            if key_f != key_f:
                if skip_key_nan:
                    continue
                ret = -1
                break
            if key_f < 0.0 or key_f > 9223372036854774784.0:
                ret = -1
                break
            key_i = <int64_t>key_f
            if key_f != <double>key_i:
                ret = -1
                break
            if key_i > max_key:
                max_key = key_i
        if ret == 0:
            if max_key < 0:
                ret = 0
            elif max_key >= nstates:
                if max_key > 2147483646:
                    ret = -1
                else:
                    ret = <int>max_key + 1
            else:
                for i in range(n):
                    if not valid_view[i]:
                        continue
                    key_f = keys_view[i]
                    if key_f != key_f:
                        if skip_key_nan:
                            continue
                        ret = -1
                        break
                    key_i = <int64_t>key_f
                    value = values_view[i]
                    if skip_value_nan and value != value:
                        continue
                    sums_view[key_i] += value
                    present_view[key_i] = 1
    return ret


def groupby_dense_f32_integral_key_f64_sum_checked(
    np.ndarray keys,
    np.ndarray values,
    np.ndarray valid,
    np.ndarray sums,
    np.ndarray present,
    bint skip_key_nan=True,
    bint skip_value_nan=False,
):
    """Checked dense float32-integral-key/float64 sum kernel."""
    if keys.ndim != 1 or values.ndim != 1 or valid.ndim != 1:
        raise ValueError("keys, values and valid must be 1-D arrays")
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.ndim != 1 or present.ndim != 1:
        raise ValueError("sums and present must be 1-D arrays")
    if keys.dtype != np.dtype(np.float32):
        raise TypeError("keys must have dtype float32")
    if values.dtype != np.dtype(np.float64):
        raise TypeError("values must have dtype float64")
    if valid.dtype != np.dtype(np.bool_):
        raise TypeError("valid must have dtype bool")
    if sums.dtype != np.dtype(np.float64):
        raise TypeError("sums must have dtype float64")
    if present.dtype != np.dtype(np.bool_):
        raise TypeError("present must have dtype bool")
    if present.shape[0] != sums.shape[0]:
        raise ValueError("present and sums must have the same length")

    cdef float[:] keys_view = keys
    cdef double[:] values_view = values
    cdef np.npy_bool[:] valid_view = valid
    cdef double[:] sums_view = sums
    cdef np.npy_bool[:] present_view = present
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef float key_f
    cdef int64_t key_i
    cdef int64_t max_key = -1
    cdef int ret = 0
    cdef double value

    with nogil:
        for i in range(n):
            if not valid_view[i]:
                continue
            key_f = keys_view[i]
            if key_f != key_f:
                if skip_key_nan:
                    continue
                ret = -1
                break
            if key_f < 0.0 or key_f > 16777216.0:
                ret = -1
                break
            key_i = <int64_t>key_f
            if key_f != <float>key_i:
                ret = -1
                break
            if key_i > max_key:
                max_key = key_i
        if ret == 0:
            if max_key < 0:
                ret = 0
            elif max_key >= nstates:
                if max_key > 2147483646:
                    ret = -1
                else:
                    ret = <int>max_key + 1
            else:
                for i in range(n):
                    if not valid_view[i]:
                        continue
                    key_f = keys_view[i]
                    if key_f != key_f:
                        if skip_key_nan:
                            continue
                        ret = -1
                        break
                    key_i = <int64_t>key_f
                    value = values_view[i]
                    if skip_value_nan and value != value:
                        continue
                    sums_view[key_i] += value
                    present_view[key_i] = 1
    return ret


# ----------------------------------------------------------------------
# Fused integer-key dense kernels
# ----------------------------------------------------------------------

ctypedef fused dense_int_key_t:
    int8_t
    uint8_t
    int16_t
    uint16_t
    int32_t
    uint32_t
    int64_t
    uint64_t


cdef inline int _dense_int_key_scan(
    dense_int_key_t[:] keys_view,
    np.npy_bool[:] valid_view,
    Py_ssize_t n,
    Py_ssize_t nstates,
    bint skip_key_null,
    int64_t key_null,
    int* ret,
) noexcept nogil:
    cdef Py_ssize_t i
    cdef int64_t key
    cdef int64_t max_key = -1
    ret[0] = 0
    for i in range(n):
        if not valid_view[i]:
            continue
        key = <int64_t>keys_view[i]
        if skip_key_null and key == key_null:
            continue
        if key < 0:
            ret[0] = -1
            return 0
        if key > max_key:
            max_key = key
    if max_key < 0:
        ret[0] = 0
    elif max_key >= nstates:
        if max_key > 2147483646:
            ret[0] = -1
        else:
            ret[0] = <int>max_key + 1
    return 0


def groupby_dense_int_size_checked(
    dense_int_key_t[:] keys,
    np.npy_bool[:] valid,
    int64_t[:] counts,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
):
    """Checked dense integer-key ``size`` kernel for all integer key widths."""
    if keys.shape[0] != valid.shape[0]:
        raise ValueError("keys and valid must have the same length")
    if counts.shape[0] != keys_present.shape[0]:
        raise ValueError("counts and keys_present must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = counts.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                counts[key] += 1
                keys_present[key] = 1
    return ret


def groupby_dense_int_count_checked(
    dense_int_key_t[:] keys,
    np.npy_bool[:] valid,
    np.npy_bool[:] values_valid,
    int64_t[:] counts,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
):
    """Checked dense integer-key non-null count kernel."""
    if keys.shape[0] != valid.shape[0] or keys.shape[0] != values_valid.shape[0]:
        raise ValueError("keys, valid and values_valid must have the same length")
    if counts.shape[0] != keys_present.shape[0]:
        raise ValueError("counts and keys_present must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = counts.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                if values_valid[i]:
                    counts[key] += 1
    return ret


def groupby_dense_int_f64_sum_checked(
    dense_int_key_t[:] keys,
    double[:] values,
    np.npy_bool[:] valid,
    double[:] sums,
    np.npy_bool[:] value_present,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
    bint skip_value_nan=False,
):
    """Checked dense integer-key float64 sum kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.shape[0] != value_present.shape[0] or sums.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef double value
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                value = values[i]
                if skip_value_nan and value != value:
                    continue
                sums[key] += value
                value_present[key] = 1
    return ret


def groupby_dense_int_i64_sum_checked(
    dense_int_key_t[:] keys,
    int64_t[:] values,
    np.npy_bool[:] valid,
    int64_t[:] sums,
    np.npy_bool[:] value_present,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
):
    """Checked dense integer-key int64 sum kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.shape[0] != value_present.shape[0] or sums.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                sums[key] += values[i]
                value_present[key] = 1
    return ret


def groupby_dense_int_f64_mean_checked(
    dense_int_key_t[:] keys,
    double[:] values,
    np.npy_bool[:] valid,
    double[:] sums,
    int64_t[:] counts,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
    bint skip_value_nan=False,
):
    """Checked dense integer-key float64 mean state kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if sums.shape[0] != counts.shape[0] or sums.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = sums.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef double value
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                value = values[i]
                if skip_value_nan and value != value:
                    continue
                sums[key] += value
                counts[key] += 1
    return ret


def groupby_dense_int_f64_min_checked(
    dense_int_key_t[:] keys,
    double[:] values,
    np.npy_bool[:] valid,
    double[:] mins,
    np.npy_bool[:] has_value,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
    bint skip_value_nan=False,
):
    """Checked dense integer-key float64 min kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if mins.shape[0] != has_value.shape[0] or mins.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = mins.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef double value
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                value = values[i]
                if skip_value_nan and value != value:
                    continue
                if not has_value[key] or value < mins[key]:
                    mins[key] = value
                has_value[key] = 1
    return ret


def groupby_dense_int_f64_max_checked(
    dense_int_key_t[:] keys,
    double[:] values,
    np.npy_bool[:] valid,
    double[:] maxs,
    np.npy_bool[:] has_value,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
    bint skip_value_nan=False,
):
    """Checked dense integer-key float64 max kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if maxs.shape[0] != has_value.shape[0] or maxs.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = maxs.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef double value
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                value = values[i]
                if skip_value_nan and value != value:
                    continue
                if not has_value[key] or value > maxs[key]:
                    maxs[key] = value
                has_value[key] = 1
    return ret


def groupby_dense_int_i64_min_checked(
    dense_int_key_t[:] keys,
    int64_t[:] values,
    np.npy_bool[:] valid,
    int64_t[:] mins,
    np.npy_bool[:] has_value,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
):
    """Checked dense integer-key int64 min kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if mins.shape[0] != has_value.shape[0] or mins.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = mins.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef int64_t value
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                value = values[i]
                if not has_value[key] or value < mins[key]:
                    mins[key] = value
                has_value[key] = 1
    return ret


def groupby_dense_int_i64_max_checked(
    dense_int_key_t[:] keys,
    int64_t[:] values,
    np.npy_bool[:] valid,
    int64_t[:] maxs,
    np.npy_bool[:] has_value,
    np.npy_bool[:] keys_present,
    bint skip_key_null=False,
    int64_t key_null=0,
):
    """Checked dense integer-key int64 max kernel."""
    if keys.shape[0] != values.shape[0] or keys.shape[0] != valid.shape[0]:
        raise ValueError("keys, values and valid must have the same length")
    if maxs.shape[0] != has_value.shape[0] or maxs.shape[0] != keys_present.shape[0]:
        raise ValueError("state arrays must have the same length")
    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t nstates = maxs.shape[0]
    cdef Py_ssize_t i
    cdef int64_t key
    cdef int64_t value
    cdef int ret
    with nogil:
        _dense_int_key_scan(keys, valid, n, nstates, skip_key_null, key_null, &ret)
        if ret == 0:
            for i in range(n):
                if not valid[i]:
                    continue
                key = <int64_t>keys[i]
                if skip_key_null and key == key_null:
                    continue
                keys_present[key] = 1
                value = values[i]
                if not has_value[key] or value > maxs[key]:
                    maxs[key] = value
                has_value[key] = 1
    return ret


# ----------------------------------------------------------------------
# Arbitrary float-key hash kernels
# ----------------------------------------------------------------------

cdef inline uint64_t _f64_bits(double value) noexcept:
    cdef uint64_t bits
    memcpy(&bits, &value, sizeof(double))
    return bits


cdef inline uint64_t _mix_u64(uint64_t x) noexcept:
    x ^= x >> 30
    x *= <uint64_t>0xbf58476d1ce4e5b9
    x ^= x >> 27
    x *= <uint64_t>0x94d049bb133111eb
    x ^= x >> 31
    return x


def groupby_hash_f64_f64(
    double[:] keys,
    double[:] values,
    np.npy_bool[:] valid,
    np.npy_bool[:] values_valid,
    bint has_values,
    bint dropna=True,
):
    """Hash arbitrary float64 keys and accumulate float64 group states.

    Returns ``(keys, row_counts, value_counts, sums, mins, maxs, has_value)``.
    NaN keys are skipped when ``dropna`` is true; otherwise all NaN bit-patterns
    are normalized into one NaN group.  ``+0.0`` and ``-0.0`` are normalized into
    the same zero group.
    """
    if keys.shape[0] != valid.shape[0]:
        raise ValueError("keys and valid must have the same length")
    if has_values and (values.shape[0] != keys.shape[0] or values_valid.shape[0] != keys.shape[0]):
        raise ValueError("values, values_valid and keys must have the same length")

    cdef Py_ssize_t n = keys.shape[0]
    cdef Py_ssize_t cap = 1024
    cdef Py_ssize_t used_count = 0
    cdef Py_ssize_t i, pos, old_pos, out_pos
    cdef uint64_t mask = <uint64_t>cap - 1
    cdef uint64_t bits, h, old_bits
    cdef double key, value
    cdef double nan_value = float("nan")
    cdef uint64_t nan_bits = <uint64_t>0x7ff8000000000000
    cdef bint value_ok

    cdef uint64_t* table_bits = <uint64_t*>malloc(cap * sizeof(uint64_t))
    cdef np.npy_bool* table_used = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))
    cdef double* table_keys = <double*>malloc(cap * sizeof(double))
    cdef int64_t* row_counts = <int64_t*>malloc(cap * sizeof(int64_t))
    cdef int64_t* value_counts = <int64_t*>malloc(cap * sizeof(int64_t))
    cdef double* sums = <double*>malloc(cap * sizeof(double))
    cdef double* mins = <double*>malloc(cap * sizeof(double))
    cdef double* maxs = <double*>malloc(cap * sizeof(double))
    cdef np.npy_bool* has_value = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))

    cdef uint64_t* new_bits
    cdef np.npy_bool* new_used
    cdef double* new_keys
    cdef int64_t* new_row_counts
    cdef int64_t* new_value_counts
    cdef double* new_sums
    cdef double* new_mins
    cdef double* new_maxs
    cdef np.npy_bool* new_has_value
    cdef Py_ssize_t old_cap
    cdef uint64_t new_mask

    if (
        table_bits == NULL
        or table_used == NULL
        or table_keys == NULL
        or row_counts == NULL
        or value_counts == NULL
        or sums == NULL
        or mins == NULL
        or maxs == NULL
        or has_value == NULL
    ):
        free(table_bits); free(table_used); free(table_keys); free(row_counts); free(value_counts)
        free(sums); free(mins); free(maxs); free(has_value)
        raise MemoryError()

    for i in range(cap):
        table_used[i] = 0

    try:
        for i in range(n):
            if not valid[i]:
                continue
            key = keys[i]
            if key != key:
                if dropna:
                    continue
                bits = nan_bits
                key = nan_value
            elif key == 0.0:
                key = 0.0
                bits = 0
            else:
                bits = _f64_bits(key)

            if (used_count + 1) * 2 >= cap:
                old_cap = cap
                cap *= 2
                mask = <uint64_t>cap - 1
                new_bits = <uint64_t*>malloc(cap * sizeof(uint64_t))
                new_used = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))
                new_keys = <double*>malloc(cap * sizeof(double))
                new_row_counts = <int64_t*>malloc(cap * sizeof(int64_t))
                new_value_counts = <int64_t*>malloc(cap * sizeof(int64_t))
                new_sums = <double*>malloc(cap * sizeof(double))
                new_mins = <double*>malloc(cap * sizeof(double))
                new_maxs = <double*>malloc(cap * sizeof(double))
                new_has_value = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))
                if (
                    new_bits == NULL
                    or new_used == NULL
                    or new_keys == NULL
                    or new_row_counts == NULL
                    or new_value_counts == NULL
                    or new_sums == NULL
                    or new_mins == NULL
                    or new_maxs == NULL
                    or new_has_value == NULL
                ):
                    free(new_bits); free(new_used); free(new_keys); free(new_row_counts); free(new_value_counts)
                    free(new_sums); free(new_mins); free(new_maxs); free(new_has_value)
                    raise MemoryError()
                for pos in range(cap):
                    new_used[pos] = 0
                for old_pos in range(old_cap):
                    if not table_used[old_pos]:
                        continue
                    old_bits = table_bits[old_pos]
                    h = _mix_u64(old_bits)
                    pos = <Py_ssize_t>(h & mask)
                    while new_used[pos]:
                        pos = <Py_ssize_t>((pos + 1) & mask)
                    new_used[pos] = 1
                    new_bits[pos] = old_bits
                    new_keys[pos] = table_keys[old_pos]
                    new_row_counts[pos] = row_counts[old_pos]
                    new_value_counts[pos] = value_counts[old_pos]
                    new_sums[pos] = sums[old_pos]
                    new_mins[pos] = mins[old_pos]
                    new_maxs[pos] = maxs[old_pos]
                    new_has_value[pos] = has_value[old_pos]
                free(table_bits); free(table_used); free(table_keys); free(row_counts); free(value_counts)
                free(sums); free(mins); free(maxs); free(has_value)
                table_bits = new_bits
                table_used = new_used
                table_keys = new_keys
                row_counts = new_row_counts
                value_counts = new_value_counts
                sums = new_sums
                mins = new_mins
                maxs = new_maxs
                has_value = new_has_value

            h = _mix_u64(bits)
            pos = <Py_ssize_t>(h & mask)
            while table_used[pos] and table_bits[pos] != bits:
                pos = <Py_ssize_t>((pos + 1) & mask)
            if not table_used[pos]:
                table_used[pos] = 1
                table_bits[pos] = bits
                table_keys[pos] = key
                row_counts[pos] = 0
                value_counts[pos] = 0
                sums[pos] = 0.0
                mins[pos] = 0.0
                maxs[pos] = 0.0
                has_value[pos] = 0
                used_count += 1

            row_counts[pos] += 1
            if has_values:
                value_ok = values_valid[i]
                if value_ok:
                    value = values[i]
                    value_counts[pos] += 1
                    sums[pos] += value
                    if not has_value[pos] or value < mins[pos]:
                        mins[pos] = value
                    if not has_value[pos] or value > maxs[pos]:
                        maxs[pos] = value
                    has_value[pos] = 1

        out_keys = np.empty(used_count, dtype=np.float64)
        out_row_counts = np.empty(used_count, dtype=np.int64)
        out_value_counts = np.empty(used_count, dtype=np.int64)
        out_sums = np.empty(used_count, dtype=np.float64)
        out_mins = np.empty(used_count, dtype=np.float64)
        out_maxs = np.empty(used_count, dtype=np.float64)
        out_has_value = np.empty(used_count, dtype=bool)

        out_pos = 0
        for pos in range(cap):
            if not table_used[pos]:
                continue
            out_keys[out_pos] = table_keys[pos]
            out_row_counts[out_pos] = row_counts[pos]
            out_value_counts[out_pos] = value_counts[pos]
            out_sums[out_pos] = sums[pos]
            out_mins[out_pos] = mins[pos]
            out_maxs[out_pos] = maxs[pos]
            out_has_value[out_pos] = has_value[pos]
            out_pos += 1
        return out_keys, out_row_counts, out_value_counts, out_sums, out_mins, out_maxs, out_has_value
    finally:
        free(table_bits); free(table_used); free(table_keys); free(row_counts); free(value_counts)
        free(sums); free(mins); free(maxs); free(has_value)


def groupby_hash_i64x2_f64(
    int64_t[:] key0,
    int64_t[:] key1,
    double[:] values,
    np.npy_bool[:] valid,
    np.npy_bool[:] values_valid,
    bint has_values,
):
    """Hash two int64-normalized keys and accumulate float64 group states."""
    if key0.shape[0] != key1.shape[0] or key0.shape[0] != valid.shape[0]:
        raise ValueError("key0, key1 and valid must have the same length")
    if has_values and (values.shape[0] != key0.shape[0] or values_valid.shape[0] != key0.shape[0]):
        raise ValueError("values, values_valid and keys must have the same length")

    cdef Py_ssize_t n = key0.shape[0]
    cdef Py_ssize_t cap = 1024
    cdef Py_ssize_t used_count = 0
    cdef Py_ssize_t i, pos, old_pos, out_pos
    cdef uint64_t mask = <uint64_t>cap - 1
    cdef uint64_t h
    cdef int64_t k0, k1
    cdef double value
    cdef bint value_ok

    cdef int64_t* table_k0 = <int64_t*>malloc(cap * sizeof(int64_t))
    cdef int64_t* table_k1 = <int64_t*>malloc(cap * sizeof(int64_t))
    cdef np.npy_bool* table_used = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))
    cdef int64_t* row_counts = <int64_t*>malloc(cap * sizeof(int64_t))
    cdef int64_t* value_counts = <int64_t*>malloc(cap * sizeof(int64_t))
    cdef double* sums = <double*>malloc(cap * sizeof(double))
    cdef double* mins = <double*>malloc(cap * sizeof(double))
    cdef double* maxs = <double*>malloc(cap * sizeof(double))
    cdef np.npy_bool* has_value = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))

    cdef int64_t* new_k0
    cdef int64_t* new_k1
    cdef np.npy_bool* new_used
    cdef int64_t* new_row_counts
    cdef int64_t* new_value_counts
    cdef double* new_sums
    cdef double* new_mins
    cdef double* new_maxs
    cdef np.npy_bool* new_has_value
    cdef Py_ssize_t old_cap

    if (
        table_k0 == NULL
        or table_k1 == NULL
        or table_used == NULL
        or row_counts == NULL
        or value_counts == NULL
        or sums == NULL
        or mins == NULL
        or maxs == NULL
        or has_value == NULL
    ):
        free(table_k0); free(table_k1); free(table_used); free(row_counts); free(value_counts)
        free(sums); free(mins); free(maxs); free(has_value)
        raise MemoryError()

    for i in range(cap):
        table_used[i] = 0

    try:
        for i in range(n):
            if not valid[i]:
                continue
            k0 = key0[i]
            k1 = key1[i]

            if (used_count + 1) * 2 >= cap:
                old_cap = cap
                cap *= 2
                mask = <uint64_t>cap - 1
                new_k0 = <int64_t*>malloc(cap * sizeof(int64_t))
                new_k1 = <int64_t*>malloc(cap * sizeof(int64_t))
                new_used = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))
                new_row_counts = <int64_t*>malloc(cap * sizeof(int64_t))
                new_value_counts = <int64_t*>malloc(cap * sizeof(int64_t))
                new_sums = <double*>malloc(cap * sizeof(double))
                new_mins = <double*>malloc(cap * sizeof(double))
                new_maxs = <double*>malloc(cap * sizeof(double))
                new_has_value = <np.npy_bool*>malloc(cap * sizeof(np.npy_bool))
                if (
                    new_k0 == NULL
                    or new_k1 == NULL
                    or new_used == NULL
                    or new_row_counts == NULL
                    or new_value_counts == NULL
                    or new_sums == NULL
                    or new_mins == NULL
                    or new_maxs == NULL
                    or new_has_value == NULL
                ):
                    free(new_k0); free(new_k1); free(new_used); free(new_row_counts); free(new_value_counts)
                    free(new_sums); free(new_mins); free(new_maxs); free(new_has_value)
                    raise MemoryError()
                for pos in range(cap):
                    new_used[pos] = 0
                for old_pos in range(old_cap):
                    if not table_used[old_pos]:
                        continue
                    h = _mix_u64(<uint64_t>table_k0[old_pos]) ^ _mix_u64(<uint64_t>table_k1[old_pos] + <uint64_t>0x9e3779b97f4a7c15)
                    pos = <Py_ssize_t>(h & mask)
                    while new_used[pos]:
                        pos = <Py_ssize_t>((pos + 1) & mask)
                    new_used[pos] = 1
                    new_k0[pos] = table_k0[old_pos]
                    new_k1[pos] = table_k1[old_pos]
                    new_row_counts[pos] = row_counts[old_pos]
                    new_value_counts[pos] = value_counts[old_pos]
                    new_sums[pos] = sums[old_pos]
                    new_mins[pos] = mins[old_pos]
                    new_maxs[pos] = maxs[old_pos]
                    new_has_value[pos] = has_value[old_pos]
                free(table_k0); free(table_k1); free(table_used); free(row_counts); free(value_counts)
                free(sums); free(mins); free(maxs); free(has_value)
                table_k0 = new_k0
                table_k1 = new_k1
                table_used = new_used
                row_counts = new_row_counts
                value_counts = new_value_counts
                sums = new_sums
                mins = new_mins
                maxs = new_maxs
                has_value = new_has_value

            h = _mix_u64(<uint64_t>k0) ^ _mix_u64(<uint64_t>k1 + <uint64_t>0x9e3779b97f4a7c15)
            pos = <Py_ssize_t>(h & mask)
            while table_used[pos] and (table_k0[pos] != k0 or table_k1[pos] != k1):
                pos = <Py_ssize_t>((pos + 1) & mask)
            if not table_used[pos]:
                table_used[pos] = 1
                table_k0[pos] = k0
                table_k1[pos] = k1
                row_counts[pos] = 0
                value_counts[pos] = 0
                sums[pos] = 0.0
                mins[pos] = 0.0
                maxs[pos] = 0.0
                has_value[pos] = 0
                used_count += 1

            row_counts[pos] += 1
            if has_values:
                value_ok = values_valid[i]
                if value_ok:
                    value = values[i]
                    value_counts[pos] += 1
                    sums[pos] += value
                    if not has_value[pos] or value < mins[pos]:
                        mins[pos] = value
                    if not has_value[pos] or value > maxs[pos]:
                        maxs[pos] = value
                    has_value[pos] = 1

        out_k0 = np.empty(used_count, dtype=np.int64)
        out_k1 = np.empty(used_count, dtype=np.int64)
        out_row_counts = np.empty(used_count, dtype=np.int64)
        out_value_counts = np.empty(used_count, dtype=np.int64)
        out_sums = np.empty(used_count, dtype=np.float64)
        out_mins = np.empty(used_count, dtype=np.float64)
        out_maxs = np.empty(used_count, dtype=np.float64)
        out_has_value = np.empty(used_count, dtype=bool)

        out_pos = 0
        for pos in range(cap):
            if not table_used[pos]:
                continue
            out_k0[out_pos] = table_k0[pos]
            out_k1[out_pos] = table_k1[pos]
            out_row_counts[out_pos] = row_counts[pos]
            out_value_counts[out_pos] = value_counts[pos]
            out_sums[out_pos] = sums[pos]
            out_mins[out_pos] = mins[pos]
            out_maxs[out_pos] = maxs[pos]
            out_has_value[out_pos] = has_value[pos]
            out_pos += 1
        return out_k0, out_k1, out_row_counts, out_value_counts, out_sums, out_mins, out_maxs, out_has_value
    finally:
        free(table_k0); free(table_k1); free(table_used); free(row_counts); free(value_counts)
        free(sums); free(mins); free(maxs); free(has_value)
