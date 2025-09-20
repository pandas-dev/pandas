from __future__ import annotations

import warnings

import numpy as np

from xarray.core import dtypes, duck_array_ops, nputils, utils
from xarray.core.duck_array_ops import (
    astype,
    count,
    fillna,
    isnull,
    sum_where,
    where,
    where_method,
)


def _maybe_null_out(result, axis, mask, min_count=1):
    """
    xarray version of pandas.core.nanops._maybe_null_out
    """
    if axis is not None and getattr(result, "ndim", False):
        null_mask = (
            np.take(mask.shape, axis).prod()
            - duck_array_ops.sum(mask, axis)
            - min_count
        ) < 0
        dtype, fill_value = dtypes.maybe_promote(result.dtype)
        result = where(null_mask, fill_value, astype(result, dtype))

    elif getattr(result, "dtype", None) not in dtypes.NAT_TYPES:
        null_mask = mask.size - duck_array_ops.sum(mask)
        result = where(null_mask < min_count, np.nan, result)

    return result


def _nan_argminmax_object(func, fill_value, value, axis=None, **kwargs):
    """In house nanargmin, nanargmax for object arrays. Always return integer
    type
    """
    valid_count = count(value, axis=axis)
    value = fillna(value, fill_value)
    data = getattr(np, func)(value, axis=axis, **kwargs)

    # TODO This will evaluate dask arrays and might be costly.
    if duck_array_ops.array_any(valid_count == 0):
        raise ValueError("All-NaN slice encountered")

    return data


def _nan_minmax_object(func, fill_value, value, axis=None, **kwargs):
    """In house nanmin and nanmax for object array"""
    valid_count = count(value, axis=axis)
    filled_value = fillna(value, fill_value)
    data = getattr(np, func)(filled_value, axis=axis, **kwargs)
    if not hasattr(data, "dtype"):  # scalar case
        data = fill_value if valid_count == 0 else data
        # we've computed a single min, max value of type object.
        # don't let np.array turn a tuple back into an array
        return utils.to_0d_object_array(data)
    return where_method(data, valid_count != 0)


def nanmin(a, axis=None, out=None):
    if a.dtype.kind == "O":
        return _nan_minmax_object("min", dtypes.get_pos_infinity(a.dtype), a, axis)

    return nputils.nanmin(a, axis=axis)


def nanmax(a, axis=None, out=None):
    if a.dtype.kind == "O":
        return _nan_minmax_object("max", dtypes.get_neg_infinity(a.dtype), a, axis)

    return nputils.nanmax(a, axis=axis)


def nanargmin(a, axis=None):
    if a.dtype.kind == "O":
        fill_value = dtypes.get_pos_infinity(a.dtype)
        return _nan_argminmax_object("argmin", fill_value, a, axis=axis)

    return nputils.nanargmin(a, axis=axis)


def nanargmax(a, axis=None):
    if a.dtype.kind == "O":
        fill_value = dtypes.get_neg_infinity(a.dtype)
        return _nan_argminmax_object("argmax", fill_value, a, axis=axis)

    return nputils.nanargmax(a, axis=axis)


def nansum(a, axis=None, dtype=None, out=None, min_count=None):
    mask = isnull(a)
    result = sum_where(a, axis=axis, dtype=dtype, where=mask)
    if min_count is not None:
        return _maybe_null_out(result, axis, mask, min_count)
    else:
        return result


def _nanmean_ddof_object(ddof, value, axis=None, dtype=None, **kwargs):
    """In house nanmean. ddof argument will be used in _nanvar method"""
    valid_count = count(value, axis=axis)
    value = fillna(value, 0)
    # As dtype inference is impossible for object dtype, we assume float
    # https://github.com/dask/dask/issues/3162
    if dtype is None and value.dtype.kind == "O":
        dtype = float

    data = np.sum(value, axis=axis, dtype=dtype, **kwargs)
    data = data / (valid_count - ddof)
    return where_method(data, valid_count != 0)


def nanmean(a, axis=None, dtype=None, out=None):
    if a.dtype.kind == "O":
        return _nanmean_ddof_object(0, a, axis=axis, dtype=dtype)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", r"Mean of empty slice", category=RuntimeWarning
        )

        return nputils.nanmean(a, axis=axis, dtype=dtype)


def nanmedian(a, axis=None, out=None):
    # The dask algorithm works by rechunking to one chunk along axis
    # Make sure we trigger the dask error when passing all dimensions
    # so that we don't rechunk the entire array to one chunk and
    # possibly blow memory
    if axis is not None and len(np.atleast_1d(axis)) == a.ndim:
        axis = None
    return nputils.nanmedian(a, axis=axis)


def _nanvar_object(value, axis=None, ddof=0, keepdims=False, **kwargs):
    value_mean = _nanmean_ddof_object(
        ddof=0, value=value, axis=axis, keepdims=True, **kwargs
    )
    squared = (astype(value, value_mean.dtype) - value_mean) ** 2
    return _nanmean_ddof_object(ddof, squared, axis=axis, keepdims=keepdims, **kwargs)


def nanvar(a, axis=None, dtype=None, out=None, ddof=0):
    if a.dtype.kind == "O":
        return _nanvar_object(a, axis=axis, dtype=dtype, ddof=ddof)

    return nputils.nanvar(a, axis=axis, dtype=dtype, ddof=ddof)


def nanstd(a, axis=None, dtype=None, out=None, ddof=0):
    return nputils.nanstd(a, axis=axis, dtype=dtype, ddof=ddof)


def nanprod(a, axis=None, dtype=None, out=None, min_count=None):
    mask = isnull(a)
    result = nputils.nanprod(a, axis=axis, dtype=dtype)
    if min_count is not None:
        return _maybe_null_out(result, axis, mask, min_count)
    else:
        return result


def nancumsum(a, axis=None, dtype=None, out=None):
    return nputils.nancumsum(a, axis=axis, dtype=dtype)


def nancumprod(a, axis=None, dtype=None, out=None):
    return nputils.nancumprod(a, axis=axis, dtype=dtype)
