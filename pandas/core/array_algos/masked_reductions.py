"""
masked_reductions.py is for reduction algorithms using a mask-based approach
for missing values.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
import warnings

import numpy as np

from pandas._libs import (
    lib,
    missing as libmissing,
)

from pandas.core.nanops import check_below_min_count

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import (
        AxisInt,
        npt,
    )


def _get_min_value(dtype: np.dtype) -> int | float:
    """Get the minimum value for a given numpy dtype."""
    if dtype.kind == "f":
        return -np.inf
    elif dtype.kind == "u":
        return 0
    else:
        return np.iinfo(dtype).min


def _get_max_value(dtype: np.dtype) -> int | float:
    """Get the maximum value for a given numpy dtype."""
    if dtype.kind == "f":
        return np.inf
    else:
        return np.iinfo(dtype).max


def _reductions(
    func: Callable,
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    min_count: int = 0,
    axis: AxisInt | None = None,
    initial: object | lib.NoDefault = lib.no_default,
    **kwargs,
):
    """
    Sum, mean or product for 1D masked array.

    Parameters
    ----------
    func : np.sum or np.prod
    values : np.ndarray
        Numpy array with the values (can be of any dtype that support the
        operation).
    mask : np.ndarray[bool]
        Boolean numpy array (True values indicate missing values).
    skipna : bool, default True
        Whether to skip NA.
    min_count : int, default 0
        The required number of valid values to perform the operation. If fewer than
        ``min_count`` non-NA values are present the result will be NA.
    axis : int, optional, default None
    initial : scalar, optional
        Starting value for the reduction. NumPy has a default value for most
        data types, but for object-dtype arrays we need to specify it explicitly
    """
    if initial is not lib.no_default:
        kwargs["initial"] = initial

    if not skipna:
        if axis is not None and values.ndim > 1:
            # For 2D with axis, compute the reduction and let the caller
            # handle per-row/column masking via _wrap_reduction_result.
            return func(values, axis=axis, **kwargs)
        if mask.any() or check_below_min_count(values.shape, None, min_count):
            return libmissing.NA
        else:
            return func(values, axis=axis, **kwargs)
    else:
        if check_below_min_count(values.shape, mask, min_count) and (
            axis is None or values.ndim == 1
        ):
            return libmissing.NA

        return func(values, where=~mask, axis=axis, **kwargs)


def sum(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    min_count: int = 0,
    axis: AxisInt | None = None,
    initial: object | lib.NoDefault = lib.no_default,
):
    return _reductions(
        np.sum,
        values=values,
        mask=mask,
        skipna=skipna,
        min_count=min_count,
        axis=axis,
        initial=initial,
    )


def prod(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    min_count: int = 0,
    axis: AxisInt | None = None,
):
    return _reductions(
        np.prod, values=values, mask=mask, skipna=skipna, min_count=min_count, axis=axis
    )


def _minmax(
    func: Callable,
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    axis: AxisInt | None = None,
):
    """
    Reduction for 1D masked array.

    Parameters
    ----------
    func : np.min or np.max
    values : np.ndarray
        Numpy array with the values (can be of any dtype that support the
        operation).
    mask : np.ndarray[bool]
        Boolean numpy array (True values indicate missing values).
    skipna : bool, default True
        Whether to skip NA.
    axis : int, optional, default None
    """
    if not skipna:
        if axis is not None and values.ndim > 1:
            # For 2D with axis, compute the reduction and let the caller
            # handle per-row/column masking via _wrap_reduction_result.
            return func(values, axis=axis)
        if mask.any() or not values.size:
            # min/max with empty array raise in numpy, pandas returns NA
            return libmissing.NA
        else:
            return func(values, axis=axis)
    else:
        if axis is not None and values.ndim > 1:
            # Can't use values[~mask] which would flatten the array.
            # Instead, replace masked values with the identity element
            # and then reduce.
            if func is np.min:
                fill_value = _get_max_value(values.dtype)
            else:
                fill_value = _get_min_value(values.dtype)
            filled = np.where(mask, fill_value, values)
            return func(filled, axis=axis)
        subset = values[~mask]
        if subset.size:
            return func(subset, axis=axis)
        else:
            # min/max with empty array raise in numpy, pandas returns NA
            return libmissing.NA


def min(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    axis: AxisInt | None = None,
):
    return _minmax(np.min, values=values, mask=mask, skipna=skipna, axis=axis)


def max(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    axis: AxisInt | None = None,
):
    return _minmax(np.max, values=values, mask=mask, skipna=skipna, axis=axis)


def mean(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    axis: AxisInt | None = None,
):
    if not values.size or (mask.all() and (axis is None or values.ndim == 1)):
        return libmissing.NA
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _reductions(np.mean, values=values, mask=mask, skipna=skipna, axis=axis)


def var(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    axis: AxisInt | None = None,
    ddof: int = 1,
):
    if not values.size or (mask.all() and (axis is None or values.ndim == 1)):
        return libmissing.NA

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _reductions(
            np.var, values=values, mask=mask, skipna=skipna, axis=axis, ddof=ddof
        )


def std(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    *,
    skipna: bool = True,
    axis: AxisInt | None = None,
    ddof: int = 1,
):
    if not values.size or (mask.all() and (axis is None or values.ndim == 1)):
        return libmissing.NA

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _reductions(
            np.std, values=values, mask=mask, skipna=skipna, axis=axis, ddof=ddof
        )
