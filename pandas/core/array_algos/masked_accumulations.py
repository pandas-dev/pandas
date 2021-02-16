from typing import Callable

import numpy as np

from pandas.core.dtypes.common import (
    is_float_dtype,
    is_integer_dtype,
)

"""
masked_accumulations.py is for accumulation algorithms using a mask-based approach
for missing values.
"""


def _cum_func(
    func: Callable,
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
):
    """
    Accumulations for 1D masked array.

    Parameters
    ----------
    func : np.cumsum, np.cumprod, np.maximum.accumulate, np.minimum.accumulate
    values : np.ndarray
        Numpy array with the values (can be of any dtype that support the
        operation).
    mask : np.ndarray
        Boolean numpy array (True values indicate missing values).
    skipna : bool, default True
        Whether to skip NA.
    """
    dtype_info = None
    if is_float_dtype(values):
        dtype_info = np.finfo(values.dtype.type)
    elif is_integer_dtype(values):
        dtype_info = np.iinfo(values.dtype.type)
    else:
        raise NotImplementedError(
            f"No masked accumulation defined for dtype {values.dtype.type}"
        )
    try:
        fill_value = {
            np.cumprod: 1,
            np.maximum.accumulate: dtype_info.min,
            np.cumsum: 0,
            np.minimum.accumulate: dtype_info.max,
        }[func]
    except KeyError:
        raise ValueError(f"No accumulation for {func} implemented on BaseMaskedArray")

    values[mask] = fill_value

    if not skipna:
        mask = np.maximum.accumulate(mask)

    values = func(values)
    return values, mask


def cumsum(values: np.ndarray, mask: np.ndarray, *, skipna: bool = True):
    return _cum_func(np.cumsum, values, mask, skipna=skipna)


def cumprod(values: np.ndarray, mask: np.ndarray, *, skipna: bool = True):
    return _cum_func(np.cumprod, values, mask, skipna=skipna)


def cummin(values: np.ndarray, mask: np.ndarray, *, skipna: bool = True):
    return _cum_func(np.minimum.accumulate, values, mask, skipna=skipna)


def cummax(values: np.ndarray, mask: np.ndarray, *, skipna: bool = True):
    return _cum_func(np.maximum.accumulate, values, mask, skipna=skipna)
