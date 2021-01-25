from typing import Callable

import numpy as np

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
    min_count: int = 0,
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
    try:
        fill_value = {
            np.cumprod: 1,
            np.maximum.accumulate: values.min(),
            np.cumsum: 0,
            np.minimum.accumulate: values.max(),
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
