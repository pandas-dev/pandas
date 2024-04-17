"""
categorical_accumulations.py is for accumulation algorithms using a mask-based
approach for missing values.
"""

from __future__ import annotations

from typing import Callable

import numpy as np


def _cum_func(
    func: Callable,
    values: np.ndarray,
    *,
    skipna: bool = True,
) -> np.ndarray:
    """
    Accumulations for 1D categorical arrays.

    We will modify values in place to replace NAs with the appropriate fill value.

    Parameters
    ----------
    func : np.maximum.accumulate, np.minimum.accumulate
    values : np.ndarray
        Numpy integer array with the values and with NAs being -1.
    skipna : bool, default True
        Whether to skip NA.
    """
    dtype_info = np.iinfo(values.dtype.type)
    try:
        fill_value = {
            np.maximum.accumulate: dtype_info.min,
            np.minimum.accumulate: dtype_info.max,
        }[func]
    except KeyError as err:
        raise NotImplementedError(
            f"No accumulation for {func} implemented on BaseMaskedArray"
        ) from err

    mask = values == -1
    values[mask] = fill_value

    if not skipna:
        mask = np.maximum.accumulate(mask)

    values = func(values)
    values[mask] = -1

    return values


def cummin(values: np.ndarray, *, skipna: bool = True) -> np.ndarray:
    return _cum_func(np.minimum.accumulate, values, skipna=skipna)


def cummax(values: np.ndarray, *, skipna: bool = True) -> np.ndarray:
    return _cum_func(np.maximum.accumulate, values, skipna=skipna)
