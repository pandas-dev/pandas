"""
datetimelke_accumulations.py is for accumulations of datetimelike extension arrays
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas._libs import iNaT
from pandas._libs.tslibs import OutOfBoundsTimedelta

from pandas.core.dtypes.missing import isna

if TYPE_CHECKING:
    from collections.abc import Callable


def _check_cumsum_overflow(
    values: np.ndarray, result: np.ndarray, mask: np.ndarray
) -> None:
    """
    Raise if a running total wrapped around or landed on the NaT sentinel.

    Parameters
    ----------
    values : np.ndarray[int64]
        The addends, with NA positions already zeroed out.
    result : np.ndarray[int64]
        The running totals, before NA positions are set to NaT.
    mask : np.ndarray[bool]
        Positions whose result will be NaT regardless, and so are exempt.
    """
    # Whether each running total is greater than the one before it, taking the
    #  total before the first entry to be zero.
    stepped_up = np.empty(result.shape, dtype=bool)
    np.greater(result[:1], 0, out=stepped_up[:1])
    np.greater(result[1:], result[:-1], out=stepped_up[1:])

    # Absent a signed wrap, the total goes up exactly when the addend is
    #  positive; a wrap flips the direction of the step.
    invalid = (values > 0) != stepped_up
    # A total of exactly int64.min does not wrap, but is indistinguishable
    #  from NaT once stored.
    invalid |= result == iNaT
    invalid &= ~mask

    if invalid.any():
        raise OutOfBoundsTimedelta("overflow in timedelta operation")


def _cum_func(
    func: Callable,
    values: np.ndarray,
    *,
    skipna: bool = True,
) -> np.ndarray:
    """
    Accumulations for 1D datetimelike arrays.

    Parameters
    ----------
    func : np.cumsum, np.maximum.accumulate, np.minimum.accumulate
    values : np.ndarray
        Numpy array with the values (can be of any dtype that support the
        operation). Values is changed is modified inplace.
    skipna : bool, default True
        Whether to skip NA.
    """
    try:
        fill_value = {
            np.maximum.accumulate: np.iinfo(np.int64).min,
            np.cumsum: 0,
            np.minimum.accumulate: np.iinfo(np.int64).max,
        }[func]
    except KeyError as err:
        raise ValueError(
            f"No accumulation for {func} implemented on BaseMaskedArray"
        ) from err

    mask = isna(values)
    y = values.view("i8")
    y[mask] = fill_value

    if not skipna:
        mask = np.maximum.accumulate(mask)

    # GH 57956
    result = func(y, axis=0)
    if func is np.cumsum:
        # GH#66551: cummin/cummax cannot leave the range, cumsum can
        _check_cumsum_overflow(y, result, mask)
    result[mask] = iNaT

    if values.dtype.kind in "mM":
        return result.view(values.dtype.base)
    return result


def cumsum(values: np.ndarray, *, skipna: bool = True) -> np.ndarray:
    return _cum_func(np.cumsum, values, skipna=skipna)


def cummin(values: np.ndarray, *, skipna: bool = True) -> np.ndarray:
    return _cum_func(np.minimum.accumulate, values, skipna=skipna)


def cummax(values: np.ndarray, *, skipna: bool = True) -> np.ndarray:
    return _cum_func(np.maximum.accumulate, values, skipna=skipna)
