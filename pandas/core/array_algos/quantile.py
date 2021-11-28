from __future__ import annotations

import numpy as np

from pandas._typing import (
    ArrayLike,
    npt,
)

from pandas.core.dtypes.missing import (
    isna,
    na_value_for_dtype,
)

from pandas.core.nanops import nanpercentile


def quantile_compat(
    values: ArrayLike, qs: npt.NDArray[np.float64], interpolation: str
) -> ArrayLike:
    """
    Compute the quantiles of the given values for each quantile in `qs`.

    Parameters
    ----------
    values : np.ndarray or ExtensionArray
    qs : np.ndarray[float64]
    interpolation : str

    Returns
    -------
    np.ndarray or ExtensionArray
    """
    if isinstance(values, np.ndarray):
        fill_value = na_value_for_dtype(values.dtype, compat=False)
        mask = isna(values)
        return quantile_with_mask(values, mask, fill_value, qs, interpolation)
    else:
        return values._quantile(qs, interpolation)


def quantile_with_mask(
    values: np.ndarray,
    mask: np.ndarray,
    fill_value,
    qs: npt.NDArray[np.float64],
    interpolation: str,
) -> np.ndarray:
    """
    Compute the quantiles of the given values for each quantile in `qs`.

    Parameters
    ----------
    values : np.ndarray
        For ExtensionArray, this is _values_for_factorize()[0]
    mask : np.ndarray[bool]
        mask = isna(values)
        For ExtensionArray, this is computed before calling _value_for_factorize
    fill_value : Scalar
        The value to interpret fill NA entries with
        For ExtensionArray, this is _values_for_factorize()[1]
    qs : np.ndarray[float64]
    interpolation : str
        Type of interpolation

    Returns
    -------
    np.ndarray

    Notes
    -----
    Assumes values is already 2D.  For ExtensionArray this means np.atleast_2d
    has been called on _values_for_factorize()[0]

    Quantile is computed along axis=1.
    """
    assert values.ndim == 2

    is_empty = values.shape[1] == 0

    if is_empty:
        # create the array of na_values
        # 2d len(values) * len(qs)
        flat = np.array([fill_value] * len(qs))
        result = np.repeat(flat, len(values)).reshape(len(values), len(qs))
    else:
        # asarray needed for Sparse, see GH#24600
        result = nanpercentile(
            values,
            np.array(qs) * 100,
            na_value=fill_value,
            mask=mask,
            interpolation=interpolation,
        )

        result = np.array(result, copy=False)
        result = result.T

    return result
