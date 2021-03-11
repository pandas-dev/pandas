from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas._libs import lib
from pandas._typing import ArrayLike

from pandas.core.dtypes.common import (
    is_list_like,
    is_sparse,
)
from pandas.core.dtypes.missing import (
    isna,
    na_value_for_dtype,
)

from pandas.core.nanops import nanpercentile

if TYPE_CHECKING:
    from pandas.core.arrays import ExtensionArray


def quantile_compat(values: ArrayLike, qs, interpolation: str, axis: int) -> ArrayLike:
    """
    Compute the quantiles of the given values for each quantile in `qs`.

    Parameters
    ----------
    values : np.ndarray or ExtensionArray
    qs : a scalar or list of the quantiles to be computed
    interpolation : str
    axis : int

    Returns
    -------
    np.ndarray or ExtensionArray
    """
    if isinstance(values, np.ndarray):
        fill_value = na_value_for_dtype(values.dtype, compat=False)
        mask = isna(values)
        return quantile_with_mask(values, mask, fill_value, qs, interpolation, axis)
    else:
        return quantile_ea_compat(values, qs, interpolation, axis)


def quantile_with_mask(
    values: np.ndarray,
    mask: np.ndarray,
    fill_value,
    qs,
    interpolation: str,
    axis: int,
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
    qs : a scalar or list of the quantiles to be computed
    interpolation : str
        Type of interpolation
    axis : int
        Axis along which to compute quantiles.

    Returns
    -------
    np.ndarray

    Notes
    -----
    Assumes values is already 2D.  For ExtensionArray this means np.atleast_2d
    has been called on _values_for_factorize()[0]
    """
    is_empty = values.shape[axis] == 0
    orig_scalar = not is_list_like(qs)
    if orig_scalar:
        # make list-like, unpack later
        qs = [qs]

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
            axis=axis,
            na_value=fill_value,
            mask=mask,
            ndim=values.ndim,
            interpolation=interpolation,
        )

        result = np.array(result, copy=False)
        result = result.T

    if orig_scalar:
        assert result.shape[-1] == 1, result.shape
        result = result[..., 0]
        result = lib.item_from_zerodim(result)

    return result


def quantile_ea_compat(
    values: ExtensionArray, qs, interpolation: str, axis: int
) -> ExtensionArray:
    """
    ExtensionArray compatibility layer for quantile_with_mask.

    We pretend that an ExtensionArray with shape (N,) is actually (1, N,)
    for compatibility with non-EA code.

    Parameters
    ----------
    values : ExtensionArray
    qs : a scalar or list of the quantiles to be computed
    interpolation: str
    axis : int

    Returns
    -------
    ExtensionArray
    """
    # TODO(EA2D): make-believe not needed with 2D EAs
    orig = values

    # asarray needed for Sparse, see GH#24600
    mask = np.asarray(values.isna())
    mask = np.atleast_2d(mask)

    # error: Incompatible types in assignment (expression has type "ndarray", variable
    # has type "ExtensionArray")
    values, fill_value = values._values_for_factorize()  # type: ignore[assignment]
    # error: No overload variant of "atleast_2d" matches argument type "ExtensionArray"
    values = np.atleast_2d(values)  # type: ignore[call-overload]

    # error: Argument 1 to "quantile_with_mask" has incompatible type "ExtensionArray";
    # expected "ndarray"
    result = quantile_with_mask(
        values, mask, fill_value, qs, interpolation, axis  # type: ignore[arg-type]
    )

    if not is_sparse(orig.dtype):
        # shape[0] should be 1 as long as EAs are 1D

        if result.ndim == 1:
            # i.e. qs was originally a scalar
            assert result.shape == (1,), result.shape
            result = type(orig)._from_factorized(result, orig)

        else:
            assert result.shape == (1, len(qs)), result.shape
            result = type(orig)._from_factorized(result[0], orig)

    # error: Incompatible return value type (got "ndarray", expected "ExtensionArray")
    return result  # type: ignore[return-value]
