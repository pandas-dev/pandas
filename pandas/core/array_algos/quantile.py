from __future__ import annotations

from functools import partial
from typing import Literal

import numpy as np

from pandas._libs import groupby as libgroupby
from pandas._typing import (
    ArrayLike,
    Scalar,
    npt,
)
from pandas.compat.numpy import np_percentile_argname

from pandas.core.dtypes.missing import (
    isna,
    na_value_for_dtype,
)


def groupby_quantile_ndim_compat(
    *,
    qs: npt.NDArray[np.float64],
    interpolation: Literal["linear", "lower", "higher", "nearest", "midpoint"],
    ngroups: int,
    ids: npt.NDArray[np.intp],
    labels_for_lexsort: npt.NDArray[np.intp],
    npy_vals: np.ndarray,
    mask: npt.NDArray[np.bool_],
    result_mask: npt.NDArray[np.bool_] | None,
) -> np.ndarray:
    """
    Compatibility layer to handle either 1D arrays or 2D ndarrays in
    GroupBy.quantile. Located here to be available to
    ExtensionArray.groupby_quantile for dispatching after casting to numpy.

    Parameters
    ----------
    qs : np.ndarray[float64]
        Values between 0 and 1 providing the quantile(s) to compute.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        Method to use when the desired quantile falls between two points.
    ngroups : int
        The number of groupby groups.
    ids : np.ndarray[intp]
        Group labels.
    labels_for_lexsort : np.ndarray[intp]
        Group labels, but with -1s moved moved to the end to sort NAs last.
    npy_vals : np.ndarray
        The values for which we are computing quantiles.
    mask : np.ndarray[bool]
        Locations to treat as NA.
    result_mask : np.ndarray[bool] or None
        If present, set to True for locations that should be treated as missing
        a result.  Modified in-place.

    Returns
    -------
    np.ndarray
    """
    nqs = len(qs)

    ncols = 1
    if npy_vals.ndim == 2:
        ncols = npy_vals.shape[0]
        shaped_labels = np.broadcast_to(
            labels_for_lexsort, (ncols, len(labels_for_lexsort))
        )
    else:
        shaped_labels = labels_for_lexsort

    npy_out = np.empty((ncols, ngroups, nqs), dtype=np.float64)

    # Get an index of values sorted by values and then labels
    order = (npy_vals, shaped_labels)
    sort_arr = np.lexsort(order).astype(np.intp, copy=False)

    func = partial(
        libgroupby.group_quantile, labels=ids, qs=qs, interpolation=interpolation
    )

    if npy_vals.ndim == 1:
        func(
            npy_out[0],
            values=npy_vals,
            mask=mask,
            sort_indexer=sort_arr,
            result_mask=result_mask,
        )
    else:
        # if we ever did get here with non-None result_mask, we'd pass result_mask[i]
        assert result_mask is None
        for i in range(ncols):
            func(
                npy_out[i],
                values=npy_vals[i],
                mask=mask[i],
                sort_indexer=sort_arr[i],
            )

    if npy_vals.ndim == 1:
        npy_out = npy_out.reshape(ngroups * nqs)
    else:
        npy_out = npy_out.reshape(ncols, ngroups * nqs)

    return npy_out


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
    mask: npt.NDArray[np.bool_],
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
    assert values.shape == mask.shape
    if values.ndim == 1:
        # unsqueeze, operate, re-squeeze
        values = np.atleast_2d(values)
        mask = np.atleast_2d(mask)
        res_values = quantile_with_mask(values, mask, fill_value, qs, interpolation)
        return res_values[0]

    assert values.ndim == 2

    is_empty = values.shape[1] == 0

    if is_empty:
        # create the array of na_values
        # 2d len(values) * len(qs)
        flat = np.array([fill_value] * len(qs))
        result = np.repeat(flat, len(values)).reshape(len(values), len(qs))
    else:
        result = _nanpercentile(
            values,
            qs * 100.0,
            na_value=fill_value,
            mask=mask,
            interpolation=interpolation,
        )

        result = np.array(result, copy=False)
        result = result.T

    return result


def _nanpercentile_1d(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    qs: npt.NDArray[np.float64],
    na_value: Scalar,
    interpolation: str,
) -> Scalar | np.ndarray:
    """
    Wrapper for np.percentile that skips missing values, specialized to
    1-dimensional case.

    Parameters
    ----------
    values : array over which to find quantiles
    mask : ndarray[bool]
        locations in values that should be considered missing
    qs : np.ndarray[float64] of quantile indices to find
    na_value : scalar
        value to return for empty or all-null values
    interpolation : str

    Returns
    -------
    quantiles : scalar or array
    """
    # mask is Union[ExtensionArray, ndarray]
    values = values[~mask]

    if len(values) == 0:
        # Can't pass dtype=values.dtype here bc we might have na_value=np.nan
        #  with values.dtype=int64 see test_quantile_empty
        # equiv: 'np.array([na_value] * len(qs))' but much faster
        return np.full(len(qs), na_value)

    return np.percentile(
        values,
        qs,
        # error: No overload variant of "percentile" matches argument
        # types "ndarray[Any, Any]", "ndarray[Any, dtype[floating[_64Bit]]]"
        # , "Dict[str, str]"  [call-overload]
        **{np_percentile_argname: interpolation},  # type: ignore[call-overload]
    )


def _nanpercentile(
    values: np.ndarray,
    qs: npt.NDArray[np.float64],
    *,
    na_value,
    mask: npt.NDArray[np.bool_],
    interpolation: str,
):
    """
    Wrapper for np.percentile that skips missing values.

    Parameters
    ----------
    values : np.ndarray[ndim=2]  over which to find quantiles
    qs : np.ndarray[float64] of quantile indices to find
    na_value : scalar
        value to return for empty or all-null values
    mask : np.ndarray[bool]
        locations in values that should be considered missing
    interpolation : str

    Returns
    -------
    quantiles : scalar or array
    """

    if values.dtype.kind in ["m", "M"]:
        # need to cast to integer to avoid rounding errors in numpy
        result = _nanpercentile(
            values.view("i8"),
            qs=qs,
            na_value=na_value.view("i8"),
            mask=mask,
            interpolation=interpolation,
        )

        # Note: we have to do `astype` and not view because in general we
        #  have float result at this point, not i8
        return result.astype(values.dtype)

    if mask.any():
        # Caller is responsible for ensuring mask shape match
        assert mask.shape == values.shape
        result = [
            _nanpercentile_1d(val, m, qs, na_value, interpolation=interpolation)
            for (val, m) in zip(list(values), list(mask))
        ]
        if values.dtype.kind == "f":
            # preserve itemsize
            result = np.array(result, dtype=values.dtype, copy=False).T
        else:
            result = np.array(result, copy=False).T
            if (
                result.dtype != values.dtype
                and not mask.all()
                and (result == result.astype(values.dtype, copy=False)).all()
            ):
                # mask.all() will never get cast back to int
                # e.g. values id integer dtype and result is floating dtype,
                #  only cast back to integer dtype if result values are all-integer.
                result = result.astype(values.dtype, copy=False)
        return result
    else:
        return np.percentile(
            values,
            qs,
            axis=1,
            # error: No overload variant of "percentile" matches argument types
            # "ndarray[Any, Any]", "ndarray[Any, dtype[floating[_64Bit]]]",
            # "int", "Dict[str, str]"  [call-overload]
            **{np_percentile_argname: interpolation},  # type: ignore[call-overload]
        )
