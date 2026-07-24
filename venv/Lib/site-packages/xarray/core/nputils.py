from __future__ import annotations

import warnings
from collections.abc import Callable

import numpy as np
import pandas as pd
from packaging.version import Version

from xarray.compat.array_api_compat import get_array_namespace
from xarray.core.utils import is_duck_array, module_available
from xarray.namedarray import pycompat

# remove once numpy 2.0 is the oldest supported version
if module_available("numpy", minversion="2.0.0.dev0"):
    from numpy.lib.array_utils import (  # type: ignore[import-not-found,unused-ignore]
        normalize_axis_index,
    )
else:
    from numpy.core.multiarray import (  # type: ignore[attr-defined,no-redef,unused-ignore]
        normalize_axis_index,
    )

# remove once numpy 2.0 is the oldest supported version
try:
    from numpy.exceptions import RankWarning  # type: ignore[attr-defined,unused-ignore]
except ImportError:
    from numpy import RankWarning  # type: ignore[attr-defined,no-redef,unused-ignore]

from xarray.core.options import OPTIONS

try:
    import bottleneck as bn

    _BOTTLENECK_AVAILABLE = True
except ImportError:
    # use numpy methods instead
    bn = np
    _BOTTLENECK_AVAILABLE = False


def _select_along_axis(values, idx, axis):
    other_ind = np.ix_(*[np.arange(s) for s in idx.shape])
    sl = other_ind[:axis] + (idx,) + other_ind[axis:]
    return values[sl]


def nanfirst(values, axis, keepdims=False):
    if isinstance(axis, tuple):
        (axis,) = axis
    axis = normalize_axis_index(axis, values.ndim)
    idx_first = np.argmax(~pd.isnull(values), axis=axis)
    result = _select_along_axis(values, idx_first, axis)
    if keepdims:
        return np.expand_dims(result, axis=axis)
    else:
        return result


def nanlast(values, axis, keepdims=False):
    if isinstance(axis, tuple):
        (axis,) = axis
    axis = normalize_axis_index(axis, values.ndim)
    rev = (slice(None),) * axis + (slice(None, None, -1),)
    idx_last = -1 - np.argmax(~pd.isnull(values)[rev], axis=axis)
    result = _select_along_axis(values, idx_last, axis)
    if keepdims:
        return np.expand_dims(result, axis=axis)
    else:
        return result


def inverse_permutation(indices: np.ndarray, N: int | None = None) -> np.ndarray:
    """Return indices for an inverse permutation.

    Parameters
    ----------
    indices : 1D np.ndarray with dtype=int
        Integer positions to assign elements to.
    N : int, optional
        Size of the array

    Returns
    -------
    inverse_permutation : 1D np.ndarray with dtype=int
        Integer indices to take from the original array to create the
        permutation.
    """
    if N is None:
        N = len(indices)
    # use intp instead of int64 because of windows :(
    inverse_permutation = np.full(N, -1, dtype=np.intp)
    inverse_permutation[indices] = np.arange(len(indices), dtype=np.intp)
    return inverse_permutation


def _ensure_bool_is_ndarray(result, *args):
    # numpy will sometimes return a scalar value from binary comparisons if it
    # can't handle the comparison instead of broadcasting, e.g.,
    # In [10]: 1 == np.array(['a', 'b'])
    # Out[10]: False
    # This function ensures that the result is the appropriate shape in these
    # cases
    if isinstance(result, bool):
        shape = np.broadcast(*args).shape
        constructor = np.ones if result else np.zeros
        result = constructor(shape, dtype=bool)
    return result


def array_eq(self, other):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"elementwise comparison failed")
        return _ensure_bool_is_ndarray(self == other, self, other)


def array_ne(self, other):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"elementwise comparison failed")
        return _ensure_bool_is_ndarray(self != other, self, other)


def _is_contiguous(positions):
    """Given a non-empty list, does it consist of contiguous integers?"""
    previous = positions[0]
    for current in positions[1:]:
        if current != previous + 1:
            return False
        previous = current
    return True


def _advanced_indexer_subspaces(key):
    """Indices of the advanced indexes subspaces for mixed indexing and vindex."""
    if not isinstance(key, tuple):
        key = (key,)
    advanced_index_positions = [
        i for i, k in enumerate(key) if not isinstance(k, slice)
    ]

    if not advanced_index_positions or not _is_contiguous(advanced_index_positions):
        # Nothing to reorder: dimensions on the indexing result are already
        # ordered like vindex. See NumPy's rule for "Combining advanced and
        # basic indexing":
        # https://numpy.org/doc/stable/reference/arrays.indexing.html#combining-advanced-and-basic-indexing
        return (), ()

    non_slices = [k for k in key if not isinstance(k, slice)]
    broadcasted_shape = np.broadcast_shapes(
        *[item.shape if is_duck_array(item) else (0,) for item in non_slices]
    )
    ndim = len(broadcasted_shape)
    mixed_positions = advanced_index_positions[0] + np.arange(ndim)
    vindex_positions = np.arange(ndim)
    return mixed_positions, vindex_positions


class NumpyVIndexAdapter:
    """Object that implements indexing like vindex on an np.ndarray.

    This is a pure Python implementation of (some of) the logic in this NumPy
    proposal: https://github.com/numpy/numpy/pull/6256
    """

    def __init__(self, array):
        self._array = array

    def __getitem__(self, key):
        mixed_positions, vindex_positions = _advanced_indexer_subspaces(key)
        return np.moveaxis(self._array[key], mixed_positions, vindex_positions)

    def __setitem__(self, key, value):
        """Value must have dimensionality matching the key."""
        mixed_positions, vindex_positions = _advanced_indexer_subspaces(key)
        self._array[key] = np.moveaxis(value, vindex_positions, mixed_positions)


def _create_method(name, npmodule=np) -> Callable:
    def f(values, axis=None, **kwargs):
        dtype = kwargs.get("dtype")
        bn_func = getattr(bn, name, None)

        xp = get_array_namespace(values)
        if xp is not np:
            func = getattr(xp, name, None)
            if func is not None:
                return func(values, axis=axis, **kwargs)
        if (
            module_available("numbagg")
            and OPTIONS["use_numbagg"]
            and isinstance(values, np.ndarray)
            # numbagg<0.7.0 uses ddof=1 only, but numpy uses ddof=0 by default
            and (
                pycompat.mod_version("numbagg") >= Version("0.7.0")
                or ("var" not in name and "std" not in name)
                or kwargs.get("ddof", 0) == 1
            )
            # TODO: bool?
            and values.dtype.kind in "uif"
            # and values.dtype.isnative
            and (dtype is None or np.dtype(dtype) == values.dtype)
            # numbagg.nanquantile only available after 0.8.0 and with linear method
            and (
                name != "nanquantile"
                or (
                    pycompat.mod_version("numbagg") >= Version("0.8.0")
                    and kwargs.get("method", "linear") == "linear"
                )
            )
        ):
            import numbagg  # type: ignore[import-not-found, unused-ignore]

            nba_func = getattr(numbagg, name, None)
            if nba_func is not None:
                # numbagg does not use dtype
                kwargs.pop("dtype", None)
                # prior to 0.7.0, numbagg did not support ddof; we ensure it's limited
                # to ddof=1 above.
                if pycompat.mod_version("numbagg") < Version("0.7.0"):
                    kwargs.pop("ddof", None)
                if name == "nanquantile":
                    kwargs["quantiles"] = kwargs.pop("q")
                    kwargs.pop("method", None)
                return nba_func(values, axis=axis, **kwargs)
        if (
            _BOTTLENECK_AVAILABLE
            and OPTIONS["use_bottleneck"]
            and isinstance(values, np.ndarray)
            and bn_func is not None
            and not isinstance(axis, tuple)
            and values.dtype.kind in "uifc"
            and values.dtype.isnative
            and (dtype is None or np.dtype(dtype) == values.dtype)
        ):
            # bottleneck does not take care dtype, min_count
            kwargs.pop("dtype", None)
            result = bn_func(values, axis=axis, **kwargs)
            # bottleneck returns python scalars for reduction over all axes
            if isinstance(result, float):
                result = np.float64(result)
        else:
            result = getattr(npmodule, name)(values, axis=axis, **kwargs)

        return result

    f.__name__ = name
    return f


def _nanpolyfit_1d(arr, x, rcond=None):
    out = np.full((x.shape[1] + 1,), np.nan)
    mask = np.isnan(arr)
    if not np.all(mask):
        out[:-1], resid, rank, _ = np.linalg.lstsq(x[~mask, :], arr[~mask], rcond=rcond)
        out[-1] = resid[0] if resid.size > 0 else np.nan
        warn_on_deficient_rank(rank, x.shape[1])
    return out


def warn_on_deficient_rank(rank, order):
    if rank != order:
        warnings.warn("Polyfit may be poorly conditioned", RankWarning, stacklevel=2)


def least_squares(lhs, rhs, rcond=None, skipna=False):
    if rhs.ndim > 2:
        out_shape = rhs.shape
        rhs = rhs.reshape(rhs.shape[0], -1)
    else:
        out_shape = None

    if skipna:
        added_dim = rhs.ndim == 1
        if added_dim:
            rhs = rhs.reshape(rhs.shape[0], 1)
        nan_cols = np.any(np.isnan(rhs), axis=0)
        out = np.empty((lhs.shape[1] + 1, rhs.shape[1]))
        if np.any(nan_cols):
            out[:, nan_cols] = np.apply_along_axis(
                _nanpolyfit_1d, 0, rhs[:, nan_cols], lhs
            )
        if np.any(~nan_cols):
            out[:-1, ~nan_cols], resids, rank, _ = np.linalg.lstsq(
                lhs, rhs[:, ~nan_cols], rcond=rcond
            )
            out[-1, ~nan_cols] = resids if resids.size > 0 else np.nan
            warn_on_deficient_rank(rank, lhs.shape[1])
        coeffs = out[:-1, :]
        residuals = out[-1, :]
        if added_dim:
            coeffs = coeffs.reshape(coeffs.shape[0])
            residuals = residuals.reshape(residuals.shape[0])
    else:
        coeffs, residuals, rank, _ = np.linalg.lstsq(lhs, rhs, rcond=rcond)
        if residuals.size == 0:
            residuals = coeffs[0] * np.nan
        warn_on_deficient_rank(rank, lhs.shape[1])

    if out_shape is not None:
        coeffs = coeffs.reshape(-1, *out_shape[1:])
        residuals = residuals.reshape(*out_shape[1:])
    return coeffs, residuals


nanmin = _create_method("nanmin")
nanmax = _create_method("nanmax")
nanmean = _create_method("nanmean")
nanmedian = _create_method("nanmedian")
nanvar = _create_method("nanvar")
nanstd = _create_method("nanstd")
nanprod = _create_method("nanprod")
nancumsum = _create_method("nancumsum")
nancumprod = _create_method("nancumprod")
nanargmin = _create_method("nanargmin")
nanargmax = _create_method("nanargmax")
nanquantile = _create_method("nanquantile")
