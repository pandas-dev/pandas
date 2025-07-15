from __future__ import annotations

import functools
from typing import Any

import numpy as np
from pandas.api.types import is_extension_array_dtype

from xarray.compat import array_api_compat, npcompat
from xarray.compat.npcompat import HAS_STRING_DTYPE
from xarray.core import utils

# Use as a sentinel value to indicate a dtype appropriate NA value.
NA = utils.ReprObject("<NA>")


@functools.total_ordering
class AlwaysGreaterThan:
    def __gt__(self, other):
        return True

    def __eq__(self, other):
        return isinstance(other, type(self))


@functools.total_ordering
class AlwaysLessThan:
    def __lt__(self, other):
        return True

    def __eq__(self, other):
        return isinstance(other, type(self))


# Equivalence to np.inf (-np.inf) for object-type
INF = AlwaysGreaterThan()
NINF = AlwaysLessThan()


# Pairs of types that, if both found, should be promoted to object dtype
# instead of following NumPy's own type-promotion rules. These type promotion
# rules match pandas instead. For reference, see the NumPy type hierarchy:
# https://numpy.org/doc/stable/reference/arrays.scalars.html
PROMOTE_TO_OBJECT: tuple[tuple[type[np.generic], type[np.generic]], ...] = (
    (np.number, np.character),  # numpy promotes to character
    (np.bool_, np.character),  # numpy promotes to character
    (np.bytes_, np.str_),  # numpy promotes to unicode
)


def maybe_promote(dtype: np.dtype) -> tuple[np.dtype, Any]:
    """Simpler equivalent of pandas.core.common._maybe_promote

    Parameters
    ----------
    dtype : np.dtype

    Returns
    -------
    dtype : Promoted dtype that can hold missing values.
    fill_value : Valid missing value for the promoted dtype.
    """
    # N.B. these casting rules should match pandas
    dtype_: np.typing.DTypeLike
    fill_value: Any
    if HAS_STRING_DTYPE and np.issubdtype(dtype, np.dtypes.StringDType()):
        # for now, we always promote string dtypes to object for consistency with existing behavior
        # TODO: refactor this once we have a better way to handle numpy vlen-string dtypes
        dtype_ = object
        fill_value = np.nan
    elif isdtype(dtype, "real floating"):
        dtype_ = dtype
        fill_value = np.nan
    elif np.issubdtype(dtype, np.timedelta64):
        # See https://github.com/numpy/numpy/issues/10685
        # np.timedelta64 is a subclass of np.integer
        # Check np.timedelta64 before np.integer
        fill_value = np.timedelta64("NaT")
        dtype_ = dtype
    elif isdtype(dtype, "integral"):
        dtype_ = np.float32 if dtype.itemsize <= 2 else np.float64
        fill_value = np.nan
    elif isdtype(dtype, "complex floating"):
        dtype_ = dtype
        fill_value = np.nan + np.nan * 1j
    elif np.issubdtype(dtype, np.datetime64):
        dtype_ = dtype
        fill_value = np.datetime64("NaT")
    else:
        dtype_ = object
        fill_value = np.nan

    dtype_out = np.dtype(dtype_)
    fill_value = dtype_out.type(fill_value)
    return dtype_out, fill_value


NAT_TYPES = {np.datetime64("NaT").dtype, np.timedelta64("NaT").dtype}


def get_fill_value(dtype):
    """Return an appropriate fill value for this dtype.

    Parameters
    ----------
    dtype : np.dtype

    Returns
    -------
    fill_value : Missing value corresponding to this dtype.
    """
    _, fill_value = maybe_promote(dtype)
    return fill_value


def get_pos_infinity(dtype, max_for_int=False):
    """Return an appropriate positive infinity for this dtype.

    Parameters
    ----------
    dtype : np.dtype
    max_for_int : bool
        Return np.iinfo(dtype).max instead of np.inf

    Returns
    -------
    fill_value : positive infinity value corresponding to this dtype.
    """
    if isdtype(dtype, "real floating"):
        return np.inf

    if isdtype(dtype, "integral"):
        if max_for_int:
            return np.iinfo(dtype).max
        else:
            return np.inf

    if isdtype(dtype, "complex floating"):
        return np.inf + 1j * np.inf

    if isdtype(dtype, "bool"):
        return True

    return np.array(INF, dtype=object)


def get_neg_infinity(dtype, min_for_int=False):
    """Return an appropriate positive infinity for this dtype.

    Parameters
    ----------
    dtype : np.dtype
    min_for_int : bool
        Return np.iinfo(dtype).min instead of -np.inf

    Returns
    -------
    fill_value : positive infinity value corresponding to this dtype.
    """
    if isdtype(dtype, "real floating"):
        return -np.inf

    if isdtype(dtype, "integral"):
        if min_for_int:
            return np.iinfo(dtype).min
        else:
            return -np.inf

    if isdtype(dtype, "complex floating"):
        return -np.inf - 1j * np.inf

    if isdtype(dtype, "bool"):
        return False

    return np.array(NINF, dtype=object)


def is_datetime_like(dtype) -> bool:
    """Check if a dtype is a subclass of the numpy datetime types"""
    return _is_numpy_subdtype(dtype, (np.datetime64, np.timedelta64))


def is_object(dtype) -> bool:
    """Check if a dtype is object"""
    return _is_numpy_subdtype(dtype, object)


def is_string(dtype) -> bool:
    """Check if a dtype is a string dtype"""
    return _is_numpy_subdtype(dtype, (np.str_, np.character))


def _is_numpy_subdtype(dtype, kind) -> bool:
    if not isinstance(dtype, np.dtype):
        return False

    kinds = kind if isinstance(kind, tuple) else (kind,)
    return any(np.issubdtype(dtype, kind) for kind in kinds)


def isdtype(dtype, kind: str | tuple[str, ...], xp=None) -> bool:
    """Compatibility wrapper for isdtype() from the array API standard.

    Unlike xp.isdtype(), kind must be a string.
    """
    # TODO(shoyer): remove this wrapper when Xarray requires
    # numpy>=2 and pandas extensions arrays are implemented in
    # Xarray via the array API
    if not isinstance(kind, str) and not (
        isinstance(kind, tuple) and all(isinstance(k, str) for k in kind)  # type: ignore[redundant-expr]
    ):
        raise TypeError(f"kind must be a string or a tuple of strings: {kind!r}")

    if isinstance(dtype, np.dtype):
        return npcompat.isdtype(dtype, kind)
    elif is_extension_array_dtype(dtype):
        # we never want to match pandas extension array dtypes
        return False
    else:
        if xp is None:
            xp = np
        return xp.isdtype(dtype, kind)


def preprocess_types(t):
    if isinstance(t, str | bytes):
        return type(t)
    elif isinstance(dtype := getattr(t, "dtype", t), np.dtype) and (
        np.issubdtype(dtype, np.str_) or np.issubdtype(dtype, np.bytes_)
    ):
        # drop the length from numpy's fixed-width string dtypes, it is better to
        # recalculate
        # TODO(keewis): remove once the minimum version of `numpy.result_type` does this
        # for us
        return dtype.type
    else:
        return t


def result_type(
    *arrays_and_dtypes: np.typing.ArrayLike | np.typing.DTypeLike,
    xp=None,
) -> np.dtype:
    """Like np.result_type, but with type promotion rules matching pandas.

    Examples of changed behavior:
    number + string -> object (not string)
    bytes + unicode -> object (not unicode)

    Parameters
    ----------
    *arrays_and_dtypes : list of arrays and dtypes
        The dtype is extracted from both numpy and dask arrays.

    Returns
    -------
    numpy.dtype for the result.
    """
    # TODO (keewis): replace `array_api_compat.result_type` with `xp.result_type` once we
    # can require a version of the Array API that supports passing scalars to it.
    from xarray.core.duck_array_ops import get_array_namespace

    if xp is None:
        xp = get_array_namespace(arrays_and_dtypes)

    types = {
        array_api_compat.result_type(preprocess_types(t), xp=xp)
        for t in arrays_and_dtypes
    }
    if any(isinstance(t, np.dtype) for t in types):
        # only check if there's numpy dtypes – the array API does not
        # define the types we're checking for
        for left, right in PROMOTE_TO_OBJECT:
            if any(np.issubdtype(t, left) for t in types) and any(
                np.issubdtype(t, right) for t in types
            ):
                return np.dtype(object)

    return array_api_compat.result_type(
        *map(preprocess_types, arrays_and_dtypes), xp=xp
    )
