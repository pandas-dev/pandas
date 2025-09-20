from __future__ import annotations

import functools
from typing import Any, Literal, TypeGuard

import numpy as np

from xarray.namedarray import utils

# Use as a sentinel value to indicate a dtype appropriate NA value.
NA = utils.ReprObject("<NA>")


@functools.total_ordering
class AlwaysGreaterThan:
    def __gt__(self, other: object) -> Literal[True]:
        return True

    def __eq__(self, other: object) -> bool:
        return isinstance(other, type(self))


@functools.total_ordering
class AlwaysLessThan:
    def __lt__(self, other: object) -> Literal[True]:
        return True

    def __eq__(self, other: object) -> bool:
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


def maybe_promote(dtype: np.dtype[np.generic]) -> tuple[np.dtype[np.generic], Any]:
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
    if np.issubdtype(dtype, np.floating):
        dtype_ = dtype
        fill_value = np.nan
    elif np.issubdtype(dtype, np.timedelta64):
        # See https://github.com/numpy/numpy/issues/10685
        # np.timedelta64 is a subclass of np.integer
        # Check np.timedelta64 before np.integer
        fill_value = np.timedelta64("NaT")
        dtype_ = dtype
    elif np.issubdtype(dtype, np.integer):
        dtype_ = np.float32 if dtype.itemsize <= 2 else np.float64
        fill_value = np.nan
    elif np.issubdtype(dtype, np.complexfloating):
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


def get_fill_value(dtype: np.dtype[np.generic]) -> Any:
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


def get_pos_infinity(
    dtype: np.dtype[np.generic], max_for_int: bool = False
) -> float | complex | AlwaysGreaterThan:
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
    if issubclass(dtype.type, np.floating):
        return np.inf

    if issubclass(dtype.type, np.integer):
        return np.iinfo(dtype.type).max if max_for_int else np.inf
    if issubclass(dtype.type, np.complexfloating):
        return np.inf + 1j * np.inf

    return INF


def get_neg_infinity(
    dtype: np.dtype[np.generic], min_for_int: bool = False
) -> float | complex | AlwaysLessThan:
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
    if issubclass(dtype.type, np.floating):
        return -np.inf

    if issubclass(dtype.type, np.integer):
        return np.iinfo(dtype.type).min if min_for_int else -np.inf
    if issubclass(dtype.type, np.complexfloating):
        return -np.inf - 1j * np.inf

    return NINF


def is_datetime_like(
    dtype: np.dtype[np.generic],
) -> TypeGuard[np.datetime64 | np.timedelta64]:
    """Check if a dtype is a subclass of the numpy datetime types"""
    return np.issubdtype(dtype, np.datetime64) or np.issubdtype(dtype, np.timedelta64)


def result_type(
    *arrays_and_dtypes: np.typing.ArrayLike | np.typing.DTypeLike,
) -> np.dtype[np.generic]:
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
    types = {np.result_type(t).type for t in arrays_and_dtypes}

    for left, right in PROMOTE_TO_OBJECT:
        if any(issubclass(t, left) for t in types) and any(
            issubclass(t, right) for t in types
        ):
            return np.dtype(object)

    return np.result_type(*arrays_and_dtypes)
