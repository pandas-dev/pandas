"""
Functions for implementing 'astype' methods according to pandas conventions,
particularly ones that differ from numpy.
"""
from __future__ import annotations

import inspect
from typing import (
    TYPE_CHECKING,
    overload,
)

import numpy as np

from pandas._libs import lib
from pandas._libs.tslibs.timedeltas import array_to_timedelta64
from pandas._typing import (
    ArrayLike,
    DtypeObj,
    IgnoreRaise,
)
from pandas.errors import IntCastingNaNError

from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_dtype_equal,
    is_integer_dtype,
    is_object_dtype,
    is_timedelta64_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import (
    ExtensionDtype,
    PandasDtype,
)
from pandas.core.dtypes.missing import isna

if TYPE_CHECKING:
    from pandas.core.arrays import ExtensionArray


_dtype_obj = np.dtype(object)


@overload
def astype_nansafe(
    arr: np.ndarray, dtype: np.dtype, copy: bool = ..., skipna: bool = ...
) -> np.ndarray:
    ...


@overload
def astype_nansafe(
    arr: np.ndarray, dtype: ExtensionDtype, copy: bool = ..., skipna: bool = ...
) -> ExtensionArray:
    ...


def astype_nansafe(
    arr: np.ndarray, dtype: DtypeObj, copy: bool = True, skipna: bool = False
) -> ArrayLike:
    """
    Cast the elements of an array to a given dtype a nan-safe manner.

    Parameters
    ----------
    arr : ndarray
    dtype : np.dtype or ExtensionDtype
    copy : bool, default True
        If False, a view will be attempted but may fail, if
        e.g. the item sizes don't align.
    skipna: bool, default False
        Whether or not we should skip NaN when casting as a string-type.

    Raises
    ------
    ValueError
        The dtype was a datetime64/timedelta64 dtype, but it had no unit.
    """

    # dispatch on extension dtype if needed
    if isinstance(dtype, ExtensionDtype):
        return dtype.construct_array_type()._from_sequence(arr, dtype=dtype, copy=copy)

    elif not isinstance(dtype, np.dtype):  # pragma: no cover
        raise ValueError("dtype must be np.dtype or ExtensionDtype")

    if arr.dtype.kind in ["m", "M"] and (
        issubclass(dtype.type, str) or dtype == _dtype_obj
    ):
        from pandas.core.construction import ensure_wrapped_if_datetimelike

        arr = ensure_wrapped_if_datetimelike(arr)
        return arr.astype(dtype, copy=copy)

    if issubclass(dtype.type, str):
        shape = arr.shape
        if arr.ndim > 1:
            arr = arr.ravel()
        return lib.ensure_string_array(
            arr, skipna=skipna, convert_na_value=False
        ).reshape(shape)

    elif is_datetime64_dtype(arr.dtype):
        if dtype == np.int64:
            if isna(arr).any():
                raise ValueError("Cannot convert NaT values to integer")
            return arr.view(dtype)

        # allow frequency conversions
        if dtype.kind == "M":
            return arr.astype(dtype)

        raise TypeError(f"cannot astype a datetimelike from [{arr.dtype}] to [{dtype}]")

    elif is_timedelta64_dtype(arr.dtype):
        if dtype == np.int64:
            if isna(arr).any():
                raise ValueError("Cannot convert NaT values to integer")
            return arr.view(dtype)

        elif dtype.kind == "m":
            # give the requested dtype for supported units (s, ms, us, ns)
            #  and doing the old convert-to-float behavior otherwise.
            from pandas.core.construction import ensure_wrapped_if_datetimelike

            arr = ensure_wrapped_if_datetimelike(arr)
            return arr.astype(dtype, copy=copy)

        raise TypeError(f"cannot astype a timedelta from [{arr.dtype}] to [{dtype}]")

    elif np.issubdtype(arr.dtype, np.floating) and is_integer_dtype(dtype):
        return _astype_float_to_int_nansafe(arr, dtype, copy)

    elif is_object_dtype(arr.dtype):

        # if we have a datetime/timedelta array of objects
        # then coerce to datetime64[ns] and use DatetimeArray.astype

        if is_datetime64_dtype(dtype):
            from pandas import to_datetime

            dti = to_datetime(arr.ravel())
            dta = dti._data.reshape(arr.shape)
            return dta.astype(dtype, copy=False)._ndarray

        elif is_timedelta64_dtype(dtype):
            # bc we know arr.dtype == object, this is equivalent to
            #  `np.asarray(to_timedelta(arr))`, but using a lower-level API that
            #  does not require a circular import.
            return array_to_timedelta64(arr).view("m8[ns]").astype(dtype, copy=False)

    if dtype.name in ("datetime64", "timedelta64"):
        msg = (
            f"The '{dtype.name}' dtype has no unit. Please pass in "
            f"'{dtype.name}[ns]' instead."
        )
        raise ValueError(msg)

    if copy or is_object_dtype(arr.dtype) or is_object_dtype(dtype):
        # Explicit copy, or required since NumPy can't view from / to object.
        return arr.astype(dtype, copy=True)

    return arr.astype(dtype, copy=copy)


def _astype_float_to_int_nansafe(
    values: np.ndarray, dtype: np.dtype, copy: bool
) -> np.ndarray:
    """
    astype with a check preventing converting NaN to an meaningless integer value.
    """
    if not np.isfinite(values).all():
        raise IntCastingNaNError(
            "Cannot convert non-finite values (NA or inf) to integer"
        )
    if dtype.kind == "u":
        # GH#45151
        if not (values >= 0).all():
            raise ValueError(f"Cannot losslessly cast from {values.dtype} to {dtype}")
    return values.astype(dtype, copy=copy)


def astype_array(values: ArrayLike, dtype: DtypeObj, copy: bool = False) -> ArrayLike:
    """
    Cast array (ndarray or ExtensionArray) to the new dtype.

    Parameters
    ----------
    values : ndarray or ExtensionArray
    dtype : dtype object
    copy : bool, default False
        copy if indicated

    Returns
    -------
    ndarray or ExtensionArray
    """
    if (
        values.dtype.kind in ["m", "M"]
        and dtype.kind in ["i", "u"]
        and isinstance(dtype, np.dtype)
        and dtype.itemsize != 8
    ):
        # TODO(2.0) remove special case once deprecation on DTA/TDA is enforced
        msg = rf"cannot astype a datetimelike from [{values.dtype}] to [{dtype}]"
        raise TypeError(msg)

    if is_dtype_equal(values.dtype, dtype):
        if copy:
            return values.copy()
        return values

    if not isinstance(values, np.ndarray):
        # i.e. ExtensionArray
        values = values.astype(dtype, copy=copy)

    else:
        values = astype_nansafe(values, dtype, copy=copy)

    # in pandas we don't store numpy str dtypes, so convert to object
    if isinstance(dtype, np.dtype) and issubclass(values.dtype.type, str):
        values = np.array(values, dtype=object)

    return values


def astype_array_safe(
    values: ArrayLike, dtype, copy: bool = False, errors: IgnoreRaise = "raise"
) -> ArrayLike:
    """
    Cast array (ndarray or ExtensionArray) to the new dtype.

    This basically is the implementation for DataFrame/Series.astype and
    includes all custom logic for pandas (NaN-safety, converting str to object,
    not allowing )

    Parameters
    ----------
    values : ndarray or ExtensionArray
    dtype : str, dtype convertible
    copy : bool, default False
        copy if indicated
    errors : str, {'raise', 'ignore'}, default 'raise'
        - ``raise`` : allow exceptions to be raised
        - ``ignore`` : suppress exceptions. On error return original object

    Returns
    -------
    ndarray or ExtensionArray
    """
    errors_legal_values = ("raise", "ignore")

    if errors not in errors_legal_values:
        invalid_arg = (
            "Expected value of kwarg 'errors' to be one of "
            f"{list(errors_legal_values)}. Supplied value is '{errors}'"
        )
        raise ValueError(invalid_arg)

    if inspect.isclass(dtype) and issubclass(dtype, ExtensionDtype):
        msg = (
            f"Expected an instance of {dtype.__name__}, "
            "but got the class instead. Try instantiating 'dtype'."
        )
        raise TypeError(msg)

    dtype = pandas_dtype(dtype)
    if isinstance(dtype, PandasDtype):
        # Ensure we don't end up with a PandasArray
        dtype = dtype.numpy_dtype

    try:
        new_values = astype_array(values, dtype, copy=copy)
    except (ValueError, TypeError):
        # e.g. astype_nansafe can fail on object-dtype of strings
        #  trying to convert to float
        if errors == "ignore":
            new_values = values
        else:
            raise

    return new_values
