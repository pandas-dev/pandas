"""Compatibility module defining operations on duck numpy-arrays.

Currently, this means Dask or NumPy arrays. None of these functions should
accept or return xarray objects.
"""

from __future__ import annotations

import contextlib
import datetime
import inspect
import warnings
from collections.abc import Callable
from functools import partial
from importlib import import_module
from typing import Any

import numpy as np
import pandas as pd
from numpy import (
    isclose,
    isnat,
    take,
    unravel_index,  # noqa: F401
)

from xarray.compat import dask_array_compat, dask_array_ops
from xarray.compat.array_api_compat import get_array_namespace
from xarray.core import dtypes, nputils
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.options import OPTIONS
from xarray.core.utils import is_duck_array, is_duck_dask_array, module_available
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import array_type, is_chunked_array

# remove once numpy 2.0 is the oldest supported version
if module_available("numpy", minversion="2.0.0.dev0"):
    from numpy.lib.array_utils import (  # type: ignore[import-not-found,unused-ignore]
        normalize_axis_index,
    )
else:
    from numpy.core.multiarray import (  # type: ignore[attr-defined,no-redef,unused-ignore]
        normalize_axis_index,
    )


dask_available = module_available("dask")


def einsum(*args, **kwargs):
    if OPTIONS["use_opt_einsum"] and module_available("opt_einsum"):
        import opt_einsum

        return opt_einsum.contract(*args, **kwargs)
    else:
        xp = get_array_namespace(*args)
        return xp.einsum(*args, **kwargs)


def tensordot(*args, **kwargs):
    xp = get_array_namespace(*args)
    return xp.tensordot(*args, **kwargs)


def cross(*args, **kwargs):
    xp = get_array_namespace(*args)
    return xp.cross(*args, **kwargs)


def gradient(f, *varargs, axis=None, edge_order=1):
    xp = get_array_namespace(f)
    return xp.gradient(f, *varargs, axis=axis, edge_order=edge_order)


def _dask_or_eager_func(
    name,
    eager_module=np,
    dask_module="dask.array",
    dask_only_kwargs=tuple(),
    numpy_only_kwargs=tuple(),
):
    """Create a function that dispatches to dask for dask array inputs."""

    def f(*args, **kwargs):
        if dask_available and any(is_duck_dask_array(a) for a in args):
            mod = (
                import_module(dask_module)
                if isinstance(dask_module, str)
                else dask_module
            )
            wrapped = getattr(mod, name)
            for kwarg in numpy_only_kwargs:
                kwargs.pop(kwarg, None)
        else:
            wrapped = getattr(eager_module, name)
            for kwarg in dask_only_kwargs:
                kwargs.pop(kwarg, None)
        return wrapped(*args, **kwargs)

    return f


def fail_on_dask_array_input(values, msg=None, func_name=None):
    if is_duck_dask_array(values):
        if msg is None:
            msg = "%r is not yet a valid method on dask arrays"
        if func_name is None:
            func_name = inspect.stack()[1][3]
        raise NotImplementedError(msg % func_name)


# Requires special-casing because pandas won't automatically dispatch to dask.isnull via NEP-18
pandas_isnull = _dask_or_eager_func("isnull", eager_module=pd, dask_module="dask.array")

# TODO replace with simply np.ma.masked_invalid once numpy/numpy#16022 is fixed
# TODO: replacing breaks iris + dask tests
masked_invalid = _dask_or_eager_func(
    "masked_invalid", eager_module=np.ma, dask_module="dask.array.ma"
)


def sliding_window_view(array, window_shape, axis=None, **kwargs):
    # TODO: some libraries (e.g. jax) don't have this, implement an alternative?
    xp = get_array_namespace(array)
    # sliding_window_view will not dispatch arbitrary kwargs (automatic_rechunk),
    # so we need to hand-code this.
    func = _dask_or_eager_func(
        "sliding_window_view",
        eager_module=xp.lib.stride_tricks,
        dask_module=dask_array_compat,
        dask_only_kwargs=("automatic_rechunk",),
        numpy_only_kwargs=("subok", "writeable"),
    )
    return func(array, window_shape, axis=axis, **kwargs)


def round(array):
    xp = get_array_namespace(array)
    return xp.round(array)


around: Callable = round


def isna(data: Any) -> bool:
    """Checks if data is literally np.nan or pd.NA.

    Parameters
    ----------
    data
        Any python object

    Returns
    -------
        Whether or not the data is np.nan or pd.NA
    """
    return data is pd.NA or data is np.nan  # noqa: PLW0177


def isnull(data):
    data = asarray(data)

    xp = get_array_namespace(data)
    scalar_type = data.dtype
    if dtypes.is_datetime_like(scalar_type):
        # datetime types use NaT for null
        # note: must check timedelta64 before integers, because currently
        # timedelta64 inherits from np.integer
        return isnat(data)
    elif dtypes.isdtype(scalar_type, ("real floating", "complex floating"), xp=xp):
        # float types use NaN for null
        xp = get_array_namespace(data)
        return xp.isnan(data)
    elif dtypes.isdtype(scalar_type, ("bool", "integral"), xp=xp) or (
        isinstance(scalar_type, np.dtype)
        and (
            np.issubdtype(scalar_type, np.character)
            or np.issubdtype(scalar_type, np.void)
        )
    ):
        # these types cannot represent missing values
        # bool_ is for backwards compat with numpy<2, and cupy
        dtype = xp.bool_ if hasattr(xp, "bool_") else xp.bool
        return full_like(data, dtype=dtype, fill_value=False)
    # at this point, array should have dtype=object
    elif isinstance(data, np.ndarray) or pd.api.types.is_extension_array_dtype(data):  # noqa: TID251
        return pandas_isnull(data)
    else:
        # Not reachable yet, but intended for use with other duck array
        # types. For full consistency with pandas, we should accept None as
        # a null value as well as NaN, but it isn't clear how to do this
        # with duck typing.
        return data != data  # noqa: PLR0124


def notnull(data):
    return ~isnull(data)


def trapz(y, x, axis):
    if axis < 0:
        axis = y.ndim + axis
    x_sl1 = (slice(1, None),) + (None,) * (y.ndim - axis - 1)
    x_sl2 = (slice(None, -1),) + (None,) * (y.ndim - axis - 1)
    slice1 = (slice(None),) * axis + (slice(1, None),)
    slice2 = (slice(None),) * axis + (slice(None, -1),)
    dx = x[x_sl1] - x[x_sl2]
    integrand = dx * 0.5 * (y[tuple(slice1)] + y[tuple(slice2)])
    return sum(integrand, axis=axis, skipna=False)


def cumulative_trapezoid(y, x, axis):
    if axis < 0:
        axis = y.ndim + axis
    x_sl1 = (slice(1, None),) + (None,) * (y.ndim - axis - 1)
    x_sl2 = (slice(None, -1),) + (None,) * (y.ndim - axis - 1)
    slice1 = (slice(None),) * axis + (slice(1, None),)
    slice2 = (slice(None),) * axis + (slice(None, -1),)
    dx = x[x_sl1] - x[x_sl2]
    integrand = dx * 0.5 * (y[tuple(slice1)] + y[tuple(slice2)])

    # Pad so that 'axis' has same length in result as it did in y
    pads = [(1, 0) if i == axis else (0, 0) for i in range(y.ndim)]

    xp = get_array_namespace(y, x)
    integrand = xp.pad(integrand, pads, mode="constant", constant_values=0.0)

    return cumsum(integrand, axis=axis, skipna=False)


def full_like(a, fill_value, **kwargs):
    xp = get_array_namespace(a)
    return xp.full_like(a, fill_value, **kwargs)


def empty_like(a, **kwargs):
    xp = get_array_namespace(a)
    return xp.empty_like(a, **kwargs)


def astype(data, dtype, *, xp=None, **kwargs):
    if not hasattr(data, "__array_namespace__") and xp is None:
        return data.astype(dtype, **kwargs)

    if xp is None:
        xp = get_array_namespace(data)

    if xp == np:
        # numpy currently doesn't have a astype:
        return data.astype(dtype, **kwargs)
    return xp.astype(data, dtype, **kwargs)


def asarray(data, xp=np, dtype=None):
    converted = data if is_duck_array(data) else xp.asarray(data)

    if dtype is None or converted.dtype == dtype:
        return converted

    if xp is np or not hasattr(xp, "astype"):
        return converted.astype(dtype)
    else:
        return xp.astype(converted, dtype)


def as_shared_dtype(scalars_or_arrays, xp=None):
    """Cast arrays to a shared dtype using xarray's type promotion rules."""
    extension_array_types = [
        x.dtype
        for x in scalars_or_arrays
        if pd.api.types.is_extension_array_dtype(x)  # noqa: TID251
    ]
    if len(extension_array_types) >= 1:
        non_nans = [x for x in scalars_or_arrays if not isna(x)]
        if len(extension_array_types) == len(non_nans) and all(
            isinstance(x, type(extension_array_types[0])) for x in extension_array_types
        ):
            return [
                x
                if not isna(x)
                else PandasExtensionArray(
                    type(non_nans[0].array)._from_sequence([x], dtype=non_nans[0].dtype)
                )
                for x in scalars_or_arrays
            ]
        raise ValueError(
            f"Cannot cast values to shared type, found values: {scalars_or_arrays}"
        )

    # Avoid calling array_type("cupy") repeatidely in the any check
    array_type_cupy = array_type("cupy")
    if any(isinstance(x, array_type_cupy) for x in scalars_or_arrays):
        import cupy as cp

        xp = cp
    elif xp is None:
        xp = get_array_namespace(scalars_or_arrays)

    # Pass arrays directly instead of dtypes to result_type so scalars
    # get handled properly.
    # Note that result_type() safely gets the dtype from dask arrays without
    # evaluating them.
    dtype = dtypes.result_type(*scalars_or_arrays, xp=xp)

    return [asarray(x, dtype=dtype, xp=xp) for x in scalars_or_arrays]


def broadcast_to(array, shape):
    xp = get_array_namespace(array)
    return xp.broadcast_to(array, shape)


def lazy_array_equiv(arr1, arr2):
    """Like array_equal, but doesn't actually compare values.
    Returns True when arr1, arr2 identical or their dask tokens are equal.
    Returns False when shapes are not equal.
    Returns None when equality cannot determined: one or both of arr1, arr2 are numpy arrays;
    or their dask tokens are not equal
    """
    if arr1 is arr2:
        return True
    arr1 = asarray(arr1)
    arr2 = asarray(arr2)
    if arr1.shape != arr2.shape:
        return False
    if dask_available and is_duck_dask_array(arr1) and is_duck_dask_array(arr2):
        from dask.base import tokenize

        # GH3068, GH4221
        if tokenize(arr1) == tokenize(arr2):
            return True
        else:
            return None
    return None


def allclose_or_equiv(arr1, arr2, rtol=1e-5, atol=1e-8):
    """Like np.allclose, but also allows values to be NaN in both arrays"""
    arr1 = asarray(arr1)
    arr2 = asarray(arr2)

    lazy_equiv = lazy_array_equiv(arr1, arr2)
    if lazy_equiv is None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
            return bool(
                array_all(isclose(arr1, arr2, rtol=rtol, atol=atol, equal_nan=True))
            )
    else:
        return lazy_equiv


def array_equiv(arr1, arr2):
    """Like np.array_equal, but also allows values to be NaN in both arrays"""
    arr1 = asarray(arr1)
    arr2 = asarray(arr2)
    lazy_equiv = lazy_array_equiv(arr1, arr2)
    if lazy_equiv is None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "In the future, 'NAT == x'")
            flag_array = (arr1 == arr2) | (isnull(arr1) & isnull(arr2))
            return bool(array_all(flag_array))
    else:
        return lazy_equiv


def array_notnull_equiv(arr1, arr2):
    """Like np.array_equal, but also allows values to be NaN in either or both
    arrays
    """
    arr1 = asarray(arr1)
    arr2 = asarray(arr2)
    lazy_equiv = lazy_array_equiv(arr1, arr2)
    if lazy_equiv is None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "In the future, 'NAT == x'")
            flag_array = (arr1 == arr2) | isnull(arr1) | isnull(arr2)
            return bool(array_all(flag_array))
    else:
        return lazy_equiv


def count(data, axis=None):
    """Count the number of non-NA in this array along the given axis or axes"""
    xp = get_array_namespace(data)
    return xp.sum(xp.logical_not(isnull(data)), axis=axis)


def sum_where(data, axis=None, dtype=None, where=None):
    xp = get_array_namespace(data)
    if where is not None:
        a = where_method(xp.zeros_like(data), where, data)
    else:
        a = data
    result = xp.sum(a, axis=axis, dtype=dtype)
    return result


def where(condition, x, y):
    """Three argument where() with better dtype promotion rules."""
    xp = get_array_namespace(condition, x, y)

    dtype = xp.bool_ if hasattr(xp, "bool_") else xp.bool
    if not is_duck_array(condition):
        condition = asarray(condition, dtype=dtype, xp=xp)
    else:
        condition = astype(condition, dtype=dtype, xp=xp)

    return xp.where(condition, *as_shared_dtype([x, y], xp=xp))


def where_method(data, cond, other=dtypes.NA):
    if other is dtypes.NA:
        other = dtypes.get_fill_value(data.dtype)
    return where(cond, data, other)


def fillna(data, other):
    # we need to pass data first so pint has a chance of returning the
    # correct unit
    # TODO: revert after https://github.com/hgrecco/pint/issues/1019 is fixed
    return where(notnull(data), data, other)


def logical_not(data):
    xp = get_array_namespace(data)
    return xp.logical_not(data)


def clip(data, min=None, max=None):
    xp = get_array_namespace(data)
    return xp.clip(data, min, max)


def concatenate(arrays, axis=0):
    """concatenate() with better dtype promotion rules."""
    # TODO: `concat` is the xp compliant name, but fallback to concatenate for
    # older numpy and for cupy
    xp = get_array_namespace(*arrays)
    if hasattr(xp, "concat"):
        return xp.concat(as_shared_dtype(arrays, xp=xp), axis=axis)
    else:
        return xp.concatenate(as_shared_dtype(arrays, xp=xp), axis=axis)


def stack(arrays, axis=0):
    """stack() with better dtype promotion rules."""
    xp = get_array_namespace(arrays[0])
    return xp.stack(as_shared_dtype(arrays, xp=xp), axis=axis)


def reshape(array, shape):
    xp = get_array_namespace(array)
    return xp.reshape(array, shape)


def ravel(array):
    return reshape(array, (-1,))


def transpose(array, axes=None):
    xp = get_array_namespace(array)
    return xp.transpose(array, axes)


def moveaxis(array, source, destination):
    xp = get_array_namespace(array)
    return xp.moveaxis(array, source, destination)


def pad(array, pad_width, **kwargs):
    xp = get_array_namespace(array)
    return xp.pad(array, pad_width, **kwargs)


def quantile(array, q, axis=None, **kwargs):
    xp = get_array_namespace(array)
    return xp.quantile(array, q, axis=axis, **kwargs)


@contextlib.contextmanager
def _ignore_warnings_if(condition):
    if condition:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            yield
    else:
        yield


def _create_nan_agg_method(name, coerce_strings=False, invariant_0d=False):
    def f(values, axis=None, skipna=None, **kwargs):
        if kwargs.pop("out", None) is not None:
            raise TypeError(f"`out` is not valid for {name}")

        # The data is invariant in the case of 0d data, so do not
        # change the data (and dtype)
        # See https://github.com/pydata/xarray/issues/4885
        if invariant_0d and axis == ():
            return values

        xp = get_array_namespace(values)
        values = asarray(values, xp=xp)

        if coerce_strings and dtypes.is_string(values.dtype):
            values = astype(values, object)

        func = None
        if skipna or (
            skipna is None
            and (
                dtypes.isdtype(
                    values.dtype, ("complex floating", "real floating"), xp=xp
                )
                or dtypes.is_object(values.dtype)
            )
        ):
            from xarray.computation import nanops

            nanname = "nan" + name
            func = getattr(nanops, nanname)
        else:
            if name in ["sum", "prod"]:
                kwargs.pop("min_count", None)

            xp = get_array_namespace(values)
            func = getattr(xp, name)

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "All-NaN slice encountered")
                return func(values, axis=axis, **kwargs)
        except AttributeError:
            if not is_duck_dask_array(values):
                raise
            try:  # dask/dask#3133 dask sometimes needs dtype argument
                # if func does not accept dtype, then raises TypeError
                return func(values, axis=axis, dtype=values.dtype, **kwargs)
            except (AttributeError, TypeError) as err:
                raise NotImplementedError(
                    f"{name} is not yet implemented on dask arrays"
                ) from err

    f.__name__ = name
    return f


# Attributes `numeric_only`, `available_min_count` is used for docs.
# See ops.inject_reduce_methods
argmax = _create_nan_agg_method("argmax", coerce_strings=True)
argmin = _create_nan_agg_method("argmin", coerce_strings=True)
max = _create_nan_agg_method("max", coerce_strings=True, invariant_0d=True)
min = _create_nan_agg_method("min", coerce_strings=True, invariant_0d=True)
sum = _create_nan_agg_method("sum", invariant_0d=True)
sum.numeric_only = True
sum.available_min_count = True
std = _create_nan_agg_method("std")
std.numeric_only = True
var = _create_nan_agg_method("var")
var.numeric_only = True
median = _create_nan_agg_method("median", invariant_0d=True)
median.numeric_only = True
prod = _create_nan_agg_method("prod", invariant_0d=True)
prod.numeric_only = True
prod.available_min_count = True
cumprod_1d = _create_nan_agg_method("cumprod", invariant_0d=True)
cumprod_1d.numeric_only = True
cumsum_1d = _create_nan_agg_method("cumsum", invariant_0d=True)
cumsum_1d.numeric_only = True


def array_all(array, axis=None, keepdims=False, **kwargs):
    xp = get_array_namespace(array)
    return xp.all(array, axis=axis, keepdims=keepdims, **kwargs)


def array_any(array, axis=None, keepdims=False, **kwargs):
    xp = get_array_namespace(array)
    return xp.any(array, axis=axis, keepdims=keepdims, **kwargs)


_mean = _create_nan_agg_method("mean", invariant_0d=True)


def _datetime_nanmin(array):
    return _datetime_nanreduce(array, min)


def _datetime_nanreduce(array, func):
    """nanreduce() function for datetime64.

    Caveats that this function deals with:

    - In numpy < 1.18, min() on datetime64 incorrectly ignores NaT
    - numpy nanmin() don't work on datetime64 (all versions at the moment of writing)
    - dask min() does not work on datetime64 (all versions at the moment of writing)
    """
    dtype = array.dtype
    assert dtypes.is_datetime_like(dtype)
    # (NaT).astype(float) does not produce NaN...
    array = where(pandas_isnull(array), np.nan, array.astype(float))
    array = func(array, skipna=True)
    if isinstance(array, float):
        array = np.array(array)
    # ...but (NaN).astype("M8") does produce NaT
    return array.astype(dtype)


def datetime_to_numeric(array, offset=None, datetime_unit=None, dtype=float):
    """Convert an array containing datetime-like data to numerical values.
    Convert the datetime array to a timedelta relative to an offset.
    Parameters
    ----------
    array : array-like
        Input data
    offset : None, datetime or cftime.datetime
        Datetime offset. If None, this is set by default to the array's minimum
        value to reduce round off errors.
    datetime_unit : {None, Y, M, W, D, h, m, s, ms, us, ns, ps, fs, as}
        If not None, convert output to a given datetime unit. Note that some
        conversions are not allowed due to non-linear relationships between units.
    dtype : dtype
        Output dtype.
    Returns
    -------
    array
        Numerical representation of datetime object relative to an offset.
    Notes
    -----
    Some datetime unit conversions won't work, for example from days to years, even
    though some calendars would allow for them (e.g. no_leap). This is because there
    is no `cftime.timedelta` object.
    """
    # Set offset to minimum if not given
    if offset is None:
        if dtypes.is_datetime_like(array.dtype):
            offset = _datetime_nanreduce(array, min)
        else:
            offset = min(array)

    # Compute timedelta object.
    # For np.datetime64, this can silently yield garbage due to overflow.
    # One option is to enforce 1970-01-01 as the universal offset.

    # This map_blocks call is for backwards compatibility.
    # dask == 2021.04.1 does not support subtracting object arrays
    # which is required for cftime
    if is_duck_dask_array(array) and dtypes.is_object(array.dtype):
        array = array.map_blocks(lambda a, b: a - b, offset, meta=array._meta)
    else:
        array = array - offset

    # Scalar is converted to 0d-array
    if not hasattr(array, "dtype"):
        array = np.array(array)

    # Convert timedelta objects to float by first converting to microseconds.
    if dtypes.is_object(array.dtype):
        return py_timedelta_to_float(array, datetime_unit or "ns").astype(dtype)

    # Convert np.NaT to np.nan
    elif dtypes.is_datetime_like(array.dtype):
        # Convert to specified timedelta units.
        if datetime_unit:
            array = array / np.timedelta64(1, datetime_unit)
        return np.where(isnull(array), np.nan, array.astype(dtype))


def timedelta_to_numeric(value, datetime_unit="ns", dtype=float):
    """Convert a timedelta-like object to numerical values.

    Parameters
    ----------
    value : datetime.timedelta, numpy.timedelta64, pandas.Timedelta, str
        Time delta representation.
    datetime_unit : {Y, M, W, D, h, m, s, ms, us, ns, ps, fs, as}
        The time units of the output values. Note that some conversions are not allowed due to
        non-linear relationships between units.
    dtype : type
        The output data type.

    """
    if isinstance(value, datetime.timedelta):
        out = py_timedelta_to_float(value, datetime_unit)
    elif isinstance(value, np.timedelta64):
        out = np_timedelta64_to_float(value, datetime_unit)
    elif isinstance(value, pd.Timedelta):
        out = pd_timedelta_to_float(value, datetime_unit)
    elif isinstance(value, str):
        try:
            a = pd.to_timedelta(value)
        except ValueError as err:
            raise ValueError(
                f"Could not convert {value!r} to timedelta64 using pandas.to_timedelta"
            ) from err
        return py_timedelta_to_float(a, datetime_unit)
    else:
        raise TypeError(
            f"Expected value of type str, pandas.Timedelta, datetime.timedelta "
            f"or numpy.timedelta64, but received {type(value).__name__}"
        )
    return out.astype(dtype)


def _to_pytimedelta(array, unit="us"):
    return array.astype(f"timedelta64[{unit}]").astype(datetime.timedelta)


def np_timedelta64_to_float(array, datetime_unit):
    """Convert numpy.timedelta64 to float, possibly at a loss of resolution."""
    unit, _ = np.datetime_data(array.dtype)
    conversion_factor = np.timedelta64(1, unit) / np.timedelta64(1, datetime_unit)
    return conversion_factor * array.astype(np.float64)


def pd_timedelta_to_float(value, datetime_unit):
    """Convert pandas.Timedelta to float.

    Notes
    -----
    Built on the assumption that pandas timedelta values are in nanoseconds,
    which is also the numpy default resolution.
    """
    value = value.to_timedelta64()
    return np_timedelta64_to_float(value, datetime_unit)


def _timedelta_to_seconds(array):
    if isinstance(array, datetime.timedelta):
        return array.total_seconds() * 1e6
    else:
        return np.reshape([a.total_seconds() for a in array.ravel()], array.shape) * 1e6


def py_timedelta_to_float(array, datetime_unit):
    """Convert a timedelta object to a float, possibly at a loss of resolution."""
    array = asarray(array)
    if is_duck_dask_array(array):
        array = array.map_blocks(
            _timedelta_to_seconds, meta=np.array([], dtype=np.float64)
        )
    else:
        array = _timedelta_to_seconds(array)
    conversion_factor = np.timedelta64(1, "us") / np.timedelta64(1, datetime_unit)
    return conversion_factor * array


def mean(array, axis=None, skipna=None, **kwargs):
    """inhouse mean that can handle np.datetime64 or cftime.datetime
    dtypes"""
    from xarray.core.common import _contains_cftime_datetimes

    array = asarray(array)
    if dtypes.is_datetime_like(array.dtype):
        dmin = _datetime_nanreduce(array, min).astype("datetime64[Y]").astype(int)
        dmax = _datetime_nanreduce(array, max).astype("datetime64[Y]").astype(int)
        offset = (
            np.array((dmin + dmax) // 2).astype("datetime64[Y]").astype(array.dtype)
        )
        # From version 2025.01.2 xarray uses np.datetime64[unit], where unit
        # is one of "s", "ms", "us", "ns".
        # To not have to worry about the resolution, we just convert the output
        # to "timedelta64" (without unit) and let the dtype of offset take precedence.
        # This is fully backwards compatible with datetime64[ns].
        return (
            _mean(
                datetime_to_numeric(array, offset), axis=axis, skipna=skipna, **kwargs
            ).astype("timedelta64")
            + offset
        )
    elif _contains_cftime_datetimes(array):
        offset = min(array)
        timedeltas = datetime_to_numeric(array, offset, datetime_unit="us")
        mean_timedeltas = _mean(timedeltas, axis=axis, skipna=skipna, **kwargs)
        return _to_pytimedelta(mean_timedeltas, unit="us") + offset
    else:
        return _mean(array, axis=axis, skipna=skipna, **kwargs)


mean.numeric_only = True  # type: ignore[attr-defined]


def _nd_cum_func(cum_func, array, axis, **kwargs):
    array = asarray(array)
    if axis is None:
        axis = tuple(range(array.ndim))
    if isinstance(axis, int):
        axis = (axis,)

    out = array
    for ax in axis:
        out = cum_func(out, axis=ax, **kwargs)
    return out


def ndim(array) -> int:
    # Required part of the duck array and the array-api, but we fall back in case
    # https://docs.xarray.dev/en/latest/internals/duck-arrays-integration.html#duck-array-requirements
    return array.ndim if hasattr(array, "ndim") else np.ndim(array)


def cumprod(array, axis=None, **kwargs):
    """N-dimensional version of cumprod."""
    return _nd_cum_func(cumprod_1d, array, axis, **kwargs)


def cumsum(array, axis=None, **kwargs):
    """N-dimensional version of cumsum."""
    return _nd_cum_func(cumsum_1d, array, axis, **kwargs)


def first(values, axis, skipna=None):
    """Return the first non-NA elements in this array along the given axis"""
    if (skipna or skipna is None) and not (
        dtypes.isdtype(values.dtype, "signed integer") or dtypes.is_string(values.dtype)
    ):
        # only bother for dtypes that can hold NaN
        if is_chunked_array(values):
            return chunked_nanfirst(values, axis)
        else:
            return nputils.nanfirst(values, axis)

    return take(values, 0, axis=axis)


def last(values, axis, skipna=None):
    """Return the last non-NA elements in this array along the given axis"""
    if (skipna or skipna is None) and not (
        dtypes.isdtype(values.dtype, "signed integer") or dtypes.is_string(values.dtype)
    ):
        # only bother for dtypes that can hold NaN
        if is_chunked_array(values):
            return chunked_nanlast(values, axis)
        else:
            return nputils.nanlast(values, axis)

    return take(values, -1, axis=axis)


def isin(element, test_elements, **kwargs):
    xp = get_array_namespace(element, test_elements)
    return xp.isin(element, test_elements, **kwargs)


def least_squares(lhs, rhs, rcond=None, skipna=False):
    """Return the coefficients and residuals of a least-squares fit."""
    if is_duck_dask_array(rhs):
        return dask_array_ops.least_squares(lhs, rhs, rcond=rcond, skipna=skipna)
    else:
        return nputils.least_squares(lhs, rhs, rcond=rcond, skipna=skipna)


def _push(array, n: int | None = None, axis: int = -1):
    """
    Use either bottleneck or numbagg depending on options & what's available
    """

    if not OPTIONS["use_bottleneck"] and not OPTIONS["use_numbagg"]:
        raise RuntimeError(
            "ffill & bfill requires bottleneck or numbagg to be enabled."
            " Call `xr.set_options(use_bottleneck=True)` or `xr.set_options(use_numbagg=True)` to enable one."
        )
    if OPTIONS["use_numbagg"] and module_available("numbagg"):
        import numbagg

        return numbagg.ffill(array, limit=n, axis=axis)

    # work around for bottleneck 178
    limit = n if n is not None else array.shape[axis]

    import bottleneck as bn

    return bn.push(array, limit, axis)


def push(array, n, axis, method="blelloch"):
    if not OPTIONS["use_bottleneck"] and not OPTIONS["use_numbagg"]:
        raise RuntimeError(
            "ffill & bfill requires bottleneck or numbagg to be enabled."
            " Call `xr.set_options(use_bottleneck=True)` or `xr.set_options(use_numbagg=True)` to enable one."
        )
    if is_duck_dask_array(array):
        return dask_array_ops.push(array, n, axis, method=method)
    else:
        return _push(array, n, axis)


def _first_last_wrapper(array, *, axis, op, keepdims):
    return op(array, axis, keepdims=keepdims)


def _chunked_first_or_last(darray, axis, op):
    chunkmanager = get_chunked_array_type(darray)

    # This will raise the same error message seen for numpy
    axis = normalize_axis_index(axis, darray.ndim)

    wrapped_op = partial(_first_last_wrapper, op=op)
    return chunkmanager.reduction(
        darray,
        func=wrapped_op,
        aggregate_func=wrapped_op,
        axis=axis,
        dtype=darray.dtype,
        keepdims=False,  # match numpy version
    )


def chunked_nanfirst(darray, axis):
    return _chunked_first_or_last(darray, axis, op=nputils.nanfirst)


def chunked_nanlast(darray, axis):
    return _chunked_first_or_last(darray, axis, op=nputils.nanlast)
