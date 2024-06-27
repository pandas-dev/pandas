"""Compatibility module defining operations on duck numpy-arrays.

Currently, this means Dask or NumPy arrays. None of these functions should
accept or return xarray objects.
"""

from __future__ import annotations

import contextlib
import datetime
import inspect
import warnings
from functools import partial
from importlib import import_module

import numpy as np
import pandas as pd
from numpy import all as array_all  # noqa
from numpy import any as array_any  # noqa
from numpy import (  # noqa
    around,  # noqa
    full_like,
    gradient,
    isclose,
    isin,
    isnat,
    take,
    tensordot,
    transpose,
    unravel_index,
)
from numpy import concatenate as _concatenate
from numpy.lib.stride_tricks import sliding_window_view  # noqa
from packaging.version import Version
from pandas.api.types import is_extension_array_dtype

from xarray.core import dask_array_ops, dtypes, nputils
from xarray.core.options import OPTIONS
from xarray.core.utils import is_duck_array, is_duck_dask_array, module_available
from xarray.namedarray import pycompat
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


def get_array_namespace(*values):
    def _get_array_namespace(x):
        if hasattr(x, "__array_namespace__"):
            return x.__array_namespace__()
        else:
            return np

    namespaces = {_get_array_namespace(t) for t in values}
    non_numpy = namespaces - {np}

    if len(non_numpy) > 1:
        raise TypeError(
            "cannot deal with more than one type supporting the array API at the same time"
        )
    elif non_numpy:
        [xp] = non_numpy
    else:
        xp = np

    return xp


def einsum(*args, **kwargs):
    from xarray.core.options import OPTIONS

    if OPTIONS["use_opt_einsum"] and module_available("opt_einsum"):
        import opt_einsum

        return opt_einsum.contract(*args, **kwargs)
    else:
        return np.einsum(*args, **kwargs)


def _dask_or_eager_func(
    name,
    eager_module=np,
    dask_module="dask.array",
):
    """Create a function that dispatches to dask for dask array inputs."""

    def f(*args, **kwargs):
        if any(is_duck_dask_array(a) for a in args):
            mod = (
                import_module(dask_module)
                if isinstance(dask_module, str)
                else dask_module
            )
            wrapped = getattr(mod, name)
        else:
            wrapped = getattr(eager_module, name)
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

# np.around has failing doctests, overwrite it so they pass:
# https://github.com/numpy/numpy/issues/19759
around.__doc__ = str.replace(
    around.__doc__ or "",
    "array([0.,  2.])",
    "array([0., 2.])",
)
around.__doc__ = str.replace(
    around.__doc__ or "",
    "array([0.,  2.])",
    "array([0., 2.])",
)
around.__doc__ = str.replace(
    around.__doc__ or "",
    "array([0.4,  1.6])",
    "array([0.4, 1.6])",
)
around.__doc__ = str.replace(
    around.__doc__ or "",
    "array([0.,  2.,  2.,  4.,  4.])",
    "array([0., 2., 2., 4., 4.])",
)
around.__doc__ = str.replace(
    around.__doc__ or "",
    (
        '    .. [2] "How Futile are Mindless Assessments of\n'
        '           Roundoff in Floating-Point Computation?", William Kahan,\n'
        "           https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf\n"
    ),
    "",
)


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
        return full_like(data, dtype=bool, fill_value=False)
    else:
        # at this point, array should have dtype=object
        if isinstance(data, np.ndarray) or is_extension_array_dtype(data):
            return pandas_isnull(data)
        else:
            # Not reachable yet, but intended for use with other duck array
            # types. For full consistency with pandas, we should accept None as
            # a null value as well as NaN, but it isn't clear how to do this
            # with duck typing.
            return data != data


def notnull(data):
    return ~isnull(data)


# TODO replace with simply np.ma.masked_invalid once numpy/numpy#16022 is fixed
masked_invalid = _dask_or_eager_func(
    "masked_invalid", eager_module=np.ma, dask_module="dask.array.ma"
)


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
    integrand = np.pad(integrand, pads, mode="constant", constant_values=0.0)

    return cumsum(integrand, axis=axis, skipna=False)


def astype(data, dtype, **kwargs):
    if hasattr(data, "__array_namespace__"):
        xp = get_array_namespace(data)
        if xp == np:
            # numpy currently doesn't have a astype:
            return data.astype(dtype, **kwargs)
        return xp.astype(data, dtype, **kwargs)
    return data.astype(dtype, **kwargs)


def asarray(data, xp=np, dtype=None):
    converted = data if is_duck_array(data) else xp.asarray(data)

    if dtype is None or converted.dtype == dtype:
        return converted

    if xp is np or not hasattr(xp, "astype"):
        return converted.astype(dtype)
    else:
        return xp.astype(converted, dtype)


def as_shared_dtype(scalars_or_arrays, xp=None):
    """Cast a arrays to a shared dtype using xarray's type promotion rules."""
    if any(is_extension_array_dtype(x) for x in scalars_or_arrays):
        extension_array_types = [
            x.dtype for x in scalars_or_arrays if is_extension_array_dtype(x)
        ]
        if len(extension_array_types) == len(scalars_or_arrays) and all(
            isinstance(x, type(extension_array_types[0])) for x in extension_array_types
        ):
            return scalars_or_arrays
        raise ValueError(
            "Cannot cast arrays to shared type, found"
            f" array types {[x.dtype for x in scalars_or_arrays]}"
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
            return bool(isclose(arr1, arr2, rtol=rtol, atol=atol, equal_nan=True).all())
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
            return bool(flag_array.all())
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
            return bool(flag_array.all())
    else:
        return lazy_equiv


def count(data, axis=None):
    """Count the number of non-NA in this array along the given axis or axes"""
    return np.sum(np.logical_not(isnull(data)), axis=axis)


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
    xp = get_array_namespace(condition)
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


def concatenate(arrays, axis=0):
    """concatenate() with better dtype promotion rules."""
    # TODO: remove the additional check once `numpy` adds `concat` to its array namespace
    if hasattr(arrays[0], "__array_namespace__") and not isinstance(
        arrays[0], np.ndarray
    ):
        xp = get_array_namespace(arrays[0])
        return xp.concat(as_shared_dtype(arrays, xp=xp), axis=axis)
    return _concatenate(as_shared_dtype(arrays), axis=axis)


def stack(arrays, axis=0):
    """stack() with better dtype promotion rules."""
    xp = get_array_namespace(arrays[0])
    return xp.stack(as_shared_dtype(arrays, xp=xp), axis=axis)


def reshape(array, shape):
    xp = get_array_namespace(array)
    return xp.reshape(array, shape)


def ravel(array):
    return reshape(array, (-1,))


@contextlib.contextmanager
def _ignore_warnings_if(condition):
    if condition:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            yield
    else:
        yield


def _create_nan_agg_method(name, coerce_strings=False, invariant_0d=False):
    from xarray.core import nanops

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
            except (AttributeError, TypeError):
                raise NotImplementedError(
                    f"{name} is not yet implemented on dask arrays"
                )

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


_mean = _create_nan_agg_method("mean", invariant_0d=True)


def _datetime_nanmin(array):
    """nanmin() function for datetime64.

    Caveats that this function deals with:

    - In numpy < 1.18, min() on datetime64 incorrectly ignores NaT
    - numpy nanmin() don't work on datetime64 (all versions at the moment of writing)
    - dask min() does not work on datetime64 (all versions at the moment of writing)
    """
    dtype = array.dtype
    assert dtypes.is_datetime_like(dtype)
    # (NaT).astype(float) does not produce NaN...
    array = where(pandas_isnull(array), np.nan, array.astype(float))
    array = min(array, skipna=True)
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
            offset = _datetime_nanmin(array)
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
    import datetime as dt

    if isinstance(value, dt.timedelta):
        out = py_timedelta_to_float(value, datetime_unit)
    elif isinstance(value, np.timedelta64):
        out = np_timedelta64_to_float(value, datetime_unit)
    elif isinstance(value, pd.Timedelta):
        out = pd_timedelta_to_float(value, datetime_unit)
    elif isinstance(value, str):
        try:
            a = pd.to_timedelta(value)
        except ValueError:
            raise ValueError(
                f"Could not convert {value!r} to timedelta64 using pandas.to_timedelta"
            )
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
    """Convert numpy.timedelta64 to float.

    Notes
    -----
    The array is first converted to microseconds, which is less likely to
    cause overflow errors.
    """
    array = array.astype("timedelta64[ns]").astype(np.float64)
    conversion_factor = np.timedelta64(1, "ns") / np.timedelta64(1, datetime_unit)
    return conversion_factor * array


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
        offset = _datetime_nanmin(array)

        # xarray always uses np.datetime64[ns] for np.datetime64 data
        dtype = "timedelta64[ns]"
        return (
            _mean(
                datetime_to_numeric(array, offset), axis=axis, skipna=skipna, **kwargs
            ).astype(dtype)
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

        if pycompat.mod_version("numbagg") < Version("0.6.2"):
            warnings.warn(
                f"numbagg >= 0.6.2 is required for bfill & ffill; {pycompat.mod_version('numbagg')} is installed. We'll attempt with bottleneck instead."
            )
        else:
            return numbagg.ffill(array, limit=n, axis=axis)

    # work around for bottleneck 178
    limit = n if n is not None else array.shape[axis]

    import bottleneck as bn

    return bn.push(array, limit, axis)


def push(array, n, axis):
    if not OPTIONS["use_bottleneck"] and not OPTIONS["use_numbagg"]:
        raise RuntimeError(
            "ffill & bfill requires bottleneck or numbagg to be enabled."
            " Call `xr.set_options(use_bottleneck=True)` or `xr.set_options(use_numbagg=True)` to enable one."
        )
    if is_duck_dask_array(array):
        return dask_array_ops.push(array, n, axis)
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
