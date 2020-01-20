import functools
import itertools
import operator
from typing import Any, Optional, Tuple, Union

import numpy as np

from pandas._config import get_option

from pandas._libs import NaT, Timedelta, Timestamp, iNaT, lib
from pandas._typing import Dtype, Scalar
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.cast import _int64_max, maybe_upcast_putmask
from pandas.core.dtypes.common import (
    _get_dtype,
    is_any_int_dtype,
    is_bool_dtype,
    is_complex,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_datetime_or_timedelta_dtype,
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_numeric_dtype,
    is_object_dtype,
    is_scalar,
    is_timedelta64_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.missing import isna, na_value_for_dtype, notna

bn = import_optional_dependency("bottleneck", raise_on_missing=False, on_version="warn")
_BOTTLENECK_INSTALLED = bn is not None
_USE_BOTTLENECK = False


def set_use_bottleneck(v: bool = True) -> None:
    # set/unset to use bottleneck
    global _USE_BOTTLENECK
    if _BOTTLENECK_INSTALLED:
        _USE_BOTTLENECK = v


set_use_bottleneck(get_option("compute.use_bottleneck"))


class disallow:
    def __init__(self, *dtypes):
        super().__init__()
        self.dtypes = tuple(pandas_dtype(dtype).type for dtype in dtypes)

    def check(self, obj) -> bool:
        return hasattr(obj, "dtype") and issubclass(obj.dtype.type, self.dtypes)

    def __call__(self, f):
        @functools.wraps(f)
        def _f(*args, **kwargs):
            obj_iter = itertools.chain(args, kwargs.values())
            if any(self.check(obj) for obj in obj_iter):
                f_name = f.__name__.replace("nan", "")
                raise TypeError(
                    f"reduction operation '{f_name}' not allowed for this dtype"
                )
            try:
                with np.errstate(invalid="ignore"):
                    return f(*args, **kwargs)
            except ValueError as e:
                # we want to transform an object array
                # ValueError message to the more typical TypeError
                # e.g. this is normally a disallowed function on
                # object arrays that contain strings
                if is_object_dtype(args[0]):
                    raise TypeError(e)
                raise

        return _f


class bottleneck_switch:
    def __init__(self, name=None, **kwargs):
        self.name = name
        self.kwargs = kwargs

    def __call__(self, alt):
        bn_name = self.name or alt.__name__

        try:
            bn_func = getattr(bn, bn_name)
        except (AttributeError, NameError):  # pragma: no cover
            bn_func = None

        @functools.wraps(alt)
        def f(
            values: np.ndarray, axis: Optional[int] = None, skipna: bool = True, **kwds
        ):
            if len(self.kwargs) > 0:
                for k, v in self.kwargs.items():
                    if k not in kwds:
                        kwds[k] = v

            if values.size == 0 and kwds.get("min_count") is None:
                # We are empty, returning NA for our type
                # Only applies for the default `min_count` of None
                # since that affects how empty arrays are handled.
                # TODO(GH-18976) update all the nanops methods to
                # correctly handle empty inputs and remove this check.
                # It *may* just be `var`
                return _na_for_min_count(values, axis)

            if _USE_BOTTLENECK and skipna and _bn_ok_dtype(values.dtype, bn_name):
                if kwds.get("mask", None) is None:
                    # `mask` is not recognised by bottleneck, would raise
                    #  TypeError if called
                    kwds.pop("mask", None)
                    result = bn_func(values, axis=axis, **kwds)

                    # prefer to treat inf/-inf as NA, but must compute the func
                    # twice :(
                    if _has_infs(result):
                        result = alt(values, axis=axis, skipna=skipna, **kwds)
                else:
                    result = alt(values, axis=axis, skipna=skipna, **kwds)
            else:
                result = alt(values, axis=axis, skipna=skipna, **kwds)

            return result

        return f


def _bn_ok_dtype(dtype: Dtype, name: str) -> bool:
    # Bottleneck chokes on datetime64
    if not is_object_dtype(dtype) and not (
        is_datetime_or_timedelta_dtype(dtype) or is_datetime64tz_dtype(dtype)
    ):

        # GH 15507
        # bottleneck does not properly upcast during the sum
        # so can overflow

        # GH 9422
        # further we also want to preserve NaN when all elements
        # are NaN, unlinke bottleneck/numpy which consider this
        # to be 0
        if name in ["nansum", "nanprod"]:
            return False

        return True
    return False


def _has_infs(result) -> bool:
    if isinstance(result, np.ndarray):
        if result.dtype == "f8":
            return lib.has_infs_f8(result.ravel())
        elif result.dtype == "f4":
            return lib.has_infs_f4(result.ravel())
    try:
        return np.isinf(result).any()
    except (TypeError, NotImplementedError):
        # if it doesn't support infs, then it can't have infs
        return False


def _get_fill_value(
    dtype: Dtype, fill_value: Optional[Scalar] = None, fill_value_typ=None
):
    """ return the correct fill value for the dtype of the values """
    if fill_value is not None:
        return fill_value
    if _na_ok_dtype(dtype):
        if fill_value_typ is None:
            return np.nan
        else:
            if fill_value_typ == "+inf":
                return np.inf
            else:
                return -np.inf
    else:
        if fill_value_typ is None:
            return iNaT
        else:
            if fill_value_typ == "+inf":
                # need the max int here
                return _int64_max
            else:
                return iNaT


def _maybe_get_mask(
    values: np.ndarray, skipna: bool, mask: Optional[np.ndarray]
) -> Optional[np.ndarray]:
    """
    Compute a mask if and only if necessary.

    This function will compute a mask iff it is necessary. Otherwise,
    return the provided mask (potentially None) when a mask does not need to be
    computed.

    A mask is never necessary if the values array is of boolean or integer
    dtypes, as these are incapable of storing NaNs. If passing a NaN-capable
    dtype that is interpretable as either boolean or integer data (eg,
    timedelta64), a mask must be provided.

    If the skipna parameter is False, a new mask will not be computed.

    The mask is computed using isna() by default. Setting invert=True selects
    notna() as the masking function.

    Parameters
    ----------
    values : ndarray
        input array to potentially compute mask for
    skipna : bool
        boolean for whether NaNs should be skipped
    mask : Optional[ndarray]
        nan-mask if known

    Returns
    -------
    Optional[np.ndarray]
    """

    if mask is None:
        if is_bool_dtype(values.dtype) or is_integer_dtype(values.dtype):
            # Boolean data cannot contain nulls, so signal via mask being None
            return None

        if skipna:
            mask = isna(values)

    return mask


def _get_values(
    values: np.ndarray,
    skipna: bool,
    fill_value: Any = None,
    fill_value_typ: Optional[str] = None,
    mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, Optional[np.ndarray], np.dtype, np.dtype, Any]:
    """
    Utility to get the values view, mask, dtype, dtype_max, and fill_value.

    If both mask and fill_value/fill_value_typ are not None and skipna is True,
    the values array will be copied.

    For input arrays of boolean or integer dtypes, copies will only occur if a
    precomputed mask, a fill_value/fill_value_typ, and skipna=True are
    provided.

    Parameters
    ----------
    values : ndarray
        input array to potentially compute mask for
    skipna : bool
        boolean for whether NaNs should be skipped
    fill_value : Any
        value to fill NaNs with
    fill_value_typ : str
        Set to '+inf' or '-inf' to handle dtype-specific infinities
    mask : Optional[np.ndarray]
        nan-mask if known

    Returns
    -------
    values : ndarray
        Potential copy of input value array
    mask : Optional[ndarray[bool]]
        Mask for values, if deemed necessary to compute
    dtype : dtype
        dtype for values
    dtype_max : dtype
        platform independent dtype
    fill_value : Any
        fill value used
    """

    # In _get_values is only called from within nanops, and in all cases
    #  with scalar fill_value.  This guarantee is important for the
    #  maybe_upcast_putmask call below
    assert is_scalar(fill_value)

    mask = _maybe_get_mask(values, skipna, mask)

    if is_datetime64tz_dtype(values):
        # lib.values_from_object returns M8[ns] dtype instead of tz-aware,
        #  so this case must be handled separately from the rest
        dtype = values.dtype
        values = getattr(values, "_values", values)
    else:
        values = lib.values_from_object(values)
        dtype = values.dtype

    if is_datetime_or_timedelta_dtype(values) or is_datetime64tz_dtype(values):
        # changing timedelta64/datetime64 to int64 needs to happen after
        #  finding `mask` above
        values = getattr(values, "asi8", values)
        values = values.view(np.int64)

    dtype_ok = _na_ok_dtype(dtype)

    # get our fill value (in case we need to provide an alternative
    # dtype for it)
    fill_value = _get_fill_value(
        dtype, fill_value=fill_value, fill_value_typ=fill_value_typ
    )

    copy = (mask is not None) and (fill_value is not None)

    if skipna and copy:
        values = values.copy()
        if dtype_ok:
            np.putmask(values, mask, fill_value)

        # promote if needed
        else:
            values, _ = maybe_upcast_putmask(values, mask, fill_value)

    # return a platform independent precision dtype
    dtype_max = dtype
    if is_integer_dtype(dtype) or is_bool_dtype(dtype):
        dtype_max = np.int64
    elif is_float_dtype(dtype):
        dtype_max = np.float64

    return values, mask, dtype, dtype_max, fill_value


def _na_ok_dtype(dtype) -> bool:
    # TODO: what about datetime64tz?  PeriodDtype?
    return not issubclass(dtype.type, (np.integer, np.timedelta64, np.datetime64))


def _wrap_results(result, dtype: Dtype, fill_value=None):
    """ wrap our results if needed """

    if is_datetime64_dtype(dtype) or is_datetime64tz_dtype(dtype):
        if fill_value is None:
            # GH#24293
            fill_value = iNaT
        if not isinstance(result, np.ndarray):
            tz = getattr(dtype, "tz", None)
            assert not isna(fill_value), "Expected non-null fill_value"
            if result == fill_value:
                result = np.nan
            result = Timestamp(result, tz=tz)
        else:
            result = result.view(dtype)
    elif is_timedelta64_dtype(dtype):
        if not isinstance(result, np.ndarray):
            if result == fill_value:
                result = np.nan

            # raise if we have a timedelta64[ns] which is too large
            if np.fabs(result) > _int64_max:
                raise ValueError("overflow in timedelta operation")

            result = Timedelta(result, unit="ns")
        else:
            result = result.astype("m8[ns]").view(dtype)

    return result


def _na_for_min_count(
    values: np.ndarray, axis: Optional[int]
) -> Union[Scalar, np.ndarray]:
    """
    Return the missing value for `values`.

    Parameters
    ----------
    values : ndarray
    axis : int or None
        axis for the reduction, required if values.ndim > 1.

    Returns
    -------
    result : scalar or ndarray
        For 1-D values, returns a scalar of the correct missing type.
        For 2-D values, returns a 1-D array where each element is missing.
    """
    # we either return np.nan or pd.NaT
    if is_numeric_dtype(values):
        values = values.astype("float64")
    fill_value = na_value_for_dtype(values.dtype)

    if values.ndim == 1:
        return fill_value
    else:
        assert axis is not None  # assertion to make mypy happy
        result_shape = values.shape[:axis] + values.shape[axis + 1 :]
        result = np.empty(result_shape, dtype=values.dtype)
        result.fill(fill_value)
        return result


def nanany(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    mask: Optional[np.ndarray] = None,
) -> bool:
    """
    Check if any elements along an axis evaluate to True.

    Parameters
    ----------
    values : ndarray
    axis : int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : bool

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2])
    >>> nanops.nanany(s)
    True

    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([np.nan])
    >>> nanops.nanany(s)
    False
    """
    values, _, _, _, _ = _get_values(values, skipna, fill_value=False, mask=mask)
    return values.any(axis)


def nanall(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    mask: Optional[np.ndarray] = None,
) -> bool:
    """
    Check if all elements along an axis evaluate to True.

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : bool

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2, np.nan])
    >>> nanops.nanall(s)
    True

    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 0])
    >>> nanops.nanall(s)
    False
    """
    values, _, _, _, _ = _get_values(values, skipna, fill_value=True, mask=mask)
    return values.all(axis)


@disallow("M8")
def nansum(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    min_count: int = 0,
    mask: Optional[np.ndarray] = None,
) -> float:
    """
    Sum the elements along an axis ignoring NaNs

    Parameters
    ----------
    values : ndarray[dtype]
    axis: int, optional
    skipna : bool, default True
    min_count: int, default 0
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : dtype

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2, np.nan])
    >>> nanops.nansum(s)
    3.0
    """
    values, mask, dtype, dtype_max, _ = _get_values(
        values, skipna, fill_value=0, mask=mask
    )
    dtype_sum = dtype_max
    if is_float_dtype(dtype):
        dtype_sum = dtype
    elif is_timedelta64_dtype(dtype):
        dtype_sum = np.float64
    the_sum = values.sum(axis, dtype=dtype_sum)
    the_sum = _maybe_null_out(the_sum, axis, mask, values.shape, min_count=min_count)

    return _wrap_results(the_sum, dtype)


@disallow("M8", DatetimeTZDtype)
@bottleneck_switch()
def nanmean(values, axis=None, skipna=True, mask=None):
    """
    Compute the mean of the element along an axis ignoring NaNs

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2, np.nan])
    >>> nanops.nanmean(s)
    1.5
    """
    values, mask, dtype, dtype_max, _ = _get_values(
        values, skipna, fill_value=0, mask=mask
    )
    dtype_sum = dtype_max
    dtype_count = np.float64
    if (
        is_integer_dtype(dtype)
        or is_timedelta64_dtype(dtype)
        or is_datetime64_dtype(dtype)
        or is_datetime64tz_dtype(dtype)
    ):
        dtype_sum = np.float64
    elif is_float_dtype(dtype):
        dtype_sum = dtype
        dtype_count = dtype
    count = _get_counts(values.shape, mask, axis, dtype=dtype_count)
    the_sum = _ensure_numeric(values.sum(axis, dtype=dtype_sum))

    if axis is not None and getattr(the_sum, "ndim", False):
        with np.errstate(all="ignore"):
            # suppress division by zero warnings
            the_mean = the_sum / count
        ct_mask = count == 0
        if ct_mask.any():
            the_mean[ct_mask] = np.nan
    else:
        the_mean = the_sum / count if count > 0 else np.nan

    return _wrap_results(the_mean, dtype)


@disallow("M8")
@bottleneck_switch()
def nanmedian(values, axis=None, skipna=True, mask=None):
    """
    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, np.nan, 2, 2])
    >>> nanops.nanmedian(s)
    2.0
    """

    def get_median(x):
        mask = notna(x)
        if not skipna and not mask.all():
            return np.nan
        return np.nanmedian(x[mask])

    values, mask, dtype, dtype_max, _ = _get_values(values, skipna, mask=mask)
    if not is_float_dtype(values):
        values = values.astype("f8")
        if mask is not None:
            values[mask] = np.nan

    if axis is None:
        values = values.ravel()

    notempty = values.size

    # an array from a frame
    if values.ndim > 1:

        # there's a non-empty array to apply over otherwise numpy raises
        if notempty:
            if not skipna:
                return _wrap_results(
                    np.apply_along_axis(get_median, axis, values), dtype
                )

            # fastpath for the skipna case
            return _wrap_results(np.nanmedian(values, axis), dtype)

        # must return the correct shape, but median is not defined for the
        # empty set so return nans of shape "everything but the passed axis"
        # since "axis" is where the reduction would occur if we had a nonempty
        # array
        shp = np.array(values.shape)
        dims = np.arange(values.ndim)
        ret = np.empty(shp[dims != axis])
        ret.fill(np.nan)
        return _wrap_results(ret, dtype)

    # otherwise return a scalar value
    return _wrap_results(get_median(values) if notempty else np.nan, dtype)


def _get_counts_nanvar(
    value_counts: Tuple[int],
    mask: Optional[np.ndarray],
    axis: Optional[int],
    ddof: int,
    dtype: Dtype = float,
) -> Tuple[Union[int, np.ndarray], Union[int, np.ndarray]]:
    """ Get the count of non-null values along an axis, accounting
    for degrees of freedom.

    Parameters
    ----------
    values_shape : Tuple[int]
        shape tuple from values ndarray, used if mask is None
    mask : Optional[ndarray[bool]]
        locations in values that should be considered missing
    axis : Optional[int]
        axis to count along
    ddof : int
        degrees of freedom
    dtype : type, optional
        type to use for count

    Returns
    -------
    count : scalar or array
    d : scalar or array
    """
    dtype = _get_dtype(dtype)
    count = _get_counts(value_counts, mask, axis, dtype=dtype)
    d = count - dtype.type(ddof)

    # always return NaN, never inf
    if is_scalar(count):
        if count <= ddof:
            count = np.nan
            d = np.nan
    else:
        mask2: np.ndarray = count <= ddof
        if mask2.any():
            np.putmask(d, mask2, np.nan)
            np.putmask(count, mask2, np.nan)
    return count, d


@disallow("M8")
@bottleneck_switch(ddof=1)
def nanstd(values, axis=None, skipna=True, ddof=1, mask=None):
    """
    Compute the standard deviation along given axis while ignoring NaNs

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    ddof : int, default 1
        Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
        where N represents the number of elements.
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, np.nan, 2, 3])
    >>> nanops.nanstd(s)
    1.0
    """
    orig_dtype = values.dtype
    values, mask, dtype, dtype_max, fill_value = _get_values(values, skipna, mask=mask)

    result = np.sqrt(nanvar(values, axis=axis, skipna=skipna, ddof=ddof, mask=mask))
    return _wrap_results(result, orig_dtype)


@disallow("M8", "m8")
@bottleneck_switch(ddof=1)
def nanvar(values, axis=None, skipna=True, ddof=1, mask=None):
    """
    Compute the variance along given axis while ignoring NaNs

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    ddof : int, default 1
        Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
        where N represents the number of elements.
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, np.nan, 2, 3])
    >>> nanops.nanvar(s)
    1.0
    """
    values = lib.values_from_object(values)
    dtype = values.dtype
    mask = _maybe_get_mask(values, skipna, mask)
    if is_any_int_dtype(values):
        values = values.astype("f8")
        if mask is not None:
            values[mask] = np.nan

    if is_float_dtype(values):
        count, d = _get_counts_nanvar(values.shape, mask, axis, ddof, values.dtype)
    else:
        count, d = _get_counts_nanvar(values.shape, mask, axis, ddof)

    if skipna and mask is not None:
        values = values.copy()
        np.putmask(values, mask, 0)

    # xref GH10242
    # Compute variance via two-pass algorithm, which is stable against
    # cancellation errors and relatively accurate for small numbers of
    # observations.
    #
    # See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    avg = _ensure_numeric(values.sum(axis=axis, dtype=np.float64)) / count
    if axis is not None:
        avg = np.expand_dims(avg, axis)
    sqr = _ensure_numeric((avg - values) ** 2)
    if mask is not None:
        np.putmask(sqr, mask, 0)
    result = sqr.sum(axis=axis, dtype=np.float64) / d

    # Return variance as np.float64 (the datatype used in the accumulator),
    # unless we were dealing with a float array, in which case use the same
    # precision as the original values array.
    if is_float_dtype(dtype):
        result = result.astype(dtype)
    return _wrap_results(result, values.dtype)


@disallow("M8", "m8")
def nansem(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    ddof: int = 1,
    mask: Optional[np.ndarray] = None,
) -> float:
    """
    Compute the standard error in the mean along given axis while ignoring NaNs

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    ddof : int, default 1
        Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
        where N represents the number of elements.
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float64
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, np.nan, 2, 3])
    >>> nanops.nansem(s)
     0.5773502691896258
    """

    # This checks if non-numeric-like data is passed with numeric_only=False
    # and raises a TypeError otherwise
    nanvar(values, axis, skipna, ddof=ddof, mask=mask)

    mask = _maybe_get_mask(values, skipna, mask)
    if not is_float_dtype(values.dtype):
        values = values.astype("f8")

    count, _ = _get_counts_nanvar(values.shape, mask, axis, ddof, values.dtype)
    var = nanvar(values, axis, skipna, ddof=ddof)

    return np.sqrt(var) / np.sqrt(count)


def _nanminmax(meth, fill_value_typ):
    @bottleneck_switch(name="nan" + meth)
    def reduction(
        values: np.ndarray,
        axis: Optional[int] = None,
        skipna: bool = True,
        mask: Optional[np.ndarray] = None,
    ) -> Dtype:

        values, mask, dtype, dtype_max, fill_value = _get_values(
            values, skipna, fill_value_typ=fill_value_typ, mask=mask
        )

        if (axis is not None and values.shape[axis] == 0) or values.size == 0:
            try:
                result = getattr(values, meth)(axis, dtype=dtype_max)
                result.fill(np.nan)
            except (AttributeError, TypeError, ValueError):
                result = np.nan
        else:
            result = getattr(values, meth)(axis)

        result = _wrap_results(result, dtype, fill_value)
        return _maybe_null_out(result, axis, mask, values.shape)

    return reduction


nanmin = _nanminmax("min", fill_value_typ="+inf")
nanmax = _nanminmax("max", fill_value_typ="-inf")


@disallow("O")
def nanargmax(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    mask: Optional[np.ndarray] = None,
) -> int:
    """
    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : int
        The index of max value in specified axis or -1 in the NA case

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2, 3, np.nan, 4])
    >>> nanops.nanargmax(s)
    4
    """
    values, mask, dtype, _, _ = _get_values(
        values, True, fill_value_typ="-inf", mask=mask
    )
    result = values.argmax(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result


@disallow("O")
def nanargmin(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    mask: Optional[np.ndarray] = None,
) -> int:
    """
    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : int
        The index of min value in specified axis or -1 in the NA case

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2, 3, np.nan, 4])
    >>> nanops.nanargmin(s)
    0
    """
    values, mask, dtype, _, _ = _get_values(
        values, True, fill_value_typ="+inf", mask=mask
    )
    result = values.argmin(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result


@disallow("M8", "m8")
def nanskew(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    mask: Optional[np.ndarray] = None,
) -> float:
    """ Compute the sample skewness.

    The statistic computed here is the adjusted Fisher-Pearson standardized
    moment coefficient G1. The algorithm computes this coefficient directly
    from the second and third central moment.

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float64
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1,np.nan, 1, 2])
    >>> nanops.nanskew(s)
    1.7320508075688787
    """
    values = lib.values_from_object(values)
    mask = _maybe_get_mask(values, skipna, mask)
    if not is_float_dtype(values.dtype):
        values = values.astype("f8")
        count = _get_counts(values.shape, mask, axis)
    else:
        count = _get_counts(values.shape, mask, axis, dtype=values.dtype)

    if skipna and mask is not None:
        values = values.copy()
        np.putmask(values, mask, 0)

    mean = values.sum(axis, dtype=np.float64) / count
    if axis is not None:
        mean = np.expand_dims(mean, axis)

    adjusted = values - mean
    if skipna and mask is not None:
        np.putmask(adjusted, mask, 0)
    adjusted2 = adjusted ** 2
    adjusted3 = adjusted2 * adjusted
    m2 = adjusted2.sum(axis, dtype=np.float64)
    m3 = adjusted3.sum(axis, dtype=np.float64)

    # floating point error
    #
    # #18044 in _libs/windows.pyx calc_skew follow this behavior
    # to fix the fperr to treat m2 <1e-14 as zero
    m2 = _zero_out_fperr(m2)
    m3 = _zero_out_fperr(m3)

    with np.errstate(invalid="ignore", divide="ignore"):
        result = (count * (count - 1) ** 0.5 / (count - 2)) * (m3 / m2 ** 1.5)

    dtype = values.dtype
    if is_float_dtype(dtype):
        result = result.astype(dtype)

    if isinstance(result, np.ndarray):
        result = np.where(m2 == 0, 0, result)
        result[count < 3] = np.nan
        return result
    else:
        result = 0 if m2 == 0 else result
        if count < 3:
            return np.nan
        return result


@disallow("M8", "m8")
def nankurt(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    mask: Optional[np.ndarray] = None,
) -> float:
    """
    Compute the sample excess kurtosis

    The statistic computed here is the adjusted Fisher-Pearson standardized
    moment coefficient G2, computed directly from the second and fourth
    central moment.

    Parameters
    ----------
    values : ndarray
    axis: int, optional
    skipna : bool, default True
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    result : float64
        Unless input is a float array, in which case use the same
        precision as the input array.

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1,np.nan, 1, 3, 2])
    >>> nanops.nankurt(s)
    -1.2892561983471076
    """
    values = lib.values_from_object(values)
    mask = _maybe_get_mask(values, skipna, mask)
    if not is_float_dtype(values.dtype):
        values = values.astype("f8")
        count = _get_counts(values.shape, mask, axis)
    else:
        count = _get_counts(values.shape, mask, axis, dtype=values.dtype)

    if skipna and mask is not None:
        values = values.copy()
        np.putmask(values, mask, 0)

    mean = values.sum(axis, dtype=np.float64) / count
    if axis is not None:
        mean = np.expand_dims(mean, axis)

    adjusted = values - mean
    if skipna and mask is not None:
        np.putmask(adjusted, mask, 0)
    adjusted2 = adjusted ** 2
    adjusted4 = adjusted2 ** 2
    m2 = adjusted2.sum(axis, dtype=np.float64)
    m4 = adjusted4.sum(axis, dtype=np.float64)

    with np.errstate(invalid="ignore", divide="ignore"):
        adj = 3 * (count - 1) ** 2 / ((count - 2) * (count - 3))
        numer = count * (count + 1) * (count - 1) * m4
        denom = (count - 2) * (count - 3) * m2 ** 2

    # floating point error
    #
    # #18044 in _libs/windows.pyx calc_kurt follow this behavior
    # to fix the fperr to treat denom <1e-14 as zero
    numer = _zero_out_fperr(numer)
    denom = _zero_out_fperr(denom)

    if not isinstance(denom, np.ndarray):
        # if ``denom`` is a scalar, check these corner cases first before
        # doing division
        if count < 4:
            return np.nan
        if denom == 0:
            return 0

    with np.errstate(invalid="ignore", divide="ignore"):
        result = numer / denom - adj

    dtype = values.dtype
    if is_float_dtype(dtype):
        result = result.astype(dtype)

    if isinstance(result, np.ndarray):
        result = np.where(denom == 0, 0, result)
        result[count < 4] = np.nan

    return result


@disallow("M8", "m8")
def nanprod(
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    min_count: int = 0,
    mask: Optional[np.ndarray] = None,
) -> float:
    """
    Parameters
    ----------
    values : ndarray[dtype]
    axis: int, optional
    skipna : bool, default True
    min_count: int, default 0
    mask : ndarray[bool], optional
        nan-mask if known

    Returns
    -------
    Dtype
        The product of all elements on a given axis. ( NaNs are treated as 1)

    Examples
    --------
    >>> import pandas.core.nanops as nanops
    >>> s = pd.Series([1, 2, 3, np.nan])
    >>> nanops.nanprod(s)
    6.0
    """
    mask = _maybe_get_mask(values, skipna, mask)

    if skipna and mask is not None:
        values = values.copy()
        values[mask] = 1
    result = values.prod(axis)
    return _maybe_null_out(result, axis, mask, values.shape, min_count=min_count)


def _maybe_arg_null_out(
    result: np.ndarray, axis: Optional[int], mask: Optional[np.ndarray], skipna: bool
) -> Union[np.ndarray, int]:
    # helper function for nanargmin/nanargmax
    if mask is None:
        return result

    if axis is None or not getattr(result, "ndim", False):
        if skipna:
            if mask.all():
                result = -1
        else:
            if mask.any():
                result = -1
    else:
        if skipna:
            na_mask = mask.all(axis)
        else:
            na_mask = mask.any(axis)
        if na_mask.any():
            result[na_mask] = -1
    return result


def _get_counts(
    values_shape: Tuple[int],
    mask: Optional[np.ndarray],
    axis: Optional[int],
    dtype: Dtype = float,
) -> Union[int, np.ndarray]:
    """ Get the count of non-null values along an axis

    Parameters
    ----------
    values_shape : Tuple[int]
        shape tuple from values ndarray, used if mask is None
    mask : Optional[ndarray[bool]]
        locations in values that should be considered missing
    axis : Optional[int]
        axis to count along
    dtype : type, optional
        type to use for count

    Returns
    -------
    count : scalar or array
    """
    dtype = _get_dtype(dtype)
    if axis is None:
        if mask is not None:
            n = mask.size - mask.sum()
        else:
            n = np.prod(values_shape)
        return dtype.type(n)

    if mask is not None:
        count = mask.shape[axis] - mask.sum(axis)
    else:
        count = values_shape[axis]

    if is_scalar(count):
        return dtype.type(count)
    try:
        return count.astype(dtype)
    except AttributeError:
        return np.array(count, dtype=dtype)


def _maybe_null_out(
    result: np.ndarray,
    axis: Optional[int],
    mask: Optional[np.ndarray],
    shape: Tuple,
    min_count: int = 1,
) -> float:
    """
    Returns
    -------
    Dtype
        The product of all elements on a given axis. ( NaNs are treated as 1)
    """
    if mask is not None and axis is not None and getattr(result, "ndim", False):
        null_mask = (mask.shape[axis] - mask.sum(axis) - min_count) < 0
        if np.any(null_mask):
            if is_numeric_dtype(result):
                if np.iscomplexobj(result):
                    result = result.astype("c16")
                else:
                    result = result.astype("f8")
                result[null_mask] = np.nan
            else:
                # GH12941, use None to auto cast null
                result[null_mask] = None
    elif result is not NaT:
        if mask is not None:
            null_mask = mask.size - mask.sum()
        else:
            null_mask = np.prod(shape)
        if null_mask < min_count:
            result = np.nan

    return result


def _zero_out_fperr(arg):
    # #18044 reference this behavior to fix rolling skew/kurt issue
    if isinstance(arg, np.ndarray):
        with np.errstate(invalid="ignore"):
            return np.where(np.abs(arg) < 1e-14, 0, arg)
    else:
        return arg.dtype.type(0) if np.abs(arg) < 1e-14 else arg


@disallow("M8", "m8")
def nancorr(
    a: np.ndarray, b: np.ndarray, method="pearson", min_periods: Optional[int] = None,
):
    """
    a, b: ndarrays
    """
    if len(a) != len(b):
        raise AssertionError("Operands to nancorr must have same size")

    if min_periods is None:
        min_periods = 1

    valid = notna(a) & notna(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

    if len(a) < min_periods:
        return np.nan

    f = get_corr_func(method)
    return f(a, b)


def get_corr_func(method):
    if method in ["kendall", "spearman"]:
        from scipy.stats import kendalltau, spearmanr
    elif method in ["pearson"]:
        pass
    elif callable(method):
        return method
    else:
        raise ValueError(
            f"Unkown method '{method}', expected one of 'kendall', 'spearman'"
        )

    def _pearson(a, b):
        return np.corrcoef(a, b)[0, 1]

    def _kendall(a, b):
        # kendallttau returns a tuple of the tau statistic and pvalue
        rs = kendalltau(a, b)
        return rs[0]

    def _spearman(a, b):
        return spearmanr(a, b)[0]

    _cor_methods = {"pearson": _pearson, "kendall": _kendall, "spearman": _spearman}
    return _cor_methods[method]


@disallow("M8", "m8")
def nancov(a: np.ndarray, b: np.ndarray, min_periods: Optional[int] = None):
    if len(a) != len(b):
        raise AssertionError("Operands to nancov must have same size")

    if min_periods is None:
        min_periods = 1

    valid = notna(a) & notna(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

    if len(a) < min_periods:
        return np.nan

    return np.cov(a, b)[0, 1]


def _ensure_numeric(x):
    if isinstance(x, np.ndarray):
        if is_integer_dtype(x) or is_bool_dtype(x):
            x = x.astype(np.float64)
        elif is_object_dtype(x):
            try:
                x = x.astype(np.complex128)
            except (TypeError, ValueError):
                x = x.astype(np.float64)
            else:
                if not np.any(np.imag(x)):
                    x = x.real
    elif not (is_float(x) or is_integer(x) or is_complex(x)):
        try:
            x = float(x)
        except ValueError:
            # e.g. "1+1j" or "foo"
            try:
                x = complex(x)
            except ValueError:
                # e.g. "foo"
                raise TypeError(f"Could not convert {x} to numeric")
    return x


# NA-friendly array comparisons


def make_nancomp(op):
    def f(x, y):
        xmask = isna(x)
        ymask = isna(y)
        mask = xmask | ymask

        with np.errstate(all="ignore"):
            result = op(x, y)

        if mask.any():
            if is_bool_dtype(result):
                result = result.astype("O")
            np.putmask(result, mask, np.nan)

        return result

    return f


nangt = make_nancomp(operator.gt)
nange = make_nancomp(operator.ge)
nanlt = make_nancomp(operator.lt)
nanle = make_nancomp(operator.le)
naneq = make_nancomp(operator.eq)
nanne = make_nancomp(operator.ne)


def _nanpercentile_1d(
    values: np.ndarray, mask: np.ndarray, q, na_value: Scalar, interpolation
) -> Union[Scalar, np.ndarray]:
    """
    Wrapper for np.percentile that skips missing values, specialized to
    1-dimensional case.

    Parameters
    ----------
    values : array over which to find quantiles
    mask : ndarray[bool]
        locations in values that should be considered missing
    q : scalar or array of quantile indices to find
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
        if lib.is_scalar(q):
            return na_value
        else:
            return np.array([na_value] * len(q), dtype=values.dtype)

    return np.percentile(values, q, interpolation=interpolation)


def nanpercentile(
    values: np.ndarray,
    q,
    axis: int,
    na_value,
    mask: np.ndarray,
    ndim: int,
    interpolation,
):
    """
    Wrapper for np.percentile that skips missing values.

    Parameters
    ----------
    values : array over which to find quantiles
    q : scalar or array of quantile indices to find
    axis : {0, 1}
    na_value : scalar
        value to return for empty or all-null values
    mask : ndarray[bool]
        locations in values that should be considered missing
    ndim : {1, 2}
    interpolation : str

    Returns
    -------
    quantiles : scalar or array
    """
    if values.dtype.kind in ["m", "M"]:
        # need to cast to integer to avoid rounding errors in numpy
        result = nanpercentile(
            values.view("i8"), q, axis, na_value.view("i8"), mask, ndim, interpolation
        )

        # Note: we have to do do `astype` and not view because in general we
        #  have float result at this point, not i8
        return result.astype(values.dtype)

    if not lib.is_scalar(mask) and mask.any():
        if ndim == 1:
            return _nanpercentile_1d(
                values, mask, q, na_value, interpolation=interpolation
            )
        else:
            # for nonconsolidatable blocks mask is 1D, but values 2D
            if mask.ndim < values.ndim:
                mask = mask.reshape(values.shape)
            if axis == 0:
                values = values.T
                mask = mask.T
            result = [
                _nanpercentile_1d(val, m, q, na_value, interpolation=interpolation)
                for (val, m) in zip(list(values), list(mask))
            ]
            result = np.array(result, dtype=values.dtype, copy=False).T
            return result
    else:
        return np.percentile(values, q, axis=axis, interpolation=interpolation)
