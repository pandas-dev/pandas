""" routings for casting """

from datetime import datetime, timedelta

import numpy as np

from pandas._libs import lib, tslib, tslibs
from pandas._libs.tslibs import NaT, OutOfBoundsDatetime, Period, iNaT
from pandas._libs.tslibs.timezones import tz_compare
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.common import (
    _INT64_DTYPE,
    _NS_DTYPE,
    _POSSIBLY_CAST_DTYPES,
    _TD_DTYPE,
    ensure_int8,
    ensure_int16,
    ensure_int32,
    ensure_int64,
    ensure_object,
    ensure_str,
    is_bool,
    is_bool_dtype,
    is_complex,
    is_complex_dtype,
    is_datetime64_dtype,
    is_datetime64_ns_dtype,
    is_datetime64tz_dtype,
    is_datetime_or_timedelta_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_object_dtype,
    is_scalar,
    is_string_dtype,
    is_timedelta64_dtype,
    is_timedelta64_ns_dtype,
    is_unsigned_integer_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import (
    DatetimeTZDtype,
    ExtensionDtype,
    IntervalDtype,
    PeriodDtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeArray,
    ABCDatetimeIndex,
    ABCPeriodArray,
    ABCPeriodIndex,
    ABCSeries,
)
from pandas.core.dtypes.inference import is_list_like
from pandas.core.dtypes.missing import isna, notna

_int8_max = np.iinfo(np.int8).max
_int16_max = np.iinfo(np.int16).max
_int32_max = np.iinfo(np.int32).max
_int64_max = np.iinfo(np.int64).max


def maybe_convert_platform(values):
    """ try to do platform conversion, allow ndarray or list here """

    if isinstance(values, (list, tuple, range)):
        values = construct_1d_object_array_from_listlike(values)
    if getattr(values, "dtype", None) == np.object_:
        if hasattr(values, "_values"):
            values = values._values
        values = lib.maybe_convert_objects(values)

    return values


def is_nested_object(obj) -> bool:
    """
    return a boolean if we have a nested object, e.g. a Series with 1 or
    more Series elements

    This may not be necessarily be performant.

    """

    if isinstance(obj, ABCSeries) and is_object_dtype(obj):

        if any(isinstance(v, ABCSeries) for v in obj.values):
            return True

    return False


def maybe_downcast_to_dtype(result, dtype):
    """ try to cast to the specified dtype (e.g. convert back to bool/int
    or could be an astype of float64->float32
    """
    do_round = False

    if is_scalar(result):
        return result
    elif isinstance(result, ABCDataFrame):
        # occurs in pivot_table doctest
        return result

    if isinstance(dtype, str):
        if dtype == "infer":
            inferred_type = lib.infer_dtype(ensure_object(result.ravel()), skipna=False)
            if inferred_type == "boolean":
                dtype = "bool"
            elif inferred_type == "integer":
                dtype = "int64"
            elif inferred_type == "datetime64":
                dtype = "datetime64[ns]"
            elif inferred_type == "timedelta64":
                dtype = "timedelta64[ns]"

            # try to upcast here
            elif inferred_type == "floating":
                dtype = "int64"
                if issubclass(result.dtype.type, np.number):
                    do_round = True

            else:
                dtype = "object"

        dtype = np.dtype(dtype)

    converted = maybe_downcast_numeric(result, dtype, do_round)
    if converted is not result:
        return converted

    # a datetimelike
    # GH12821, iNaT is casted to float
    if dtype.kind in ["M", "m"] and result.dtype.kind in ["i", "f"]:
        if hasattr(dtype, "tz"):
            # not a numpy dtype
            if dtype.tz:
                # convert to datetime and change timezone
                from pandas import to_datetime

                result = to_datetime(result).tz_localize("utc")
                result = result.tz_convert(dtype.tz)
        else:
            result = result.astype(dtype)

    elif dtype.type is Period:
        # TODO(DatetimeArray): merge with previous elif
        from pandas.core.arrays import PeriodArray

        try:
            return PeriodArray(result, freq=dtype.freq)
        except TypeError:
            # e.g. TypeError: int() argument must be a string, a
            #  bytes-like object or a number, not 'Period
            pass

    return result


def maybe_downcast_numeric(result, dtype, do_round: bool = False):
    """
    Subset of maybe_downcast_to_dtype restricted to numeric dtypes.

    Parameters
    ----------
    result : ndarray or ExtensionArray
    dtype : np.dtype or ExtensionDtype
    do_round : bool

    Returns
    -------
    ndarray or ExtensionArray
    """
    if not isinstance(dtype, np.dtype):
        # e.g. SparseDtype has no itemsize attr
        return result

    if isinstance(result, list):
        # reached via groupoby.agg _ohlc; really this should be handled
        #  earlier
        result = np.array(result)

    def trans(x):
        if do_round:
            return x.round()
        return x

    if dtype.kind == result.dtype.kind:
        # don't allow upcasts here (except if empty)
        if result.dtype.itemsize <= dtype.itemsize and result.size:
            return result

    if is_bool_dtype(dtype) or is_integer_dtype(dtype):

        if not result.size:
            # if we don't have any elements, just astype it
            return trans(result).astype(dtype)

        # do a test on the first element, if it fails then we are done
        r = result.ravel()
        arr = np.array([r[0]])

        if isna(arr).any():
            # if we have any nulls, then we are done
            return result

        elif not isinstance(r[0], (np.integer, np.floating, np.bool, int, float, bool)):
            # a comparable, e.g. a Decimal may slip in here
            return result

        if (
            issubclass(result.dtype.type, (np.object_, np.number))
            and notna(result).all()
        ):
            new_result = trans(result).astype(dtype)
            if new_result.dtype.kind == "O" or result.dtype.kind == "O":
                # np.allclose may raise TypeError on object-dtype
                if (new_result == result).all():
                    return new_result
            else:
                if np.allclose(new_result, result, rtol=0):
                    return new_result

    elif (
        issubclass(dtype.type, np.floating)
        and not is_bool_dtype(result.dtype)
        and not is_string_dtype(result.dtype)
    ):
        return result.astype(dtype)

    return result


def maybe_upcast_putmask(result: np.ndarray, mask: np.ndarray, other):
    """
    A safe version of putmask that potentially upcasts the result.
    The result is replaced with the first N elements of other,
    where N is the number of True values in mask.
    If the length of other is shorter than N, other will be repeated.

    Parameters
    ----------
    result : ndarray
        The destination array. This will be mutated in-place if no upcasting is
        necessary.
    mask : boolean ndarray
    other : scalar
        The source value.

    Returns
    -------
    result : ndarray
    changed : bool
        Set to true if the result array was upcasted.

    Examples
    --------
    >>> result, _ = maybe_upcast_putmask(np.arange(1,6),
    np.array([False, True, False, True, True]), np.arange(21,23))
    >>> result
    array([1, 21, 3, 22, 21])
    """

    if not isinstance(result, np.ndarray):
        raise ValueError("The result input must be a ndarray.")
    if not is_scalar(other):
        # We _could_ support non-scalar other, but until we have a compelling
        #  use case, we assume away the possibility.
        raise ValueError("other must be a scalar")

    if mask.any():
        # Two conversions for date-like dtypes that can't be done automatically
        # in np.place:
        #   NaN -> NaT
        #   integer or integer array -> date-like array
        if result.dtype.kind in ["m", "M"]:
            if is_scalar(other):
                if isna(other):
                    other = result.dtype.type("nat")
                elif is_integer(other):
                    other = np.array(other, dtype=result.dtype)
            elif is_integer_dtype(other):
                other = np.array(other, dtype=result.dtype)

        def changeit():

            # try to directly set by expanding our array to full
            # length of the boolean
            try:
                om = other[mask]
            except (IndexError, TypeError):
                # IndexError occurs in test_upcast when we have a boolean
                #  mask of the wrong shape
                # TypeError occurs in test_upcast when `other` is a bool
                pass
            else:
                om_at = om.astype(result.dtype)
                if (om == om_at).all():
                    new_result = result.values.copy()
                    new_result[mask] = om_at
                    result[:] = new_result
                    return result, False

            # we are forced to change the dtype of the result as the input
            # isn't compatible
            r, _ = maybe_upcast(result, fill_value=other, copy=True)
            np.place(r, mask, other)

            return r, True

        # we want to decide whether place will work
        # if we have nans in the False portion of our mask then we need to
        # upcast (possibly), otherwise we DON't want to upcast (e.g. if we
        # have values, say integers, in the success portion then it's ok to not
        # upcast)
        new_dtype, _ = maybe_promote(result.dtype, other)
        if new_dtype != result.dtype:

            # we have a scalar or len 0 ndarray
            # and its nan and we are changing some values
            if is_scalar(other) or (isinstance(other, np.ndarray) and other.ndim < 1):
                if isna(other):
                    return changeit()

            # we have an ndarray and the masking has nans in it
            else:

                if isna(other).any():
                    return changeit()

        try:
            np.place(result, mask, other)
        except TypeError:
            # e.g. int-dtype result and float-dtype other
            return changeit()

    return result, False


def maybe_promote(dtype, fill_value=np.nan):
    """
    Find the minimal dtype that can hold both the given dtype and fill_value.

    Parameters
    ----------
    dtype : np.dtype or ExtensionDtype
    fill_value : scalar, default np.nan

    Returns
    -------
    dtype
        Upcasted from dtype argument if necessary.
    fill_value
        Upcasted from fill_value argument if necessary.
    """
    if not is_scalar(fill_value) and not is_object_dtype(dtype):
        # with object dtype there is nothing to promote, and the user can
        #  pass pretty much any weird fill_value they like
        raise ValueError("fill_value must be a scalar")

    # if we passed an array here, determine the fill value by dtype
    if isinstance(fill_value, np.ndarray):
        if issubclass(fill_value.dtype.type, (np.datetime64, np.timedelta64)):
            fill_value = fill_value.dtype.type("NaT", "ns")
        else:

            # we need to change to object type as our
            # fill_value is of object type
            if fill_value.dtype == np.object_:
                dtype = np.dtype(np.object_)
            fill_value = np.nan

        if dtype == np.object_ or dtype.kind in ["U", "S"]:
            # We treat string-like dtypes as object, and _always_ fill
            #  with np.nan
            fill_value = np.nan
            dtype = np.dtype(np.object_)

    # returns tuple of (dtype, fill_value)
    if issubclass(dtype.type, np.datetime64):
        if isinstance(fill_value, datetime) and fill_value.tzinfo is not None:
            # Trying to insert tzaware into tznaive, have to cast to object
            dtype = np.dtype(np.object_)
        elif is_integer(fill_value) or (is_float(fill_value) and not isna(fill_value)):
            dtype = np.dtype(np.object_)
        else:
            try:
                fill_value = tslibs.Timestamp(fill_value).to_datetime64()
            except (TypeError, ValueError):
                dtype = np.dtype(np.object_)
    elif issubclass(dtype.type, np.timedelta64):
        if (
            is_integer(fill_value)
            or (is_float(fill_value) and not np.isnan(fill_value))
            or isinstance(fill_value, str)
        ):
            # TODO: What about str that can be a timedelta?
            dtype = np.dtype(np.object_)
        else:
            try:
                fv = tslibs.Timedelta(fill_value)
            except ValueError:
                dtype = np.dtype(np.object_)
            else:
                if fv is NaT:
                    # NaT has no `to_timedelta64` method
                    fill_value = np.timedelta64("NaT", "ns")
                else:
                    fill_value = fv.to_timedelta64()
    elif is_datetime64tz_dtype(dtype):
        if isna(fill_value):
            fill_value = NaT
        elif not isinstance(fill_value, datetime):
            dtype = np.dtype(np.object_)
        elif fill_value.tzinfo is None:
            dtype = np.dtype(np.object_)
        elif not tz_compare(fill_value.tzinfo, dtype.tz):
            # TODO: sure we want to cast here?
            dtype = np.dtype(np.object_)

    elif is_extension_array_dtype(dtype) and isna(fill_value):
        fill_value = dtype.na_value

    elif is_float(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.dtype(np.object_)

        elif issubclass(dtype.type, np.integer):
            dtype = np.dtype(np.float64)

        elif dtype.kind == "f":
            mst = np.min_scalar_type(fill_value)
            if mst > dtype:
                # e.g. mst is np.float64 and dtype is np.float32
                dtype = mst

        elif dtype.kind == "c":
            mst = np.min_scalar_type(fill_value)
            dtype = np.promote_types(dtype, mst)

    elif is_bool(fill_value):
        if not issubclass(dtype.type, np.bool_):
            dtype = np.dtype(np.object_)

    elif is_integer(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.dtype(np.object_)

        elif issubclass(dtype.type, np.integer):
            if not np.can_cast(fill_value, dtype):
                # upcast to prevent overflow
                mst = np.min_scalar_type(fill_value)
                dtype = np.promote_types(dtype, mst)
                if dtype.kind == "f":
                    # Case where we disagree with numpy
                    dtype = np.dtype(np.object_)

    elif is_complex(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.dtype(np.object_)

        elif issubclass(dtype.type, (np.integer, np.floating)):
            mst = np.min_scalar_type(fill_value)
            dtype = np.promote_types(dtype, mst)

        elif dtype.kind == "c":
            mst = np.min_scalar_type(fill_value)
            if mst > dtype:
                # e.g. mst is np.complex128 and dtype is np.complex64
                dtype = mst

    elif fill_value is None:
        if is_float_dtype(dtype) or is_complex_dtype(dtype):
            fill_value = np.nan
        elif is_integer_dtype(dtype):
            dtype = np.float64
            fill_value = np.nan
        elif is_datetime_or_timedelta_dtype(dtype):
            fill_value = dtype.type("NaT", "ns")
        else:
            dtype = np.dtype(np.object_)
            fill_value = np.nan
    else:
        dtype = np.dtype(np.object_)

    # in case we have a string that looked like a number
    if is_extension_array_dtype(dtype):
        pass
    elif issubclass(np.dtype(dtype).type, (bytes, str)):
        dtype = np.dtype(np.object_)

    fill_value = _ensure_dtype_type(fill_value, dtype)
    return dtype, fill_value


def _ensure_dtype_type(value, dtype):
    """
    Ensure that the given value is an instance of the given dtype.

    e.g. if out dtype is np.complex64, we should have an instance of that
    as opposed to a python complex object.

    Parameters
    ----------
    value : object
    dtype : np.dtype or ExtensionDtype

    Returns
    -------
    object
    """

    # Start with exceptions in which we do _not_ cast to numpy types
    if is_extension_array_dtype(dtype):
        return value
    elif dtype == np.object_:
        return value
    elif isna(value):
        # e.g. keep np.nan rather than try to cast to np.float32(np.nan)
        return value

    return dtype.type(value)


def infer_dtype_from(val, pandas_dtype: bool = False):
    """
    Interpret the dtype from a scalar or array.

    Parameters
    ----------
    val : object
    pandas_dtype : bool, default False
        whether to infer dtype including pandas extension types.
        If False, scalar/array belongs to pandas extension types is inferred as
        object
    """
    if is_scalar(val):
        return infer_dtype_from_scalar(val, pandas_dtype=pandas_dtype)
    return infer_dtype_from_array(val, pandas_dtype=pandas_dtype)


def infer_dtype_from_scalar(val, pandas_dtype: bool = False):
    """
    Interpret the dtype from a scalar.

    Parameters
    ----------
    pandas_dtype : bool, default False
        whether to infer dtype including pandas extension types.
        If False, scalar belongs to pandas extension types is inferred as
        object
    """

    dtype = np.object_

    # a 1-element ndarray
    if isinstance(val, np.ndarray):
        msg = "invalid ndarray passed to infer_dtype_from_scalar"
        if val.ndim != 0:
            raise ValueError(msg)

        dtype = val.dtype
        val = val.item()

    elif isinstance(val, str):

        # If we create an empty array using a string to infer
        # the dtype, NumPy will only allocate one character per entry
        # so this is kind of bad. Alternately we could use np.repeat
        # instead of np.empty (but then you still don't want things
        # coming out as np.str_!

        dtype = np.object_

    elif isinstance(val, (np.datetime64, datetime)):
        val = tslibs.Timestamp(val)
        if val is tslibs.NaT or val.tz is None:
            dtype = np.dtype("M8[ns]")
        else:
            if pandas_dtype:
                dtype = DatetimeTZDtype(unit="ns", tz=val.tz)
            else:
                # return datetimetz as object
                return np.object_, val
        val = val.value

    elif isinstance(val, (np.timedelta64, timedelta)):
        val = tslibs.Timedelta(val).value
        dtype = np.dtype("m8[ns]")

    elif is_bool(val):
        dtype = np.bool_

    elif is_integer(val):
        if isinstance(val, np.integer):
            dtype = type(val)
        else:
            dtype = np.int64

    elif is_float(val):
        if isinstance(val, np.floating):
            dtype = type(val)
        else:
            dtype = np.float64

    elif is_complex(val):
        dtype = np.complex_

    elif pandas_dtype:
        if lib.is_period(val):
            dtype = PeriodDtype(freq=val.freq)
            val = val.ordinal
        elif lib.is_interval(val):
            subtype = infer_dtype_from_scalar(val.left, pandas_dtype=True)[0]
            dtype = IntervalDtype(subtype=subtype)

    return dtype, val


def infer_dtype_from_array(arr, pandas_dtype: bool = False):
    """
    Infer the dtype from an array.

    Parameters
    ----------
    arr : array
    pandas_dtype : bool, default False
        whether to infer dtype including pandas extension types.
        If False, array belongs to pandas extension types
        is inferred as object

    Returns
    -------
    tuple (numpy-compat/pandas-compat dtype, array)

    Notes
    -----
    if pandas_dtype=False. these infer to numpy dtypes
    exactly with the exception that mixed / object dtypes
    are not coerced by stringifying or conversion

    if pandas_dtype=True. datetime64tz-aware/categorical
    types will retain there character.

    Examples
    --------
    >>> np.asarray([1, '1'])
    array(['1', '1'], dtype='<U21')

    >>> infer_dtype_from_array([1, '1'])
    (numpy.object_, [1, '1'])
    """

    if isinstance(arr, np.ndarray):
        return arr.dtype, arr

    if not is_list_like(arr):
        arr = [arr]

    if pandas_dtype and is_extension_array_dtype(arr):
        return arr.dtype, arr

    elif isinstance(arr, ABCSeries):
        return arr.dtype, np.asarray(arr)

    # don't force numpy coerce with nan's
    inferred = lib.infer_dtype(arr, skipna=False)
    if inferred in ["string", "bytes", "mixed", "mixed-integer"]:
        return (np.object_, arr)

    arr = np.asarray(arr)
    return arr.dtype, arr


def maybe_infer_dtype_type(element):
    """
    Try to infer an object's dtype, for use in arithmetic ops.

    Uses `element.dtype` if that's available.
    Objects implementing the iterator protocol are cast to a NumPy array,
    and from there the array's type is used.

    Parameters
    ----------
    element : object
        Possibly has a `.dtype` attribute, and possibly the iterator
        protocol.

    Returns
    -------
    tipo : type

    Examples
    --------
    >>> from collections import namedtuple
    >>> Foo = namedtuple("Foo", "dtype")
    >>> maybe_infer_dtype_type(Foo(np.dtype("i8")))
    numpy.int64
    """
    tipo = None
    if hasattr(element, "dtype"):
        tipo = element.dtype
    elif is_list_like(element):
        element = np.asarray(element)
        tipo = element.dtype
    return tipo


def maybe_upcast(values, fill_value=np.nan, dtype=None, copy: bool = False):
    """
    Provide explicit type promotion and coercion.

    Parameters
    ----------
    values : ndarray or ExtensionArray
        The array that we want to maybe upcast.
    fill_value : what we want to fill with
    dtype : if None, then use the dtype of the values, else coerce to this type
    copy : bool, default True
        If True always make a copy even if no upcast is required.
    """
    if not is_scalar(fill_value) and not is_object_dtype(values.dtype):
        # We allow arbitrary fill values for object dtype
        raise ValueError("fill_value must be a scalar")

    if is_extension_array_dtype(values):
        if copy:
            values = values.copy()
    else:
        if dtype is None:
            dtype = values.dtype
        new_dtype, fill_value = maybe_promote(dtype, fill_value)
        if new_dtype != values.dtype:
            values = values.astype(new_dtype)
        elif copy:
            values = values.copy()

    return values, fill_value


def invalidate_string_dtypes(dtype_set):
    """Change string like dtypes to object for
    ``DataFrame.select_dtypes()``.
    """
    non_string_dtypes = dtype_set - {np.dtype("S").type, np.dtype("<U").type}
    if non_string_dtypes != dtype_set:
        raise TypeError("string dtypes are not allowed, use 'object' instead")


def coerce_indexer_dtype(indexer, categories):
    """ coerce the indexer input array to the smallest dtype possible """
    length = len(categories)
    if length < _int8_max:
        return ensure_int8(indexer)
    elif length < _int16_max:
        return ensure_int16(indexer)
    elif length < _int32_max:
        return ensure_int32(indexer)
    return ensure_int64(indexer)


def coerce_to_dtypes(result, dtypes):
    """
    given a dtypes and a result set, coerce the result elements to the
    dtypes
    """
    if len(result) != len(dtypes):
        raise AssertionError("_coerce_to_dtypes requires equal len arrays")

    def conv(r, dtype):
        if np.any(isna(r)):
            pass
        elif dtype == _NS_DTYPE:
            r = tslibs.Timestamp(r)
        elif dtype == _TD_DTYPE:
            r = tslibs.Timedelta(r)
        elif dtype == np.bool_:
            # messy. non 0/1 integers do not get converted.
            if is_integer(r) and r not in [0, 1]:
                return int(r)
            r = bool(r)
        elif dtype.kind == "f":
            r = float(r)
        elif dtype.kind == "i":
            r = int(r)

        return r

    return [conv(r, dtype) for r, dtype in zip(result, dtypes)]


def astype_nansafe(arr, dtype, copy: bool = True, skipna: bool = False):
    """
    Cast the elements of an array to a given dtype a nan-safe manner.

    Parameters
    ----------
    arr : ndarray
    dtype : np.dtype
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
    if is_extension_array_dtype(dtype):
        return dtype.construct_array_type()._from_sequence(arr, dtype=dtype, copy=copy)

    if not isinstance(dtype, np.dtype):
        dtype = pandas_dtype(dtype)

    if issubclass(dtype.type, str):
        return lib.astype_str(arr.ravel(), skipna=skipna).reshape(arr.shape)

    elif is_datetime64_dtype(arr):
        if is_object_dtype(dtype):
            return tslib.ints_to_pydatetime(arr.view(np.int64))
        elif dtype == np.int64:
            if isna(arr).any():
                raise ValueError("Cannot convert NaT values to integer")
            return arr.view(dtype)

        # allow frequency conversions
        if dtype.kind == "M":
            return arr.astype(dtype)

        raise TypeError(f"cannot astype a datetimelike from [{arr.dtype}] to [{dtype}]")

    elif is_timedelta64_dtype(arr):
        if is_object_dtype(dtype):
            return tslibs.ints_to_pytimedelta(arr.view(np.int64))
        elif dtype == np.int64:
            if isna(arr).any():
                raise ValueError("Cannot convert NaT values to integer")
            return arr.view(dtype)

        if dtype not in [_INT64_DTYPE, _TD_DTYPE]:

            # allow frequency conversions
            # we return a float here!
            if dtype.kind == "m":
                mask = isna(arr)
                result = arr.astype(dtype).astype(np.float64)
                result[mask] = np.nan
                return result
        elif dtype == _TD_DTYPE:
            return arr.astype(_TD_DTYPE, copy=copy)

        raise TypeError(f"cannot astype a timedelta from [{arr.dtype}] to [{dtype}]")

    elif np.issubdtype(arr.dtype, np.floating) and np.issubdtype(dtype, np.integer):

        if not np.isfinite(arr).all():
            raise ValueError("Cannot convert non-finite values (NA or inf) to integer")

    elif is_object_dtype(arr):

        # work around NumPy brokenness, #1987
        if np.issubdtype(dtype.type, np.integer):
            return lib.astype_intsafe(arr.ravel(), dtype).reshape(arr.shape)

        # if we have a datetime/timedelta array of objects
        # then coerce to a proper dtype and recall astype_nansafe

        elif is_datetime64_dtype(dtype):
            from pandas import to_datetime

            return astype_nansafe(to_datetime(arr).values, dtype, copy=copy)
        elif is_timedelta64_dtype(dtype):
            from pandas import to_timedelta

            return astype_nansafe(to_timedelta(arr).values, dtype, copy=copy)

    if dtype.name in ("datetime64", "timedelta64"):
        msg = (
            f"The '{dtype.name}' dtype has no unit. Please pass in "
            f"'{dtype.name}[ns]' instead."
        )
        raise ValueError(msg)

    if copy or is_object_dtype(arr) or is_object_dtype(dtype):
        # Explicit copy, or required since NumPy can't view from / to object.
        return arr.astype(dtype, copy=True)

    return arr.view(dtype)


def maybe_convert_objects(values: np.ndarray, convert_numeric: bool = True):
    """
    If we have an object dtype array, try to coerce dates and/or numbers.

    Parameters
    ----------
    values : ndarray
    convert_numeric : bool, default True

    Returns
    -------
    ndarray or DatetimeIndex
    """
    validate_bool_kwarg(convert_numeric, "convert_numeric")

    orig_values = values

    # convert dates
    if is_object_dtype(values.dtype):
        values = lib.maybe_convert_objects(values, convert_datetime=True)

    # convert timedeltas
    if is_object_dtype(values.dtype):
        values = lib.maybe_convert_objects(values, convert_timedelta=True)

    # convert to numeric
    if is_object_dtype(values.dtype):
        if convert_numeric:
            try:
                new_values = lib.maybe_convert_numeric(
                    values, set(), coerce_numeric=True
                )
            except (ValueError, TypeError):
                pass
            else:
                # if we are all nans then leave me alone
                if not isna(new_values).all():
                    values = new_values

        else:
            # soft-conversion
            values = lib.maybe_convert_objects(values)

    if values is orig_values:
        values = values.copy()

    return values


def soft_convert_objects(
    values: np.ndarray,
    datetime: bool = True,
    numeric: bool = True,
    timedelta: bool = True,
    coerce: bool = False,
    copy: bool = True,
):
    """ if we have an object dtype, try to coerce dates and/or numbers """

    validate_bool_kwarg(datetime, "datetime")
    validate_bool_kwarg(numeric, "numeric")
    validate_bool_kwarg(timedelta, "timedelta")
    validate_bool_kwarg(coerce, "coerce")
    validate_bool_kwarg(copy, "copy")

    conversion_count = sum((datetime, numeric, timedelta))
    if conversion_count == 0:
        raise ValueError("At least one of datetime, numeric or timedelta must be True.")
    elif conversion_count > 1 and coerce:
        raise ValueError(
            "Only one of 'datetime', 'numeric' or "
            "'timedelta' can be True when when coerce=True."
        )

    if not is_object_dtype(values.dtype):
        # If not object, do not attempt conversion
        values = values.copy() if copy else values
        return values

    # If 1 flag is coerce, ensure 2 others are False
    if coerce:
        # Immediate return if coerce
        if datetime:
            from pandas import to_datetime

            return to_datetime(values, errors="coerce").to_numpy()
        elif timedelta:
            from pandas import to_timedelta

            return to_timedelta(values, errors="coerce").to_numpy()
        elif numeric:
            from pandas import to_numeric

            return to_numeric(values, errors="coerce")

    # Soft conversions
    if datetime:
        # GH 20380, when datetime is beyond year 2262, hence outside
        # bound of nanosecond-resolution 64-bit integers.
        try:
            values = lib.maybe_convert_objects(values, convert_datetime=True)
        except OutOfBoundsDatetime:
            pass

    if timedelta and is_object_dtype(values.dtype):
        # Object check to ensure only run if previous did not convert
        values = lib.maybe_convert_objects(values, convert_timedelta=True)

    if numeric and is_object_dtype(values.dtype):
        try:
            converted = lib.maybe_convert_numeric(values, set(), coerce_numeric=True)
        except (ValueError, TypeError):
            pass
        else:
            # If all NaNs, then do not-alter
            values = converted if not isna(converted).all() else values
            values = values.copy() if copy else values

    return values


def maybe_castable(arr) -> bool:
    # return False to force a non-fastpath

    # check datetime64[ns]/timedelta64[ns] are valid
    # otherwise try to coerce
    kind = arr.dtype.kind
    if kind == "M":
        return is_datetime64_ns_dtype(arr.dtype)
    elif kind == "m":
        return is_timedelta64_ns_dtype(arr.dtype)

    return arr.dtype.name not in _POSSIBLY_CAST_DTYPES


def maybe_infer_to_datetimelike(value, convert_dates: bool = False):
    """
    we might have a array (or single object) that is datetime like,
    and no dtype is passed don't change the value unless we find a
    datetime/timedelta set

    this is pretty strict in that a datetime/timedelta is REQUIRED
    in addition to possible nulls/string likes

    Parameters
    ----------
    value : np.array / Series / Index / list-like
    convert_dates : bool, default False
       if True try really hard to convert dates (such as datetime.date), other
       leave inferred dtype 'date' alone

    """

    # TODO: why not timedelta?
    if isinstance(
        value, (ABCDatetimeIndex, ABCPeriodIndex, ABCDatetimeArray, ABCPeriodArray)
    ):
        return value
    elif isinstance(value, ABCSeries):
        if isinstance(value._values, ABCDatetimeIndex):
            return value._values

    v = value

    if not is_list_like(v):
        v = [v]
    v = np.array(v, copy=False)

    # we only care about object dtypes
    if not is_object_dtype(v):
        return value

    shape = v.shape
    if not v.ndim == 1:
        v = v.ravel()

    if not len(v):
        return value

    def try_datetime(v):
        # safe coerce to datetime64
        try:
            # GH19671
            v = tslib.array_to_datetime(v, require_iso8601=True, errors="raise")[0]
        except ValueError:

            # we might have a sequence of the same-datetimes with tz's
            # if so coerce to a DatetimeIndex; if they are not the same,
            # then these stay as object dtype, xref GH19671
            from pandas._libs.tslibs import conversion
            from pandas import DatetimeIndex

            try:

                values, tz = conversion.datetime_to_datetime64(v)
                return DatetimeIndex(values).tz_localize("UTC").tz_convert(tz=tz)
            except (ValueError, TypeError):
                pass

        except Exception:
            pass

        return v.reshape(shape)

    def try_timedelta(v):
        # safe coerce to timedelta64

        # will try first with a string & object conversion
        from pandas import to_timedelta

        try:
            return to_timedelta(v)._ndarray_values.reshape(shape)
        except ValueError:
            return v.reshape(shape)

    inferred_type = lib.infer_datetimelike_array(ensure_object(v))

    if inferred_type == "date" and convert_dates:
        value = try_datetime(v)
    elif inferred_type == "datetime":
        value = try_datetime(v)
    elif inferred_type == "timedelta":
        value = try_timedelta(v)
    elif inferred_type == "nat":

        # if all NaT, return as datetime
        if isna(v).all():
            value = try_datetime(v)
        else:

            # We have at least a NaT and a string
            # try timedelta first to avoid spurious datetime conversions
            # e.g. '00:00:01' is a timedelta but technically is also a datetime
            value = try_timedelta(v)
            if lib.infer_dtype(value, skipna=False) in ["mixed"]:
                # cannot skip missing values, as NaT implies that the string
                # is actually a datetime
                value = try_datetime(v)

    return value


def maybe_cast_to_datetime(value, dtype, errors: str = "raise"):
    """ try to cast the array/value to a datetimelike dtype, converting float
    nan to iNaT
    """
    from pandas.core.tools.timedeltas import to_timedelta
    from pandas.core.tools.datetimes import to_datetime

    if dtype is not None:
        if isinstance(dtype, str):
            dtype = np.dtype(dtype)

        is_datetime64 = is_datetime64_dtype(dtype)
        is_datetime64tz = is_datetime64tz_dtype(dtype)
        is_timedelta64 = is_timedelta64_dtype(dtype)

        if is_datetime64 or is_datetime64tz or is_timedelta64:

            # Force the dtype if needed.
            msg = (
                f"The '{dtype.name}' dtype has no unit. "
                f"Please pass in '{dtype.name}[ns]' instead."
            )

            if is_datetime64 and not is_dtype_equal(dtype, _NS_DTYPE):

                # pandas supports dtype whose granularity is less than [ns]
                # e.g., [ps], [fs], [as]
                if dtype <= np.dtype("M8[ns]"):
                    if dtype.name == "datetime64":
                        raise ValueError(msg)
                    dtype = _NS_DTYPE
                else:
                    raise TypeError(f"cannot convert datetimelike to dtype [{dtype}]")
            elif is_datetime64tz:

                # our NaT doesn't support tz's
                # this will coerce to DatetimeIndex with
                # a matching dtype below
                if is_scalar(value) and isna(value):
                    value = [value]

            elif is_timedelta64 and not is_dtype_equal(dtype, _TD_DTYPE):

                # pandas supports dtype whose granularity is less than [ns]
                # e.g., [ps], [fs], [as]
                if dtype <= np.dtype("m8[ns]"):
                    if dtype.name == "timedelta64":
                        raise ValueError(msg)
                    dtype = _TD_DTYPE
                else:
                    raise TypeError(f"cannot convert timedeltalike to dtype [{dtype}]")

            if is_scalar(value):
                if value == iNaT or isna(value):
                    value = iNaT
            else:
                value = np.array(value, copy=False)

                # have a scalar array-like (e.g. NaT)
                if value.ndim == 0:
                    value = iNaT

                # we have an array of datetime or timedeltas & nulls
                elif np.prod(value.shape) or not is_dtype_equal(value.dtype, dtype):
                    try:
                        if is_datetime64:
                            value = to_datetime(value, errors=errors)
                            # GH 25843: Remove tz information since the dtype
                            # didn't specify one
                            if value.tz is not None:
                                value = value.tz_localize(None)
                            value = value._values
                        elif is_datetime64tz:
                            # The string check can be removed once issue #13712
                            # is solved. String data that is passed with a
                            # datetime64tz is assumed to be naive which should
                            # be localized to the timezone.
                            is_dt_string = is_string_dtype(value)
                            value = to_datetime(value, errors=errors).array
                            if is_dt_string:
                                # Strings here are naive, so directly localize
                                value = value.tz_localize(dtype.tz)
                            else:
                                # Numeric values are UTC at this point,
                                # so localize and convert
                                value = value.tz_localize("UTC").tz_convert(dtype.tz)
                        elif is_timedelta64:
                            value = to_timedelta(value, errors=errors)._values
                    except OutOfBoundsDatetime:
                        raise
                    except (AttributeError, ValueError, TypeError):
                        pass

        # coerce datetimelike to object
        elif is_datetime64_dtype(value) and not is_datetime64_dtype(dtype):
            if is_object_dtype(dtype):
                if value.dtype != _NS_DTYPE:
                    value = value.astype(_NS_DTYPE)
                ints = np.asarray(value).view("i8")
                return tslib.ints_to_pydatetime(ints)

            # we have a non-castable dtype that was passed
            raise TypeError(f"Cannot cast datetime64 to {dtype}")

    else:

        is_array = isinstance(value, np.ndarray)

        # catch a datetime/timedelta that is not of ns variety
        # and no coercion specified
        if is_array and value.dtype.kind in ["M", "m"]:
            dtype = value.dtype

            if dtype.kind == "M" and dtype != _NS_DTYPE:
                value = tslibs.conversion.ensure_datetime64ns(value)

            elif dtype.kind == "m" and dtype != _TD_DTYPE:
                value = to_timedelta(value)

        # only do this if we have an array and the dtype of the array is not
        # setup already we are not an integer/object, so don't bother with this
        # conversion
        elif not (
            is_array
            and not (
                issubclass(value.dtype.type, np.integer) or value.dtype == np.object_
            )
        ):
            value = maybe_infer_to_datetimelike(value)

    return value


def find_common_type(types):
    """
    Find a common data type among the given dtypes.

    Parameters
    ----------
    types : list of dtypes

    Returns
    -------
    pandas extension or numpy dtype

    See Also
    --------
    numpy.find_common_type

    """

    if len(types) == 0:
        raise ValueError("no types given")

    first = types[0]

    # workaround for find_common_type([np.dtype('datetime64[ns]')] * 2)
    # => object
    if all(is_dtype_equal(first, t) for t in types[1:]):
        return first

    if any(isinstance(t, ExtensionDtype) for t in types):
        return np.object

    # take lowest unit
    if all(is_datetime64_dtype(t) for t in types):
        return np.dtype("datetime64[ns]")
    if all(is_timedelta64_dtype(t) for t in types):
        return np.dtype("timedelta64[ns]")

    # don't mix bool / int or float or complex
    # this is different from numpy, which casts bool with float/int as int
    has_bools = any(is_bool_dtype(t) for t in types)
    if has_bools:
        for t in types:
            if is_integer_dtype(t) or is_float_dtype(t) or is_complex_dtype(t):
                return np.object

    return np.find_common_type(types, [])


def cast_scalar_to_array(shape, value, dtype=None):
    """
    Create np.ndarray of specified shape and dtype, filled with values.

    Parameters
    ----------
    shape : tuple
    value : scalar value
    dtype : np.dtype, optional
        dtype to coerce

    Returns
    -------
    ndarray of shape, filled with value, of specified / inferred dtype

    """

    if dtype is None:
        dtype, fill_value = infer_dtype_from_scalar(value)
    else:
        fill_value = value

    values = np.empty(shape, dtype=dtype)
    values.fill(fill_value)

    return values


def construct_1d_arraylike_from_scalar(value, length: int, dtype):
    """
    create a np.ndarray / pandas type of specified shape and dtype
    filled with values

    Parameters
    ----------
    value : scalar value
    length : int
    dtype : pandas_dtype / np.dtype

    Returns
    -------
    np.ndarray / pandas type of length, filled with value

    """
    if is_extension_array_dtype(dtype):
        cls = dtype.construct_array_type()
        subarr = cls._from_sequence([value] * length, dtype=dtype)

    else:
        if not isinstance(dtype, (np.dtype, type(np.dtype))):
            dtype = dtype.dtype

        if length and is_integer_dtype(dtype) and isna(value):
            # coerce if we have nan for an integer dtype
            dtype = np.dtype("float64")
        elif isinstance(dtype, np.dtype) and dtype.kind in ("U", "S"):
            # we need to coerce to object dtype to avoid
            # to allow numpy to take our string as a scalar value
            dtype = object
            if not isna(value):
                value = ensure_str(value)

        subarr = np.empty(length, dtype=dtype)
        subarr.fill(value)

    return subarr


def construct_1d_object_array_from_listlike(values):
    """
    Transform any list-like object in a 1-dimensional numpy array of object
    dtype.

    Parameters
    ----------
    values : any iterable which has a len()

    Raises
    ------
    TypeError
        * If `values` does not have a len()

    Returns
    -------
    1-dimensional numpy array of dtype object
    """
    # numpy will try to interpret nested lists as further dimensions, hence
    # making a 1D array that contains list-likes is a bit tricky:
    result = np.empty(len(values), dtype="object")
    result[:] = values
    return result


def construct_1d_ndarray_preserving_na(values, dtype=None, copy: bool = False):
    """
    Construct a new ndarray, coercing `values` to `dtype`, preserving NA.

    Parameters
    ----------
    values : Sequence
    dtype : numpy.dtype, optional
    copy : bool, default False
        Note that copies may still be made with ``copy=False`` if casting
        is required.

    Returns
    -------
    arr : ndarray[dtype]

    Examples
    --------
    >>> np.array([1.0, 2.0, None], dtype='str')
    array(['1.0', '2.0', 'None'], dtype='<U4')

    >>> construct_1d_ndarray_preserving_na([1.0, 2.0, None], dtype=np.dtype('str'))
    array(['1.0', '2.0', None], dtype=object)
    """
    subarr = np.array(values, dtype=dtype, copy=copy)

    if dtype is not None and dtype.kind in ("U", "S"):
        # GH-21083
        # We can't just return np.array(subarr, dtype='str') since
        # NumPy will convert the non-string objects into strings
        # Including NA values. Se we have to go
        # string -> object -> update NA, which requires an
        # additional pass over the data.
        na_values = isna(values)
        subarr2 = subarr.astype(object)
        subarr2[na_values] = np.asarray(values, dtype=object)[na_values]
        subarr = subarr2

    return subarr


def maybe_cast_to_integer_array(arr, dtype, copy: bool = False):
    """
    Takes any dtype and returns the casted version, raising for when data is
    incompatible with integer/unsigned integer dtypes.

    .. versionadded:: 0.24.0

    Parameters
    ----------
    arr : array-like
        The array to cast.
    dtype : str, np.dtype
        The integer dtype to cast the array to.
    copy: bool, default False
        Whether to make a copy of the array before returning.

    Returns
    -------
    int_arr : ndarray
        An array of integer or unsigned integer dtype

    Raises
    ------
    OverflowError : the dtype is incompatible with the data
    ValueError : loss of precision has occurred during casting

    Examples
    --------
    If you try to coerce negative values to unsigned integers, it raises:

    >>> Series([-1], dtype="uint64")
    Traceback (most recent call last):
        ...
    OverflowError: Trying to coerce negative values to unsigned integers

    Also, if you try to coerce float values to integers, it raises:

    >>> Series([1, 2, 3.5], dtype="int64")
    Traceback (most recent call last):
        ...
    ValueError: Trying to coerce float values to integers
    """

    try:
        if not hasattr(arr, "astype"):
            casted = np.array(arr, dtype=dtype, copy=copy)
        else:
            casted = arr.astype(dtype, copy=copy)
    except OverflowError:
        raise OverflowError(
            "The elements provided in the data cannot all be "
            f"casted to the dtype {dtype}"
        )

    if np.array_equal(arr, casted):
        return casted

    # We do this casting to allow for proper
    # data and dtype checking.
    #
    # We didn't do this earlier because NumPy
    # doesn't handle `uint64` correctly.
    arr = np.asarray(arr)

    if is_unsigned_integer_dtype(dtype) and (arr < 0).any():
        raise OverflowError("Trying to coerce negative values to unsigned integers")

    if is_integer_dtype(dtype) and (is_float_dtype(arr) or is_object_dtype(arr)):
        raise ValueError("Trying to coerce float values to integers")
