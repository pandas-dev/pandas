""" routings for casting """

from datetime import datetime, timedelta
import warnings

import numpy as np

from pandas._libs import lib, tslib, tslibs
from pandas._libs.tslibs import NaT, OutOfBoundsDatetime, Period, iNaT
from pandas.util._validators import validate_bool_kwarg

from .common import (
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
    is_categorical_dtype,
    is_complex,
    is_complex_dtype,
    is_datetime64_dtype,
    is_datetime64_ns_dtype,
    is_datetime64tz_dtype,
    is_datetime_or_timedelta_dtype,
    is_datetimelike,
    is_dtype_equal,
    is_extension_array_dtype,
    is_extension_type,
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
from .dtypes import DatetimeTZDtype, ExtensionDtype, PeriodDtype
from .generic import (
    ABCDataFrame,
    ABCDatetimeArray,
    ABCDatetimeIndex,
    ABCIndexClass,
    ABCPeriodArray,
    ABCPeriodIndex,
    ABCSeries,
)
from .inference import is_list_like
from .missing import isna, notna

_int8_max = np.iinfo(np.int8).max
_int16_max = np.iinfo(np.int16).max
_int32_max = np.iinfo(np.int32).max
_int64_max = np.iinfo(np.int64).max
_int64_min = np.iinfo(np.int64).min
_uint64_max = np.iinfo(np.uint64).max
_float32_max = np.finfo(np.float32).max


def _is_iNaT(x):
    """
    Helper function to circumvent numpy bug for timedeltas

    Specifically, comparing a scalar timedelta against another scalar value may
    raise a spurious DeprecationWarning, see numpy/numpy#10095
    """
    if not is_scalar(x):
        return False
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        result = x == iNaT
    return result


def maybe_convert_platform(values):
    """ try to do platform conversion, allow ndarray or list here """

    if isinstance(values, (list, tuple, range)):
        values = construct_1d_object_array_from_listlike(values)
    if getattr(values, "dtype", None) == np.object_:
        if hasattr(values, "_values"):
            values = values._values
        values = lib.maybe_convert_objects(values)

    return values


def is_nested_object(obj):
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

        if isna(arr).any() or not np.allclose(arr, trans(arr).astype(dtype), rtol=0):
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


def maybe_upcast_putmask(result, mask, other):
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
    other : ndarray or scalar
        The source array or value

    Returns
    -------
    result : ndarray
    changed : boolean
        Set to true if the result array was upcasted

    Examples
    --------
    >>> result, _ = maybe_upcast_putmask(np.arange(1,6),
    np.array([False, True, False, True, True]), np.arange(21,23))
    >>> result
    array([1, 21, 3, 22, 21])
    """

    if not isinstance(result, np.ndarray):
        raise ValueError("The result input must be a ndarray.")

    if mask.any():
        # Two conversions for date-like dtypes that can't be done automatically
        # in np.place:
        #   NaN -> NaT
        #   integer or integer array -> date-like array
        if is_datetimelike(result.dtype):
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
    Determine minimal dtype to hold fill_value, when starting from dtype

    Parameters
    ----------
    dtype : DType
        The dtype to start from.
    fill_value : scalar or np.ndarray / Series / Index
        The value that the output dtype needs to be able to hold.

    Returns
    -------
    dtype : DType
        The updated dtype.
    fill_value : scalar
        The type of this value depends on the type of the passed fill_value

        * If fill_value is a scalar, the method returns that scalar, but
          modified to fit the updated dtype. For example, a datetime fill_value
          will be returned as an integer (representing ns) for M8[ns], and
          values considered missing (see pd.isna) will be returned as the
          corresponding missing value marker for the updated dtype.
        * If fill_value is an ndarray/Series/Index, this method will always
          return the missing value marker for the updated dtype. This value
          will be None for dtypes that cannot hold missing values (integers,
          booleans, bytes).

    See Also
    --------
    _maybe_promote_with_scalar : underlying method for scalar case
    _maybe_promote_with_array : underlying method for array case
    """
    if is_scalar(fill_value) or isinstance(fill_value, tuple):
        return _maybe_promote_with_scalar(dtype, fill_value)
    elif isinstance(fill_value, (np.ndarray, ABCSeries, ABCIndexClass)):
        return _maybe_promote_with_array(dtype, fill_value)
    else:
        fill_type = type(fill_value).__name__
        raise ValueError(
            "fill_value must either be scalar, or a Series / "
            "Index / np.ndarray; received {}".format(fill_type)
        )


def _maybe_promote_with_scalar(dtype, fill_value=np.nan):
    """
    Determine minimal dtype to hold fill_value, when starting from dtype

    Parameters
    ----------
    dtype : DType
        The dtype to start from.
    fill_value : scalar or np.ndarray / Series / Index
        The value that the output dtype needs to be able to hold.

    Returns
    -------
    dtype : DType
        The updated dtype.
    fill_value : scalar
        The passed fill_value, potentially modified to fit the updated dtype.
        For example, assuming a datetime dtype, a datetime.datetime fill_value
        will be returned as an integer (representing ns) for M8[ns]. Similarly,
        values considered missing (see pd.isna) will be returned as the
        corresponding missing value marker for the updated dtype.

    See Also
    --------
    _maybe_promote_with_array : similar method for array case
        This method contains the actual promotion logic for both cases.

    Examples
    --------
    >>> maybe_promote(np.dtype('int'), fill_value=np.nan)
    (dtype('float64'), nan)
    >>> maybe_promote(np.dtype('float'), fill_value='abcd')
    (dtype('O'), 'abcd')

    For datetimes, timedeltas and datetimes with a timezone, the missing value
    marker is pandas._libs.tslibs.iNaT (== np.iinfo('int64').min):

    >>> maybe_promote(np.dtype('datetime64[ns]'), fill_value=np.nan)
    (dtype('<M8[ns]'), -9223372036854775808)

    The method will generally infer as conservatively as possible, also for
    subtypes of integers / floats / complex:

    >>> maybe_promote(np.dtype('uint8'), fill_value=np.iinfo('uint8').max + 1)
    (dtype('uint16'), 256)
    >>> maybe_promote(np.dtype('uint8'), fill_value=-1)
    (dtype('int16'), -1)
    """
    from pandas import Series

    if not (is_scalar(fill_value) or isinstance(fill_value, tuple)):
        raise ValueError(
            "fill_value must be a scalar, received " "{}".format(type(fill_value))
        )

    # unify handling of scalar and array values to simplify actual
    # promotion logic in _maybe_promote_with_array;
    if is_object_dtype(dtype) and fill_value is not None:
        # inserting into object does not cast (except for None -> np.nan)
        return np.dtype(object), fill_value

    # use Series to construct, since np.array cannot deal with pandas-internal
    # dtypes (e.g. DatetimeTZDtype); furthermore, we want to treat tuples as
    # scalar, but numpy casts those to a new dimension
    fill_array = Series([fill_value], dtype=object)
    dtype, na_value = _maybe_promote_with_array(dtype, fill_array)

    # _maybe_promote_with_array returns the na-marker for the new dtype;
    # _maybe_promote_with_scalar always casts fill_value to the new dtype
    if is_integer_dtype(dtype) and _is_iNaT(fill_value):
        # _maybe_promote_with_array considers iNaT a missing value, and since
        # int dtypes cannot hold missing values, that method returns None as
        # the na_value. For scalars, we need to keep it however, to ensure
        # correct operations for datetime/timedelta code.
        fill_value = iNaT
    elif fill_value is NaT and is_object_dtype(dtype):
        # the presence of pd.NaT forced upcasting to object, and therefore
        # fill_value does not get cast to na-marker of object (cf. below)
        pass
    elif isna(fill_value) or _is_iNaT(fill_value):
        # cast missing values (incl. iNaT) to correct missing value marker for
        # the updated dtype
        fill_value = na_value
    # otherwise casts fill_value (= only entry of fill_array) to new dtype
    elif is_datetime_or_timedelta_dtype(dtype):
        # for datetime/timedelta, we need to return the underlying ints
        fill_value = fill_array.astype(dtype)[0].value
    else:
        fill_value = fill_array.astype(dtype)[0]

    return dtype, fill_value


def _maybe_promote_with_array(dtype, fill_value=np.nan):
    """
    Determine minimal dtype to hold fill_value, when starting from dtype

    This will also return the default missing value for the resulting dtype, if
    necessary (e.g. for datetime / timedelta, the missing value will be `iNaT`)

    Parameters
    ----------
    dtype : DType
        The dtype to start from.
    fill_value : np.ndarray / Series / Index
        Array-like of values that the output dtype needs to be able to hold.

    Returns
    -------
    dtype : DType
        The updated dtype.
    na_value : scalar
        The missing value for the new dtype. Returns None or dtypes that
        cannot hold missing values (integers, booleans, bytes).

    See Also
    --------
    _maybe_promote_with_scalar : similar method for scalar case

    Examples
    --------
    >>> maybe_promote(np.dtype('int'), fill_value=np.array([None]))
    (dtype('float64'), nan)
    >>> maybe_promote(np.dtype('float'), fill_value=np.array(['abcd']))
    (dtype('O'), nan)

    For datetimes, timedeltas and datetimes with a timezone, the missing value
    marker is pandas._libs.tslibs.iNaT (== np.iinfo('int64').min):

    >>> maybe_promote(np.dtype('datetime64[ns]'), fill_value=np.array([None]))
    (dtype('<M8[ns]'), -9223372036854775808)

    String values do not get cast to datetime/timedelta automatically, but
    force an upcast to object (with corresponding missing value marker nan).

    >>> maybe_promote(np.dtype('datetime64[ns]'),
    ...               fill_value=np.array(['2018-01-01']))
    (dtype('O'), nan)

    The method will infer as conservatively as possible for integer types:

    >>> maybe_promote(np.dtype('uint8'),
    ...               fill_value=np.array([np.iinfo('uint8').max + 1]))
    (dtype('uint16'), None)
    >>> maybe_promote(np.dtype('uint8'), fill_value=np.array([-1]))
    (dtype('int16'), None)
    """

    if isinstance(fill_value, np.ndarray):
        if fill_value.ndim == 0:
            # zero-dimensional arrays cannot be iterated over
            fill_value = np.expand_dims(fill_value, 0)
        elif fill_value.ndim > 1:
            # ndarray, but too high-dimensional
            fill_value = fill_value.ravel()
    elif not isinstance(fill_value, (ABCSeries, ABCIndexClass)):
        fill_type = type(fill_value).__name__
        raise ValueError(
            "fill_value must either be a Series / Index / "
            "np.ndarray, received {}".format(fill_type)
        )

    if all(isna(x) or _is_iNaT(x) for x in fill_value):
        # only missing values (or no values at all)

        if is_datetime_or_timedelta_dtype(dtype):
            return dtype, iNaT
        elif is_datetime64tz_dtype(dtype):
            # DatetimeTZDtype does not use iNaT as missing value marker
            return dtype, NaT

        na_value = np.nan
        if len(fill_value) == 0:
            # empty array; no values to force change
            if is_integer_dtype(dtype) or dtype in (bool, bytes):
                # these types do not have a missing value marker
                na_value = None
            # otherwise nothing changes
        elif any(x is NaT for x in fill_value):
            # presence of pd.NaT upcasts everything that's not
            # datetime/timedelta (see above) to object
            dtype = np.dtype(object)
        elif (
            is_integer_dtype(dtype)
            and dtype == "uint64"
            and all(x == iNaT for x in fill_value)
        ):
            # uint64 + negative int casts to object
            dtype = np.dtype(object)
        elif is_integer_dtype(dtype) and all(x == iNaT for x in fill_value):
            # integer + iNaT casts to int64
            dtype = np.dtype("int64")
            na_value = None
        elif is_integer_dtype(dtype):
            # integer + other missing value (np.nan / None) casts to float
            dtype = np.dtype("float64")
        elif is_extension_array_dtype(dtype):
            na_value = dtype.na_value
        elif is_string_dtype(dtype) or dtype in (bool, bytes):
            # original dtype cannot hold nans
            dtype = np.dtype(object)

        return dtype, na_value

    fill_dtype = fill_value.dtype
    if fill_dtype == object:
        # for object dtype, we determine if we actually need to upcast
        # by inferring the dtype of fill_value
        inferred_dtype = lib.infer_dtype(fill_value, skipna=True)

        # cases that would yield 'empty' have been treated in branch above
        if inferred_dtype in ["period", "interval", "datetime64tz"]:
            # TODO: handle & test pandas-dtypes
            # TODO: lib.infer_dtype does not support datetime64tz yet
            pass
        else:
            # rest can be mapped to numpy dtypes
            map_inferred_to_numpy = {
                "floating": float,
                "mixed-integer-float": float,
                "decimal": float,
                "integer": int,
                "boolean": bool,
                "complex": complex,
                "bytes": bytes,
                "datetime64": "datetime64[ns]",
                "datetime": "datetime64[ns]",
                "date": "datetime64[ns]",
                "timedelta64": "timedelta64[ns]",
                "timedelta": "timedelta64[ns]",
                "time": object,  # time cannot be cast to datetime/timedelta
                "string": object,
                "mixed-integer": object,
                "mixed": object,
            }
            fill_dtype = np.dtype(map_inferred_to_numpy[inferred_dtype])

    # now that we have the correct dtype; check how we must upcast
    # * extension arrays
    # * int vs int
    # * int vs float / complex
    # * float vs float
    # * float vs complex (and vice versa)
    # * bool
    # * bytes
    # * datetimetz
    # * datetime
    # * timedelta
    # * string/object

    # if (is_extension_array_dtype(dtype)
    #         or is_extension_array_dtype(fill_dtype)):
    #     # TODO: dispatch to ExtensionDType.maybe_promote? GH 24246
    if is_integer_dtype(dtype) and is_integer_dtype(fill_dtype):
        if is_unsigned_integer_dtype(dtype) and all(fill_value >= 0):
            # can stay unsigned
            fill_max = fill_value.max()
            if fill_max > _uint64_max:
                return np.dtype(object), np.nan

            while fill_max > np.iinfo(dtype).max:
                # itemsize is the number of bytes; times eight is number of
                # bits, which is used in the string identifier of the dtype;
                # if fill_max is above the max for that dtype,
                # we double the number of bytes/bits.
                dtype = np.dtype("uint{}".format(dtype.itemsize * 8 * 2))
            return dtype, None
        else:
            # cannot stay unsigned
            if dtype == "uint64":
                # need to hold negative values, but int64 cannot hold
                # maximum of uint64 -> needs object
                return np.dtype(object), np.nan
            elif is_unsigned_integer_dtype(dtype):
                # need to turn into signed integers to hold negative values
                # int8 cannot hold maximum of uint8; similar for 16/32
                # therefore, upcast at least to next higher int-type
                dtype = np.dtype("int{}".format(dtype.itemsize * 8 * 2))

            fill_max = fill_value.max()
            fill_min = fill_value.min()
            if isinstance(fill_max, np.uint64):
                # numpy comparator is broken for uint64;
                # see https://github.com/numpy/numpy/issues/12525
                # use .item to get int object
                fill_max = fill_max.item()

            # comparison mechanics are broken above _int64_max;
            # use greater equal instead of equal
            if fill_max >= _int64_max + 1 or fill_min <= _int64_min - 1:
                return np.dtype(object), np.nan

            while fill_max > np.iinfo(dtype).max or fill_min < np.iinfo(dtype).min:
                # same mechanism as above, but for int instead of uint
                dtype = np.dtype("int{}".format(dtype.itemsize * 8 * 2))
            return dtype, None
    elif is_integer_dtype(dtype) and is_float_dtype(fill_dtype):
        # int with float: always upcasts to float64
        return np.dtype("float64"), np.nan
    elif is_integer_dtype(dtype) and is_complex_dtype(fill_dtype):
        # int with complex: always upcasts to complex128
        return np.dtype("complex128"), np.nan
    elif (is_float_dtype(dtype) or is_complex_dtype(dtype)) and is_integer_dtype(
        fill_dtype
    ):
        # float/complex with int: always stays original float/complex dtype
        return dtype, np.nan
    elif is_float_dtype(dtype) and is_float_dtype(fill_dtype):
        # float with float; upcasts depending on absolute max of fill_value
        if dtype == "float32" and np.abs(fill_value).max() <= _float32_max:
            return dtype, np.nan
        # all other cases return float64
        return np.dtype("float64"), np.nan
    elif (is_float_dtype(dtype) or is_complex_dtype(dtype)) and (
        is_float_dtype(fill_dtype) or is_complex_dtype(fill_dtype)
    ):
        # at least one is complex; otherwise we'd have hit float/float above
        with warnings.catch_warnings():
            # work around GH 27610
            warnings.filterwarnings("ignore", category=FutureWarning)
            if (
                dtype in ["float32", "complex64"]
                and max(
                    np.abs(np.real(fill_value)).max(),  # also works for float
                    np.abs(np.imag(fill_value)).max(),
                )
                <= _float32_max
            ):
                return np.complex64, np.nan
            # all other cases return complex128
            return np.dtype("complex128"), np.nan
    elif is_bool_dtype(dtype) and is_bool_dtype(fill_dtype):
        # bool with bool is the only combination that stays bool; any other
        # combination involving bool upcasts to object, see else-clause below
        return dtype, None
    elif issubclass(dtype.type, np.bytes_) and issubclass(fill_dtype.type, np.bytes_):
        # bytes with bytes is the only combination that stays bytes; any other
        # combination involving bytes upcasts to object, see else-clause below
        return dtype, None
    elif (
        is_datetime64tz_dtype(dtype)
        and is_datetime64tz_dtype(fill_dtype)
        and (dtype.tz == fill_dtype.tz)
    ):
        # datetimetz with datetimetz with the same timezone is the only
        # combination that stays datetimetz (in particular, mixing timezones or
        # tz-aware and tz-naive datetimes will cast to object);  any other
        # combination involving datetimetz upcasts to object, see below
        return dtype, iNaT
    elif (is_timedelta64_dtype(dtype) and is_timedelta64_dtype(fill_dtype)) or (
        is_datetime64_dtype(dtype) and is_datetime64_dtype(fill_dtype)
    ):
        # datetime and timedelta try to cast; if successful, keep dtype,
        # otherwise upcast to object
        try:
            with warnings.catch_warnings():
                msg = (
                    "parsing timezone aware datetimes is deprecated; "
                    "this will raise an error in the future"
                )
                warnings.filterwarnings(
                    "ignore", message=msg, category=DeprecationWarning
                )
                fill_value.astype(dtype)
            na_value = iNaT
        except (ValueError, TypeError):
            dtype = np.dtype(object)
            na_value = np.nan
        return dtype, na_value
    else:
        # anything else (e.g. strings, objects, or unmatched
        # bool / bytes / datetime / datetimetz / timedelta)
        return np.dtype(object), np.nan


def infer_dtype_from(val, pandas_dtype=False):
    """
    interpret the dtype from a scalar or array. This is a convenience
    routines to infer dtype from a scalar or an array

    Parameters
    ----------
    pandas_dtype : bool, default False
        whether to infer dtype including pandas extension types.
        If False, scalar/array belongs to pandas extension types is inferred as
        object
    """
    if is_scalar(val):
        return infer_dtype_from_scalar(val, pandas_dtype=pandas_dtype)
    return infer_dtype_from_array(val, pandas_dtype=pandas_dtype)


def infer_dtype_from_scalar(val, pandas_dtype=False):
    """
    interpret the dtype from a scalar

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

    return dtype, val


def infer_dtype_from_array(arr, pandas_dtype=False):
    """
    infer the dtype from a scalar or array

    Parameters
    ----------
    arr : scalar or array
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

    if pandas_dtype and is_extension_type(arr):
        return arr.dtype, arr

    elif isinstance(arr, ABCSeries):
        return arr.dtype, np.asarray(arr)

    # don't force numpy coerce with nan's
    inferred = lib.infer_dtype(arr, skipna=False)
    if inferred in ["string", "bytes", "unicode", "mixed", "mixed-integer"]:
        return (np.object_, arr)

    arr = np.asarray(arr)
    return arr.dtype, arr


def maybe_infer_dtype_type(element):
    """Try to infer an object's dtype, for use in arithmetic ops

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


def maybe_upcast(values, fill_value=np.nan, dtype=None, copy=False):
    """ provide explicit type promotion and coercion

    Parameters
    ----------
    values : the ndarray that we want to maybe upcast
    fill_value : what we want to fill with
    dtype : if None, then use the dtype of the values, else coerce to this type
    copy : if True always make a copy even if no upcast is required
    """

    if is_extension_type(values):
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


def astype_nansafe(arr, dtype, copy=True, skipna=False):
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
            return arr.view(dtype)

        # allow frequency conversions
        if dtype.kind == "M":
            return arr.astype(dtype)

        raise TypeError(
            "cannot astype a datetimelike from [{from_dtype}] "
            "to [{to_dtype}]".format(from_dtype=arr.dtype, to_dtype=dtype)
        )

    elif is_timedelta64_dtype(arr):
        if is_object_dtype(dtype):
            return tslibs.ints_to_pytimedelta(arr.view(np.int64))
        elif dtype == np.int64:
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

        raise TypeError(
            "cannot astype a timedelta from [{from_dtype}] "
            "to [{to_dtype}]".format(from_dtype=arr.dtype, to_dtype=dtype)
        )

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
        msg = "The '{dtype}' dtype has no unit. Please pass in '{dtype}[ns]' instead."
        raise ValueError(msg.format(dtype=dtype.name))

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


def maybe_castable(arr):
    # return False to force a non-fastpath

    # check datetime64[ns]/timedelta64[ns] are valid
    # otherwise try to coerce
    kind = arr.dtype.kind
    if kind == "M":
        return is_datetime64_ns_dtype(arr.dtype)
    elif kind == "m":
        return is_timedelta64_ns_dtype(arr.dtype)

    return arr.dtype.name not in _POSSIBLY_CAST_DTYPES


def maybe_infer_to_datetimelike(value, convert_dates=False):
    """
    we might have a array (or single object) that is datetime like,
    and no dtype is passed don't change the value unless we find a
    datetime/timedelta set

    this is pretty strict in that a datetime/timedelta is REQUIRED
    in addition to possible nulls/string likes

    Parameters
    ----------
    value : np.array / Series / Index / list-like
    convert_dates : boolean, default False
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


def maybe_cast_to_datetime(value, dtype, errors="raise"):
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
                "The '{dtype}' dtype has no unit. "
                "Please pass in '{dtype}[ns]' instead."
            )

            if is_datetime64 and not is_dtype_equal(dtype, _NS_DTYPE):

                # pandas supports dtype whose granularity is less than [ns]
                # e.g., [ps], [fs], [as]
                if dtype <= np.dtype("M8[ns]"):
                    if dtype.name == "datetime64":
                        raise ValueError(msg.format(dtype=dtype.name))
                    dtype = _NS_DTYPE
                else:
                    raise TypeError(
                        "cannot convert datetimelike to "
                        "dtype [{dtype}]".format(dtype=dtype)
                    )
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
                        raise ValueError(msg.format(dtype=dtype.name))
                    dtype = _TD_DTYPE
                else:
                    raise TypeError(
                        "cannot convert timedeltalike to "
                        "dtype [{dtype}]".format(dtype=dtype)
                    )

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
            raise TypeError("Cannot cast datetime64 to {dtype}".format(dtype=dtype))

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
    create np.ndarray of specified shape and dtype, filled with values

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


def construct_1d_arraylike_from_scalar(value, length, dtype):
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
    if is_datetime64tz_dtype(dtype):
        from pandas import DatetimeIndex

        subarr = DatetimeIndex([value] * length, dtype=dtype)
    elif is_categorical_dtype(dtype):
        from pandas import Categorical

        subarr = Categorical([value] * length, dtype=dtype)
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


def construct_1d_ndarray_preserving_na(values, dtype=None, copy=False):
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

    >>> construct_1d_ndarray_preserving_na([1.0, 2.0, None], dtype='str')


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


def maybe_cast_to_integer_array(arr, dtype, copy=False):
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
    copy: boolean, default False
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
            "casted to the dtype {dtype}".format(dtype=dtype)
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
