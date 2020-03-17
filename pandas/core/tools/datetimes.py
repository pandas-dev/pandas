from collections import abc
from datetime import datetime, time
from functools import partial
from itertools import islice
from typing import List, Optional, TypeVar, Union

import numpy as np

from pandas._libs import tslib, tslibs
from pandas._libs.tslibs import Timestamp, conversion, parsing
from pandas._libs.tslibs.parsing import (  # noqa
    DateParseError,
    _format_is_iso,
    _guess_datetime_format,
)
from pandas._libs.tslibs.strptime import array_strptime
from pandas._typing import ArrayLike

from pandas.core.dtypes.common import (
    ensure_object,
    is_datetime64_dtype,
    is_datetime64_ns_dtype,
    is_datetime64tz_dtype,
    is_float,
    is_integer,
    is_integer_dtype,
    is_list_like,
    is_numeric_dtype,
    is_scalar,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeIndex,
    ABCIndex,
    ABCIndexClass,
    ABCSeries,
)
from pandas.core.dtypes.missing import notna

from pandas.arrays import DatetimeArray, IntegerArray
from pandas.core import algorithms
from pandas.core.algorithms import unique
from pandas.core.arrays.datetimes import tz_to_dtype

# ---------------------------------------------------------------------
# types used in annotations

ArrayConvertible = Union[list, tuple, ArrayLike, ABCSeries]
Scalar = Union[int, float, str]
DatetimeScalar = TypeVar("DatetimeScalar", Scalar, datetime)
DatetimeScalarOrArrayConvertible = Union[
    DatetimeScalar, list, tuple, ArrayLike, ABCSeries
]


# ---------------------------------------------------------------------


def _guess_datetime_format_for_array(arr, **kwargs):
    # Try to guess the format based on the first non-NaN element
    non_nan_elements = notna(arr).nonzero()[0]
    if len(non_nan_elements):
        return _guess_datetime_format(arr[non_nan_elements[0]], **kwargs)


def should_cache(
    arg: ArrayConvertible, unique_share: float = 0.7, check_count: Optional[int] = None
) -> bool:
    """
    Decides whether to do caching.

    If the percent of unique elements among `check_count` elements less
    than `unique_share * 100` then we can do caching.

    Parameters
    ----------
    arg: listlike, tuple, 1-d array, Series
    unique_share: float, default=0.7, optional
        0 < unique_share < 1
    check_count: int, optional
        0 <= check_count <= len(arg)

    Returns
    -------
    do_caching: bool

    Notes
    -----
    By default for a sequence of less than 50 items in size, we don't do
    caching; for the number of elements less than 5000, we take ten percent of
    all elements to check for a uniqueness share; if the sequence size is more
    than 5000, then we check only the first 500 elements.
    All constants were chosen empirically by.
    """
    do_caching = True

    # default realization
    if check_count is None:
        # in this case, the gain from caching is negligible
        if len(arg) <= 50:
            return False

        if len(arg) <= 5000:
            check_count = int(len(arg) * 0.1)
        else:
            check_count = 500
    else:
        assert (
            0 <= check_count <= len(arg)
        ), "check_count must be in next bounds: [0; len(arg)]"
        if check_count == 0:
            return False

    assert 0 < unique_share < 1, "unique_share must be in next bounds: (0; 1)"

    unique_elements = set(islice(arg, check_count))
    if len(unique_elements) > check_count * unique_share:
        do_caching = False
    return do_caching


def _maybe_cache(arg, format, cache, convert_listlike):
    """
    Create a cache of unique dates from an array of dates

    Parameters
    ----------
    arg : listlike, tuple, 1-d array, Series
    format : string
        Strftime format to parse time
    cache : boolean
        True attempts to create a cache of converted values
    convert_listlike : function
        Conversion function to apply on dates

    Returns
    -------
    cache_array : Series
        Cache of converted, unique dates. Can be empty
    """
    from pandas import Series

    cache_array = Series(dtype=object)

    if cache:
        # Perform a quicker unique check
        if not should_cache(arg):
            return cache_array

        unique_dates = unique(arg)
        if len(unique_dates) < len(arg):
            cache_dates = convert_listlike(unique_dates, format)
            cache_array = Series(cache_dates, index=unique_dates)
    return cache_array


def _box_as_indexlike(
    dt_array: ArrayLike, utc: Optional[bool] = None, name: Optional[str] = None
) -> Union[ABCIndex, ABCDatetimeIndex]:
    """
    Properly boxes the ndarray of datetimes to DatetimeIndex
    if it is possible or to generic Index instead

    Parameters
    ----------
    dt_array: 1-d array
        Array of datetimes to be wrapped in an Index.
    tz : object
        None or 'utc'
    name : string, default None
        Name for a resulting index

    Returns
    -------
    result : datetime of converted dates
        - DatetimeIndex if convertible to sole datetime64 type
        - general Index otherwise
    """
    from pandas import DatetimeIndex, Index

    if is_datetime64_dtype(dt_array):
        tz = "utc" if utc else None
        return DatetimeIndex(dt_array, tz=tz, name=name)
    return Index(dt_array, name=name)


def _convert_and_box_cache(
    arg: DatetimeScalarOrArrayConvertible,
    cache_array: ABCSeries,
    name: Optional[str] = None,
) -> ABCIndexClass:
    """
    Convert array of dates with a cache and wrap the result in an Index.

    Parameters
    ----------
    arg : integer, float, string, datetime, list, tuple, 1-d array, Series
    cache_array : Series
        Cache of converted, unique dates
    name : string, default None
        Name for a DatetimeIndex

    Returns
    -------
    result : Index-like of converted dates
    """
    from pandas import Series

    result = Series(arg).map(cache_array)
    return _box_as_indexlike(result, utc=None, name=name)


def _return_parsed_timezone_results(result, timezones, tz, name):
    """
    Return results from array_strptime if a %z or %Z directive was passed.

    Parameters
    ----------
    result : ndarray
        int64 date representations of the dates
    timezones : ndarray
        pytz timezone objects
    tz : object
        None or pytz timezone object
    name : string, default None
        Name for a DatetimeIndex

    Returns
    -------
    tz_result : Index-like of parsed dates with timezone
    """
    if tz is not None:
        raise ValueError(
            "Cannot pass a tz argument when parsing strings with timezone information."
        )
    tz_results = np.array(
        [Timestamp(res).tz_localize(zone) for res, zone in zip(result, timezones)]
    )
    from pandas import Index

    return Index(tz_results, name=name)


def _convert_listlike_datetimes(
    arg,
    format,
    name=None,
    tz=None,
    unit=None,
    errors=None,
    infer_datetime_format=None,
    dayfirst=None,
    yearfirst=None,
    exact=None,
):
    """
    Helper function for to_datetime. Performs the conversions of 1D listlike
    of dates

    Parameters
    ----------
    arg : list, tuple, ndarray, Series, Index
        date to be parced
    name : object
        None or string for the Index name
    tz : object
        None or 'utc'
    unit : string
        None or string of the frequency of the passed data
    errors : string
        error handing behaviors from to_datetime, 'raise', 'coerce', 'ignore'
    infer_datetime_format : boolean
        inferring format behavior from to_datetime
    dayfirst : boolean
        dayfirst parsing behavior from to_datetime
    yearfirst : boolean
        yearfirst parsing behavior from to_datetime
    exact : boolean
        exact format matching behavior from to_datetime

    Returns
    -------
    Index-like of parsed dates
    """
    from pandas import DatetimeIndex
    from pandas.core.arrays.datetimes import (
        maybe_convert_dtype,
        objects_to_datetime64ns,
    )

    if isinstance(arg, (list, tuple)):
        arg = np.array(arg, dtype="O")

    # these are shortcutable
    if is_datetime64tz_dtype(arg):
        if not isinstance(arg, (DatetimeArray, DatetimeIndex)):
            return DatetimeIndex(arg, tz=tz, name=name)
        if tz == "utc":
            # error: Item "DatetimeIndex" of "Union[DatetimeArray, DatetimeIndex]" has
            # no attribute "tz_convert"
            arg = arg.tz_convert(None).tz_localize(tz)  # type: ignore
        return arg

    elif is_datetime64_ns_dtype(arg):
        if not isinstance(arg, (DatetimeArray, DatetimeIndex)):
            try:
                return DatetimeIndex(arg, tz=tz, name=name)
            except ValueError:
                pass
        elif tz:
            # DatetimeArray, DatetimeIndex
            # error: Item "DatetimeIndex" of "Union[DatetimeArray, DatetimeIndex]" has
            # no attribute "tz_localize"
            return arg.tz_localize(tz)  # type: ignore

        return arg

    elif unit is not None:
        if format is not None:
            raise ValueError("cannot specify both format and unit")
        arg = getattr(arg, "_values", arg)

        # GH 30050 pass an ndarray to tslib.array_with_unit_to_datetime
        # because it expects an ndarray argument
        if isinstance(arg, IntegerArray):
            result = arg.astype(f"datetime64[{unit}]")
            tz_parsed = None
        else:

            result, tz_parsed = tslib.array_with_unit_to_datetime(
                arg, unit, errors=errors
            )

        if errors == "ignore":
            from pandas import Index

            result = Index(result, name=name)
        else:
            result = DatetimeIndex(result, name=name)
        # GH 23758: We may still need to localize the result with tz
        # GH 25546: Apply tz_parsed first (from arg), then tz (from caller)
        # result will be naive but in UTC
        try:
            result = result.tz_localize("UTC").tz_convert(tz_parsed)
        except AttributeError:
            # Regular Index from 'ignore' path
            return result
        if tz is not None:
            if result.tz is None:
                result = result.tz_localize(tz)
            else:
                result = result.tz_convert(tz)
        return result
    elif getattr(arg, "ndim", 1) > 1:
        raise TypeError(
            "arg must be a string, datetime, list, tuple, 1-d array, or Series"
        )

    # warn if passing timedelta64, raise for PeriodDtype
    # NB: this must come after unit transformation
    orig_arg = arg
    try:
        arg, _ = maybe_convert_dtype(arg, copy=False)
    except TypeError:
        if errors == "coerce":
            result = np.array(["NaT"], dtype="datetime64[ns]").repeat(len(arg))
            return DatetimeIndex(result, name=name)
        elif errors == "ignore":
            from pandas import Index

            result = Index(arg, name=name)
            return result
        raise

    arg = ensure_object(arg)
    require_iso8601 = False

    if infer_datetime_format and format is None:
        format = _guess_datetime_format_for_array(arg, dayfirst=dayfirst)

    if format is not None:
        # There is a special fast-path for iso8601 formatted
        # datetime strings, so in those cases don't use the inferred
        # format because this path makes process slower in this
        # special case
        format_is_iso8601 = _format_is_iso(format)
        if format_is_iso8601:
            require_iso8601 = not infer_datetime_format
            format = None

    tz_parsed = None
    result = None

    if format is not None:
        try:
            # shortcut formatting here
            if format == "%Y%m%d":
                try:
                    # pass orig_arg as float-dtype may have been converted to
                    # datetime64[ns]
                    orig_arg = ensure_object(orig_arg)
                    result = _attempt_YYYYMMDD(orig_arg, errors=errors)
                except (ValueError, TypeError, tslibs.OutOfBoundsDatetime) as err:
                    raise ValueError(
                        "cannot convert the input to '%Y%m%d' date format"
                    ) from err

            # fallback
            if result is None:
                try:
                    result, timezones = array_strptime(
                        arg, format, exact=exact, errors=errors
                    )
                    if "%Z" in format or "%z" in format:
                        return _return_parsed_timezone_results(
                            result, timezones, tz, name
                        )
                except tslibs.OutOfBoundsDatetime:
                    if errors == "raise":
                        raise
                    elif errors == "coerce":
                        result = np.empty(arg.shape, dtype="M8[ns]")
                        iresult = result.view("i8")
                        iresult.fill(tslibs.iNaT)
                    else:
                        result = arg
                except ValueError:
                    # if format was inferred, try falling back
                    # to array_to_datetime - terminate here
                    # for specified formats
                    if not infer_datetime_format:
                        if errors == "raise":
                            raise
                        elif errors == "coerce":
                            result = np.empty(arg.shape, dtype="M8[ns]")
                            iresult = result.view("i8")
                            iresult.fill(tslibs.iNaT)
                        else:
                            result = arg
        except ValueError as e:
            # Fallback to try to convert datetime objects if timezone-aware
            #  datetime objects are found without passing `utc=True`
            try:
                values, tz = conversion.datetime_to_datetime64(arg)
                dta = DatetimeArray(values, dtype=tz_to_dtype(tz))
                return DatetimeIndex._simple_new(dta, name=name)
            except (ValueError, TypeError):
                raise e

    if result is None:
        assert format is None or infer_datetime_format
        utc = tz == "utc"
        result, tz_parsed = objects_to_datetime64ns(
            arg,
            dayfirst=dayfirst,
            yearfirst=yearfirst,
            utc=utc,
            errors=errors,
            require_iso8601=require_iso8601,
            allow_object=True,
        )

    if tz_parsed is not None:
        # We can take a shortcut since the datetime64 numpy array
        # is in UTC
        dta = DatetimeArray(result, dtype=tz_to_dtype(tz_parsed))
        return DatetimeIndex._simple_new(dta, name=name)

    utc = tz == "utc"
    return _box_as_indexlike(result, utc=utc, name=name)


def _adjust_to_origin(arg, origin, unit):
    """
    Helper function for to_datetime.
    Adjust input argument to the specified origin

    Parameters
    ----------
    arg : list, tuple, ndarray, Series, Index
        date to be adjusted
    origin : 'julian' or Timestamp
        origin offset for the arg
    unit : string
        passed unit from to_datetime, must be 'D'

    Returns
    -------
    ndarray or scalar of adjusted date(s)
    """
    if origin == "julian":
        original = arg
        j0 = Timestamp(0).to_julian_date()
        if unit != "D":
            raise ValueError("unit must be 'D' for origin='julian'")
        try:
            arg = arg - j0
        except TypeError as err:
            raise ValueError(
                "incompatible 'arg' type for given 'origin'='julian'"
            ) from err

        # preemptively check this for a nice range
        j_max = Timestamp.max.to_julian_date() - j0
        j_min = Timestamp.min.to_julian_date() - j0
        if np.any(arg > j_max) or np.any(arg < j_min):
            raise tslibs.OutOfBoundsDatetime(
                f"{original} is Out of Bounds for origin='julian'"
            )
    else:
        # arg must be numeric
        if not (
            (is_scalar(arg) and (is_integer(arg) or is_float(arg)))
            or is_numeric_dtype(np.asarray(arg))
        ):
            raise ValueError(
                f"'{arg}' is not compatible with origin='{origin}'; "
                "it must be numeric with a unit specified"
            )

        # we are going to offset back to unix / epoch time
        try:
            offset = Timestamp(origin)
        except tslibs.OutOfBoundsDatetime as err:
            raise tslibs.OutOfBoundsDatetime(
                f"origin {origin} is Out of Bounds"
            ) from err
        except ValueError as err:
            raise ValueError(
                f"origin {origin} cannot be converted to a Timestamp"
            ) from err

        if offset.tz is not None:
            raise ValueError(f"origin offset {offset} must be tz-naive")
        offset -= Timestamp(0)

        # convert the offset to the unit of the arg
        # this should be lossless in terms of precision
        offset = offset // tslibs.Timedelta(1, unit=unit)

        # scalars & ndarray-like can handle the addition
        if is_list_like(arg) and not isinstance(
            arg, (ABCSeries, ABCIndexClass, np.ndarray)
        ):
            arg = np.asarray(arg)
        arg = arg + offset
    return arg


def to_datetime(
    arg,
    errors="raise",
    dayfirst=False,
    yearfirst=False,
    utc=None,
    format=None,
    exact=True,
    unit=None,
    infer_datetime_format=False,
    origin="unix",
    cache=True,
):
    """
    Convert argument to datetime.

    Parameters
    ----------
    arg : int, float, str, datetime, list, tuple, 1-d array, Series, DataFrame/dict-like
        The object to convert to a datetime.
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception.
        - If 'coerce', then invalid parsing will be set as NaT.
        - If 'ignore', then invalid parsing will return the input.
    dayfirst : bool, default False
        Specify a date parse order if `arg` is str or its list-likes.
        If True, parses dates with the day first, eg 10/11/12 is parsed as
        2012-11-10.
        Warning: dayfirst=True is not strict, but will prefer to parse
        with day first (this is a known bug, based on dateutil behavior).
    yearfirst : bool, default False
        Specify a date parse order if `arg` is str or its list-likes.

        - If True parses dates with the year first, eg 10/11/12 is parsed as
          2010-11-12.
        - If both dayfirst and yearfirst are True, yearfirst is preceded (same
          as dateutil).

        Warning: yearfirst=True is not strict, but will prefer to parse
        with year first (this is a known bug, based on dateutil behavior).
    utc : bool, default None
        Return UTC DatetimeIndex if True (converting any tz-aware
        datetime.datetime objects as well).
    format : str, default None
        The strftime to parse time, eg "%d/%m/%Y", note that "%f" will parse
        all the way up to nanoseconds.
        See strftime documentation for more information on choices:
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior.
    exact : bool, True by default
        Behaves as:
        - If True, require an exact format match.
        - If False, allow the format to match anywhere in the target string.

    unit : str, default 'ns'
        The unit of the arg (D,s,ms,us,ns) denote the unit, which is an
        integer or float number. This will be based off the origin.
        Example, with unit='ms' and origin='unix' (the default), this
        would calculate the number of milliseconds to the unix epoch start.
    infer_datetime_format : bool, default False
        If True and no `format` is given, attempt to infer the format of the
        datetime strings, and if it can be inferred, switch to a faster
        method of parsing them. In some cases this can increase the parsing
        speed by ~5-10x.
    origin : scalar, default 'unix'
        Define the reference date. The numeric values would be parsed as number
        of units (defined by `unit`) since this reference date.

        - If 'unix' (or POSIX) time; origin is set to 1970-01-01.
        - If 'julian', unit must be 'D', and origin is set to beginning of
          Julian Calendar. Julian day number 0 is assigned to the day starting
          at noon on January 1, 4713 BC.
        - If Timestamp convertible, origin is set to Timestamp identified by
          origin.
    cache : bool, default True
        If True, use a cache of unique, converted dates to apply the datetime
        conversion. May produce significant speed-up when parsing duplicate
        date strings, especially ones with timezone offsets. The cache is only
        used when there are at least 50 values. The presence of out-of-bounds
        values will render the cache unusable and may slow down parsing.

        .. versionadded:: 0.23.0

        .. versionchanged:: 0.25.0
            - changed default value from False to True.

    Returns
    -------
    datetime
        If parsing succeeded.
        Return type depends on input:

        - list-like: DatetimeIndex
        - Series: Series of datetime64 dtype
        - scalar: Timestamp

        In case when it is not possible to return designated types (e.g. when
        any element of input is before Timestamp.min or after Timestamp.max)
        return will have datetime.datetime type (or corresponding
        array/Series).

    See Also
    --------
    DataFrame.astype : Cast argument to a specified dtype.
    to_timedelta : Convert argument to timedelta.
    convert_dtypes : Convert dtypes.

    Examples
    --------
    Assembling a datetime from multiple columns of a DataFrame. The keys can be
    common abbreviations like ['year', 'month', 'day', 'minute', 'second',
    'ms', 'us', 'ns']) or plurals of the same

    >>> df = pd.DataFrame({'year': [2015, 2016],
    ...                    'month': [2, 3],
    ...                    'day': [4, 5]})
    >>> pd.to_datetime(df)
    0   2015-02-04
    1   2016-03-05
    dtype: datetime64[ns]

    If a date does not meet the `timestamp limitations
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html
    #timeseries-timestamp-limits>`_, passing errors='ignore'
    will return the original input instead of raising any exception.

    Passing errors='coerce' will force an out-of-bounds date to NaT,
    in addition to forcing non-dates (or non-parseable dates) to NaT.

    >>> pd.to_datetime('13000101', format='%Y%m%d', errors='ignore')
    datetime.datetime(1300, 1, 1, 0, 0)
    >>> pd.to_datetime('13000101', format='%Y%m%d', errors='coerce')
    NaT

    Passing infer_datetime_format=True can often-times speedup a parsing
    if its not an ISO8601 format exactly, but in a regular format.

    >>> s = pd.Series(['3/11/2000', '3/12/2000', '3/13/2000'] * 1000)
    >>> s.head()
    0    3/11/2000
    1    3/12/2000
    2    3/13/2000
    3    3/11/2000
    4    3/12/2000
    dtype: object

    >>> %timeit pd.to_datetime(s, infer_datetime_format=True)  # doctest: +SKIP
    100 loops, best of 3: 10.4 ms per loop

    >>> %timeit pd.to_datetime(s, infer_datetime_format=False)  # doctest: +SKIP
    1 loop, best of 3: 471 ms per loop

    Using a unix epoch time

    >>> pd.to_datetime(1490195805, unit='s')
    Timestamp('2017-03-22 15:16:45')
    >>> pd.to_datetime(1490195805433502912, unit='ns')
    Timestamp('2017-03-22 15:16:45.433502912')

    .. warning:: For float arg, precision rounding might happen. To prevent
        unexpected behavior use a fixed-width exact type.

    Using a non-unix epoch origin

    >>> pd.to_datetime([1, 2, 3], unit='D',
    ...                origin=pd.Timestamp('1960-01-01'))
    DatetimeIndex(['1960-01-02', '1960-01-03', '1960-01-04'], \
dtype='datetime64[ns]', freq=None)
    """
    if arg is None:
        return None

    if origin != "unix":
        arg = _adjust_to_origin(arg, origin, unit)

    tz = "utc" if utc else None
    convert_listlike = partial(
        _convert_listlike_datetimes,
        tz=tz,
        unit=unit,
        dayfirst=dayfirst,
        yearfirst=yearfirst,
        errors=errors,
        exact=exact,
        infer_datetime_format=infer_datetime_format,
    )

    if isinstance(arg, Timestamp):
        result = arg
        if tz is not None:
            if arg.tz is not None:
                result = result.tz_convert(tz)
            else:
                result = result.tz_localize(tz)
    elif isinstance(arg, ABCSeries):
        cache_array = _maybe_cache(arg, format, cache, convert_listlike)
        if not cache_array.empty:
            result = arg.map(cache_array)
        else:
            values = convert_listlike(arg._values, format)
            result = arg._constructor(values, index=arg.index, name=arg.name)
    elif isinstance(arg, (ABCDataFrame, abc.MutableMapping)):
        result = _assemble_from_unit_mappings(arg, errors, tz)
    elif isinstance(arg, ABCIndexClass):
        cache_array = _maybe_cache(arg, format, cache, convert_listlike)
        if not cache_array.empty:
            result = _convert_and_box_cache(arg, cache_array, name=arg.name)
        else:
            convert_listlike = partial(convert_listlike, name=arg.name)
            result = convert_listlike(arg, format)
    elif is_list_like(arg):
        try:
            cache_array = _maybe_cache(arg, format, cache, convert_listlike)
        except tslibs.OutOfBoundsDatetime:
            # caching attempts to create a DatetimeIndex, which may raise
            # an OOB. If that's the desired behavior, then just reraise...
            if errors == "raise":
                raise
            # ... otherwise, continue without the cache.
            from pandas import Series

            cache_array = Series([], dtype=object)  # just an empty array
        if not cache_array.empty:
            result = _convert_and_box_cache(arg, cache_array)
        else:
            result = convert_listlike(arg, format)
    else:
        result = convert_listlike(np.array([arg]), format)[0]

    return result


# mappings for assembling units
_unit_map = {
    "year": "year",
    "years": "year",
    "month": "month",
    "months": "month",
    "day": "day",
    "days": "day",
    "hour": "h",
    "hours": "h",
    "minute": "m",
    "minutes": "m",
    "second": "s",
    "seconds": "s",
    "ms": "ms",
    "millisecond": "ms",
    "milliseconds": "ms",
    "us": "us",
    "microsecond": "us",
    "microseconds": "us",
    "ns": "ns",
    "nanosecond": "ns",
    "nanoseconds": "ns",
}


def _assemble_from_unit_mappings(arg, errors, tz):
    """
    assemble the unit specified fields from the arg (DataFrame)
    Return a Series for actual parsing

    Parameters
    ----------
    arg : DataFrame
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'

        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaT
        - If 'ignore', then invalid parsing will return the input
    tz : None or 'utc'

    Returns
    -------
    Series
    """
    from pandas import to_timedelta, to_numeric, DataFrame

    arg = DataFrame(arg)
    if not arg.columns.is_unique:
        raise ValueError("cannot assemble with duplicate keys")

    # replace passed unit with _unit_map
    def f(value):
        if value in _unit_map:
            return _unit_map[value]

        # m is case significant
        if value.lower() in _unit_map:
            return _unit_map[value.lower()]

        return value

    unit = {k: f(k) for k in arg.keys()}
    unit_rev = {v: k for k, v in unit.items()}

    # we require at least Ymd
    required = ["year", "month", "day"]
    req = sorted(set(required) - set(unit_rev.keys()))
    if len(req):
        _required = ",".join(req)
        raise ValueError(
            "to assemble mappings requires at least that "
            f"[year, month, day] be specified: [{_required}] is missing"
        )

    # keys we don't recognize
    excess = sorted(set(unit_rev.keys()) - set(_unit_map.values()))
    if len(excess):
        _excess = ",".join(excess)
        raise ValueError(
            f"extra keys have been passed to the datetime assemblage: [{_excess}]"
        )

    def coerce(values):
        # we allow coercion to if errors allows
        values = to_numeric(values, errors=errors)

        # prevent overflow in case of int8 or int16
        if is_integer_dtype(values):
            values = values.astype("int64", copy=False)
        return values

    values = (
        coerce(arg[unit_rev["year"]]) * 10000
        + coerce(arg[unit_rev["month"]]) * 100
        + coerce(arg[unit_rev["day"]])
    )
    try:
        values = to_datetime(values, format="%Y%m%d", errors=errors, utc=tz)
    except (TypeError, ValueError) as err:
        raise ValueError(f"cannot assemble the datetimes: {err}") from err

    for u in ["h", "m", "s", "ms", "us", "ns"]:
        value = unit_rev.get(u)
        if value is not None and value in arg:
            try:
                values += to_timedelta(coerce(arg[value]), unit=u, errors=errors)
            except (TypeError, ValueError) as err:
                raise ValueError(
                    f"cannot assemble the datetimes [{value}]: {err}"
                ) from err
    return values


def _attempt_YYYYMMDD(arg, errors):
    """
    try to parse the YYYYMMDD/%Y%m%d format, try to deal with NaT-like,
    arg is a passed in as an object dtype, but could really be ints/strings
    with nan-like/or floats (e.g. with nan)

    Parameters
    ----------
    arg : passed value
    errors : 'raise','ignore','coerce'
    """

    def calc(carg):
        # calculate the actual result
        carg = carg.astype(object)
        parsed = parsing.try_parse_year_month_day(
            carg / 10000, carg / 100 % 100, carg % 100
        )
        return tslib.array_to_datetime(parsed, errors=errors)[0]

    def calc_with_mask(carg, mask):
        result = np.empty(carg.shape, dtype="M8[ns]")
        iresult = result.view("i8")
        iresult[~mask] = tslibs.iNaT

        masked_result = calc(carg[mask].astype(np.float64).astype(np.int64))
        result[mask] = masked_result.astype("M8[ns]")
        return result

    # try intlike / strings that are ints
    try:
        return calc(arg.astype(np.int64))
    except (ValueError, OverflowError, TypeError):
        pass

    # a float with actual np.nan
    try:
        carg = arg.astype(np.float64)
        return calc_with_mask(carg, notna(carg))
    except (ValueError, OverflowError, TypeError):
        pass

    # string with NaN-like
    try:
        mask = ~algorithms.isin(arg, list(tslib.nat_strings))
        return calc_with_mask(arg, mask)
    except (ValueError, OverflowError, TypeError):
        pass

    return None


# Fixed time formats for time parsing
_time_formats = [
    "%H:%M",
    "%H%M",
    "%I:%M%p",
    "%I%M%p",
    "%H:%M:%S",
    "%H%M%S",
    "%I:%M:%S%p",
    "%I%M%S%p",
]


def _guess_time_format_for_array(arr):
    # Try to guess the format based on the first non-NaN element
    non_nan_elements = notna(arr).nonzero()[0]
    if len(non_nan_elements):
        element = arr[non_nan_elements[0]]
        for time_format in _time_formats:
            try:
                datetime.strptime(element, time_format)
                return time_format
            except ValueError:
                pass

    return None


def to_time(arg, format=None, infer_time_format=False, errors="raise"):
    """
    Parse time strings to time objects using fixed strptime formats ("%H:%M",
    "%H%M", "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
    "%I%M%S%p")

    Use infer_time_format if all the strings are in the same format to speed
    up conversion.

    Parameters
    ----------
    arg : string in time format, datetime.time, list, tuple, 1-d array,  Series
    format : str, default None
        Format used to convert arg into a time object.  If None, fixed formats
        are used.
    infer_time_format: bool, default False
        Infer the time format based on the first non-NaN element.  If all
        strings are in the same format, this will speed up conversion.
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as None
        - If 'ignore', then invalid parsing will return the input

    Returns
    -------
    datetime.time
    """

    def _convert_listlike(arg, format):

        if isinstance(arg, (list, tuple)):
            arg = np.array(arg, dtype="O")

        elif getattr(arg, "ndim", 1) > 1:
            raise TypeError(
                "arg must be a string, datetime, list, tuple, 1-d array, or Series"
            )

        arg = ensure_object(arg)

        if infer_time_format and format is None:
            format = _guess_time_format_for_array(arg)

        times: List[Optional[time]] = []
        if format is not None:
            for element in arg:
                try:
                    times.append(datetime.strptime(element, format).time())
                except (ValueError, TypeError) as err:
                    if errors == "raise":
                        msg = (
                            f"Cannot convert {element} to a time with given "
                            f"format {format}"
                        )
                        raise ValueError(msg) from err
                    elif errors == "ignore":
                        return arg
                    else:
                        times.append(None)
        else:
            formats = _time_formats[:]
            format_found = False
            for element in arg:
                time_object = None
                for time_format in formats:
                    try:
                        time_object = datetime.strptime(element, time_format).time()
                        if not format_found:
                            # Put the found format in front
                            fmt = formats.pop(formats.index(time_format))
                            formats.insert(0, fmt)
                            format_found = True
                        break
                    except (ValueError, TypeError):
                        continue

                if time_object is not None:
                    times.append(time_object)
                elif errors == "raise":
                    raise ValueError(f"Cannot convert arg {arg} to a time")
                elif errors == "ignore":
                    return arg
                else:
                    times.append(None)

        return times

    if arg is None:
        return arg
    elif isinstance(arg, time):
        return arg
    elif isinstance(arg, ABCSeries):
        values = _convert_listlike(arg._values, format)
        return arg._constructor(values, index=arg.index, name=arg.name)
    elif isinstance(arg, ABCIndexClass):
        return _convert_listlike(arg, format)
    elif is_list_like(arg):
        return _convert_listlike(arg, format)

    return _convert_listlike(np.array([arg]), format)[0]
