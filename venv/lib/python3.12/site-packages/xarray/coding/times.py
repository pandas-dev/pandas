from __future__ import annotations

import contextlib
import re
import warnings
from collections.abc import Callable, Hashable
from datetime import datetime, timedelta
from functools import partial
from typing import TYPE_CHECKING, Union, cast

import numpy as np
import pandas as pd
from pandas.errors import OutOfBoundsDatetime, OutOfBoundsTimedelta

from xarray.coding.common import (
    SerializationWarning,
    VariableCoder,
    lazy_elemwise_func,
    pop_to,
    safe_setitem,
    unpack_for_decoding,
    unpack_for_encoding,
)
from xarray.compat.pdcompat import default_precision_timestamp, timestamp_as_unit
from xarray.core import indexing
from xarray.core.common import contains_cftime_datetimes, is_np_datetime_like
from xarray.core.duck_array_ops import array_all, asarray, ravel, reshape
from xarray.core.formatting import first_n_items, format_timestamp, last_item
from xarray.core.utils import attempt_import, emit_user_level_warning
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import T_ChunkedArray, get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array, to_numpy
from xarray.namedarray.utils import is_duck_dask_array

try:
    import cftime
except ImportError:
    cftime = None

from xarray.core.types import (
    CFCalendar,
    CFTimeDatetime,
    NPDatetimeUnitOptions,
    PDDatetimeUnitOptions,
    T_DuckArray,
)

T_Name = Union[Hashable, None]

# standard calendars recognized by cftime
_STANDARD_CALENDARS = {"standard", "gregorian", "proleptic_gregorian"}

_NS_PER_TIME_DELTA = {
    "ns": 1,
    "us": int(1e3),
    "ms": int(1e6),
    "s": int(1e9),
    "m": int(1e9) * 60,
    "h": int(1e9) * 60 * 60,
    "D": int(1e9) * 60 * 60 * 24,
}

_US_PER_TIME_DELTA = {
    "microseconds": 1,
    "milliseconds": 1_000,
    "seconds": 1_000_000,
    "minutes": 60 * 1_000_000,
    "hours": 60 * 60 * 1_000_000,
    "days": 24 * 60 * 60 * 1_000_000,
}

_NETCDF_TIME_UNITS_CFTIME = [
    "days",
    "hours",
    "minutes",
    "seconds",
    "milliseconds",
    "microseconds",
]

_NETCDF_TIME_UNITS_NUMPY = _NETCDF_TIME_UNITS_CFTIME + ["nanoseconds"]

TIME_UNITS = frozenset(
    [
        "days",
        "hours",
        "minutes",
        "seconds",
        "milliseconds",
        "microseconds",
        "nanoseconds",
    ]
)


_INVALID_LITERAL_TIMEDELTA64_ENCODING_KEYS = [
    "add_offset",
    "scale_factor",
]


def _is_standard_calendar(calendar: str) -> bool:
    return calendar.lower() in _STANDARD_CALENDARS


def _is_numpy_compatible_time_range(times):
    if is_np_datetime_like(times.dtype):
        return True
    # times array contains cftime objects
    times = np.asarray(times)
    tmin = times.min()
    tmax = times.max()
    try:
        # before relaxing the nanosecond constrained
        # this raised OutOfBoundsDatetime for
        # times < 1678 and times > 2262
        # this isn't the case anymore for other resolutions like "s"
        # now, we raise for dates before 1582-10-15
        _check_date_is_after_shift(tmin, "standard")
        _check_date_is_after_shift(tmax, "standard")
        convert_time_or_go_back(tmin, pd.Timestamp)
        convert_time_or_go_back(tmax, pd.Timestamp)
    except pd.errors.OutOfBoundsDatetime:
        return False
    except ValueError as err:
        if err.args[0] == "year 0 is out of range":
            return False
        raise
    else:
        return True


def _netcdf_to_numpy_timeunit(units: str) -> NPDatetimeUnitOptions:
    units = units.lower()
    if not units.endswith("s"):
        units = f"{units}s"
    return cast(
        NPDatetimeUnitOptions,
        {
            "nanoseconds": "ns",
            "microseconds": "us",
            "milliseconds": "ms",
            "seconds": "s",
            "minutes": "m",
            "hours": "h",
            "days": "D",
        }[units],
    )


def _numpy_to_netcdf_timeunit(units: NPDatetimeUnitOptions) -> str:
    return {
        "ns": "nanoseconds",
        "us": "microseconds",
        "ms": "milliseconds",
        "s": "seconds",
        "m": "minutes",
        "h": "hours",
        "D": "days",
    }[units]


def _numpy_dtype_to_netcdf_timeunit(dtype: np.dtype) -> str:
    unit, _ = np.datetime_data(dtype)
    unit = cast(NPDatetimeUnitOptions, unit)
    return _numpy_to_netcdf_timeunit(unit)


def _ensure_padded_year(ref_date: str) -> str:
    # Reference dates without a padded year (e.g. since 1-1-1 or since 2-3-4)
    # are ambiguous (is it YMD or DMY?). This can lead to some very odd
    # behaviour e.g. pandas (via dateutil) passes '1-1-1 00:00:0.0' as
    # '2001-01-01 00:00:00' (because it assumes a) DMY and b) that year 1 is
    # shorthand for 2001 (like 02 would be shorthand for year 2002)).

    # Here we ensure that there is always a four-digit year, with the
    # assumption being that year comes first if we get something ambiguous.
    matches_year = re.match(r".*\d{4}.*", ref_date)
    if matches_year:
        # all good, return
        return ref_date

    # No four-digit strings, assume the first digits are the year and pad
    # appropriately
    matches_start_digits = re.match(r"(\d+)(.*)", ref_date)
    if not matches_start_digits:
        raise ValueError(f"invalid reference date for time units: {ref_date}")
    ref_year, everything_else = (s for s in matches_start_digits.groups())
    ref_date_padded = f"{int(ref_year):04d}{everything_else}"

    warning_msg = (
        f"Ambiguous reference date string: {ref_date}. The first value is "
        "assumed to be the year hence will be padded with zeros to remove "
        f"the ambiguity (the padded reference date string is: {ref_date_padded}). "
        "To remove this message, remove the ambiguity by padding your reference "
        "date strings with zeros."
    )
    warnings.warn(warning_msg, SerializationWarning, stacklevel=2)

    return ref_date_padded


def _unpack_netcdf_time_units(units: str) -> tuple[str, str]:
    # CF datetime units follow the format: "UNIT since DATE"
    # this parses out the unit and date allowing for extraneous
    # whitespace. It also ensures that the year is padded with zeros
    # so it will be correctly understood by pandas (via dateutil).
    matches = re.match(r"(.+) since (.+)", units)
    if not matches:
        raise ValueError(f"invalid time units: {units}")

    delta_units, ref_date = (s.strip() for s in matches.groups())
    ref_date = _ensure_padded_year(ref_date)

    return delta_units, ref_date


def named(name: str, pattern: str) -> str:
    return "(?P<" + name + ">" + pattern + ")"


def optional(x: str) -> str:
    return "(?:" + x + ")?"


def trailing_optional(xs: list[str]) -> str:
    if not xs:
        return ""
    return xs[0] + optional(trailing_optional(xs[1:]))


def build_pattern(
    date_sep: str = r"\-",
    datetime_sep: str = r"T",
    time_sep: str = r"\:",
    micro_sep: str = r".",
) -> str:
    pieces = [
        (None, "year", r"[+-]?\d{4,5}"),
        (date_sep, "month", r"\d{2}"),
        (date_sep, "day", r"\d{2}"),
        (datetime_sep, "hour", r"\d{2}"),
        (time_sep, "minute", r"\d{2}"),
        (time_sep, "second", r"\d{2}"),
        (micro_sep, "microsecond", r"\d{1,6}"),
    ]
    pattern_list = []
    for sep, name, sub_pattern in pieces:
        pattern_list.append((sep or "") + named(name, sub_pattern))
        # TODO: allow timezone offsets?
    return "^" + trailing_optional(pattern_list) + "$"


_BASIC_PATTERN = build_pattern(date_sep="", time_sep="")
_EXTENDED_PATTERN = build_pattern()
_CFTIME_PATTERN = build_pattern(datetime_sep=" ")
_PATTERNS = [_BASIC_PATTERN, _EXTENDED_PATTERN, _CFTIME_PATTERN]


def parse_iso8601_like(datetime_string: str) -> dict[str, str | None]:
    for pattern in _PATTERNS:
        match = re.match(pattern, datetime_string)
        if match:
            return match.groupdict()
    raise ValueError(
        f"no ISO-8601 or cftime-string-like match for string: {datetime_string}"
    )


def _parse_iso8601(date_type, timestr):
    default = date_type(1, 1, 1)
    result = parse_iso8601_like(timestr)
    replace = {}

    for attr in ["year", "month", "day", "hour", "minute", "second", "microsecond"]:
        value = result.get(attr, None)
        if value is not None:
            resolution = attr
            if attr == "microsecond":
                if len(value) <= 3:
                    resolution = "millisecond"
                # convert match string into valid microsecond value
                value = 10 ** (6 - len(value)) * int(value)
            replace[attr] = int(value)
    return default.replace(**replace), resolution


def _maybe_strip_tz_from_timestamp(date: pd.Timestamp) -> pd.Timestamp:
    # If the ref_date Timestamp is timezone-aware, convert to UTC and
    # make it timezone-naive (GH 2649).
    if date.tz is not None:
        return date.tz_convert("UTC").tz_convert(None)
    return date


def _unpack_time_unit_and_ref_date(
    units: str,
) -> tuple[NPDatetimeUnitOptions, pd.Timestamp]:
    # same us _unpack_netcdf_time_units but finalizes ref_date for
    # processing in encode_cf_datetime
    time_unit, _ref_date = _unpack_netcdf_time_units(units)
    time_unit = _netcdf_to_numpy_timeunit(time_unit)
    ref_date = pd.Timestamp(_ref_date)
    ref_date = _maybe_strip_tz_from_timestamp(ref_date)
    return time_unit, ref_date


def _unpack_time_units_and_ref_date_cftime(units: str, calendar: str):
    # same as _unpack_netcdf_time_units but finalizes ref_date for
    # processing in encode_cf_datetime
    time_units, ref_date = _unpack_netcdf_time_units(units)
    ref_date = cftime.num2date(
        0,
        units=f"microseconds since {ref_date}",
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )
    return time_units, ref_date


def _decode_cf_datetime_dtype(
    data,
    units: str,
    calendar: str | None,
    use_cftime: bool | None,
    time_unit: PDDatetimeUnitOptions = "ns",
) -> np.dtype:
    # Verify that at least the first and last date can be decoded
    # successfully. Otherwise, tracebacks end up swallowed by
    # Dataset.__repr__ when users try to view their lazily decoded array.
    values = indexing.ImplicitToExplicitIndexingAdapter(indexing.as_indexable(data))
    example_value = np.concatenate(
        [to_numpy(first_n_items(values, 1)), to_numpy(last_item(values))]
    )

    try:
        result = decode_cf_datetime(
            example_value, units, calendar, use_cftime, time_unit
        )
    except Exception as err:
        calendar_msg = (
            "the default calendar" if calendar is None else f"calendar {calendar!r}"
        )
        msg = (
            f"unable to decode time units {units!r} with {calendar_msg!r}. Try "
            "opening your dataset with decode_times=False or installing cftime "
            "if it is not installed."
        )
        raise ValueError(msg) from err
    else:
        dtype = getattr(result, "dtype", np.dtype("object"))

    return dtype


def _decode_datetime_with_cftime(
    num_dates: np.ndarray, units: str, calendar: str
) -> np.ndarray:
    if TYPE_CHECKING:
        import cftime
    else:
        cftime = attempt_import("cftime")
    if num_dates.size > 0:
        return np.asarray(
            cftime.num2date(num_dates, units, calendar, only_use_cftime_datetimes=True)
        )
    else:
        return np.array([], dtype=object)


def _check_date_for_units_since_refdate(
    date, unit: NPDatetimeUnitOptions, ref_date: pd.Timestamp
) -> pd.Timestamp:
    # check for out-of-bounds floats and raise
    if date > np.iinfo("int64").max or date < np.iinfo("int64").min:
        raise OutOfBoundsTimedelta(
            f"Value {date} can't be represented as Datetime/Timedelta."
        )
    delta = date * np.timedelta64(1, unit)
    if not np.isnan(delta):
        # this will raise on dtype overflow for integer dtypes
        if date.dtype.kind in "u" and not np.int64(delta) == date:
            raise OutOfBoundsTimedelta(
                "DType overflow in Datetime/Timedelta calculation."
            )
        # this will raise on overflow if ref_date + delta
        # can't be represented in the current ref_date resolution
        return timestamp_as_unit(ref_date + delta, ref_date.unit)
    else:
        # if date is exactly NaT (np.iinfo("int64").min) return NaT
        # to make follow-up checks work
        return pd.Timestamp("NaT")


def _check_timedelta_range(value, data_unit, time_unit):
    if value > np.iinfo("int64").max or value < np.iinfo("int64").min:
        OutOfBoundsTimedelta(f"Value {value} can't be represented as Timedelta.")
    # on windows multiplying nan leads to RuntimeWarning
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "invalid value encountered in multiply", RuntimeWarning
        )
        delta = value * np.timedelta64(1, data_unit)
    if not np.isnan(delta):
        # this will raise on dtype overflow for integer dtypes
        if value.dtype.kind in "u" and not np.int64(delta) == value:
            raise OutOfBoundsTimedelta(
                "DType overflow in Datetime/Timedelta calculation."
            )
        # this will raise on overflow if delta cannot be represented with the
        # resolutions supported by pandas.
        pd.to_timedelta(delta)


def _align_reference_date_and_unit(
    ref_date: pd.Timestamp, unit: NPDatetimeUnitOptions
) -> pd.Timestamp:
    # align to the highest needed resolution of ref_date or unit
    if np.timedelta64(1, ref_date.unit) > np.timedelta64(1, unit):
        # this will raise accordingly
        # if data can't be represented in the higher resolution
        return timestamp_as_unit(ref_date, cast(PDDatetimeUnitOptions, unit))
    return ref_date


def _check_date_is_after_shift(
    date: pd.Timestamp | datetime | CFTimeDatetime, calendar: str
) -> None:
    # if we have gregorian/standard we need to raise
    # if we are outside the well-defined date range
    # proleptic_gregorian and standard/gregorian are only equivalent
    # if reference date and date range is >= 1582-10-15
    if calendar != "proleptic_gregorian" and date < type(date)(1582, 10, 15):
        raise OutOfBoundsDatetime(
            f"Dates before 1582-10-15 cannot be decoded "
            f"with pandas using {calendar!r} calendar: {date}"
        )


def _check_higher_resolution(
    flat_num_dates: np.ndarray,
    time_unit: PDDatetimeUnitOptions,
) -> tuple[np.ndarray, PDDatetimeUnitOptions]:
    """Iterate until fitting resolution found."""
    res: list[PDDatetimeUnitOptions] = ["s", "ms", "us", "ns"]
    new_units = res[res.index(time_unit) :]
    for new_time_unit in new_units:
        if not ((np.unique(flat_num_dates % 1) > 0).any() and new_time_unit != "ns"):
            break
        flat_num_dates *= 1000
    return flat_num_dates, new_time_unit


def _decode_datetime_with_pandas(
    flat_num_dates: np.ndarray,
    units: str,
    calendar: str,
    time_resolution: PDDatetimeUnitOptions = "ns",
) -> np.ndarray:
    if not _is_standard_calendar(calendar):
        raise OutOfBoundsDatetime(
            f"Cannot decode times from a non-standard calendar, {calendar!r}, using "
            "pandas."
        )

    # Work around pandas.to_timedelta issue with dtypes smaller than int64 and
    # NumPy 2.0 by casting all int and uint data to int64 and uint64,
    # respectively. See https://github.com/pandas-dev/pandas/issues/56996 for
    # more details.
    if flat_num_dates.dtype.kind == "i":
        flat_num_dates = flat_num_dates.astype(np.int64)
    elif flat_num_dates.dtype.kind == "u":
        flat_num_dates = flat_num_dates.astype(np.uint64)

    try:
        time_unit, ref_date = _unpack_time_unit_and_ref_date(units)
        ref_date = _align_reference_date_and_unit(ref_date, time_unit)
        # here the highest wanted resolution is set
        ref_date = _align_reference_date_and_unit(ref_date, time_resolution)
    except ValueError as err:
        # ValueError is raised by pd.Timestamp for non-ISO timestamp
        # strings, in which case we fall back to using cftime
        raise OutOfBoundsDatetime from err

    _check_date_is_after_shift(ref_date, calendar)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "invalid value encountered", RuntimeWarning)
        if flat_num_dates.size > 0:
            # avoid size 0 datetimes GH1329
            _check_date_for_units_since_refdate(
                flat_num_dates.min(), time_unit, ref_date
            )
            _check_date_for_units_since_refdate(
                flat_num_dates.max(), time_unit, ref_date
            )

    # To avoid integer overflow when converting to nanosecond units for integer
    # dtypes smaller than np.int64 cast all integer and unsigned integer dtype
    # arrays to np.int64 (GH 2002, GH 6589).  Note this is safe even in the case
    # of np.uint64 values, because any np.uint64 value that would lead to
    # overflow when converting to np.int64 would not be representable with a
    # timedelta64 value, and therefore would raise an error in the lines above.
    if flat_num_dates.dtype.kind in "iu":
        flat_num_dates = flat_num_dates.astype(np.int64)
    elif flat_num_dates.dtype.kind in "f":
        flat_num_dates = flat_num_dates.astype(np.float64)

    timedeltas = _numbers_to_timedelta(
        flat_num_dates, time_unit, ref_date.unit, "datetimes"
    )

    # add timedeltas to ref_date
    return ref_date + timedeltas


def decode_cf_datetime(
    num_dates,
    units: str,
    calendar: str | None = None,
    use_cftime: bool | None = None,
    time_unit: PDDatetimeUnitOptions = "ns",
) -> np.ndarray:
    """Given an array of numeric dates in netCDF format, convert it into a
    numpy array of date time objects.

    For standard (Gregorian) calendars, this function uses vectorized
    operations, which makes it much faster than cftime.num2date. In such a
    case, the returned array will be of type np.datetime64.

    Note that time unit in `units` must not be smaller than microseconds and
    not larger than days.

    See Also
    --------
    cftime.num2date
    """
    num_dates = to_numpy(num_dates)
    flat_num_dates = ravel(num_dates)
    if calendar is None:
        calendar = "standard"

    if use_cftime is None:
        try:
            dates = _decode_datetime_with_pandas(
                flat_num_dates, units, calendar, time_unit
            )
        except (KeyError, OutOfBoundsDatetime, OutOfBoundsTimedelta, OverflowError):
            dates = _decode_datetime_with_cftime(
                flat_num_dates.astype(float), units, calendar
            )
            # retrieve cftype
            dates_min = dates[np.nanargmin(num_dates)]
            dates_max = dates[np.nanargmax(num_dates)]
            cftype = type(dates_min)
            # create first day of gregorian calendar in current cf calendar type
            border = cftype(1582, 10, 15)
            # "ns" borders
            # between ['1677-09-21T00:12:43.145224193', '2262-04-11T23:47:16.854775807']
            lower = cftype(1677, 9, 21, 0, 12, 43, 145224)
            upper = cftype(2262, 4, 11, 23, 47, 16, 854775)

            if dates_min < border:
                if _is_standard_calendar(calendar):
                    emit_user_level_warning(
                        "Unable to decode time axis into full "
                        "numpy.datetime64 objects, continuing using "
                        "cftime.datetime objects instead, reason: dates prior "
                        "reform date (1582-10-15). To silence this warning specify "
                        "'use_cftime=True'.",
                        SerializationWarning,
                    )
            elif time_unit == "ns" and (dates_min < lower or dates_max > upper):
                emit_user_level_warning(
                    "Unable to decode time axis into full "
                    "numpy.datetime64[ns] objects, continuing using "
                    "cftime.datetime objects instead, reason: dates out "
                    "of range. To silence this warning use a coarser resolution "
                    "'time_unit' or specify 'use_cftime=True'.",
                    SerializationWarning,
                )
            elif _is_standard_calendar(calendar):
                dates = cftime_to_nptime(dates, time_unit=time_unit)
    elif use_cftime:
        dates = _decode_datetime_with_cftime(flat_num_dates, units, calendar)
    else:
        dates = _decode_datetime_with_pandas(flat_num_dates, units, calendar, time_unit)

    return reshape(dates, num_dates.shape)


def to_datetime_unboxed(value, **kwargs):
    result = pd.to_datetime(value, **kwargs).to_numpy()
    assert np.issubdtype(result.dtype, "datetime64")
    return result


def _numbers_to_timedelta(
    flat_num: np.ndarray,
    time_unit: NPDatetimeUnitOptions,
    ref_unit: PDDatetimeUnitOptions,
    datatype: str,
    target_unit: PDDatetimeUnitOptions | None = None,
) -> np.ndarray:
    """Transform numbers to np.timedelta64."""
    # keep NaT/nan mask
    if flat_num.dtype.kind == "f":
        nan = np.asarray(np.isnan(flat_num))
    elif flat_num.dtype.kind == "i":
        nan = np.asarray(flat_num == np.iinfo(np.int64).min)

    # in case we need to change the unit, we fix the numbers here
    # this should be safe, as errors would have been raised above
    ns_time_unit = _NS_PER_TIME_DELTA[time_unit]
    ns_ref_date_unit = _NS_PER_TIME_DELTA[ref_unit]
    if ns_time_unit > ns_ref_date_unit:
        flat_num = np.asarray(flat_num * np.int64(ns_time_unit / ns_ref_date_unit))
        time_unit = ref_unit

    # estimate fitting resolution for floating point values
    # this iterates until all floats are fractionless or time_unit == "ns"
    if flat_num.dtype.kind == "f" and time_unit != "ns":
        flat_num, new_time_unit = _check_higher_resolution(
            flat_num, cast(PDDatetimeUnitOptions, time_unit)
        )
        if time_unit != new_time_unit:
            if target_unit is None or np.timedelta64(1, target_unit) > np.timedelta64(
                1, new_time_unit
            ):
                if datatype == "datetimes":
                    kwarg = "decode_times"
                    coder = "CFDatetimeCoder"
                else:
                    kwarg = "decode_timedelta"
                    coder = "CFTimedeltaCoder"
                formatted_kwarg = f"{kwarg}={coder}(time_unit={new_time_unit!r})"
                message = (
                    f"Can't decode floating point {datatype} to {time_unit!r} "
                    f"without precision loss; decoding to {new_time_unit!r} "
                    f"instead. To silence this warning pass {formatted_kwarg} "
                    f"to your opening function."
                )
                emit_user_level_warning(message, SerializationWarning)
            time_unit = new_time_unit

    # Cast input ordinals to integers and properly handle NaN/NaT
    # to prevent casting NaN to int
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        flat_num = flat_num.astype(np.int64)
    if nan.any():
        flat_num[nan] = np.iinfo(np.int64).min

    # cast to wanted type
    return flat_num.astype(f"timedelta64[{time_unit}]")


def decode_cf_timedelta(
    num_timedeltas, units: str, time_unit: PDDatetimeUnitOptions = "ns"
) -> np.ndarray:
    """Given an array of numeric timedeltas in netCDF format, convert it into a
    numpy timedelta64 ["s", "ms", "us", "ns"] array.
    """
    num_timedeltas = to_numpy(num_timedeltas)
    unit = _netcdf_to_numpy_timeunit(units)

    # special case empty arrays
    is_empty_array = num_timedeltas.size == 0

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN slice encountered", RuntimeWarning)
        if not is_empty_array:
            _check_timedelta_range(np.nanmin(num_timedeltas), unit, time_unit)
            _check_timedelta_range(np.nanmax(num_timedeltas), unit, time_unit)

    timedeltas = _numbers_to_timedelta(
        num_timedeltas, unit, "s", "timedeltas", target_unit=time_unit
    )
    pd_timedeltas = pd.to_timedelta(ravel(timedeltas))

    if not is_empty_array and np.isnat(timedeltas).all():
        empirical_unit = time_unit
    else:
        empirical_unit = pd_timedeltas.unit

    if is_empty_array or np.timedelta64(1, time_unit) > np.timedelta64(
        1, empirical_unit
    ):
        time_unit = empirical_unit

    if time_unit not in {"s", "ms", "us", "ns"}:
        raise ValueError(
            f"time_unit must be one of 's', 'ms', 'us', or 'ns'. Got: {time_unit}"
        )

    result = pd_timedeltas.as_unit(time_unit).to_numpy()
    return reshape(result, num_timedeltas.shape)


def _unit_timedelta_cftime(units: str) -> timedelta:
    return timedelta(microseconds=_US_PER_TIME_DELTA[units])


def _unit_timedelta_numpy(units: str) -> np.timedelta64:
    numpy_units = _netcdf_to_numpy_timeunit(units)
    return np.timedelta64(1, numpy_units)


def _infer_time_units_from_diff(unique_timedeltas) -> str:
    # todo: check, if this function works correctly wrt np.timedelta64
    unit_timedelta: Callable[[str], timedelta] | Callable[[str], np.timedelta64]
    zero_timedelta: timedelta | np.timedelta64
    unique_timedeltas = asarray(unique_timedeltas)
    if unique_timedeltas.dtype == np.dtype("O"):
        time_units = _NETCDF_TIME_UNITS_CFTIME
        unit_timedelta = _unit_timedelta_cftime
        zero_timedelta = timedelta(microseconds=0)
    else:
        time_units = _NETCDF_TIME_UNITS_NUMPY
        unit_timedelta = _unit_timedelta_numpy
        zero_timedelta = np.timedelta64(0, "ns")
    for time_unit in time_units:
        if array_all(unique_timedeltas % unit_timedelta(time_unit) == zero_timedelta):
            return time_unit
    return "seconds"


def _time_units_to_timedelta(units: str) -> timedelta:
    return timedelta(microseconds=_US_PER_TIME_DELTA[units])


def infer_calendar_name(dates) -> CFCalendar:
    """Given an array of datetimes, infer the CF calendar name"""
    if is_np_datetime_like(dates.dtype):
        return "proleptic_gregorian"
    elif dates.dtype == np.dtype("O") and dates.size > 0:
        # Logic copied from core.common.contains_cftime_datetimes.
        if cftime is not None:
            sample = np.asarray(dates).flat[0]
            if is_duck_dask_array(sample):
                sample = sample.compute()
                if isinstance(sample, np.ndarray):
                    sample = sample.item()
            if isinstance(sample, cftime.datetime):
                return sample.calendar

    # Error raise if dtype is neither datetime or "O", if cftime is not importable, and if element of 'O' dtype is not cftime.
    raise ValueError("Array does not contain datetime objects.")


def infer_datetime_units(dates) -> str:
    """Given an array of datetimes, returns a CF compatible time-unit string of
    the form "{time_unit} since {date[0]}", where `time_unit` is 'days',
    'hours', 'minutes' or 'seconds' (the first one that can evenly divide all
    unique time deltas in `dates`)
    """
    dates = ravel(np.asarray(dates))
    if np.issubdtype(np.asarray(dates).dtype, "datetime64"):
        dates = to_datetime_unboxed(dates)
        dates = dates[pd.notnull(dates)]
        reference_date = dates[0] if len(dates) > 0 else "1970-01-01"
        reference_date = pd.Timestamp(reference_date)
    else:
        reference_date = dates[0] if len(dates) > 0 else "1970-01-01"
        reference_date = format_cftime_datetime(reference_date)
    unique_timedeltas = np.unique(np.diff(dates))
    units = _infer_time_units_from_diff(unique_timedeltas)
    return f"{units} since {reference_date}"


def format_cftime_datetime(date) -> str:
    """Converts a cftime.datetime object to a string with the format:
    YYYY-MM-DD HH:MM:SS.UUUUUU
    """
    return f"{date.year:04d}-{date.month:02d}-{date.day:02d} {date.hour:02d}:{date.minute:02d}:{date.second:02d}.{date.microsecond:06d}"


def infer_timedelta_units(deltas) -> str:
    """Given an array of timedeltas, returns a CF compatible time-unit from
    {'days', 'hours', 'minutes' 'seconds'} (the first one that can evenly
    divide all unique time deltas in `deltas`)
    """
    deltas = ravel(deltas)
    unique_timedeltas = np.unique(deltas[pd.notnull(deltas)])
    return _infer_time_units_from_diff(unique_timedeltas)


def cftime_to_nptime(
    times, raise_on_invalid: bool = True, time_unit: PDDatetimeUnitOptions = "ns"
) -> np.ndarray:
    """Given an array of cftime.datetime objects, return an array of
    numpy.datetime64 objects of the same size

    If raise_on_invalid is True (default), invalid dates trigger a ValueError.
    Otherwise, the invalid element is replaced by np.NaT."""
    times = np.asarray(times)
    new = []
    dt: np.datetime64
    for _i, t in np.ndenumerate(times):
        try:
            # We expect either "us" resolution or "s" resolution depending on
            # whether 'microseconds' are defined for the input or not.
            dt = (
                pd.Timestamp(np.datetime64(t.isoformat())).as_unit(time_unit).to_numpy()
            )
        except ValueError as e:
            if raise_on_invalid:
                raise ValueError(
                    f"Cannot convert date {t} to a date in the "
                    f"standard calendar.  Reason: {e}."
                ) from e
            else:
                dt = np.datetime64("NaT")
        new.append(dt)
    return np.asarray(new).reshape(times.shape)


def convert_times(times, date_type, raise_on_invalid: bool = True) -> np.ndarray:
    """Given an array of datetimes, return the same dates in another cftime or numpy date type.

    Useful to convert between calendars in numpy and cftime or between cftime calendars.

    If raise_on_valid is True (default), invalid dates trigger a ValueError.
    Otherwise, the invalid element is replaced by np.nan for cftime types and np.NaT for np.datetime64.
    """
    if date_type in (pd.Timestamp, np.datetime64) and not is_np_datetime_like(
        times.dtype
    ):
        return cftime_to_nptime(times, raise_on_invalid=raise_on_invalid)
    if is_np_datetime_like(times.dtype):
        # Convert datetime64 objects to Timestamps since those have year, month, day, etc. attributes
        times = pd.DatetimeIndex(times)
    new = np.empty(times.shape, dtype="O")
    for i, t in enumerate(times):
        try:
            dt = date_type(
                t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond
            )
        except ValueError as e:
            if raise_on_invalid:
                raise ValueError(
                    f"Cannot convert date {t} to a date in the "
                    f"{date_type(2000, 1, 1).calendar} calendar.  Reason: {e}."
                ) from e
            else:
                dt = np.nan

        new[i] = dt
    return new


def convert_time_or_go_back(date, date_type):
    """Convert a single date to a new date_type (cftime.datetime or pd.Timestamp).

    If the new date is invalid, it goes back a day and tries again. If it is still
    invalid, goes back a second day.

    This is meant to convert end-of-month dates into a new calendar.
    """
    if date_type == pd.Timestamp:
        date_type = default_precision_timestamp
    try:
        return date_type(
            date.year,
            date.month,
            date.day,
            date.hour,
            date.minute,
            date.second,
            date.microsecond,
        )
    except OutOfBoundsDatetime:
        raise
    except ValueError:
        # Day is invalid, happens at the end of months, try again the day before
        try:
            return date_type(
                date.year,
                date.month,
                date.day - 1,
                date.hour,
                date.minute,
                date.second,
                date.microsecond,
            )
        except ValueError:
            # Still invalid, happens for 360_day to non-leap february. Try again 2 days before date.
            return date_type(
                date.year,
                date.month,
                date.day - 2,
                date.hour,
                date.minute,
                date.second,
                date.microsecond,
            )


def _should_cftime_be_used(
    source, target_calendar: str, use_cftime: bool | None
) -> bool:
    """Return whether conversion of the source to the target calendar should
    result in a cftime-backed array.

    Source is a 1D datetime array, target_cal a string (calendar name) and
    use_cftime is a boolean or None. If use_cftime is None, this returns True
    if the source's range and target calendar are convertible to np.datetime64 objects.
    """
    # Arguments Checks for target
    if use_cftime is not True:
        if _is_standard_calendar(target_calendar):
            if _is_numpy_compatible_time_range(source):
                # Conversion is possible with pandas, force False if it was None
                return False
            elif use_cftime is False:
                raise ValueError(
                    "Source time range is not valid for numpy datetimes. Try using `use_cftime=True`."
                )
        elif use_cftime is False:
            raise ValueError(
                f"Calendar '{target_calendar}' is only valid with cftime. Try using `use_cftime=True`."
            )
    return True


def _cleanup_netcdf_time_units(units: str) -> str:
    time_units, ref_date = _unpack_netcdf_time_units(units)
    time_units = time_units.lower()
    if not time_units.endswith("s"):
        time_units = f"{time_units}s"
    # don't worry about reifying the units if they're out of bounds or
    # formatted badly
    with contextlib.suppress(OutOfBoundsDatetime, ValueError):
        units = f"{time_units} since {format_timestamp(ref_date)}"
    return units


def _encode_datetime_with_cftime(dates, units: str, calendar: str) -> np.ndarray:
    """Fallback method for encoding dates using cftime.

    This method is more flexible than xarray's parsing using datetime64[ns]
    arrays but also slower because it loops over each element.
    """
    if TYPE_CHECKING:
        import cftime
    else:
        cftime = attempt_import("cftime")

    dates = np.asarray(dates)
    original_shape = dates.shape

    if np.issubdtype(dates.dtype, np.datetime64):
        # numpy's broken datetime conversion only works for us precision
        dates = dates.astype("M8[us]").astype(datetime)

    dates = np.atleast_1d(dates)

    # Find all the None position
    none_position = dates == None  # noqa: E711
    filtered_dates = dates[~none_position]

    # Since netCDF files do not support storing float128 values, we ensure
    # that float64 values are used by setting longdouble=False in num2date.
    # This try except logic can be removed when xarray's minimum version of
    # cftime is at least 1.6.2.
    try:
        encoded_nums = cftime.date2num(
            filtered_dates, units, calendar, longdouble=False
        )
    except TypeError:
        encoded_nums = cftime.date2num(filtered_dates, units, calendar)

    if filtered_dates.size == none_position.size:
        return encoded_nums.reshape(original_shape)

    # Create a full matrix of NaN
    # And fill the num dates in the not NaN or None position
    result = np.full(dates.shape, np.nan)
    result[np.nonzero(~none_position)] = encoded_nums
    return result.reshape(original_shape)


def cast_to_int_if_safe(num) -> np.ndarray:
    int_num = np.asarray(num, dtype=np.int64)
    if array_all(num == int_num):
        num = int_num
    return num


def _division(deltas, delta, floor):
    if floor:
        # calculate int64 floor division
        # to preserve integer dtype if possible (GH 4045, GH7817).
        num = deltas // delta.astype(np.int64)
        num = num.astype(np.int64, copy=False)
    else:
        num = deltas / delta
    return num


def encode_cf_datetime(
    dates: T_DuckArray,  # type: ignore[misc]
    units: str | None = None,
    calendar: str | None = None,
    dtype: np.dtype | None = None,
) -> tuple[T_DuckArray, str, str]:
    """Given an array of datetime objects, returns the tuple `(num, units,
    calendar)` suitable for a CF compliant time variable.

    Unlike `date2num`, this function can handle datetime64 arrays.

    See Also
    --------
    cftime.date2num
    """
    dates = asarray(dates)
    if is_chunked_array(dates):
        return _lazily_encode_cf_datetime(dates, units, calendar, dtype)
    else:
        return _eagerly_encode_cf_datetime(dates, units, calendar, dtype)


def _infer_needed_units_numpy(ref_date, data_units):
    needed_units, data_ref_date = _unpack_time_unit_and_ref_date(data_units)
    needed_units = _numpy_to_netcdf_timeunit(needed_units)
    ref_delta = abs(data_ref_date - ref_date).to_timedelta64()
    data_delta = _unit_timedelta_numpy(needed_units)
    if (ref_delta % data_delta) > np.timedelta64(0, "ns"):
        needed_units = _infer_time_units_from_diff(ref_delta)
    return needed_units


def _infer_needed_units_cftime(ref_date, data_units, calendar):
    needed_units, data_ref_date = _unpack_time_units_and_ref_date_cftime(
        data_units, calendar
    )
    ref_delta = abs(data_ref_date - ref_date)
    data_delta = _time_units_to_timedelta(needed_units)
    if (ref_delta % data_delta) > timedelta(seconds=0):
        needed_units = _infer_time_units_from_diff(ref_delta)
    return needed_units


def _eagerly_encode_cf_datetime(
    dates: T_DuckArray,  # type: ignore[misc]
    units: str | None = None,
    calendar: str | None = None,
    dtype: np.dtype | None = None,
    allow_units_modification: bool = True,
) -> tuple[T_DuckArray, str, str]:
    dates = asarray(dates)
    data_units = infer_datetime_units(dates)
    if units is None:
        units = data_units
    else:
        units = _cleanup_netcdf_time_units(units)

    if calendar is None:
        calendar = infer_calendar_name(dates)

    raise_incompatible_units_error = False
    raise_gregorian_proleptic_gregorian_mismatch_error = False
    try:
        if not _is_standard_calendar(calendar) or dates.dtype.kind == "O":
            # parse with cftime instead
            raise OutOfBoundsDatetime
        assert np.issubdtype(dates.dtype, "datetime64")
        if calendar in ["standard", "gregorian"] and np.nanmin(dates).astype(
            "=M8[us]"
        ).astype(datetime) < datetime(1582, 10, 15):
            raise_gregorian_proleptic_gregorian_mismatch_error = True

        time_unit, ref_date = _unpack_time_unit_and_ref_date(units)
        # calendar equivalence only for days after the reform
        _check_date_is_after_shift(ref_date, calendar)
        time_delta = np.timedelta64(1, time_unit)

        # Wrap the dates in a DatetimeIndex to do the subtraction to ensure
        # an OverflowError is raised if the ref_date is too far away from
        # dates to be encoded (GH 2272).
        # DatetimeIndex will convert to units of ["s", "ms", "us", "ns"]
        dates_as_index = pd.DatetimeIndex(ravel(dates))
        time_deltas = dates_as_index - ref_date

        # retrieve needed units to faithfully encode to int64
        needed_units = _infer_needed_units_numpy(ref_date, data_units)
        needed_time_delta = _unit_timedelta_numpy(needed_units)

        floor_division = np.issubdtype(dtype, np.integer) or dtype is None
        if time_delta > needed_time_delta:
            floor_division = False
            if dtype is None:
                emit_user_level_warning(
                    f"Times can't be serialized faithfully to int64 with requested units {units!r}. "
                    f"Resolution of {needed_units!r} needed. Serializing times to floating point instead. "
                    f"Set encoding['dtype'] to integer dtype to serialize to int64. "
                    f"Set encoding['dtype'] to floating point dtype to silence this warning."
                )
            elif np.issubdtype(dtype, np.integer) and allow_units_modification:
                new_units = f"{needed_units} since {format_timestamp(ref_date)}"
                emit_user_level_warning(
                    f"Times can't be serialized faithfully to int64 with requested units {units!r}. "
                    f"Serializing with units {new_units!r} instead. "
                    f"Set encoding['dtype'] to floating point dtype to serialize with units {units!r}. "
                    f"Set encoding['units'] to {new_units!r} to silence this warning ."
                )
                units = new_units
                time_delta = needed_time_delta
                floor_division = True
            elif np.issubdtype(dtype, np.integer) and not allow_units_modification:
                new_units = f"{needed_units} since {format_timestamp(ref_date)}"
                raise_incompatible_units_error = True

        # get resolution of TimedeltaIndex and align time_delta
        # todo: check, if this works in any case
        num = _division(
            time_deltas, time_delta.astype(f"=m8[{time_deltas.unit}]"), floor_division
        )
        num = reshape(num.values, dates.shape)

    except (OutOfBoundsDatetime, OverflowError, ValueError):
        time_units, ref_date = _unpack_time_units_and_ref_date_cftime(units, calendar)
        time_delta_cftime = _time_units_to_timedelta(time_units)
        needed_units = _infer_needed_units_cftime(ref_date, data_units, calendar)
        needed_time_delta_cftime = _time_units_to_timedelta(needed_units)

        if (
            np.issubdtype(dtype, np.integer)
            and time_delta_cftime > needed_time_delta_cftime
        ):
            new_units = f"{needed_units} since {format_cftime_datetime(ref_date)}"
            if allow_units_modification:
                emit_user_level_warning(
                    f"Times can't be serialized faithfully to int64 with requested units {units!r}. "
                    f"Serializing with units {new_units!r} instead. "
                    f"Set encoding['dtype'] to floating point dtype to serialize with units {units!r}. "
                    f"Set encoding['units'] to {new_units!r} to silence this warning ."
                )
                units = new_units
            else:
                raise_incompatible_units_error = True

        num = _encode_datetime_with_cftime(dates, units, calendar)
        # do it now only for cftime-based flow
        # we already covered for this in pandas-based flow
        num = cast_to_int_if_safe(num)

    if raise_incompatible_units_error:
        raise ValueError(
            f"Times can't be serialized faithfully to int64 with requested units {units!r}. "
            f"Consider setting encoding['dtype'] to a floating point dtype to serialize with "
            f"units {units!r}. Consider setting encoding['units'] to {new_units!r} to "
            f"serialize with an integer dtype."
        )
    if raise_gregorian_proleptic_gregorian_mismatch_error:
        raise ValueError(
            f"Unable to encode np.datetime64 values with {calendar} "
            f"calendar, because some or all values are prior to the reform "
            f"date of 1582-10-15. To encode these times, set "
            f"encoding['calendar'] to 'proleptic_gregorian' instead, which "
            f"is the true calendar that np.datetime64 values use. The "
            f"'standard' or 'gregorian' calendar is only equivalent to the "
            f"'proleptic_gregorian' calendar after the reform date."
        )

    return num, units, calendar


def _encode_cf_datetime_within_map_blocks(
    dates: T_DuckArray,  # type: ignore[misc]
    units: str,
    calendar: str,
    dtype: np.dtype,
) -> T_DuckArray:
    num, *_ = _eagerly_encode_cf_datetime(
        dates, units, calendar, dtype, allow_units_modification=False
    )
    return num


def _lazily_encode_cf_datetime(
    dates: T_ChunkedArray,
    units: str | None = None,
    calendar: str | None = None,
    dtype: np.dtype | None = None,
) -> tuple[T_ChunkedArray, str, str]:
    if calendar is None:
        # This will only trigger minor compute if dates is an object dtype array.
        calendar = infer_calendar_name(dates)

    if units is None and dtype is None:
        if dates.dtype == "O":
            units = "microseconds since 1970-01-01"
            dtype = np.dtype("int64")
        else:
            netcdf_unit = _numpy_dtype_to_netcdf_timeunit(dates.dtype)
            units = f"{netcdf_unit} since 1970-01-01"
            dtype = np.dtype("int64")

    if units is None or dtype is None:
        raise ValueError(
            f"When encoding chunked arrays of datetime values, both the units "
            f"and dtype must be prescribed or both must be unprescribed. "
            f"Prescribing only one or the other is not currently supported. "
            f"Got a units encoding of {units} and a dtype encoding of {dtype}."
        )

    chunkmanager = get_chunked_array_type(dates)
    num = chunkmanager.map_blocks(
        _encode_cf_datetime_within_map_blocks,
        dates,
        units,
        calendar,
        dtype,
        dtype=dtype,
    )
    return num, units, calendar


def encode_cf_timedelta(
    timedeltas: T_DuckArray,  # type: ignore[misc]
    units: str | None = None,
    dtype: np.dtype | None = None,
) -> tuple[T_DuckArray, str]:
    timedeltas = asarray(timedeltas)
    if is_chunked_array(timedeltas):
        return _lazily_encode_cf_timedelta(timedeltas, units, dtype)
    else:
        return _eagerly_encode_cf_timedelta(timedeltas, units, dtype)


def _eagerly_encode_cf_timedelta(
    timedeltas: T_DuckArray,  # type: ignore[misc]
    units: str | None = None,
    dtype: np.dtype | None = None,
    allow_units_modification: bool = True,
) -> tuple[T_DuckArray, str]:
    data_units = infer_timedelta_units(timedeltas)
    if units is None:
        units = data_units
    # units take precedence in the case of zero-size array
    if timedeltas.size == 0:
        data_units = units

    time_delta = _unit_timedelta_numpy(units)
    time_deltas = pd.TimedeltaIndex(ravel(timedeltas))
    # get resolution of TimedeltaIndex and align time_delta
    deltas_unit = time_deltas.unit
    time_delta = time_delta.astype(f"=m8[{deltas_unit}]")

    # retrieve needed units to faithfully encode to int64
    needed_units = data_units
    if data_units != units:
        needed_units = _infer_time_units_from_diff(np.unique(time_deltas.dropna()))

    # needed time delta to encode faithfully to int64
    needed_time_delta = _unit_timedelta_numpy(needed_units)

    floor_division = np.issubdtype(dtype, np.integer) or dtype is None
    if time_delta > needed_time_delta:
        floor_division = False
        if dtype is None:
            emit_user_level_warning(
                f"Timedeltas can't be serialized faithfully to int64 with requested units {units!r}. "
                f"Resolution of {needed_units!r} needed. Serializing timeseries to floating point instead. "
                f"Set encoding['dtype'] to integer dtype to serialize to int64. "
                f"Set encoding['dtype'] to floating point dtype to silence this warning."
            )
        elif np.issubdtype(dtype, np.integer) and allow_units_modification:
            emit_user_level_warning(
                f"Timedeltas can't be serialized faithfully with requested units {units!r}. "
                f"Serializing with units {needed_units!r} instead. "
                f"Set encoding['dtype'] to floating point dtype to serialize with units {units!r}. "
                f"Set encoding['units'] to {needed_units!r} to silence this warning ."
            )
            units = needed_units
            time_delta = needed_time_delta
            time_delta = time_delta.astype(f"=m8[{deltas_unit}]")
            floor_division = True
        elif np.issubdtype(dtype, np.integer) and not allow_units_modification:
            raise ValueError(
                f"Timedeltas can't be serialized faithfully to int64 with requested units {units!r}. "
                f"Consider setting encoding['dtype'] to a floating point dtype to serialize with "
                f"units {units!r}. Consider setting encoding['units'] to {needed_units!r} to "
                f"serialize with an integer dtype."
            )

    num = _division(time_deltas, time_delta, floor_division)
    num = reshape(num.values, timedeltas.shape)

    return num, units


def _encode_cf_timedelta_within_map_blocks(
    timedeltas: T_DuckArray,  # type: ignore[misc]
    units: str,
    dtype: np.dtype,
) -> T_DuckArray:
    num, _ = _eagerly_encode_cf_timedelta(
        timedeltas, units, dtype, allow_units_modification=False
    )
    return num


def _lazily_encode_cf_timedelta(
    timedeltas: T_ChunkedArray, units: str | None = None, dtype: np.dtype | None = None
) -> tuple[T_ChunkedArray, str]:
    if units is None and dtype is None:
        units = _numpy_dtype_to_netcdf_timeunit(timedeltas.dtype)
        dtype = np.dtype("int64")

    if units is None or dtype is None:
        raise ValueError(
            f"When encoding chunked arrays of timedelta values, both the "
            f"units and dtype must be prescribed or both must be "
            f"unprescribed. Prescribing only one or the other is not "
            f"currently supported. Got a units encoding of {units} and a "
            f"dtype encoding of {dtype}."
        )

    chunkmanager = get_chunked_array_type(timedeltas)
    num = chunkmanager.map_blocks(
        _encode_cf_timedelta_within_map_blocks,
        timedeltas,
        units,
        dtype,
        dtype=dtype,
    )
    return num, units


class CFDatetimeCoder(VariableCoder):
    """Coder for CF Datetime coding.

    Parameters
    ----------
    use_cftime : bool, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64`` objects.  If False, always
        decode times to ``np.datetime64`` objects; if this is not possible
        raise an error.
        May not be supported by all the backends.
    time_unit : PDDatetimeUnitOptions
          Target resolution when decoding dates. Defaults to "ns".
    """

    def __init__(
        self,
        use_cftime: bool | None = None,
        time_unit: PDDatetimeUnitOptions = "ns",
    ) -> None:
        self.use_cftime = use_cftime
        self.time_unit = time_unit

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        if np.issubdtype(
            variable.data.dtype, np.datetime64
        ) or contains_cftime_datetimes(variable):
            dims, data, attrs, encoding = unpack_for_encoding(variable)

            units = encoding.pop("units", None)
            calendar = encoding.pop("calendar", None)
            dtype = encoding.get("dtype", None)

            # in the case of packed data we need to encode into
            # float first, the correct dtype will be established
            # via CFScaleOffsetCoder/CFMaskCoder
            if "add_offset" in encoding or "scale_factor" in encoding:
                dtype = data.dtype if data.dtype.kind == "f" else "float64"
            (data, units, calendar) = encode_cf_datetime(data, units, calendar, dtype)

            safe_setitem(attrs, "units", units, name=name)
            safe_setitem(attrs, "calendar", calendar, name=name)

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        units = variable.attrs.get("units", None)
        if isinstance(units, str) and "since" in units:
            dims, data, attrs, encoding = unpack_for_decoding(variable)

            units = pop_to(attrs, encoding, "units")
            calendar = pop_to(attrs, encoding, "calendar")
            dtype = _decode_cf_datetime_dtype(
                data, units, calendar, self.use_cftime, self.time_unit
            )
            transform = partial(
                decode_cf_datetime,
                units=units,
                calendar=calendar,
                use_cftime=self.use_cftime,
                time_unit=self.time_unit,
            )
            data = lazy_elemwise_func(data, transform, dtype)

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable


def has_timedelta64_encoding_dtype(attrs_or_encoding: dict) -> bool:
    dtype = attrs_or_encoding.get("dtype")
    return isinstance(dtype, str) and dtype.startswith("timedelta64")


def resolve_time_unit_from_attrs_dtype(
    attrs_dtype: str, name: T_Name
) -> PDDatetimeUnitOptions:
    dtype = np.dtype(attrs_dtype)
    resolution, _ = np.datetime_data(dtype)
    resolution = cast(NPDatetimeUnitOptions, resolution)
    time_unit: PDDatetimeUnitOptions
    if np.timedelta64(1, resolution) > np.timedelta64(1, "s"):
        time_unit = "s"
        message = (
            f"Following pandas, xarray only supports decoding to timedelta64 "
            f"values with a resolution of 's', 'ms', 'us', or 'ns'. Encoded "
            f"values for variable {name!r} have a resolution of "
            f"{resolution!r}. Attempting to decode to a resolution of 's'. "
            f"Note, depending on the encoded values, this may lead to an "
            f"OverflowError. Additionally, data will not be identically round "
            f"tripped; xarray will choose an encoding dtype of "
            f"'timedelta64[s]' when re-encoding."
        )
        emit_user_level_warning(message)
    elif np.timedelta64(1, resolution) < np.timedelta64(1, "ns"):
        time_unit = "ns"
        message = (
            f"Following pandas, xarray only supports decoding to timedelta64 "
            f"values with a resolution of 's', 'ms', 'us', or 'ns'. Encoded "
            f"values for variable {name!r} have a resolution of "
            f"{resolution!r}. Attempting to decode to a resolution of 'ns'. "
            f"Note, depending on the encoded values, this may lead to loss of "
            f"precision. Additionally, data will not be identically round "
            f"tripped; xarray will choose an encoding dtype of "
            f"'timedelta64[ns]' when re-encoding."
        )
        emit_user_level_warning(message)
    else:
        time_unit = cast(PDDatetimeUnitOptions, resolution)
    return time_unit


class CFTimedeltaCoder(VariableCoder):
    """Coder for CF Timedelta coding.

    Parameters
    ----------
    time_unit : PDDatetimeUnitOptions
        Target resolution when decoding timedeltas via units. Defaults to "ns".
        When decoding via dtype, the resolution is specified in the dtype
        attribute, so this parameter is ignored.
    decode_via_units : bool
        Whether to decode timedeltas based on the presence of a timedelta-like
        units attribute, e.g. "seconds". Defaults to True, but in the future
        will default to False.
    decode_via_dtype : bool
        Whether to decode timedeltas based on the presence of a np.timedelta64
        dtype attribute, e.g. "timedelta64[s]". Defaults to True.
    """

    def __init__(
        self,
        time_unit: PDDatetimeUnitOptions | None = None,
        decode_via_units: bool = True,
        decode_via_dtype: bool = True,
    ) -> None:
        self.time_unit = time_unit
        self.decode_via_units = decode_via_units
        self.decode_via_dtype = decode_via_dtype
        self._emit_decode_timedelta_future_warning = False

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        if np.issubdtype(variable.data.dtype, np.timedelta64):
            dims, data, attrs, encoding = unpack_for_encoding(variable)
            dtype = encoding.get("dtype", None)
            units = encoding.pop("units", None)

            # in the case of packed data we need to encode into
            # float first, the correct dtype will be established
            # via CFScaleOffsetCoder/CFMaskCoder
            if "add_offset" in encoding or "scale_factor" in encoding:
                dtype = data.dtype if data.dtype.kind == "f" else "float64"

            resolution, _ = np.datetime_data(variable.dtype)
            attrs_dtype = f"timedelta64[{resolution}]"
            safe_setitem(attrs, "dtype", attrs_dtype, name=name)

            data, units = encode_cf_timedelta(data, units, dtype)
            safe_setitem(attrs, "units", units, name=name)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        units = variable.attrs.get("units", None)
        has_timedelta_units = isinstance(units, str) and units in TIME_UNITS
        has_timedelta_dtype = has_timedelta64_encoding_dtype(variable.attrs)
        is_dtype_decodable = has_timedelta_units and has_timedelta_dtype
        is_units_decodable = has_timedelta_units
        if (is_dtype_decodable and self.decode_via_dtype) or (
            is_units_decodable and self.decode_via_units
        ):
            dims, data, attrs, encoding = unpack_for_decoding(variable)
            units = pop_to(attrs, encoding, "units")
            if is_dtype_decodable:
                attrs_dtype = attrs.pop("dtype")
                if self.time_unit is None:
                    time_unit = resolve_time_unit_from_attrs_dtype(attrs_dtype, name)
                else:
                    time_unit = self.time_unit
            else:
                if self._emit_decode_timedelta_future_warning:
                    var_string = f"the variable {name!r}" if name else ""
                    emit_user_level_warning(
                        "In a future version, xarray will not decode "
                        f"{var_string} into a timedelta64 dtype based on the "
                        "presence of a timedelta-like 'units' attribute by "
                        "default. Instead it will rely on the presence of a "
                        "timedelta64 'dtype' attribute, which is now xarray's "
                        "default way of encoding timedelta64 values.\n"
                        "To continue decoding into a timedelta64 dtype, either "
                        "set `decode_timedelta=True` when opening this "
                        "dataset, or add the attribute "
                        "`dtype='timedelta64[ns]'` to this variable on disk.\n"
                        "To opt-in to future behavior, set "
                        "`decode_timedelta=False`.",
                        FutureWarning,
                    )
                if self.time_unit is None:
                    time_unit = "ns"
                else:
                    time_unit = self.time_unit

                # Handle edge case that decode_via_dtype=False and
                # decode_via_units=True, and timedeltas were encoded with a
                # dtype attribute. We need to remove the dtype attribute
                # to prevent an error during round tripping.
                if has_timedelta_dtype:
                    attrs.pop("dtype")

            dtype = np.dtype(f"timedelta64[{time_unit}]")
            transform = partial(decode_cf_timedelta, units=units, time_unit=time_unit)
            data = lazy_elemwise_func(data, transform, dtype=dtype)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable
