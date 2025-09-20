"""DatetimeIndex analog for cftime.datetime objects"""

# The pandas.Index subclass defined here was copied and adapted for
# use with cftime.datetime objects based on the source code defining
# pandas.DatetimeIndex.

# For reference, here is a copy of the pandas copyright notice:

# (c) 2011-2012, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2008-2011 AQR Capital Management, LLC
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the copyright holder nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from __future__ import annotations

import math
from datetime import timedelta
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from packaging.version import Version

from xarray.coding.times import (
    _STANDARD_CALENDARS,
    _parse_iso8601,
    cftime_to_nptime,
    infer_calendar_name,
)
from xarray.core.common import _contains_cftime_datetimes
from xarray.core.options import OPTIONS
from xarray.core.types import PDDatetimeUnitOptions
from xarray.core.utils import attempt_import, emit_user_level_warning, is_scalar

if TYPE_CHECKING:
    from xarray.coding.cftime_offsets import BaseCFTimeOffset
    from xarray.core.types import Self


# constants for cftimeindex.repr
CFTIME_REPR_LENGTH = 19
ITEMS_IN_REPR_MAX_ELSE_ELLIPSIS = 100
REPR_ELLIPSIS_SHOW_ITEMS_FRONT_END = 10


OUT_OF_BOUNDS_TIMEDELTA_ERRORS: tuple[type[Exception], ...]
try:
    OUT_OF_BOUNDS_TIMEDELTA_ERRORS = (pd.errors.OutOfBoundsTimedelta, OverflowError)
except AttributeError:
    OUT_OF_BOUNDS_TIMEDELTA_ERRORS = (OverflowError,)


def _parsed_string_to_bounds(date_type, resolution, parsed):
    """Generalization of
    pandas.tseries.index.DatetimeIndex._parsed_string_to_bounds
    for use with non-standard calendars and cftime.datetime
    objects.
    """
    if resolution == "year":
        return (
            date_type(parsed.year, 1, 1),
            date_type(parsed.year + 1, 1, 1) - timedelta(microseconds=1),
        )
    elif resolution == "month":
        if parsed.month == 12:
            end = date_type(parsed.year + 1, 1, 1) - timedelta(microseconds=1)
        else:
            end = date_type(parsed.year, parsed.month + 1, 1) - timedelta(
                microseconds=1
            )
        return date_type(parsed.year, parsed.month, 1), end
    elif resolution == "day":
        start = date_type(parsed.year, parsed.month, parsed.day)
        return start, start + timedelta(days=1, microseconds=-1)
    elif resolution == "hour":
        start = date_type(parsed.year, parsed.month, parsed.day, parsed.hour)
        return start, start + timedelta(hours=1, microseconds=-1)
    elif resolution == "minute":
        start = date_type(
            parsed.year, parsed.month, parsed.day, parsed.hour, parsed.minute
        )
        return start, start + timedelta(minutes=1, microseconds=-1)
    elif resolution == "second":
        start = date_type(
            parsed.year,
            parsed.month,
            parsed.day,
            parsed.hour,
            parsed.minute,
            parsed.second,
        )
        return start, start + timedelta(seconds=1, microseconds=-1)
    else:
        raise KeyError


def get_date_field(datetimes, field):
    """Adapted from pandas.tslib.get_date_field"""
    return np.array([getattr(date, field) for date in datetimes], dtype=np.int64)


def _field_accessor(name, docstring=None, min_cftime_version="0.0"):
    """Adapted from pandas.tseries.index._field_accessor"""

    def f(self, min_cftime_version=min_cftime_version):
        if TYPE_CHECKING:
            import cftime
        else:
            cftime = attempt_import("cftime")

        if Version(cftime.__version__) >= Version(min_cftime_version):
            return get_date_field(self._data, name)
        else:
            raise ImportError(
                f"The {name:!r} accessor requires a minimum "
                f"version of cftime of {min_cftime_version}. Found an "
                f"installed version of {cftime.__version__}."
            )

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def get_date_type(self):
    if self._data.size:
        return type(self._data[0])
    else:
        return None


def assert_all_valid_date_type(data):
    if TYPE_CHECKING:
        import cftime
    else:
        cftime = attempt_import("cftime")

    if len(data) > 0:
        sample = data[0]
        date_type = type(sample)
        if not isinstance(sample, cftime.datetime):
            raise TypeError(
                "CFTimeIndex requires cftime.datetime "
                f"objects. Got object of {date_type}."
            )
        if not all(isinstance(value, date_type) for value in data):
            raise TypeError(
                "CFTimeIndex requires using datetime "
                f"objects of all the same type.  Got\n{data}."
            )


def format_row(times, indent=0, separator=", ", row_end=",\n"):
    """Format a single row from format_times."""
    return indent * " " + separator.join(map(str, times)) + row_end


def format_times(
    index,
    max_width,
    offset,
    separator=", ",
    first_row_offset=0,
    intermediate_row_end=",\n",
    last_row_end="",
):
    """Format values of cftimeindex as pd.Index."""
    n_per_row = max(max_width // (CFTIME_REPR_LENGTH + len(separator)), 1)
    n_rows = math.ceil(len(index) / n_per_row)

    representation = ""
    for row in range(n_rows):
        indent = first_row_offset if row == 0 else offset
        row_end = last_row_end if row == n_rows - 1 else intermediate_row_end
        times_for_row = index[row * n_per_row : (row + 1) * n_per_row]
        representation += format_row(
            times_for_row, indent=indent, separator=separator, row_end=row_end
        )

    return representation


def format_attrs(index, separator=", "):
    """Format attributes of CFTimeIndex for __repr__."""
    attrs = {
        "dtype": f"'{index.dtype}'",
        "length": f"{len(index)}",
        "calendar": f"{index.calendar!r}",
        "freq": f"{index.freq!r}",
    }

    attrs_str = [f"{k}={v}" for k, v in attrs.items()]
    attrs_str = f"{separator}".join(attrs_str)
    return attrs_str


class CFTimeIndex(pd.Index):
    """Custom Index for working with CF calendars and dates

    All elements of a CFTimeIndex must be cftime.datetime objects.

    Parameters
    ----------
    data : array or CFTimeIndex
        Sequence of cftime.datetime objects to use in index
    name : str, default: None
        Name of the resulting index

    See Also
    --------
    date_range
    """

    _data: np.ndarray

    year = _field_accessor("year", "The year of the datetime")
    month = _field_accessor("month", "The month of the datetime")
    day = _field_accessor("day", "The days of the datetime")
    hour = _field_accessor("hour", "The hours of the datetime")
    minute = _field_accessor("minute", "The minutes of the datetime")
    second = _field_accessor("second", "The seconds of the datetime")
    microsecond = _field_accessor("microsecond", "The microseconds of the datetime")
    dayofyear = _field_accessor(
        "dayofyr", "The ordinal day of year of the datetime", "1.0.2.1"
    )
    dayofweek = _field_accessor("dayofwk", "The day of week of the datetime", "1.0.2.1")
    days_in_month = _field_accessor(
        "daysinmonth", "The number of days in the month of the datetime", "1.1.0.0"
    )
    date_type = property(get_date_type)

    def __new__(cls, data, name=None, **kwargs):
        assert_all_valid_date_type(data)
        if name is None and hasattr(data, "name"):
            name = data.name

        result = object.__new__(cls)
        result._data = np.array(data, dtype="O")
        result.name = name
        result._cache = {}
        return result

    def __repr__(self):
        """
        Return a string representation for this object.
        """
        klass_name = type(self).__name__
        display_width = OPTIONS["display_width"]
        offset = len(klass_name) + 2

        if len(self) <= ITEMS_IN_REPR_MAX_ELSE_ELLIPSIS:
            datastr = format_times(
                self.values, display_width, offset=offset, first_row_offset=0
            )
        else:
            front_str = format_times(
                self.values[:REPR_ELLIPSIS_SHOW_ITEMS_FRONT_END],
                display_width,
                offset=offset,
                first_row_offset=0,
                last_row_end=",",
            )
            end_str = format_times(
                self.values[-REPR_ELLIPSIS_SHOW_ITEMS_FRONT_END:],
                display_width,
                offset=offset,
                first_row_offset=offset,
            )
            datastr = "\n".join([front_str, f"{' ' * offset}...", end_str])

        attrs_str = format_attrs(self)
        # oneliner only if smaller than display_width
        full_repr_str = f"{klass_name}([{datastr}], {attrs_str})"
        if len(full_repr_str) > display_width:
            # if attrs_str too long, one per line
            if len(attrs_str) >= display_width - offset:
                attrs_str = attrs_str.replace(",", f",\n{' ' * (offset - 2)}")
            full_repr_str = (
                f"{klass_name}([{datastr}],\n{' ' * (offset - 1)}{attrs_str})"
            )

        return full_repr_str

    def _partial_date_slice(self, resolution, parsed):
        """Adapted from
        pandas.tseries.index.DatetimeIndex._partial_date_slice

        Note that when using a CFTimeIndex, if a partial-date selection
        returns a single element, it will never be converted to a scalar
        coordinate; this is in slight contrast to the behavior when using
        a DatetimeIndex, which sometimes will return a DataArray with a scalar
        coordinate depending on the resolution of the datetimes used in
        defining the index.  For example:

        >>> from cftime import DatetimeNoLeap
        >>> da = xr.DataArray(
        ...     [1, 2],
        ...     coords=[[DatetimeNoLeap(2001, 1, 1), DatetimeNoLeap(2001, 2, 1)]],
        ...     dims=["time"],
        ... )
        >>> da.sel(time="2001-01-01")
        <xarray.DataArray (time: 1)> Size: 8B
        array([1])
        Coordinates:
          * time     (time) object 8B 2001-01-01 00:00:00
        >>> da = xr.DataArray(
        ...     [1, 2],
        ...     coords=[[pd.Timestamp(2001, 1, 1), pd.Timestamp(2001, 2, 1)]],
        ...     dims=["time"],
        ... )
        >>> da.sel(time="2001-01-01")
        <xarray.DataArray ()> Size: 8B
        array(1)
        Coordinates:
            time     datetime64[ns] 8B 2001-01-01
        >>> da = xr.DataArray(
        ...     [1, 2],
        ...     coords=[[pd.Timestamp(2001, 1, 1, 1), pd.Timestamp(2001, 2, 1)]],
        ...     dims=["time"],
        ... )
        >>> da.sel(time="2001-01-01")
        <xarray.DataArray (time: 1)> Size: 8B
        array([1])
        Coordinates:
          * time     (time) datetime64[ns] 8B 2001-01-01T01:00:00
        """
        start, end = _parsed_string_to_bounds(self.date_type, resolution, parsed)

        times = self._data

        if self.is_monotonic_increasing:
            if len(times) and (
                (start < times[0] and end < times[0])
                or (start > times[-1] and end > times[-1])
            ):
                # we are out of range
                raise KeyError

            # a monotonic (sorted) series can be sliced
            left = times.searchsorted(start, side="left")
            right = times.searchsorted(end, side="right")
            return slice(left, right)

        lhs_mask = times >= start
        rhs_mask = times <= end
        return np.flatnonzero(lhs_mask & rhs_mask)

    def _get_string_slice(self, key):
        """Adapted from pandas.tseries.index.DatetimeIndex._get_string_slice"""
        parsed, resolution = _parse_iso8601(self.date_type, key)
        try:
            loc = self._partial_date_slice(resolution, parsed)
        except KeyError as err:
            raise KeyError(key) from err
        return loc

    def _get_nearest_indexer(self, target, limit, tolerance):
        """Adapted from pandas.Index._get_nearest_indexer"""
        left_indexer = self.get_indexer(target, "pad", limit=limit)
        right_indexer = self.get_indexer(target, "backfill", limit=limit)
        left_distances = abs(self.values[left_indexer] - target.values)
        right_distances = abs(self.values[right_indexer] - target.values)

        if self.is_monotonic_increasing:
            condition = (left_distances < right_distances) | (right_indexer == -1)
        else:
            condition = (left_distances <= right_distances) | (right_indexer == -1)
        indexer = np.where(condition, left_indexer, right_indexer)

        if tolerance is not None:
            indexer = self._filter_indexer_tolerance(target, indexer, tolerance)
        return indexer

    def _filter_indexer_tolerance(self, target, indexer, tolerance):
        """Adapted from pandas.Index._filter_indexer_tolerance"""
        if isinstance(target, pd.Index):
            distance = abs(self.values[indexer] - target.values)
        else:
            distance = abs(self.values[indexer] - target)
        indexer = np.where(distance <= tolerance, indexer, -1)
        return indexer

    def get_loc(self, key):
        """Adapted from pandas.tseries.index.DatetimeIndex.get_loc"""
        if isinstance(key, str):
            return self._get_string_slice(key)
        else:
            return super().get_loc(key)

    def _maybe_cast_slice_bound(self, label, side):
        """Adapted from
        pandas.tseries.index.DatetimeIndex._maybe_cast_slice_bound
        """
        if not isinstance(label, str):
            return label

        parsed, resolution = _parse_iso8601(self.date_type, label)
        start, end = _parsed_string_to_bounds(self.date_type, resolution, parsed)
        if self.is_monotonic_decreasing and len(self) > 1:
            return end if side == "left" else start
        return start if side == "left" else end

    # TODO: Add ability to use integer range outside of iloc?
    # e.g. series[1:5].
    def get_value(self, series, key):
        """Adapted from pandas.tseries.index.DatetimeIndex.get_value"""
        if np.asarray(key).dtype == np.dtype(bool):
            return series.iloc[key]
        elif isinstance(key, slice):
            return series.iloc[self.slice_indexer(key.start, key.stop, key.step)]
        else:
            return series.iloc[self.get_loc(key)]

    def __contains__(self, key: Any) -> bool:
        """Adapted from
        pandas.tseries.base.DatetimeIndexOpsMixin.__contains__"""
        try:
            result = self.get_loc(key)
            return (
                is_scalar(result)
                or isinstance(result, slice)
                or (isinstance(result, np.ndarray) and result.size > 0)
            )
        except (KeyError, TypeError, ValueError):
            return False

    def contains(self, key: Any) -> bool:
        """Needed for .loc based partial-string indexing"""
        return self.__contains__(key)

    def shift(  # type: ignore[override,unused-ignore]
        self,
        periods: int | float,
        freq: str | timedelta | BaseCFTimeOffset | None = None,
    ) -> Self:
        """Shift the CFTimeIndex a multiple of the given frequency.

        See the documentation for :py:func:`~xarray.date_range` for a
        complete listing of valid frequency strings.

        Parameters
        ----------
        periods : int, float if freq of days or below
            Periods to shift by
        freq : str, datetime.timedelta or BaseCFTimeOffset
            A frequency string or datetime.timedelta object to shift by

        Returns
        -------
        CFTimeIndex

        See Also
        --------
        pandas.DatetimeIndex.shift

        Examples
        --------
        >>> index = xr.date_range("2000", periods=1, freq="ME", use_cftime=True)
        >>> index
        CFTimeIndex([2000-01-31 00:00:00],
                    dtype='object', length=1, calendar='standard', freq=None)
        >>> index.shift(1, "ME")
        CFTimeIndex([2000-02-29 00:00:00],
                    dtype='object', length=1, calendar='standard', freq=None)
        >>> index.shift(1.5, "24h")
        CFTimeIndex([2000-02-01 12:00:00],
                    dtype='object', length=1, calendar='standard', freq=None)
        """
        from xarray.coding.cftime_offsets import BaseCFTimeOffset

        if freq is None:
            # None type is required to be compatible with base pd.Index class
            raise TypeError(
                f"`freq` argument cannot be None for {type(self).__name__}.shift"
            )

        if isinstance(freq, timedelta):
            return self + periods * freq

        if isinstance(freq, str | BaseCFTimeOffset):
            from xarray.coding.cftime_offsets import to_offset

            return self + periods * to_offset(freq)

        raise TypeError(
            f"'freq' must be of type str or datetime.timedelta, got {type(freq)}."
        )

    def __add__(self, other) -> Self:
        if isinstance(other, pd.TimedeltaIndex):
            other = other.to_pytimedelta()
        return type(self)(np.array(self) + other)

    def __radd__(self, other) -> Self:
        if isinstance(other, pd.TimedeltaIndex):
            other = other.to_pytimedelta()
        return type(self)(other + np.array(self))

    def __sub__(self, other):
        if _contains_datetime_timedeltas(other):
            return type(self)(np.array(self) - other)
        if isinstance(other, pd.TimedeltaIndex):
            return type(self)(np.array(self) - other.to_pytimedelta())
        if _contains_cftime_datetimes(np.array(other)):
            try:
                return pd.TimedeltaIndex(np.array(self) - np.array(other))
            except OUT_OF_BOUNDS_TIMEDELTA_ERRORS as err:
                raise ValueError(
                    "The time difference exceeds the range of values "
                    "that can be expressed at the nanosecond resolution."
                ) from err
        return NotImplemented

    def __rsub__(self, other):
        try:
            return pd.TimedeltaIndex(other - np.array(self))
        except OUT_OF_BOUNDS_TIMEDELTA_ERRORS as err:
            raise ValueError(
                "The time difference exceeds the range of values "
                "that can be expressed at the nanosecond resolution."
            ) from err

    def to_datetimeindex(
        self, unsafe: bool = False, time_unit: PDDatetimeUnitOptions | None = None
    ) -> pd.DatetimeIndex:
        """If possible, convert this index to a pandas.DatetimeIndex.

        Parameters
        ----------
        unsafe : bool
            Flag to turn off calendar mismatch warnings (default ``False``).
        time_unit : str
            Time resolution of resulting DatetimeIndex. Can be one of `"s"`,
            ``"ms"``, ``"us"``, or ``"ns"`` (default ``"ns"``).

        Returns
        -------
        pandas.DatetimeIndex

        Raises
        ------
        ValueError
            If the CFTimeIndex contains dates that are not possible in the
            standard calendar or outside the range representable by the
            specified ``time_unit``.

        Warns
        -----
        RuntimeWarning
            If converting from a non-standard calendar, or a Gregorian
            calendar with dates prior to the reform (1582-10-15).

        Warnings
        --------
        Note that for non-proleptic Gregorian calendars, this will change the
        calendar type of the index. In that case the result of this method
        should be used with caution.

        Examples
        --------
        >>> times = xr.date_range(
        ...     "2000", periods=2, calendar="gregorian", use_cftime=True
        ... )
        >>> times
        CFTimeIndex([2000-01-01 00:00:00, 2000-01-02 00:00:00],
                    dtype='object', length=2, calendar='standard', freq=None)
        >>> times.to_datetimeindex(time_unit="ns")
        DatetimeIndex(['2000-01-01', '2000-01-02'], dtype='datetime64[ns]', freq=None)
        """

        if not self._data.size:
            return pd.DatetimeIndex([])

        if time_unit is None:
            emit_user_level_warning(
                "In a future version of xarray to_datetimeindex will default "
                "to returning a 'us'-resolution DatetimeIndex instead of a "
                "'ns'-resolution DatetimeIndex. This warning can be silenced "
                "by explicitly passing the `time_unit` keyword argument.",
                FutureWarning,
            )
            time_unit = "ns"

        nptimes = cftime_to_nptime(self, time_unit=time_unit)
        calendar = infer_calendar_name(self)
        if calendar not in _STANDARD_CALENDARS and not unsafe:
            emit_user_level_warning(
                "Converting a CFTimeIndex with dates from a non-standard "
                f"calendar, {calendar!r}, to a pandas.DatetimeIndex, which "
                "uses dates from the standard calendar.  This may lead to "
                "subtle errors in operations that depend on the length of "
                "time between dates.",
                RuntimeWarning,
            )
        if calendar == "standard" and not unsafe:
            reform_date = self.date_type(1582, 10, 15)
            if self.min() < reform_date:
                emit_user_level_warning(
                    "Converting a CFTimeIndex with dates from a Gregorian "
                    "calendar that fall before the reform date of 1582-10-15 "
                    "to a pandas.DatetimeIndex. During this time period the "
                    "Gregorian calendar and the proleptic Gregorian calendar "
                    "of the DatetimeIndex do not exactly align. This warning "
                    "can be silenced by setting unsafe=True.",
                    RuntimeWarning,
                )

        return pd.DatetimeIndex(nptimes)

    def strftime(self, date_format):
        """
        Return an Index of formatted strings specified by date_format, which
        supports the same string format as the python standard library. Details
        of the string format can be found in `python string format doc
        <https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior>`__

        Parameters
        ----------
        date_format : str
            Date format string (e.g. "%Y-%m-%d")

        Returns
        -------
        pandas.Index
            Index of formatted strings

        Examples
        --------
        >>> rng = xr.date_range(
        ...     start="2000",
        ...     periods=5,
        ...     freq="2MS",
        ...     calendar="noleap",
        ...     use_cftime=True,
        ... )
        >>> rng.strftime("%B %d, %Y, %r")
        Index(['January 01, 2000, 12:00:00 AM', 'March 01, 2000, 12:00:00 AM',
               'May 01, 2000, 12:00:00 AM', 'July 01, 2000, 12:00:00 AM',
               'September 01, 2000, 12:00:00 AM'],
              dtype='object')
        """
        return pd.Index([date.strftime(date_format) for date in self._data])

    @property
    def asi8(self):
        """Convert to integers with units of microseconds since 1970-01-01."""
        from xarray.core.resample_cftime import exact_cftime_datetime_difference

        if not self._data.size:
            return np.array([], dtype=np.int64)

        epoch = self.date_type(1970, 1, 1)
        return np.array(
            [
                _total_microseconds(exact_cftime_datetime_difference(epoch, date))
                for date in self.values
            ],
            dtype=np.int64,
        )

    @property
    def calendar(self):
        """The calendar used by the datetimes in the index."""
        if not self._data.size:
            return None

        return infer_calendar_name(self)

    @property
    def freq(self):
        """The frequency used by the dates in the index."""
        from xarray.coding.frequencies import infer_freq

        # min 3 elemtents required to determine freq
        if self._data.size < 3:
            return None

        return infer_freq(self)

    def _round_via_method(self, freq, method):
        """Round dates using a specified method."""
        from xarray.coding.cftime_offsets import CFTIME_TICKS, Day, to_offset

        if not self._data.size:
            return CFTimeIndex(np.array(self))

        offset = to_offset(freq)
        if isinstance(offset, Day):
            # Following pandas, "In the 'round' context, Day unambiguously
            # means 24h, not calendar-day"
            offset_as_timedelta = timedelta(days=offset.n)
        elif isinstance(offset, CFTIME_TICKS):
            offset_as_timedelta = offset.as_timedelta()
        else:
            raise ValueError(f"{offset} is a non-fixed frequency")

        unit = _total_microseconds(offset_as_timedelta)
        values = self.asi8
        rounded = method(values, unit)
        return _cftimeindex_from_i8(rounded, self.date_type, self.name)

    def floor(self, freq):
        """Round dates down to fixed frequency.

        Parameters
        ----------
        freq : str
            The frequency level to round the index to.  Must be a fixed
            frequency like 'S' (second) not 'ME' (month end).  See `frequency
            aliases <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`_
            for a list of possible values.

        Returns
        -------
        CFTimeIndex
        """
        return self._round_via_method(freq, _floor_int)

    def ceil(self, freq):
        """Round dates up to fixed frequency.

        Parameters
        ----------
        freq : str
            The frequency level to round the index to.  Must be a fixed
            frequency like 'S' (second) not 'ME' (month end).  See `frequency
            aliases <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`_
            for a list of possible values.

        Returns
        -------
        CFTimeIndex
        """
        return self._round_via_method(freq, _ceil_int)

    def round(self, freq):
        """Round dates to a fixed frequency.

        Parameters
        ----------
        freq : str
            The frequency level to round the index to.  Must be a fixed
            frequency like 'S' (second) not 'ME' (month end).  See `frequency
            aliases <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`_
            for a list of possible values.

        Returns
        -------
        CFTimeIndex
        """
        return self._round_via_method(freq, _round_to_nearest_half_even)

    @property
    def is_leap_year(self):
        if TYPE_CHECKING:
            import cftime
        else:
            cftime = attempt_import("cftime")
        func = np.vectorize(cftime.is_leap_year)
        return func(self.year, calendar=self.calendar)


def _parse_array_of_cftime_strings(strings, date_type):
    """Create a numpy array from an array of strings.

    For use in generating dates from strings for use with interp.  Assumes the
    array is either 0-dimensional or 1-dimensional.

    Parameters
    ----------
    strings : array of strings
        Strings to convert to dates
    date_type : cftime.datetime type
        Calendar type to use for dates

    Returns
    -------
    np.array
    """
    return np.array([_parse_iso8601(date_type, s)[0] for s in strings.ravel()]).reshape(
        strings.shape
    )


def _contains_datetime_timedeltas(array):
    """Check if an input array contains datetime.timedelta objects."""
    array = np.atleast_1d(array)
    return isinstance(array[0], timedelta)


def _cftimeindex_from_i8(values, date_type, name):
    """Construct a CFTimeIndex from an array of integers.

    Parameters
    ----------
    values : np.array
        Integers representing microseconds since 1970-01-01.
    date_type : cftime.datetime
        Type of date for the index.
    name : str
        Name of the index.

    Returns
    -------
    CFTimeIndex
    """
    epoch = date_type(1970, 1, 1)
    dates = np.array([epoch + timedelta(microseconds=int(value)) for value in values])
    return CFTimeIndex(dates, name=name)


def _total_microseconds(delta):
    """Compute the total number of microseconds of a datetime.timedelta.

    Parameters
    ----------
    delta : datetime.timedelta
        Input timedelta.

    Returns
    -------
    int
    """
    return delta / timedelta(microseconds=1)


def _floor_int(values, unit):
    """Copied from pandas."""
    return values - np.remainder(values, unit)


def _ceil_int(values, unit):
    """Copied from pandas."""
    return values + np.remainder(-values, unit)


def _round_to_nearest_half_even(values, unit):
    """Copied from pandas."""
    if unit % 2:
        return _ceil_int(values - unit // 2, unit)
    quotient, remainder = np.divmod(values, unit)
    mask = np.logical_or(
        remainder > (unit // 2), np.logical_and(remainder == (unit // 2), quotient % 2)
    )
    quotient[mask] += 1
    return quotient * unit
