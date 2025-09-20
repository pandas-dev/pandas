from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Generic

import numpy as np
import pandas as pd

from xarray.coding.calendar_ops import _decimal_year
from xarray.coding.times import infer_calendar_name
from xarray.core import duck_array_ops
from xarray.core.common import (
    _contains_datetime_like_objects,
    full_like,
    is_np_datetime_like,
    is_np_timedelta_like,
)
from xarray.core.types import T_DataArray
from xarray.core.variable import IndexVariable, Variable
from xarray.namedarray.utils import is_duck_dask_array

if TYPE_CHECKING:
    from typing import Self

    from numpy.typing import DTypeLike

    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import CFCalendar


def _season_from_months(months):
    """Compute season (DJF, MAM, JJA, SON) from month ordinal"""
    # TODO: Move "season" accessor upstream into pandas
    seasons = np.array(["DJF", "MAM", "JJA", "SON", "nan"])
    months = np.asarray(months)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="invalid value encountered in floor_divide"
        )
        warnings.filterwarnings(
            "ignore", message="invalid value encountered in remainder"
        )
        idx = (months // 3) % 4

    idx[np.isnan(idx)] = 4
    return seasons[idx.astype(int)]


def _access_through_cftimeindex(values, name):
    """Coerce an array of datetime-like values to a CFTimeIndex
    and access requested datetime component
    """
    from xarray.coding.cftimeindex import CFTimeIndex

    if not isinstance(values, CFTimeIndex):
        values_as_cftimeindex = CFTimeIndex(duck_array_ops.ravel(values))
    else:
        values_as_cftimeindex = values
    if name == "season":
        months = values_as_cftimeindex.month
        field_values = _season_from_months(months)
    elif name == "date":
        raise AttributeError(
            "'CFTimeIndex' object has no attribute `date`. Consider using the floor method "
            "instead, for instance: `.time.dt.floor('D')`."
        )
    else:
        field_values = getattr(values_as_cftimeindex, name)
    return field_values.reshape(values.shape)


def _access_through_series(values, name):
    """Coerce an array of datetime-like values to a pandas Series and
    access requested datetime component
    """
    values_as_series = pd.Series(duck_array_ops.ravel(values), copy=False)
    if name == "season":
        months = values_as_series.dt.month.values
        field_values = _season_from_months(months)
    elif name == "total_seconds":
        field_values = values_as_series.dt.total_seconds().values
    elif name == "isocalendar":
        # special NaT-handling can be removed when
        # https://github.com/pandas-dev/pandas/issues/54657 is resolved
        field_values = values_as_series.dt.isocalendar()
        # test for <NA> and apply needed dtype
        hasna = any(field_values.year.isnull())
        if hasna:
            field_values = np.dstack(
                [
                    getattr(field_values, name).astype(np.float64, copy=False).values
                    for name in ["year", "week", "day"]
                ]
            )
        else:
            field_values = np.array(field_values, dtype=np.int64)
        # isocalendar returns iso- year, week, and weekday -> reshape
        return field_values.T.reshape(3, *values.shape)
    else:
        field_values = getattr(values_as_series.dt, name).values

    return field_values.reshape(values.shape)


def _get_date_field(values, name, dtype):
    """Indirectly access pandas' libts.get_date_field by wrapping data
    as a Series and calling through `.dt` attribute.

    Parameters
    ----------
    values : np.ndarray or dask.array-like
        Array-like container of datetime-like values
    name : str
        Name of datetime field to access
    dtype : dtype-like
        dtype for output date field values

    Returns
    -------
    datetime_fields : same type as values
        Array-like of datetime fields accessed for each element in values

    """
    if is_np_datetime_like(values.dtype):
        access_method = _access_through_series
    else:
        access_method = _access_through_cftimeindex

    if is_duck_dask_array(values):
        from dask.array import map_blocks

        new_axis = chunks = None
        # isocalendar adds an axis
        if name == "isocalendar":
            chunks = (3,) + values.chunksize
            new_axis = 0

        return map_blocks(
            access_method, values, name, dtype=dtype, new_axis=new_axis, chunks=chunks
        )
    else:
        out = access_method(values, name)
        # cast only for integer types to keep float64 in presence of NaT
        # see https://github.com/pydata/xarray/issues/7928
        if np.issubdtype(out.dtype, np.integer):
            out = out.astype(dtype, copy=False)
        return out


def _round_through_series_or_index(values, name, freq):
    """Coerce an array of datetime-like values to a pandas Series or xarray
    CFTimeIndex and apply requested rounding
    """
    from xarray.coding.cftimeindex import CFTimeIndex

    if is_np_datetime_like(values.dtype):
        values_as_series = pd.Series(duck_array_ops.ravel(values), copy=False)
        method = getattr(values_as_series.dt, name)
    else:
        values_as_cftimeindex = CFTimeIndex(duck_array_ops.ravel(values))
        method = getattr(values_as_cftimeindex, name)

    field_values = method(freq=freq).values

    return field_values.reshape(values.shape)


def _round_field(values, name, freq):
    """Indirectly access rounding functions by wrapping data
    as a Series or CFTimeIndex

    Parameters
    ----------
    values : np.ndarray or dask.array-like
        Array-like container of datetime-like values
    name : {"ceil", "floor", "round"}
        Name of rounding function
    freq : str
        a freq string indicating the rounding resolution

    Returns
    -------
    rounded timestamps : same type as values
        Array-like of datetime fields accessed for each element in values

    """
    if is_duck_dask_array(values):
        from dask.array import map_blocks

        dtype = np.datetime64 if is_np_datetime_like(values.dtype) else np.dtype("O")
        return map_blocks(
            _round_through_series_or_index, values, name, freq=freq, dtype=dtype
        )
    else:
        return _round_through_series_or_index(values, name, freq)


def _strftime_through_cftimeindex(values, date_format: str):
    """Coerce an array of cftime-like values to a CFTimeIndex
    and access requested datetime component
    """
    from xarray.coding.cftimeindex import CFTimeIndex

    values_as_cftimeindex = CFTimeIndex(duck_array_ops.ravel(values))

    field_values = values_as_cftimeindex.strftime(date_format)
    return field_values.to_numpy().reshape(values.shape)


def _strftime_through_series(values, date_format: str):
    """Coerce an array of datetime-like values to a pandas Series and
    apply string formatting
    """
    values_as_series = pd.Series(duck_array_ops.ravel(values), copy=False)
    strs = values_as_series.dt.strftime(date_format)
    return strs.to_numpy().reshape(values.shape)


def _strftime(values, date_format):
    if is_np_datetime_like(values.dtype):
        access_method = _strftime_through_series
    else:
        access_method = _strftime_through_cftimeindex
    if is_duck_dask_array(values):
        from dask.array import map_blocks

        return map_blocks(access_method, values, date_format)
    else:
        return access_method(values, date_format)


def _index_or_data(obj):
    if isinstance(obj.variable, IndexVariable):
        return obj.to_index()
    else:
        return obj.data


class TimeAccessor(Generic[T_DataArray]):
    __slots__ = ("_obj",)

    def __init__(self, obj: T_DataArray) -> None:
        self._obj = obj

    def _date_field(self, name: str, dtype: DTypeLike) -> T_DataArray:
        if dtype is None:
            dtype = self._obj.dtype
        result = _get_date_field(_index_or_data(self._obj), name, dtype)
        newvar = Variable(
            dims=self._obj.dims,
            attrs=self._obj.attrs,
            encoding=self._obj.encoding,
            data=result,
        )
        return self._obj._replace(newvar, name=name)

    def _tslib_round_accessor(self, name: str, freq: str) -> T_DataArray:
        result = _round_field(_index_or_data(self._obj), name, freq)
        newvar = Variable(
            dims=self._obj.dims,
            attrs=self._obj.attrs,
            encoding=self._obj.encoding,
            data=result,
        )
        return self._obj._replace(newvar, name=name)

    def floor(self, freq: str) -> T_DataArray:
        """
        Round timestamps downward to specified frequency resolution.

        Parameters
        ----------
        freq : str
            a freq string indicating the rounding resolution e.g. "D" for daily resolution

        Returns
        -------
        floor-ed timestamps : same type as values
            Array-like of datetime fields accessed for each element in values
        """

        return self._tslib_round_accessor("floor", freq)

    def ceil(self, freq: str) -> T_DataArray:
        """
        Round timestamps upward to specified frequency resolution.

        Parameters
        ----------
        freq : str
            a freq string indicating the rounding resolution e.g. "D" for daily resolution

        Returns
        -------
        ceil-ed timestamps : same type as values
            Array-like of datetime fields accessed for each element in values
        """
        return self._tslib_round_accessor("ceil", freq)

    def round(self, freq: str) -> T_DataArray:
        """
        Round timestamps to specified frequency resolution.

        Parameters
        ----------
        freq : str
            a freq string indicating the rounding resolution e.g. "D" for daily resolution

        Returns
        -------
        rounded timestamps : same type as values
            Array-like of datetime fields accessed for each element in values
        """
        return self._tslib_round_accessor("round", freq)


class DatetimeAccessor(TimeAccessor[T_DataArray]):
    """Access datetime fields for DataArrays with datetime-like dtypes.

    Fields can be accessed through the `.dt` attribute
    for applicable DataArrays.

    Examples
    ---------
    >>> dates = pd.date_range(start="2000/01/01", freq="D", periods=10)
    >>> ts = xr.DataArray(dates, dims=("time"))
    >>> ts
    <xarray.DataArray (time: 10)> Size: 80B
    array(['2000-01-01T00:00:00.000000000', '2000-01-02T00:00:00.000000000',
           '2000-01-03T00:00:00.000000000', '2000-01-04T00:00:00.000000000',
           '2000-01-05T00:00:00.000000000', '2000-01-06T00:00:00.000000000',
           '2000-01-07T00:00:00.000000000', '2000-01-08T00:00:00.000000000',
           '2000-01-09T00:00:00.000000000', '2000-01-10T00:00:00.000000000'],
          dtype='datetime64[ns]')
    Coordinates:
      * time     (time) datetime64[ns] 80B 2000-01-01 2000-01-02 ... 2000-01-10
    >>> ts.dt  # doctest: +ELLIPSIS
    <xarray.core.accessor_dt.DatetimeAccessor object at 0x...>
    >>> ts.dt.dayofyear
    <xarray.DataArray 'dayofyear' (time: 10)> Size: 80B
    array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10])
    Coordinates:
      * time     (time) datetime64[ns] 80B 2000-01-01 2000-01-02 ... 2000-01-10
    >>> ts.dt.quarter
    <xarray.DataArray 'quarter' (time: 10)> Size: 80B
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    Coordinates:
      * time     (time) datetime64[ns] 80B 2000-01-01 2000-01-02 ... 2000-01-10

    """

    def strftime(self, date_format: str) -> T_DataArray:
        """
        Return an array of formatted strings specified by date_format, which
        supports the same string format as the python standard library. Details
        of the string format can be found in `python string format doc
        <https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior>`__

        Parameters
        ----------
        date_format : str
            date format string (e.g. "%Y-%m-%d")

        Returns
        -------
        formatted strings : same type as values
            Array-like of strings formatted for each element in values

        Examples
        --------
        >>> import datetime
        >>> rng = xr.Dataset({"time": datetime.datetime(2000, 1, 1)})
        >>> rng["time"].dt.strftime("%B %d, %Y, %r")
        <xarray.DataArray 'strftime' ()> Size: 8B
        array('January 01, 2000, 12:00:00 AM', dtype=object)
        """
        obj_type = type(self._obj)

        result = _strftime(self._obj.data, date_format)

        return obj_type(
            result, name="strftime", coords=self._obj.coords, dims=self._obj.dims
        )

    def isocalendar(self) -> Dataset:
        """Dataset containing ISO year, week number, and weekday.

        Notes
        -----
        The iso year and weekday differ from the nominal year and weekday.
        """

        from xarray.core.dataset import Dataset

        if not is_np_datetime_like(self._obj.data.dtype):
            raise AttributeError("'CFTimeIndex' object has no attribute 'isocalendar'")

        values = _get_date_field(self._obj.data, "isocalendar", np.int64)

        obj_type = type(self._obj)
        data_vars = {}
        for i, name in enumerate(["year", "week", "weekday"]):
            data_vars[name] = obj_type(
                values[i], name=name, coords=self._obj.coords, dims=self._obj.dims
            )

        return Dataset(data_vars)

    @property
    def year(self) -> T_DataArray:
        """The year of the datetime"""
        return self._date_field("year", np.int64)

    @property
    def month(self) -> T_DataArray:
        """The month as January=1, December=12"""
        return self._date_field("month", np.int64)

    @property
    def day(self) -> T_DataArray:
        """The days of the datetime"""
        return self._date_field("day", np.int64)

    @property
    def hour(self) -> T_DataArray:
        """The hours of the datetime"""
        return self._date_field("hour", np.int64)

    @property
    def minute(self) -> T_DataArray:
        """The minutes of the datetime"""
        return self._date_field("minute", np.int64)

    @property
    def second(self) -> T_DataArray:
        """The seconds of the datetime"""
        return self._date_field("second", np.int64)

    @property
    def microsecond(self) -> T_DataArray:
        """The microseconds of the datetime"""
        return self._date_field("microsecond", np.int64)

    @property
    def nanosecond(self) -> T_DataArray:
        """The nanoseconds of the datetime"""
        return self._date_field("nanosecond", np.int64)

    @property
    def weekofyear(self) -> DataArray:
        "The week ordinal of the year"

        warnings.warn(
            "dt.weekofyear and dt.week have been deprecated. Please use "
            "dt.isocalendar().week instead.",
            FutureWarning,
            stacklevel=2,
        )

        weekofyear = self.isocalendar().week

        return weekofyear

    week = weekofyear

    @property
    def dayofweek(self) -> T_DataArray:
        """The day of the week with Monday=0, Sunday=6"""
        return self._date_field("dayofweek", np.int64)

    weekday = dayofweek

    @property
    def dayofyear(self) -> T_DataArray:
        """The ordinal day of the year"""
        return self._date_field("dayofyear", np.int64)

    @property
    def quarter(self) -> T_DataArray:
        """The quarter of the date"""
        return self._date_field("quarter", np.int64)

    @property
    def days_in_month(self) -> T_DataArray:
        """The number of days in the month"""
        return self._date_field("days_in_month", np.int64)

    daysinmonth = days_in_month

    @property
    def season(self) -> T_DataArray:
        """Season of the year"""
        return self._date_field("season", object)

    @property
    def time(self) -> T_DataArray:
        """Timestamps corresponding to datetimes"""
        return self._date_field("time", object)

    @property
    def date(self) -> T_DataArray:
        """Date corresponding to datetimes"""
        return self._date_field("date", object)

    @property
    def is_month_start(self) -> T_DataArray:
        """Indicate whether the date is the first day of the month"""
        return self._date_field("is_month_start", bool)

    @property
    def is_month_end(self) -> T_DataArray:
        """Indicate whether the date is the last day of the month"""
        return self._date_field("is_month_end", bool)

    @property
    def is_quarter_start(self) -> T_DataArray:
        """Indicate whether the date is the first day of a quarter"""
        return self._date_field("is_quarter_start", bool)

    @property
    def is_quarter_end(self) -> T_DataArray:
        """Indicate whether the date is the last day of a quarter"""
        return self._date_field("is_quarter_end", bool)

    @property
    def is_year_start(self) -> T_DataArray:
        """Indicate whether the date is the first day of a year"""
        return self._date_field("is_year_start", bool)

    @property
    def is_year_end(self) -> T_DataArray:
        """Indicate whether the date is the last day of the year"""
        return self._date_field("is_year_end", bool)

    @property
    def is_leap_year(self) -> T_DataArray:
        """Indicate if the date belongs to a leap year"""
        return self._date_field("is_leap_year", bool)

    @property
    def calendar(self) -> CFCalendar:
        """The name of the calendar of the dates.

        Only relevant for arrays of :py:class:`cftime.datetime` objects,
        returns "proleptic_gregorian" for arrays of :py:class:`numpy.datetime64` values.
        """
        return infer_calendar_name(self._obj.data)

    @property
    def days_in_year(self) -> T_DataArray:
        """Each datetime as the year plus the fraction of the year elapsed."""
        if self.calendar == "360_day":
            result = full_like(self.year, 360)
        else:
            result = self.is_leap_year.astype(int) + 365
        newvar = Variable(
            dims=self._obj.dims,
            attrs=self._obj.attrs,
            encoding=self._obj.encoding,
            data=result,
        )
        return self._obj._replace(newvar, name="days_in_year")

    @property
    def decimal_year(self) -> T_DataArray:
        """Convert the dates as a fractional year."""
        result = _decimal_year(self._obj)
        newvar = Variable(
            dims=self._obj.dims,
            attrs=self._obj.attrs,
            encoding=self._obj.encoding,
            data=result,
        )
        return self._obj._replace(newvar, name="decimal_year")


class TimedeltaAccessor(TimeAccessor[T_DataArray]):
    """Access Timedelta fields for DataArrays with Timedelta-like dtypes.

    Fields can be accessed through the `.dt` attribute for applicable DataArrays.

    Examples
    --------
    >>> dates = pd.timedelta_range(start="1 day", freq="6h", periods=20)
    >>> ts = xr.DataArray(dates, dims=("time"))
    >>> ts
    <xarray.DataArray (time: 20)> Size: 160B
    array([ 86400000000000, 108000000000000, 129600000000000, 151200000000000,
           172800000000000, 194400000000000, 216000000000000, 237600000000000,
           259200000000000, 280800000000000, 302400000000000, 324000000000000,
           345600000000000, 367200000000000, 388800000000000, 410400000000000,
           432000000000000, 453600000000000, 475200000000000, 496800000000000],
          dtype='timedelta64[ns]')
    Coordinates:
      * time     (time) timedelta64[ns] 160B 1 days 00:00:00 ... 5 days 18:00:00
    >>> ts.dt  # doctest: +ELLIPSIS
    <xarray.core.accessor_dt.TimedeltaAccessor object at 0x...>
    >>> ts.dt.days
    <xarray.DataArray 'days' (time: 20)> Size: 160B
    array([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5])
    Coordinates:
      * time     (time) timedelta64[ns] 160B 1 days 00:00:00 ... 5 days 18:00:00
    >>> ts.dt.microseconds
    <xarray.DataArray 'microseconds' (time: 20)> Size: 160B
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    Coordinates:
      * time     (time) timedelta64[ns] 160B 1 days 00:00:00 ... 5 days 18:00:00
    >>> ts.dt.seconds
    <xarray.DataArray 'seconds' (time: 20)> Size: 160B
    array([    0, 21600, 43200, 64800,     0, 21600, 43200, 64800,     0,
           21600, 43200, 64800,     0, 21600, 43200, 64800,     0, 21600,
           43200, 64800])
    Coordinates:
      * time     (time) timedelta64[ns] 160B 1 days 00:00:00 ... 5 days 18:00:00
    >>> ts.dt.total_seconds()
    <xarray.DataArray 'total_seconds' (time: 20)> Size: 160B
    array([ 86400., 108000., 129600., 151200., 172800., 194400., 216000.,
           237600., 259200., 280800., 302400., 324000., 345600., 367200.,
           388800., 410400., 432000., 453600., 475200., 496800.])
    Coordinates:
      * time     (time) timedelta64[ns] 160B 1 days 00:00:00 ... 5 days 18:00:00
    """

    @property
    def days(self) -> T_DataArray:
        """Number of days for each element"""
        return self._date_field("days", np.int64)

    @property
    def seconds(self) -> T_DataArray:
        """Number of seconds (>= 0 and less than 1 day) for each element"""
        return self._date_field("seconds", np.int64)

    @property
    def microseconds(self) -> T_DataArray:
        """Number of microseconds (>= 0 and less than 1 second) for each element"""
        return self._date_field("microseconds", np.int64)

    @property
    def nanoseconds(self) -> T_DataArray:
        """Number of nanoseconds (>= 0 and less than 1 microsecond) for each element"""
        return self._date_field("nanoseconds", np.int64)

    # Not defined as a property in order to match the Pandas API
    def total_seconds(self) -> T_DataArray:
        """Total duration of each element expressed in seconds."""
        return self._date_field("total_seconds", np.float64)


class CombinedDatetimelikeAccessor(
    DatetimeAccessor[T_DataArray], TimedeltaAccessor[T_DataArray]
):
    def __new__(cls, obj: T_DataArray) -> Self:
        # CombinedDatetimelikeAccessor isn't really instantiated. Instead
        # we need to choose which parent (datetime or timedelta) is
        # appropriate. Since we're checking the dtypes anyway, we'll just
        # do all the validation here.
        if not _contains_datetime_like_objects(obj.variable):
            # We use an AttributeError here so that `obj.dt` raises an error that
            # `getattr` expects; https://github.com/pydata/xarray/issues/8718. It's a
            # bit unusual in a `__new__`, but that's the only case where we use this
            # class.
            raise AttributeError(
                "'.dt' accessor only available for "
                "DataArray with datetime64 timedelta64 dtype or "
                "for arrays containing cftime datetime "
                "objects."
            )

        if is_np_timedelta_like(obj.dtype):
            return TimedeltaAccessor(obj)  # type: ignore[return-value]
        else:
            return DatetimeAccessor(obj)  # type: ignore[return-value]
