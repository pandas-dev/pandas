from __future__ import annotations

from datetime import (
    datetime,
    timedelta,
)
from typing import (
    TYPE_CHECKING,
    Self,
)
import warnings

import numpy as np

from pandas._libs import index as libindex
from pandas._libs.tslibs import (
    BaseOffset,
    Day,
    NaT,
    Period,
    Resolution,
    Tick,
    to_offset,
)
from pandas._libs.tslibs.dtypes import OFFSET_TO_PERIOD_FREQSTR
from pandas.errors import Pandas4Warning
from pandas.util._decorators import (
    set_module,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_integer,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import is_valid_na_for_dtype

from pandas.core.arrays.period import (
    PeriodArray,
    period_array,
    raise_on_incompatible,
    validate_dtype_freq,
)
import pandas.core.common as com
from pandas.core.indexes.base import maybe_extract_name
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin
from pandas.core.indexes.datetimes import (
    DatetimeIndex,
    Index,
)

if TYPE_CHECKING:
    from collections.abc import Hashable

    from pandas._typing import (
        Dtype,
        DtypeObj,
        npt,
    )


def _new_PeriodIndex(cls, **d):
    # GH13277 for unpickling
    values = d.pop("data")
    if values.dtype == "int64":
        freq = d.pop("freq", None)
        dtype = PeriodDtype(freq)
        values = PeriodArray(values, dtype=dtype)
        return cls._simple_new(values, **d)
    else:
        return cls(values, **d)


@set_module("pandas")
class PeriodIndex(DatetimeIndexOpsMixin):
    """
    Immutable ndarray holding ordinal values indicating regular periods in time.

    Index keys are boxed to Period objects which carries the metadata (eg,
    frequency information).

    Parameters
    ----------
    data : array-like (1d int np.ndarray or PeriodArray), optional
        Optional period-like data to construct index with.
    freq : str or period object, optional
        One of pandas period strings or corresponding objects.
    dtype : str or PeriodDtype, default None
        A dtype from which to extract a freq.
    copy : bool, default None
        Whether to copy input data, only relevant for array, Series, and Index
        inputs (for other input, e.g. a list, a new array is created anyway).
        Defaults to True for array input and False for Index/Series.
        Set to False to avoid copying array input at your own risk (if you
        know the input data won't be modified elsewhere).
        Set to True to force copying Series/Index input up front.
    name : str, default None
        Name of the resulting PeriodIndex.

    Attributes
    ----------
    day
    dayofweek
    day_of_week
    dayofyear
    day_of_year
    days_in_month
    daysinmonth
    end_time
    freq
    freqstr
    hour
    is_leap_year
    minute
    month
    quarter
    qyear
    second
    start_time
    week
    weekday
    weekofyear
    year

    Methods
    -------
    asfreq
    strftime
    to_timestamp
    from_fields
    from_ordinals

    Raises
    ------
    ValueError
        Passing the parameter data as a list without specifying either freq or
        dtype will raise a ValueError: "freq not specified and cannot be inferred"

    See Also
    --------
    Index : The base pandas Index type.
    Period : Represents a period of time.
    DatetimeIndex : Index with datetime64 data.
    TimedeltaIndex : Index of timedelta64 data.
    period_range : Create a fixed-frequency PeriodIndex.

    Examples
    --------
    >>> idx = pd.PeriodIndex(data=["2000Q1", "2002Q3"], freq="Q")
    >>> idx
    PeriodIndex(['2000Q1', '2002Q3'], dtype='period[Q-DEC]')
    """

    _typ = "periodindex"

    _data: PeriodArray
    freq: BaseOffset
    dtype: PeriodDtype

    _data_cls = PeriodArray
    _supports_partial_string_indexing = True

    _warn_quarter: bool = False

    @property
    def _engine_type(self) -> type[libindex.PeriodEngine]:
        return libindex.PeriodEngine

    # --------------------------------------------------------------------
    # methods that dispatch to array and wrap result in Index

    def asfreq(self, freq=None, how: str = "E") -> Self:
        """
        Convert the PeriodIndex to the specified frequency `freq`.

        Equivalent to applying :meth:`pandas.Period.asfreq` with the given arguments
        to each :class:`~pandas.Period` in this PeriodIndex.

        Parameters
        ----------
        freq : str
            A frequency.
        how : {'end', 'start', 'e', 's'}, default 'end'
            Whether the elements should be aligned to the end or start of
            each period, e.g. January 31st vs. January 1st. Case-insensitive.

        Returns
        -------
        PeriodIndex
            The transformed PeriodIndex with the new frequency.

        See Also
        --------
        arrays.PeriodArray.asfreq: Convert each Period in a PeriodArray to
            the given frequency.
        Period.asfreq : Convert a :class:`~pandas.Period` object to the given frequency.

        Examples
        --------
        >>> pidx = pd.period_range("2010-01-01", "2015-01-01", freq="Y")
        >>> pidx
        PeriodIndex(['2010', '2011', '2012', '2013', '2014', '2015'],
        dtype='period[Y-DEC]')

        >>> pidx.asfreq("M")
        PeriodIndex(['2010-12', '2011-12', '2012-12', '2013-12', '2014-12',
        '2015-12'], dtype='period[M]')

        >>> pidx.asfreq("M", how="S")
        PeriodIndex(['2010-01', '2011-01', '2012-01', '2013-01', '2014-01',
        '2015-01'], dtype='period[M]')
        """
        arr = self._data.asfreq(freq, how)
        return type(self)._simple_new(arr, name=self.name)

    def to_timestamp(self, freq=None, how: str = "start") -> DatetimeIndex:
        """
        Cast to DatetimeIndex.

        If possible, gives microsecond-unit DatetimeIndex. Otherwise
        gives nanosecond unit.

        Parameters
        ----------
        freq : str or DateOffset, optional
            Target frequency. The default is 'D' for week or longer,
            's' otherwise.
        how : {'start', 'end', 's', 'e'}, default 'start'
            Whether to use the start or end of the time period being converted.
            Case-insensitive.

        Returns
        -------
        DatetimeIndex
            Timestamp representation of given Period-like object.

        See Also
        --------
        PeriodIndex.day : The days of the period.
        PeriodIndex.from_fields : Construct a PeriodIndex from fields
            (year, month, day, etc.).
        PeriodIndex.from_ordinals : Construct a PeriodIndex from ordinals.
        PeriodIndex.hour : The hour of the period.
        PeriodIndex.minute : The minute of the period.
        PeriodIndex.month : The month as January=1, December=12.
        PeriodIndex.second : The second of the period.
        PeriodIndex.year : The year of the period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.to_timestamp()
        DatetimeIndex(['2023-01-01', '2023-02-01', '2023-03-01'],
        dtype='datetime64[us]', freq='MS')

        The frequency will not be inferred if the index contains less than
        three elements, or if the values of index are not strictly monotonic:

        >>> idx = pd.PeriodIndex(["2023-01", "2023-02"], freq="M")
        >>> idx.to_timestamp()
        DatetimeIndex(['2023-01-01', '2023-02-01'], dtype='datetime64[us]', freq=None)

        >>> idx = pd.PeriodIndex(
        ...     ["2023-01", "2023-02", "2023-02", "2023-03"], freq="2M"
        ... )
        >>> idx.to_timestamp()
        DatetimeIndex(['2023-01-01', '2023-02-01', '2023-02-01', '2023-03-01'],
        dtype='datetime64[us]', freq=None)
        """
        parr = self._data
        arr = parr.to_timestamp(freq, how)
        result = DatetimeIndex._simple_new(arr, name=self.name)
        result._freq = parr._to_timestamp_freq(arr, target_freq=freq, how=how)
        return result

    def strftime(self, date_format: str) -> Index:
        """
        Convert to Index using specified date_format.

        Return an Index of formatted strings specified by date_format, which
        supports the same string format as the python standard library. Details
        of the string format can be found in `python string format
        doc <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior>`__.
        :class:`Period` objects additionally support several directives not
        covered there, detailed in :meth:`Period.strftime`.

        Parameters
        ----------
        date_format : str
            Date format string (e.g. "%Y-%m-%d").

        Returns
        -------
        Index
            Index of formatted strings.

        See Also
        --------
        to_datetime : Convert the given argument to datetime.
        Period.strftime : Format a single Period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.strftime("%B %d, %Y")
        Index(['January 31, 2023', 'February 28, 2023', 'March 31, 2023'], dtype='str')
        """
        arr = self._data.strftime(date_format)
        return Index(arr, name=self.name, dtype=arr.dtype, copy=False)

    def _wrap_field(self, name: str) -> Index:
        result = getattr(self._data, name)
        return Index(result, name=self.name, dtype=result.dtype, copy=False)

    @property
    def start_time(self) -> DatetimeIndex:
        """
        Get the Timestamp for the start of the period.

        Returns a DatetimeIndex with the exact start Timestamp of each
        period in the index.

        Returns
        -------
        DatetimeIndex

        See Also
        --------
        PeriodIndex.end_time : Return the end Timestamp.
        PeriodIndex.to_timestamp : Cast to DatetimeIndex.
        Period.start_time : Return the start Timestamp for a single Period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02"], freq="M")
        >>> idx.start_time
        DatetimeIndex(['2023-01-01', '2023-02-01'], dtype='datetime64[us]', freq=None)
        """
        return DatetimeIndex(self._data.start_time, name=self.name, copy=False)

    @property
    def end_time(self) -> DatetimeIndex:
        """
        Get the Timestamp for the end of the period.

        Returns a DatetimeIndex with the last possible moment within each
        period in the index (e.g. 23:59:59.999999 for a daily period).

        Returns
        -------
        DatetimeIndex

        See Also
        --------
        PeriodIndex.start_time : Return the start Timestamp.
        PeriodIndex.to_timestamp : Cast to DatetimeIndex.
        Period.end_time : Return the end Timestamp for a single Period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02"], freq="M")
        >>> idx.end_time
        DatetimeIndex(['2023-01-31 23:59:59.999999', '2023-02-28 23:59:59.999999'],
                      dtype='datetime64[us]', freq=None)
        """
        return DatetimeIndex(self._data.end_time, name=self.name, copy=False)

    @property
    def year(self) -> Index:
        """
        The year of the period.

        Returns the year component for each period in the index.

        See Also
        --------
        PeriodIndex.day_of_year : The ordinal day of the year.
        PeriodIndex.is_leap_year : Logical indicating if the date belongs to a
            leap year.
        PeriodIndex.weekofyear : The week ordinal of the year.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023", "2024", "2025"], freq="Y")
        >>> idx.year
        Index([2023, 2024, 2025], dtype='int64')
        """
        return self._wrap_field("year")

    @property
    def month(self) -> Index:
        """
        The month as January=1, December=12.

        Returns the month component for each period in the index as an
        integer, where January is 1 and December is 12.

        See Also
        --------
        PeriodIndex.days_in_month : The number of days in the month.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.month
        Index([1, 2, 3], dtype='int64')
        """
        return self._wrap_field("month")

    @property
    def day(self) -> Index:
        """
        The days of the period.

        Returns the day-of-month component for each period in the index.

        See Also
        --------
        PeriodIndex.day_of_week : The day of the week with Monday=0, Sunday=6.
        PeriodIndex.day_of_year : The ordinal day of the year.
        PeriodIndex.days_in_month : The number of days in the month.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2020-01-31", "2020-02-28"], freq="D")
        >>> idx.day
        Index([31, 28], dtype='int64')
        """
        return self._wrap_field("day")

    @property
    def hour(self) -> Index:
        """
        The hour of the period.

        Returns the hour component for each period in the index.

        See Also
        --------
        PeriodIndex.minute : The minute of the period.
        PeriodIndex.second : The second of the period.
        PeriodIndex.to_timestamp : Cast to DatetimeArray/Index.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01-01 10:00", "2023-01-01 11:00"], freq="h")
        >>> idx.hour
        Index([10, 11], dtype='int64')
        """
        return self._wrap_field("hour")

    @property
    def minute(self) -> Index:
        """
        The minute of the period.

        Returns the minute component for each period in the index.

        See Also
        --------
        PeriodIndex.hour : The hour of the period.
        PeriodIndex.second : The second of the period.
        PeriodIndex.to_timestamp : Cast to DatetimeArray/Index.

        Examples
        --------
        >>> idx = pd.PeriodIndex(
        ...     ["2023-01-01 10:30:00", "2023-01-01 11:50:00"], freq="min"
        ... )
        >>> idx.minute
        Index([30, 50], dtype='int64')
        """
        return self._wrap_field("minute")

    @property
    def second(self) -> Index:
        """
        The second of the period.

        Returns the second component for each period in the index.

        See Also
        --------
        PeriodIndex.hour : The hour of the period.
        PeriodIndex.minute : The minute of the period.
        PeriodIndex.to_timestamp : Cast to DatetimeArray/Index.

        Examples
        --------
        >>> idx = pd.PeriodIndex(
        ...     ["2023-01-01 10:00:30", "2023-01-01 10:00:31"], freq="s"
        ... )
        >>> idx.second
        Index([30, 31], dtype='int64')
        """
        return self._wrap_field("second")

    @property
    def weekofyear(self) -> Index:
        """
        The week ordinal of the year.

        Returns the week number (1 through 53) for each period in the index.

        See Also
        --------
        PeriodIndex.day_of_week : The day of the week with Monday=0, Sunday=6.
        PeriodIndex.week : The week ordinal of the year.
        PeriodIndex.year : The year of the period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.week  # It can be written `weekofyear`
        Index([5, 9, 13], dtype='int64')
        """
        return self._wrap_field("weekofyear")

    week = weekofyear

    @property
    def day_of_week(self) -> Index:
        """
        The day of the week with Monday=0, Sunday=6.

        Returns the day-of-week component for each period, following the
        Python convention where Monday is 0 and Sunday is 6.

        See Also
        --------
        PeriodIndex.day : The days of the period.
        PeriodIndex.day_of_year : The ordinal day of the year.
        PeriodIndex.week : The week ordinal of the year.
        PeriodIndex.weekofyear : The week ordinal of the year.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01-01", "2023-01-02", "2023-01-03"], freq="D")
        >>> idx.day_of_week
        Index([6, 0, 1], dtype='int64')
        """
        return self._wrap_field("day_of_week")

    @property
    def weekday(self) -> Index:
        """
        The day of the week with Monday=0, Sunday=6.

        .. deprecated:: 3.1.0
            Use :attr:`PeriodIndex.day_of_week` instead.
        """
        return self._wrap_field("weekday")

    @property
    def dayofweek(self) -> Index:
        """
        The day of the week with Monday=0, Sunday=6.

        .. deprecated:: 3.1.0
            Use :attr:`PeriodIndex.day_of_week` instead.
        """
        return self._wrap_field("dayofweek")

    @property
    def day_of_year(self) -> Index:
        """
        The ordinal day of the year.

        Returns the day-of-year component for each period, ranging from
        1 (January 1st) to 365 or 366 for leap years.

        See Also
        --------
        PeriodIndex.day : The days of the period.
        PeriodIndex.day_of_week : The day of the week with Monday=0, Sunday=6.
        PeriodIndex.weekofyear : The week ordinal of the year.
        PeriodIndex.year : The year of the period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01-10", "2023-02-01", "2023-03-01"], freq="D")
        >>> idx.day_of_year
        Index([10, 32, 60], dtype='int64')

        >>> idx = pd.PeriodIndex(["2023", "2024", "2025"], freq="Y")
        >>> idx
        PeriodIndex(['2023', '2024', '2025'], dtype='period[Y-DEC]')
        >>> idx.day_of_year
        Index([365, 366, 365], dtype='int64')
        """
        return self._wrap_field("day_of_year")

    @property
    def dayofyear(self) -> Index:
        """
        The ordinal day of the year.

        .. deprecated:: 3.1.0
            Use :attr:`PeriodIndex.day_of_year` instead.
        """
        return self._wrap_field("dayofyear")

    @property
    def quarter(self) -> Index:
        """
        The quarter of the date.

        Returns the quarter (1 through 4) for each period in the index.

        See Also
        --------
        PeriodIndex.qyear : Fiscal year the Period lies in according to its
            starting-quarter.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.quarter
        Index([1, 1, 1], dtype='int64')
        """
        return self._wrap_field("quarter")

    @property
    def qyear(self) -> Index:
        """
        Fiscal year the Period lies in according to its starting-quarter.

        The `year` and the `qyear` of the period will be the same if the fiscal
        and calendar years are the same. When they are not, the fiscal year
        can be different from the calendar year of the period.

        Returns
        -------
        Index
            The fiscal year of each period.

        See Also
        --------
        PeriodIndex.quarter : The quarter of the date.
        PeriodIndex.year : The year of the period.

        Examples
        --------
        If the natural and fiscal year are the same, `qyear` and `year` will
        be the same.

        >>> idx = pd.PeriodIndex(["2018Q1"], freq="Q")
        >>> idx.qyear
        Index([2018], dtype='int64')
        >>> idx.year
        Index([2018], dtype='int64')

        If the fiscal year starts in April (`Q-MAR`), the first quarter of
        2018 will start in April 2017. `year` will then be 2017, but `qyear`
        will be the fiscal year, 2018.

        >>> idx = pd.PeriodIndex(["2018Q1"], freq="Q-MAR")
        >>> idx.start_time
        DatetimeIndex(['2017-04-01'], dtype='datetime64[us]', freq=None)
        >>> idx.qyear
        Index([2018], dtype='int64')
        >>> idx.year
        Index([2017], dtype='int64')
        """
        return self._wrap_field("qyear")

    @property
    def days_in_month(self) -> Index:
        """
        The number of days in the month.

        Returns the total number of days in the month of each period,
        accounting for leap years.

        See Also
        --------
        PeriodIndex.day : The days of the period.
        PeriodIndex.month : The month as January=1, December=12.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.days_in_month
        Index([31, 28, 31], dtype='int64')
        """
        return self._wrap_field("days_in_month")

    @property
    def daysinmonth(self) -> Index:
        """
        The number of days in the month.

        .. deprecated:: 3.1.0
            Use :attr:`PeriodIndex.days_in_month` instead.
        """
        return self._wrap_field("daysinmonth")

    @property
    def is_leap_year(self) -> npt.NDArray[np.bool_]:
        """
        Logical indicating if the date belongs to a leap year.

        Returns a boolean array where ``True`` indicates the period's year
        is a leap year.

        See Also
        --------
        PeriodIndex.qyear : Fiscal year the Period lies in according to its
            starting-quarter.
        PeriodIndex.year : The year of the period.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023", "2024", "2025"], freq="Y")
        >>> idx.is_leap_year
        array([False,  True, False])
        """
        return self._data.is_leap_year

    # ------------------------------------------------------------------------
    # Index Constructors

    def __new__(
        cls,
        data=None,
        freq=None,
        dtype: Dtype | None = None,
        copy: bool | None = None,
        name: Hashable | None = None,
    ) -> Self:
        refs = None
        if not copy and isinstance(data, (Index, ABCSeries)):
            refs = data._references
        if dtype is not None:
            dtype = pandas_dtype(dtype)

        name = maybe_extract_name(name, data, cls)

        if freq is not None:
            freq = to_offset(freq, is_period=True)
        dtype2 = PeriodDtype(freq) if freq is not None else None
        dtype = validate_dtype_freq(dtype, dtype2)
        if dtype is not None:
            freq = dtype._freq

        # GH#63388
        data, copy = cls._maybe_copy_array_input(data, copy, dtype)

        data = period_array(data=data, dtype=dtype)

        if copy:
            data = data.copy()

        return cls._simple_new(data, name=name, refs=refs)

    @classmethod
    def from_fields(
        cls,
        *,
        year=None,
        quarter=None,
        month=None,
        day=None,
        hour=None,
        minute=None,
        second=None,
        freq=None,
    ) -> Self:
        """
        Construct a PeriodIndex from fields (year, month, day, etc.).

        Each field (year, quarter, month, day, hour, minute, second) can be
        specified as a scalar or list-like. At least one field must be
        list-like; scalar fields are broadcast to its length. The frequency
        is inferred from the fields provided or can be given explicitly.

        Parameters
        ----------
        year : int or list-like, default None
            Year for the PeriodIndex.
        quarter : int or list-like, default None
            Quarter for the PeriodIndex.
        month : int or list-like, default None
            Month for the PeriodIndex.
        day : int or list-like, default None
            Day for the PeriodIndex.
        hour : int or list-like, default None
            Hour for the PeriodIndex.
        minute : int or list-like, default None
            Minute for the PeriodIndex.
        second : int or list-like, default None
            Second for the PeriodIndex.
        freq : str or period object, optional
            One of pandas period strings or corresponding objects.

        Returns
        -------
        PeriodIndex

        See Also
        --------
        PeriodIndex.from_ordinals : Construct a PeriodIndex from ordinals.
        PeriodIndex.to_timestamp : Cast to DatetimeArray/Index.

        Notes
        -----
        A list-like field must be a ``list``, ``tuple``, ``np.ndarray``, or
        ``Series``. Other array-like inputs, such as a ``range`` or an
        ``Index``, are not accepted.

        Examples
        --------
        >>> idx = pd.PeriodIndex.from_fields(year=[2000, 2002], quarter=[1, 3])
        >>> idx
        PeriodIndex(['2000Q1', '2002Q3'], dtype='period[Q-DEC]')
        """
        fields = {
            "year": year,
            "quarter": quarter,
            "month": month,
            "day": day,
            "hour": hour,
            "minute": minute,
            "second": second,
        }
        fields = {key: value for key, value in fields.items() if value is not None}
        arr = PeriodArray._from_fields(fields=fields, freq=freq)
        return cls._simple_new(arr)

    @classmethod
    def from_ordinals(cls, ordinals, *, freq, name=None) -> Self:
        """
        Construct a PeriodIndex from ordinals.

        Ordinals are integer offsets from the proleptic Gregorian epoch,
        interpreted according to the given frequency.

        Parameters
        ----------
        ordinals : array-like of int
            The period offsets from the proleptic Gregorian epoch.
        freq : str or period object
            One of pandas period strings or corresponding objects.
        name : str, default None
            Name of the resulting PeriodIndex.

        Returns
        -------
        PeriodIndex

        See Also
        --------
        PeriodIndex.from_fields : Construct a PeriodIndex from fields
            (year, month, day, etc.).
        PeriodIndex.to_timestamp : Cast to DatetimeArray/Index.

        Examples
        --------
        >>> idx = pd.PeriodIndex.from_ordinals([-1, 0, 1], freq="Q")
        >>> idx
        PeriodIndex(['1969Q4', '1970Q1', '1970Q2'], dtype='period[Q-DEC]')
        """
        ordinals = np.asarray(ordinals, dtype=np.int64)
        dtype = PeriodDtype(freq)
        data = PeriodArray._simple_new(ordinals, dtype=dtype)
        return cls._simple_new(data, name=name)

    # ------------------------------------------------------------------------
    # Data

    @property
    def values(self) -> npt.NDArray[np.object_]:
        warnings.warn(
            "PeriodIndex.values returning an object-dtype ndarray is "
            "deprecated. In a future version, this will return the "
            "underlying PeriodArray instead. Use 'PeriodIndex.to_numpy()' "
            "to get a NumPy array, or 'PeriodIndex.array' to get the "
            "ExtensionArray.",
            Pandas4Warning,
            stacklevel=2,
        )
        return np.asarray(self, dtype=object)

    def _mpl_repr(self) -> np.ndarray:
        # Return ordinals directly so matplotlib receives numeric x-values,
        # bypassing a round-trip through Period scalar objects.  GH#10578
        return self.asi8

    def _maybe_convert_timedelta(self, other) -> int | npt.NDArray[np.int64]:
        """
        Convert timedelta-like input to an integer multiple of self.freq

        Parameters
        ----------
        other : timedelta, np.timedelta64, DateOffset, int, np.ndarray

        Returns
        -------
        converted : int, np.ndarray[int64]

        Raises
        ------
        IncompatibleFrequency : if the input cannot be written as a multiple
            of self.freq.  Note IncompatibleFrequency subclasses ValueError.
        """
        if isinstance(other, (timedelta, np.timedelta64, Tick, np.ndarray)):
            if isinstance(self.freq, (Tick, Day)):
                # _check_timedeltalike_freq_compat will raise if incompatible
                delta = self._data._check_timedeltalike_freq_compat(other)
                return delta
        elif isinstance(other, BaseOffset):
            if other.base == self.freq.base:
                return other.n

            raise raise_on_incompatible(self, other)
        elif is_integer(other):
            assert isinstance(other, int)
            return other

        # raise when input doesn't have freq
        raise raise_on_incompatible(self, None)

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        """
        Can we compare values of the given dtype to our own?
        """
        return self.dtype == dtype

    # ------------------------------------------------------------------------
    # Index Methods

    def asof_locs(self, where: Index, mask: npt.NDArray[np.bool_]) -> np.ndarray:
        """
        where : array of timestamps
        mask : np.ndarray[bool]
            Array of booleans where data is not NA.
        """
        if isinstance(where, DatetimeIndex):
            where = PeriodIndex(where._values, freq=self.freq, copy=False)
        elif not isinstance(where, PeriodIndex):
            raise TypeError("asof_locs `where` must be DatetimeIndex or PeriodIndex")

        return super().asof_locs(where, mask)

    @property
    def is_full(self) -> bool:
        """
        Return True if the index contains all periods from start to end
        (inclusive) with no gaps.

        Requires monotonic increasing order. Duplicate periods are allowed.

        .. deprecated:: 3.1.0
            ``PeriodIndex.is_full`` is deprecated and will be removed in
            a future version. Use
            ``index.empty or len(index.unique()) ==
            len(period_range(index.min(), index.max(), freq=index.freq))``
            instead. Unlike ``is_full``, this does not raise on a
            non-monotonic index.
        """
        warnings.warn(
            "PeriodIndex.is_full is deprecated and will be removed in a "
            "future version. Use index.empty or len(index.unique()) == "
            "len(period_range(index.min(), index.max(), freq=index.freq)) "
            "instead.",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )
        if len(self) == 0:
            return True
        if not self.is_monotonic_increasing:
            raise ValueError("Index is not monotonic")
        values = self.asi8
        return bool(((values[1:] - values[:-1]) < 2).all())

    @property
    def inferred_type(self) -> str:
        # b/c data is represented as ints make sure we can't have ambiguous
        # indexing
        return "period"

    # ------------------------------------------------------------------------
    # Indexing Methods

    def _convert_tolerance(self, tolerance, target):
        # Returned tolerance must be in dtype/units so that
        #  `|self._get_engine_target() - target._engine_target()| <= tolerance`
        #  is meaningful.  Since PeriodIndex returns int64 for engine_target,
        #  we may need to convert timedelta64 tolerance to int64.
        tolerance = super()._convert_tolerance(tolerance, target)

        if self.dtype == target.dtype:
            # convert tolerance to i8
            tolerance = self._maybe_convert_timedelta(tolerance)

        return tolerance

    def get_loc(self, key):
        """
        Get integer location for requested label.

        Parameters
        ----------
        key : Period, NaT, str, or datetime
            String or datetime key must be parsable as Period.

        Returns
        -------
        loc : int or ndarray[int64]

        Raises
        ------
        KeyError
            Key is not present in the index.
        TypeError
            If key is listlike or otherwise not hashable.
        """
        orig_key = key

        self._check_indexing_error(key)

        if is_valid_na_for_dtype(key, self.dtype):
            key = NaT

        elif isinstance(key, str):
            try:
                parsed, reso = self._parse_with_reso(key)
            except ValueError as err:
                # A string with invalid format
                raise KeyError(f"Cannot interpret '{key}' as period") from err

            if self._can_partial_date_slice(reso):
                try:
                    return self._partial_date_slice(reso, parsed)
                except KeyError as err:
                    raise KeyError(key) from err

            if reso == self._resolution_obj:
                # the reso < self._resolution_obj case goes
                #  through _get_string_slice
                key = self._cast_partial_indexing_scalar(parsed)
            else:
                raise KeyError(key)

        elif isinstance(key, Period):
            self._disallow_mismatched_indexing(key)

        elif isinstance(key, datetime):
            key = self._cast_partial_indexing_scalar(key)

        else:
            # in particular integer, which Period constructor would cast to string
            raise KeyError(key)

        try:
            return Index.get_loc(self, key)
        except KeyError as err:
            raise KeyError(orig_key) from err

    def _disallow_mismatched_indexing(self, key: Period) -> None:
        if key._dtype != self.dtype:
            raise KeyError(key)

    def _cast_partial_indexing_scalar(self, label: datetime) -> Period:
        try:
            period = Period(label, freq=self.freq)
        except ValueError as err:
            # we cannot construct the Period
            raise KeyError(label) from err
        return period

    def _maybe_cast_slice_bound(self, label, side: str):
        """
        If label is a string, cast it to scalar type according to resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}

        Returns
        -------
        label : object

        Notes
        -----
        Value of `side` parameter should be validated in caller.
        """
        if isinstance(label, datetime):
            label = self._cast_partial_indexing_scalar(label)

        return super()._maybe_cast_slice_bound(label, side)

    def _parsed_string_to_bounds(self, reso: Resolution, parsed: datetime):
        freq = OFFSET_TO_PERIOD_FREQSTR.get(reso.attr_abbrev, reso.attr_abbrev)
        iv = Period(parsed, freq=freq)
        return (iv.asfreq(self.freq, how="start"), iv.asfreq(self.freq, how="end"))

    def shift(self, periods: int = 1, freq=None) -> Self:
        """
        Shift index by desired number of time frequency increments.

        This method is for shifting the values of datetime-like indexes
        by a specified time increment a given number of times.

        Parameters
        ----------
        periods : int, default 1
            Number of periods (or increments) to shift by,
            can be positive or negative.
        freq : pandas.DateOffset, pandas.Timedelta or string, optional
            Frequency increment to shift by.
            If None, the index is shifted by its own `freq` attribute.
            Offset aliases are valid strings, e.g., 'D', 'W', 'M' etc.

        Returns
        -------
        pandas.DatetimeIndex
            Shifted index.

        See Also
        --------
        Index.shift : Shift values of Index.
        PeriodIndex.shift : Shift values of PeriodIndex.
        """
        if freq is not None:
            raise TypeError(
                f"`freq` argument is not supported for {type(self).__name__}.shift"
            )
        return self + periods


@set_module("pandas")
def period_range(
    start=None,
    end=None,
    periods: int | None = None,
    freq=None,
    name: Hashable | None = None,
) -> PeriodIndex:
    """
    Return a fixed frequency PeriodIndex.

    The day (calendar) is the default frequency.

    Parameters
    ----------
    start : str, datetime, date, pandas.Timestamp, or period-like, default None
        Left bound for generating periods.
    end : str, datetime, date, pandas.Timestamp, or period-like, default None
        Right bound for generating periods.
    periods : int, default None
        Number of periods to generate.
    freq : str or DateOffset, optional
        Frequency alias. By default the freq is taken from `start` or `end`
        if those are Period objects. Otherwise, the default is ``"D"`` for
        daily frequency.
    name : str, default None
        Name of the resulting PeriodIndex.

    Returns
    -------
    PeriodIndex
        A PeriodIndex of fixed frequency periods.

    See Also
    --------
    date_range : Returns a fixed frequency DatetimeIndex.
    Period : Represents a period of time.
    PeriodIndex : Immutable ndarray holding ordinal values indicating regular periods
        in time.

    Notes
    -----
    Of the three parameters: ``start``, ``end``, and ``periods``, exactly two
    must be specified.

    To learn more about the frequency strings, please see
    :ref:`this link<timeseries.offset_aliases>`.

    Examples
    --------
    >>> pd.period_range(start="2017-01-01", end="2018-01-01", freq="M")
    PeriodIndex(['2017-01', '2017-02', '2017-03', '2017-04', '2017-05', '2017-06',
             '2017-07', '2017-08', '2017-09', '2017-10', '2017-11', '2017-12',
             '2018-01'],
            dtype='period[M]')

    If ``start`` or ``end`` are ``Period`` objects, they will be used as anchor
    endpoints for a ``PeriodIndex`` with frequency matching that of the
    ``period_range`` constructor.

    >>> pd.period_range(
    ...     start=pd.Period("2017Q1", freq="Q"),
    ...     end=pd.Period("2017Q2", freq="Q"),
    ...     freq="M",
    ... )
    PeriodIndex(['2017-03', '2017-04', '2017-05', '2017-06'],
                dtype='period[M]')
    """
    if com.count_not_none(start, end, periods) != 2:
        raise ValueError(
            "Of the three parameters: start, end, and periods, "
            "exactly two must be specified"
        )
    if freq is None and (not isinstance(start, Period) and not isinstance(end, Period)):
        freq = "D"

    data, freq = PeriodArray._generate_range(start, end, periods, freq)
    dtype = PeriodDtype(freq)
    data = PeriodArray(data, dtype=dtype)
    return PeriodIndex(data, name=name, copy=False)
