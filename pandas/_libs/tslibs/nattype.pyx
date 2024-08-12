from cpython.datetime cimport (
    PyDate_Check,
    PyDateTime_Check,
    PyDelta_Check,
    datetime,
    import_datetime,
)

import_datetime()
from cpython.object cimport (
    Py_EQ,
    Py_NE,
    PyObject_RichCompare,
)

import numpy as np

cimport numpy as cnp
from numpy cimport int64_t

cnp.import_array()

cimport pandas._libs.tslibs.util as util

# ----------------------------------------------------------------------
# Constants
nat_strings = {"NaT", "nat", "NAT", "nan", "NaN", "NAN"}
cdef set c_nat_strings = nat_strings

cdef int64_t NPY_NAT = util.get_nat()
iNaT = NPY_NAT  # python-visible constant

# ----------------------------------------------------------------------


def _make_nan_func(func_name: str, doc: str):
    def f(*args, **kwargs):
        return np.nan
    f.__name__ = func_name
    f.__doc__ = doc
    return f


def _make_nat_func(func_name: str, doc: str):
    def f(*args, **kwargs):
        return c_NaT
    f.__name__ = func_name
    f.__doc__ = doc
    return f


def _make_error_func(func_name: str, cls):
    def f(*args, **kwargs):
        raise ValueError(f"NaTType does not support {func_name}")

    f.__name__ = func_name
    if isinstance(cls, str):
        # passed the literal docstring directly
        f.__doc__ = cls
    elif cls is not None:
        f.__doc__ = getattr(cls, func_name).__doc__
    return f


cdef _nat_divide_op(self, other):
    if PyDelta_Check(other) or cnp.is_timedelta64_object(other) or other is c_NaT:
        return np.nan
    if util.is_integer_object(other) or util.is_float_object(other):
        return c_NaT
    return NotImplemented


cdef _nat_rdivide_op(self, other):
    if PyDelta_Check(other):
        return np.nan
    return NotImplemented


def _nat_unpickle(*args):
    # return constant defined in the module
    return c_NaT

# ----------------------------------------------------------------------


cdef class _NaT(datetime):
    # cdef readonly:
    #    int64_t value

    # higher than np.ndarray and np.matrix
    __array_priority__ = 100

    def __richcmp__(_NaT self, object other, int op):
        if cnp.is_datetime64_object(other) or PyDateTime_Check(other):
            # We treat NaT as datetime-like for this comparison
            return op == Py_NE

        elif cnp.is_timedelta64_object(other) or PyDelta_Check(other):
            # We treat NaT as timedelta-like for this comparison
            return op == Py_NE

        elif util.is_array(other):
            if other.dtype.kind in "mM":
                result = np.empty(other.shape, dtype=np.bool_)
                result.fill(op == Py_NE)
            elif other.dtype.kind == "O":
                result = np.array([PyObject_RichCompare(self, x, op) for x in other])
            elif op == Py_EQ:
                result = np.zeros(other.shape, dtype=bool)
            elif op == Py_NE:
                result = np.ones(other.shape, dtype=bool)
            else:
                return NotImplemented
            return result

        elif PyDate_Check(other):
            # GH#39151 don't defer to datetime.date object
            if op == Py_EQ:
                return False
            if op == Py_NE:
                return True
            raise TypeError("Cannot compare NaT with datetime.date object")

        return NotImplemented

    def __add__(self, other):
        if PyDateTime_Check(other):
            return c_NaT
        elif PyDelta_Check(other):
            return c_NaT
        elif cnp.is_datetime64_object(other) or cnp.is_timedelta64_object(other):
            return c_NaT

        elif util.is_integer_object(other):
            # For Period compat
            return c_NaT

        elif util.is_array(other):
            if other.dtype.kind in "mM":
                # If we are adding to datetime64, we treat NaT as timedelta
                #  Either way, result dtype is datetime64
                result = np.empty(other.shape, dtype="datetime64[ns]")
                result.fill("NaT")
                return result
            raise TypeError(f"Cannot add NaT to ndarray with dtype {other.dtype}")

        # Includes Period, DateOffset going through here
        return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        # Duplicate some logic from _Timestamp.__sub__ to avoid needing
        # to subclass; allows us to @final(_Timestamp.__sub__)

        if PyDateTime_Check(other):
            return c_NaT
        elif PyDelta_Check(other):
            return c_NaT
        elif cnp.is_datetime64_object(other) or cnp.is_timedelta64_object(other):
            return c_NaT

        elif util.is_integer_object(other):
            # For Period compat
            return c_NaT

        elif util.is_array(other):
            if other.dtype.kind == "m":
                # NaT - timedelta64 we treat NaT as datetime64, so result
                #  is datetime64
                result = np.empty(other.shape, dtype="datetime64[ns]")
                result.fill("NaT")
                return result

            elif other.dtype.kind == "M":
                # We treat NaT as a datetime, so regardless of whether this is
                #  NaT - other or other - NaT, the result is timedelta64
                result = np.empty(other.shape, dtype="timedelta64[ns]")
                result.fill("NaT")
                return result

            raise TypeError(
                f"Cannot subtract NaT from ndarray with dtype {other.dtype}"
            )

        # Includes Period, DateOffset going through here
        return NotImplemented

    def __rsub__(self, other):
        if util.is_array(other):
            if other.dtype.kind == "m":
                # timedelta64 - NaT we have to treat NaT as timedelta64
                #  for this to be meaningful, and the result is timedelta64
                result = np.empty(other.shape, dtype="timedelta64[ns]")
                result.fill("NaT")
                return result

            elif other.dtype.kind == "M":
                # We treat NaT as a datetime, so regardless of whether this is
                #  NaT - other or other - NaT, the result is timedelta64
                result = np.empty(other.shape, dtype="timedelta64[ns]")
                result.fill("NaT")
                return result
        # other cases are same, swap operands is allowed even though we subtract
        # because this is NaT
        return self.__sub__(other)

    def __pos__(self):
        return NaT

    def __neg__(self):
        return NaT

    def __truediv__(self, other):
        return _nat_divide_op(self, other)

    def __floordiv__(self, other):
        return _nat_divide_op(self, other)

    def __mul__(self, other):
        if util.is_integer_object(other) or util.is_float_object(other):
            return NaT
        return NotImplemented

    @property
    def asm8(self) -> np.datetime64:
        return np.datetime64(NPY_NAT, "ns")

    def to_datetime64(self) -> np.datetime64:
        """
        Return a numpy.datetime64 object with same precision.

        Examples
        --------
        >>> ts = pd.Timestamp(year=2023, month=1, day=1,
        ...                   hour=10, second=15)
        >>> ts
        Timestamp('2023-01-01 10:00:15')
        >>> ts.to_datetime64()
        numpy.datetime64('2023-01-01T10:00:15.000000')
        """
        return np.datetime64("NaT", "ns")

    def to_numpy(self, dtype=None, copy=False) -> np.datetime64 | np.timedelta64:
        """
        Convert the Timestamp to a NumPy datetime64 or timedelta64.

        With the default 'dtype', this is an alias method for `NaT.to_datetime64()`.

        The copy parameter is available here only for compatibility. Its value
        will not affect the return value.

        Returns
        -------
        numpy.datetime64 or numpy.timedelta64

        See Also
        --------
        DatetimeIndex.to_numpy : Similar method for DatetimeIndex.

        Examples
        --------
        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')
        >>> ts.to_numpy()
        numpy.datetime64('2020-03-14T15:32:52.192548651')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.to_numpy()
        numpy.datetime64('NaT')

        >>> pd.NaT.to_numpy("m8[ns]")
        numpy.timedelta64('NaT','ns')
        """
        if dtype is not None:
            # GH#44460
            dtype = np.dtype(dtype)
            if dtype.kind == "M":
                return np.datetime64("NaT").astype(dtype)
            elif dtype.kind == "m":
                return np.timedelta64("NaT").astype(dtype)
            else:
                raise ValueError(
                    "NaT.to_numpy dtype must be a datetime64 dtype, timedelta64 "
                    "dtype, or None."
                )
        return self.to_datetime64()

    def __repr__(self) -> str:
        return "NaT"

    def __str__(self) -> str:
        return "NaT"

    def isoformat(self, sep: str = "T", timespec: str = "auto") -> str:
        # This allows Timestamp(ts.isoformat()) to always correctly roundtrip.
        return "NaT"

    def __hash__(self) -> int:
        return NPY_NAT

    @property
    def is_leap_year(self) -> bool:
        return False

    @property
    def is_month_start(self) -> bool:
        return False

    @property
    def is_quarter_start(self) -> bool:
        return False

    @property
    def is_year_start(self) -> bool:
        return False

    @property
    def is_month_end(self) -> bool:
        return False

    @property
    def is_quarter_end(self) -> bool:
        return False

    @property
    def is_year_end(self) -> bool:
        return False


class NaTType(_NaT):
    """
    (N)ot-(A)-(T)ime, the time equivalent of NaN.

    Examples
    --------
    >>> pd.DataFrame([pd.Timestamp("2023"), np.nan], columns=["col_1"])
            col_1
    0  2023-01-01
    1         NaT
    """

    def __new__(cls):
        cdef _NaT base

        base = _NaT.__new__(cls, 1, 1, 1)
        base._value= NPY_NAT

        return base

    @property
    def value(self) -> int:
        return self._value

    def __reduce_ex__(self, protocol):
        # python 3.6 compat
        # https://bugs.python.org/issue28730
        # now __reduce_ex__ is defined and higher priority than __reduce__
        return self.__reduce__()

    def __reduce__(self):
        return (_nat_unpickle, (None, ))

    def __rtruediv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rfloordiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rmul__(self, other):
        if util.is_integer_object(other) or util.is_float_object(other):
            return c_NaT
        return NotImplemented

    # ----------------------------------------------------------------------
    # inject the Timestamp field properties
    # these by definition return np.nan

    year = property(fget=lambda self: np.nan)
    quarter = property(fget=lambda self: np.nan)
    month = property(fget=lambda self: np.nan)
    day = property(fget=lambda self: np.nan)
    hour = property(fget=lambda self: np.nan)
    minute = property(fget=lambda self: np.nan)
    second = property(fget=lambda self: np.nan)
    millisecond = property(fget=lambda self: np.nan)
    microsecond = property(fget=lambda self: np.nan)
    nanosecond = property(fget=lambda self: np.nan)

    week = property(fget=lambda self: np.nan)
    dayofyear = property(fget=lambda self: np.nan)
    day_of_year = property(fget=lambda self: np.nan)
    weekofyear = property(fget=lambda self: np.nan)
    days_in_month = property(fget=lambda self: np.nan)
    daysinmonth = property(fget=lambda self: np.nan)
    dayofweek = property(fget=lambda self: np.nan)
    day_of_week = property(fget=lambda self: np.nan)

    # inject Timedelta properties
    days = property(fget=lambda self: np.nan)
    seconds = property(fget=lambda self: np.nan)
    microseconds = property(fget=lambda self: np.nan)
    nanoseconds = property(fget=lambda self: np.nan)

    # inject pd.Period properties
    qyear = property(fget=lambda self: np.nan)

    # ----------------------------------------------------------------------
    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    # These are the ones that can get their docstrings from datetime.

    # nan methods
    weekday = _make_nan_func(
        "weekday",
        """
        Return the day of the week represented by the date.

        Monday == 0 ... Sunday == 6.

        See Also
        --------
        Timestamp.dayofweek : Return the day of the week with Monday=0, Sunday=6.
        Timestamp.isoweekday : Return the day of the week with Monday=1, Sunday=7.
        datetime.date.weekday : Equivalent method in datetime module.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01')
        >>> ts
        Timestamp('2023-01-01  00:00:00')
        >>> ts.weekday()
        6
        """,
    )
    isoweekday = _make_nan_func(
        "isoweekday",
        """
        Return the day of the week represented by the date.

        Monday == 1 ... Sunday == 7.

        See Also
        --------
        Timestamp.weekday : Return the day of the week with Monday=0, Sunday=6.
        Timestamp.isocalendar : Return a tuple containing ISO year, week number
            and weekday.
        datetime.date.isoweekday : Equivalent method in datetime module.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00')
        >>> ts
        Timestamp('2023-01-01 10:00:00')
        >>> ts.isoweekday()
        7
        """,
    )
    total_seconds = _make_nan_func(
        "total_seconds",
        """
        Total seconds in the duration.

        Examples
        --------
        >>> td = pd.Timedelta('1min')
        >>> td
        Timedelta('0 days 00:01:00')
        >>> td.total_seconds()
        60.0
        """,
    )
    month_name = _make_nan_func(
        "month_name",
        """
        Return the month name of the Timestamp with specified locale.

        Parameters
        ----------
        locale : str, default None (English locale)
            Locale determining the language in which to return the month name.

        Returns
        -------
        str

        Examples
        --------
        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')
        >>> ts.month_name()
        'March'

        Analogous for ``pd.NaT``:

        >>> pd.NaT.month_name()
        nan
        """,
    )
    day_name = _make_nan_func(
        "day_name",
        """
        Return the day name of the Timestamp with specified locale.

        Parameters
        ----------
        locale : str, default None (English locale)
            Locale determining the language in which to return the day name.

        Returns
        -------
        str

        See Also
        --------
        Timestamp.day_of_week : Return day of the week.
        Timestamp.day_of_year : Return day of the year.

        Examples
        --------
        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')
        >>> ts.day_name()
        'Saturday'

        Analogous for ``pd.NaT``:

        >>> pd.NaT.day_name()
        nan
        """,
    )
    # _nat_methods

    # "fromisocalendar" was introduced in 3.8
    fromisocalendar = _make_error_func("fromisocalendar", datetime)

    # ----------------------------------------------------------------------
    # The remaining methods have docstrings copy/pasted from the analogous
    # Timestamp methods.
    isocalendar = _make_error_func(
        "isocalendar",
        """
        Return a named tuple containing ISO year, week number, and weekday.

        See Also
        --------
        DatetimeIndex.isocalendar : Return a 3-tuple containing ISO year,
            week number, and weekday for the given DatetimeIndex object.
        datetime.date.isocalendar : The equivalent method for `datetime.date` objects.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00')
        >>> ts
        Timestamp('2023-01-01 10:00:00')
        >>> ts.isocalendar()
        datetime.IsoCalendarDate(year=2022, week=52, weekday=7)
        """
        )
    dst = _make_error_func(
        "dst",
        """
        Return the daylight saving time (DST) adjustment.

        This method returns the DST adjustment as a `datetime.timedelta` object
        if the Timestamp is timezone-aware and DST is applicable.

        See Also
        --------
        Timestamp.tz_localize : Localize the Timestamp to a timezone.
        Timestamp.tz_convert : Convert timezone-aware Timestamp to another time zone.

        Examples
        --------
        >>> ts = pd.Timestamp('2000-06-01 00:00:00', tz='Europe/Brussels')
        >>> ts
        Timestamp('2000-06-01 00:00:00+0200', tz='Europe/Brussels')
        >>> ts.dst()
        datetime.timedelta(seconds=3600)
        """
        )
    date = _make_nat_func(
        "date",
        """
        Return date object with same year, month and day.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00.00')
        >>> ts
        Timestamp('2023-01-01 10:00:00')
        >>> ts.date()
        datetime.date(2023, 1, 1)
        """
        )
    utctimetuple = _make_error_func(
        "utctimetuple",
        """
        Return UTC time tuple, compatible with time.localtime().

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00', tz='Europe/Brussels')
        >>> ts
        Timestamp('2023-01-01 10:00:00+0100', tz='Europe/Brussels')
        >>> ts.utctimetuple()
        time.struct_time(tm_year=2023, tm_mon=1, tm_mday=1, tm_hour=9,
        tm_min=0, tm_sec=0, tm_wday=6, tm_yday=1, tm_isdst=0)
        """
        )
    utcoffset = _make_error_func(
        "utcoffset",
        """
        Return utc offset.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00', tz='Europe/Brussels')
        >>> ts
        Timestamp('2023-01-01 10:00:00+0100', tz='Europe/Brussels')
        >>> ts.utcoffset()
        datetime.timedelta(seconds=3600)
        """
        )
    tzname = _make_error_func(
        "tzname",
        """
        Return time zone name.

        This method returns the name of the Timestamp's time zone as a string.

        See Also
        --------
        Timestamp.tzinfo : Returns the timezone information of the Timestamp.
        Timestamp.tz_convert : Convert timezone-aware Timestamp to another time zone.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00', tz='Europe/Brussels')
        >>> ts
        Timestamp('2023-01-01 10:00:00+0100', tz='Europe/Brussels')
        >>> ts.tzname()
        'CET'
        """
        )
    time = _make_error_func(
        "time",
        """
        Return time object with same time but with tzinfo=None.

        This method extracts the time part of the `Timestamp` object, excluding any
        timezone information. It returns a `datetime.time` object which only represents
        the time (hours, minutes, seconds, and microseconds).

        See Also
        --------
        Timestamp.date : Return date object with same year, month and day.
        Timestamp.tz_convert : Convert timezone-aware Timestamp to another time zone.
        Timestamp.tz_localize : Localize the Timestamp to a timezone.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00')
        >>> ts
        Timestamp('2023-01-01 10:00:00')
        >>> ts.time()
        datetime.time(10, 0)
        """,
        )
    timetuple = _make_error_func(
        "timetuple",
        """
        Return time tuple, compatible with time.localtime().

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00')
        >>> ts
        Timestamp('2023-01-01 10:00:00')
        >>> ts.timetuple()
        time.struct_time(tm_year=2023, tm_mon=1, tm_mday=1,
        tm_hour=10, tm_min=0, tm_sec=0, tm_wday=6, tm_yday=1, tm_isdst=-1)
        """
        )
    timetz = _make_error_func(
        "timetz",
        """
        Return time object with same time and tzinfo.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00', tz='Europe/Brussels')
        >>> ts
        Timestamp('2023-01-01 10:00:00+0100', tz='Europe/Brussels')
        >>> ts.timetz()
        datetime.time(10, 0, tzinfo=<DstTzInfo 'Europe/Brussels' CET+1:00:00 STD>)
        """
        )
    toordinal = _make_error_func(
        "toordinal",
        """
        Return proleptic Gregorian ordinal. January 1 of year 1 is day 1.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:50')
        >>> ts
        Timestamp('2023-01-01 10:00:50')
        >>> ts.toordinal()
        738521
        """
        )
    ctime = _make_error_func(
        "ctime",
        """
        Return ctime() style string.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 10:00:00.00')
        >>> ts
        Timestamp('2023-01-01 10:00:00')
        >>> ts.ctime()
        'Sun Jan  1 10:00:00 2023'
        """,
    )

    strftime = _make_error_func(
        "strftime",
        """
        Return a formatted string of the Timestamp.

        Parameters
        ----------
        format : str
            Format string to convert Timestamp to string.
            See strftime documentation for more information on the format string:
            https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior.

        See Also
        --------
        Timestamp.isoformat : Return the time formatted according to ISO 8601.
        pd.to_datetime : Convert argument to datetime.
        Period.strftime : Format a single Period.

        Examples
        --------
        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')
        >>> ts.strftime('%Y-%m-%d %X')
        '2020-03-14 15:32:52'
        """,
    )

    strptime = _make_error_func(
        "strptime",
        """
        Timestamp.strptime(string, format)

        Function is not implemented. Use pd.to_datetime().

        Examples
        --------
        >>> pd.Timestamp.strptime("2023-01-01", "%d/%m/%y")
        Traceback (most recent call last):
        NotImplementedError
        """,
    )

    utcfromtimestamp = _make_error_func(
        "utcfromtimestamp",
        """
        Timestamp.utcfromtimestamp(ts)

        Construct a timezone-aware UTC datetime from a POSIX timestamp.

        This method creates a datetime object from a POSIX timestamp, keeping the
        Timestamp object's timezone.

        Parameters
        ----------
        ts : float
            POSIX timestamp.

        See Also
        --------
        Timezone.tzname : Return time zone name.
        Timestamp.utcnow : Return a new Timestamp representing UTC day and time.
        Timestamp.fromtimestamp : Transform timestamp[, tz] to tz's local
            time from POSIX timestamp.

        Notes
        -----
        Timestamp.utcfromtimestamp behavior differs from datetime.utcfromtimestamp
        in returning a timezone-aware object.

        Examples
        --------
        >>> pd.Timestamp.utcfromtimestamp(1584199972)
        Timestamp('2020-03-14 15:32:52+0000', tz='UTC')
        """,
    )
    fromtimestamp = _make_error_func(
        "fromtimestamp",
        """
        Timestamp.fromtimestamp(ts)

        Transform timestamp[, tz] to tz's local time from POSIX timestamp.

        Examples
        --------
        >>> pd.Timestamp.fromtimestamp(1584199972)  # doctest: +SKIP
        Timestamp('2020-03-14 15:32:52')

        Note that the output may change depending on your local time.
        """,
    )
    combine = _make_error_func(
        "combine",
        """
        Timestamp.combine(date, time)

        Combine date, time into datetime with same date and time fields.

        Examples
        --------
        >>> from datetime import date, time
        >>> pd.Timestamp.combine(date(2020, 3, 14), time(15, 30, 15))
        Timestamp('2020-03-14 15:30:15')
        """,
    )
    utcnow = _make_error_func(
        "utcnow",
        """
        Timestamp.utcnow()

        Return a new Timestamp representing UTC day and time.

        See Also
        --------
        Timestamp : Constructs an arbitrary datetime.
        Timestamp.now : Return the current local date and time, which
            can be timezone-aware.
        Timestamp.today : Return the current local date and time with
            timezone information set to None.
        to_datetime : Convert argument to datetime.
        date_range : Return a fixed frequency DatetimeIndex.
        Timestamp.utctimetuple : Return UTC time tuple, compatible with
            time.localtime().

        Examples
        --------
        >>> pd.Timestamp.utcnow()   # doctest: +SKIP
        Timestamp('2020-11-16 22:50:18.092888+0000', tz='UTC')
        """,
    )

    timestamp = _make_error_func(
        "timestamp",
        """
        Return POSIX timestamp as float.

        Examples
        --------
        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548')
        >>> ts.timestamp()
        1584199972.192548
        """
    )

    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    astimezone = _make_error_func(
        "astimezone",
        """
        Convert timezone-aware Timestamp to another time zone.

        This method is used to convert a timezone-aware Timestamp object to a
        different time zone. The original UTC time remains the same; only the
        time zone information is changed. If the Timestamp is timezone-naive, a
        TypeError is raised.

        Parameters
        ----------
        tz : str, zoneinfo.ZoneInfo, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding UTC time.

        Returns
        -------
        converted : Timestamp

        Raises
        ------
        TypeError
            If Timestamp is tz-naive.

        See Also
        --------
        Timestamp.tz_localize : Localize the Timestamp to a timezone.
        DatetimeIndex.tz_convert : Convert a DatetimeIndex to another time zone.
        DatetimeIndex.tz_localize : Localize a DatetimeIndex to a specific time zone.
        datetime.datetime.astimezone : Convert a datetime object to another time zone.

        Examples
        --------
        Create a timestamp object with UTC timezone:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651', tz='UTC')
        >>> ts
        Timestamp('2020-03-14 15:32:52.192548651+0000', tz='UTC')

        Change to Tokyo timezone:

        >>> ts.tz_convert(tz='Asia/Tokyo')
        Timestamp('2020-03-15 00:32:52.192548651+0900', tz='Asia/Tokyo')

        Can also use ``astimezone``:

        >>> ts.astimezone(tz='Asia/Tokyo')
        Timestamp('2020-03-15 00:32:52.192548651+0900', tz='Asia/Tokyo')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.tz_convert(tz='Asia/Tokyo')
        NaT
        """,
    )
    fromordinal = _make_error_func(
        "fromordinal",
        """
        Construct a timestamp from a a proleptic Gregorian ordinal.

        Parameters
        ----------
        ordinal : int
            Date corresponding to a proleptic Gregorian ordinal.
        tz : str, zoneinfo.ZoneInfo, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for the Timestamp.

        Notes
        -----
        By definition there cannot be any tz info on the ordinal itself.

        Examples
        --------
        >>> pd.Timestamp.fromordinal(737425)
        Timestamp('2020-01-01 00:00:00')
        """,
    )

    # _nat_methods
    to_pydatetime = _make_nat_func(
        "to_pydatetime",
        """
        Convert a Timestamp object to a native Python datetime object.

        This method is useful for when you need to utilize a pandas Timestamp
        object in contexts where native Python datetime objects are expected
        or required. The conversion discards the nanoseconds component, and a
        warning can be issued in such cases if desired.

        Parameters
        ----------
        warn : bool, default True
            If True, issues a warning when the timestamp includes nonzero
            nanoseconds, as these will be discarded during the conversion.

        Returns
        -------
        datetime.datetime or NaT
            Returns a datetime.datetime object representing the timestamp,
            with year, month, day, hour, minute, second, and microsecond components.
            If the timestamp is NaT (Not a Time), returns NaT.

        See Also
        --------
        datetime.datetime : The standard Python datetime class that this method
            returns.
        Timestamp.timestamp : Convert a Timestamp object to POSIX timestamp.
        Timestamp.to_datetime64 : Convert a Timestamp object to numpy.datetime64.

        Examples
        --------
        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548')
        >>> ts.to_pydatetime()
        datetime.datetime(2020, 3, 14, 15, 32, 52, 192548)

        Analogous for ``pd.NaT``:

        >>> pd.NaT.to_pydatetime()
        NaT
        """,
    )

    now = _make_nat_func(
        "now",
        """
        Return new Timestamp object representing current time local to tz.

        This method returns a new `Timestamp` object that represents the current time.
        If a timezone is provided, the current time will be localized to that timezone.
        Otherwise, it returns the current local time.

        Parameters
        ----------
        tz : str or timezone object, default None
            Timezone to localize to.

        See Also
        --------
        to_datetime : Convert argument to datetime.
        Timestamp.utcnow : Return a new Timestamp representing UTC day and time.
        Timestamp.today : Return the current time in the local timezone.

        Examples
        --------
        >>> pd.Timestamp.now()  # doctest: +SKIP
        Timestamp('2020-11-16 22:06:16.378782')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.now()
        NaT
        """,
    )
    today = _make_nat_func(
        "today",
        """
        Return the current time in the local timezone.

        This differs from datetime.today() in that it can be localized to a
        passed timezone.

        Parameters
        ----------
        tz : str or timezone object, default None
            Timezone to localize to.

        Examples
        --------
        >>> pd.Timestamp.today()    # doctest: +SKIP
        Timestamp('2020-11-16 22:37:39.969883')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.today()
        NaT
        """,
    )
    round = _make_nat_func(
        "round",
        """
        Round the Timestamp to the specified resolution.

        This method rounds the given Timestamp down to a specified frequency
        level. It is particularly useful in data analysis to normalize timestamps
        to regular frequency intervals. For instance, rounding to the nearest
        minute, hour, or day can help in time series comparisons or resampling
        operations.

        Parameters
        ----------
        freq : str
            Frequency string indicating the rounding resolution.
        ambiguous : bool or {'raise', 'NaT'}, default 'raise'
            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise a ValueError for an ambiguous time.

        nonexistent : {'raise', 'shift_forward', 'shift_backward', 'NaT', \
timedelta}, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise a ValueError if there are
              nonexistent times.

        Returns
        -------
        a new Timestamp rounded to the given resolution of `freq`

        Raises
        ------
        ValueError if the freq cannot be converted

        See Also
        --------
        datetime.round : Similar behavior in native Python datetime module.
        Timestamp.floor : Round the Timestamp downward to the nearest multiple
            of the specified frequency.
        Timestamp.ceil : Round the Timestamp upward to the nearest multiple of
            the specified frequency.

        Notes
        -----
        If the Timestamp has a timezone, rounding will take place relative to the
        local ("wall") time and re-localized to the same timezone. When rounding
        near daylight savings time, use ``nonexistent`` and ``ambiguous`` to
        control the re-localization behavior.

        Examples
        --------
        Create a timestamp object:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')

        A timestamp can be rounded using multiple frequency units:

        >>> ts.round(freq='h')  # hour
        Timestamp('2020-03-14 16:00:00')

        >>> ts.round(freq='min')  # minute
        Timestamp('2020-03-14 15:33:00')

        >>> ts.round(freq='s')  # seconds
        Timestamp('2020-03-14 15:32:52')

        >>> ts.round(freq='ms')  # milliseconds
        Timestamp('2020-03-14 15:32:52.193000')

        ``freq`` can also be a multiple of a single unit, like '5min' (i.e.  5 minutes):

        >>> ts.round(freq='5min')
        Timestamp('2020-03-14 15:35:00')

        or a combination of multiple units, like '1h30min' (i.e. 1 hour and 30 minutes):

        >>> ts.round(freq='1h30min')
        Timestamp('2020-03-14 15:00:00')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.round()
        NaT

        When rounding near a daylight savings time transition, use ``ambiguous`` or
        ``nonexistent`` to control how the timestamp should be re-localized.

        >>> ts_tz = pd.Timestamp("2021-10-31 01:30:00").tz_localize("Europe/Amsterdam")

        >>> ts_tz.round("h", ambiguous=False)
        Timestamp('2021-10-31 02:00:00+0100', tz='Europe/Amsterdam')

        >>> ts_tz.round("h", ambiguous=True)
        Timestamp('2021-10-31 02:00:00+0200', tz='Europe/Amsterdam')
        """,
    )
    floor = _make_nat_func(
        "floor",
        """
        Return a new Timestamp floored to this resolution.

        Parameters
        ----------
        freq : str
            Frequency string indicating the flooring resolution.
        ambiguous : bool or {'raise', 'NaT'}, default 'raise'
            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise a ValueError for an ambiguous time.

        nonexistent : {'raise', 'shift_forward', 'shift_backward', 'NaT', \
timedelta}, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise a ValueError if there are
              nonexistent times.

        Raises
        ------
        ValueError if the freq cannot be converted.

        See Also
        --------
        Timestamp.ceil : Round up a Timestamp to the specified resolution.
        Timestamp.round : Round a Timestamp to the specified resolution.
        Series.dt.floor : Round down the datetime values in a Series.

        Notes
        -----
        If the Timestamp has a timezone, flooring will take place relative to the
        local ("wall") time and re-localized to the same timezone. When flooring
        near daylight savings time, use ``nonexistent`` and ``ambiguous`` to
        control the re-localization behavior.

        Examples
        --------
        Create a timestamp object:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')

        A timestamp can be floored using multiple frequency units:

        >>> ts.floor(freq='h')  # hour
        Timestamp('2020-03-14 15:00:00')

        >>> ts.floor(freq='min')  # minute
        Timestamp('2020-03-14 15:32:00')

        >>> ts.floor(freq='s')  # seconds
        Timestamp('2020-03-14 15:32:52')

        >>> ts.floor(freq='ns')  # nanoseconds
        Timestamp('2020-03-14 15:32:52.192548651')

        ``freq`` can also be a multiple of a single unit, like '5min' (i.e.  5 minutes):

        >>> ts.floor(freq='5min')
        Timestamp('2020-03-14 15:30:00')

        or a combination of multiple units, like '1h30min' (i.e. 1 hour and 30 minutes):

        >>> ts.floor(freq='1h30min')
        Timestamp('2020-03-14 15:00:00')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.floor()
        NaT

        When rounding near a daylight savings time transition, use ``ambiguous`` or
        ``nonexistent`` to control how the timestamp should be re-localized.

        >>> ts_tz = pd.Timestamp("2021-10-31 03:30:00").tz_localize("Europe/Amsterdam")

        >>> ts_tz.floor("2h", ambiguous=False)
        Timestamp('2021-10-31 02:00:00+0100', tz='Europe/Amsterdam')

        >>> ts_tz.floor("2h", ambiguous=True)
        Timestamp('2021-10-31 02:00:00+0200', tz='Europe/Amsterdam')
        """,
    )
    ceil = _make_nat_func(
        "ceil",
        """
        Return a new Timestamp ceiled to this resolution.

        Parameters
        ----------
        freq : str
            Frequency string indicating the ceiling resolution.
        ambiguous : bool or {'raise', 'NaT'}, default 'raise'
            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise a ValueError for an ambiguous time.

        nonexistent : {'raise', 'shift_forward', 'shift_backward', 'NaT', \
timedelta}, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise a ValueError if there are
              nonexistent times.

        Raises
        ------
        ValueError if the freq cannot be converted.

        See Also
        --------
        Timestamp.floor : Round down a Timestamp to the specified resolution.
        Timestamp.round : Round a Timestamp to the specified resolution.
        Series.dt.ceil : Ceil the datetime values in a Series.

        Notes
        -----
        If the Timestamp has a timezone, ceiling will take place relative to the
        local ("wall") time and re-localized to the same timezone. When ceiling
        near daylight savings time, use ``nonexistent`` and ``ambiguous`` to
        control the re-localization behavior.

        Examples
        --------
        Create a timestamp object:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')

        A timestamp can be ceiled using multiple frequency units:

        >>> ts.ceil(freq='h')  # hour
        Timestamp('2020-03-14 16:00:00')

        >>> ts.ceil(freq='min')  # minute
        Timestamp('2020-03-14 15:33:00')

        >>> ts.ceil(freq='s')  # seconds
        Timestamp('2020-03-14 15:32:53')

        >>> ts.ceil(freq='us')  # microseconds
        Timestamp('2020-03-14 15:32:52.192549')

        ``freq`` can also be a multiple of a single unit, like '5min' (i.e.  5 minutes):

        >>> ts.ceil(freq='5min')
        Timestamp('2020-03-14 15:35:00')

        or a combination of multiple units, like '1h30min' (i.e. 1 hour and 30 minutes):

        >>> ts.ceil(freq='1h30min')
        Timestamp('2020-03-14 16:30:00')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.ceil()
        NaT

        When rounding near a daylight savings time transition, use ``ambiguous`` or
        ``nonexistent`` to control how the timestamp should be re-localized.

        >>> ts_tz = pd.Timestamp("2021-10-31 01:30:00").tz_localize("Europe/Amsterdam")

        >>> ts_tz.ceil("h", ambiguous=False)
        Timestamp('2021-10-31 02:00:00+0100', tz='Europe/Amsterdam')

        >>> ts_tz.ceil("h", ambiguous=True)
        Timestamp('2021-10-31 02:00:00+0200', tz='Europe/Amsterdam')
        """,
    )

    tz_convert = _make_nat_func(
        "tz_convert",
        """
        Convert timezone-aware Timestamp to another time zone.

        This method is used to convert a timezone-aware Timestamp object to a
        different time zone. The original UTC time remains the same; only the
        time zone information is changed. If the Timestamp is timezone-naive, a
        TypeError is raised.

        Parameters
        ----------
        tz : str, zoneinfo.ZoneInfo, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding UTC time.

        Returns
        -------
        converted : Timestamp

        Raises
        ------
        TypeError
            If Timestamp is tz-naive.

        See Also
        --------
        Timestamp.tz_localize : Localize the Timestamp to a timezone.
        DatetimeIndex.tz_convert : Convert a DatetimeIndex to another time zone.
        DatetimeIndex.tz_localize : Localize a DatetimeIndex to a specific time zone.
        datetime.datetime.astimezone : Convert a datetime object to another time zone.

        Examples
        --------
        Create a timestamp object with UTC timezone:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651', tz='UTC')
        >>> ts
        Timestamp('2020-03-14 15:32:52.192548651+0000', tz='UTC')

        Change to Tokyo timezone:

        >>> ts.tz_convert(tz='Asia/Tokyo')
        Timestamp('2020-03-15 00:32:52.192548651+0900', tz='Asia/Tokyo')

        Can also use ``astimezone``:

        >>> ts.astimezone(tz='Asia/Tokyo')
        Timestamp('2020-03-15 00:32:52.192548651+0900', tz='Asia/Tokyo')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.tz_convert(tz='Asia/Tokyo')
        NaT
        """,
    )
    tz_localize = _make_nat_func(
        "tz_localize",
        """
        Localize the Timestamp to a timezone.

        Convert naive Timestamp to local time zone or remove
        timezone from timezone-aware Timestamp.

        Parameters
        ----------
        tz : str, zoneinfo.ZoneInfo, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding local time.

        ambiguous : bool, 'NaT', default 'raise'
            When clocks moved backward due to DST, ambiguous times may arise.
            For example in Central European Time (UTC+01), when going from
            03:00 DST to 02:00 non-DST, 02:30:00 local time occurs both at
            00:30:00 UTC and at 01:30:00 UTC. In such a situation, the
            `ambiguous` parameter dictates how ambiguous times should be
            handled.

            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise a ValueError for an ambiguous time.

        nonexistent : 'shift_forward', 'shift_backward', 'NaT', timedelta, \
default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            The behavior is as follows:

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise a ValueError if there are
              nonexistent times.

        Returns
        -------
        localized : Timestamp

        Raises
        ------
        TypeError
            If the Timestamp is tz-aware and tz is not None.

        See Also
        --------
        Timestamp.tzinfo : Returns the timezone information of the Timestamp.
        Timestamp.tz_convert : Convert timezone-aware Timestamp to another time zone.
        DatetimeIndex.tz_localize : Localize a DatetimeIndex to a specific time zone.
        datetime.datetime.astimezone : Convert a datetime object to another time zone.

        Examples
        --------
        Create a naive timestamp object:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651')
        >>> ts
        Timestamp('2020-03-14 15:32:52.192548651')

        Add 'Europe/Stockholm' as timezone:

        >>> ts.tz_localize(tz='Europe/Stockholm')
        Timestamp('2020-03-14 15:32:52.192548651+0100', tz='Europe/Stockholm')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.tz_localize()
        NaT
        """,
    )
    replace = _make_nat_func(
        "replace",
        """
        Implements datetime.replace, handles nanoseconds.

        Parameters
        ----------
        year : int, optional
        month : int, optional
        day : int, optional
        hour : int, optional
        minute : int, optional
        second : int, optional
        microsecond : int, optional
        nanosecond : int, optional
        tzinfo : tz-convertible, optional
        fold : int, optional

        Returns
        -------
        Timestamp with fields replaced

        Examples
        --------
        Create a timestamp object:

        >>> ts = pd.Timestamp('2020-03-14T15:32:52.192548651', tz='UTC')
        >>> ts
        Timestamp('2020-03-14 15:32:52.192548651+0000', tz='UTC')

        Replace year and the hour:

        >>> ts.replace(year=1999, hour=10)
        Timestamp('1999-03-14 10:32:52.192548651+0000', tz='UTC')

        Replace timezone (not a conversion):

        >>> import zoneinfo
        >>> ts.replace(tzinfo=zoneinfo.ZoneInfo('US/Pacific'))
        Timestamp('2020-03-14 15:32:52.192548651-0700', tz='US/Pacific')

        Analogous for ``pd.NaT``:

        >>> pd.NaT.replace(tzinfo=zoneinfo.ZoneInfo('US/Pacific'))
        NaT
        """,
    )

    @property
    def tz(self) -> None:
        return None

    @property
    def tzinfo(self) -> None:
        return None

    def as_unit(self, str unit, bint round_ok=True) -> "NaTType":
        """
        Convert the underlying int64 representaton to the given unit.

        Parameters
        ----------
        unit : {"ns", "us", "ms", "s"}
        round_ok : bool, default True
            If False and the conversion requires rounding, raise.

        Returns
        -------
        Timestamp

        See Also
        --------
        Timestamp.asm8 : Return numpy datetime64 format with same precision.
        Timestamp.to_pydatetime : Convert Timestamp object to a native
            Python datetime object.
        to_timedelta : Convert argument into timedelta object,
            which can represent differences in times.

        Examples
        --------
        >>> ts = pd.Timestamp('2023-01-01 00:00:00.01')
        >>> ts
        Timestamp('2023-01-01 00:00:00.010000')
        >>> ts.unit
        'ms'
        >>> ts = ts.as_unit('s')
        >>> ts
        Timestamp('2023-01-01 00:00:00')
        >>> ts.unit
        's'
        """
        return c_NaT


c_NaT = NaTType()  # C-visible
NaT = c_NaT        # Python-visible


# ----------------------------------------------------------------------

cdef bint checknull_with_nat(object val):
    """
    Utility to check if a value is a nat or not.
    """
    return val is None or util.is_nan(val) or val is c_NaT


cdef bint is_dt64nat(object val):
    """
    Is this a np.datetime64 object np.datetime64("NaT").
    """
    if cnp.is_datetime64_object(val):
        return cnp.get_datetime64_value(val) == NPY_NAT
    return False


cdef bint is_td64nat(object val):
    """
    Is this a np.timedelta64 object np.timedelta64("NaT").
    """
    if cnp.is_timedelta64_object(val):
        return cnp.get_timedelta64_value(val) == NPY_NAT
    return False
