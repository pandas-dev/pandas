from cpython.object cimport (
    PyObject_RichCompare,
    Py_GT, Py_GE, Py_EQ, Py_NE, Py_LT, Py_LE)

from cpython.datetime cimport (datetime, timedelta,
                               PyDateTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)

from cpython.version cimport PY_MINOR_VERSION

PyDateTime_IMPORT

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t
cnp.import_array()

from pandas._libs.tslibs.np_datetime cimport (
    get_datetime64_value, get_timedelta64_value)
cimport pandas._libs.tslibs.util as util
from pandas._libs.tslibs.util cimport (
    get_nat, is_integer_object, is_float_object, is_datetime64_object,
    is_timedelta64_object)


# ----------------------------------------------------------------------
# Constants
nat_strings = {'NaT', 'nat', 'NAT', 'nan', 'NaN', 'NAN'}

cdef int64_t NPY_NAT = get_nat()
iNaT = NPY_NAT  # python-visible constant

cdef bint _nat_scalar_rules[6]
_nat_scalar_rules[Py_EQ] = False
_nat_scalar_rules[Py_NE] = True
_nat_scalar_rules[Py_LT] = False
_nat_scalar_rules[Py_LE] = False
_nat_scalar_rules[Py_GT] = False
_nat_scalar_rules[Py_GE] = False

# ----------------------------------------------------------------------


def _make_nan_func(func_name, doc):
    def f(*args, **kwargs):
        return np.nan
    f.__name__ = func_name
    f.__doc__ = doc
    return f


def _make_nat_func(func_name, doc):
    def f(*args, **kwargs):
        return c_NaT
    f.__name__ = func_name
    f.__doc__ = doc
    return f


def _make_error_func(func_name, cls):
    def f(*args, **kwargs):
        raise ValueError("NaTType does not support " + func_name)

    f.__name__ = func_name
    if isinstance(cls, str):
        # passed the literal docstring directly
        f.__doc__ = cls
    elif cls is not None:
        f.__doc__ = getattr(cls, func_name).__doc__
    return f


cdef _nat_divide_op(self, other):
    if PyDelta_Check(other) or is_timedelta64_object(other) or other is c_NaT:
        return np.nan
    if is_integer_object(other) or is_float_object(other):
        return c_NaT
    return NotImplemented


cdef _nat_rdivide_op(self, other):
    if PyDelta_Check(other):
        return np.nan
    return NotImplemented


def __nat_unpickle(*args):
    # return constant defined in the module
    return c_NaT

# ----------------------------------------------------------------------


cdef class _NaT(datetime):
    # cdef readonly:
    #    int64_t value
    #    object freq

    # higher than np.ndarray and np.matrix
    __array_priority__ = 100

    def __richcmp__(_NaT self, object other, int op):
        cdef:
            int ndim = getattr(other, 'ndim', -1)

        if ndim == -1:
            return _nat_scalar_rules[op]

        elif util.is_array(other):
            result = np.empty(other.shape, dtype=np.bool_)
            result.fill(_nat_scalar_rules[op])
            return result

        elif ndim == 0:
            if is_datetime64_object(other):
                return _nat_scalar_rules[op]
            else:
                raise TypeError(f'Cannot compare type {type(self).__name__} '
                                f'with type {type(other).__name__}')

        # Note: instead of passing "other, self, _reverse_ops[op]", we observe
        # that `_nat_scalar_rules` is invariant under `_reverse_ops`,
        # rendering it unnecessary.
        return PyObject_RichCompare(other, self, op)

    def __add__(self, other):
        if self is not c_NaT:
            # cython __radd__ semantics
            self, other = other, self

        if PyDateTime_Check(other):
            return c_NaT
        elif PyDelta_Check(other):
            return c_NaT
        elif is_datetime64_object(other) or is_timedelta64_object(other):
            return c_NaT
        elif hasattr(other, 'delta'):
            # Timedelta, offsets.Tick, offsets.Week
            return c_NaT

        elif is_integer_object(other) or util.is_period_object(other):
            # For Period compat
            # TODO: the integer behavior is deprecated, remove it
            return c_NaT

        elif util.is_array(other):
            if other.dtype.kind in 'mM':
                # If we are adding to datetime64, we treat NaT as timedelta
                #  Either way, result dtype is datetime64
                result = np.empty(other.shape, dtype="datetime64[ns]")
                result.fill("NaT")
                return result
            raise TypeError(f"Cannot add NaT to ndarray with dtype {other.dtype}")

        return NotImplemented

    def __sub__(self, other):
        # Duplicate some logic from _Timestamp.__sub__ to avoid needing
        # to subclass; allows us to @final(_Timestamp.__sub__)
        cdef:
            bint is_rsub = False

        if self is not c_NaT:
            # cython __rsub__ semantics
            self, other = other, self
            is_rsub = True

        if PyDateTime_Check(other):
            return c_NaT
        elif PyDelta_Check(other):
            return c_NaT
        elif is_datetime64_object(other) or is_timedelta64_object(other):
            return c_NaT
        elif hasattr(other, 'delta'):
            # offsets.Tick, offsets.Week
            return c_NaT

        elif is_integer_object(other) or util.is_period_object(other):
            # For Period compat
            # TODO: the integer behavior is deprecated, remove it
            return c_NaT

        elif util.is_array(other):
            if other.dtype.kind == 'm':
                if not is_rsub:
                    # NaT - timedelta64 we treat NaT as datetime64, so result
                    #  is datetime64
                    result = np.empty(other.shape, dtype="datetime64[ns]")
                    result.fill("NaT")
                    return result

                # timedelta64 - NaT we have to treat NaT as timedelta64
                #  for this to be meaningful, and the result is timedelta64
                result = np.empty(other.shape, dtype="timedelta64[ns]")
                result.fill("NaT")
                return result

            elif other.dtype.kind == 'M':
                # We treat NaT as a datetime, so regardless of whether this is
                #  NaT - other or other - NaT, the result is timedelta64
                result = np.empty(other.shape, dtype="timedelta64[ns]")
                result.fill("NaT")
                return result

            raise TypeError(f"Cannot subtract NaT from ndarray with "
                            f"dtype {other.dtype}")

        return NotImplemented

    def __pos__(self):
        return NaT

    def __neg__(self):
        return NaT

    def __div__(self, other):
        return _nat_divide_op(self, other)

    def __truediv__(self, other):
        return _nat_divide_op(self, other)

    def __floordiv__(self, other):
        return _nat_divide_op(self, other)

    def __mul__(self, other):
        if is_integer_object(other) or is_float_object(other):
            return NaT
        return NotImplemented

    @property
    def asm8(self) -> np.datetime64:
        return np.datetime64(NPY_NAT, 'ns')

    def to_datetime64(self) -> np.datetime64:
        """
        Return a numpy.datetime64 object with 'ns' precision.
        """
        return np.datetime64('NaT', 'ns')

    def to_numpy(self, dtype=None, copy=False) -> np.datetime64:
        """
        Convert the Timestamp to a NumPy datetime64.

        .. versionadded:: 0.25.0

        This is an alias method for `Timestamp.to_datetime64()`. The dtype and
        copy parameters are available here only for compatibility. Their values
        will not affect the return value.

        Returns
        -------
        numpy.datetime64

        See Also
        --------
        DatetimeIndex.to_numpy : Similar method for DatetimeIndex.
        """
        return self.to_datetime64()

    def __repr__(self) -> str:
        return 'NaT'

    def __str__(self) -> str:
        return 'NaT'

    def isoformat(self, sep='T') -> str:
        # This allows Timestamp(ts.isoformat()) to always correctly roundtrip.
        return 'NaT'

    def __hash__(self):
        return NPY_NAT

    def __int__(self):
        return NPY_NAT

    def __long__(self):
        return NPY_NAT

    @property
    def is_leap_year(self):
        return False

    @property
    def is_month_start(self):
        return False

    @property
    def is_quarter_start(self):
        return False

    @property
    def is_year_start(self):
        return False

    @property
    def is_month_end(self):
        return False

    @property
    def is_quarter_end(self):
        return False

    @property
    def is_year_end(self):
        return False


class NaTType(_NaT):
    """(N)ot-(A)-(T)ime, the time equivalent of NaN"""

    def __new__(cls):
        cdef _NaT base

        base = _NaT.__new__(cls, 1, 1, 1)
        base.value = NPY_NAT
        base.freq = None

        return base

    def __reduce_ex__(self, protocol):
        # python 3.6 compat
        # https://bugs.python.org/issue28730
        # now __reduce_ex__ is defined and higher priority than __reduce__
        return self.__reduce__()

    def __reduce__(self):
        return (__nat_unpickle, (None, ))

    def __rdiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rtruediv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rfloordiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rmul__(self, other):
        if is_integer_object(other) or is_float_object(other):
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
    weekofyear = property(fget=lambda self: np.nan)
    days_in_month = property(fget=lambda self: np.nan)
    daysinmonth = property(fget=lambda self: np.nan)
    dayofweek = property(fget=lambda self: np.nan)

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
    weekday = _make_nan_func('weekday', datetime.weekday.__doc__)
    isoweekday = _make_nan_func('isoweekday', datetime.isoweekday.__doc__)
    total_seconds = _make_nan_func('total_seconds', timedelta.total_seconds.__doc__)
    month_name = _make_nan_func('month_name',  # noqa:E128
        """
        Return the month name of the Timestamp with specified locale.

        Parameters
        ----------
        locale : string, default None (English locale)
            Locale determining the language in which to return the month name.

        Returns
        -------
        month_name : string

        .. versionadded:: 0.23.0
        """)
    day_name = _make_nan_func('day_name',  # noqa:E128
        """
        Return the day name of the Timestamp with specified locale.

        Parameters
        ----------
        locale : string, default None (English locale)
            Locale determining the language in which to return the day name.

        Returns
        -------
        day_name : string

        .. versionadded:: 0.23.0
        """)
    # _nat_methods
    date = _make_nat_func('date', datetime.date.__doc__)

    utctimetuple = _make_error_func('utctimetuple', datetime)
    timetz = _make_error_func('timetz', datetime)
    timetuple = _make_error_func('timetuple', datetime)
    strftime = _make_error_func('strftime', datetime)
    isocalendar = _make_error_func('isocalendar', datetime)
    dst = _make_error_func('dst', datetime)
    ctime = _make_error_func('ctime', datetime)
    time = _make_error_func('time', datetime)
    toordinal = _make_error_func('toordinal', datetime)
    tzname = _make_error_func('tzname', datetime)
    utcoffset = _make_error_func('utcoffset', datetime)

    # "fromisocalendar" was introduced in 3.8
    if PY_MINOR_VERSION >= 8:
        fromisocalendar = _make_error_func('fromisocalendar', datetime)

    # ----------------------------------------------------------------------
    # The remaining methods have docstrings copy/pasted from the analogous
    # Timestamp methods.

    strptime = _make_error_func('strptime',  # noqa:E128
        """
        Timestamp.strptime(string, format)

        Function is not implemented. Use pd.to_datetime().
        """
    )

    utcfromtimestamp = _make_error_func('utcfromtimestamp',  # noqa:E128
        """
        Timestamp.utcfromtimestamp(ts)

        Construct a naive UTC datetime from a POSIX timestamp.
        """
    )
    fromtimestamp = _make_error_func('fromtimestamp',  # noqa:E128
        """
        Timestamp.fromtimestamp(ts)

        timestamp[, tz] -> tz's local time from POSIX timestamp.
        """
    )
    combine = _make_error_func('combine',  # noqa:E128
        """
        Timestamp.combine(date, time)

        date, time -> datetime with same date and time fields.
        """
    )
    utcnow = _make_error_func('utcnow',  # noqa:E128
        """
        Timestamp.utcnow()

        Return a new Timestamp representing UTC day and time.
        """
    )

    timestamp = _make_error_func('timestamp',  # noqa:E128
        """Return POSIX timestamp as float.""")

    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    astimezone = _make_error_func('astimezone',  # noqa:E128
        """
        Convert tz-aware Timestamp to another time zone.

        Parameters
        ----------
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding UTC time.

        Returns
        -------
        converted : Timestamp

        Raises
        ------
        TypeError
            If Timestamp is tz-naive.
        """)
    fromordinal = _make_error_func('fromordinal',  # noqa:E128
        """
        Timestamp.fromordinal(ordinal, freq=None, tz=None)

        Passed an ordinal, translate and convert to a ts.
        Note: by definition there cannot be any tz info on the ordinal itself.

        Parameters
        ----------
        ordinal : int
            Date corresponding to a proleptic Gregorian ordinal.
        freq : str, DateOffset
            Offset to apply to the Timestamp.
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for the Timestamp.
        """)

    # _nat_methods
    to_pydatetime = _make_nat_func('to_pydatetime',  # noqa:E128
        """
        Convert a Timestamp object to a native Python datetime object.

        If warn=True, issue a warning if nanoseconds is nonzero.
        """)

    now = _make_nat_func('now',  # noqa:E128
        """
        Timestamp.now(tz=None)

        Return new Timestamp object representing current time local to
        tz.

        Parameters
        ----------
        tz : str or timezone object, default None
            Timezone to localize to.
        """)
    today = _make_nat_func('today',  # noqa:E128
        """
        Timestamp.today(cls, tz=None)

        Return the current time in the local timezone.  This differs
        from datetime.today() in that it can be localized to a
        passed timezone.

        Parameters
        ----------
        tz : str or timezone object, default None
            Timezone to localize to.
        """)
    round = _make_nat_func('round',  # noqa:E128
        """
        Round the Timestamp to the specified resolution.

        Parameters
        ----------
        freq : str
            Frequency string indicating the rounding resolution.
        ambiguous : bool or {'raise', 'NaT'}, default 'raise'
            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise an AmbiguousTimeError for an ambiguous time.

            .. versionadded:: 0.24.0
        nonexistent : {'raise', 'shift_forward', 'shift_backward, 'NaT', \
timedelta}, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise an NonExistentTimeError if there are
              nonexistent times.

            .. versionadded:: 0.24.0

        Returns
        -------
        a new Timestamp rounded to the given resolution of `freq`

        Raises
        ------
        ValueError if the freq cannot be converted
        """)
    floor = _make_nat_func('floor',  # noqa:E128
        """
        return a new Timestamp floored to this resolution.

        Parameters
        ----------
        freq : str
            Frequency string indicating the flooring resolution.
        ambiguous : bool or {'raise', 'NaT'}, default 'raise'
            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise an AmbiguousTimeError for an ambiguous time.

            .. versionadded:: 0.24.0
        nonexistent : {'raise', 'shift_forward', 'shift_backward, 'NaT', \
timedelta}, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise an NonExistentTimeError if there are
              nonexistent times.

            .. versionadded:: 0.24.0

        Raises
        ------
        ValueError if the freq cannot be converted.
        """)
    ceil = _make_nat_func('ceil',  # noqa:E128
        """
        return a new Timestamp ceiled to this resolution.

        Parameters
        ----------
        freq : str
            Frequency string indicating the ceiling resolution.
        ambiguous : bool or {'raise', 'NaT'}, default 'raise'
            The behavior is as follows:

            * bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates).
            * 'NaT' will return NaT for an ambiguous time.
            * 'raise' will raise an AmbiguousTimeError for an ambiguous time.

            .. versionadded:: 0.24.0
        nonexistent : {'raise', 'shift_forward', 'shift_backward, 'NaT', \
timedelta}, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            * 'shift_forward' will shift the nonexistent time forward to the
              closest existing time.
            * 'shift_backward' will shift the nonexistent time backward to the
              closest existing time.
            * 'NaT' will return NaT where there are nonexistent times.
            * timedelta objects will shift nonexistent times by the timedelta.
            * 'raise' will raise an NonExistentTimeError if there are
              nonexistent times.

            .. versionadded:: 0.24.0

        Raises
        ------
        ValueError if the freq cannot be converted.
        """)

    tz_convert = _make_nat_func('tz_convert',  # noqa:E128
        """
        Convert tz-aware Timestamp to another time zone.

        Parameters
        ----------
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding UTC time.

        Returns
        -------
        converted : Timestamp

        Raises
        ------
        TypeError
            If Timestamp is tz-naive.
        """)
    tz_localize = _make_nat_func('tz_localize',  # noqa:E128
        """
        Convert naive Timestamp to local time zone, or remove
        timezone from tz-aware Timestamp.

        Parameters
        ----------
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
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
            * 'raise' will raise an AmbiguousTimeError for an ambiguous time.

        nonexistent : 'shift_forward', 'shift_backward, 'NaT', timedelta, \
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
            * 'raise' will raise an NonExistentTimeError if there are
              nonexistent times.

            .. versionadded:: 0.24.0

        Returns
        -------
        localized : Timestamp

        Raises
        ------
        TypeError
            If the Timestamp is tz-aware and tz is not None.
        """)
    replace = _make_nat_func('replace',  # noqa:E128
        """
        implements datetime.replace, handles nanoseconds.

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
        fold : int, optional, default is 0

        Returns
        -------
        Timestamp with fields replaced
        """)


c_NaT = NaTType()  # C-visible
NaT = c_NaT        # Python-visible


# ----------------------------------------------------------------------

cdef inline bint checknull_with_nat(object val):
    """ utility to check if a value is a nat or not """
    return val is None or util.is_nan(val) or val is c_NaT


cpdef bint is_null_datetimelike(object val, bint inat_is_null=True):
    """
    Determine if we have a null for a timedelta/datetime (or integer versions)

    Parameters
    ----------
    val : object
    inat_is_null : bool, default True
        Whether to treat integer iNaT value as null

    Returns
    -------
    null_datetimelike : bool
    """
    if val is None:
        return True
    elif val is c_NaT:
        return True
    elif util.is_float_object(val) or util.is_complex_object(val):
        return val != val
    elif util.is_timedelta64_object(val):
        return get_timedelta64_value(val) == NPY_NAT
    elif util.is_datetime64_object(val):
        return get_datetime64_value(val) == NPY_NAT
    elif inat_is_null and util.is_integer_object(val):
        return val == NPY_NAT
    return False
