"""
_Timestamp is a c-defined subclass of datetime.datetime

It is separate from timestamps.pyx to prevent circular cimports

This allows _Timestamp to be imported in other modules
so that isinstance(obj, _Timestamp) checks can be performed

_Timestamp is PITA. Because we inherit from datetime, which has very specific
construction requirements, we need to do object instantiation in python
(see Timestamp class below). This will serve as a C extension type that
shadows the python class, where we do any heavy lifting.
"""

import warnings

from cpython.object cimport (PyObject_RichCompareBool, PyObject_RichCompare,
                             Py_GT, Py_GE, Py_EQ, Py_NE, Py_LT, Py_LE)

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t, int8_t, uint8_t, ndarray
cnp.import_array()

from cpython.datetime cimport (datetime,
                               PyDateTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT

from pandas._libs.tslibs.util cimport (
    is_datetime64_object, is_timedelta64_object, is_integer_object,
    is_array)

from pandas._libs.tslibs.fields import get_start_end_field, get_date_name_field
from pandas._libs.tslibs.nattype cimport c_NaT as NaT
from pandas._libs.tslibs.np_datetime import OutOfBoundsDatetime
from pandas._libs.tslibs.np_datetime cimport (
    reverse_ops, cmp_scalar)
from pandas._libs.tslibs.timezones cimport (
    get_timezone, is_utc, tz_compare)
from pandas._libs.tslibs.timezones import UTC
from pandas._libs.tslibs.tzconversion cimport tz_convert_single


class NullFrequencyError(ValueError):
    """
    Error raised when a null `freq` attribute is used in an operation
    that needs a non-null frequency, particularly `DatetimeIndex.shift`,
    `TimedeltaIndex.shift`, `PeriodIndex.shift`.
    """
    pass


def integer_op_not_supported(obj):
    # GH#22535 add/sub of integers and int-arrays is no longer allowed
    # Note we return rather than raise the exception so we can raise in
    #  the caller; mypy finds this more palatable.
    cls = type(obj).__name__

    # GH#30886 using an fstring raises SystemError
    int_addsub_msg = (
        "Addition/subtraction of integers and integer-arrays with {cls} is "
        "no longer supported.  Instead of adding/subtracting `n`, "
        "use `n * obj.freq`"
    ).format(cls=cls)
    return TypeError(int_addsub_msg)


cdef class _Timestamp(datetime):

    # higher than np.ndarray and np.matrix
    __array_priority__ = 100

    def __hash__(_Timestamp self):
        if self.nanosecond:
            return hash(self.value)
        return datetime.__hash__(self)

    def __richcmp__(_Timestamp self, object other, int op):
        cdef:
            _Timestamp ots
            int ndim

        if isinstance(other, _Timestamp):
            ots = other
        elif other is NaT:
            return op == Py_NE
        elif PyDateTime_Check(other):
            if self.nanosecond == 0:
                val = self.to_pydatetime()
                return PyObject_RichCompareBool(val, other, op)

            try:
                ots = type(self)(other)
            except ValueError:
                return self._compare_outside_nanorange(other, op)
        else:
            ndim = getattr(other, "ndim", -1)

            if ndim != -1:
                if ndim == 0:
                    if is_datetime64_object(other):
                        other = type(self)(other)
                    elif is_array(other):
                        # zero-dim array, occurs if try comparison with
                        #  datetime64 scalar on the left hand side
                        # Unfortunately, for datetime64 values, other.item()
                        #  incorrectly returns an integer, so we need to use
                        #  the numpy C api to extract it.
                        other = cnp.PyArray_ToScalar(cnp.PyArray_DATA(other),
                                                     other)
                        other = type(self)(other)
                    else:
                        return NotImplemented
                elif is_array(other):
                    # avoid recursion error GH#15183
                    return PyObject_RichCompare(np.array([self]), other, op)
                return PyObject_RichCompare(other, self, reverse_ops[op])
            else:
                return NotImplemented

        self._assert_tzawareness_compat(other)
        return cmp_scalar(self.value, ots.value, op)

    def __reduce_ex__(self, protocol):
        # python 3.6 compat
        # http://bugs.python.org/issue28730
        # now __reduce_ex__ is defined and higher priority than __reduce__
        return self.__reduce__()

    def __repr__(self) -> str:
        stamp = self._repr_base
        zone = None

        try:
            stamp += self.strftime('%z')
            if self.tzinfo:
                zone = get_timezone(self.tzinfo)
        except ValueError:
            year2000 = self.replace(year=2000)
            stamp += year2000.strftime('%z')
            if self.tzinfo:
                zone = get_timezone(self.tzinfo)

        try:
            stamp += zone.strftime(' %%Z')
        except AttributeError:
            # e.g. tzlocal has no `strftime`
            pass

        tz = f", tz='{zone}'" if zone is not None else ""
        freq = "" if self.freq is None else f", freq='{self.freqstr}'"

        return f"Timestamp('{stamp}'{tz}{freq})"

    cdef bint _compare_outside_nanorange(_Timestamp self, datetime other,
                                         int op) except -1:
        cdef:
            datetime dtval = self.to_pydatetime()

        self._assert_tzawareness_compat(other)

        if self.nanosecond == 0:
            return PyObject_RichCompareBool(dtval, other, op)
        else:
            if op == Py_EQ:
                return False
            elif op == Py_NE:
                return True
            elif op == Py_LT:
                return dtval < other
            elif op == Py_LE:
                return dtval < other
            elif op == Py_GT:
                return dtval >= other
            elif op == Py_GE:
                return dtval >= other

    cdef _assert_tzawareness_compat(_Timestamp self, datetime other):
        if self.tzinfo is None:
            if other.tzinfo is not None:
                raise TypeError('Cannot compare tz-naive and tz-aware '
                                'timestamps')
        elif other.tzinfo is None:
            raise TypeError('Cannot compare tz-naive and tz-aware timestamps')

    cpdef datetime to_pydatetime(_Timestamp self, bint warn=True):
        """
        Convert a Timestamp object to a native Python datetime object.

        If warn=True, issue a warning if nanoseconds is nonzero.
        """
        if self.nanosecond != 0 and warn:
            warnings.warn("Discarding nonzero nanoseconds in conversion",
                          UserWarning, stacklevel=2)

        return datetime(self.year, self.month, self.day,
                        self.hour, self.minute, self.second,
                        self.microsecond, self.tzinfo)

    cpdef to_datetime64(self):
        """
        Return a numpy.datetime64 object with 'ns' precision.
        """
        return np.datetime64(self.value, 'ns')

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

    def __add__(self, other):
        cdef:
            int64_t other_int, nanos = 0

        if is_timedelta64_object(other):
            other_int = other.astype('timedelta64[ns]').view('i8')
            return type(self)(self.value + other_int, tz=self.tzinfo, freq=self.freq)

        elif is_integer_object(other):
            raise integer_op_not_supported(self)

        elif PyDelta_Check(other) or hasattr(other, 'delta'):
            # delta --> offsets.Tick
            # logic copied from delta_to_nanoseconds to prevent circular import
            if hasattr(other, 'nanos'):
                nanos = other.nanos
            elif hasattr(other, 'delta'):
                nanos = other.delta
            elif PyDelta_Check(other):
                nanos = (other.days * 24 * 60 * 60 * 1000000 +
                         other.seconds * 1000000 +
                         other.microseconds) * 1000

            result = type(self)(self.value + nanos, tz=self.tzinfo, freq=self.freq)
            return result

        elif is_array(other):
            if other.dtype.kind in ['i', 'u']:
                raise integer_op_not_supported(self)

        # index/series like
        elif hasattr(other, '_typ'):
            return NotImplemented

        result = datetime.__add__(self, other)
        if PyDateTime_Check(result):
            result = type(self)(result)
            result.nanosecond = self.nanosecond
        return result

    def __sub__(self, other):

        if (is_timedelta64_object(other) or is_integer_object(other) or
                PyDelta_Check(other) or hasattr(other, 'delta')):
            # `delta` attribute is for offsets.Tick or offsets.Week obj
            neg_other = -other
            return self + neg_other

        elif is_array(other):
            if other.dtype.kind in ['i', 'u']:
                raise integer_op_not_supported(self)

        typ = getattr(other, '_typ', None)
        if typ is not None:
            return NotImplemented

        if other is NaT:
            return NaT

        # coerce if necessary if we are a Timestamp-like
        if (PyDateTime_Check(self)
                and (PyDateTime_Check(other) or is_datetime64_object(other))):
            if isinstance(self, _Timestamp):
                other = type(self)(other)
            else:
                self = type(other)(self)

            # validate tz's
            if not tz_compare(self.tzinfo, other.tzinfo):
                raise TypeError("Timestamp subtraction must have the "
                                "same timezones or no timezones")

            # scalar Timestamp/datetime - Timestamp/datetime -> yields a
            # Timedelta
            from pandas._libs.tslibs.timedeltas import Timedelta
            try:
                return Timedelta(self.value - other.value)
            except (OverflowError, OutOfBoundsDatetime):
                pass
        elif is_datetime64_object(self):
            # GH#28286 cython semantics for __rsub__, `other` is actually
            #  the Timestamp
            return type(other)(self) - other

        return NotImplemented

    cdef int64_t _maybe_convert_value_to_local(self):
        """Convert UTC i8 value to local i8 value if tz exists"""
        cdef:
            int64_t val
        val = self.value
        if self.tz is not None and not is_utc(self.tz):
            val = tz_convert_single(self.value, UTC, self.tz)
        return val

    cpdef bint _get_start_end_field(self, str field):
        cdef:
            int64_t val
            dict kwds
            ndarray[uint8_t, cast=True] out
            int month_kw

        freq = self.freq
        if freq:
            kwds = freq.kwds
            month_kw = kwds.get('startingMonth', kwds.get('month', 12))
            freqstr = self.freqstr
        else:
            month_kw = 12
            freqstr = None

        val = self._maybe_convert_value_to_local()
        out = get_start_end_field(np.array([val], dtype=np.int64),
                                  field, freqstr, month_kw)
        return out[0]

    cpdef _get_date_name_field(self, object field, object locale):
        cdef:
            int64_t val
            object[:] out

        val = self._maybe_convert_value_to_local()
        out = get_date_name_field(np.array([val], dtype=np.int64),
                                  field, locale=locale)
        return out[0]

    @property
    def _repr_base(self) -> str:
        return f"{self._date_repr} {self._time_repr}"

    @property
    def _date_repr(self) -> str:
        # Ideal here would be self.strftime("%Y-%m-%d"), but
        # the datetime strftime() methods require year >= 1900
        return f'{self.year}-{self.month:02d}-{self.day:02d}'

    @property
    def _time_repr(self) -> str:
        result = f'{self.hour:02d}:{self.minute:02d}:{self.second:02d}'

        if self.nanosecond != 0:
            result += f'.{self.nanosecond + 1000 * self.microsecond:09d}'
        elif self.microsecond != 0:
            result += f'.{self.microsecond:06d}'

        return result

    @property
    def _short_repr(self) -> str:
        # format a Timestamp with only _date_repr if possible
        # otherwise _repr_base
        if (self.hour == 0 and
                self.minute == 0 and
                self.second == 0 and
                self.microsecond == 0 and
                self.nanosecond == 0):
            return self._date_repr
        return self._repr_base

    @property
    def asm8(self) -> np.datetime64:
        """
        Return numpy datetime64 format in nanoseconds.
        """
        return np.datetime64(self.value, 'ns')

    def timestamp(self):
        """Return POSIX timestamp as float."""
        # GH 17329
        # Note: Naive timestamps will not match datetime.stdlib
        return round(self.value / 1e9, 6)
