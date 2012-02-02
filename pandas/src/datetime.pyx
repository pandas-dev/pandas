cimport numpy as np
import numpy as np

from numpy cimport int32_t, int64_t, import_array, ndarray
from cpython cimport *

from libc.stdlib cimport malloc, free
from libc.math cimport floor

# this is our datetime.pxd
from datetime cimport *
from util cimport is_integer_object

# initialize numpy
np.import_array()
np.import_ufunc()

# import datetime C API
PyDateTime_IMPORT

# in numpy 1.7, will prob need this
# numpy_pydatetime_import

ctypedef enum time_res:
    r_min = 0
    r_microsecond
    r_second
    r_minute
    r_hour
    r_day
    r_month
    r_year
    r_max = 98
    r_invalid = 99

# Objects to support date/time arithmetic
# --------------------------------------------------------------------------------

cdef class Timestamp:
    """
    A timestamp (absolute moment in time) to microsecond resolution; number of
    microseconds since the POSIX epoch, ignoring leap seconds (thereby different
    from UTC).
    """
    cdef:
        int64_t value
        npy_datetimestruct dts

    def __init__(self, object ts):
        """
        Construct a timestamp that is datetime64-compatible from any of:
            - int64 pyarray scalar object
            - python int or long object
            - iso8601 string object
            - python datetime object
            - another timestamp object
        """
        cdef:
            Py_ssize_t strlen
            npy_bool islocal, special
            NPY_DATETIMEUNIT out_bestunit

        if is_integer_object(ts) or PyInt_Check(ts) or PyLong_Check(ts):
            self.value = ts
            PyArray_DatetimeToDatetimeStruct(self.value, NPY_FR_us, &self.dts)
        elif PyString_Check(ts):
            parse_iso_8601_datetime(ts, len(ts), NPY_FR_us, NPY_UNSAFE_CASTING,
                                    &self.dts, &islocal, &out_bestunit, &special)
            self.value = PyArray_DatetimeStructToDatetime(NPY_FR_us, &self.dts)
        elif PyDateTime_Check(ts):
            convert_pydatetime_to_datetimestruct(<PyObject *>ts, &self.dts,
                                                 &out_bestunit, 1)
            self.value = PyArray_DatetimeStructToDatetime(out_bestunit, &self.dts)
        else:
            raise ValueError("Could not construct Timestamp from argument")

    def __repr__(self):
        return "Timestamp(%d)" % self.value

    def __hash__(self):
        return hash(self.value)

    def __sub__(self, object other):
        """
        Subtract two timestamps, results in an duration.
        """
        if isinstance(other, Timestamp):
            return self.value - other.asint
        elif isinstance(other, Delta):
            return other.__sub__(self)
        else:
            raise NotImplementedError("Sub operation not supported")

    def __richcmp__(Timestamp self, object other, int op):
        if not isinstance(other, Timestamp):
            raise ValueError("Cannot compare to non-Timestamp")

        if op == 0:
            return self.value < other.asint
        if op == 2:
            return self.value == other.asint
        if op == 4:
            return self.value > other.asint
        if op == 1:
            return self.value <= other.asint
        if op == 3:
            return self.value != other.asint
        if op == 5:
            return self.value >= other.asint

        raise NotImplementedError("Op %d not recognized" % op)

    def __add__(self, object other):
        """
        Add an Interval, Duration, or Period to the Timestamp, resulting in
        new Timestamp.
        """
        if isinstance(other, Duration):
            return Timestamp(self.asint + other.length)
        elif isinstance(other, Delta):
            return other.__add__(self)
        else:
            raise NotImplementedError("Add operation not supported")

    def __str__(self):
        """
        Output ISO8601 format string representation of timestamp.
        """
        cdef:
            int outlen
            char *isostr
            bytes py_str

        outlen = get_datetime_iso_8601_strlen(0, NPY_FR_us)

        isostr = <char *>malloc(outlen)
        make_iso_8601_datetime(&self.dts, isostr, outlen, 0, NPY_FR_us,
                               0, NPY_UNSAFE_CASTING)
        py_str = isostr
        free(isostr)

        return py_str

    def replace(Timestamp self, int year=-1, int month=-1, int day=-1,
                                int hour=-1, int minute=-1, int second=-1,
                                int microsecond=-1):
        cdef:
            npy_datetimestruct dts

        dts = self.dts

        if year >= 0:
            dts.year = year
        if month >= 1:
            dts.month = month
        if day >= 1:
            dts.day = day
        if hour >= 0:
            dts.hour = hour
        if minute >= 0:
            dts.min = minute
        if second >= 0:
            dts.sec = second
        if microsecond >= 0:
            dts.us = microsecond

        return Timestamp(PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts))

    cdef normalize(Timestamp self, time_res res):
        cdef:
            npy_datetimestruct dts

        dts = self.dts

        if res > r_microsecond:
            dts.us = 0
        if res > r_second:
            dts.sec = 0
        if res > r_minute:
            dts.min = 0
        if res > r_hour:
            dts.hour = 0
        if res > r_day:
            dts.day = 1
        if res > r_month:
            dts.month = 1
        if res > r_year:
            raise ValueError("Invalid resolution")

        return Timestamp(PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts))

    property asint:
        def __get__(Timestamp self):
            return self.value

    property year:
        def __get__(Timestamp self):
            return self.dts.year

    property month:
        def __get__(Timestamp self):
            return self.dts.month

    property day:
        def __get__(Timestamp self):
            return self.dts.day

    property hour:
        def __get__(Timestamp self):
            return self.dts.hour

    property minute:
        def __get__(Timestamp self):
            return self.dts.min

    property second:
        def __get__(Timestamp self):
            return self.dts.sec

    property microsecond:
        def __get__(Timestamp self):
            return self.dts.us

    def weekday(Timestamp self):
        return dayofweek(self.dts.year, self.dts.month, self.dts.day)


cdef class Duration:
    """
    Absolute length of time, similar to timedelta
    """
    cdef:
        int64_t length, days, seconds, microseconds

    def __init__(self, int64_t days = 0,
                       int64_t seconds = 0,
                       int64_t microseconds = 0,
                       int64_t milliseconds = 0,
                       int64_t minutes = 0,
                       int64_t hours = 0,
                       int64_t weeks = 0):

        self.days = days
        self.seconds = seconds
        self.microseconds = microseconds
        self.length = (microseconds + 1000 * (milliseconds
                                    + 1000 * (seconds
                                    + 60   * (minutes
                                    + 60   * (hours
                                    + 24   * (days
                                    +  7   * weeks))))))

    @staticmethod
    def from_micros(int64_t length):
        return Duration(microseconds = length)

    property length:
        def __get__(self):
            return self.length

    property microseconds:
        def __get__(self):
            return self.microseconds # length % 1000000

    property seconds:
        def __get__(self):
            return self.seconds      # (length // 1000000) % 86400

    property days:
        def __get__(self):
            return self.days         # (self.length // 1000000) // 86400

    def __repr__(self):
        return "Duration(%d, %d, %d)" % (self.days, self.seconds, self.microseconds)


cdef convert_to_res(object res):
    if res == 'microsecond':
        return r_microsecond
    if res == 'second':
        return r_second
    if res == 'minute':
        return r_minute
    if res == 'hour':
        return r_hour
    if res == 'day':
        return r_day
    if res == 'month':
        return r_month
    if res == 'year':
        return r_year
    return r_invalid

cdef conversion_factor(time_res res1, time_res res2):
    cdef:
        time_res min_res, max_res
        int64_t factor

    min_res = min(res1, res2)
    max_res = max(res1, res2)
    factor = 1

    if min_res == max_res:
        return factor

    while min_res < max_res:
        if min_res < r_microsecond:
            raise "Cannot convert from less than us"
        elif min_res == r_microsecond:
            factor *= 1000000
            min_res = r_second
        elif min_res == r_second:
            factor *= 60
            min_res = r_minute
        elif min_res == r_minute:
            factor *= 60
            min_res = r_hour
        elif min_res == r_hour:
            factor *= 24
            min_res = r_day
        else:
            raise "Cannot convert to month or year"

    return factor

# The following is derived from relativedelta.py in dateutil package
# ------------------------------------------------------------------------------
# Copyright (c) 2003-2010  Gustavo Niemeyer <gustavo@niemeyer.net>
# under Simplified BSD

cdef class Weekday:
    cdef:
        int64_t weekday, n

    def __init__(self, int64_t weekday, int64_t n = INT64_MIN):
        if weekday < 0 or weekday > 6:
            raise ValueError("Invalid weekday: %d", weekday)

        self.weekday = weekday
        self.n = n

    def __call__(self, int n):
        if n == self.n:
            return self
        else:
            return self.__class__(self.weekday, n)

    def __richcmp__(self, other, int op):
        isequal = False

        if not isinstance(other, Weekday):
            isequal = False
        else:
            isequal = (self.weekday == other.weekday and self.n == other.n)

        if op == 2: # equals
            return isequal
        if op == 3: # not equals
            return not isequal

        raise NotImplementedError("Comparison not supported")

    property weekday:
        def __get__(self):
            return self.weekday

    property n:
        def __get__(self):
            return self.n if self.n != INT64_MIN else None

    def __repr__(self):
        s = ("MO", "TU", "WE", "TH", "FR", "SA", "SU")[self.weekday]
        if self.n == INT64_MIN:
            return s
        else:
            return "%s(%+d)" % (s, self.n)

MO, TU, WE, TH, FR, SA, SU = weekdays = tuple(Weekday(x) for x in range(7))

cdef class Delta:
    """
    There's two different ways to build a Delta instance. The
    first one is passing it two Timestamp classes:

        Delta(Timestamp1, Timestamp1)

    In which case the following holds:

        Timestamp1 + Delta(Timestamp1, Timestamp2) == TimeStamp2

    And the other way is to use the following keyword arguments:

        year, month, day, hour, minute, second, microsecond:
            Absolute information.

        years, months, weeks, days, hours, minutes, seconds, microseconds:
            Relative information, may be negative.

        weekday:
            One of the weekday instances (MO, TU, etc). These instances may
            receive a parameter N, specifying the Nth weekday, which could
            be positive or negative (like MO(+1) or MO(-2). Not specifying
            it is the same as specifying +1. You can also use an integer,
            where 0=MO.

        leapdays:
            Will add given days to the date found, if year is a leap
            year, and the date found is post 28 of february.

        yearday, nlyearday:
            Set the yearday or the non-leap year day (jump leap days).
            These are converted to day/month/leapdays information.

    Here is the behavior of operations with Delta:

    1) Calculate the absolute year, using the 'year' argument, or the
    original datetime year, if the argument is not present.

    2) Add the relative 'years' argument to the absolute year.

    3) Do steps 1 and 2 for month/months.

    4) Calculate the absolute day, using the 'day' argument, or the
    original datetime day, if the argument is not present. Then,
    subtract from the day until it fits in the year and month
    found after their operations.

    5) Add the relative 'days' argument to the absolute day. Notice
    that the 'weeks' argument is multiplied by 7 and added to
    'days'.

    6) Do steps 1 and 2 for hour/hours, minute/minutes, second/seconds,
    microsecond/microseconds.

    7) If the 'weekday' argument is present, calculate the weekday,
    with the given (wday, nth) tuple. wday is the index of the
    weekday (0-6, 0=Mon), and nth is the number of weeks to add
    forward or backward, depending on its signal. Notice that if
    the calculated date is already Monday, for example, using
    (0, 1) or (0, -1) won't change the day.
    """

    cdef:
        int64_t years, months, days, leapdays, hours, minutes, seconds, microseconds
        int64_t year, month, day, hour, minute, second, microsecond
        object weekday

    def __init__(self,

                 object ts1=None,
                 object ts2=None,

                 int64_t years=0,
                 int64_t months=0,
                 int64_t days=0,
                 int64_t leapdays=0,
                 int64_t weeks=0,
                 int64_t hours=0,
                 int64_t minutes=0,
                 int64_t seconds=0,
                 int64_t microseconds=0,

                 int64_t year=-1,
                 int64_t month=-1,
                 int64_t day=-1,
                 int64_t yearday=-1,
                 int64_t nlyearday=-1,
                 int64_t hour=-1,
                 int64_t minute=-1,
                 int64_t second=-1,
                 int64_t microsecond=-1,

                 object weekday=None):

        if ts1 and ts2:
            if not (isinstance(ts1, Timestamp) and isinstance(ts2, Timestamp)):
                raise TypeError("Delta only diffs Timestamp")

            self.years = 0
            self.months = 0
            self.days = 0
            self.leapdays = 0
            self.hours = 0
            self.minutes = 0
            self.seconds = 0
            self.microseconds = 0

            self.year = -1
            self.month = -1
            self.day = -1
            self.hour = -1
            self.minute = -1
            self.second = -1
            self.microsecond = -1
            self.weekday = None

            # difference in months
            months = (ts1.year * 12 + ts1.month) - (ts2.year * 12 + ts2.month)
            self._set_months(months)

            # add ourself (delta) to ts2
            dtm = self.__add__(ts2)

            if ts1 < ts2:
                while ts1 > dtm:
                    months += 1
                    self._set_months(months)
                    dtm = self.__add__(ts2)
            else:
                while ts1 < dtm:
                    months -= 1
                    self._set_months(months)
                    dtm = self.__add__(ts2)
            delta = ts1 - dtm
            self.seconds = delta.seconds + delta.days * 86400
            self.microseconds = delta.microseconds
        else:
            self.years = years
            self.months = months
            self.days = days + weeks * 7
            self.leapdays = leapdays
            self.hours = hours
            self.minutes = minutes
            self.seconds = seconds
            self.microseconds = microseconds

            self.year = year
            self.month = month
            self.day = day
            self.hour = hour
            self.minute = minute
            self.second = second
            self.microsecond = microsecond

            if isinstance(weekday, Weekday):
                self.weekday = weekday
            elif isinstance(weekday, type(None)):
                self.weekday = None
            else:
                self.weekday = weekdays[weekday]

            yday = 0
            if nlyearday != -1:
                yday = nlyearday
            elif yearday != -1:
                yday = yearday
                if yearday > 59:
                    self.leapdays = -1
            if yday:
                ydayidx = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334,
                           366]
                for idx, ydays in enumerate(ydayidx):
                    if yday <= ydays:
                        self.month = idx + 1
                        if idx == 0:
                            self.day = yday
                        else:
                            self.day = yday - ydayidx[idx-1]
                        break
                else:
                    raise ValueError("invalid year day (%d)" % yday)

        self._fix()

    def _fix(self):
        if abs(self.microseconds) > 999999:
            s = self.microseconds//abs(self.microseconds)
            div, mod = divmod(self.microseconds*s, 1000000)
            self.microseconds = mod*s
            self.seconds += div*s
        if abs(self.seconds) > 59:
            s = self.seconds//abs(self.seconds)
            div, mod = divmod(self.seconds*s, 60)
            self.seconds = mod*s
            self.minutes += div*s
        if abs(self.minutes) > 59:
            s = self.minutes//abs(self.minutes)
            div, mod = divmod(self.minutes*s, 60)
            self.minutes = mod*s
            self.hours += div*s
        if abs(self.hours) > 23:
            s = self.hours//abs(self.hours)
            div, mod = divmod(self.hours*s, 24)
            self.hours = mod*s
            self.days += div*s
        if abs(self.months) > 11:
            s = self.months//abs(self.months)
            div, mod = divmod(self.months*s, 12)
            self.months = mod*s
            self.years += div*s

    def _set_months(self, months):
        self.months = months
        if abs(self.months) > 11:
            s = self.months//abs(self.months)
            div, mod = divmod(self.months*s, 12)
            self.months = mod*s
            self.years = div*s
        else:
            self.years = 0

    def __add__(self, other):
        if isinstance(other, Delta):
            return self._add_delta(other)

        if isinstance(other, Timestamp):
            return self._add_timestamp(other)

        raise ValueError("Cannot add to Delta")

    def _add_timestamp(self, other):
        year = (self.year if self.year != -1 else other.year) + self.years
        month = (self.month if self.month != -1 else other.month)

        if self.months:
            assert 1 <= abs(self.months) <= 12
            month += self.months
            if month > 12:
                year += 1
                month -= 12
            elif month < 1:
                year -= 1
                month += 12
        day = min(monthrange(year, month)[1],
                  self.day if self.day != -1 else other.day)
        repl = {"year": year, "month": month, "day": day}
        for attr in ["hour", "minute", "second", "microsecond"]:
            value = getattr(self, attr)
            if value != -1:
                repl[attr] = value
        days = self.days
        if self.leapdays and month > 2 and isleapyear(year):
            days += self.leapdays
        ret = (other.replace(**repl)
               + Duration(days=days,
                          hours=self.hours,
                          minutes=self.minutes,
                          seconds=self.seconds,
                          microseconds=self.microseconds))
        if self.weekday:
            weekday, nth = self.weekday.weekday, (self.weekday.n or 1)

            jumpdays = (abs(nth)-1)*7
            if nth > 0:
                jumpdays += (7-ret.weekday()+weekday)%7
            else:
                jumpdays += (ret.weekday()-weekday)%7
                jumpdays *= -1
            ret += Duration(days=jumpdays)

        return ret

    def _add_delta(self, other):
        return Delta(years=other.years+self.years,
                    months=other.months+self.months,
                    days=other.days+self.days,
                    hours=other.hours+self.hours,
                    minutes=other.minutes+self.minutes,
                    seconds=other.seconds+self.seconds,
                    microseconds=other.microseconds+self.microseconds,
                    leapdays=other.leapdays if other.leapdays != -1 else self.leapdays,
                    year=other.year if other.year != -1 else self.year,
                    month=other.month if other.month != -1 else self.month,
                    day=other.day if other.day != -1 else self.day,
                    weekday=other.weekday or self.weekday,
                    hour=other.hour if other.hour != -1 else self.hour,
                    minute=other.minute if other.minute != -1 else self.minute,
                    second=other.second if other.second != -1 else self.second,
                    microsecond=(other.microsecond if other.microsecond != -1
                                                    else self.microsecond))


    def __sub__(self, other):
        if isinstance(other, Delta):
            return self._sub_delta(other)
        else:
            return self.__neg__().__add__(other)

    def _sub_delta(self, other):
        return Delta(years=other.years-self.years,
                    months=other.months-self.months,
                    days=other.days-self.days,
                    hours=other.hours-self.hours,
                    minutes=other.minutes-self.minutes,
                    seconds=other.seconds-self.seconds,
                    microseconds=other.microseconds-self.microseconds,
                    leapdays=other.leapdays if other.leapdays != -1 else self.leapdays,
                    year=other.year if other.year != -1 else self.year,
                    month=other.month if other.month != -1 else self.month,
                    day=other.day if other.day != -1 else self.day,
                    weekday=other.weekday or self.weekday,
                    hour=other.hour if other.hour != -1 else self.hour,
                    minute=other.minute if other.minute != -1 else self.minute,
                    second=other.second if other.second != -1 else self.second,
                    microsecond=(other.microsecond if other.microsecond != -1
                                                    else self.microsecond))

    def __neg__(self):
        return Delta(years=-self.years,
                     months=-self.months,
                     days=-self.days,
                     hours=-self.hours,
                     minutes=-self.minutes,
                     seconds=-self.seconds,
                     microseconds=-self.microseconds,
                     leapdays=self.leapdays,
                     year=self.year,
                     month=self.month,
                     day=self.day,
                     weekday=self.weekday,
                     hour=self.hour,
                     minute=self.minute,
                     second=self.second,
                     microsecond=self.microsecond)


    def __mul__(self, int f):
        return Delta(years=self.years*f,
                     months=self.months*f,
                     days=self.days*f,
                     hours=self.hours*f,
                     minutes=self.minutes*f,
                     seconds=self.seconds*f,
                     microseconds=self.microseconds*f,
                     leapdays=self.leapdays,
                     year=self.year,
                     month=self.month,
                     day=self.day,
                     weekday=self.weekday,
                     hour=self.hour,
                     minute=self.minute,
                     second=self.second,
                     microsecond=self.microsecond)

    def __repr__(self):
        l = []
        for attr in ["years", "months", "days", "leapdays",
                     "hours", "minutes", "seconds", "microseconds"]:
            value = getattr(self, attr)
            if value:
                l.append("%s=%+d" % (attr, value))
        for attr in ["year", "month", "day", "weekday",
                     "hour", "minute", "second", "microsecond"]:
            value = getattr(self, attr)
            if value != -1:
                l.append("%s=%s" % (attr, value))
        return "%s(%s)" % (self.__class__.__name__, ", ".join(l))

    property year:
        def __get__(self):
            return self.year

    property month:
        def __get__(self):
            return self.month

    property day:
        def __get__(self):
            return self.day

    property weekday:
        def __get__(self):
            return self.weekday

    property hour:
        def __get__(self):
            return self.hour

    property minute:
        def __get__(self):
            return self.minute

    property second:
        def __get__(self):
            return self.second

    property microsecond:
        def __get__(self):
            return self.microsecond

    property years:
        def __get__(self):
            return self.years

    property months:
        def __get__(self):
            return self.months

    property days:
        def __get__(self):
            return self.days

    property leapdays:
        def __get__(self):
            return self.leapdays

    property hours:
        def __get__(self):
            return self.hours

    property minutes:
        def __get__(self):
            return self.minutes

    property seconds:
        def __get__(self):
            return self.seconds

    property microseconds:
        def __get__(self):
            return self.microseconds

# End derivation from dateutil

# Conversion routines
# ------------------------------------------------------------------------------

def pydt_to_i8(object pydt):
    '''
    Convert from python datetime object to int64 representation compatible with
    numpy datetime64; converts to UTC
    '''
    cdef:
        npy_datetimestruct dts
        NPY_DATETIMEUNIT out_bestunit

    if PyDateTime_Check(pydt):
        # TODO: this function can prob be optimized
        convert_pydatetime_to_datetimestruct(<PyObject *>pydt, &dts,
                                             &out_bestunit, 1)

        return PyArray_DatetimeStructToDatetime(out_bestunit, &dts)

    raise ValueError("Expected a datetime, received a %s" % type(pydt))

def i8_to_pydt(int64_t i8, object tzinfo = None):
    '''
    Inverse of pydt_to_i8
    '''
    cdef:
        npy_datetimestruct dts
        object result

    PyArray_DatetimeToDatetimeStruct(i8, NPY_FR_us, &dts)

    result = <object>PyDateTime_FromDateAndTime(dts.year, dts.month, dts.day,
                                                dts.hour, dts.min, dts.sec, dts.us)

    return result


# Accessors
# ------------------------------------------------------------------------------

def fast_field_accessor(ndarray[int64_t] dtindex, object field):
    '''
    Given a int64-based datetime index, extract the year, month, etc.,
    field and return an array of these values.
    '''
    cdef:
        npy_datetimestruct dts
        Py_ssize_t i, count = 0
        ndarray[int32_t] out

    count = len(dtindex)
    out = np.empty(count, dtype='i4')

    if field == 'Y':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.year
        return out

    elif field == 'M':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.month
        return out

    elif field == 'D':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.day
        return out

    elif field == 'h':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.hour
        return out

    elif field == 'm':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.min
        return out

    elif field == 's':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.sec
        return out

    elif field == 'us':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.us
        return out

    else:
        raise ValueError("Field %s not supported; not in (Y,M,D,h,m,s,us)" % field)

# Another accessor, for datetime object -
# ------------------------------------------------------------------------------

def fast_field_accessor2(ndarray[object] dtindex, object field):
    '''
    Given a int64-based datetime index, extract the year, month, etc.,
    field and return an array of these values.
    '''
    cdef:
        npy_datetimestruct dts
        Py_ssize_t i, count = 0
        ndarray[int32_t] out

    count = len(dtindex)
    out = np.empty(count, dtype='i4')

    if field == 'Y':
        for i in range(count):
                out[i] = PyDateTime_GET_YEAR(dtindex[i])
        return out
    else:
        raise ValueError("Field %s not supported; not in (Y,M,D,h,m,s,us)" % field)

# Some general helper functions
# ------------------------------------------------------------------------------

def isleapyear(int64_t year):
    return is_leapyear(year)

def monthrange(int64_t year, int64_t month):
    cdef:
        int64_t days
        int64_t day_of_week

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = _days_per_month_table[is_leapyear(year)][month-1]

    return (dayofweek(year, month, 1), days)

