cimport numpy as np
import numpy as np

from numpy cimport int32_t, int64_t, import_array, ndarray
from cpython cimport *

# this is our datetime.pxd
from datetime cimport *

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()

# in numpy 1.7, will prob need this
# numpy_pydatetime_import

# Objects to support date/time arithmetic, inspired by the architecture of the
# lubridate R package, to eventually replace most of pandas/core/datetools.py
# --------------------------------------------------------------------------------

cdef class Instant:
    '''
    A timestamp (absolute moment in time) to microsecond resolution, in UTC.
    '''
    cdef:
        int64_t timestamp
        npy_datetimestruct dts

    def __init__(self, int64_t ts):
        self.timestamp = ts

        # decompose datetime64 to components
        PyArray_DatetimeToDatetimeStruct(ts, NPY_FR_us, &self.dts)

    def __sub__(self, object other):
        if isinstance(other, Instant):
            if other.timestamp > self.timestamp:
                return Interval(self.timestamp, other.timestamp)
            else:
                return Interval(other.timestamp, self.timestamp)

    def __add__(self, object other):
        if isinstance(other, Interval):
            return Instant(self.timestamp + other.dur.length)
        elif isinstance(other, Duration):
            return Instant(self.timestamp + other.length)

    property timestamp:
        def __get__(self):
            return self.timestamp

    property year:
        def __get__(self):
            return self.dts.year

    property month:
        def __get__(self):
            return self.dts.month

    property day:
        def __get__(self):
            return self.dts.day

    property hour:
        def __get__(self):
            return self.dts.hour

    property min:
        def __get__(self):
            return self.dts.min

    property sec:
        def __get__(self):
            return self.dts.sec

    property ms:
        def __get__(self):
            return self.dts.us / 1000.

    property us:
        def __get__(self):
            return self.dts.us

cdef class Interval:
    '''
    An absolute time span, from one instant to another
    '''
    cdef:
        int64_t start
        int64_t end

    def __init__(self, int64_t start, int64_t end):
        self.start = start
        self.end = end

    property start:
        def __get__(self):
            return self.start

    property end:
        def __get__(self):
            return self.end

    property dur:
        def __get__(self):
            return Duration(self.end - self.start)

cdef class Duration:
    '''
    Absolute length of time, in microseconds
    '''
    cdef int64_t length

    def __init__(self, int64_t length = 1):
        self.length = length

    def __str__(self):
        return "Duration (%d)" % self.length

    property us:
        def __get__(self):
            return self.length

    property ms:
        def __get__(self):
            return self.length / 1000.

    property secs:
        def __get__(self):
            return self.length / 1000000.

cdef class Period:
    '''
    Relative length of time
    '''
    cdef:
        npy_datetimestruct dts
        int isbiz

    def __init__(self, int years = 0, int months = 0, int days = 0,
                       int hours = 0, int mins = 0, int secs = 0):
        self.dts.year = years
        self.dts.month = months
        self.dts.day = days
        self.dts.hour = hours
        self.dts.min = mins
        self.dts.sec = secs

    def __add__(self, other):
        if issubclass(other, Period):
            return Period(years = self.years + other.years,
                          months = self.months + other.months,
                          days = self.days + other.days,
                          hours = self.hours + other.hours,
                          mins = self.mins + other.mins,
                          secs = self.secs + other.secs)
        raise ValueError("Could not add operand to Period")

    def __str__(self):
        strbuf = ""
        if self.dts.year > 0:
            strbuf += "%d " % self.dts.year
            strbuf += "year"
            if self.dts.year > 1:
                strbuf += "s"
        if self.dts.month > 0:
            if len(strbuf):
                strbuf += ", "
            strbuf += "%d " % self.dts.month
            strbuf += "month"
            if self.dts.month > 1:
                strbuf += "s"
        if self.dts.day > 0:
            if len(strbuf):
                strbuf += ", "
            strbuf += "%d " % self.dts.day
            strbuf += "day"
            if self.dts.day > 1:
                strbuf += "s"
        if self.dts.hour > 0:
            if len(strbuf):
                strbuf += ", "
            strbuf += "%d " % self.dts.hour
            strbuf += "hour"
            if self.dts.hour > 1:
                strbuf += "s"
        if self.dts.min > 0:
            if len(strbuf):
                strbuf += ", "
            strbuf += "%d " % self.dts.min
            strbuf += "min"
            if self.dts.min > 1:
                strbuf += "s"
        if self.dts.sec > 0:
            if len(strbuf):
                strbuf += ", "
            strbuf += "%d " % self.dts.sec
            strbuf += "sec"
            if self.dts.sec > 1:
                strbuf += "s"

        return "Period: %s" % strbuf

    property years:
        def __get__(self):
            return self.dts.year
    property months:
        def __get__(self):
            return self.dts.month
    property days:
        def __get__(self):
            return self.dts.day
    property hours:
        def __get__(self):
            return self.dts.hour
    property mins:
        def __get__(self):
            return self.dts.min
    property secs:
        def __get__(self):
            return self.dts.sec

def seconds(int count):
    return Period(secs = count)

def minutes(int count):
    return Period(mins = count)

def hours(int count):
    return Period(hours = count)

def days(int count):
    return Period(days = count)

# TODO: fixme
def bdays(int count):
    return Period(days = count)

def weeks(int count):
    return Period(days = 7 * count)

def months(int count):
    return Period(months = count)

def quarters(int count):
    return Period(months = 3 * count)

def years(int count):
    return Period(years = count)

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

# Some general helper functions we need
# ------------------------------------------------------------------------------

def isleapyear(int64_t year):
    return is_leapyear(year)
