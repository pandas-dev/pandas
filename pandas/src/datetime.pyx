cimport numpy as np
import numpy as np

from numpy cimport int32_t, int64_t, import_array, ndarray
from cpython cimport *
from libc.stdlib cimport malloc, free

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

# Objects to support date/time arithmetic, inspired by the architecture of the
# lubridate R package, to eventually handle all datetime logic in pandas.
# --------------------------------------------------------------------------------

cdef class Timestamp:
    '''
    A timestamp (absolute moment in time) to microsecond resolution, in UTC.
    '''
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

    def __sub__(self, object other):
        ''' 
        Subtract two timestamps, results in an interval with the start being
        the earlier of the two timestamps. 
        '''
        if isinstance(other, Timestamp):
            if other.value > self.value:
                return Interval(self.value, other.value)
            else:
                return Interval(other.value, self.value)

    def __add__(self, object other):
        '''
        Add an Interval, Duration, or Period to the Timestamp, resulting in 
        new Timestamp.
        '''
        if isinstance(other, Interval):
            return Timestamp(self.value + other.dur.length)
        elif isinstance(other, Duration):
            return Timestamp(self.value + other.length)
        elif isinstance(other, Period):
            # TODO: fix me
            raise ValueError("TODO: Period needs to be implemented")

    def __str__(self):
        '''
        Output ISO8601 format string representation of timestamp.
        '''
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

    property value:
        def __get__(self):
            return self.value

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
            return self.dts.us // 1000.

    property us:
        def __get__(self):
            return self.dts.us

cdef class Interval:
    '''
    An absolute time span, from one timestamp to another
    '''
    cdef:
        Timestamp start
        Timestamp end

    def __init__(self, Timestamp start, Timestamp end):
        self.start = start
        self.end = end

    property start:
        def __get__(self):
            return self.start

    property end:
        def __get__(self):
            return self.end

    property us:
        def __get__(self):
            return Duration(self.end - self.start).us

    property ms:
        def __get__(self):
            return Duration(self.end - self.start).ms

    property secs:
        def __get__(self):
            return Duration(self.end - self.start).secs

cdef class Duration:
    '''
    Absolute length of time, in microseconds
    '''
    cdef int64_t length

    def __init__(self, int64_t us = 1):
        self.length = us

    def __str__(self):
        return "Duration (%d)" % self.length

    property us:
        def __get__(self):
            return self.length

    property ms:
        def __get__(self):
            return self.length // 1000

    property secs:
        def __get__(self):
            return self.length // 1000000

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
        raise ValueError("Could not add Period to operand")

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

# Some general helper functions
# ------------------------------------------------------------------------------

def isleapyear(int64_t year):
    return is_leapyear(year)
