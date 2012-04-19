# cython: profile=False

cimport numpy as np
import numpy as np

from numpy cimport int32_t, int64_t, import_array, ndarray
from cpython cimport *

from libc.stdlib cimport malloc, free, abs
from libc.math cimport floor

# this is our datetime.pxd
from datetime cimport *
from util cimport is_integer_object, is_datetime64_object

from dateutil.parser import parse as parse_date

# initialize numpy
import_array()
import_ufunc()

# import datetime C API
PyDateTime_IMPORT

# in numpy 1.7, will prob need the following:
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

try:
    basestring
except NameError: # py3
    basestring = str


# Python front end to C extension type _Timestamp
# This serves as the box for datetime64
class Timestamp(_Timestamp):

    def __new__(cls, object ts_input, object offset=None, tz=None):
        if isinstance(ts_input, float):
            # to do, do we want to support this, ie with fractional seconds?
            raise TypeError("Cannot convert a float to datetime")

        if isinstance(ts_input, basestring):
            try:
                ts_input = parse_date(ts_input)
            except Exception:
                pass

        ts = convert_to_tsobject(ts_input, tz)

        # make datetime happy
        ts_base = _Timestamp.__new__(
            cls,
            ts.dtval.year,
            ts.dtval.month,
            ts.dtval.day,
            ts.dtval.hour,
            ts.dtval.minute,
            ts.dtval.second,
            ts.dtval.microsecond,
            ts.dtval.tzinfo)

        # fill out rest of data
        ts_base.value = ts.value
        ts_base.offset = offset

        return ts_base

    @property
    def tz(self):
        """
        Alias for tzinfo
        """
        return self.tzinfo

    @property
    def freq(self):
        return self.offset

    def __setstate__(self, state):
        self.value = state[0]
        self.offset = state[1]
        self.tzinfo = state[2]

    def __reduce__(self):
        object_state = self.value, self.offset, self.tzinfo
        return (Timestamp, object_state)

    def to_period(self, freq=None):
        """
        Return an period of which this timestamp is an observation.
        """
        from pandas.tseries.period import Period

        if freq == None:
            freq = self.freq

        return Period(self, freq=freq)

def apply_offset(ndarray[object] values, object offset):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] new_values
        object boxed

    result = np.empty(n, dtype='M8[us]')
    new_values = result.view('i8')
    pass

# This is PITA. Because we inherit from datetime, which has very specific
# construction requirements, we need to do object instantiation in python
# (see Timestamp class above). This will serve as a C extension type that
# shadows the python class, where we do any heavy lifting.
cdef class _Timestamp(datetime):
    cdef:
        int64_t value       # numpy int64
        object offset       # frequency reference

    def __add__(self, other):
        if is_integer_object(other):
            if self.offset is None:
                msg = ("Cannot add integral value to Timestamp "
                       "without offset.")
                raise ValueError(msg)
            else:
                return Timestamp((self.offset.__mul__(other)).apply(self))
        else:
            return datetime.__add__(self, other)

    def __sub__(self, other):
        if is_integer_object(other):
            return self.__add__(-other)
        else:
            return datetime.__sub__(self, other)

# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    cdef:
        datetime dtval      # python datetime
        int64_t value       # numpy dt64

    property dtval:
        def __get__(self):
            return self.dtval

    property value:
        def __get__(self):
            return self.value

# helper to extract datetime and int64 from several different possibilities
cdef convert_to_tsobject(object ts, object tzinfo=None):
    """
    Extract datetime and int64 from any of:
        - np.int64
        - np.datetime64
        - python int or long object
        - iso8601 string object
        - python datetime object
        - another timestamp object
    """
    cdef:
        Py_ssize_t strlen
        npy_bool islocal, special
        NPY_DATETIMEUNIT out_bestunit
        npy_datetimestruct dts
        _Timestamp tmp
        _TSObject retval

    if isinstance(ts, _TSObject) or ts is None:
        return ts

    retval = _TSObject()

    # pretty expensive - faster way to access as i8?
    if is_datetime64_object(ts):
        retval.value = ts.view('i8')
        PyArray_DatetimeToDatetimeStruct(retval.value, NPY_FR_us, &dts)
        retval.dtval = <object>PyDateTime_FromDateAndTime(
                            dts.year, dts.month,
                            dts.day, dts.hour,
                            dts.min, dts.sec, dts.us)
    # this is cheap
    elif is_integer_object(ts):
        retval.value = ts
        PyArray_DatetimeToDatetimeStruct(retval.value, NPY_FR_us, &dts)
        retval.dtval = <object>PyDateTime_FromDateAndTime(
                            dts.year, dts.month,
                            dts.day, dts.hour,
                            dts.min, dts.sec, dts.us)
    # this is pretty cheap
    elif PyString_Check(ts):
        # we might want to fall back on dateutil parser?
        parse_iso_8601_datetime(ts, len(ts), NPY_FR_us, NPY_UNSAFE_CASTING,
                                &dts, &islocal, &out_bestunit, &special)
        retval.value = PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts)
        retval.dtval = <object>PyDateTime_FromDateAndTime(
                            dts.year, dts.month,
                            dts.day, dts.hour,
                            dts.min, dts.sec, dts.us)
    # pretty cheap
    elif PyDateTime_Check(ts):
        retval.dtval = ts
        # to do this is expensive (10x other constructors)
        # convert_pydatetime_to_datetimestruct(<PyObject *>ts, &dts,
        #                                     &out_bestunit, 0)
        dts.year = PyDateTime_GET_YEAR(ts)
        dts.month = PyDateTime_GET_MONTH(ts)
        dts.day = PyDateTime_GET_DAY(ts)
        dts.hour = PyDateTime_DATE_GET_HOUR(ts)
        dts.min = PyDateTime_DATE_GET_MINUTE(ts)
        dts.sec = PyDateTime_DATE_GET_SECOND(ts)
        dts.us = PyDateTime_DATE_GET_MICROSECOND(ts)
        retval.value = PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts)
    elif PyDate_Check(ts):
        dts.year = PyDateTime_GET_YEAR(ts)
        dts.month = PyDateTime_GET_MONTH(ts)
        dts.day = PyDateTime_GET_DAY(ts)
        retval.dtval = PyDateTime_FromDateAndTime(dts.year, dts.month, dts.day,
                                                  0, 0, 0, 0)
        dts.hour = 0
        dts.min = 0
        dts.sec = 0
        dts.us = 0
        retval.value = PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts)
    # pretty cheap
    elif isinstance(ts, _Timestamp):
        tmp = ts
        retval.value = tmp.value
        retval.dtval = tmp
    # fallback, does it at least have the right fields?
    elif isinstance(ts, object):
        dts.year = ts.year
        dts.month = ts.month
        dts.day = ts.day
        dts.hour = ts.hour
        dts.min = ts.minute
        dts.sec = ts.second
        dts.us = ts.microsecond
        retval.dtval = <object>PyDateTime_FromDateAndTime(
                            dts.year, dts.month,
                            dts.day, dts.hour,
                            dts.min, dts.sec, dts.us)
        retval.value = PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts)
    else:
        raise ValueError("Could not construct Timestamp from argument %s" % type(ts))

    if tzinfo is not None:
        retval.dtval = retval.dtval.replace(tzinfo=tzinfo)

    return retval

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

# Logic to generate ranges
# -----------------------------------------------------------------------------

cdef inline int64_t weekend_adjustment(int64_t dow, int bkwd):
    if dow > 4:                         # sat or sun?
        if bkwd:                        # roll back 1 or 2 days
            return (4 - dow)
        else:                           # roll forward 2 or 1 days
            return (7 - dow)
    return 0

cdef int64_t us_in_day = conversion_factor(r_microsecond, r_day)

cdef class _Offset:
    """
    Base class to generate timestamps. Set the anchor, and then move offsets
    with next & prev. Retrieve timestamp with ts attribute.
    """
    cdef:
        int64_t t, dow, biz, dayoffset
        object start
        _TSObject ts

    def __cinit__(self):
        self.t=0
        self.dow=0
        self.biz=0
        self.dayoffset=0

    cpdef anchor(self, object start=None):
        if start is not None:
            self.start = start
        self.ts = convert_to_tsobject(self.start)
        self._setup()

    cdef _setup(self):
        pass

    cpdef next(self):
        pass

    cpdef prev(self):
        pass

    cdef int64_t _ts(self):
        """
        Access the current timestamp value, with a possible weekday
        adjustment.
        """
        cdef int64_t adj

        if self.biz != 0:
            adj = weekend_adjustment(self.dow, self.biz < 0)
            return self.t + us_in_day * adj
        else:
            return self.t

    cdef int64_t _get_anchor(self):
        """
        Retrieve an anchor relating to current offset we're on.
        """
        return self.t - self.dayoffset * us_in_day

    property ts:
        def __get__(self):
            return self._ts()

cdef class YearOffset(_Offset):
    """
    Generate annual timestamps from provided start time; apply dayoffset to
    each timestamp. If biz > 0, we choose the next business day at each time;
    previous if < 0.

    Parameters
    ----------
    dayoffset : int
    biz : int
    """
    cdef:
        int64_t y, ly

    def __init__(self, int64_t dayoffset=0, int64_t biz=0, object anchor=None):
        self.dayoffset = dayoffset
        self.biz = biz

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts

        self.t = ts.value + self.dayoffset * us_in_day
        self.y = ts.dtval.year

        self.ly = (ts.dtval.month > 2 or
                   ts.dtval.month == 2 and ts.dtval.day == 29)

        if self.biz != 0:
            self.dow = (ts_dayofweek(ts) + self.dayoffset) % 7

    cpdef next(self):
        cdef int64_t days

        days = 365 + is_leapyear(self.y + self.ly)

        self.t += days * us_in_day
        self.y += 1

        if self.biz != 0:
            self.dow = (self.dow + days) % 7

    cpdef prev(self):
        cdef int64_t days

        days = 365 + is_leapyear(self.y - (1-self.ly))

        self.t -= days * us_in_day
        self.y -= 1

        if self.biz != 0:
            self.dow = (self.dow - days) % 7

cdef class MonthOffset(_Offset):
    """
    Generate monthly timestamps from provided start time, and apply dayoffset
    to each timestamp.  Stride to construct strided timestamps (eg quarterly).
    If biz > 0, we choose the next business day at each time; previous if < 0.

    Parameters
    ----------
    dayoffset : int
    stride : int, > 0
    biz : int
    """
    cdef:
        Py_ssize_t stride, ly, m
        int64_t y

    def __init__(self, int64_t dayoffset=0, Py_ssize_t stride=1,
                 int64_t biz=0, object anchor=None):
        self.dayoffset = dayoffset
        self.stride = stride
        self.biz = biz

        if stride <= 0:
            raise ValueError("Stride must be positive")

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts

        self.t = ts.value + (self.dayoffset * us_in_day)

        # for day counting
        self.m  = ts.dtval.month - 1
        self.y  = ts.dtval.year
        self.ly = is_leapyear(self.y)

        if self.biz != 0:
            self.dow = (ts_dayofweek(ts) + self.dayoffset) % 7

    cpdef next(self):
        cdef:
            int64_t tmp, days
            Py_ssize_t j

        days = 0
        for j in range(0, self.stride):
            if self.m >= 12:
                self.m -= 12
                self.y += 1
                self.ly = is_leapyear(self.y)
            days += _days_per_month_table[self.ly][self.m]
            self.m += 1

        self.t += days * us_in_day

        if self.biz != 0:
            self.dow = (self.dow + days) % 7

    cpdef prev(self):
        cdef:
            int64_t tmp, days
            Py_ssize_t j

        days = 0
        for j in range(0, self.stride):
            self.m -= 1
            if self.m < 0:
                self.m += 12
                self.y -= 1
                self.ly = is_leapyear(self.y)
            days += _days_per_month_table[self.ly][self.m]

        self.t -= days * us_in_day

        if self.biz != 0:
            self.dow = (self.dow - days) % 7

cdef class DayOfMonthOffset(_Offset):
    """
    Generate relative monthly timestamps from month & year of provided start
    time. For example, fridays of the third week of each month (week=3, day=4);
    or, thursdays of the last week of each month (week=-1, day=3).

    Parameters
    ----------
    week : int
    day : int, 0 to 6
    """
    cdef:
        Py_ssize_t ly, m
        int64_t y, day, week

    def __init__(self, int64_t week=0, int64_t day=0, object anchor=None):
        self.week = week
        self.day = day

        if self.day < 0 or self.day > 6:
            raise ValueError("Day offset must be 0 to 6")

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts

        # rewind to beginning of month
        self.t = ts.value - (ts.dtval.day - 1) * us_in_day
        self.dow = dayofweek(ts.dtval.year, ts.dtval.month, 1)

        # for day counting
        self.m = ts.dtval.month - 1
        self.y = ts.dtval.year
        self.ly = is_leapyear(self.y)

    cpdef next(self):
        cdef:
            int64_t tmp, days

        days = _days_per_month_table[self.ly][self.m]
        self.t += days * us_in_day
        self.dow = (self.dow + days) % 7

        self.m += 1
        if self.m >= 12:
            self.m -= 12
            self.y += 1
            self.ly = is_leapyear(self.y)

    cpdef prev(self):
        cdef:
            int64_t tmp, days

        days = _days_per_month_table[self.ly][(self.m - 1) % 12]
        self.t -= days * us_in_day
        self.dow = (self.dow - days) % 7

        self.m -= 1
        if self.m < 0:
            self.m += 12
            self.y -= 1
            self.ly = is_leapyear(self.y)

    cdef int64_t _ts(self):
        """
        Overwrite default adjustment
        """
        cdef int64_t adj = (self.week * 7) + (self.day - self.dow) % 7
        return self.t + us_in_day * adj

cdef class DayOffset(_Offset):
    """
    Generate daily timestamps beginning with first valid time >= start time. If
    biz != 0, we skip weekends. Stride, to construct weekly timestamps.

    Parameters
    ----------
    stride : int, > 0
    biz : boolean
    """
    cdef:
        Py_ssize_t stride

    def __init__(self, int64_t stride=1, int64_t biz=0, object anchor=None):
        self.stride = stride
        self.biz = biz

        if self.stride <= 0:
            raise ValueError("Stride must be positive")

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts
        self.t = ts.value
        if self.biz != 0:
            self.dow = ts_dayofweek(ts)

    cpdef next(self):
        self.t += (self.stride * us_in_day)
        if self.biz != 0:
            self.dow = (self.dow + self.stride) % 7
            if self.dow >= 5:
                self.t += (7 - self.dow) * us_in_day
                self.dow = 0

    cpdef prev(self):
        self.t -= (self.stride * us_in_day)
        if self.biz != 0:
            self.dow = (self.dow - self.stride) % 7
            if self.dow >= 5:
                self.t += (4 - self.dow) * us_in_day
                self.dow = 4

#cdef ndarray[int64_t] _generate_range(_Offset offset, Py_ssize_t periods):
#    """
#    Generate timestamps according to offset.
#    """
#    cdef:
#        Py_ssize_t i
#        ndarray[int64_t] dtindex

#    dtindex = np.empty(periods, np.int64)
#    for i in range(periods):
#        dtindex[i] = offset._ts()
#        offset.next()
#    return dtindex

#cdef int64_t _count_range(_Offset offset, object end):
#    """
#    Count timestamps in range according to offset up to (and including)
#    end time.
#    """
#    cdef:
#        Py_ssize_t i=0
#        _TSObject e

#    e = convert_to_tsobject(end)
#    while offset._ts() <= e.value:
#        i += 1
#        offset.next()
#    return i

def string_to_datetime(ndarray[object] strings, raise_=False, dayfirst=False):
    cdef:
        Py_ssize_t i, n = len(strings)
        object val
        ndarray[int64_t] iresult
        ndarray[object] oresult

    from dateutil.parser import parse


    try:
        result = np.empty(n, dtype='M8[us]')
        iresult = result.view('i8')
        for i in range(n):
            val = strings[i]
            if util._checknull(val):
                result[i] = NaT
            elif PyDateTime_Check(val):
                result[i] = val
            else:
                if len(val) == 0:
                    result[i] = NaT
                    continue
                try:
                    result[i] = parse(val, dayfirst=dayfirst)
                except Exception:
                    raise TypeError
        return result
    except TypeError:
        oresult = np.empty(n, dtype=object)

        for i in range(n):
            val = strings[i]
            if util._checknull(val):
                oresult[i] = val
            else:
                if len(val) == 0:
                    oresult[i] = NaT
                    continue
                try:
                    oresult[i] = parse(val, dayfirst=dayfirst)
                except Exception:
                    if raise_:
                        raise
                    oresult[i] = val

        return oresult


# Conversion routines
# ------------------------------------------------------------------------------

def pydt_to_i8(object pydt):
    '''
    Convert to int64 representation compatible with numpy datetime64; converts
    to UTC
    '''
    cdef:
        _TSObject ts

    ts = convert_to_tsobject(pydt)

    return ts.value

def i8_to_pydt(int64_t i8, object tzinfo = None):
    '''
    Inverse of pydt_to_i8
    '''
    return Timestamp(i8)

# time zone conversion helpers
# ------------------------------------------------------------------------------

try:
    import pytz
    have_pytz = True
except:
    have_pytz = False

trans_cache = {}
utc_offset_cache = {}

cdef ndarray[int64_t] _get_transitions(object tz):
    """
    Get UTC times of DST transitions
    """
    if tz not in trans_cache:
        arr = np.array(tz._utc_transition_times, dtype='M8[us]')
        trans_cache[tz] = np.array(arr.view('i8'))
    return trans_cache[tz]

cdef ndarray[int64_t] _unbox_utcoffsets(object transinfo):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] arr

    sz = len(transinfo)
    arr = np.empty(sz, dtype='i8')

    for i in range(sz):
        arr[i] = int(transinfo[i][0].total_seconds()) * 1000000

    return arr

cdef int64_t get_utcoffset(object tz, Py_ssize_t idx):
    """
    Get UTC offsets in microseconds corresponding to DST transitions
    """
    cdef:
        ndarray[int64_t] arr
    if tz not in utc_offset_cache:
        utc_offset_cache[tz] = _unbox_utcoffsets(tz._transition_info)
    arr = utc_offset_cache[tz]
    return arr[idx]

def tz_normalize_array(ndarray[int64_t] vals, object tz1, object tz2):
    """
    Convert DateRange from one time zone to another (using pytz)

    Returns
    -------
    normalized : DateRange
    """
    cdef:
        ndarray[int64_t] result
        ndarray[int64_t] trans
        Py_ssize_t i, sz, tzidx
        int64_t v, tz1offset, tz2offset

    if not have_pytz:
        raise Exception("Could not find pytz module")

    sz = len(vals)

    if sz == 0:
        return np.empty(0, dtype=np.int64)

    result = np.empty(sz, dtype=np.int64)
    trans = _get_transitions(tz1)

    tzidx = np.searchsorted(trans, vals[0])

    tz1offset = get_utcoffset(tz1, tzidx)
    tz2offset = get_utcoffset(tz2, tzidx)

    for i in range(sz):
        v = vals[i]
        if v >= trans[tzidx + 1]:
            tzidx += 1
            tz1offset = get_utcoffset(tz1, tzidx)
            tz2offset = get_utcoffset(tz2, tzidx)

        result[i] = (v - tz1offset) + tz2offset

    return result

def tz_localize_array(ndarray[int64_t] vals, object tz):
    """
    Localize tzinfo-naive DateRange to given time zone (using pytz). If
    there are ambiguities in the values, raise AmbiguousTimeError.

    Returns
    -------
    localized : DatetimeIndex
    """
    cdef:
        ndarray[int64_t] trans
        Py_ssize_t i, sz, tzidx
        int64_t v, t1, t2, currtrans, tmp

    if not have_pytz:
        raise Exception("Could not find pytz module")

    if tz == pytz.utc or tz is None:
        return vals

    sz = len(vals)

    if sz == 0:
        return np.empty(0, dtype=np.int64)

    result = np.empty(sz, dtype=np.int64)
    trans = _get_transitions(tz)
    tzidx = np.searchsorted(trans, vals[0])

    currtrans = trans[tzidx]
    t1 = currtrans + get_utcoffset(tz, tzidx-1)
    t2 = currtrans + get_utcoffset(tz, tzidx)

    for i in range(sz):
        v = vals[i]
        if v >= trans[tzidx + 1]:
            tzidx += 1
            currtrans = trans[tzidx]
            t1 = currtrans + get_utcoffset(tz, tzidx-1)
            t2 = currtrans + get_utcoffset(tz, tzidx)

        if t1 > t2:
            tmp = t1
            t1 = t2
            t2 = tmp

        if t1 <= v and v <= t2:
            msg = "Cannot localize, ambiguous time %s found" % Timestamp(v)
            raise pytz.AmbiguousTimeError(msg)

    return vals

# Accessors
# ------------------------------------------------------------------------------

def build_field_sarray(ndarray[int64_t] dtindex):
    '''
    Datetime as int64 representation to a structured array of fields
    '''
    cdef:
        Py_ssize_t i, count = 0
        int isleap
        npy_datetimestruct dts

    count = len(dtindex)

    sa_dtype = [('Y', '>i4'), # year
                ('M', '>i4'), # month
                ('D', '>i4'), # day
                ('h', '>i4'), # hour
                ('m', '>i4'), # min
                ('s', '>i4'), # second
                ('u', '>i4')] # microsecond

    out = np.empty(count, dtype=sa_dtype)

    for i in range(count):
        PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
        out[i] = (dts.year, dts.month, dts.day, dts.hour, dts.min, dts.sec,
                  dts.us)
    return out

@cython.wraparound(False)
def fast_field_accessor(ndarray[int64_t] dtindex, object field):
    '''
    Given a int64-based datetime index, extract the year, month, etc.,
    field and return an array of these values.
    '''
    cdef:
        _TSObject ts
        Py_ssize_t i, count = 0
        ndarray[int32_t] out
        ndarray[int32_t, ndim=2] _month_offset
        int isleap
        npy_datetimestruct dts

    _month_offset = np.array(
        [[ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 ],
         [ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 ]],
         dtype=np.int32 )

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

    elif field == 'doy':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            isleap = is_leapyear(dts.year)
            out[i] = _month_offset[isleap, dts.month-1] + dts.day
        return out

    elif field == 'dow':
        for i in range(count):
            ts = convert_to_tsobject(dtindex[i])
            out[i] = ts_dayofweek(ts)
        return out

    elif field == 'woy':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            isleap = is_leapyear(dts.year)
            out[i] = _month_offset[isleap, dts.month - 1] + dts.day
            out[i] = ((out[i] - 1) / 7) + 1
        return out

    elif field == 'q':
        for i in range(count):
            PyArray_DatetimeToDatetimeStruct(dtindex[i], NPY_FR_us, &dts)
            out[i] = dts.month
            out[i] = ((out[i] - 1) / 3) + 1
        return out

    raise ValueError("Field %s not supported" % field)

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

cdef inline int64_t ts_dayofweek(_TSObject ts):
    return dayofweek(ts.dtval.year, ts.dtval.month, ts.dtval.day)

# Period logic
# ------------------------------------------------------------------------------

cdef long apply_mult(long period_ord, long mult):
    """
    Get base+multiple ordinal value from corresponding base-only ordinal value.
    For example, 5min ordinal will be 1/5th the 1min ordinal (rounding down to
    integer).
    """
    if mult == 1:
        return period_ord

    return (period_ord - 1) // mult

cdef long remove_mult(long period_ord_w_mult, long mult):
    """
    Get base-only ordinal value from corresponding base+multiple ordinal.
    """
    if mult == 1:
        return period_ord_w_mult

    return period_ord_w_mult * mult + 1;

def dt64arr_to_periodarr(ndarray[int64_t] dtarr, int base, long mult):
    """
    Convert array of datetime64 values (passed in as 'i8' dtype) to a set of
    periods corresponding to desired frequency, per period convention.
    """
    cdef:
        ndarray[int64_t] out
        Py_ssize_t i, l
        npy_datetimestruct dts

    l = len(dtarr)

    out = np.empty(l, dtype='i8')

    for i in range(l):
        PyArray_DatetimeToDatetimeStruct(dtarr[i], NPY_FR_us, &dts)
        out[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                  dts.hour, dts.min, dts.sec, base)
        out[i] = apply_mult(out[i], mult)
    return out

def periodarr_to_dt64arr(ndarray[int64_t] periodarr, int base, long mult):
    """
    Convert array to datetime64 values from a set of ordinals corresponding to
    periods per period convention.
    """
    cdef:
        ndarray[int64_t] out
        Py_ssize_t i, l

    l = len(periodarr)

    out = np.empty(l, dtype='i8')

    for i in range(l):
        out[i] = period_ordinal_to_dt64(periodarr[i], base, mult)

    return out

cpdef long period_asfreq(long period_ordinal, int base1, long mult1,
                           int base2, long mult2, object relation='E'):
    """
    Convert period ordinal from one frequency to another, and if upsampling,
    choose to use start ('S') or end ('E') of period.
    """
    cdef:
        long retval

    if relation not in ('S', 'E'):
        raise ValueError('relation argument must be one of S or E')

    period_ordinal = remove_mult(period_ordinal, mult1)

    if mult1 != 1 and relation == 'E':
        period_ordinal += (mult1 - 1)

    retval = asfreq(period_ordinal, base1, base2, (<char*>relation)[0])
    retval = apply_mult(retval, mult2)

    return retval

def period_asfreq_arr(ndarray[int64_t] arr, int base1, long mult1, int base2,
                        long mult2, object relation='E'):
    """
    Convert int64-array of period ordinals from one frequency to another, and if
    upsampling, choose to use start ('S') or end ('E') of period.
    """
    cdef:
        ndarray[int64_t] new_arr
        Py_ssize_t i, sz

    if relation not in ('S', 'E'):
        raise ValueError('relation argument must be one of S or E')

    sz = len(arr)
    new_arr = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        new_arr[i] = period_asfreq(arr[i], base1, mult1, base2, mult2, relation)

    return new_arr

def period_ordinal(int y, int m, int d, int h, int min, int s,
                   int base, long mult):
    cdef:
        long ordinal

    ordinal = get_period_ordinal(y, m, d, h, min, s, base)

    return apply_mult(ordinal, mult)

cpdef int64_t period_ordinal_to_dt64(long period_ordinal, int base, long mult):
    cdef:
        long ordinal
        npy_datetimestruct dts
        date_info dinfo

    ordinal = remove_mult(period_ordinal, mult)

    get_date_info(ordinal, base, &dinfo)

    dts.year = dinfo.year
    dts.month = dinfo.month
    dts.day = dinfo.day
    dts.hour = dinfo.hour
    dts.min = dinfo.minute
    dts.sec = int(dinfo.second)
    dts.us = 0

    return PyArray_DatetimeStructToDatetime(NPY_FR_us, &dts)

def period_ordinal_to_string(long value, int base, long mult):
    cdef:
        char *ptr

    ptr = period_to_string(remove_mult(value, mult), base)

    if ptr == NULL:
        raise ValueError("Could not create string from ordinal '%d'" % value)

    return <object>ptr

def period_strftime(long value, int freq, long mult, object fmt):
    cdef:
        char *ptr

    value = remove_mult(value, mult)
    ptr = period_to_string2(value, freq, <char*>fmt)

    if ptr == NULL:
        raise ValueError("Could not create string with fmt '%s'" % fmt)

    return <object>ptr

# period accessors

ctypedef int (*accessor)(long ordinal, int base) except -1

cdef int apply_accessor(accessor func, long value, int base,
                        long mult) except -1:
    value = remove_mult(value, mult)
    return func(value, base)

cpdef int get_period_year(long value, int base, long mult) except -1:
    return apply_accessor(pyear, value, base, mult)

cpdef int get_period_qyear(long value, int base, long mult) except -1:
    return apply_accessor(pqyear, value, base, mult)

cpdef int get_period_quarter(long value, int base, long mult) except -1:
    return apply_accessor(pquarter, value, base, mult)

cpdef int get_period_month(long value, int base, long mult) except -1:
    return apply_accessor(pmonth, value, base, mult)

cpdef int get_period_day(long value, int base, long mult) except -1:
    return apply_accessor(pday, value, base, mult)

cpdef int get_period_hour(long value, int base, long mult) except -1:
    return apply_accessor(phour, value, base, mult)

cpdef int get_period_minute(long value, int base, long mult) except -1:
    return apply_accessor(pminute, value, base, mult)

cpdef int get_period_second(long value, int base, long mult) except -1:
    return apply_accessor(psecond, value, base, mult)

cpdef int get_period_dow(long value, int base, long mult) except -1:
    return apply_accessor(pday_of_week, value, base, mult)

cpdef int get_period_week(long value, int base, long mult) except -1:
    return apply_accessor(pweek, value, base, mult)

cpdef int get_period_weekday(long value, int base, long mult) except -1:
    return apply_accessor(pweekday, value, base, mult)

cpdef int get_period_doy(long value, int base, long mult) except -1:
    return apply_accessor(pday_of_year, value, base, mult)

# same but for arrays

cdef ndarray[int64_t] apply_accessor_arr(accessor func,
                                         ndarray[int64_t] arr,
                                         int base, long mult):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] out

    sz = len(arr)
    out = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        out[i] = remove_mult(arr[i], mult)
        out[i] = func(out[i], base)

    return out

def get_period_year_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pyear, arr, base, mult)

def get_period_qyear_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pqyear, arr, base, mult)

def get_period_quarter_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pquarter, arr, base, mult)

def get_period_month_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pmonth, arr, base, mult)

def get_period_day_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pday, arr, base, mult)

def get_period_hour_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(phour, arr, base, mult)

def get_period_minute_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pminute, arr, base, mult)

def get_period_second_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(psecond, arr, base, mult)

def get_period_dow_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pday_of_week, arr, base, mult)

def get_period_week_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pweek, arr, base, mult)

def get_period_weekday_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pweekday, arr, base, mult)

def get_period_doy_arr(ndarray[int64_t] arr, int base, long mult):
    return apply_accessor_arr(pday_of_year, arr, base, mult)

def get_abs_time(freq, dailyDate, originalDate):
    return getAbsTime(freq, dailyDate, originalDate)
