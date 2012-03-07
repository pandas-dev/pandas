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

# Python front end to C extension type _Timestamp
# This serves as the box for datetime64
class Timestamp(_Timestamp):
    def __new__(cls, object ts_input, object offset=None, tzinfo=None):
        ts = convert_to_tsobject(ts_input, tzinfo)

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

    def __setstate__(self, state):
        self.value = state[0]
        self.offset = state[1]
        self.tzinfo = state[2]

    def __reduce__(self):
        object_state = self.value, self.offset, self.tzinfo
        return (Timestamp, object_state)


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
            return super(_Timestamp, self).__add__(other)

    def __sub__(self, other):
        if is_integer_object(other):
            return self.__add__(-other)
        else:
            return super(_Timestamp, self).__sub__(other)

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
    # pretty cheap
    elif isinstance(ts, _Timestamp):
        tmp = ts
        retval.value = tmp.value
        retval.dtval = tmp
    elif isinstance(ts, Delta):
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

#cdef convert_to_res(object res):
#    if res == 'microsecond':
#        return r_microsecond
#    if res == 'second':
#        return r_second
#    if res == 'minute':
#        return r_minute
#    if res == 'hour':
#        return r_hour
#    if res == 'day':
#        return r_day
#    if res == 'month':
#        return r_month
#    if res == 'year':
#        return r_year
#    return r_invalid

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

cdef ndarray[int64_t] _generate_range(_Offset offset, Py_ssize_t periods):
    """
    Generate timestamps according to offset.
    """
    cdef:
        Py_ssize_t i
        ndarray[int64_t] dtindex

    dtindex = np.empty(periods, np.int64)
    for i in range(periods):
        dtindex[i] = offset._ts()
        offset.next()
    return dtindex

cdef int64_t _count_range(_Offset offset, object end):
    """
    Count timestamps in range according to offset up to (and including)
    end time.
    """
    cdef:
        Py_ssize_t i=0
        _TSObject e

    e = convert_to_tsobject(end)
    while offset._ts() <= e.value:
        i += 1
        offset.next()
    return i

# Here's some frequency caching logic
# ----------------------------------------------------------------------------------

# cache size...
# daily 1850-2050 takes up 73049 x 8bytes / 2**20bytes => 0.56Mb of cache. seems ok
# double this at least, for indexer

cdef class DatetimeCache:
    """
    Holds a array of datetimes according to some regular offset rule, along
    with a int64 => Py_ssize_t hashtable to lookup array index values.
    """
    cdef:
        object start
        object end
        object buf
        object periods

        _Offset offset
        Int64HashTable indexer
        object is_dirty

        int64_t end_anchor

    def __init__(self, _Offset offset, object end=None, object periods=None):
        """
        Note, prefer 'periods' argument over 'end' for generating range.
        """

        self.offset = offset
        self.start = offset.ts
        self.end = convert_to_tsobject(end)
        self.buf = None
        self.periods = periods
        self.is_dirty = True

    cdef rebuild(self):
        """
        Rebuild cache so that start and end fall within the range of
        times generated.
        """
        cdef:
            int64_t periods
            Py_ssize_t i
            ndarray[int64_t] buf

        # TODO: minimize memory allocation by copying existing data

        if self.periods is not None:
            periods = self.periods
            if self.buf is not None and periods < len(self.buf):
                periods = len(self.buf)
        else:
            periods = self.count()

        self.offset.anchor(self.start)
        buf = _generate_range(self.offset, periods)

        self.end = convert_to_tsobject(buf[-1])
        self.buf = buf

        self.offset.prev()
        self.end_anchor = self.offset._get_anchor()

        self.indexer = Int64HashTable(size_hint=periods)
        for i in range(periods):
            self.indexer.set_item(buf[i], i)

        self.is_dirty = False

    cdef int search(self, object ts, object side='left'):
        cdef _TSObject t = convert_to_tsobject(ts)
        return np.searchsorted(self.buf, t, side=side)

    cdef set_start(self, _Offset off):
        self.start = convert_to_tsobject(off._get_anchor())
        self.is_dirty = True

    cdef set_end(self, _Offset off):
        self.end = off._ts()
        self.is_dirty = True

    cdef set_periods(self, object periods):
        self.periods = periods
        self.is_dirty = True

    cdef int _lookup(self, int64_t val):
        cdef:
            kh_int64_t *table = self.indexer.table
            cdef khiter_t k

        k = kh_get_int64(table, val)
        if k != table.n_buckets:
            return table.vals[k]
        else:
            return -1

    cpdef ndarray[int64_t] extend(self, int64_t first, int64_t last, int n=0):
        """
        Extend cache to at least n periods beyond first and last
        """
        cdef:
            _Offset offset
            Py_ssize_t i, j
            int an

        an = abs(n)
        offset = self.offset

        offset.anchor()
        # if first is before current start
        if offset._ts() > first:
            # move back until on or just past first
            while offset._ts() > first:
                offset.prev()
            # move back an additional n periods
            for i in range(an):
                offset.prev()
            self.set_start(offset)
        # if first is after current start
        else:
            # move forward up to n periods until on or just past first
            for i in range(an):
                if offset._ts() >= first:
                    break
                offset.next()
            # now move back n periods
            for i in range(an):
                offset.prev()
            # are we earlier than start?
            if offset._ts() < self.start.value:
                self.set_start(offset)

        offset.anchor(self.end_anchor)
        # if last is after current end
        if offset._ts() < last:
            # move forward until on or just past last
            while offset._ts() < last:
                offset.next()
            # move forward an additional n periods
            for i in range(an):
                offset.next()
            self.set_end(offset)
        # if last is before current end
        else:
            # move back up to n periods until on or just past last
            for i in range(an):
                if offset._ts() <= last:
                    break
                offset.prev()
            # move forward n periods
            for j in range(an):
                offset.next()
            # are we further than end?
            if offset._ts() > self.end.value:
                self.set_end(offset)

        if self.is_dirty:
            self.rebuild()

        return self.buf

    # user/python-accessible methods

    cpdef Py_ssize_t count(self):
        if not self.is_dirty:
            return len(self.buf)

        self.offset.anchor(self.start)
        return _count_range(self.offset, self.end)

    cpdef lookup(self, object tslike):
        cdef:
            _TSObject ts = convert_to_tsobject(tslike)
            int64_t idx

        if ts.value < self.start.value:
            self.extend(ts.value, self.end.value)

        if ts.value > self.end.value:
            self.extend(self.start.value, ts.value)

        idx = self._lookup(ts.value)
        if idx < 0:
            raise KeyError(ts.value)
        else:
            return idx

    cpdef ndarray[int64_t] cache(self):
        return self.buf

_DEFAULT_BEGIN = Timestamp(datetime(1850, 1, 1))
_DEFAULT_END = Timestamp(datetime(2050, 1, 1))

_tcaches = {}

_months = {
    'JAN' : 1,
    'FEB' : 2,
    'MAR' : 3,
    'APR' : 4,
    'MAY' : 5,
    'JUN' : 6,
    'JUL' : 7,
    'AUG' : 8,
    'SEP' : 9,
    'OCT' : 10,
    'NOV' : 11,
    'DEC' : 12
}

_quarters = {
    'JAN' : 0,
    'FEB' : 1,
    'MAR' : 2
}

_weekdays = {
    'MON' : 0,
    'TUE' : 1,
    'WED' : 2,
    'THU' : 3,
    'FRI' : 4
}

# first two letters of frequency determines major rank for considering up- or
# down-sampling
_freqrank = {
    'DA' : 5,  # upsampling
    'WE' : 4,
    'W@' : 3,
    'EO' : 2,
    'Q@' : 1,
    'A@' : 0,  # downsampling
}

def flush_tcache(object freq):
    if freq in _tcaches:
        del _tcaches[freq]

# TODO: missing feature, user-provided striding

cpdef DatetimeCache get_tcache(freq, object first=None, object last=None):
    """
    Retrieve from cache (or generate, first time through) times that correspond
    to the frequency we care about.

    First and last allow overriding default cache range.
    """
    cdef:
        _Offset offset
        int64_t s, e
        DatetimeCache tc

    if first is None:
        first = _DEFAULT_BEGIN

    if last is None:
        last  = _DEFAULT_END

    if freq not in _tcaches:
        if freq == 'WEEKDAY':
            offset = DayOffset(stride=1, biz=1)

        elif freq == 'DAILY':
            offset = DayOffset(stride=1, biz=0)

        elif freq == 'EOM':
            offset = MonthOffset(dayoffset=-1, biz=-1)

        elif freq.startswith('W@'):
            wkd = first.weekday()
            dow = freq[-3:]
            if dow not in _weekdays:
                raise ValueError('Bad weekday %s' % freq)
            first += timedelta(days=(_weekdays[dow] - wkd) % 7)
            offset = DayOffset(stride=7, biz=1)

        elif freq.startswith('Q@'):
            mo = freq[-3:]
            if mo not in _quarters:
                raise ValueError('Bad month quarter %s' % freq)
            first += Delta(months=_quarters[mo])
            offset = MonthOffset(dayoffset=-1, stride=3, biz=-1)

        elif freq.startswith('A@'):
            mo = freq[-3:]
            if mo not in _months:
                raise ValueError('Bad annual month in %s' % freq)
            if _months[mo] < 12:
                first += Delta(months=_months[mo])
            offset = YearOffset(dayoffset=-1, biz=-1)

        elif freq.startswith('WOM@'):
            offset = freq.split('@')[1]
            week = int(offset[:-3])
            dow = offset[-3:]
            if dow not in _weekdays:
                raise ValueError('Bad weekday in %s' % freq)
            offset = DayOfMonthOffset(week=week, day=_weekdays[dow])

        else:
            raise ValueError('Supplied frequency %s not implemented' % freq)

        offset.anchor(first)
        last = Timestamp(last)

        _tcaches[freq] = DatetimeCache(offset, last.value)

    tc = _tcaches[freq]

    if tc.is_dirty:
        tc.rebuild()

    return _tcaches[freq]

# Helper methods for frequency-based analysis
# -------------------------------------------------------------

@cython.wraparound(False)
def conformity_check(ndarray[int64_t] data, object freq):
    """
    Return first non-conforming time, otherwise None. Also, return
    whether all times are concecutive according to the frequency.
    """
    cdef:
        Py_ssize_t i, ld, lc
        int idx, previdx
        object regular
        ndarray[int64_t] cache
        DatetimeCache tc

    regular = True
    ld = len(data)

    if ld == 0:
        return None, regular

    tc = get_tcache(freq)
    cache = tc.cache()

    lc = len(cache)

    # make sure cache is large enough to handle data
    if data[0] < cache[0] or data[ld-1] > cache[lc-1]:
        tc.extend(data[0], data[ld-1])

    previdx = tc._lookup(data[0])
    if previdx == -1:
        return data[0], False

    i = 1
    while i < ld:
        idx = tc._lookup(data[i])
        if idx == -1:
            return data[i], False
        if previdx != (idx - 1):
            regular = False
        i += 1

    return None, regular

@cython.wraparound(False)
cpdef ndarray[int64_t] fast_shift(ndarray[int64_t] data, object freq, int64_t n):
    """
    Shift times n periods according to the frequency.
    """
    cdef:
        DatetimeCache tc
        ndarray[int64_t] result, cache
        Py_ssize_t i, ld, lc
        int idx, s, e

    ld = len(data)

    tc = get_tcache(freq)
    cache = tc.cache()

    lc = len(cache)

    # make sure cache is large enough to handle data
    # plus a shift of N
    s = cache.searchsorted(data[0])
    e = cache.searchsorted(data[ld-1])

    if (data[0] < cache[0] or s + n < 0 or
        data[ld-1] > cache[lc-1] or e + n > lc):
        cache = tc.extend(data[0], data[ld-1], n)

    result = np.empty(ld, dtype='i8')
    for i in range(ld):
        idx = tc._lookup(data[i]) + n
        if idx == -1:
            raise ValueError("Shift hit nonconforming time %s" 
                             % np.datetime64(data[i]))
        result[i] = cache[idx]

    return result


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
        if not isinstance(self, Delta):
            tmp = self
            self = other
            other = tmp

        if isinstance(other, Delta):
            return self._add_delta(other)

        if isinstance(other, Timestamp):
            return self._add_timestamp(other)

        if isinstance(other, datetime):
            return self._add_timestamp(Timestamp(other))

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
               + timedelta(days=days,
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
            ret += timedelta(days=jumpdays)

        return Timestamp(ret)

    def __richcmp__(self, other, op):
        if op == 3:
            return not (self == other)
        elif op == 2:
            if not isinstance(other, Delta):
                return False
            if self.weekday or other.weekday:
                if not self.weekday or not other.weekday:
                    return False
                if self.weekday.weekday != other.weekday.weekday:
                    return False
                n1, n2 = self.weekday.n, other.weekday.n
                if n1 != n2 and not ((not n1 or n1 == 1) and (not n2 or n2 == 1)):
                    return False
            return (self.years == other.years and
                    self.months == other.months and
                    self.days == other.days and
                    self.hours == other.hours and
                    self.minutes == other.minutes and
                    self.seconds == other.seconds and
                    self.leapdays == other.leapdays and
                    self.year == other.year and
                    self.month == other.month and
                    self.day == other.day and
                    self.hour == other.hour and
                    self.minute == other.minute and
                    self.second == other.second and
                    self.microsecond == other.microsecond)
        else:
            raise NotImplementedError("Delta doesn't implement that comparison")

    def _add_delta(self, other):
        if not isinstance(self, Delta):
            tmp = self
            self = other
            other = tmp

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
        if not isinstance(self, Delta):
            tmp = self
            self = other
            other = tmp

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


    def __mul__(self, v):
        cdef int64_t f

        if not isinstance(self, Delta):
            tmp = self
            self = v
            v = tmp

        f = v

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
            out[i] = _month_offset[isleap][dts.month-1] + dts.day
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
            out[i] = _month_offset[isleap][dts.month-1] + dts.day
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

# Interval logic
# ------------------------------------------------------------------------------

def interval_freq_conv(int dtordinal, int freq1, int freq2, char relation='E'):
    cdef:
        int retval

    retval = frequency_conversion(dtordinal, freq1, freq2, relation)

    return retval

def skts_ordinal(int y, int m, int d, int h, int min, int s, int freq):
    cdef:
        long ordinal

    ordinal = get_skts_ordinal(y, m, d, h, min, s, freq)

    return ordinal

def skts_ordinal_to_dt(long skts_ordinal, int freq):
    return datetime.fromordinal(get_python_ordinal(skts_ordinal, freq))
