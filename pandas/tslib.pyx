# cython: profile=False

cimport numpy as np
from numpy cimport (int32_t, int64_t, import_array, ndarray,
                    NPY_INT64, NPY_DATETIME)
import numpy as np

from cpython cimport *

# Cython < 0.17 doesn't have this in cpython
cdef extern from "Python.h":
    cdef PyTypeObject *Py_TYPE(object)


from libc.stdlib cimport free

from util cimport is_integer_object, is_datetime64_object
cimport util

from datetime cimport *
from khash cimport *
cimport cython

from datetime import timedelta, datetime
from dateutil.parser import parse as parse_date

cdef extern from "Python.h":
    int PySlice_Check(object)

# initialize numpy
import_array()
#import_ufunc()

# import datetime C API
PyDateTime_IMPORT

# in numpy 1.7, will prob need the following:
# numpy_pydatetime_import

cdef int64_t NPY_NAT = util.get_nat()


try:
    basestring
except NameError: # py3
    basestring = str

def ints_to_pydatetime(ndarray[int64_t] arr, tz=None):
    cdef:
        Py_ssize_t i, n = len(arr)
        pandas_datetimestruct dts
        ndarray[object] result = np.empty(n, dtype=object)

    if tz is not None:
        if _is_utc(tz):
            for i in range(n):
                pandas_datetime_to_datetimestruct(arr[i], PANDAS_FR_ns, &dts)
                result[i] = datetime(dts.year, dts.month, dts.day, dts.hour,
                                     dts.min, dts.sec, dts.us, tz)
        elif _is_tzlocal(tz) or _is_fixed_offset(tz):
            for i in range(n):
                pandas_datetime_to_datetimestruct(arr[i], PANDAS_FR_ns, &dts)
                dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                              dts.min, dts.sec, dts.us, tz)
                result[i] = dt + tz.utcoffset(dt)
        else:
            trans = _get_transitions(tz)
            deltas = _get_deltas(tz)
            for i in range(n):
                # Adjust datetime64 timestamp, recompute datetimestruct
                pos = trans.searchsorted(arr[i]) - 1
                inf = tz._transition_info[pos]

                pandas_datetime_to_datetimestruct(arr[i] + deltas[pos],
                                                  PANDAS_FR_ns, &dts)
                result[i] = datetime(dts.year, dts.month, dts.day, dts.hour,
                                     dts.min, dts.sec, dts.us,
                                     tz._tzinfos[inf])
    else:
        for i in range(n):
            pandas_datetime_to_datetimestruct(arr[i], PANDAS_FR_ns, &dts)
            result[i] = datetime(dts.year, dts.month, dts.day, dts.hour,
                                 dts.min, dts.sec, dts.us)

    return result

from dateutil.tz import tzlocal

def _is_tzlocal(tz):
    return isinstance(tz, tzlocal)

def _is_fixed_offset(tz):
    try:
        tz._transition_info
        return False
    except AttributeError:
        return True

# Python front end to C extension type _Timestamp
# This serves as the box for datetime64
class Timestamp(_Timestamp):

    def __new__(cls, object ts_input, object offset=None, tz=None):
        cdef _TSObject ts
        cdef _Timestamp ts_base

        if PyFloat_Check(ts_input):
            # to do, do we want to support this, ie with fractional seconds?
            raise TypeError("Cannot convert a float to datetime")

        if util.is_string_object(ts_input):
            try:
                ts_input = parse_date(ts_input)
            except Exception:
                pass

        ts = convert_to_tsobject(ts_input, tz)

        if ts.value == NPY_NAT:
            return NaT

        # make datetime happy
        ts_base = _Timestamp.__new__(cls, ts.dts.year, ts.dts.month,
                                     ts.dts.day, ts.dts.hour, ts.dts.min,
                                     ts.dts.sec, ts.dts.us, ts.tzinfo)

        # fill out rest of data
        ts_base.value = ts.value
        ts_base.offset = offset
        ts_base.nanosecond = ts.dts.ps / 1000

        return ts_base

    def __repr__(self):
        result = self._repr_base

        try:
            result += self.strftime('%z')
            if self.tzinfo:
                zone = _get_zone(self.tzinfo)
                result += _tz_format(self, zone)
        except ValueError:
            year2000 = self.replace(year=2000)
            result += year2000.strftime('%z')
            if self.tzinfo:
                zone = _get_zone(self.tzinfo)
                result += _tz_format(year2000, zone)

        return '<Timestamp: %s>' % result

    @property
    def _repr_base(self):
        result = '%d-%.2d-%.2d %.2d:%.2d:%.2d' % (self.year, self.month,
                                                  self.day, self.hour,
                                                  self.minute, self.second)

        if self.nanosecond != 0:
            nanos = self.nanosecond + 1000 * self.microsecond
            result += '.%.9d' % nanos
        elif self.microsecond != 0:
            result += '.%.6d' % self.microsecond

        return result

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

        if freq is None:
            freq = self.freq

        return Period(self, freq=freq)

    @property
    def dayofweek(self):
        return self.weekday()

    @property
    def dayofyear(self):
        return self._get_field('doy')

    @property
    def week(self):
        return self._get_field('woy')

    weekofyear = week

    @property
    def quarter(self):
        return self._get_field('q')

    @property
    def freqstr(self):
        return getattr(self.offset, 'freqstr', self.offset)

    @property
    def asm8(self):
        return np.int64(self.value).view('M8[ns]')

    def tz_localize(self, tz):
        """
        Convert naive Timestamp to local time zone

        Parameters
        ----------
        tz : pytz.timezone

        Returns
        -------
        localized : Timestamp
        """
        if self.tzinfo is None:
            # tz naive, localize
            return Timestamp(self.to_pydatetime(), tz=tz)
        else:
            raise Exception('Cannot localize tz-aware Timestamp, use '
                            'tz_convert for conversions')

    def tz_convert(self, tz):
        """
        Convert Timestamp to another time zone or localize to requested time
        zone

        Parameters
        ----------
        tz : pytz.timezone

        Returns
        -------
        converted : Timestamp
        """
        if self.tzinfo is None:
            # tz naive, use tz_localize
            raise Exception('Cannot convert tz-naive Timestamp, use '
                            'tz_localize to localize')
        else:
            # Same UTC timestamp, different time zone
            return Timestamp(self.value, tz=tz)

    astimezone = tz_convert

    def replace(self, **kwds):
        return Timestamp(datetime.replace(self, **kwds),
                         offset=self.offset)

    def to_pydatetime(self, warn=True):
        """
        If warn=True, issue warning if nanoseconds is nonzero
        """
        cdef:
            pandas_datetimestruct dts
            _TSObject ts

        if self.nanosecond != 0 and warn:
            print 'Warning: discarding nonzero nanoseconds'
        ts = convert_to_tsobject(self, self.tzinfo)

        return datetime(ts.dts.year, ts.dts.month, ts.dts.day,
                        ts.dts.hour, ts.dts.min, ts.dts.sec,
                        ts.dts.us, ts.tzinfo)


class NaTType(_NaT):

    def __new__(cls):
        cdef _NaT base

        base = _NaT.__new__(cls, 1, 1, 1)
        mangle_nat(base)
        base.value = NPY_NAT

        return base

    def __repr__(self):
        return 'NaT'

    def weekday(self):
        return -1

    def toordinal(self):
        return -1

fields = ['year', 'quarter', 'month', 'day', 'hour',
          'minute', 'second', 'microsecond', 'nanosecond',
          'week', 'dayofyear']
for field in fields:
    prop = property(fget=lambda self: -1)
    setattr(NaTType, field, prop)


NaT = NaTType()

iNaT = util.get_nat()

cdef _tz_format(object obj, object zone):
    try:
        return obj.strftime(' %%Z, tz=%s' % zone)
    except:
        return ', tz=%s' % zone

def is_timestamp_array(ndarray[object] values):
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not is_timestamp(values[i]):
            return False
    return True


cpdef object get_value_box(ndarray arr, object loc):
    cdef:
        Py_ssize_t i, sz
        void* data_ptr
    if util.is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    sz = np.PyArray_SIZE(arr)

    if i < 0 and sz > 0:
        i += sz
    elif i >= sz or sz == 0:
        raise IndexError('index out of bounds')

    if arr.descr.type_num == NPY_DATETIME:
        return Timestamp(util.get_value_1d(arr, i))
    else:
        return util.get_value_1d(arr, i)


#----------------------------------------------------------------------
# Frequency inference

def unique_deltas(ndarray[int64_t] arr):
    cdef:
        Py_ssize_t i, n = len(arr)
        int64_t val
        khiter_t k
        kh_int64_t *table
        int ret = 0
        list uniques = []

    table = kh_init_int64()
    kh_resize_int64(table, 10)
    for i in range(n - 1):
        val = arr[i + 1] - arr[i]
        k = kh_get_int64(table, val)
        if k == table.n_buckets:
            kh_put_int64(table, val, &ret)
            uniques.append(val)
    kh_destroy_int64(table)

    result = np.array(uniques, dtype=np.int64)
    result.sort()
    return result


cdef inline bint _is_multiple(int64_t us, int64_t mult):
    return us % mult == 0


def apply_offset(ndarray[object] values, object offset):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] new_values
        object boxed

    result = np.empty(n, dtype='M8[ns]')
    new_values = result.view('i8')
    pass


# This is PITA. Because we inherit from datetime, which has very specific
# construction requirements, we need to do object instantiation in python
# (see Timestamp class above). This will serve as a C extension type that
# shadows the python class, where we do any heavy lifting.
cdef class _Timestamp(datetime):
    cdef readonly:
        int64_t value, nanosecond
        object offset       # frequency reference

    def __hash__(self):
        if self.nanosecond:
            return hash(self.value)
        else:
            return datetime.__hash__(self)

    def __richcmp__(_Timestamp self, object other, int op):
        cdef _Timestamp ots

        if isinstance(other, _Timestamp):
            ots = other
        elif type(other) is datetime:
            if self.nanosecond == 0:
                val = self.to_datetime()
                return PyObject_RichCompareBool(val, other, op)

            try:
                ots = Timestamp(other)
            except ValueError:
                return self._compare_outside_nanorange(other, op)
        else:
            if op == 2:
                return False
            elif op == 3:
                return True
            else:
                raise TypeError('Cannot compare Timestamp with %s' % str(other))

        self._assert_tzawareness_compat(other)

        if op == 2: # ==
            return self.value == ots.value
        elif op == 3: # !=
            return self.value != ots.value
        elif op == 0: # <
            return self.value < ots.value
        elif op == 1: # <=
            return self.value <= ots.value
        elif op == 4: # >
            return self.value > ots.value
        elif op == 5: # >=
            return self.value >= ots.value

    cdef _compare_outside_nanorange(self, object other, int op):
        dtval = self.to_datetime()

        self._assert_tzawareness_compat(other)

        if self.nanosecond == 0:
            if op == 2: # ==
                return dtval == other
            elif op == 3: # !=
                return dtval != other
            elif op == 0: # <
                return dtval < other
            elif op == 1: # <=
                return dtval <= other
            elif op == 4: # >
                return dtval > other
            elif op == 5: # >=
                return dtval >= other
        else:
            if op == 2: # ==
                return False
            elif op == 3: # !=
                return True
            elif op == 0: # <
                return dtval < other
            elif op == 1: # <=
                return dtval < other
            elif op == 4: # >
                return dtval >= other
            elif op == 5: # >=
                return dtval >= other

    cdef _assert_tzawareness_compat(self, object other):
        if self.tzinfo is None:
            if other.tzinfo is not None:
                raise Exception('Cannot compare tz-naive and '
                                'tz-aware timestamps')
        elif other.tzinfo is None:
            raise Exception('Cannot compare tz-naive and tz-aware timestamps')

    cpdef to_datetime(self):
        cdef:
            pandas_datetimestruct dts
            _TSObject ts
        ts = convert_to_tsobject(self, self.tzinfo)
        dts = ts.dts
        return datetime(dts.year, dts.month, dts.day,
                        dts.hour, dts.min, dts.sec,
                        dts.us, ts.tzinfo)

    def __add__(self, other):
        if is_integer_object(other):
            if self.offset is None:
                msg = ("Cannot add integral value to Timestamp "
                       "without offset.")
                raise ValueError(msg)
            else:
                return Timestamp((self.offset.__mul__(other)).apply(self))
        else:
            if isinstance(other, timedelta) or hasattr(other, 'delta'):
                nanos = _delta_to_nanoseconds(other)
                return Timestamp(self.value + nanos, tz=self.tzinfo)
            else:
                result = datetime.__add__(self, other)
                if isinstance(result, datetime):
                    result = Timestamp(result)
                    result.nanosecond = self.nanosecond
                return result

    def __sub__(self, other):
        if is_integer_object(other):
            return self.__add__(-other)
        else:
            return datetime.__sub__(self, other)

    cpdef _get_field(self, field):
        out = get_date_field(np.array([self.value], dtype=np.int64), field)
        return out[0]


cdef PyTypeObject* ts_type = <PyTypeObject*> Timestamp


cdef inline bint is_timestamp(object o):
    return Py_TYPE(o) == ts_type # isinstance(o, Timestamp)


cdef class _NaT(_Timestamp):

    def __richcmp__(_NaT self, object other, int op):
        # if not isinstance(other, (_NaT, _Timestamp)):
        #     raise TypeError('Cannot compare %s with NaT' % type(other))

        if op == 2: # ==
            return False
        elif op == 3: # !=
            return True
        elif op == 0: # <
            return False
        elif op == 1: # <=
            return False
        elif op == 4: # >
            return False
        elif op == 5: # >=
            return False




def _delta_to_nanoseconds(delta):
    try:
        delta = delta.delta
    except:
        pass
    return (delta.days * 24 * 60 * 60 * 1000000
            + delta.seconds * 1000000
            + delta.microseconds) * 1000


# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    cdef:
        pandas_datetimestruct dts      # pandas_datetimestruct
        int64_t value               # numpy dt64
        object tzinfo

    property value:
        def __get__(self):
            return self.value

cpdef _get_utcoffset(tzinfo, obj):
    try:
        return tzinfo._utcoffset
    except AttributeError:
        return tzinfo.utcoffset(obj)

# helper to extract datetime and int64 from several different possibilities
cdef convert_to_tsobject(object ts, object tz):
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
        _TSObject obj
        bint utc_convert = 1

    if tz is not None:
        if isinstance(tz, basestring):
            tz = pytz.timezone(tz)

    obj = _TSObject()

    if is_datetime64_object(ts):
        obj.value = _get_datetime64_nanos(ts)
        pandas_datetime_to_datetimestruct(obj.value, PANDAS_FR_ns, &obj.dts)
    elif is_integer_object(ts):
        obj.value = ts
        pandas_datetime_to_datetimestruct(ts, PANDAS_FR_ns, &obj.dts)
    elif util.is_string_object(ts):
        _string_to_dts(ts, &obj.dts)
        obj.value = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &obj.dts)
    elif PyDateTime_Check(ts):
        if tz is not None:
            # sort of a temporary hack
            if ts.tzinfo is not None:
                if hasattr(tz, 'normalize'):
                    ts = tz.normalize(ts)
                    obj.value = _pydatetime_to_dts(ts, &obj.dts)
                    obj.tzinfo = ts.tzinfo
                else: #tzoffset
                    ts_offset = _get_utcoffset(ts.tzinfo, ts)
                    obj.value = _pydatetime_to_dts(ts, &obj.dts)
                    obj.value -= _delta_to_nanoseconds(ts_offset)
                    tz_offset = _get_utcoffset(tz, ts)
                    obj.value += _delta_to_nanoseconds(tz_offset)

                    obj.tzinfo = tz
            elif not _is_utc(tz):
                try:
                    ts = tz.localize(ts)
                except AttributeError:
                    ts = ts.replace(tzinfo=tz)

                obj.value = _pydatetime_to_dts(ts, &obj.dts)
                offset = _get_utcoffset(ts.tzinfo, ts)
                obj.value -= _delta_to_nanoseconds(offset)
                obj.tzinfo = ts.tzinfo
            else:
                # UTC
                obj.value = _pydatetime_to_dts(ts, &obj.dts)
                obj.tzinfo = pytz.utc
        else:
            obj.value = _pydatetime_to_dts(ts, &obj.dts)
            obj.tzinfo = ts.tzinfo
            if obj.tzinfo is not None and not _is_utc(obj.tzinfo):
                offset = _get_utcoffset(obj.tzinfo, ts)
                obj.value -= _delta_to_nanoseconds(offset)

        if is_timestamp(ts):
            obj.value += ts.nanosecond
        _check_dts_bounds(obj.value, &obj.dts)
        return obj
    elif PyDate_Check(ts):
        obj.value  = _date_to_datetime64(ts, &obj.dts)
    else:
        raise ValueError("Could not construct Timestamp from argument %s" %
                         type(ts))

    if obj.value != NPY_NAT:
        _check_dts_bounds(obj.value, &obj.dts)

    if tz is not None:
        _localize_tso(obj, tz)

    return obj

cdef inline void _localize_tso(_TSObject obj, object tz):
    if _is_utc(tz):
        obj.tzinfo = tz
    elif _is_tzlocal(tz):
        pandas_datetime_to_datetimestruct(obj.value, PANDAS_FR_ns, &obj.dts)
        dt = datetime(obj.dts.year, obj.dts.month, obj.dts.day, obj.dts.hour,
                      obj.dts.min, obj.dts.sec, obj.dts.us, tz)
        delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
        pandas_datetime_to_datetimestruct(obj.value + delta,
                                          PANDAS_FR_ns, &obj.dts)
        obj.tzinfo = tz
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans = _get_transitions(tz)
        deltas = _get_deltas(tz)
        pos = trans.searchsorted(obj.value, side='right') - 1

        # statictzinfo
        if not hasattr(tz, '_transition_info'):
            pandas_datetime_to_datetimestruct(obj.value + deltas[0],
                                              PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz
        else:
            inf = tz._transition_info[pos]
            pandas_datetime_to_datetimestruct(obj.value + deltas[pos],
                                              PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz._tzinfos[inf]


def get_timezone(tz):
    return _get_zone(tz)

cdef inline bint _is_utc(object tz):
    return tz is UTC or isinstance(tz, _du_utc)

cdef inline object _get_zone(object tz):
    if _is_utc(tz):
        return 'UTC'
    else:
        try:
            zone = tz.zone
            if zone is None:
                return tz
            return zone
        except AttributeError:
            return tz

# cdef int64_t _NS_LOWER_BOUND = -9223285636854775809LL
# cdef int64_t _NS_UPPER_BOUND = -9223372036854775807LL

cdef inline _check_dts_bounds(int64_t value, pandas_datetimestruct *dts):
    cdef pandas_datetimestruct dts2
    if dts.year <= 1677 or dts.year >= 2262:
        pandas_datetime_to_datetimestruct(value, PANDAS_FR_ns, &dts2)
        if dts2.year != dts.year:
            fmt = '%d-%.2d-%.2d %.2d:%.2d:%.2d' % (dts.year, dts.month,
                                                   dts.day, dts.hour,
                                                   dts.min, dts.sec)

            raise ValueError('Out of bounds nanosecond timestamp: %s' % fmt)

# elif isinstance(ts, _Timestamp):
#     tmp = ts
#     obj.value = (<_Timestamp> ts).value
#     obj.dtval =
# elif isinstance(ts, object):
#     # If all else fails
#     obj.value = _dtlike_to_datetime64(ts, &obj.dts)
#     obj.dtval = _dts_to_pydatetime(&obj.dts)

def datetime_to_datetime64(ndarray[object] values):
    cdef:
        Py_ssize_t i, n = len(values)
        object val, inferred_tz = None
        ndarray[int64_t] iresult
        pandas_datetimestruct dts
        _TSObject _ts

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')
    for i in range(n):
        val = values[i]
        if util._checknull(val):
            iresult[i] = iNaT
        elif PyDateTime_Check(val):
            if val.tzinfo is not None:
                if inferred_tz is not None:
                    if _get_zone(val.tzinfo) != inferred_tz:
                        raise ValueError('Array must be all same time zone')
                else:
                    inferred_tz = _get_zone(val.tzinfo)

                _ts = convert_to_tsobject(val, None)
                iresult[i] = _ts.value
                _check_dts_bounds(iresult[i], &_ts.dts)
            else:
                if inferred_tz is not None:
                    raise ValueError('Cannot mix tz-aware with tz-naive values')
                iresult[i] = _pydatetime_to_dts(val, &dts)
                _check_dts_bounds(iresult[i], &dts)
        else:
            raise TypeError('Unrecognized value type: %s' % type(val))

    return result, inferred_tz


def array_to_datetime(ndarray[object] values, raise_=False, dayfirst=False,
                      format=None, utc=None):
    cdef:
        Py_ssize_t i, n = len(values)
        object val
        ndarray[int64_t] iresult
        ndarray[object] oresult
        pandas_datetimestruct dts
        bint utc_convert = bool(utc)
        _TSObject _ts

    from dateutil.parser import parse

    try:
        result = np.empty(n, dtype='M8[ns]')
        iresult = result.view('i8')
        for i in range(n):
            val = values[i]
            if util._checknull(val) or val is NaT:
                iresult[i] = iNaT
            elif PyDateTime_Check(val):
                if val.tzinfo is not None:
                    if utc_convert:
                        _ts = convert_to_tsobject(val, None)
                        iresult[i] = _ts.value
                        _check_dts_bounds(iresult[i], &_ts.dts)
                    else:
                        raise ValueError('Tz-aware datetime.datetime cannot '
                                         'be converted to datetime64 unless '
                                         'utc=True')
                else:
                    iresult[i] = _pydatetime_to_dts(val, &dts)
                    if is_timestamp(val):
                        iresult[i] += (<_Timestamp>val).nanosecond
                    _check_dts_bounds(iresult[i], &dts)
            elif PyDate_Check(val):
                iresult[i] = _date_to_datetime64(val, &dts)
                _check_dts_bounds(iresult[i], &dts)
            elif util.is_datetime64_object(val):
                iresult[i] = _get_datetime64_nanos(val)
            elif util.is_integer_object(val):
                iresult[i] = val
            else:
                if len(val) == 0:
                    iresult[i] = iNaT
                    continue

                try:
                    _string_to_dts(val, &dts)
                    iresult[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns,
                                                                   &dts)
                    _check_dts_bounds(iresult[i], &dts)
                except ValueError:
                    try:
                        result[i] = parse(val, dayfirst=dayfirst)
                    except Exception:
                        raise TypeError
                    pandas_datetime_to_datetimestruct(iresult[i], PANDAS_FR_ns,
                                                      &dts)
                    _check_dts_bounds(iresult[i], &dts)
        return result
    except TypeError:
        oresult = np.empty(n, dtype=object)

        for i in range(n):
            val = values[i]
            if util._checknull(val):
                oresult[i] = val
            else:
                if len(val) == 0:
                    # TODO: ??
                    oresult[i] = 'NaT'
                    continue
                try:
                    oresult[i] = parse(val, dayfirst=dayfirst)
                except Exception:
                    if raise_:
                        raise
                    return values
                    # oresult[i] = val

        return oresult

cdef inline _get_datetime64_nanos(object val):
    cdef:
        pandas_datetimestruct dts
        PANDAS_DATETIMEUNIT unit
        npy_datetime ival

    unit = get_datetime64_unit(val)
    if unit == 3:
        raise ValueError('NumPy 1.6.1 business freq not supported')

    ival = get_datetime64_value(val)

    if unit != PANDAS_FR_ns:
        pandas_datetime_to_datetimestruct(ival, unit, &dts)
        return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
    else:
        return ival


def cast_to_nanoseconds(ndarray arr):
    cdef:
        Py_ssize_t i, n = arr.size
        ndarray[int64_t] ivalues, iresult
        PANDAS_DATETIMEUNIT unit
        pandas_datetimestruct dts

    shape = (<object> arr).shape

    ivalues = arr.view(np.int64).ravel()

    result = np.empty(shape, dtype='M8[ns]')
    iresult = result.ravel().view(np.int64)

    unit = get_datetime64_unit(arr.flat[0])
    if unit == 3:
        raise ValueError('NumPy 1.6.1 business freq not supported')

    for i in range(n):
        pandas_datetime_to_datetimestruct(ivalues[i], unit, &dts)
        iresult[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)

    return result

#----------------------------------------------------------------------
# Conversion routines


def pydt_to_i8(object pydt):
    '''
    Convert to int64 representation compatible with numpy datetime64; converts
    to UTC
    '''
    cdef:
        _TSObject ts

    ts = convert_to_tsobject(pydt, None)

    return ts.value

def i8_to_pydt(int64_t i8, object tzinfo = None):
    '''
    Inverse of pydt_to_i8
    '''
    return Timestamp(i8)

#----------------------------------------------------------------------
# time zone conversion helpers

try:
    from dateutil.tz import tzutc as _du_utc
    import pytz
    UTC = pytz.utc
    have_pytz = True
except:
    have_pytz = False

def tz_convert(ndarray[int64_t] vals, object tz1, object tz2):
    cdef:
        ndarray[int64_t] utc_dates, result, trans, deltas
        Py_ssize_t i, pos, n = len(vals)
        int64_t v, offset
        pandas_datetimestruct dts

    if not have_pytz:
        import pytz

    # Convert to UTC

    if _get_zone(tz1) != 'UTC':
        utc_dates = np.empty(n, dtype=np.int64)
        if _is_tzlocal(tz1):
            for i in range(n):
                v = vals[i]
                pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
                dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                              dts.min, dts.sec, dts.us, tz1)
                delta = int(total_seconds(_get_utcoffset(tz1, dt))) * 1000000000
                utc_dates[i] = v - delta
        else:
            deltas = _get_deltas(tz1)
            trans = _get_transitions(tz1)
            pos = trans.searchsorted(vals[0]) - 1
            if pos < 0:
                raise ValueError('First time before start of DST info')

            offset = deltas[pos]
            for i in range(n):
                v = vals[i]
                if v >= [pos + 1]:
                    pos += 1
                    offset = deltas[pos]
                utc_dates[i] = v - offset
    else:
        utc_dates = vals

    if _get_zone(tz2) == 'UTC':
        return utc_dates

    result = np.empty(n, dtype=np.int64)
    if _is_tzlocal(tz2):
        for i in range(n):
            v = utc_dates[i]
            pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz2)
            delta = int(total_seconds(_get_utcoffset(tz2, dt))) * 1000000000
            result[i] = v + delta
            return result

    # Convert UTC to other timezone
    trans = _get_transitions(tz2)
    deltas = _get_deltas(tz2)
    pos = trans.searchsorted(utc_dates[0])
    if pos == 0:
        raise ValueError('First time before start of DST info')
    elif pos == len(trans):
        return utc_dates + deltas[-1]

    # TODO: this assumed sortedness :/
    pos -= 1

    offset = deltas[pos]
    for i in range(n):
        v = utc_dates[i]
        if v >= trans[pos + 1]:
            pos += 1
            offset = deltas[pos]
        result[i] = v + offset

    return result

def tz_convert_single(int64_t val, object tz1, object tz2):
    cdef:
        ndarray[int64_t] trans, deltas
        Py_ssize_t pos
        int64_t v, offset, utc_date
        pandas_datetimestruct dts

    if not have_pytz:
        import pytz

    # Convert to UTC
    if _is_tzlocal(tz1):
        pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz1)
        delta = int(total_seconds(_get_utcoffset(tz1, dt))) * 1000000000
        utc_date = val - delta
    elif _get_zone(tz1) != 'UTC':
        deltas = _get_deltas(tz1)
        trans = _get_transitions(tz1)
        pos = trans.searchsorted(val) - 1
        if pos < 0:
            raise ValueError('First time before start of DST info')
        offset = deltas[pos]
        utc_date = val - offset
    else:
        utc_date = val

    if _get_zone(tz2) == 'UTC':
        return utc_date
    if _is_tzlocal(tz2):
        pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz2)
        delta = int(total_seconds(_get_utcoffset(tz2, dt))) * 1000000000
        return utc_date + delta
    # Convert UTC to other timezone
    trans = _get_transitions(tz2)
    deltas = _get_deltas(tz2)
    pos = trans.searchsorted(utc_date) - 1
    if pos < 0:
        raise ValueError('First time before start of DST info')

    offset = deltas[pos]
    return utc_date + offset


trans_cache = {}
utc_offset_cache = {}

def _get_transitions(tz):
    """
    Get UTC times of DST transitions
    """
    try:
        # tzoffset not hashable in Python 3
        hash(tz)
    except TypeError:
        return np.array([NPY_NAT + 1], dtype=np.int64)

    if tz not in trans_cache:
        if hasattr(tz, '_utc_transition_times'):
            arr = np.array(tz._utc_transition_times, dtype='M8[ns]')
            arr = arr.view('i8')
            try:
                if tz._utc_transition_times[0].year == 1:
                    arr[0] = NPY_NAT + 1
            except Exception:
                pass
        else:
            arr = np.array([NPY_NAT + 1], dtype=np.int64)
        trans_cache[tz] = arr
    return trans_cache[tz]

def _get_deltas(tz):
    """
    Get UTC offsets in microseconds corresponding to DST transitions
    """
    try:
        # tzoffset not hashable in Python 3
        hash(tz)
    except TypeError:
        num = int(total_seconds(_get_utcoffset(tz, None))) * 1000000000
        return np.array([num], dtype=np.int64)

    if tz not in utc_offset_cache:
        if hasattr(tz, '_utc_transition_times'):
            utc_offset_cache[tz] = _unbox_utcoffsets(tz._transition_info)
        else:
            # static tzinfo
            num = int(total_seconds(_get_utcoffset(tz, None))) * 1000000000
            utc_offset_cache[tz] = np.array([num], dtype=np.int64)

    return utc_offset_cache[tz]

cdef double total_seconds(object td): # Python 2.6 compat
    return ((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) //
            10**6)

def tot_seconds(td):
    return total_seconds(td)

cpdef ndarray _unbox_utcoffsets(object transinfo):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] arr

    sz = len(transinfo)
    arr = np.empty(sz, dtype='i8')

    for i in range(sz):
        arr[i] = int(total_seconds(transinfo[i][0])) * 1000000000

    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, object tz):
    """
    Localize tzinfo-naive DateRange to given time zone (using pytz). If
    there are ambiguities in the values, raise AmbiguousTimeError.

    Returns
    -------
    localized : DatetimeIndex
    """
    cdef:
        ndarray[int64_t] trans, deltas, idx_shifted
        Py_ssize_t i, idx, pos, ntrans, n = len(vals)
        int64_t *tdata
        int64_t v, left, right
        ndarray[int64_t] result, result_a, result_b
        pandas_datetimestruct dts

    # Vectorized version of DstTzInfo.localize

    if not have_pytz:
        raise Exception("Could not find pytz module")

    if tz == UTC or tz is None:
        return vals

    result = np.empty(n, dtype=np.int64)

    if _is_tzlocal(tz):
        for i in range(n):
            v = vals[i]
            pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
            result[i] = v - delta
        return result

    trans = _get_transitions(tz)  # transition dates
    deltas = _get_deltas(tz)      # utc offsets

    tdata = <int64_t*> trans.data
    ntrans = len(trans)

    result_a = np.empty(n, dtype=np.int64)
    result_b = np.empty(n, dtype=np.int64)
    result_a.fill(NPY_NAT)
    result_b.fill(NPY_NAT)

    # left side
    idx_shifted = _ensure_int64(
        np.maximum(0, trans.searchsorted(vals - DAY_NS, side='right') - 1))

    for i in range(n):
        v = vals[i] - deltas[idx_shifted[i]]
        pos = bisect_right_i8(tdata, v, ntrans) - 1

        # timestamp falls to the left side of the DST transition
        if v + deltas[pos] == vals[i]:
            result_a[i] = v

    # right side
    idx_shifted = _ensure_int64(
        np.maximum(0, trans.searchsorted(vals + DAY_NS, side='right') - 1))

    for i in range(n):
        v = vals[i] - deltas[idx_shifted[i]]
        pos = bisect_right_i8(tdata, v, ntrans) - 1

        # timestamp falls to the right side of the DST transition
        if v + deltas[pos] == vals[i]:
            result_b[i] = v

    for i in range(n):
        left = result_a[i]
        right = result_b[i]
        if left != NPY_NAT and right != NPY_NAT:
            if left == right:
                result[i] = left
            else:
                stamp = Timestamp(vals[i])
                raise pytz.AmbiguousTimeError(stamp)
        elif left != NPY_NAT:
            result[i] = left
        elif right != NPY_NAT:
            result[i] = right
        else:
            stamp = Timestamp(vals[i])
            raise pytz.NonExistentTimeError(stamp)

    return result

cdef _ensure_int64(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_INT64:
            return arr
        else:
            return arr.astype(np.int64)
    else:
        return np.array(arr, dtype=np.int64)


cdef inline bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n):
    cdef Py_ssize_t pivot, left = 0, right = n

    # edge cases
    if val > data[n - 1]:
        return n

    if val < data[0]:
        return 0

    while left < right:
        pivot = left + (right - left) // 2

        if data[pivot] <= val:
            left = pivot + 1
        else:
            right = pivot

    return left


# Accessors
#----------------------------------------------------------------------

def build_field_sarray(ndarray[int64_t] dtindex):
    '''
    Datetime as int64 representation to a structured array of fields
    '''
    cdef:
        Py_ssize_t i, count = 0
        int isleap
        pandas_datetimestruct dts
        ndarray[int32_t] years, months, days, hours, minutes, seconds, mus

    count = len(dtindex)

    sa_dtype = [('Y', 'i4'), # year
                ('M', 'i4'), # month
                ('D', 'i4'), # day
                ('h', 'i4'), # hour
                ('m', 'i4'), # min
                ('s', 'i4'), # second
                ('u', 'i4')] # microsecond

    out = np.empty(count, dtype=sa_dtype)

    years = out['Y']
    months = out['M']
    days = out['D']
    hours = out['h']
    minutes = out['m']
    seconds = out['s']
    mus = out['u']

    for i in range(count):
        pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
        years[i] = dts.year
        months[i] = dts.month
        days[i] = dts.day
        hours[i] = dts.hour
        minutes[i] = dts.min
        seconds[i] = dts.sec
        mus[i] = dts.us

    return out

def get_time_micros(ndarray[int64_t] dtindex):
    '''
    Datetime as int64 representation to a structured array of fields
    '''
    cdef:
        Py_ssize_t i, n = len(dtindex)
        pandas_datetimestruct dts
        ndarray[int64_t] micros

    micros = np.empty(n, dtype=np.int64)

    for i in range(n):
        pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
        micros[i] = 1000000LL * (dts.hour * 60 * 60 +
                                 60 * dts.min + dts.sec) + dts.us

    return micros

@cython.wraparound(False)
def get_date_field(ndarray[int64_t] dtindex, object field):
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
        pandas_datetimestruct dts

    _month_offset = np.array(
        [[ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 ],
         [ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 ]],
         dtype=np.int32 )

    count = len(dtindex)
    out = np.empty(count, dtype='i4')

    if field == 'Y':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.year
        return out

    elif field == 'M':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.month
        return out

    elif field == 'D':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.day
        return out

    elif field == 'h':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.hour
        return out

    elif field == 'm':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.min
        return out

    elif field == 's':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.sec
        return out

    elif field == 'us':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.us
        return out
    elif field == 'ns':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.ps / 1000
        return out
    elif field == 'doy':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            isleap = is_leapyear(dts.year)
            out[i] = _month_offset[isleap, dts.month-1] + dts.day
        return out

    elif field == 'dow':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            ts = convert_to_tsobject(dtindex[i], None)
            out[i] = ts_dayofweek(ts)
        return out

    elif field == 'woy':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            isleap = is_leapyear(dts.year)
            out[i] = _month_offset[isleap, dts.month - 1] + dts.day
            out[i] = ((out[i] - 1) / 7) + 1
        return out

    elif field == 'q':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = -1; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            out[i] = dts.month
            out[i] = ((out[i] - 1) / 3) + 1
        return out

    raise ValueError("Field %s not supported" % field)


cdef inline int m8_weekday(int64_t val):
    ts = convert_to_tsobject(val, None)
    return ts_dayofweek(ts)

cdef int64_t DAY_NS = 86400000000000LL


def date_normalize(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts
        _TSObject tso
        ndarray[int64_t] result = np.empty(n, dtype=np.int64)

    if tz is not None:
        tso = _TSObject()
        if isinstance(tz, basestring):
            tz = pytz.timezone(tz)
        result = _normalize_local(stamps, tz)
    else:
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            result[i] = _normalized_stamp(&dts)

    return result

cdef _normalize_local(ndarray[int64_t] stamps, object tz):
    cdef:
        Py_ssize_t n = len(stamps)
        ndarray[int64_t] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans, deltas, pos
        pandas_datetimestruct dts

    if _is_utc(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            result[i] = _normalized_stamp(&dts)
    elif _is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns,
                                              &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
            pandas_datetime_to_datetimestruct(stamps[i] + delta,
                                              PANDAS_FR_ns, &dts)
            result[i] = _normalized_stamp(&dts)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans = _get_transitions(tz)
        deltas = _get_deltas(tz)
        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if not hasattr(tz, '_transition_info'):
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[0],
                                                  PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)
        else:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos[i]],
                                                  PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)

    return result

cdef inline int64_t _normalized_stamp(pandas_datetimestruct *dts):
    dts.hour = 0
    dts.min = 0
    dts.sec = 0
    dts.us = 0
    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, dts)


cdef inline void m8_populate_tsobject(int64_t stamp, _TSObject tso, object tz):
    tso.value = stamp
    pandas_datetime_to_datetimestruct(tso.value, PANDAS_FR_ns, &tso.dts)

    if tz is not None:
        _localize_tso(tso, tz)


def dates_normalized(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts

    if tz is None or _is_utc(tz):
        for i in range(n):
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False
    elif _is_tzlocal(tz):
        for i in range(n):
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            if (dts.min + dts.sec + dts.us) > 0:
                return False
            dt = datetime(dts.year, dts.month, dts.day, dts.hour, dts.min,
                          dts.sec, dts.us, tz)
            dt = dt + tz.utcoffset(dt)
            if dt.hour > 0:
                return False
    else:
        trans = _get_transitions(tz)
        deltas = _get_deltas(tz)
        for i in range(n):
            # Adjust datetime64 timestamp, recompute datetimestruct
            pos = trans.searchsorted(stamps[i]) - 1
            inf = tz._transition_info[pos]

            pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos],
                                              PANDAS_FR_ns, &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False

    return True

# Some general helper functions
#----------------------------------------------------------------------

def isleapyear(int64_t year):
    return is_leapyear(year)

def monthrange(int64_t year, int64_t month):
    cdef:
        int64_t days
        int64_t day_of_week

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = days_per_month_table[is_leapyear(year)][month-1]

    return (dayofweek(year, month, 1), days)

cdef inline int64_t ts_dayofweek(_TSObject ts):
    return dayofweek(ts.dts.year, ts.dts.month, ts.dts.day)


cpdef normalize_date(object dt):
    '''
    Normalize datetime.datetime value to midnight. Returns datetime.date as a
    datetime.datetime at midnight

    Returns
    -------
    normalized : datetime.datetime or Timestamp
    '''
    if PyDateTime_Check(dt):
        return dt.replace(hour=0, minute=0, second=0, microsecond=0)
    elif PyDate_Check(dt):
        return datetime(dt.year, dt.month, dt.day)
    else:
        raise TypeError('Unrecognized type: %s' % type(dt))

cdef ndarray[int64_t] localize_dt64arr_to_period(ndarray[int64_t] stamps,
                                                 int freq, object tz):
    cdef:
        Py_ssize_t n = len(stamps)
        ndarray[int64_t] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans, deltas, pos
        pandas_datetimestruct dts

    if not have_pytz:
        raise Exception('Could not find pytz module')

    if _is_utc(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            result[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                           dts.hour, dts.min, dts.sec, freq)

    elif _is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns,
                                              &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
            pandas_datetime_to_datetimestruct(stamps[i] + delta,
                                              PANDAS_FR_ns, &dts)
            result[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                           dts.hour, dts.min, dts.sec, freq)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans = _get_transitions(tz)
        deltas = _get_deltas(tz)
        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if not hasattr(tz, '_transition_info'):
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[0],
                                                  PANDAS_FR_ns, &dts)
                result[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                               dts.hour, dts.min, dts.sec, freq)
        else:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos[i]],
                                                  PANDAS_FR_ns, &dts)
                result[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                               dts.hour, dts.min, dts.sec, freq)

    return result


cdef extern from "period.h":
    ctypedef struct date_info:
        int64_t absdate
        double abstime
        double second
        int minute
        int hour
        int day
        int month
        int quarter
        int year
        int day_of_week
        int day_of_year
        int calendar

    ctypedef struct asfreq_info:
        int from_week_end
        int to_week_end

        int from_a_year_end
        int to_a_year_end

        int from_q_year_end
        int to_q_year_end

    ctypedef int64_t (*freq_conv_func)(int64_t, char, asfreq_info*)

    int64_t asfreq(int64_t dtordinal, int freq1, int freq2, char relation) except INT32_MIN
    freq_conv_func get_asfreq_func(int fromFreq, int toFreq)
    void get_asfreq_info(int fromFreq, int toFreq, asfreq_info *af_info)

    int64_t get_period_ordinal(int year, int month, int day,
                          int hour, int minute, int second,
                          int freq) except INT32_MIN

    int64_t get_python_ordinal(int64_t period_ordinal, int freq) except INT32_MIN

    int get_date_info(int64_t ordinal, int freq, date_info *dinfo) except INT32_MIN
    double getAbsTime(int, int64_t, int64_t)

    int pyear(int64_t ordinal, int freq) except INT32_MIN
    int pqyear(int64_t ordinal, int freq) except INT32_MIN
    int pquarter(int64_t ordinal, int freq) except INT32_MIN
    int pmonth(int64_t ordinal, int freq) except INT32_MIN
    int pday(int64_t ordinal, int freq) except INT32_MIN
    int pweekday(int64_t ordinal, int freq) except INT32_MIN
    int pday_of_week(int64_t ordinal, int freq) except INT32_MIN
    int pday_of_year(int64_t ordinal, int freq) except INT32_MIN
    int pweek(int64_t ordinal, int freq) except INT32_MIN
    int phour(int64_t ordinal, int freq) except INT32_MIN
    int pminute(int64_t ordinal, int freq) except INT32_MIN
    int psecond(int64_t ordinal, int freq) except INT32_MIN
    char *c_strftime(date_info *dinfo, char *fmt)
    int get_yq(int64_t ordinal, int freq, int *quarter, int *year)

# Period logic
#----------------------------------------------------------------------

cdef inline int64_t apply_mult(int64_t period_ord, int64_t mult):
    """
    Get freq+multiple ordinal value from corresponding freq-only ordinal value.
    For example, 5min ordinal will be 1/5th the 1min ordinal (rounding down to
    integer).
    """
    if mult == 1:
        return period_ord

    return (period_ord - 1) // mult

cdef inline int64_t remove_mult(int64_t period_ord_w_mult, int64_t mult):
    """
    Get freq-only ordinal value from corresponding freq+multiple ordinal.
    """
    if mult == 1:
        return period_ord_w_mult

    return period_ord_w_mult * mult + 1;

def dt64arr_to_periodarr(ndarray[int64_t] dtarr, int freq, tz=None):
    """
    Convert array of datetime64 values (passed in as 'i8' dtype) to a set of
    periods corresponding to desired frequency, per period convention.
    """
    cdef:
        ndarray[int64_t] out
        Py_ssize_t i, l
        pandas_datetimestruct dts

    l = len(dtarr)

    out = np.empty(l, dtype='i8')

    if tz is None:
        for i in range(l):
            pandas_datetime_to_datetimestruct(dtarr[i], PANDAS_FR_ns, &dts)
            out[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                        dts.hour, dts.min, dts.sec, freq)
    else:
        out = localize_dt64arr_to_period(dtarr, freq, tz)
    return out

def periodarr_to_dt64arr(ndarray[int64_t] periodarr, int freq):
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
        out[i] = period_ordinal_to_dt64(periodarr[i], freq)

    return out

cdef char START = 'S'
cdef char END = 'E'

cpdef int64_t period_asfreq(int64_t period_ordinal, int freq1, int freq2,
                            bint end):
    """
    Convert period ordinal from one frequency to another, and if upsampling,
    choose to use start ('S') or end ('E') of period.
    """
    cdef:
        int64_t retval

    if end:
        retval = asfreq(period_ordinal, freq1, freq2, END)
    else:
        retval = asfreq(period_ordinal, freq1, freq2, START)

    if retval == INT32_MIN:
        raise ValueError('Frequency conversion failed')

    return retval

def period_asfreq_arr(ndarray[int64_t] arr, int freq1, int freq2, bint end):
    """
    Convert int64-array of period ordinals from one frequency to another, and
    if upsampling, choose to use start ('S') or end ('E') of period.
    """
    cdef:
        ndarray[int64_t] result
        Py_ssize_t i, n
        freq_conv_func func
        asfreq_info finfo
        int64_t val, ordinal
        char relation

    n = len(arr)
    result = np.empty(n, dtype=np.int64)

    func = get_asfreq_func(freq1, freq2)
    get_asfreq_info(freq1, freq2, &finfo)

    if end:
        relation = END
    else:
        relation = START

    for i in range(n):
        val = func(arr[i], relation, &finfo)
        if val == INT32_MIN:
            raise ValueError("Unable to convert to desired frequency.")
        result[i] = val

    return result

def period_ordinal(int y, int m, int d, int h, int min, int s, int freq):
    cdef:
        int64_t ordinal

    return get_period_ordinal(y, m, d, h, min, s, freq)


cpdef int64_t period_ordinal_to_dt64(int64_t ordinal, int freq):
    cdef:
        pandas_datetimestruct dts
        date_info dinfo

    get_date_info(ordinal, freq, &dinfo)

    dts.year = dinfo.year
    dts.month = dinfo.month
    dts.day = dinfo.day
    dts.hour = dinfo.hour
    dts.min = dinfo.minute
    dts.sec = int(dinfo.second)
    dts.us = dts.ps = 0

    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)

def period_format(int64_t value, int freq, object fmt=None):
    cdef:
        int freq_group

    if fmt is None:
        freq_group = (freq // 1000) * 1000
        if freq_group == 1000: # FR_ANN
            fmt = b'%Y'
        elif freq_group == 2000: # FR_QTR
            fmt = b'%FQ%q'
        elif freq_group == 3000: # FR_MTH
            fmt = b'%Y-%m'
        elif freq_group == 4000: # WK
            left = period_asfreq(value, freq, 6000, 0)
            right = period_asfreq(value, freq, 6000, 1)
            return '%s/%s' % (period_format(left, 6000),
                              period_format(right, 6000))
        elif (freq_group == 5000 # BUS
              or freq_group == 6000): # DAY
            fmt = b'%Y-%m-%d'
        elif freq_group == 7000: # HR
            fmt = b'%Y-%m-%d %H:00'
        elif freq_group == 8000: # MIN
            fmt = b'%Y-%m-%d %H:%M'
        elif freq_group == 9000: # SEC
            fmt = b'%Y-%m-%d %H:%M:%S'
        else:
            raise ValueError('Unknown freq: %d' % freq)

    return _period_strftime(value, freq, fmt)


cdef list extra_fmts = [(b"%q", b"^`AB`^"),
                        (b"%f", b"^`CD`^"),
                        (b"%F", b"^`EF`^")]

cdef list str_extra_fmts = ["^`AB`^", "^`CD`^", "^`EF`^"]

cdef _period_strftime(int64_t value, int freq, object fmt):
    cdef:
        Py_ssize_t i
        date_info dinfo
        char *formatted
        object pat, repl, result
        list found_pat = [False] * len(extra_fmts)
        int year, quarter

    if PyUnicode_Check(fmt):
        fmt = fmt.encode('utf-8')

    get_date_info(value, freq, &dinfo)
    for i in range(len(extra_fmts)):
        pat = extra_fmts[i][0]
        repl = extra_fmts[i][1]
        if pat in fmt:
            fmt = fmt.replace(pat, repl)
            found_pat[i] = True

    formatted = c_strftime(&dinfo, <char*> fmt)

    result = util.char_to_string(formatted)
    free(formatted)

    for i in range(len(extra_fmts)):
        if found_pat[i]:
            if get_yq(value, freq, &quarter, &year) < 0:
                raise ValueError('Unable to get quarter and year')

            if i == 0:
                repl = '%d' % quarter
            elif i == 1:  # %f, 2-digit year
                repl = '%.2d' % (year % 100)
            elif i == 2:
                repl = '%d' % year

            result = result.replace(str_extra_fmts[i], repl)

    # Py3?
    if not PyString_Check(result):
        result = str(result)

    return result

# period accessors

ctypedef int (*accessor)(int64_t ordinal, int freq) except INT32_MIN

def get_period_field(int code, int64_t value, int freq):
    cdef accessor f = _get_accessor_func(code)
    return f(value, freq)

def get_period_field_arr(int code, ndarray[int64_t] arr, int freq):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] out
        accessor f

    f = _get_accessor_func(code)

    sz = len(arr)
    out = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        out[i] = f(arr[i], freq)

    return out



cdef accessor _get_accessor_func(int code):
    if code == 0:
        return &pyear
    elif code == 1:
        return &pqyear
    elif code == 2:
        return &pquarter
    elif code == 3:
        return &pmonth
    elif code == 4:
        return &pday
    elif code == 5:
        return &phour
    elif code == 6:
        return &pminute
    elif code == 7:
        return &psecond
    elif code == 8:
        return &pweek
    elif code == 9:
        return &pday_of_year
    elif code == 10:
        return &pweekday
    else:
        raise ValueError('Unrecognized code: %s' % code)


def extract_ordinals(ndarray[object] values, freq):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] ordinals = np.empty(n, dtype=np.int64)
        object p

    for i in range(n):
        p = values[i]
        ordinals[i] = p.ordinal
        if p.freq != freq:
            raise ValueError("%s is wrong freq" % p)

    return ordinals

cpdef resolution(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts
        int reso = D_RESO, curr_reso

    if tz is not None:
        if isinstance(tz, basestring):
            tz = pytz.timezone(tz)
        return _reso_local(stamps, tz)
    else:
        for i in range(n):
            if stamps[i] == NPY_NAT:
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
        return reso

US_RESO = 0
S_RESO = 1
T_RESO = 2
H_RESO = 3
D_RESO = 4

cdef inline int _reso_stamp(pandas_datetimestruct *dts):
    if dts.us != 0:
        return US_RESO
    elif dts.sec != 0:
        return S_RESO
    elif dts.min != 0:
        return T_RESO
    elif dts.hour != 0:
        return H_RESO
    return D_RESO

cdef _reso_local(ndarray[int64_t] stamps, object tz):
    cdef:
        Py_ssize_t n = len(stamps)
        int reso = D_RESO, curr_reso
        ndarray[int64_t] trans, deltas, pos
        pandas_datetimestruct dts

    if _is_utc(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
    elif _is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns,
                                              &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
            pandas_datetime_to_datetimestruct(stamps[i] + delta,
                                              PANDAS_FR_ns, &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans = _get_transitions(tz)
        deltas = _get_deltas(tz)
        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if not hasattr(tz, '_transition_info'):
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[0],
                                                  PANDAS_FR_ns, &dts)
                curr_reso = _reso_stamp(&dts)
                if curr_reso < reso:
                    reso = curr_reso
        else:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos[i]],
                                                  PANDAS_FR_ns, &dts)
                curr_reso = _reso_stamp(&dts)
                if curr_reso < reso:
                    reso = curr_reso

    return reso
