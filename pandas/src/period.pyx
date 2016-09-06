from datetime import datetime, date, timedelta
import operator

from cpython cimport (
    PyObject_RichCompareBool,
    Py_EQ, Py_NE,
)

from numpy cimport (int8_t, int32_t, int64_t, import_array, ndarray,
                    NPY_INT64, NPY_DATETIME, NPY_TIMEDELTA)
import numpy as np

cdef extern from "datetime_helper.h":
    double total_seconds(object)

from libc.stdlib cimport free

from pandas import compat

from pandas.tseries import offsets
from pandas.tseries.tools import parse_time_string

cimport cython
from datetime cimport *
cimport util
cimport lib
from lib cimport is_null_datetimelike, is_period
import lib
from pandas import tslib
from tslib import Timedelta, Timestamp, iNaT, NaT
from tslib import have_pytz, _get_utcoffset
from tslib cimport (
    maybe_get_tz,
    _is_utc,
    _is_tzlocal,
    _get_dst_info,
    _nat_scalar_rules,
)

from pandas.tseries import frequencies

from sys import version_info

cdef bint PY2 = version_info[0] == 2

cdef int64_t NPY_NAT = util.get_nat()

cdef int US_RESO = frequencies.US_RESO
cdef int MS_RESO = frequencies.MS_RESO
cdef int S_RESO = frequencies.S_RESO
cdef int T_RESO = frequencies.T_RESO
cdef int H_RESO = frequencies.H_RESO
cdef int D_RESO = frequencies.D_RESO

cdef extern from "period_helper.h":
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

    void initialize_daytime_conversion_factor_matrix()
    int64_t asfreq(int64_t dtordinal, int freq1, int freq2,
                   char relation) except INT32_MIN
    freq_conv_func get_asfreq_func(int fromFreq, int toFreq)
    void get_asfreq_info(int fromFreq, int toFreq, asfreq_info *af_info)

    int64_t get_period_ordinal(int year, int month, int day,
                               int hour, int minute, int second,
                               int microseconds, int picoseconds,
                               int freq) nogil except INT32_MIN

    int64_t get_python_ordinal(int64_t period_ordinal,
                               int freq) except INT32_MIN

    int get_date_info(int64_t ordinal, int freq,
                      date_info *dinfo) nogil except INT32_MIN
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
    int pdays_in_month(int64_t ordinal, int freq) except INT32_MIN
    char *c_strftime(date_info *dinfo, char *fmt)
    int get_yq(int64_t ordinal, int freq, int *quarter, int *year)

initialize_daytime_conversion_factor_matrix()

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


@cython.wraparound(False)
@cython.boundscheck(False)
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
        with nogil:
            for i in range(l):
                if dtarr[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(dtarr[i], PANDAS_FR_ns, &dts)
                out[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                            dts.hour, dts.min, dts.sec,
                                            dts.us, dts.ps, freq)
    else:
        out = localize_dt64arr_to_period(dtarr, freq, tz)
    return out


@cython.wraparound(False)
@cython.boundscheck(False)
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

    with nogil:
        for i in range(l):
            if periodarr[i] == NPY_NAT:
                out[i] = NPY_NAT
                continue
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

    if period_ordinal == iNaT:
        return iNaT

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

    mask = arr == iNaT
    if mask.any():      # NaT process
        for i in range(n):
            val = arr[i]
            if val != iNaT:
                val = func(val, relation, &finfo)
                if val == INT32_MIN:
                    raise ValueError("Unable to convert to desired frequency.")
            result[i] = val
    else:
        for i in range(n):
            val = func(arr[i], relation, &finfo)
            if val == INT32_MIN:
                raise ValueError("Unable to convert to desired frequency.")
            result[i] = val

    return result


def period_ordinal(int y, int m, int d, int h, int min,
                   int s, int us, int ps, int freq):
    cdef:
        int64_t ordinal

    return get_period_ordinal(y, m, d, h, min, s, us, ps, freq)


cpdef int64_t period_ordinal_to_dt64(int64_t ordinal, int freq) nogil:
    cdef:
        pandas_datetimestruct dts
        date_info dinfo
        float subsecond_fraction

    if ordinal == NPY_NAT:
        return NPY_NAT

    get_date_info(ordinal, freq, &dinfo)

    dts.year = dinfo.year
    dts.month = dinfo.month
    dts.day = dinfo.day
    dts.hour = dinfo.hour
    dts.min = dinfo.minute
    dts.sec = int(dinfo.second)
    subsecond_fraction = dinfo.second - dts.sec
    dts.us = int((subsecond_fraction) * 1e6)
    dts.ps = int(((subsecond_fraction) * 1e6 - dts.us) * 1e6)

    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)


def period_format(int64_t value, int freq, object fmt=None):
    cdef:
        int freq_group

    if value == iNaT:
        return repr(NaT)

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
        elif freq_group == 10000: # MILLISEC
            fmt = b'%Y-%m-%d %H:%M:%S.%l'
        elif freq_group == 11000: # MICROSEC
            fmt = b'%Y-%m-%d %H:%M:%S.%u'
        elif freq_group == 12000: # NANOSEC
            fmt = b'%Y-%m-%d %H:%M:%S.%n'
        else:
            raise ValueError('Unknown freq: %d' % freq)

    return _period_strftime(value, freq, fmt)


cdef list extra_fmts = [(b"%q", b"^`AB`^"),
                        (b"%f", b"^`CD`^"),
                        (b"%F", b"^`EF`^"),
                        (b"%l", b"^`GH`^"),
                        (b"%u", b"^`IJ`^"),
                        (b"%n", b"^`KL`^")]

cdef list str_extra_fmts = ["^`AB`^", "^`CD`^", "^`EF`^",
                            "^`GH`^", "^`IJ`^", "^`KL`^"]

cdef object _period_strftime(int64_t value, int freq, object fmt):
    import sys

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
            elif i == 3:
                repl = '%03d' % (value % 1000)
            elif i == 4:
                repl = '%06d' % (value % 1000000)
            elif i == 5:
                repl = '%09d' % (value % 1000000000)

            result = result.replace(str_extra_fmts[i], repl)

    if PY2:
        result = result.decode('utf-8', 'ignore')

    return result

# period accessors

ctypedef int (*accessor)(int64_t ordinal, int freq) except INT32_MIN


def get_period_field(int code, int64_t value, int freq):
    cdef accessor f = _get_accessor_func(code)
    if f is NULL:
        raise ValueError('Unrecognized period code: %d' % code)
    if value == iNaT:
        return np.nan
    return f(value, freq)


def get_period_field_arr(int code, ndarray[int64_t] arr, int freq):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] out
        accessor f

    f = _get_accessor_func(code)
    if f is NULL:
        raise ValueError('Unrecognized period code: %d' % code)

    sz = len(arr)
    out = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        if arr[i] == iNaT:
            out[i] = -1
            continue
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
    elif code == 11:
        return &pdays_in_month
    return NULL


def extract_ordinals(ndarray[object] values, freq):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] ordinals = np.empty(n, dtype=np.int64)
        object p

    freqstr = Period._maybe_convert_freq(freq).freqstr

    for i in range(n):
        p = values[i]

        if is_null_datetimelike(p):
            ordinals[i] = tslib.iNaT
        else:
            try:
                ordinals[i] = p.ordinal

                if p.freqstr != freqstr:
                    msg = _DIFFERENT_FREQ_INDEX.format(freqstr, p.freqstr)
                    raise IncompatibleFrequency(msg)

            except AttributeError:
                p = Period(p, freq=freq)
                if p is tslib.NaT:
                    # input may contain NaT-like string
                    ordinals[i] = tslib.iNaT
                else:
                    ordinals[i] = p.ordinal

    return ordinals


def extract_freq(ndarray[object] values):
    cdef:
        Py_ssize_t i, n = len(values)
        object p

    for i in range(n):
        p = values[i]

        try:
            # now Timestamp / NaT has freq attr
            if is_period(p):
                return p.freq
        except AttributeError:
            pass

    raise ValueError('freq not specified and cannot be inferred')


cpdef resolution(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts
        int reso = D_RESO, curr_reso

    if tz is not None:
        tz = maybe_get_tz(tz)
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


cdef inline int _reso_stamp(pandas_datetimestruct *dts):
    if dts.us != 0:
        if dts.us % 1000 == 0:
            return MS_RESO
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
        trans, deltas, typ = _get_dst_info(tz)

        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if typ not in ['pytz', 'dateutil']:
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


# period helpers

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
                                           dts.hour, dts.min, dts.sec,
                                           dts.us, dts.ps, freq)

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
                                           dts.hour, dts.min, dts.sec,
                                           dts.us, dts.ps, freq)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = _get_dst_info(tz)

        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if typ not in ['pytz', 'dateutil']:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[0],
                                                  PANDAS_FR_ns, &dts)
                result[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                               dts.hour, dts.min, dts.sec,
                                               dts.us, dts.ps, freq)
        else:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos[i]],
                                                  PANDAS_FR_ns, &dts)
                result[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                               dts.hour, dts.min, dts.sec,
                                               dts.us, dts.ps, freq)

    return result


_DIFFERENT_FREQ = "Input has different freq={1} from Period(freq={0})"
_DIFFERENT_FREQ_INDEX = ("Input has different freq={1} "
                         "from PeriodIndex(freq={0})")


class IncompatibleFrequency(ValueError):
    pass


cdef class _Period(object):

    cdef public:
        int64_t ordinal
        object freq

    _comparables = ['name', 'freqstr']
    _typ = 'period'

    @classmethod
    def _maybe_convert_freq(cls, object freq):

        if isinstance(freq, (int, tuple)):
            code, stride = frequencies.get_freq_code(freq)
            freq = frequencies._get_freq_str(code, stride)

        freq = frequencies.to_offset(freq)

        if freq.n <= 0:
            raise ValueError('Frequency must be positive, because it'
                             ' represents span: {0}'.format(freq.freqstr))

        return freq

    @classmethod
    def _from_ordinal(cls, ordinal, freq):
        """
        Fast creation from an ordinal and freq that are already validated!
        """
        if ordinal == tslib.iNaT:
            return tslib.NaT
        else:
            self = _Period.__new__(cls)
            self.ordinal = ordinal
            self.freq = cls._maybe_convert_freq(freq)
            return self

    def __richcmp__(self, other, op):
        if isinstance(other, Period):
            if other.freq != self.freq:
                msg = _DIFFERENT_FREQ.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)
            return PyObject_RichCompareBool(self.ordinal, other.ordinal, op)
        elif other is tslib.NaT:
            return _nat_scalar_rules[op]
        # index/series like
        elif hasattr(other, '_typ'):
            return NotImplemented
        else:
            if op == Py_EQ:
                return NotImplemented
            elif op == Py_NE:
                return NotImplemented
            raise TypeError('Cannot compare type %r with type %r' %
                            (type(self).__name__, type(other).__name__))

    def __hash__(self):
        return hash((self.ordinal, self.freqstr))

    def _add_delta(self, other):
        if isinstance(other, (timedelta, np.timedelta64,
                              offsets.Tick, Timedelta)):
            offset = frequencies.to_offset(self.freq.rule_code)
            if isinstance(offset, offsets.Tick):
                nanos = tslib._delta_to_nanoseconds(other)
                offset_nanos = tslib._delta_to_nanoseconds(offset)

                if nanos % offset_nanos == 0:
                    ordinal = self.ordinal + (nanos // offset_nanos)
                    return Period(ordinal=ordinal, freq=self.freq)
            msg = 'Input cannot be converted to Period(freq={0})'
            raise IncompatibleFrequency(msg.format(self.freqstr))
        elif isinstance(other, offsets.DateOffset):
            freqstr = other.rule_code
            base = frequencies.get_base_alias(freqstr)
            if base == self.freq.rule_code:
                ordinal = self.ordinal + other.n
                return Period(ordinal=ordinal, freq=self.freq)
            msg = _DIFFERENT_FREQ.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)
        else: # pragma no cover
            return NotImplemented

    def __add__(self, other):
        if isinstance(self, Period):
            if isinstance(other, (timedelta, np.timedelta64,
                                  offsets.Tick, offsets.DateOffset,
                                  Timedelta)):
                return self._add_delta(other)
            elif other is tslib.NaT:
                return tslib.NaT
            elif lib.is_integer(other):
                ordinal = self.ordinal + other * self.freq.n
                return Period(ordinal=ordinal, freq=self.freq)
            else:  # pragma: no cover
                return NotImplemented
        elif isinstance(other, Period):
            return other + self
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(self, Period):
            if isinstance(other, (timedelta, np.timedelta64,
                                  offsets.Tick, offsets.DateOffset,
                                  Timedelta)):
                neg_other = -other
                return self + neg_other
            elif lib.is_integer(other):
                ordinal = self.ordinal - other * self.freq.n
                return Period(ordinal=ordinal, freq=self.freq)
            elif isinstance(other, Period):
                if other.freq != self.freq:
                    msg = _DIFFERENT_FREQ.format(self.freqstr, other.freqstr)
                    raise IncompatibleFrequency(msg)
                return self.ordinal - other.ordinal
            elif getattr(other, '_typ', None) == 'periodindex':
                return -other.__sub__(self)
            else:  # pragma: no cover
                return NotImplemented
        elif isinstance(other, Period):
            if self is tslib.NaT:
                return tslib.NaT
            return NotImplemented
        else:
            return NotImplemented

    def asfreq(self, freq, how='E'):
        """
        Convert Period to desired frequency, either at the start or end of the
        interval

        Parameters
        ----------
        freq : string
        how : {'E', 'S', 'end', 'start'}, default 'end'
            Start or end of the timespan

        Returns
        -------
        resampled : Period
        """
        freq = self._maybe_convert_freq(freq)
        how = _validate_end_alias(how)
        base1, mult1 = frequencies.get_freq_code(self.freq)
        base2, mult2 = frequencies.get_freq_code(freq)

        # mult1 can't be negative or 0
        end = how == 'E'
        if end:
            ordinal = self.ordinal + mult1 - 1
        else:
            ordinal = self.ordinal
        ordinal = period_asfreq(ordinal, base1, base2, end)

        return Period(ordinal=ordinal, freq=freq)

    @property
    def start_time(self):
        return self.to_timestamp(how='S')

    @property
    def end_time(self):
        # freq.n can't be negative or 0
        # ordinal = (self + self.freq.n).start_time.value - 1
        ordinal = (self + 1).start_time.value - 1
        return Timestamp(ordinal)

    def to_timestamp(self, freq=None, how='start', tz=None):
        """
        Return the Timestamp representation of the Period at the target
        frequency at the specified end (how) of the Period

        Parameters
        ----------
        freq : string or DateOffset, default is 'D' if self.freq is week or
               longer and 'S' otherwise
            Target frequency
        how: str, default 'S' (start)
            'S', 'E'. Can be aliased as case insensitive
            'Start', 'Finish', 'Begin', 'End'

        Returns
        -------
        Timestamp
        """
        if freq is not None:
            freq = self._maybe_convert_freq(freq)
        how = _validate_end_alias(how)

        if freq is None:
            base, mult = frequencies.get_freq_code(self.freq)
            freq = frequencies.get_to_timestamp_base(base)

        base, mult = frequencies.get_freq_code(freq)
        val = self.asfreq(freq, how)

        dt64 = period_ordinal_to_dt64(val.ordinal, base)
        return Timestamp(dt64, tz=tz)

    cdef _field(self, alias):
        base, mult = frequencies.get_freq_code(self.freq)
        return get_period_field(alias, self.ordinal, base)

    property year:
        def __get__(self):
            return self._field(0)
    property month:
        def __get__(self):
            return self._field(3)
    property day:
        def __get__(self):
            return self._field(4)
    property hour:
        def __get__(self):
            return self._field(5)
    property minute:
        def __get__(self):
            return self._field(6)
    property second:
        def __get__(self):
            return self._field(7)
    property weekofyear:
        def __get__(self):
            return self._field(8)
    property week:
        def __get__(self):
            return self.weekofyear
    property dayofweek:
        def __get__(self):
            return self._field(10)
    property weekday:
        def __get__(self):
            return self.dayofweek
    property dayofyear:
        def __get__(self):
            return self._field(9)
    property quarter:
        def __get__(self):
            return self._field(2)
    property qyear:
        def __get__(self):
            return self._field(1)
    property days_in_month:
        def __get__(self):
            return self._field(11)
    property daysinmonth:
        def __get__(self):
            return self.days_in_month
    property is_leap_year:
        def __get__(self):
            return bool(is_leapyear(self._field(0)))

    @classmethod
    def now(cls, freq=None):
        return Period(datetime.now(), freq=freq)

    # HACK IT UP AND YOU BETTER FIX IT SOON
    def __str__(self):
        return self.__unicode__()

    @property
    def freqstr(self):
        return self.freq.freqstr

    def __repr__(self):
        base, mult = frequencies.get_freq_code(self.freq)
        formatted = period_format(self.ordinal, base)
        return "Period('%s', '%s')" % (formatted, self.freqstr)

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        base, mult = frequencies.get_freq_code(self.freq)
        formatted = period_format(self.ordinal, base)
        value = ("%s" % formatted)
        return value

    def __setstate__(self, state):
        self.freq=state[1]
        self.ordinal=state[2]

    def __reduce__(self):
        object_state = None, self.freq, self.ordinal
        return (Period, object_state)

    def strftime(self, fmt):
        """
        Returns the string representation of the :class:`Period`, depending
        on the selected :keyword:`format`. :keyword:`format` must be a string
        containing one or several directives.  The method recognizes the same
        directives as the :func:`time.strftime` function of the standard Python
        distribution, as well as the specific additional directives ``%f``,
        ``%F``, ``%q``. (formatting & docs originally from scikits.timeries)

        +-----------+--------------------------------+-------+
        | Directive | Meaning                        | Notes |
        +===========+================================+=======+
        | ``%a``    | Locale's abbreviated weekday   |       |
        |           | name.                          |       |
        +-----------+--------------------------------+-------+
        | ``%A``    | Locale's full weekday name.    |       |
        +-----------+--------------------------------+-------+
        | ``%b``    | Locale's abbreviated month     |       |
        |           | name.                          |       |
        +-----------+--------------------------------+-------+
        | ``%B``    | Locale's full month name.      |       |
        +-----------+--------------------------------+-------+
        | ``%c``    | Locale's appropriate date and  |       |
        |           | time representation.           |       |
        +-----------+--------------------------------+-------+
        | ``%d``    | Day of the month as a decimal  |       |
        |           | number [01,31].                |       |
        +-----------+--------------------------------+-------+
        | ``%f``    | 'Fiscal' year without a        | \(1)  |
        |           | century  as a decimal number   |       |
        |           | [00,99]                        |       |
        +-----------+--------------------------------+-------+
        | ``%F``    | 'Fiscal' year with a century   | \(2)  |
        |           | as a decimal number            |       |
        +-----------+--------------------------------+-------+
        | ``%H``    | Hour (24-hour clock) as a      |       |
        |           | decimal number [00,23].        |       |
        +-----------+--------------------------------+-------+
        | ``%I``    | Hour (12-hour clock) as a      |       |
        |           | decimal number [01,12].        |       |
        +-----------+--------------------------------+-------+
        | ``%j``    | Day of the year as a decimal   |       |
        |           | number [001,366].              |       |
        +-----------+--------------------------------+-------+
        | ``%m``    | Month as a decimal number      |       |
        |           | [01,12].                       |       |
        +-----------+--------------------------------+-------+
        | ``%M``    | Minute as a decimal number     |       |
        |           | [00,59].                       |       |
        +-----------+--------------------------------+-------+
        | ``%p``    | Locale's equivalent of either  | \(3)  |
        |           | AM or PM.                      |       |
        +-----------+--------------------------------+-------+
        | ``%q``    | Quarter as a decimal number    |       |
        |           | [01,04]                        |       |
        +-----------+--------------------------------+-------+
        | ``%S``    | Second as a decimal number     | \(4)  |
        |           | [00,61].                       |       |
        +-----------+--------------------------------+-------+
        | ``%U``    | Week number of the year        | \(5)  |
        |           | (Sunday as the first day of    |       |
        |           | the week) as a decimal number  |       |
        |           | [00,53].  All days in a new    |       |
        |           | year preceding the first       |       |
        |           | Sunday are considered to be in |       |
        |           | week 0.                        |       |
        +-----------+--------------------------------+-------+
        | ``%w``    | Weekday as a decimal number    |       |
        |           | [0(Sunday),6].                 |       |
        +-----------+--------------------------------+-------+
        | ``%W``    | Week number of the year        | \(5)  |
        |           | (Monday as the first day of    |       |
        |           | the week) as a decimal number  |       |
        |           | [00,53].  All days in a new    |       |
        |           | year preceding the first       |       |
        |           | Monday are considered to be in |       |
        |           | week 0.                        |       |
        +-----------+--------------------------------+-------+
        | ``%x``    | Locale's appropriate date      |       |
        |           | representation.                |       |
        +-----------+--------------------------------+-------+
        | ``%X``    | Locale's appropriate time      |       |
        |           | representation.                |       |
        +-----------+--------------------------------+-------+
        | ``%y``    | Year without century as a      |       |
        |           | decimal number [00,99].        |       |
        +-----------+--------------------------------+-------+
        | ``%Y``    | Year with century as a decimal |       |
        |           | number.                        |       |
        +-----------+--------------------------------+-------+
        | ``%Z``    | Time zone name (no characters  |       |
        |           | if no time zone exists).       |       |
        +-----------+--------------------------------+-------+
        | ``%%``    | A literal ``'%'`` character.   |       |
        +-----------+--------------------------------+-------+

        .. note::

            (1)
                The ``%f`` directive is the same as ``%y`` if the frequency is
                not quarterly.
                Otherwise, it corresponds to the 'fiscal' year, as defined by
                the :attr:`qyear` attribute.

            (2)
                The ``%F`` directive is the same as ``%Y`` if the frequency is
                not quarterly.
                Otherwise, it corresponds to the 'fiscal' year, as defined by
                the :attr:`qyear` attribute.

            (3)
                The ``%p`` directive only affects the output hour field
                if the ``%I`` directive is used to parse the hour.

            (4)
                The range really is ``0`` to ``61``; this accounts for leap
                seconds and the (very rare) double leap seconds.

            (5)
                The ``%U`` and ``%W`` directives are only used in calculations
                when the day of the week and the year are specified.

        .. rubric::  Examples

            >>> a = Period(freq='Q@JUL', year=2006, quarter=1)
            >>> a.strftime('%F-Q%q')
            '2006-Q1'
            >>> # Output the last month in the quarter of this date
            >>> a.strftime('%b-%Y')
            'Oct-2005'
            >>>
            >>> a = Period(freq='D', year=2001, month=1, day=1)
            >>> a.strftime('%d-%b-%Y')
            '01-Jan-2006'
            >>> a.strftime('%b. %d, %Y was a %A')
            'Jan. 01, 2001 was a Monday'
        """
        base, mult = frequencies.get_freq_code(self.freq)
        return period_format(self.ordinal, base, fmt)


class Period(_Period):
    """
    Represents an period of time

    Parameters
    ----------
    value : Period or compat.string_types, default None
        The time period represented (e.g., '4Q2005')
    freq : str, default None
        One of pandas period strings or corresponding objects
    year : int, default None
    month : int, default 1
    quarter : int, default None
    day : int, default 1
    hour : int, default 0
    minute : int, default 0
    second : int, default 0
    """

    def __new__(cls, value=None, freq=None, ordinal=None,
                year=None, month=None, quarter=None, day=None,
                hour=None, minute=None, second=None):
        # freq points to a tuple (base, mult);  base is one of the defined
        # periods such as A, Q, etc. Every five minutes would be, e.g.,
        # ('T', 5) but may be passed in as a string like '5T'

        # ordinal is the period offset from the gregorian proleptic epoch

        cdef _Period self

        if freq is not None:
            freq = cls._maybe_convert_freq(freq)

        if ordinal is not None and value is not None:
            raise ValueError(("Only value or ordinal but not both should be "
                              "given but not both"))
        elif ordinal is not None:
            if not lib.is_integer(ordinal):
                raise ValueError("Ordinal must be an integer")
            if freq is None:
                raise ValueError('Must supply freq for ordinal value')

        elif value is None:
            if (year is None and month is None and
                        quarter is None and day is None and
                        hour is None and minute is None and second is None):
                ordinal = tslib.iNaT
            else:
                if freq is None:
                    raise ValueError("If value is None, freq cannot be None")

                # set defaults
                month = 1 if month is None else month
                day = 1 if day is None else day
                hour = 0 if hour is None else hour
                minute = 0 if minute is None else minute
                second = 0 if second is None else second

                ordinal = _ordinal_from_fields(year, month, quarter, day,
                                               hour, minute, second, freq)

        elif isinstance(value, Period):
            other = value
            if freq is None or frequencies.get_freq_code(
                    freq) == frequencies.get_freq_code(other.freq):
                ordinal = other.ordinal
                freq = other.freq
            else:
                converted = other.asfreq(freq)
                ordinal = converted.ordinal

        elif is_null_datetimelike(value) or value in tslib._nat_strings:
            ordinal = tslib.iNaT

        elif isinstance(value, compat.string_types) or lib.is_integer(value):
            if lib.is_integer(value):
                value = str(value)
            value = value.upper()
            dt, _, reso = parse_time_string(value, freq)

            if freq is None:
                try:
                    freq = frequencies.Resolution.get_freq(reso)
                except KeyError:
                    raise ValueError(
                        "Invalid frequency or could not infer: %s" % reso)

        elif isinstance(value, datetime):
            dt = value
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        elif isinstance(value, np.datetime64):
            dt = Timestamp(value)
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        elif isinstance(value, date):
            dt = datetime(year=value.year, month=value.month, day=value.day)
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        else:
            msg = "Value must be Period, string, integer, or datetime"
            raise ValueError(msg)

        if ordinal is None:
            base, mult = frequencies.get_freq_code(freq)
            ordinal = get_period_ordinal(dt.year, dt.month, dt.day,
                                         dt.hour, dt.minute, dt.second,
                                         dt.microsecond, 0, base)

        return cls._from_ordinal(ordinal, freq)


def _ordinal_from_fields(year, month, quarter, day,
                         hour, minute, second, freq):
    base, mult = frequencies.get_freq_code(freq)
    if quarter is not None:
        year, month = _quarter_to_myear(year, quarter, freq)

    return get_period_ordinal(year, month, day, hour,
                              minute, second, 0, 0, base)


def _quarter_to_myear(year, quarter, freq):
    if quarter is not None:
        if quarter <= 0 or quarter > 4:
            raise ValueError('Quarter must be 1 <= q <= 4')

        mnum = frequencies._month_numbers[
            frequencies._get_rule_month(freq)] + 1
        month = (mnum + (quarter - 1) * 3) % 12 + 1
        if month > mnum:
            year -= 1

    return year, month


def _validate_end_alias(how):
    how_dict = {'S': 'S', 'E': 'E',
                'START': 'S', 'FINISH': 'E',
                'BEGIN': 'S', 'END': 'E'}
    how = how_dict.get(str(how).upper())
    if how not in set(['S', 'E']):
        raise ValueError('How must be one of S or E')
    return how
