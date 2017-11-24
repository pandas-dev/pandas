# -*- coding: utf-8 -*-
# cython: profile=False

from cython cimport Py_ssize_t

import numpy as np
cimport numpy as np
from numpy cimport ndarray, int64_t
np.import_array()

from util cimport is_string_object, get_nat

from pandas._libs.khash cimport (
    khiter_t,
    kh_destroy_int64, kh_put_int64,
    kh_init_int64, kh_int64_t,
    kh_resize_int64, kh_get_int64)

from cpython.datetime cimport datetime

from np_datetime cimport (pandas_datetimestruct,
                          dtstruct_to_dt64, dt64_to_dtstruct)
from frequencies cimport get_freq_code
from timezones cimport (
    is_utc, is_tzlocal,
    maybe_get_tz, get_dst_info, get_utcoffset)
from fields import build_field_sarray
from conversion import tz_convert

from pandas._libs.properties import cache_readonly
from pandas._libs.tslib import Timestamp

from pandas.core.algorithms import unique  # TODO: Avoid this non-cython import

# ----------------------------------------------------------------------
# Constants

cdef int64_t NPY_NAT = get_nat()

cdef int RESO_NS = 0
cdef int RESO_US = 1
cdef int RESO_MS = 2
cdef int RESO_SEC = 3
cdef int RESO_MIN = 4
cdef int RESO_HR = 5
cdef int RESO_DAY = 6

_ONE_MICRO = 1000L
_ONE_MILLI = _ONE_MICRO * 1000
_ONE_SECOND = _ONE_MILLI * 1000
_ONE_MINUTE = 60 * _ONE_SECOND
_ONE_HOUR = 60 * _ONE_MINUTE
_ONE_DAY = 24 * _ONE_HOUR

DAYS = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
_weekday_rule_aliases = {k: v for k, v in enumerate(DAYS)}

_MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
           'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
_MONTH_ALIASES = {(k + 1): v for k, v in enumerate(_MONTHS)}

# ----------------------------------------------------------------------

cpdef resolution(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts
        int reso = RESO_DAY, curr_reso

    if tz is not None:
        tz = maybe_get_tz(tz)
        return _reso_local(stamps, tz)
    else:
        for i in range(n):
            if stamps[i] == NPY_NAT:
                continue
            dt64_to_dtstruct(stamps[i], &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
        return reso


cdef _reso_local(ndarray[int64_t] stamps, object tz):
    cdef:
        Py_ssize_t n = len(stamps)
        int reso = RESO_DAY, curr_reso
        ndarray[int64_t] trans, deltas, pos
        pandas_datetimestruct dts

    if is_utc(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                continue
            dt64_to_dtstruct(stamps[i], &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
    elif is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                continue
            dt64_to_dtstruct(stamps[i], &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000
            dt64_to_dtstruct(stamps[i] + delta, &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if typ not in ['pytz', 'dateutil']:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    continue
                dt64_to_dtstruct(stamps[i] + deltas[0], &dts)
                curr_reso = _reso_stamp(&dts)
                if curr_reso < reso:
                    reso = curr_reso
        else:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    continue
                dt64_to_dtstruct(stamps[i] + deltas[pos[i]], &dts)
                curr_reso = _reso_stamp(&dts)
                if curr_reso < reso:
                    reso = curr_reso

    return reso


cdef inline int _reso_stamp(pandas_datetimestruct *dts):
    if dts.us != 0:
        if dts.us % 1000 == 0:
            return RESO_MS
        return RESO_US
    elif dts.sec != 0:
        return RESO_SEC
    elif dts.min != 0:
        return RESO_MIN
    elif dts.hour != 0:
        return RESO_HR
    return RESO_DAY


def get_freq_group(freq):
    """
    Return frequency code group of given frequency str or offset.

    Example
    -------
    >>> get_freq_group('W-MON')
    4000

    >>> get_freq_group('W-FRI')
    4000
    """
    if getattr(freq, '_typ', None) == 'dateoffset':
        freq = freq.rule_code

    if is_string_object(freq):
        base, mult = get_freq_code(freq)
        freq = base
    elif isinstance(freq, int):
        pass
    else:
        raise ValueError('input must be str, offset or int')
    return (freq // 1000) * 1000


class Resolution(object):

    # Note: cython won't allow us to reference the cdef versions at the
    # module level
    RESO_NS = 0
    RESO_US = 1
    RESO_MS = 2
    RESO_SEC = 3
    RESO_MIN = 4
    RESO_HR = 5
    RESO_DAY = 6

    _reso_str_map = {
        RESO_NS: 'nanosecond',
        RESO_US: 'microsecond',
        RESO_MS: 'millisecond',
        RESO_SEC: 'second',
        RESO_MIN: 'minute',
        RESO_HR: 'hour',
        RESO_DAY: 'day'}

    # factor to multiply a value by to convert it to the next finer grained
    # resolution
    _reso_mult_map = {
        RESO_NS: None,
        RESO_US: 1000,
        RESO_MS: 1000,
        RESO_SEC: 1000,
        RESO_MIN: 60,
        RESO_HR: 60,
        RESO_DAY: 24}

    _reso_str_bump_map = {
        'D': 'H',
        'H': 'T',
        'T': 'S',
        'S': 'L',
        'L': 'U',
        'U': 'N',
        'N': None}

    _str_reso_map = {v: k for k, v in _reso_str_map.items()}

    _reso_freq_map = {
        'year': 'A',
        'quarter': 'Q',
        'month': 'M',
        'day': 'D',
        'hour': 'H',
        'minute': 'T',
        'second': 'S',
        'millisecond': 'L',
        'microsecond': 'U',
        'nanosecond': 'N'}

    _freq_reso_map = {v: k for k, v in _reso_freq_map.items()}

    @classmethod
    def get_str(cls, reso):
        """
        Return resolution str against resolution code.

        Example
        -------
        >>> Resolution.get_str(Resolution.RESO_SEC)
        'second'
        """
        return cls._reso_str_map.get(reso, 'day')

    @classmethod
    def get_reso(cls, resostr):
        """
        Return resolution str against resolution code.

        Example
        -------
        >>> Resolution.get_reso('second')
        2

        >>> Resolution.get_reso('second') == Resolution.RESO_SEC
        True
        """
        return cls._str_reso_map.get(resostr, cls.RESO_DAY)

    @classmethod
    def get_freq_group(cls, resostr):
        """
        Return frequency str against resolution str.

        Example
        -------
        >>> f.Resolution.get_freq_group('day')
        4000
        """
        return get_freq_group(cls.get_freq(resostr))

    @classmethod
    def get_freq(cls, resostr):
        """
        Return frequency str against resolution str.

        Example
        -------
        >>> f.Resolution.get_freq('day')
        'D'
        """
        return cls._reso_freq_map[resostr]

    @classmethod
    def get_str_from_freq(cls, freq):
        """
        Return resolution str against frequency str.

        Example
        -------
        >>> Resolution.get_str_from_freq('H')
        'hour'
        """
        return cls._freq_reso_map.get(freq, 'day')

    @classmethod
    def get_reso_from_freq(cls, freq):
        """
        Return resolution code against frequency str.

        Example
        -------
        >>> Resolution.get_reso_from_freq('H')
        4

        >>> Resolution.get_reso_from_freq('H') == Resolution.RESO_HR
        True
        """
        return cls.get_reso(cls.get_str_from_freq(freq))

    @classmethod
    def get_stride_from_decimal(cls, value, freq):
        """
        Convert freq with decimal stride into a higher freq with integer stride

        Parameters
        ----------
        value : integer or float
        freq : string
            Frequency string

        Raises
        ------
        ValueError
            If the float cannot be converted to an integer at any resolution.

        Example
        -------
        >>> Resolution.get_stride_from_decimal(1.5, 'T')
        (90, 'S')

        >>> Resolution.get_stride_from_decimal(1.04, 'H')
        (3744, 'S')

        >>> Resolution.get_stride_from_decimal(1, 'D')
        (1, 'D')
        """
        if np.isclose(value % 1, 0):
            return int(value), freq
        else:
            start_reso = cls.get_reso_from_freq(freq)
            if start_reso == 0:
                raise ValueError("Could not convert to integer offset "
                                 "at any resolution")

            next_value = cls._reso_mult_map[start_reso] * value
            next_name = cls._reso_str_bump_map[freq]
            return cls.get_stride_from_decimal(next_value, next_name)


# ----------------------------------------------------------------------
# Frequency Inference


# TODO: this is non performiant logic here (and duplicative) and this
# simply should call unique_1d directly
# plus no reason to depend on khash directly
cdef unique_deltas(ndarray[int64_t] arr):
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


def _is_multiple(us, mult):
    return us % mult == 0


def _maybe_add_count(base, count):
    if count != 1:
        return '{count}{base}'.format(count=int(count), base=base)
    else:
        return base


class _FrequencyInferer(object):
    """
    Not sure if I can avoid the state machine here
    """

    def __init__(self, index, warn=True):
        self.index = index
        self.values = np.asarray(index).view('i8')

        # This moves the values, which are implicitly in UTC, to the
        # the timezone so they are in local time
        if hasattr(index, 'tz'):
            if index.tz is not None:
                self.values = tz_convert(self.values, 'UTC', index.tz)

        self.warn = warn

        if len(index) < 3:
            raise ValueError('Need at least 3 dates to infer frequency')

        self.is_monotonic = (self.index.is_monotonic_increasing or
                             self.index.is_monotonic_decreasing)

    @cache_readonly
    def deltas(self):
        return unique_deltas(self.values)

    @cache_readonly
    def deltas_asi8(self):
        return unique_deltas(self.index.asi8)

    @cache_readonly
    def is_unique(self):
        return len(self.deltas) == 1

    @cache_readonly
    def is_unique_asi8(self):
        return len(self.deltas_asi8) == 1

    def get_freq(self):
        if not self.is_monotonic or not self.index.is_unique:
            return None

        delta = self.deltas[0]
        if _is_multiple(delta, _ONE_DAY):
            return self._infer_daily_rule()
        else:
            # Business hourly, maybe. 17: one day / 65: one weekend
            if self.hour_deltas in ([1, 17], [1, 65], [1, 17, 65]):
                return 'BH'
            # Possibly intraday frequency.  Here we use the
            # original .asi8 values as the modified values
            # will not work around DST transitions.  See #8772
            elif not self.is_unique_asi8:
                return None
            delta = self.deltas_asi8[0]
            if _is_multiple(delta, _ONE_HOUR):
                # Hours
                return _maybe_add_count('H', delta / _ONE_HOUR)
            elif _is_multiple(delta, _ONE_MINUTE):
                # Minutes
                return _maybe_add_count('T', delta / _ONE_MINUTE)
            elif _is_multiple(delta, _ONE_SECOND):
                # Seconds
                return _maybe_add_count('S', delta / _ONE_SECOND)
            elif _is_multiple(delta, _ONE_MILLI):
                # Milliseconds
                return _maybe_add_count('L', delta / _ONE_MILLI)
            elif _is_multiple(delta, _ONE_MICRO):
                # Microseconds
                return _maybe_add_count('U', delta / _ONE_MICRO)
            else:
                # Nanoseconds
                return _maybe_add_count('N', delta)

    @cache_readonly
    def day_deltas(self):
        return [x / _ONE_DAY for x in self.deltas]

    @cache_readonly
    def hour_deltas(self):
        return [x / _ONE_HOUR for x in self.deltas]

    @cache_readonly
    def fields(self):
        return build_field_sarray(self.values)

    @cache_readonly
    def rep_stamp(self):
        return Timestamp(self.values[0])

    def month_position_check(self):
        # TODO: cythonize this, very slow
        calendar_end = True
        business_end = True
        calendar_start = True
        business_start = True

        years = self.fields['Y']
        months = self.fields['M']
        days = self.fields['D']
        weekdays = self.index.dayofweek

        from calendar import monthrange
        for y, m, d, wd in zip(years, months, days, weekdays):

            if calendar_start:
                calendar_start &= d == 1
            if business_start:
                business_start &= d == 1 or (d <= 3 and wd == 0)

            if calendar_end or business_end:
                _, daysinmonth = monthrange(y, m)
                cal = d == daysinmonth
                if calendar_end:
                    calendar_end &= cal
                if business_end:
                    business_end &= cal or (daysinmonth - d < 3 and wd == 4)
            elif not calendar_start and not business_start:
                break

        if calendar_end:
            return 'ce'
        elif business_end:
            return 'be'
        elif calendar_start:
            return 'cs'
        elif business_start:
            return 'bs'
        else:
            return None

    @cache_readonly
    def mdiffs(self):
        nmonths = self.fields['Y'] * 12 + self.fields['M']
        return unique_deltas(nmonths.astype('i8'))

    @cache_readonly
    def ydiffs(self):
        return unique_deltas(self.fields['Y'].astype('i8'))

    def _infer_daily_rule(self):
        annual_rule = self._get_annual_rule()
        if annual_rule:
            nyears = self.ydiffs[0]
            month = _MONTH_ALIASES[self.rep_stamp.month]
            alias = '{prefix}-{month}'.format(prefix=annual_rule, month=month)
            return _maybe_add_count(alias, nyears)

        quarterly_rule = self._get_quarterly_rule()
        if quarterly_rule:
            nquarters = self.mdiffs[0] / 3
            mod_dict = {0: 12, 2: 11, 1: 10}
            month = _MONTH_ALIASES[mod_dict[self.rep_stamp.month % 3]]
            alias = '{prefix}-{month}'.format(prefix=quarterly_rule,
                                              month=month)
            return _maybe_add_count(alias, nquarters)

        monthly_rule = self._get_monthly_rule()
        if monthly_rule:
            return _maybe_add_count(monthly_rule, self.mdiffs[0])

        if self.is_unique:
            days = self.deltas[0] / _ONE_DAY
            if days % 7 == 0:
                # Weekly
                day = _weekday_rule_aliases[self.rep_stamp.weekday()]
                return _maybe_add_count('W-{day}'.format(day=day), days / 7)
            else:
                return _maybe_add_count('D', days)

        if self._is_business_daily():
            return 'B'

        wom_rule = self._get_wom_rule()
        if wom_rule:
            return wom_rule

    def _get_annual_rule(self):
        if len(self.ydiffs) > 1:
            return None

        if len(unique(self.fields['M'])) > 1:
            return None

        pos_check = self.month_position_check()
        return {'cs': 'AS', 'bs': 'BAS',
                'ce': 'A', 'be': 'BA'}.get(pos_check)

    def _get_quarterly_rule(self):
        if len(self.mdiffs) > 1:
            return None

        if not self.mdiffs[0] % 3 == 0:
            return None

        pos_check = self.month_position_check()
        return {'cs': 'QS', 'bs': 'BQS',
                'ce': 'Q', 'be': 'BQ'}.get(pos_check)

    def _get_monthly_rule(self):
        if len(self.mdiffs) > 1:
            return None
        pos_check = self.month_position_check()
        return {'cs': 'MS', 'bs': 'BMS',
                'ce': 'M', 'be': 'BM'}.get(pos_check)

    def _is_business_daily(self):
        # quick check: cannot be business daily
        if self.day_deltas != [1, 3]:
            return False

        # probably business daily, but need to confirm
        first_weekday = self.index[0].weekday()
        shifts = np.diff(self.index.asi8)
        shifts = np.floor_divide(shifts, _ONE_DAY)
        weekdays = np.mod(first_weekday + np.cumsum(shifts), 7)
        return np.all(((weekdays == 0) & (shifts == 3)) |
                      ((weekdays > 0) & (weekdays <= 4) & (shifts == 1)))

    def _get_wom_rule(self):
        #         wdiffs = unique(np.diff(self.index.week))
        # We also need -47, -49, -48 to catch index spanning year boundary
        #     if not lib.ismember(wdiffs, set([4, 5, -47, -49, -48])).all():
        #         return None

        weekdays = unique(self.index.weekday)
        if len(weekdays) > 1:
            return None

        week_of_months = unique((self.index.day - 1) // 7)
        # Only attempt to infer up to WOM-4. See #9425
        week_of_months = week_of_months[week_of_months < 4]
        if len(week_of_months) == 0 or len(week_of_months) > 1:
            return None

        # get which week
        week = week_of_months[0] + 1
        wd = _weekday_rule_aliases[weekdays[0]]

        return 'WOM-{week}{weekday}'.format(week=week, weekday=wd)


class _TimedeltaFrequencyInferer(_FrequencyInferer):

    def _infer_daily_rule(self):
        if self.is_unique:
            days = self.deltas[0] / _ONE_DAY
            if days % 7 == 0:
                # Weekly
                wd = _weekday_rule_aliases[self.rep_stamp.weekday()]
                alias = 'W-{weekday}'.format(weekday=wd)
                return _maybe_add_count(alias, days / 7)
            else:
                return _maybe_add_count('D', days)
