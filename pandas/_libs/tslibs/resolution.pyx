# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython
from cython cimport Py_ssize_t

import numpy as np
from numpy cimport ndarray, int64_t, int32_t

from util cimport is_string_object, get_nat

from np_datetime cimport npy_datetimestruct, dt64_to_dtstruct
from frequencies cimport get_freq_code
from timezones cimport is_utc, is_tzlocal, maybe_get_tz, get_dst_info
from conversion cimport tz_convert_utc_to_tzlocal
from ccalendar cimport get_days_in_month

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

# ----------------------------------------------------------------------

cpdef resolution(int64_t[:] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        npy_datetimestruct dts
        int reso = RESO_DAY, curr_reso

    if tz is not None:
        tz = maybe_get_tz(tz)
    return _reso_local(stamps, tz)


cdef _reso_local(int64_t[:] stamps, object tz):
    cdef:
        Py_ssize_t i, n = len(stamps)
        int reso = RESO_DAY, curr_reso
        ndarray[int64_t] trans
        int64_t[:] deltas
        Py_ssize_t[:] pos
        npy_datetimestruct dts
        int64_t local_val, delta

    if is_utc(tz) or tz is None:
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
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
            dt64_to_dtstruct(local_val, &dts)
            curr_reso = _reso_stamp(&dts)
            if curr_reso < reso:
                reso = curr_reso
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        if typ not in ['pytz', 'dateutil']:
            # static/fixed; in this case we know that len(delta) == 1
            delta = deltas[0]
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    continue
                dt64_to_dtstruct(stamps[i] + delta, &dts)
                curr_reso = _reso_stamp(&dts)
                if curr_reso < reso:
                    reso = curr_reso
        else:
            pos = trans.searchsorted(stamps, side='right') - 1
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    continue
                dt64_to_dtstruct(stamps[i] + deltas[pos[i]], &dts)
                curr_reso = _reso_stamp(&dts)
                if curr_reso < reso:
                    reso = curr_reso

    return reso


cdef inline int _reso_stamp(npy_datetimestruct *dts):
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

def month_position_check(fields, weekdays):
    cdef:
        int32_t daysinmonth, y, m, d
        bint calendar_end = True
        bint business_end = True
        bint calendar_start = True
        bint business_start = True
        bint cal
        int32_t[:] years
        int32_t[:] months
        int32_t[:] days

    years = fields['Y']
    months = fields['M']
    days = fields['D']

    for y, m, d, wd in zip(years, months, days, weekdays):
        if calendar_start:
            calendar_start &= d == 1
        if business_start:
            business_start &= d == 1 or (d <= 3 and wd == 0)

        if calendar_end or business_end:
            daysinmonth = get_days_in_month(y, m)
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
