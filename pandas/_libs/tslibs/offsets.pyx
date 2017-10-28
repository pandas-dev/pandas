# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython

import time
from cpython.datetime cimport time as dt_time

import numpy as np
cimport numpy as np
np.import_array()


from util cimport is_string_object


from pandas._libs.tslib import pydt_to_i8, tz_convert_single

# ---------------------------------------------------------------------
# Constants

# Duplicated in tslib
_MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
           'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
_int_to_month = {(k + 1): v for k, v in enumerate(_MONTHS)}
_month_to_int = dict((v, k) for k, v in _int_to_month.items())


class WeekDay(object):
    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6


_int_to_weekday = {
    WeekDay.MON: 'MON',
    WeekDay.TUE: 'TUE',
    WeekDay.WED: 'WED',
    WeekDay.THU: 'THU',
    WeekDay.FRI: 'FRI',
    WeekDay.SAT: 'SAT',
    WeekDay.SUN: 'SUN'}

_weekday_to_int = {_int_to_weekday[key]: key for key in _int_to_weekday}


_offset_to_period_map = {
    'WEEKDAY': 'D',
    'EOM': 'M',
    'BM': 'M',
    'BQS': 'Q',
    'QS': 'Q',
    'BQ': 'Q',
    'BA': 'A',
    'AS': 'A',
    'BAS': 'A',
    'MS': 'M',
    'D': 'D',
    'C': 'C',
    'B': 'B',
    'T': 'T',
    'S': 'S',
    'L': 'L',
    'U': 'U',
    'N': 'N',
    'H': 'H',
    'Q': 'Q',
    'A': 'A',
    'W': 'W',
    'M': 'M',
    'Y': 'A',
    'BY': 'A',
    'YS': 'A',
    'BYS': 'A'}

need_suffix = ['QS', 'BQ', 'BQS', 'YS', 'AS', 'BY', 'BA', 'BYS', 'BAS']


for __prefix in need_suffix:
    for _m in _MONTHS:
        key = '%s-%s' % (__prefix, _m)
        _offset_to_period_map[key] = _offset_to_period_map[__prefix]

for __prefix in ['A', 'Q']:
    for _m in _MONTHS:
        _alias = '%s-%s' % (__prefix, _m)
        _offset_to_period_map[_alias] = _alias

_days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
for _d in _days:
    _offset_to_period_map['W-%s' % _d] = 'W-%s' % _d


# ---------------------------------------------------------------------
# Misc Helpers

def as_datetime(obj):
    f = getattr(obj, 'to_pydatetime', None)
    if f is not None:
        obj = f()
    return obj


def _is_normalized(dt):
    if (dt.hour != 0 or dt.minute != 0 or dt.second != 0 or
            dt.microsecond != 0 or getattr(dt, 'nanosecond', 0) != 0):
        return False
    return True


# ---------------------------------------------------------------------
# Business Helpers

def _get_firstbday(wkday):
    """
    wkday is the result of monthrange(year, month)

    If it's a saturday or sunday, increment first business day to reflect this
    """
    first = 1
    if wkday == 5:  # on Saturday
        first = 3
    elif wkday == 6:  # on Sunday
        first = 2
    return first


def _get_calendar(weekmask, holidays, calendar):
    """Generate busdaycalendar"""
    if isinstance(calendar, np.busdaycalendar):
        if not holidays:
            holidays = tuple(calendar.holidays)
        elif not isinstance(holidays, tuple):
            holidays = tuple(holidays)
        else:
            # trust that calendar.holidays and holidays are
            # consistent
            pass
        return calendar, holidays

    if holidays is None:
        holidays = []
    try:
        holidays = holidays + calendar.holidays().tolist()
    except AttributeError:
        pass
    holidays = [_to_dt64(dt, dtype='datetime64[D]') for dt in holidays]
    holidays = tuple(sorted(holidays))

    kwargs = {'weekmask': weekmask}
    if holidays:
        kwargs['holidays'] = holidays

    busdaycalendar = np.busdaycalendar(**kwargs)
    return busdaycalendar, holidays


def _to_dt64(dt, dtype='datetime64'):
    # Currently
    # > np.datetime64(dt.datetime(2013,5,1),dtype='datetime64[D]')
    # numpy.datetime64('2013-05-01T02:00:00.000000+0200')
    # Thus astype is needed to cast datetime to datetime64[D]
    if getattr(dt, 'tzinfo', None) is not None:
        i8 = pydt_to_i8(dt)
        dt = tz_convert_single(i8, 'UTC', dt.tzinfo)
        dt = np.int64(dt).astype('datetime64[ns]')
    else:
        dt = np.datetime64(dt)
    if dt.dtype.name != dtype:
        dt = dt.astype(dtype)
    return dt


# ---------------------------------------------------------------------
# Validation


def _validate_business_time(t_input):
    if is_string_object(t_input):
        try:
            t = time.strptime(t_input, '%H:%M')
            return dt_time(hour=t.tm_hour, minute=t.tm_min)
        except ValueError:
            raise ValueError("time data must match '%H:%M' format")
    elif isinstance(t_input, dt_time):
        if t_input.second != 0 or t_input.microsecond != 0:
            raise ValueError(
                "time data must be specified only with hour and minute")
        return t_input
    else:
        raise ValueError("time data must be string or datetime.time")

# ---------------------------------------------------------------------
# Mixins & Singletons


class ApplyTypeError(TypeError):
    # sentinel class for catching the apply error to return NotImplemented
    pass


# TODO: unused.  remove?
class CacheableOffset(object):
    _cacheable = True
