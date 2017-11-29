# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython
from cython cimport Py_ssize_t

import time
from cpython.datetime cimport datetime, timedelta, time as dt_time

from dateutil.relativedelta import relativedelta

import numpy as np
cimport numpy as np
from numpy cimport int64_t
np.import_array()


from util cimport is_string_object, is_integer_object

from conversion cimport tz_convert_single, pydt_to_i8
from frequencies cimport get_freq_code
from nattype cimport NPY_NAT
from np_datetime cimport (pandas_datetimestruct,
                          dtstruct_to_dt64, dt64_to_dtstruct,
                          is_leapyear, days_per_month_table, dayofweek)

# ---------------------------------------------------------------------
# Constants

# Duplicated in tslib
_MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
           'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
_int_to_month = {(k + 1): v for k, v in enumerate(_MONTHS)}
_month_to_int = {v: k for k, v in _int_to_month.items()}


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


cpdef bint _is_normalized(dt):
    if (dt.hour != 0 or dt.minute != 0 or dt.second != 0 or
            dt.microsecond != 0 or getattr(dt, 'nanosecond', 0) != 0):
        return False
    return True


def apply_index_wraps(func):
    # Note: normally we would use `@functools.wraps(func)`, but this does
    # not play nicely wtih cython class methods
    def wrapper(self, other):
        result = func(self, other)
        if self.normalize:
            result = result.to_period('D').to_timestamp()
        return result

    # do @functools.wraps(func) manually since it doesn't work on cdef funcs
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    try:
        wrapper.__module__ = func.__module__
    except AttributeError:
        # AttributeError: 'method_descriptor' object has no
        # attribute '__module__'
        pass
    return wrapper


# ---------------------------------------------------------------------
# Business Helpers

cpdef int get_lastbday(int year, int month) nogil:
    """
    Find the last day of the month that is a business day.

    Parameters
    ----------
    year : int
    month : int

    Returns
    -------
    last_bday : int
    """
    cdef:
        int wkday, days_in_month

    wkday = dayofweek(year, month, 1)
    days_in_month = get_days_in_month(year, month)
    return days_in_month - max(((wkday + days_in_month - 1) % 7) - 4, 0)


cpdef int get_firstbday(int year, int month) nogil:
    """
    Find the first day of the month that is a business day.

    Parameters
    ----------
    year : int
    month : int

    Returns
    -------
    first_bday : int
    """
    cdef:
        int first, wkday

    wkday = dayofweek(year, month, 1)
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
# Constructor Helpers

relativedelta_kwds = set([
    'years', 'months', 'weeks', 'days',
    'year', 'month', 'week', 'day', 'weekday',
    'hour', 'minute', 'second', 'microsecond',
    'nanosecond', 'nanoseconds',
    'hours', 'minutes', 'seconds', 'milliseconds', 'microseconds'])


def _determine_offset(kwds):
    # timedelta is used for sub-daily plural offsets and all singular
    # offsets relativedelta is used for plural offsets of daily length or
    # more nanosecond(s) are handled by apply_wraps
    kwds_no_nanos = dict(
        (k, v) for k, v in kwds.items()
        if k not in ('nanosecond', 'nanoseconds')
    )
    # TODO: Are nanosecond and nanoseconds allowed somewhere?

    _kwds_use_relativedelta = ('years', 'months', 'weeks', 'days',
                               'year', 'month', 'week', 'day', 'weekday',
                               'hour', 'minute', 'second', 'microsecond')

    use_relativedelta = False
    if len(kwds_no_nanos) > 0:
        if any(k in _kwds_use_relativedelta for k in kwds_no_nanos):
            offset = relativedelta(**kwds_no_nanos)
            use_relativedelta = True
        else:
            # sub-daily offset - use timedelta (tz-aware)
            offset = timedelta(**kwds_no_nanos)
    else:
        offset = timedelta(1)
    return offset, use_relativedelta


# ---------------------------------------------------------------------
# Mixins & Singletons


class ApplyTypeError(TypeError):
    # sentinel class for catching the apply error to return NotImplemented
    pass


# TODO: unused.  remove?
class CacheableOffset(object):
    _cacheable = True


class BeginMixin(object):
    # helper for vectorized offsets

    def _beg_apply_index(self, i, freq):
        """Offsets index to beginning of Period frequency"""

        off = i.to_perioddelta('D')

        base, mult = get_freq_code(freq)
        base_period = i.to_period(base)
        if self.n <= 0:
            # when subtracting, dates on start roll to prior
            roll = np.where(base_period.to_timestamp() == i - off,
                            self.n, self.n + 1)
        else:
            roll = self.n

        base = (base_period + roll).to_timestamp()
        return base + off


class EndMixin(object):
    # helper for vectorized offsets

    def _end_apply_index(self, i, freq):
        """Offsets index to end of Period frequency"""

        off = i.to_perioddelta('D')

        base, mult = get_freq_code(freq)
        base_period = i.to_period(base)
        if self.n > 0:
            # when adding, dates on end roll to next
            roll = np.where(base_period.to_timestamp(how='end') == i - off,
                            self.n, self.n - 1)
        else:
            roll = self.n

        base = (base_period + roll).to_timestamp(how='end')
        return base + off


# ---------------------------------------------------------------------
# Base Classes

class _BaseOffset(object):
    """
    Base class for DateOffset methods that are not overriden by subclasses
    and will (after pickle errors are resolved) go into a cdef class.
    """
    _typ = "dateoffset"
    _normalize_cache = True
    _cacheable = False
    _day_opt = None

    def __call__(self, other):
        return self.apply(other)

    def __mul__(self, someInt):
        return self.__class__(n=someInt * self.n, normalize=self.normalize,
                              **self.kwds)

    def __neg__(self):
        # Note: we are defering directly to __mul__ instead of __rmul__, as
        # that allows us to use methods that can go in a `cdef class`
        return self * -1

    def copy(self):
        # Note: we are defering directly to __mul__ instead of __rmul__, as
        # that allows us to use methods that can go in a `cdef class`
        return self * 1

    # TODO: this is never true.  fix it or get rid of it
    def _should_cache(self):
        return self.isAnchored() and self._cacheable

    def __repr__(self):
        className = getattr(self, '_outputName', type(self).__name__)

        if abs(self.n) != 1:
            plural = 's'
        else:
            plural = ''

        n_str = ""
        if self.n != 1:
            n_str = "%s * " % self.n

        out = '<%s' % n_str + className + plural + self._repr_attrs() + '>'
        return out

    def _get_offset_day(self, datetime other):
        # subclass must implement `_day_opt`; calling from the base class
        # will raise NotImplementedError.
        return get_day_of_month(other, self._day_opt)

    def _validate_n(self, n):
        """
        Require that `n` be a nonzero integer.

        Parameters
        ----------
        n : int

        Returns
        -------
        nint : int

        Raises
        ------
        TypeError if `int(n)` raises
        ValueError if n != int(n)
        """
        try:
            nint = int(n)
        except (ValueError, TypeError):
            raise TypeError('`n` argument must be an integer, '
                            'got {ntype}'.format(ntype=type(n)))
        if n != nint:
            raise ValueError('`n` argument must be an integer, '
                             'got {n}'.format(n=n))
        return nint


class BaseOffset(_BaseOffset):
    # Here we add __rfoo__ methods that don't play well with cdef classes
    def __rmul__(self, someInt):
        return self.__mul__(someInt)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        if getattr(other, '_typ', None) in ['datetimeindex', 'series']:
            # i.e. isinstance(other, (ABCDatetimeIndex, ABCSeries))
            return other - self
        return -self + other


# ----------------------------------------------------------------------
# RelativeDelta Arithmetic

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int get_days_in_month(int year, int month) nogil:
    return days_per_month_table[is_leapyear(year)][month - 1]


cdef inline int year_add_months(pandas_datetimestruct dts, int months) nogil:
    """new year number after shifting pandas_datetimestruct number of months"""
    return dts.year + (dts.month + months - 1) / 12


cdef inline int month_add_months(pandas_datetimestruct dts, int months) nogil:
    """
    New month number after shifting pandas_datetimestruct
    number of months.
    """
    cdef int new_month = (dts.month + months) % 12
    return 12 if new_month == 0 else new_month


@cython.wraparound(False)
@cython.boundscheck(False)
def shift_quarters(int64_t[:] dtindex, int quarters,
                   int q1start_month, object day, int modby=3):
    """
    Given an int64 array representing nanosecond timestamps, shift all elements
    by the specified number of quarters using DateOffset semantics.

    Parameters
    ----------
    dtindex : int64_t[:] timestamps for input dates
    quarters : int number of quarters to shift
    q1start_month : int month in which Q1 begins by convention
    day : {'start', 'end', 'business_start', 'business_end'}
    modby : int (3 for quarters, 12 for years)

    Returns
    -------
    out : ndarray[int64_t]
    """
    cdef:
        Py_ssize_t i
        pandas_datetimestruct dts
        int count = len(dtindex)
        int months_to_roll, months_since, n, compare_day
        bint roll_check
        int64_t[:] out = np.empty(count, dtype='int64')

    if day == 'start':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby

                # offset semantics - if on the anchor point and going backwards
                # shift to next
                if n <= 0 and (months_since != 0 or
                               (months_since == 0 and dts.day > 1)):
                    n += 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)
                dts.day = 1

                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'end':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby

                if n <= 0 and months_since != 0:
                    # The general case of this condition would be
                    # `months_since != 0 or (months_since == 0 and
                    #    dts.day > get_days_in_month(dts.year, dts.month))`
                    # but the get_days_in_month inequality would never hold.
                    n += 1
                elif n > 0 and (months_since == 0 and
                                dts.day < get_days_in_month(dts.year,
                                                            dts.month)):
                    n -= 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)
                dts.day = get_days_in_month(dts.year, dts.month)

                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_start':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby
                compare_month = dts.month - months_since
                compare_month = compare_month or 12
                # compare_day is only relevant for comparison in the case
                # where months_since == 0.
                compare_day = get_firstbday(dts.year, compare_month)

                if n <= 0 and (months_since != 0 or
                               (months_since == 0 and dts.day > compare_day)):
                    # make sure to roll forward, so negate
                    n += 1
                elif n > 0 and (months_since == 0 and dts.day < compare_day):
                    # pretend to roll back if on same month but
                    # before compare_day
                    n -= 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)

                dts.day = get_firstbday(dts.year, dts.month)

                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_end':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby
                compare_month = dts.month - months_since
                compare_month = compare_month or 12
                # compare_day is only relevant for comparison in the case
                # where months_since == 0.
                compare_day = get_lastbday(dts.year, compare_month)

                if n <= 0 and (months_since != 0 or
                               (months_since == 0 and dts.day > compare_day)):
                    # make sure to roll forward, so negate
                    n += 1
                elif n > 0 and (months_since == 0 and dts.day < compare_day):
                    # pretend to roll back if on same month but
                    # before compare_day
                    n -= 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)

                dts.day = get_lastbday(dts.year, dts.month)

                out[i] = dtstruct_to_dt64(&dts)

    else:
        raise ValueError("day must be None, 'start', 'end', "
                         "'business_start', or 'business_end'")

    return np.asarray(out)


@cython.wraparound(False)
@cython.boundscheck(False)
def shift_months(int64_t[:] dtindex, int months, object day=None):
    """
    Given an int64-based datetime index, shift all elements
    specified number of months using DateOffset semantics

    day: {None, 'start', 'end'}
       * None: day of month
       * 'start' 1st day of month
       * 'end' last day of month
    """
    cdef:
        Py_ssize_t i
        pandas_datetimestruct dts
        int count = len(dtindex)
        int months_to_roll
        bint roll_check
        int64_t[:] out = np.empty(count, dtype='int64')

    if day is None:
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                dts.year = year_add_months(dts, months)
                dts.month = month_add_months(dts, months)

                dts.day = min(dts.day, get_days_in_month(dts.year, dts.month))
                out[i] = dtstruct_to_dt64(&dts)
    elif day == 'start':
        roll_check = False
        if months <= 0:
            months += 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months

                # offset semantics - if on the anchor point and going backwards
                # shift to next
                if roll_check and dts.day == 1:
                    months_to_roll -= 1

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)
                dts.day = 1

                out[i] = dtstruct_to_dt64(&dts)
    elif day == 'end':
        roll_check = False
        if months > 0:
            months -= 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months

                # similar semantics - when adding shift forward by one
                # month if already at an end of month
                if roll_check and dts.day == get_days_in_month(dts.year,
                                                               dts.month):
                    months_to_roll += 1

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_days_in_month(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_start':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months
                compare_day = get_firstbday(dts.year, dts.month)

                if months_to_roll > 0 and dts.day < compare_day:
                    months_to_roll -= 1
                elif months_to_roll <= 0 and dts.day > compare_day:
                    # as if rolled forward already
                    months_to_roll += 1

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_firstbday(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_end':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months
                compare_day = get_lastbday(dts.year, dts.month)

                if months_to_roll > 0 and dts.day < compare_day:
                    months_to_roll -= 1
                elif months_to_roll <= 0 and dts.day > compare_day:
                    # as if rolled forward already
                    months_to_roll += 1

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_lastbday(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    else:
        raise ValueError("day must be None, 'start', 'end', "
                         "'business_start', or 'business_end'")

    return np.asarray(out)


cpdef datetime shift_month(datetime stamp, int months, object day_opt=None):
    """
    Given a datetime (or Timestamp) `stamp`, an integer `months` and an
    option `day_opt`, return a new datetimelike that many months later,
    with day determined by `day_opt` using relativedelta semantics.

    Scalar analogue of shift_months

    Parameters
    ----------
    stamp : datetime or Timestamp
    months : int
    day_opt : None, 'start', 'end', or an integer
        None: returned datetimelike has the same day as the input, or the
              last day of the month if the new month is too short
        'start': returned datetimelike has day=1
        'end': returned datetimelike has day on the last day of the month
        int: returned datetimelike has day equal to day_opt

    Returns
    -------
    shifted : datetime or Timestamp (same as input `stamp`)
    """
    cdef:
        int year, month, day
        int days_in_month, dy

    dy = (stamp.month + months) // 12
    month = (stamp.month + months) % 12

    if month == 0:
        month = 12
        dy -= 1
    year = stamp.year + dy

    if day_opt is None:
        days_in_month = get_days_in_month(year, month)
        day = min(stamp.day, days_in_month)
    elif day_opt == 'start':
        day = 1
    elif day_opt == 'end':
        day = get_days_in_month(year, month)
    elif day_opt == 'business_start':
        # first business day of month
        day = get_firstbday(year, month)
    elif day_opt == 'business_end':
        # last business day of month
        day = get_lastbday(year, month)
    elif is_integer_object(day_opt):
        days_in_month = get_days_in_month(year, month)
        day = min(day_opt, days_in_month)
    else:
        raise ValueError(day_opt)
    return stamp.replace(year=year, month=month, day=day)


cpdef int get_day_of_month(datetime other, day_opt) except? -1:
    """
    Find the day in `other`'s month that satisfies a DateOffset's onOffset
    policy, as described by the `day_opt` argument.

    Parameters
    ----------
    other : datetime or Timestamp
    day_opt : 'start', 'end'
        'start': returns 1
        'end': returns last day  of the month

    Returns
    -------
    day_of_month : int

    Examples
    -------
    >>> other = datetime(2017, 11, 14)
    >>> get_day_of_month(other, 'start')
    1
    >>> get_day_of_month(other, 'end')
    30

    """
    cdef:
        int days_in_month

    if day_opt == 'start':
        return 1
    elif day_opt == 'end':
        days_in_month = get_days_in_month(other.year, other.month)
        return days_in_month
    elif day_opt == 'business_start':
        # first business day of month
        return get_firstbday(other.year, other.month)
    elif day_opt == 'business_end':
        # last business day of month
        return get_lastbday(other.year, other.month)
    elif is_integer_object(day_opt):
        days_in_month = get_days_in_month(other.year, other.month)
        return min(day_opt, days_in_month)
    elif day_opt is None:
        # Note: unlike `shift_month`, get_day_of_month does not
        # allow day_opt = None
        raise NotImplementedError
    else:
        raise ValueError(day_opt)


cpdef int roll_yearday(other, n, month, day_opt='start') except? -1:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : datetime or Timestamp
    n : number of periods to increment, before adjusting for rolling
    day_opt : 'start', 'end'
        'start': returns 1
        'end': returns last day  of the month

    Returns
    -------
    n : int number of periods to increment

    Notes
    -----
    * Mirrors `roll_check` in tslib.shift_months

    Examples
    -------
    >>> month = 3
    >>> day_opt = 'start'              # `other` will be compared to March 1
    >>> other = datetime(2017, 2, 10)  # before March 1
    >>> roll_yearday(other, 2, month, day_opt)
    1
    >>> roll_yearday(other, -7, month, day_opt)
    -7
    >>>
    >>> other = Timestamp('2014-03-15', tz='US/Eastern')  # after March 1
    >>> roll_yearday(other, 2, month, day_opt)
    2
    >>> roll_yearday(other, -7, month, day_opt)
    -6

    >>> month = 6
    >>> day_opt = 'end'                # `other` will be compared to June 30
    >>> other = datetime(1999, 6, 29)  # before June 30
    >>> roll_yearday(other, 5, month, day_opt)
    4
    >>> roll_yearday(other, -7, month, day_opt)
    -7
    >>>
    >>> other = Timestamp(2072, 8, 24, 6, 17, 18)  # after June 30
    >>> roll_yearday(other, 5, month, day_opt)
    5
    >>> roll_yearday(other, -7, month, day_opt)
    -6

    """
    # Note: The other.day < ... condition will never hold when day_opt=='start'
    # and the other.day > ... condition will never hold when day_opt=='end'.
    # At some point these extra checks may need to be optimized away.
    # But that point isn't today.
    if n > 0:
        if other.month < month or (other.month == month and
                                   other.day < get_day_of_month(other,
                                                                day_opt)):
            n -= 1
    elif n <= 0:
        if other.month > month or (other.month == month and
                                   other.day > get_day_of_month(other,
                                                                day_opt)):
            n += 1
    return n
