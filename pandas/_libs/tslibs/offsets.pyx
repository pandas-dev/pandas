# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython

import time
from cpython.datetime cimport timedelta, time as dt_time

from dateutil.relativedelta import relativedelta

import numpy as np
cimport numpy as np
from numpy cimport int64_t
np.import_array()


from util cimport is_string_object

from pandas._libs.tslib import pydt_to_i8

from frequencies cimport get_freq_code
from conversion cimport tz_convert_single

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

cpdef int _get_firstbday(int wkday):
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
# Constructor Helpers

_rd_kwds = set([
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


def __unpickle():
    return

class PickleMixin(object):
    """Handle issues with backwards compat related to making DateOffset
    immutable.
    """
    '''
    def __reduce__(self):
        # We need to define this explicitly or else pickle tests fail
        wlist_attrs = tuple(getattr(self, name) for name in sorted(self.kwds))
        tup = (self.n, self.normalize,) + wlist_attrs
        return (self.__class__, tup)
        # FIXME: Isn't this going to screw up on DateOffset?


    def __reduce_ex__(self, protocol):
        # python 3.6 compat
        # http://bugs.python.org/issue28730
        # now __reduce_ex__ is defined and higher priority than __reduce__
        return self.__reduce__()

    '''

    def __getstate__(self):
        """Return a pickleable state"""
        state = self.__dict__.copy()

        # Add attributes from the C base class that aren't in self.__dict__
        state['n'] = self.n
        state['normalize'] = self.normalize

        # we don't want to actually pickle the calendar object
        # as its a np.busyday; we recreate on deserilization
        if 'calendar' in state:
            del state['calendar']
        try:
            state['kwds'].pop('calendar')
        except KeyError:
            pass

        return state


class BusinessMixin(object):
    """ mixin to business types to provide related functions """
    pass


# ---------------------------------------------------------------------
# Base Classes
@cython.auto_pickle(False)
cdef class _BaseOffset(object):
    """
    Base class for DateOffset methods that are not overriden by subclasses
    and will (after pickle errors are resolved) go into a cdef class.
    """
    _typ = "dateoffset"
    _normalize_cache = True
    _cacheable = False

    cdef readonly:
        int64_t n
        bint normalize

    def __init__(self, n, normalize):
        self.n = n
        self.normalize = normalize

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

    def __setstate__(self, state):
        """Reconstruct an instance from a pickled state"""
        # Note: __setstate__ needs to be defined in the cython class otherwise
        # trying to set self.n and self.normalize below will
        # raise an AttributeError.
        if 'normalize' not in state:
            # default for prior pickles
            # See GH #7748, #7789
            state['normalize'] = False

        if 'offset' in state:
            # Older versions have offset attribute instead of _offset
            if '_offset' in state:  # pragma: no cover
                raise ValueError('Unexpected key `_offset`')
            state['_offset'] = state.pop('offset')
            state['kwds']['offset'] = state['_offset']
        
        self.n = state.pop('n')
        self.normalize = state.pop('normalize')
        self.__dict__ = state

        if 'weekmask' in state and 'holidays' in state:
            calendar, holidays = _get_calendar(weekmask=self.weekmask,
                                               holidays=self.holidays,
                                               calendar=None)
            self.kwds['calendar'] = self.calendar = calendar
            self.kwds['holidays'] = self.holidays = holidays
            self.kwds['weekmask'] = state['weekmask']


class BaseOffset(_BaseOffset, PickleMixin):
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
