# -*- coding: utf-8 -*-

import cython
from cython import Py_ssize_t

import time
from cpython.datetime cimport (PyDateTime_IMPORT,
                               PyDateTime_Check,
                               datetime, timedelta,
                               time as dt_time)
PyDateTime_IMPORT

from dateutil.relativedelta import relativedelta

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t
cnp.import_array()


from util cimport is_string_object, is_integer_object

from ccalendar import MONTHS, DAYS
from ccalendar cimport get_days_in_month, dayofweek
from conversion cimport tz_convert_single, pydt_to_i8, localize_pydatetime
from nattype cimport NPY_NAT
from np_datetime cimport (npy_datetimestruct,
                          dtstruct_to_dt64, dt64_to_dtstruct)

# ---------------------------------------------------------------------
# Constants


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
    for _m in MONTHS:
        key = '%s-%s' % (__prefix, _m)
        _offset_to_period_map[key] = _offset_to_period_map[__prefix]

for __prefix in ['A', 'Q']:
    for _m in MONTHS:
        _alias = '%s-%s' % (__prefix, _m)
        _offset_to_period_map[_alias] = _alias

for _d in DAYS:
    _offset_to_period_map['W-%s' % _d] = 'W-%s' % _d


# ---------------------------------------------------------------------
# Misc Helpers

cdef to_offset(object obj):
    """
    Wrap pandas.tseries.frequencies.to_offset to keep centralize runtime
    imports
    """
    from pandas.tseries.frequencies import to_offset
    return to_offset(obj)


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
    # not play nicely with cython class methods
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

relativedelta_kwds = {'years', 'months', 'weeks', 'days', 'year', 'month',
                      'day', 'weekday', 'hour', 'minute', 'second',
                      'microsecond', 'nanosecond', 'nanoseconds', 'hours',
                      'minutes', 'seconds', 'microseconds'}


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


# ---------------------------------------------------------------------
# Base Classes

class _BaseOffset(object):
    """
    Base class for DateOffset methods that are not overridden by subclasses
    and will (after pickle errors are resolved) go into a cdef class.
    """
    _typ = "dateoffset"
    _day_opt = None
    _attributes = frozenset(['n', 'normalize'])

    def __init__(self, n=1, normalize=False):
        n = self._validate_n(n)
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "normalize", normalize)
        object.__setattr__(self, "_cache", {})

    def __setattr__(self, name, value):
        raise AttributeError("DateOffset objects are immutable.")

    def __eq__(self, other):
        if is_string_object(other):
            other = to_offset(other)

        try:
            return self._params == other._params
        except AttributeError:
            # other is not a DateOffset object
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._params)

    @property
    def _params(self):
        """
        Returns a tuple containing all of the attributes needed to evaluate
        equality between two DateOffset objects.
        """
        # NB: non-cython subclasses override property with cache_readonly
        all_paras = self.__dict__.copy()
        if 'holidays' in all_paras and not all_paras['holidays']:
            all_paras.pop('holidays')
        exclude = ['kwds', 'name', 'calendar']
        attrs = [(k, v) for k, v in all_paras.items()
                 if (k not in exclude) and (k[0] != '_')]
        attrs = sorted(set(attrs))
        params = tuple([str(self.__class__)] + attrs)
        return params

    @property
    def kwds(self):
        # for backwards-compatibility
        kwds = {name: getattr(self, name, None) for name in self._attributes
                if name not in ['n', 'normalize']}
        return {name: kwds[name] for name in kwds if kwds[name] is not None}

    def __add__(self, other):
        if getattr(other, "_typ", None) in ["datetimeindex", "periodindex",
                                            "series", "period", "dataframe"]:
            # defer to the other class's implementation
            return other + self
        try:
            return self.apply(other)
        except ApplyTypeError:
            return NotImplemented

    def __sub__(self, other):
        if PyDateTime_Check(other):
            raise TypeError('Cannot subtract datetime from offset.')
        elif type(other) == type(self):
            return type(self)(self.n - other.n, normalize=self.normalize,
                              **self.kwds)
        else:  # pragma: no cover
            return NotImplemented

    def __call__(self, other):
        return self.apply(other)

    def __mul__(self, other):
        return type(self)(n=other * self.n, normalize=self.normalize,
                          **self.kwds)

    def __neg__(self):
        # Note: we are defering directly to __mul__ instead of __rmul__, as
        # that allows us to use methods that can go in a `cdef class`
        return self * -1

    def copy(self):
        # Note: we are defering directly to __mul__ instead of __rmul__, as
        # that allows us to use methods that can go in a `cdef class`
        return self * 1

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

    def __setstate__(self, state):
        """Reconstruct an instance from a pickled state"""
        if 'offset' in state:
            # Older (<0.22.0) versions have offset attribute instead of _offset
            if '_offset' in state:  # pragma: no cover
                raise AssertionError('Unexpected key `_offset`')
            state['_offset'] = state.pop('offset')
            state['kwds']['offset'] = state['_offset']

        if '_offset' in state and not isinstance(state['_offset'], timedelta):
            # relativedelta, we need to populate using its kwds
            offset = state['_offset']
            odict = offset.__dict__
            kwds = {key: odict[key] for key in odict if odict[key]}
            state.update(kwds)

        if '_cache' not in state:
            state['_cache'] = {}

        self.__dict__.update(state)

        if 'weekmask' in state and 'holidays' in state:
            calendar, holidays = _get_calendar(weekmask=self.weekmask,
                                               holidays=self.holidays,
                                               calendar=None)
            object.__setattr__(self, "calendar", calendar)
            object.__setattr__(self, "holidays", holidays)

    def __getstate__(self):
        """Return a pickleable state"""
        state = self.__dict__.copy()

        # we don't want to actually pickle the calendar object
        # as its a np.busyday; we recreate on deserilization
        if 'calendar' in state:
            del state['calendar']
        try:
            state['kwds'].pop('calendar')
        except KeyError:
            pass

        return state


class BaseOffset(_BaseOffset):
    # Here we add __rfoo__ methods that don't play well with cdef classes
    def __rmul__(self, other):
        return self.__mul__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        if getattr(other, '_typ', None) in ['datetimeindex', 'series']:
            # i.e. isinstance(other, (ABCDatetimeIndex, ABCSeries))
            return other - self
        return -self + other


class _Tick(object):
    """
    dummy class to mix into tseries.offsets.Tick so that in tslibs.period we
    can do isinstance checks on _Tick and avoid importing tseries.offsets
    """
    pass


# ----------------------------------------------------------------------
# RelativeDelta Arithmetic

def shift_day(other: datetime, days: int) -> datetime:
    """
    Increment the datetime `other` by the given number of days, retaining
    the time-portion of the datetime.  For tz-naive datetimes this is
    equivalent to adding a timedelta.  For tz-aware datetimes it is similar to
    dateutil's relativedelta.__add__, but handles pytz tzinfo objects.

    Parameters
    ----------
    other : datetime or Timestamp
    days : int

    Returns
    -------
    shifted: datetime or Timestamp
    """
    if other.tzinfo is None:
        return other + timedelta(days=days)

    tz = other.tzinfo
    naive = other.replace(tzinfo=None)
    shifted = naive + timedelta(days=days)
    return localize_pydatetime(shifted, tz)


cdef inline int year_add_months(npy_datetimestruct dts, int months) nogil:
    """new year number after shifting npy_datetimestruct number of months"""
    return dts.year + (dts.month + months - 1) / 12


cdef inline int month_add_months(npy_datetimestruct dts, int months) nogil:
    """
    New month number after shifting npy_datetimestruct
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
        npy_datetimestruct dts
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
                # compare_day is only relevant for comparison in the case
                # where months_since == 0.
                compare_day = get_firstbday(dts.year, dts.month)

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
                # compare_day is only relevant for comparison in the case
                # where months_since == 0.
                compare_day = get_lastbday(dts.year, dts.month)

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
        npy_datetimestruct dts
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

                months_to_roll = roll_convention(dts.day, months_to_roll,
                                                 compare_day)

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

                months_to_roll = roll_convention(dts.day, months_to_roll,
                                                 compare_day)

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_lastbday(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    else:
        raise ValueError("day must be None, 'start', 'end', "
                         "'business_start', or 'business_end'")

    return np.asarray(out)


def shift_month(stamp: datetime, months: int,
                day_opt: object = None) -> datetime:
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
        'end': returns last day of the month

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


cpdef int roll_convention(int other, int n, int compare) nogil:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : int, generally the day component of a datetime
    n : number of periods to increment, before adjusting for rolling
    compare : int, generally the day component of a datetime, in the same
              month as the datetime form which `other` was taken.

    Returns
    -------
    n : int number of periods to increment
    """
    if n > 0 and other < compare:
        n -= 1
    elif n <= 0 and other > compare:
        # as if rolled forward already
        n += 1
    return n


def roll_qtrday(other: datetime, n: int, month: int,
                day_opt: object, modby: int = 3) -> int:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : datetime or Timestamp
    n : number of periods to increment, before adjusting for rolling
    month : int reference month giving the first month of the year
    day_opt : 'start', 'end', 'business_start', 'business_end'
        The convention to use in finding the day in a given month against
        which to compare for rollforward/rollbackward decisions.
    modby : int 3 for quarters, 12 for years

    Returns
    -------
    n : int number of periods to increment
    """
    cdef:
        int months_since
    # TODO: Merge this with roll_yearday by setting modby=12 there?
    #       code de-duplication versus perf hit?
    # TODO: with small adjustments this could be used in shift_quarters
    months_since = other.month % modby - month % modby

    if n > 0:
        if months_since < 0 or (months_since == 0 and
                                other.day < get_day_of_month(other,
                                                             day_opt)):
            # pretend to roll back if on same month but
            # before compare_day
            n -= 1
    else:
        if months_since > 0 or (months_since == 0 and
                                other.day > get_day_of_month(other,
                                                             day_opt)):
            # make sure to roll forward, so negate
            n += 1
    return n


def roll_yearday(other: datetime, n: int, month: int, day_opt: object) -> int:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : datetime or Timestamp
    n : number of periods to increment, before adjusting for rolling
    month : reference month giving the first month of the year
    day_opt : 'start', 'end'
        'start': returns 1
        'end': returns last day of the month

    Returns
    -------
    n : int number of periods to increment

    Notes
    -----
    * Mirrors `roll_check` in shift_months

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
    else:
        if other.month > month or (other.month == month and
                                   other.day > get_day_of_month(other,
                                                                day_opt)):
            n += 1
    return n
