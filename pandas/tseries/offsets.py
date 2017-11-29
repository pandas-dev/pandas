# -*- coding: utf-8 -*-
from datetime import date, datetime, timedelta
import functools
import operator

from pandas.compat import range
from pandas import compat
import numpy as np

from pandas.core.dtypes.generic import ABCSeries, ABCDatetimeIndex, ABCPeriod
from pandas.core.tools.datetimes import to_datetime, normalize_date
from pandas.core.common import AbstractMethodError

# import after tools, dateutil check
from dateutil.easter import easter
from pandas._libs import tslib, Timestamp, OutOfBoundsDatetime, Timedelta
from pandas.util._decorators import cache_readonly

from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds
import pandas._libs.tslibs.offsets as liboffsets
from pandas._libs.tslibs.offsets import (
    ApplyTypeError,
    as_datetime, _is_normalized,
    _get_calendar, _to_dt64, _validate_business_time,
    _int_to_weekday, _weekday_to_int,
    _determine_offset,
    apply_index_wraps,
    roll_yearday,
    shift_month,
    EndMixin,
    BaseOffset)


__all__ = ['Day', 'BusinessDay', 'BDay', 'CustomBusinessDay', 'CDay',
           'CBMonthEnd', 'CBMonthBegin',
           'MonthBegin', 'BMonthBegin', 'MonthEnd', 'BMonthEnd',
           'SemiMonthEnd', 'SemiMonthBegin',
           'BusinessHour', 'CustomBusinessHour',
           'YearBegin', 'BYearBegin', 'YearEnd', 'BYearEnd',
           'QuarterBegin', 'BQuarterBegin', 'QuarterEnd', 'BQuarterEnd',
           'LastWeekOfMonth', 'FY5253Quarter', 'FY5253',
           'Week', 'WeekOfMonth', 'Easter',
           'Hour', 'Minute', 'Second', 'Milli', 'Micro', 'Nano',
           'DateOffset']

# convert to/from datetime/timestamp to allow invalid Timestamp ranges to
# pass thru


def as_timestamp(obj):
    if isinstance(obj, Timestamp):
        return obj
    try:
        return Timestamp(obj)
    except (OutOfBoundsDatetime):
        pass
    return obj


def apply_wraps(func):
    @functools.wraps(func)
    def wrapper(self, other):
        if other is tslib.NaT:
            return tslib.NaT
        elif isinstance(other, (timedelta, Tick, DateOffset)):
            # timedelta path
            return func(self, other)
        elif isinstance(other, (np.datetime64, datetime, date)):
            other = as_timestamp(other)

        tz = getattr(other, 'tzinfo', None)
        nano = getattr(other, 'nanosecond', 0)

        try:
            if self._adjust_dst and isinstance(other, Timestamp):
                other = other.tz_localize(None)

            result = func(self, other)

            if self._adjust_dst:
                result = tslib._localize_pydatetime(result, tz)

            result = Timestamp(result)
            if self.normalize:
                result = result.normalize()

            # nanosecond may be deleted depending on offset process
            if not self.normalize and nano != 0:
                if not isinstance(self, Nano) and result.nanosecond != nano:
                    if result.tz is not None:
                        # convert to UTC
                        value = tslib.tz_convert_single(
                            result.value, 'UTC', result.tz)
                    else:
                        value = result.value
                    result = Timestamp(value + nano)

            if tz is not None and result.tzinfo is None:
                result = tslib._localize_pydatetime(result, tz)

        except OutOfBoundsDatetime:
            result = func(self, as_datetime(other))

            if self.normalize:
                # normalize_date returns normal datetime
                result = normalize_date(result)

            if tz is not None and result.tzinfo is None:
                result = tslib._localize_pydatetime(result, tz)

        return result
    return wrapper


# ---------------------------------------------------------------------
# DateOffset


class DateOffset(BaseOffset):
    """
    Standard kind of date increment used for a date range.

    Works exactly like relativedelta in terms of the keyword args you
    pass in, use of the keyword n is discouraged-- you would be better
    off specifying n in the keywords you use, but regardless it is
    there for you. n is needed for DateOffset subclasses.

    DateOffets work as follows.  Each offset specify a set of dates
    that conform to the DateOffset.  For example, Bday defines this
    set to be the set of dates that are weekdays (M-F).  To test if a
    date is in the set of a DateOffset dateOffset we can use the
    onOffset method: dateOffset.onOffset(date).

    If a date is not on a valid date, the rollback and rollforward
    methods can be used to roll the date to the nearest valid date
    before/after the date.

    DateOffsets can be created to move dates forward a given number of
    valid dates.  For example, Bday(2) can be added to a date to move
    it two business days forward.  If the date does not start on a
    valid date, first it is moved to a valid date.  Thus pseudo code
    is:

    def __add__(date):
      date = rollback(date) # does nothing if date is valid
      return date + <n number of periods>

    When a date offset is created for a negative number of periods,
    the date is first rolled forward.  The pseudo code is:

    def __add__(date):
      date = rollforward(date) # does nothing is date is valid
      return date + <n number of periods>

    Zero presents a problem.  Should it roll forward or back?  We
    arbitrarily have it rollforward:

    date + BDay(0) == BDay.rollforward(date)

    Since 0 is a bit weird, we suggest avoiding its use.
    """
    _use_relativedelta = False
    _adjust_dst = False

    # default for prior pickles
    normalize = False

    def __init__(self, n=1, normalize=False, **kwds):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.kwds = kwds

        self._offset, self._use_relativedelta = _determine_offset(kwds)

    @apply_wraps
    def apply(self, other):
        if self._use_relativedelta:
            other = as_datetime(other)

        if len(self.kwds) > 0:
            tzinfo = getattr(other, 'tzinfo', None)
            if tzinfo is not None and self._use_relativedelta:
                # perform calculation in UTC
                other = other.replace(tzinfo=None)

            if self.n > 0:
                for i in range(self.n):
                    other = other + self._offset
            else:
                for i in range(-self.n):
                    other = other - self._offset

            if tzinfo is not None and self._use_relativedelta:
                # bring tz back from UTC calculation
                other = tslib._localize_pydatetime(other, tzinfo)

            return as_timestamp(other)
        else:
            return other + timedelta(self.n)

    @apply_index_wraps
    def apply_index(self, i):
        """
        Vectorized apply of DateOffset to DatetimeIndex,
        raises NotImplentedError for offsets without a
        vectorized implementation

        Parameters
        ----------
        i : DatetimeIndex

        Returns
        -------
        y : DatetimeIndex
        """

        if not type(self) is DateOffset:
            raise NotImplementedError("DateOffset subclass {name} "
                                      "does not have a vectorized "
                                      "implementation".format(
                                          name=self.__class__.__name__))
        relativedelta_fast = set(['years', 'months', 'weeks',
                                  'days', 'hours', 'minutes',
                                  'seconds', 'microseconds'])
        # relativedelta/_offset path only valid for base DateOffset
        if (self._use_relativedelta and
                set(self.kwds).issubset(relativedelta_fast)):

            months = ((self.kwds.get('years', 0) * 12 +
                       self.kwds.get('months', 0)) * self.n)
            if months:
                shifted = liboffsets.shift_months(i.asi8, months)
                i = i._shallow_copy(shifted)

            weeks = (self.kwds.get('weeks', 0)) * self.n
            if weeks:
                i = (i.to_period('W') + weeks).to_timestamp() + \
                    i.to_perioddelta('W')

            timedelta_kwds = {k: v for k, v in self.kwds.items()
                              if k in ['days', 'hours', 'minutes',
                                       'seconds', 'microseconds']}
            if timedelta_kwds:
                delta = Timedelta(**timedelta_kwds)
                i = i + (self.n * delta)
            return i
        elif not self._use_relativedelta and hasattr(self, '_offset'):
            # timedelta
            return i + (self._offset * self.n)
        else:
            # relativedelta with other keywords
            kwd = set(self.kwds) - relativedelta_fast
            raise NotImplementedError("DateOffset with relativedelta "
                                      "keyword(s) {kwd} not able to be "
                                      "applied vectorized".format(kwd=kwd))

    def isAnchored(self):
        # TODO: Does this make sense for the general case?  It would help
        # if there were a canonical docstring for what isAnchored means.
        return (self.n == 1)

    def _params(self):
        all_paras = dict(list(vars(self).items()) + list(self.kwds.items()))
        if 'holidays' in all_paras and not all_paras['holidays']:
            all_paras.pop('holidays')
        exclude = ['kwds', 'name', 'normalize', 'calendar']
        attrs = [(k, v) for k, v in all_paras.items()
                 if (k not in exclude) and (k[0] != '_')]
        attrs = sorted(set(attrs))
        params = tuple([str(self.__class__)] + attrs)
        return params

    # TODO: Combine this with BusinessMixin version by defining a whitelisted
    # set of attributes on each object rather than the existing behavior of
    # iterating over internal ``__dict__``
    def _repr_attrs(self):
        exclude = set(['n', 'inc', 'normalize'])
        attrs = []
        for attr in sorted(self.__dict__):
            if attr.startswith('_'):
                continue
            elif attr == 'kwds':  # TODO: get rid of this
                kwds_new = {}
                for key in self.kwds:
                    if not hasattr(self, key):
                        kwds_new[key] = self.kwds[key]
                if len(kwds_new) > 0:
                    attrs.append('kwds={kwds_new}'.format(kwds_new=kwds_new))
            elif attr not in exclude:
                value = getattr(self, attr)
                attrs.append('{attr}={value}'.format(attr=attr, value=value))

        out = ''
        if attrs:
            out += ': ' + ', '.join(attrs)
        return out

    @property
    def name(self):
        return self.rule_code

    def __eq__(self, other):
        if other is None:
            return False

        if isinstance(other, compat.string_types):
            from pandas.tseries.frequencies import to_offset

            other = to_offset(other)

        if not isinstance(other, DateOffset):
            return False

        return self._params() == other._params()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._params())

    def __add__(self, other):
        if isinstance(other, (ABCDatetimeIndex, ABCSeries)):
            return other + self
        elif isinstance(other, ABCPeriod):
            return other + self
        try:
            return self.apply(other)
        except ApplyTypeError:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, datetime):
            raise TypeError('Cannot subtract datetime from offset.')
        elif type(other) == type(self):
            return self.__class__(self.n - other.n, normalize=self.normalize,
                                  **self.kwds)
        else:  # pragma: no cover
            return NotImplemented

    def rollback(self, dt):
        """Roll provided date backward to next offset only if not on offset"""
        dt = as_timestamp(dt)
        if not self.onOffset(dt):
            dt = dt - self.__class__(1, normalize=self.normalize, **self.kwds)
        return dt

    def rollforward(self, dt):
        """Roll provided date forward to next offset only if not on offset"""
        dt = as_timestamp(dt)
        if not self.onOffset(dt):
            dt = dt + self.__class__(1, normalize=self.normalize, **self.kwds)
        return dt

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        # XXX, see #1395
        if type(self) == DateOffset or isinstance(self, Tick):
            return True

        # Default (slow) method for determining if some date is a member of the
        # date range generated by this offset. Subclasses may have this
        # re-implemented in a nicer way.
        a = dt
        b = ((dt + self) - self)
        return a == b

    # way to get around weirdness with rule_code
    @property
    def _prefix(self):
        raise NotImplementedError('Prefix not defined')

    @property
    def rule_code(self):
        return self._prefix

    @property
    def freqstr(self):
        try:
            code = self.rule_code
        except NotImplementedError:
            return repr(self)

        if self.n != 1:
            fstr = '{n}{code}'.format(n=self.n, code=code)
        else:
            fstr = code

        try:
            if self._offset:
                fstr += self._offset_str()
        except AttributeError:
            # TODO: standardize `_offset` vs `offset` naming convention
            pass

        return fstr

    def _offset_str(self):
        return ''

    @property
    def nanos(self):
        raise ValueError("{name} is a non-fixed frequency".format(name=self))


class SingleConstructorOffset(DateOffset):
    @classmethod
    def _from_name(cls, suffix=None):
        # default _from_name calls cls with no args
        if suffix:
            raise ValueError("Bad freq suffix {suffix}".format(suffix=suffix))
        return cls()


class BusinessMixin(object):
    """ mixin to business types to provide related functions """

    @property
    def offset(self):
        """Alias for self._offset"""
        # Alias for backward compat
        return self._offset

    def _repr_attrs(self):
        if self.offset:
            attrs = ['offset={offset!r}'.format(offset=self.offset)]
        else:
            attrs = None
        out = ''
        if attrs:
            out += ': ' + ', '.join(attrs)
        return out

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

    def __setstate__(self, state):
        """Reconstruct an instance from a pickled state"""
        if 'offset' in state:
            # Older versions have offset attribute instead of _offset
            if '_offset' in state:  # pragma: no cover
                raise ValueError('Unexpected key `_offset`')
            state['_offset'] = state.pop('offset')
            state['kwds']['offset'] = state['_offset']
        self.__dict__ = state
        if 'weekmask' in state and 'holidays' in state:
            calendar, holidays = _get_calendar(weekmask=self.weekmask,
                                               holidays=self.holidays,
                                               calendar=None)
            self.kwds['calendar'] = self.calendar = calendar
            self.kwds['holidays'] = self.holidays = holidays
            self.kwds['weekmask'] = state['weekmask']


class BusinessDay(BusinessMixin, SingleConstructorOffset):
    """
    DateOffset subclass representing possibly n business days
    """
    _prefix = 'B'
    _adjust_dst = True

    def __init__(self, n=1, normalize=False, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.kwds = {'offset': offset}
        self._offset = offset

    def _offset_str(self):
        def get_str(td):
            off_str = ''
            if td.days > 0:
                off_str += str(td.days) + 'D'
            if td.seconds > 0:
                s = td.seconds
                hrs = int(s / 3600)
                if hrs != 0:
                    off_str += str(hrs) + 'H'
                    s -= hrs * 3600
                mts = int(s / 60)
                if mts != 0:
                    off_str += str(mts) + 'Min'
                    s -= mts * 60
                if s != 0:
                    off_str += str(s) + 's'
            if td.microseconds > 0:
                off_str += str(td.microseconds) + 'us'
            return off_str

        if isinstance(self.offset, timedelta):
            zero = timedelta(0, 0, 0)
            if self.offset >= zero:
                off_str = '+' + get_str(self.offset)
            else:
                off_str = '-' + get_str(-self.offset)
            return off_str
        else:
            return '+' + repr(self.offset)

    @apply_wraps
    def apply(self, other):
        if isinstance(other, datetime):
            n = self.n

            if n == 0 and other.weekday() > 4:
                n = 1

            result = other

            # avoid slowness below
            if abs(n) > 5:
                k = n // 5
                result = result + timedelta(7 * k)
                if n < 0 and result.weekday() > 4:
                    n += 1
                n -= 5 * k
                if n == 0 and result.weekday() > 4:
                    n -= 1

            while n != 0:
                k = n // abs(n)
                result = result + timedelta(k)
                if result.weekday() < 5:
                    n -= k

            if self.offset:
                result = result + self.offset
            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other,
                        normalize=self.normalize)
        else:
            raise ApplyTypeError('Only know how to combine business day with '
                                 'datetime or timedelta.')

    @apply_index_wraps
    def apply_index(self, i):
        time = i.to_perioddelta('D')
        # to_period rolls forward to next BDay; track and
        # reduce n where it does when rolling forward
        shifted = (i.to_perioddelta('B') - time).asi8 != 0
        if self.n > 0:
            roll = np.where(shifted, self.n - 1, self.n)
        else:
            roll = self.n

        return (i.to_period('B') + roll).to_timestamp() + time

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.weekday() < 5


class BusinessHourMixin(BusinessMixin):

    def __init__(self, start='09:00', end='17:00', offset=timedelta(0)):
        # must be validated here to equality check
        kwds = {'offset': offset}
        self.start = kwds['start'] = _validate_business_time(start)
        self.end = kwds['end'] = _validate_business_time(end)
        self.kwds = kwds
        self._offset = offset

    def _get_daytime_flag(self):
        if self.start == self.end:
            raise ValueError('start and end must not be the same')
        elif self.start < self.end:
            return True
        else:
            return False

    def _next_opening_time(self, other):
        """
        If n is positive, return tomorrow's business day opening time.
        Otherwise yesterday's business day's opening time.

        Opening time always locates on BusinessDay.
        Otherwise, closing time may not if business hour extends over midnight.
        """
        if not self.next_bday.onOffset(other):
            other = other + self.next_bday
        else:
            if self.n >= 0 and self.start < other.time():
                other = other + self.next_bday
            elif self.n < 0 and other.time() < self.start:
                other = other + self.next_bday
        return datetime(other.year, other.month, other.day,
                        self.start.hour, self.start.minute)

    def _prev_opening_time(self, other):
        """
        If n is positive, return yesterday's business day opening time.
        Otherwise yesterday business day's opening time.
        """
        if not self.next_bday.onOffset(other):
            other = other - self.next_bday
        else:
            if self.n >= 0 and other.time() < self.start:
                other = other - self.next_bday
            elif self.n < 0 and other.time() > self.start:
                other = other - self.next_bday
        return datetime(other.year, other.month, other.day,
                        self.start.hour, self.start.minute)

    def _get_business_hours_by_sec(self):
        """
        Return business hours in a day by seconds.
        """
        if self._get_daytime_flag():
            # create dummy datetime to calculate businesshours in a day
            dtstart = datetime(2014, 4, 1, self.start.hour, self.start.minute)
            until = datetime(2014, 4, 1, self.end.hour, self.end.minute)
            return (until - dtstart).total_seconds()
        else:
            self.daytime = False
            dtstart = datetime(2014, 4, 1, self.start.hour, self.start.minute)
            until = datetime(2014, 4, 2, self.end.hour, self.end.minute)
            return (until - dtstart).total_seconds()

    @apply_wraps
    def rollback(self, dt):
        """Roll provided date backward to next offset only if not on offset"""
        if not self.onOffset(dt):
            businesshours = self._get_business_hours_by_sec()
            if self.n >= 0:
                dt = self._prev_opening_time(
                    dt) + timedelta(seconds=businesshours)
            else:
                dt = self._next_opening_time(
                    dt) + timedelta(seconds=businesshours)
        return dt

    @apply_wraps
    def rollforward(self, dt):
        """Roll provided date forward to next offset only if not on offset"""
        if not self.onOffset(dt):
            if self.n >= 0:
                return self._next_opening_time(dt)
            else:
                return self._prev_opening_time(dt)
        return dt

    @apply_wraps
    def apply(self, other):
        # calculate here because offset is not immutable
        daytime = self._get_daytime_flag()
        businesshours = self._get_business_hours_by_sec()
        bhdelta = timedelta(seconds=businesshours)

        if isinstance(other, datetime):
            # used for detecting edge condition
            nanosecond = getattr(other, 'nanosecond', 0)
            # reset timezone and nanosecond
            # other may be a Timestamp, thus not use replace
            other = datetime(other.year, other.month, other.day,
                             other.hour, other.minute,
                             other.second, other.microsecond)
            n = self.n
            if n >= 0:
                if (other.time() == self.end or
                        not self._onOffset(other, businesshours)):
                    other = self._next_opening_time(other)
            else:
                if other.time() == self.start:
                    # adjustment to move to previous business day
                    other = other - timedelta(seconds=1)
                if not self._onOffset(other, businesshours):
                    other = self._next_opening_time(other)
                    other = other + bhdelta

            bd, r = divmod(abs(n * 60), businesshours // 60)
            if n < 0:
                bd, r = -bd, -r

            if bd != 0:
                skip_bd = BusinessDay(n=bd)
                # midnight business hour may not on BusinessDay
                if not self.next_bday.onOffset(other):
                    remain = other - self._prev_opening_time(other)
                    other = self._next_opening_time(other + skip_bd) + remain
                else:
                    other = other + skip_bd

            hours, minutes = divmod(r, 60)
            result = other + timedelta(hours=hours, minutes=minutes)

            # because of previous adjustment, time will be larger than start
            if ((daytime and (result.time() < self.start or
                              self.end < result.time())) or
                    not daytime and (self.end < result.time() < self.start)):
                if n >= 0:
                    bday_edge = self._prev_opening_time(other)
                    bday_edge = bday_edge + bhdelta
                    # calculate remainder
                    bday_remain = result - bday_edge
                    result = self._next_opening_time(other)
                    result += bday_remain
                else:
                    bday_edge = self._next_opening_time(other)
                    bday_remain = result - bday_edge
                    result = self._next_opening_time(result) + bhdelta
                    result += bday_remain
            # edge handling
            if n >= 0:
                if result.time() == self.end:
                    result = self._next_opening_time(result)
            else:
                if result.time() == self.start and nanosecond == 0:
                    # adjustment to move to previous business day
                    result = self._next_opening_time(
                        result - timedelta(seconds=1)) + bhdelta

            return result
        else:
            # TODO: Figure out the end of this sente
            raise ApplyTypeError(
                'Only know how to combine business hour with ')

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False

        if dt.tzinfo is not None:
            dt = datetime(dt.year, dt.month, dt.day, dt.hour,
                          dt.minute, dt.second, dt.microsecond)
        # Valid BH can be on the different BusinessDay during midnight
        # Distinguish by the time spent from previous opening time
        businesshours = self._get_business_hours_by_sec()
        return self._onOffset(dt, businesshours)

    def _onOffset(self, dt, businesshours):
        """
        Slight speedups using calculated values
        """
        # if self.normalize and not _is_normalized(dt):
        #     return False
        # Valid BH can be on the different BusinessDay during midnight
        # Distinguish by the time spent from previous opening time
        if self.n >= 0:
            op = self._prev_opening_time(dt)
        else:
            op = self._next_opening_time(dt)
        span = (dt - op).total_seconds()
        if span <= businesshours:
            return True
        else:
            return False

    def _repr_attrs(self):
        out = super(BusinessHourMixin, self)._repr_attrs()
        start = self.start.strftime('%H:%M')
        end = self.end.strftime('%H:%M')
        attrs = ['{prefix}={start}-{end}'.format(prefix=self._prefix,
                                                 start=start, end=end)]
        out += ': ' + ', '.join(attrs)
        return out


class BusinessHour(BusinessHourMixin, SingleConstructorOffset):
    """
    DateOffset subclass representing possibly n business days

    .. versionadded:: 0.16.1

    """
    _prefix = 'BH'
    _anchor = 0

    def __init__(self, n=1, normalize=False, start='09:00',
                 end='17:00', offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        super(BusinessHour, self).__init__(start=start, end=end, offset=offset)

    @cache_readonly
    def next_bday(self):
        # used for moving to next businessday
        if self.n >= 0:
            nb_offset = 1
        else:
            nb_offset = -1
        return BusinessDay(n=nb_offset)


class CustomBusinessDay(BusinessDay):
    """
    DateOffset subclass representing possibly n custom business days,
    excluding holidays

    Parameters
    ----------
    n : int, default 1
    offset : timedelta, default timedelta(0)
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``
    calendar : pd.HolidayCalendar or np.busdaycalendar
    """
    _cacheable = False
    _prefix = 'C'

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset
        self.kwds = {}

        calendar, holidays = _get_calendar(weekmask=weekmask,
                                           holidays=holidays,
                                           calendar=calendar)
        # CustomBusinessDay instances are identified by the
        # following two attributes. See DateOffset._params()
        # holidays, weekmask

        self.kwds['weekmask'] = self.weekmask = weekmask
        self.kwds['holidays'] = self.holidays = holidays
        self.kwds['calendar'] = self.calendar = calendar
        self.kwds['offset'] = offset

    @apply_wraps
    def apply(self, other):
        if self.n <= 0:
            roll = 'forward'
        else:
            roll = 'backward'

        if isinstance(other, datetime):
            date_in = other
            np_dt = np.datetime64(date_in.date())

            np_incr_dt = np.busday_offset(np_dt, self.n, roll=roll,
                                          busdaycal=self.calendar)

            dt_date = np_incr_dt.astype(datetime)
            result = datetime.combine(dt_date, date_in.time())

            if self.offset:
                result = result + self.offset
            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other,
                        normalize=self.normalize)
        else:
            raise ApplyTypeError('Only know how to combine trading day with '
                                 'datetime, datetime64 or timedelta.')

    def apply_index(self, i):
        raise NotImplementedError

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        day64 = _to_dt64(dt, 'datetime64[D]')
        return np.is_busday(day64, busdaycal=self.calendar)


class CustomBusinessHour(BusinessHourMixin, SingleConstructorOffset):
    """
    DateOffset subclass representing possibly n custom business days

    .. versionadded:: 0.18.1

    """
    _prefix = 'CBH'
    _anchor = 0

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None,
                 start='09:00', end='17:00', offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        super(CustomBusinessHour, self).__init__(start=start,
                                                 end=end, offset=offset)

        calendar, holidays = _get_calendar(weekmask=weekmask,
                                           holidays=holidays,
                                           calendar=calendar)
        self.kwds['weekmask'] = self.weekmask = weekmask
        self.kwds['holidays'] = self.holidays = holidays
        self.kwds['calendar'] = self.calendar = calendar

    @cache_readonly
    def next_bday(self):
        # used for moving to next businessday
        if self.n >= 0:
            nb_offset = 1
        else:
            nb_offset = -1
        return CustomBusinessDay(n=nb_offset,
                                 weekmask=self.weekmask,
                                 holidays=self.holidays,
                                 calendar=self.calendar)


# ---------------------------------------------------------------------
# Month-Based Offset Classes


class MonthOffset(SingleConstructorOffset):
    _adjust_dst = True

    def __init__(self, n=1, normalize=False):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.kwds = {}

    @property
    def name(self):
        if self.isAnchored:
            return self.rule_code
        else:
            month = liboffsets._int_to_month[self.n]
            return "{code}-{month}".format(code=self.rule_code,
                                           month=month)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.day == self._get_offset_day(dt)

    @apply_wraps
    def apply(self, other):
        n = self.n
        compare_day = self._get_offset_day(other)

        if n > 0 and other.day < compare_day:
            n -= 1
        elif n <= 0 and other.day > compare_day:
            # as if rolled forward already
            n += 1

        return shift_month(other, n, self._day_opt)

    @apply_index_wraps
    def apply_index(self, i):
        shifted = liboffsets.shift_months(i.asi8, self.n, self._day_opt)
        return i._shallow_copy(shifted)


class MonthEnd(MonthOffset):
    """DateOffset of one month end"""
    _prefix = 'M'
    _day_opt = 'end'


class MonthBegin(MonthOffset):
    """DateOffset of one month at beginning"""
    _prefix = 'MS'
    _day_opt = 'start'


class BusinessMonthEnd(MonthOffset):
    """DateOffset increments between business EOM dates"""
    _prefix = 'BM'
    _day_opt = 'business_end'


class BusinessMonthBegin(MonthOffset):
    """DateOffset of one business month at beginning"""
    _prefix = 'BMS'
    _day_opt = 'business_start'


class CustomBusinessMonthEnd(BusinessMixin, MonthOffset):
    """
    DateOffset subclass representing one custom business month, incrementing
    between end of month dates

    Parameters
    ----------
    n : int, default 1
    offset : timedelta, default timedelta(0)
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``
    calendar : pd.HolidayCalendar or np.busdaycalendar
    """

    _cacheable = False
    _prefix = 'CBM'

    onOffset = DateOffset.onOffset  # override MonthOffset method
    apply_index = DateOffset.apply_index  # override MonthOffset method

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset
        self.kwds = {}

        calendar, holidays = _get_calendar(weekmask=weekmask,
                                           holidays=holidays,
                                           calendar=calendar)
        self.kwds['weekmask'] = self.weekmask = weekmask
        self.kwds['holidays'] = self.holidays = holidays
        self.kwds['calendar'] = self.calendar = calendar
        self.kwds['offset'] = offset

    @cache_readonly
    def cbday(self):
        kwds = self.kwds
        return CustomBusinessDay(n=self.n, normalize=self.normalize, **kwds)

    @cache_readonly
    def m_offset(self):
        return MonthEnd(n=1, normalize=self.normalize)

    @apply_wraps
    def apply(self, other):
        n = self.n

        # First move to month offset
        cur_mend = self.m_offset.rollforward(other)

        # Find this custom month offset
        cur_cmend = self.cbday.rollback(cur_mend)

        # handle zero case. arbitrarily rollforward
        if n == 0 and other != cur_cmend:
            n += 1

        if other < cur_cmend and n >= 1:
            n -= 1
        elif other > cur_cmend and n <= -1:
            n += 1

        new = cur_mend + n * self.m_offset
        result = self.cbday.rollback(new)
        return result


class CustomBusinessMonthBegin(BusinessMixin, MonthOffset):
    """
    DateOffset subclass representing one custom business month, incrementing
    between beginning of month dates

    Parameters
    ----------
    n : int, default 1
    offset : timedelta, default timedelta(0)
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``
    calendar : pd.HolidayCalendar or np.busdaycalendar
    """

    _cacheable = False
    _prefix = 'CBMS'

    onOffset = DateOffset.onOffset  # override MonthOffset method
    apply_index = DateOffset.apply_index  # override MonthOffset method

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset
        self.kwds = {}

        # _get_calendar does validation and possible transformation
        # of calendar and holidays.
        calendar, holidays = _get_calendar(weekmask=weekmask,
                                           holidays=holidays,
                                           calendar=calendar)
        self.kwds['calendar'] = self.calendar = calendar
        self.kwds['weekmask'] = self.weekmask = weekmask
        self.kwds['holidays'] = self.holidays = holidays
        self.kwds['offset'] = offset

    @cache_readonly
    def cbday(self):
        kwds = self.kwds
        return CustomBusinessDay(n=self.n, normalize=self.normalize, **kwds)

    @cache_readonly
    def m_offset(self):
        return MonthBegin(n=1, normalize=self.normalize)

    @apply_wraps
    def apply(self, other):
        n = self.n
        dt_in = other

        # First move to month offset
        cur_mbegin = self.m_offset.rollback(dt_in)

        # Find this custom month offset
        cur_cmbegin = self.cbday.rollforward(cur_mbegin)

        # handle zero case. arbitrarily rollforward
        if n == 0 and dt_in != cur_cmbegin:
            n += 1

        if dt_in > cur_cmbegin and n <= -1:
            n += 1
        elif dt_in < cur_cmbegin and n >= 1:
            n -= 1

        new = cur_mbegin + n * self.m_offset
        result = self.cbday.rollforward(new)
        return result


# ---------------------------------------------------------------------
# Semi-Month Based Offset Classes

class SemiMonthOffset(DateOffset):
    _adjust_dst = True
    _default_day_of_month = 15
    _min_day_of_month = 2

    def __init__(self, n=1, normalize=False, day_of_month=None):
        if day_of_month is None:
            self.day_of_month = self._default_day_of_month
        else:
            self.day_of_month = int(day_of_month)
        if not self._min_day_of_month <= self.day_of_month <= 27:
            msg = 'day_of_month must be {min}<=day_of_month<=27, got {day}'
            raise ValueError(msg.format(min=self._min_day_of_month,
                                        day=self.day_of_month))

        self.n = self._validate_n(n)
        self.normalize = normalize
        self.kwds = {'day_of_month': self.day_of_month}

    @classmethod
    def _from_name(cls, suffix=None):
        return cls(day_of_month=suffix)

    @property
    def rule_code(self):
        suffix = '-{day_of_month}'.format(day_of_month=self.day_of_month)
        return self._prefix + suffix

    @apply_wraps
    def apply(self, other):
        n = self.n
        if not self.onOffset(other):
            _, days_in_month = tslib.monthrange(other.year, other.month)
            if 1 < other.day < self.day_of_month:
                other = other.replace(day=self.day_of_month)
                if n > 0:
                    # rollforward so subtract 1
                    n -= 1
            elif self.day_of_month < other.day < days_in_month:
                other = other.replace(day=self.day_of_month)
                if n < 0:
                    # rollforward in the negative direction so add 1
                    n += 1
                elif n == 0:
                    n = 1

        return self._apply(n, other)

    def _apply(self, n, other):
        """Handle specific apply logic for child classes"""
        raise AbstractMethodError(self)

    @apply_index_wraps
    def apply_index(self, i):
        # determine how many days away from the 1st of the month we are
        days_from_start = i.to_perioddelta('M').asi8
        delta = Timedelta(days=self.day_of_month - 1).value

        # get boolean array for each element before the day_of_month
        before_day_of_month = days_from_start < delta

        # get boolean array for each element after the day_of_month
        after_day_of_month = days_from_start > delta

        # determine the correct n for each date in i
        roll = self._get_roll(i, before_day_of_month, after_day_of_month)

        # isolate the time since it will be striped away one the next line
        time = i.to_perioddelta('D')

        # apply the correct number of months
        i = (i.to_period('M') + (roll // 2)).to_timestamp()

        # apply the correct day
        i = self._apply_index_days(i, roll)

        return i + time

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        """Return an array with the correct n for each date in i.

        The roll array is based on the fact that i gets rolled back to
        the first day of the month.
        """
        raise AbstractMethodError(self)

    def _apply_index_days(self, i, roll):
        """Apply the correct day for each date in i"""
        raise AbstractMethodError(self)


class SemiMonthEnd(SemiMonthOffset):
    """
    Two DateOffset's per month repeating on the last
    day of the month and day_of_month.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    n: int
    normalize : bool, default False
    day_of_month: int, {1, 3,...,27}, default 15
    """
    _prefix = 'SM'
    _min_day_of_month = 1

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        _, days_in_month = tslib.monthrange(dt.year, dt.month)
        return dt.day in (self.day_of_month, days_in_month)

    def _apply(self, n, other):
        # if other.day is not day_of_month move to day_of_month and update n
        if n > 0 and other.day < self.day_of_month:
            n -= 1
        elif other.day > self.day_of_month:
            n += 1

        months = n // 2
        day = 31 if n % 2 else self.day_of_month
        return shift_month(other, months, day)

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        n = self.n
        is_month_end = i.is_month_end
        if n > 0:
            roll_end = np.where(is_month_end, 1, 0)
            roll_before = np.where(before_day_of_month, n, n + 1)
            roll = roll_end + roll_before
        elif n == 0:
            roll_after = np.where(after_day_of_month, 2, 0)
            roll_before = np.where(~after_day_of_month, 1, 0)
            roll = roll_before + roll_after
        else:
            roll = np.where(after_day_of_month, n + 2, n + 1)
        return roll

    def _apply_index_days(self, i, roll):
        i += (roll % 2) * Timedelta(days=self.day_of_month).value
        return i + Timedelta(days=-1)


class SemiMonthBegin(SemiMonthOffset):
    """
    Two DateOffset's per month repeating on the first
    day of the month and day_of_month.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    n: int
    normalize : bool, default False
    day_of_month: int, {2, 3,...,27}, default 15
    """
    _prefix = 'SMS'

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.day in (1, self.day_of_month)

    def _apply(self, n, other):
        # if other.day is not day_of_month move to day_of_month and update n
        if other.day < self.day_of_month:
            n -= 1
        elif n <= 0 and other.day > self.day_of_month:
            n += 1

        months = n // 2 + n % 2
        day = 1 if n % 2 else self.day_of_month
        return shift_month(other, months, day)

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        n = self.n
        is_month_start = i.is_month_start
        if n > 0:
            roll = np.where(before_day_of_month, n, n + 1)
        elif n == 0:
            roll_start = np.where(is_month_start, 0, 1)
            roll_after = np.where(after_day_of_month, 1, 0)
            roll = roll_start + roll_after
        else:
            roll_after = np.where(after_day_of_month, n + 2, n + 1)
            roll_start = np.where(is_month_start, -1, 0)
            roll = roll_after + roll_start
        return roll

    def _apply_index_days(self, i, roll):
        return i + (roll % 2) * Timedelta(days=self.day_of_month - 1).value


# ---------------------------------------------------------------------
# Week-Based Offset Classes

class Week(EndMixin, DateOffset):
    """
    Weekly offset

    Parameters
    ----------
    weekday : int, default None
        Always generate specific day of week. 0 for Monday
    """
    _adjust_dst = True
    _inc = timedelta(weeks=1)
    _prefix = 'W'

    def __init__(self, n=1, normalize=False, weekday=None):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.weekday = weekday

        if self.weekday is not None:
            if self.weekday < 0 or self.weekday > 6:
                raise ValueError('Day must be 0<=day<=6, got {day}'
                                 .format(day=self.weekday))

        self.kwds = {'weekday': weekday}

    def isAnchored(self):
        return (self.n == 1 and self.weekday is not None)

    @apply_wraps
    def apply(self, other):
        if self.weekday is None:
            return other + self.n * self._inc

        k = self.n
        otherDay = other.weekday()
        if otherDay != self.weekday:
            other = other + timedelta((self.weekday - otherDay) % 7)
            if k > 0:
                k -= 1

        return other + timedelta(weeks=k)

    @apply_index_wraps
    def apply_index(self, i):
        if self.weekday is None:
            return ((i.to_period('W') + self.n).to_timestamp() +
                    i.to_perioddelta('W'))
        else:
            return self._end_apply_index(i, self.freqstr)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.weekday() == self.weekday

    @property
    def rule_code(self):
        suffix = ''
        if self.weekday is not None:
            suffix = '-{weekday}'.format(weekday=_int_to_weekday[self.weekday])
        return self._prefix + suffix

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            weekday = None
        else:
            weekday = _weekday_to_int[suffix]
        return cls(weekday=weekday)


class WeekOfMonth(DateOffset):
    """
    Describes monthly dates like "the Tuesday of the 2nd week of each month"

    Parameters
    ----------
    n : int
    week : {0, 1, 2, 3, ...}, default None
        0 is 1st week of month, 1 2nd week, etc.
    weekday : {0, 1, ..., 6}, default None
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    """
    _prefix = 'WOM'
    _adjust_dst = True

    def __init__(self, n=1, normalize=False, week=None, weekday=None):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.weekday = weekday
        self.week = week

        if self.n == 0:
            raise ValueError('N cannot be 0')

        if self.weekday < 0 or self.weekday > 6:
            raise ValueError('Day must be 0<=day<=6, got {day}'
                             .format(day=self.weekday))
        if self.week < 0 or self.week > 3:
            raise ValueError('Week must be 0<=week<=3, got {week}'
                             .format(week=self.week))

        self.kwds = {'weekday': weekday, 'week': week}

    @apply_wraps
    def apply(self, other):
        base = other
        offsetOfMonth = self.getOffsetOfMonth(other)

        months = self.n
        if months > 0 and offsetOfMonth > other:
            months -= 1
        elif months <= 0 and offsetOfMonth < other:
            months += 1

        other = self.getOffsetOfMonth(shift_month(other, months, 'start'))
        other = datetime(other.year, other.month, other.day, base.hour,
                         base.minute, base.second, base.microsecond)
        return other

    def getOffsetOfMonth(self, dt):
        w = Week(weekday=self.weekday)
        d = datetime(dt.year, dt.month, 1, tzinfo=dt.tzinfo)
        # TODO: Is this DST-safe?
        d = w.rollforward(d)
        return d + timedelta(weeks=self.week)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        d = datetime(dt.year, dt.month, dt.day, tzinfo=dt.tzinfo)
        return d == self.getOffsetOfMonth(dt)

    @property
    def rule_code(self):
        weekday = _int_to_weekday.get(self.weekday, '')
        return '{prefix}-{week}{weekday}'.format(prefix=self._prefix,
                                                 week=self.week + 1,
                                                 weekday=weekday)

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            raise ValueError("Prefix {prefix!r} requires a suffix."
                             .format(prefix=cls._prefix))
        # TODO: handle n here...
        # only one digit weeks (1 --> week 0, 2 --> week 1, etc.)
        week = int(suffix[0]) - 1
        weekday = _weekday_to_int[suffix[1:]]
        return cls(week=week, weekday=weekday)


class LastWeekOfMonth(DateOffset):
    """
    Describes monthly dates in last week of month like "the last Tuesday of
    each month"

    Parameters
    ----------
    n : int, default 1
    weekday : {0, 1, ..., 6}, default None
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays

    """
    _prefix = 'LWOM'

    def __init__(self, n=1, normalize=False, weekday=None):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.weekday = weekday

        if self.n == 0:
            raise ValueError('N cannot be 0')

        if self.weekday < 0 or self.weekday > 6:
            raise ValueError('Day must be 0<=day<=6, got {day}'
                             .format(day=self.weekday))

        self.kwds = {'weekday': weekday}

    @apply_wraps
    def apply(self, other):
        offsetOfMonth = self.getOffsetOfMonth(other)

        months = self.n
        if months > 0 and offsetOfMonth > other:
            months -= 1
        elif months <= 0 and offsetOfMonth < other:
            months += 1

        return self.getOffsetOfMonth(shift_month(other, months, 'start'))

    def getOffsetOfMonth(self, dt):
        m = MonthEnd()
        d = datetime(dt.year, dt.month, 1, dt.hour, dt.minute,
                     dt.second, dt.microsecond, tzinfo=dt.tzinfo)
        eom = m.rollforward(d)
        # TODO: Is this DST-safe?
        w = Week(weekday=self.weekday)
        return w.rollback(eom)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt == self.getOffsetOfMonth(dt)

    @property
    def rule_code(self):
        weekday = _int_to_weekday.get(self.weekday, '')
        return '{prefix}-{weekday}'.format(prefix=self._prefix,
                                           weekday=weekday)

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            raise ValueError("Prefix {prefix!r} requires a suffix."
                             .format(prefix=cls._prefix))
        # TODO: handle n here...
        weekday = _weekday_to_int[suffix]
        return cls(weekday=weekday)

# ---------------------------------------------------------------------
# Quarter-Based Offset Classes


class QuarterOffset(DateOffset):
    """Quarter representation - doesn't call super"""
    _default_startingMonth = None
    _from_name_startingMonth = None
    _adjust_dst = True
    # TODO: Consider combining QuarterOffset and YearOffset __init__ at some
    #       point

    def __init__(self, n=1, normalize=False, startingMonth=None):
        self.n = self._validate_n(n)
        self.normalize = normalize
        if startingMonth is None:
            startingMonth = self._default_startingMonth
        self.startingMonth = startingMonth

        self.kwds = {'startingMonth': startingMonth}

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    @classmethod
    def _from_name(cls, suffix=None):
        kwargs = {}
        if suffix:
            kwargs['startingMonth'] = liboffsets._month_to_int[suffix]
        else:
            if cls._from_name_startingMonth is not None:
                kwargs['startingMonth'] = cls._from_name_startingMonth
        return cls(**kwargs)

    @property
    def rule_code(self):
        month = liboffsets._int_to_month[self.startingMonth]
        return '{prefix}-{month}'.format(prefix=self._prefix, month=month)

    @apply_wraps
    def apply(self, other):
        n = self.n
        compare_day = self._get_offset_day(other)

        months_since = (other.month - self.startingMonth) % 3

        if n <= 0 and (months_since != 0 or
                       (months_since == 0 and other.day > compare_day)):
            # make sure to roll forward, so negate
            n += 1
        elif n > 0 and (months_since == 0 and other.day < compare_day):
            # pretend to roll back if on same month but before compare_day
            n -= 1

        return shift_month(other, 3 * n - months_since, self._day_opt)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        modMonth = (dt.month - self.startingMonth) % 3
        return modMonth == 0 and dt.day == self._get_offset_day(dt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = liboffsets.shift_quarters(dtindex.asi8, self.n,
                                            self.startingMonth, self._day_opt)
        return dtindex._shallow_copy(shifted)


class BQuarterEnd(QuarterOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...
    """
    _outputName = 'BusinessQuarterEnd'
    _default_startingMonth = 3
    _from_name_startingMonth = 12
    _prefix = 'BQ'
    _day_opt = 'business_end'


# TODO: This is basically the same as BQuarterEnd
class BQuarterBegin(QuarterOffset):
    _outputName = "BusinessQuarterBegin"
    # I suspect this is wrong for *all* of them.
    _default_startingMonth = 3
    _from_name_startingMonth = 1
    _prefix = 'BQS'
    _day_opt = 'business_start'


class QuarterEnd(QuarterOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/31/2007, 6/30/2007, ...
    """
    _outputName = 'QuarterEnd'
    _default_startingMonth = 3
    _prefix = 'Q'
    _day_opt = 'end'


class QuarterBegin(QuarterOffset):
    _outputName = 'QuarterBegin'
    _default_startingMonth = 3
    _from_name_startingMonth = 1
    _prefix = 'QS'
    _day_opt = 'start'


# ---------------------------------------------------------------------
# Year-Based Offset Classes

class YearOffset(DateOffset):
    """DateOffset that just needs a month"""
    _adjust_dst = True

    def _get_offset_day(self, other):
        # override BaseOffset method to use self.month instead of other.month
        # TODO: there may be a more performant way to do this
        return liboffsets.get_day_of_month(other.replace(month=self.month),
                                           self._day_opt)

    @apply_wraps
    def apply(self, other):
        years = roll_yearday(other, self.n, self.month, self._day_opt)
        months = years * 12 + (self.month - other.month)
        return shift_month(other, months, self._day_opt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = liboffsets.shift_quarters(dtindex.asi8, self.n,
                                            self.month, self._day_opt,
                                            modby=12)
        return dtindex._shallow_copy(shifted)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.month == self.month and dt.day == self._get_offset_day(dt)

    def __init__(self, n=1, normalize=False, month=None):
        month = month if month is not None else self._default_month
        self.month = month

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, normalize=normalize, month=month)

    @classmethod
    def _from_name(cls, suffix=None):
        kwargs = {}
        if suffix:
            kwargs['month'] = liboffsets._month_to_int[suffix]
        return cls(**kwargs)

    @property
    def rule_code(self):
        month = liboffsets._int_to_month[self.month]
        return '{prefix}-{month}'.format(prefix=self._prefix, month=month)


class BYearEnd(YearOffset):
    """DateOffset increments between business EOM dates"""
    _outputName = 'BusinessYearEnd'
    _default_month = 12
    _prefix = 'BA'
    _day_opt = 'business_end'


class BYearBegin(YearOffset):
    """DateOffset increments between business year begin dates"""
    _outputName = 'BusinessYearBegin'
    _default_month = 1
    _prefix = 'BAS'
    _day_opt = 'business_start'


class YearEnd(YearOffset):
    """DateOffset increments between calendar year ends"""
    _default_month = 12
    _prefix = 'A'
    _day_opt = 'end'


class YearBegin(YearOffset):
    """DateOffset increments between calendar year begin dates"""
    _default_month = 1
    _prefix = 'AS'
    _day_opt = 'start'


# ---------------------------------------------------------------------
# Special Offset Classes

class FY5253(DateOffset):
    """
    Describes 52-53 week fiscal year. This is also known as a 4-4-5 calendar.

    It is used by companies that desire that their
    fiscal year always end on the same day of the week.

    It is a method of managing accounting periods.
    It is a common calendar structure for some industries,
    such as retail, manufacturing and parking industry.

    For more information see:
    http://en.wikipedia.org/wiki/4%E2%80%934%E2%80%935_calendar


    The year may either:
    - end on the last X day of the Y month.
    - end on the last X day closest to the last day of the Y month.

    X is a specific day of the week.
    Y is a certain month of the year

    Parameters
    ----------
    n : int
    weekday : {0, 1, ..., 6}
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    startingMonth : The month in which fiscal years end. {1, 2, ... 12}
    variation : str
        {"nearest", "last"} for "LastOfMonth" or "NearestEndMonth"
    """
    _prefix = 'RE'
    _adjust_dst = True

    def __init__(self, n=1, normalize=False, weekday=0, startingMonth=1,
                 variation="nearest"):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.startingMonth = startingMonth
        self.weekday = weekday

        self.variation = variation

        self.kwds = {'weekday': weekday, 'startingMonth': startingMonth,
                     'variation': variation}

        if self.n == 0:
            raise ValueError('N cannot be 0')

        if self.variation not in ["nearest", "last"]:
            raise ValueError('{variation} is not a valid variation'
                             .format(variation=self.variation))

    @cache_readonly
    def _offset_lwom(self):
        if self.variation == "nearest":
            return None
        else:
            return LastWeekOfMonth(n=1, weekday=self.weekday)

    def isAnchored(self):
        return (self.n == 1 and
                self.startingMonth is not None and
                self.weekday is not None)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        dt = datetime(dt.year, dt.month, dt.day)
        year_end = self.get_year_end(dt)

        if self.variation == "nearest":
            # We have to check the year end of "this" cal year AND the previous
            return (year_end == dt or
                    self.get_year_end(shift_month(dt, -1, None)) == dt)
        else:
            return year_end == dt

    @apply_wraps
    def apply(self, other):
        n = self.n
        prev_year = self.get_year_end(
            datetime(other.year - 1, self.startingMonth, 1))
        cur_year = self.get_year_end(
            datetime(other.year, self.startingMonth, 1))
        next_year = self.get_year_end(
            datetime(other.year + 1, self.startingMonth, 1))

        prev_year = tslib._localize_pydatetime(prev_year, other.tzinfo)
        cur_year = tslib._localize_pydatetime(cur_year, other.tzinfo)
        next_year = tslib._localize_pydatetime(next_year, other.tzinfo)

        if other == prev_year:
            n -= 1
        elif other == cur_year:
            pass
        elif other == next_year:
            n += 1
            # TODO: Not hit in tests
        elif n > 0:
            if other < prev_year:
                n -= 2
                # TODO: Not hit in tests
            elif other < cur_year:
                n -= 1
            elif other < next_year:
                pass
            else:
                assert False
        else:
            if other > next_year:
                n += 2
                # TODO: Not hit in tests
            elif other > cur_year:
                n += 1
            elif other > prev_year:
                pass
            else:
                assert False

        shifted = datetime(other.year + n, self.startingMonth, 1)
        result = self.get_year_end(shifted)
        result = datetime(result.year, result.month, result.day,
                          other.hour, other.minute, other.second,
                          other.microsecond)
        return result

    def get_year_end(self, dt):
        if self.variation == "nearest":
            return self._get_year_end_nearest(dt)
        else:
            return self._get_year_end_last(dt)

    def get_target_month_end(self, dt):
        target_month = datetime(dt.year, self.startingMonth, 1,
                                tzinfo=dt.tzinfo)
        return shift_month(target_month, 0, 'end')
        # TODO: is this DST-safe?

    def _get_year_end_nearest(self, dt):
        target_date = self.get_target_month_end(dt)
        wkday_diff = self.weekday - target_date.weekday()
        if wkday_diff == 0:
            return target_date

        days_forward = wkday_diff % 7
        if days_forward <= 3:
            # The upcoming self.weekday is closer than the previous one
            return target_date + timedelta(days_forward)
        else:
            # The previous self.weekday is closer than the upcoming one
            return target_date + timedelta(days_forward - 7)

    def _get_year_end_last(self, dt):
        current_year = datetime(dt.year, self.startingMonth, 1,
                                tzinfo=dt.tzinfo)
        return current_year + self._offset_lwom

    @property
    def rule_code(self):
        prefix = self._prefix
        suffix = self.get_rule_code_suffix()
        return "{prefix}-{suffix}".format(prefix=prefix, suffix=suffix)

    def _get_suffix_prefix(self):
        if self.variation == "nearest":
            return 'N'
        else:
            return 'L'

    def get_rule_code_suffix(self):
        prefix = self._get_suffix_prefix()
        month = liboffsets._int_to_month[self.startingMonth]
        weekday = _int_to_weekday[self.weekday]
        return '{prefix}-{month}-{weekday}'.format(prefix=prefix, month=month,
                                                   weekday=weekday)

    @classmethod
    def _parse_suffix(cls, varion_code, startingMonth_code, weekday_code):
        if varion_code == "N":
            variation = "nearest"
        elif varion_code == "L":
            variation = "last"
        else:
            raise ValueError("Unable to parse varion_code: "
                             "{code}".format(code=varion_code))

        startingMonth = liboffsets._month_to_int[startingMonth_code]
        weekday = _weekday_to_int[weekday_code]

        return {"weekday": weekday,
                "startingMonth": startingMonth,
                "variation": variation}

    @classmethod
    def _from_name(cls, *args):
        return cls(**cls._parse_suffix(*args))


class FY5253Quarter(DateOffset):
    """
    DateOffset increments between business quarter dates
    for 52-53 week fiscal year (also known as a 4-4-5 calendar).

    It is used by companies that desire that their
    fiscal year always end on the same day of the week.

    It is a method of managing accounting periods.
    It is a common calendar structure for some industries,
    such as retail, manufacturing and parking industry.

    For more information see:
    http://en.wikipedia.org/wiki/4%E2%80%934%E2%80%935_calendar

    The year may either:
    - end on the last X day of the Y month.
    - end on the last X day closest to the last day of the Y month.

    X is a specific day of the week.
    Y is a certain month of the year

    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...

    Parameters
    ----------
    n : int
    weekday : {0, 1, ..., 6}
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    startingMonth : The month in which fiscal years end. {1, 2, ... 12}
    qtr_with_extra_week : The quarter number that has the leap
        or 14 week when needed. {1,2,3,4}
    variation : str
        {"nearest", "last"} for "LastOfMonth" or "NearestEndMonth"
    """

    _prefix = 'REQ'
    _adjust_dst = True

    def __init__(self, n=1, normalize=False, weekday=0, startingMonth=1,
                 qtr_with_extra_week=1, variation="nearest"):
        self.n = self._validate_n(n)
        self.normalize = normalize

        self.weekday = weekday
        self.startingMonth = startingMonth
        self.qtr_with_extra_week = qtr_with_extra_week
        self.variation = variation

        self.kwds = {'weekday': weekday, 'startingMonth': startingMonth,
                     'qtr_with_extra_week': qtr_with_extra_week,
                     'variation': variation}

        if self.n == 0:
            raise ValueError('N cannot be 0')

    @cache_readonly
    def _offset(self):
        return FY5253(startingMonth=self.startingMonth,
                      weekday=self.weekday,
                      variation=self.variation)

    def isAnchored(self):
        return self.n == 1 and self._offset.isAnchored()

    @apply_wraps
    def apply(self, other):
        base = other
        n = self.n

        if n > 0:
            while n > 0:
                if not self._offset.onOffset(other):
                    qtr_lens = self.get_weeks(other)
                    start = other - self._offset
                else:
                    start = other
                    qtr_lens = self.get_weeks(other + self._offset)

                for weeks in qtr_lens:
                    start += timedelta(weeks=weeks)
                    if start > other:
                        other = start
                        n -= 1
                        break

        else:
            n = -n
            while n > 0:
                if not self._offset.onOffset(other):
                    qtr_lens = self.get_weeks(other)
                    end = other + self._offset
                else:
                    end = other
                    qtr_lens = self.get_weeks(other)

                for weeks in reversed(qtr_lens):
                    end -= timedelta(weeks=weeks)
                    if end < other:
                        other = end
                        n -= 1
                        break
        other = datetime(other.year, other.month, other.day,
                         base.hour, base.minute, base.second, base.microsecond)
        return other

    def get_weeks(self, dt):
        ret = [13] * 4

        year_has_extra_week = self.year_has_extra_week(dt)

        if year_has_extra_week:
            ret[self.qtr_with_extra_week - 1] = 14

        return ret

    def year_has_extra_week(self, dt):
        if self._offset.onOffset(dt):
            prev_year_end = dt - self._offset
            next_year_end = dt
        else:
            next_year_end = dt + self._offset
            prev_year_end = dt - self._offset

        week_in_year = (next_year_end - prev_year_end).days / 7

        return week_in_year == 53

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        if self._offset.onOffset(dt):
            return True

        next_year_end = dt - self._offset

        qtr_lens = self.get_weeks(dt)

        current = next_year_end
        for qtr_len in qtr_lens[0:4]:
            current += timedelta(weeks=qtr_len)
            if dt == current:
                return True
        return False

    @property
    def rule_code(self):
        suffix = self._offset.get_rule_code_suffix()
        qtr = self.qtr_with_extra_week
        return "{prefix}-{suffix}-{qtr}".format(prefix=self._prefix,
                                                suffix=suffix, qtr=qtr)

    @classmethod
    def _from_name(cls, *args):
        return cls(**dict(FY5253._parse_suffix(*args[:-1]),
                          qtr_with_extra_week=int(args[-1])))


class Easter(DateOffset):
    """
    DateOffset for the Easter holiday using
    logic defined in dateutil.  Right now uses
    the revised method which is valid in years
    1583-4099.
    """
    _adjust_dst = True

    def __init__(self, n=1, normalize=False):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.kwds = {}

    @apply_wraps
    def apply(self, other):
        current_easter = easter(other.year)
        current_easter = datetime(current_easter.year,
                                  current_easter.month, current_easter.day)
        current_easter = tslib._localize_pydatetime(current_easter,
                                                    other.tzinfo)

        n = self.n
        if n >= 0 and other < current_easter:
            n -= 1
        elif n < 0 and other > current_easter:
            n += 1

        # NOTE: easter returns a datetime.date so we have to convert to type of
        # other
        new = easter(other.year + n)
        new = datetime(new.year, new.month, new.day, other.hour,
                       other.minute, other.second, other.microsecond)
        return new

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return date(dt.year, dt.month, dt.day) == easter(dt.year)

# ---------------------------------------------------------------------
# Ticks


def _tick_comp(op):
    def f(self, other):
        return op(self.delta, other.delta)

    return f


class Tick(SingleConstructorOffset):
    _inc = Timedelta(microseconds=1000)
    _prefix = 'undefined'

    def __init__(self, n=1, normalize=False):
        # TODO: do Tick classes with normalize=True make sense?
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.kwds = {}

    __gt__ = _tick_comp(operator.gt)
    __ge__ = _tick_comp(operator.ge)
    __lt__ = _tick_comp(operator.lt)
    __le__ = _tick_comp(operator.le)
    __eq__ = _tick_comp(operator.eq)
    __ne__ = _tick_comp(operator.ne)

    def __add__(self, other):
        if isinstance(other, Tick):
            if type(self) == type(other):
                return type(self)(self.n + other.n)
            else:
                return _delta_to_tick(self.delta + other.delta)
        elif isinstance(other, ABCPeriod):
            return other + self
        try:
            return self.apply(other)
        except ApplyTypeError:
            return NotImplemented
        except OverflowError:
            raise OverflowError("the add operation between {self} and {other} "
                                "will overflow".format(self=self, other=other))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            from pandas.tseries.frequencies import to_offset

            other = to_offset(other)

        if isinstance(other, Tick):
            return self.delta == other.delta
        else:
            # TODO: Are there cases where this should raise TypeError?
            return False

    # This is identical to DateOffset.__hash__, but has to be redefined here
    # for Python 3, because we've redefined __eq__.
    def __hash__(self):
        return hash(self._params())

    def __ne__(self, other):
        if isinstance(other, compat.string_types):
            from pandas.tseries.frequencies import to_offset

            other = to_offset(other)

        if isinstance(other, Tick):
            return self.delta != other.delta
        else:
            # TODO: Are there cases where this should raise TypeError?
            return True

    @property
    def delta(self):
        return self.n * self._inc

    @property
    def nanos(self):
        return delta_to_nanoseconds(self.delta)

    # TODO: Should Tick have its own apply_index?
    def apply(self, other):
        # Timestamp can handle tz and nano sec, thus no need to use apply_wraps
        if isinstance(other, Timestamp):

            # GH 15126
            # in order to avoid a recursive
            # call of __add__ and __radd__ if there is
            # an exception, when we call using the + operator,
            # we directly call the known method
            result = other.__add__(self)
            if result == NotImplemented:
                raise OverflowError
            return result
        elif isinstance(other, (datetime, np.datetime64, date)):
            return as_timestamp(other) + self

        if isinstance(other, timedelta):
            return other + self.delta
        elif isinstance(other, type(self)):
            return type(self)(self.n + other.n)

        raise ApplyTypeError('Unhandled type: {type_str}'
                             .format(type_str=type(other).__name__))

    def isAnchored(self):
        return False


def _delta_to_tick(delta):
    if delta.microseconds == 0:
        if delta.seconds == 0:
            return Day(delta.days)
        else:
            seconds = delta.days * 86400 + delta.seconds
            if seconds % 3600 == 0:
                return Hour(seconds / 3600)
            elif seconds % 60 == 0:
                return Minute(seconds / 60)
            else:
                return Second(seconds)
    else:
        nanos = delta_to_nanoseconds(delta)
        if nanos % 1000000 == 0:
            return Milli(nanos // 1000000)
        elif nanos % 1000 == 0:
            return Micro(nanos // 1000)
        else:  # pragma: no cover
            return Nano(nanos)


class Day(Tick):
    _inc = Timedelta(days=1)
    _prefix = 'D'


class Hour(Tick):
    _inc = Timedelta(hours=1)
    _prefix = 'H'


class Minute(Tick):
    _inc = Timedelta(minutes=1)
    _prefix = 'T'


class Second(Tick):
    _inc = Timedelta(seconds=1)
    _prefix = 'S'


class Milli(Tick):
    _inc = Timedelta(milliseconds=1)
    _prefix = 'L'


class Micro(Tick):
    _inc = Timedelta(microseconds=1)
    _prefix = 'U'


class Nano(Tick):
    _inc = Timedelta(nanoseconds=1)
    _prefix = 'N'


BDay = BusinessDay
BMonthEnd = BusinessMonthEnd
BMonthBegin = BusinessMonthBegin
CBMonthEnd = CustomBusinessMonthEnd
CBMonthBegin = CustomBusinessMonthBegin
CDay = CustomBusinessDay

# ---------------------------------------------------------------------


def generate_range(start=None, end=None, periods=None,
                   offset=BDay(), time_rule=None):
    """
    Generates a sequence of dates corresponding to the specified time
    offset. Similar to dateutil.rrule except uses pandas DateOffset
    objects to represent time increments

    Parameters
    ----------
    start : datetime (default None)
    end : datetime (default None)
    periods : int, optional
    time_rule : (legacy) name of DateOffset object to be used, optional
        Corresponds with names expected by tseries.frequencies.get_offset

    Notes
    -----
    * This method is faster for generating weekdays than dateutil.rrule
    * At least two of (start, end, periods) must be specified.
    * If both start and end are specified, the returned dates will
    satisfy start <= date <= end.
    * If both time_rule and offset are specified, time_rule supersedes offset.

    Returns
    -------
    dates : generator object

    """
    if time_rule is not None:
        from pandas.tseries.frequencies import get_offset

        offset = get_offset(time_rule)

    start = to_datetime(start)
    end = to_datetime(end)

    if start and not offset.onOffset(start):
        start = offset.rollforward(start)

    elif end and not offset.onOffset(end):
        end = offset.rollback(end)

    if periods is None and end < start:
        end = None
        periods = 0

    if end is None:
        end = start + (periods - 1) * offset

    if start is None:
        start = end - (periods - 1) * offset

    cur = start
    if offset.n >= 0:
        while cur <= end:
            yield cur

            # faster than cur + offset
            next_date = offset.apply(cur)
            if next_date <= cur:
                raise ValueError('Offset {offset} did not increment date'
                                 .format(offset=offset))
            cur = next_date
    else:
        while cur >= end:
            yield cur

            # faster than cur + offset
            next_date = offset.apply(cur)
            if next_date >= cur:
                raise ValueError('Offset {offset} did not decrement date'
                                 .format(offset=offset))
            cur = next_date


prefix_mapping = dict((offset._prefix, offset) for offset in [
    YearBegin,                 # 'AS'
    YearEnd,                   # 'A'
    BYearBegin,                # 'BAS'
    BYearEnd,                  # 'BA'
    BusinessDay,               # 'B'
    BusinessMonthBegin,        # 'BMS'
    BusinessMonthEnd,          # 'BM'
    BQuarterEnd,               # 'BQ'
    BQuarterBegin,             # 'BQS'
    BusinessHour,              # 'BH'
    CustomBusinessDay,         # 'C'
    CustomBusinessMonthEnd,    # 'CBM'
    CustomBusinessMonthBegin,  # 'CBMS'
    CustomBusinessHour,        # 'CBH'
    MonthEnd,                  # 'M'
    MonthBegin,                # 'MS'
    Nano,                      # 'N'
    SemiMonthEnd,              # 'SM'
    SemiMonthBegin,            # 'SMS'
    Week,                      # 'W'
    Second,                    # 'S'
    Minute,                    # 'T'
    Micro,                     # 'U'
    QuarterEnd,                # 'Q'
    QuarterBegin,              # 'QS'
    Milli,                     # 'L'
    Hour,                      # 'H'
    Day,                       # 'D'
    WeekOfMonth,               # 'WOM'
    FY5253,
    FY5253Quarter,
])
