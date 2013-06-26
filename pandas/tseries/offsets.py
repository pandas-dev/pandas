from datetime import date, datetime, timedelta
import numpy as np

from pandas.tseries.tools import to_datetime

# import after tools, dateutil check
from dateutil.relativedelta import relativedelta
import pandas.tslib as tslib

__all__ = ['Day', 'BusinessDay', 'BDay', 'CustomBusinessDay', 'CDay',
           'MonthBegin', 'BMonthBegin', 'MonthEnd', 'BMonthEnd',
           'YearBegin', 'BYearBegin', 'YearEnd', 'BYearEnd',
           'QuarterBegin', 'BQuarterBegin', 'QuarterEnd', 'BQuarterEnd',
           'Week', 'WeekOfMonth',
           'Hour', 'Minute', 'Second', 'Milli', 'Micro', 'Nano']

#----------------------------------------------------------------------
# DateOffset


class CacheableOffset(object):

    _cacheable = True


class DateOffset(object):
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
    valid date, first it is moved to a valid date.  Thus psedo code
    is:

    def __add__(date):
      date = rollback(date) # does nothing if date is valid
      return date + <n number of periods>

    When a date offset is created for a negitive number of periods,
    the date is first rolled forward.  The pseudo code is:

    def __add__(date):
      date = rollforward(date) # does nothing is date is valid
      return date + <n number of periods>

    Zero presents a problem.  Should it roll forward or back?  We
    arbitrarily have it rollforward:

    date + BDay(0) == BDay.rollforward(date)

    Since 0 is a bit weird, we suggest avoiding its use.
    """
    _cacheable = False
    _normalize_cache = True

    def __init__(self, n=1, **kwds):
        self.n = int(n)
        self.kwds = kwds
        if len(kwds) > 0:
            self._offset = relativedelta(**kwds)
        else:
            self._offset = timedelta(1)

    def apply(self, other):
        if len(self.kwds) > 0:
            if self.n > 0:
                for i in xrange(self.n):
                    other = other + self._offset
            else:
                for i in xrange(-self.n):
                    other = other - self._offset
            return other
        else:
            return other + timedelta(self.n)

    def isAnchored(self):
        return (self.n == 1)

    def copy(self):
        return self.__class__(self.n, **self.kwds)

    def _should_cache(self):
        return self.isAnchored() and self._cacheable

    def _params(self):
        attrs = [(k, v) for k, v in vars(self).iteritems()
                 if k not in ['kwds', '_offset', 'name', 'normalize',
                 'busdaycalendar']]
        attrs.extend(self.kwds.items())
        attrs = sorted(set(attrs))

        params = tuple([str(self.__class__)] + attrs)
        return params

    def __repr__(self):
        if hasattr(self, 'name') and len(self.name):
            return self.name

        className = getattr(self, '_outputName', type(self).__name__)
        exclude = set(['n', 'inc'])
        attrs = []
        for attr in self.__dict__:
            if ((attr == 'kwds' and len(self.kwds) == 0)
                    or attr.startswith('_')):
                continue
            if attr not in exclude:
                attrs.append('='.join((attr, repr(getattr(self, attr)))))

        if abs(self.n) != 1:
            plural = 's'
        else:
            plural = ''

        out = '<%s ' % self.n + className + plural
        if attrs:
            out += ': ' + ', '.join(attrs)
        out += '>'
        return out

    def __eq__(self, other):
        if other is None:
            return False

        if isinstance(other, basestring):
            from pandas.tseries.frequencies import to_offset
            other = to_offset(other)

        if not isinstance(other, DateOffset):
            return False

        return self._params() == other._params()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._params())

    def __call__(self, other):
        return self.apply(other)

    def __add__(self, other):
        return self.apply(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, datetime):
            raise TypeError('Cannot subtract datetime from offset!')
        elif type(other) == type(self):
            return self.__class__(self.n - other.n, **self.kwds)
        else:  # pragma: no cover
            raise TypeError('Cannot subtract %s from %s'
                            % (type(other), type(self)))

    def __rsub__(self, other):
        return self.__class__(-self.n, **self.kwds) + other

    def __mul__(self, someInt):
        return self.__class__(n=someInt * self.n, **self.kwds)

    def __rmul__(self, someInt):
        return self.__mul__(someInt)

    def __neg__(self):
        return self.__class__(-self.n, **self.kwds)

    def rollback(self, dt):
        """Roll provided date backward to next offset only if not on offset"""
        if type(dt) == date:
            dt = datetime(dt.year, dt.month, dt.day)

        if not self.onOffset(dt):
            dt = dt - self.__class__(1, **self.kwds)
        return dt

    def rollforward(self, dt):
        """Roll provided date forward to next offset only if not on offset"""
        if type(dt) == date:
            dt = datetime(dt.year, dt.month, dt.day)

        if not self.onOffset(dt):
            dt = dt + self.__class__(1, **self.kwds)
        return dt

    def onOffset(self, dt):
        # XXX, see #1395
        if type(self) == DateOffset or isinstance(self, Tick):
            return True

        # Default (slow) method for determining if some date is a member of the
        # date range generated by this offset. Subclasses may have this
        # re-implemented in a nicer way.
        a = dt
        b = ((dt + self) - self)
        return a == b

    @property
    def rule_code(self):
        raise NotImplementedError

    @property
    def freqstr(self):
        try:
            code = self.rule_code
        except NotImplementedError:
            return repr(self)

        if self.n != 1:
            fstr = '%d%s' % (self.n, code)
        else:
            fstr = code

        return fstr


class BusinessDay(CacheableOffset, DateOffset):
    """
    DateOffset subclass representing possibly n business days
    """
    def __init__(self, n=1, **kwds):
        self.n = int(n)
        self.kwds = kwds
        self.offset = kwds.get('offset', timedelta(0))
        self.normalize = kwds.get('normalize', False)

    @property
    def rule_code(self):
        return 'B'

    def __repr__(self):
        if hasattr(self, 'name') and len(self.name):
            return self.name

        className = getattr(self, '_outputName', self.__class__.__name__)
        attrs = []

        if self.offset:
            attrs = ['offset=%s' % repr(self.offset)]

        if abs(self.n) != 1:
            plural = 's'
        else:
            plural = ''

        out = '<%s ' % self.n + className + plural
        if attrs:
            out += ': ' + ', '.join(attrs)
        out += '>'
        return out

    @property
    def freqstr(self):
        try:
            code = self.rule_code
        except NotImplementedError:
            return repr(self)

        if self.n != 1:
            fstr = '%d%s' % (self.n, code)
        else:
            fstr = code

        if self.offset:
            fstr += self._offset_str()

        return fstr

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

    def isAnchored(self):
        return (self.n == 1)

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

            while n != 0:
                k = n // abs(n)
                result = result + timedelta(k)
                if result.weekday() < 5:
                    n -= k

            if self.normalize:
                result = datetime(result.year, result.month, result.day)

            if self.offset:
                result = result + self.offset

            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other,
                        normalize=self.normalize)
        else:
            raise TypeError('Only know how to combine business day with '
                            'datetime or timedelta!')

    @classmethod
    def onOffset(cls, dt):
        return dt.weekday() < 5


class CustomBusinessDay(BusinessDay):
    """
    **EXPERIMENTAL** DateOffset subclass representing possibly n business days
    excluding holidays

    .. warning:: EXPERIMENTAL

        This class is not officially supported and the API is likely to change
        in future versions. Use this at your own risk.

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
    """

    _cacheable = False

    def __init__(self, n=1, **kwds):
        # Check we have the required numpy version
        from distutils.version import LooseVersion
        if LooseVersion(np.__version__) < '1.7.0':
            raise NotImplementedError("CustomBusinessDay requires numpy >= "
                                      "1.7.0. Current version: " +
                                      np.__version__)

        self.n = int(n)
        self.kwds = kwds
        self.offset = kwds.get('offset', timedelta(0))
        self.normalize = kwds.get('normalize', False)
        self.weekmask = kwds.get('weekmask', 'Mon Tue Wed Thu Fri')

        holidays = kwds.get('holidays', [])
        holidays = [self._to_dt64(dt, dtype='datetime64[D]') for dt in
                    holidays]
        self.holidays = tuple(sorted(holidays))
        self.kwds['holidays'] = self.holidays
        self._set_busdaycalendar()

    def _set_busdaycalendar(self):
        holidays = np.array(self.holidays, dtype='datetime64[D]')
        self.busdaycalendar = np.busdaycalendar(holidays=holidays,
                                                weekmask=self.weekmask)

    def __getstate__(self):
        """"Return a pickleable state"""
        state = self.__dict__.copy()
        del state['busdaycalendar']
        return state

    def __setstate__(self, state):
        """Reconstruct an instance from a pickled state"""
        self.__dict__ = state
        self._set_busdaycalendar()

    @property
    def rule_code(self):
        return 'C'

    @staticmethod
    def _to_dt64(dt, dtype='datetime64'):
        if isinstance(dt, (datetime, basestring)):
            dt = np.datetime64(dt, dtype=dtype)
        if isinstance(dt, np.datetime64):
            dt = dt.astype(dtype)
        else:
            raise TypeError('dt must be datestring, datetime or datetime64')
        return dt

    def apply(self, other):
        if isinstance(other, datetime):
            dtype = type(other)
        elif isinstance(other, np.datetime64):
            dtype = other.dtype
        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other,
                        normalize=self.normalize)
        else:
            raise TypeError('Only know how to combine trading day with '
                            'datetime, datetime64 or timedelta!')
        dt64 = self._to_dt64(other)

        day64 = dt64.astype('datetime64[D]')
        time = dt64 - day64

        if self.n<=0:
            roll = 'forward'
        else:
            roll = 'backward'

        result = np.busday_offset(day64, self.n, roll=roll,
                                  busdaycal=self.busdaycalendar)

        if not self.normalize:
            result = result + time

        result = result.astype(dtype)

        if self.offset:
            result = result + self.offset

        return result

    def onOffset(self, dt):
        day64 = self._to_dt64(dt).astype('datetime64[D]')
        return np.is_busday(day64, busdaycal=self.busdaycalendar)


class MonthEnd(DateOffset, CacheableOffset):
    """DateOffset of one month end"""

    def apply(self, other):
        other = datetime(other.year, other.month, other.day, tzinfo=other.tzinfo)

        n = self.n
        _, days_in_month = tslib.monthrange(other.year, other.month)
        if other.day != days_in_month:
            other = other + relativedelta(months=-1, day=31)
            if n <= 0:
                n = n + 1
        other = other + relativedelta(months=n, day=31)
        return other

    @classmethod
    def onOffset(cls, dt):
        days_in_month = tslib.monthrange(dt.year, dt.month)[1]
        return dt.day == days_in_month

    @property
    def rule_code(self):
        return 'M'


class MonthBegin(DateOffset, CacheableOffset):
    """DateOffset of one month at beginning"""

    def apply(self, other):
        n = self.n

        if other.day > 1 and n <= 0:  # then roll forward if n<=0
            n += 1

        other = other + relativedelta(months=n, day=1)
        return other

    @classmethod
    def onOffset(cls, dt):
        return dt.day == 1

    @property
    def rule_code(self):
        return 'MS'


class BusinessMonthEnd(CacheableOffset, DateOffset):
    """DateOffset increments between business EOM dates"""

    def isAnchored(self):
        return (self.n == 1)

    def apply(self, other):
        other = datetime(other.year, other.month, other.day)

        n = self.n

        wkday, days_in_month = tslib.monthrange(other.year, other.month)
        lastBDay = days_in_month - max(((wkday + days_in_month - 1)
                                       % 7) - 4, 0)

        if n > 0 and not other.day >= lastBDay:
            n = n - 1
        elif n <= 0 and other.day > lastBDay:
            n = n + 1
        other = other + relativedelta(months=n, day=31)

        if other.weekday() > 4:
            other = other - BDay()
        return other

    @property
    def rule_code(self):
        return 'BM'


class BusinessMonthBegin(DateOffset, CacheableOffset):
    """DateOffset of one business month at beginning"""

    def apply(self, other):
        n = self.n

        wkday, _ = tslib.monthrange(other.year, other.month)
        first = _get_firstbday(wkday)

        if other.day > first and n <= 0:
            # as if rolled forward already
            n += 1
        elif other.day < first and n > 0:
            other = other + timedelta(days=first - other.day)
            n -= 1

        other = other + relativedelta(months=n)
        wkday, _ = tslib.monthrange(other.year, other.month)
        first = _get_firstbday(wkday)
        result = datetime(other.year, other.month, first)
        return result

    @classmethod
    def onOffset(cls, dt):
        first_weekday, _ = tslib.monthrange(dt.year, dt.month)
        if first_weekday == 5:
            return dt.day == 3
        elif first_weekday == 6:
            return dt.day == 2
        else:
            return dt.day == 1

    @property
    def rule_code(self):
        return 'BMS'


class Week(DateOffset, CacheableOffset):
    """
    Weekly offset

    Parameters
    ----------
    weekday : int, default None
        Always generate specific day of week. 0 for Monday
    """
    def __init__(self, n=1, **kwds):
        self.n = n
        self.weekday = kwds.get('weekday', None)

        if self.weekday is not None:
            if self.weekday < 0 or self.weekday > 6:
                raise ValueError('Day must be 0<=day<=6, got %d' %
                                self.weekday)

        self._inc = timedelta(weeks=1)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.weekday is not None)

    def apply(self, other):
        if self.weekday is None:
            return other + self.n * self._inc

        if self.n > 0:
            k = self.n
            otherDay = other.weekday()
            if otherDay != self.weekday:
                other = other + timedelta((self.weekday - otherDay) % 7)
                k = k - 1
            for i in xrange(k):
                other = other + self._inc
        else:
            k = self.n
            otherDay = other.weekday()
            if otherDay != self.weekday:
                other = other + timedelta((self.weekday - otherDay) % 7)
            for i in xrange(-k):
                other = other - self._inc
        return other

    def onOffset(self, dt):
        return dt.weekday() == self.weekday

    @property
    def rule_code(self):
        suffix = ''
        if self.weekday is not None:
            suffix = '-%s' % (_weekday_dict[self.weekday])
        return 'W' + suffix

_weekday_dict = {
    0: 'MON',
    1: 'TUE',
    2: 'WED',
    3: 'THU',
    4: 'FRI',
    5: 'SAT',
    6: 'SUN'
}


class WeekOfMonth(DateOffset, CacheableOffset):
    """
    Describes monthly dates like "the Tuesday of the 2nd week of each month"

    Parameters
    ----------
    n : int
    week : {0, 1, 2, 3, ...}
        0 is 1st week of month, 1 2nd week, etc.
    weekday : {0, 1, ..., 6}
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    """
    def __init__(self, n=1, **kwds):
        self.n = n
        self.weekday = kwds['weekday']
        self.week = kwds['week']

        if self.n == 0:
            raise ValueError('N cannot be 0')

        if self.weekday < 0 or self.weekday > 6:
            raise ValueError('Day must be 0<=day<=6, got %d' %
                            self.weekday)
        if self.week < 0 or self.week > 3:
            raise ValueError('Week must be 0<=day<=3, got %d' %
                            self.week)

        self.kwds = kwds

    def apply(self, other):
        offsetOfMonth = self.getOffsetOfMonth(other)

        if offsetOfMonth > other:
            if self.n > 0:
                months = self.n - 1
            else:
                months = self.n
        elif offsetOfMonth == other:
            months = self.n
        else:
            if self.n > 0:
                months = self.n
            else:
                months = self.n + 1

        return self.getOffsetOfMonth(other + relativedelta(months=months, day=1))

    def getOffsetOfMonth(self, dt):
        w = Week(weekday=self.weekday)
        d = datetime(dt.year, dt.month, 1)

        d = w.rollforward(d)

        for i in xrange(self.week):
            d = w.apply(d)

        return d

    def onOffset(self, dt):
        return dt == self.getOffsetOfMonth(dt)

    @property
    def rule_code(self):
        suffix = '-%d%s' % (self.week + 1, _weekday_dict.get(self.weekday, ''))
        return 'WOM' + suffix


class BQuarterEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...
    """
    _outputName = 'BusinessQuarterEnd'

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = BMonthEnd(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = tslib.monthrange(other.year, other.month)
        lastBDay = days_in_month - max(((wkday + days_in_month - 1)
                                       % 7) - 4, 0)

        monthsToGo = 3 - ((other.month - self.startingMonth) % 3)
        if monthsToGo == 3:
            monthsToGo = 0

        if n > 0 and not (other.day >= lastBDay and monthsToGo == 0):
            n = n - 1
        elif n <= 0 and other.day > lastBDay and monthsToGo == 0:
            n = n + 1

        other = other + relativedelta(months=monthsToGo + 3 * n, day=31)

        if other.weekday() > 4:
            other = other - BDay()

        return other

    def onOffset(self, dt):
        modMonth = (dt.month - self.startingMonth) % 3
        return BMonthEnd().onOffset(dt) and modMonth == 0

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.startingMonth]
        return 'BQ' + suffix


_month_dict = {
    1: 'JAN',
    2: 'FEB',
    3: 'MAR',
    4: 'APR',
    5: 'MAY',
    6: 'JUN',
    7: 'JUL',
    8: 'AUG',
    9: 'SEP',
    10: 'OCT',
    11: 'NOV',
    12: 'DEC'
}


class BQuarterBegin(DateOffset, CacheableOffset):
    _outputName = "BusinessQuarterBegin"

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = BMonthBegin(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, _ = tslib.monthrange(other.year, other.month)

        first = _get_firstbday(wkday)

        monthsSince = (other.month - self.startingMonth) % 3

        if n <= 0 and monthsSince != 0:  # make sure to roll forward so negate
            monthsSince = monthsSince - 3

        # roll forward if on same month later than first bday
        if n <= 0 and (monthsSince == 0 and other.day > first):
            n = n + 1
        # pretend to roll back if on same month but before firstbday
        elif n > 0 and (monthsSince == 0 and other.day < first):
            n = n - 1

        # get the first bday for result
        other = other + relativedelta(months=3 * n - monthsSince)
        wkday, _ = tslib.monthrange(other.year, other.month)
        first = _get_firstbday(wkday)
        result = datetime(other.year, other.month, first,
                          other.hour, other.minute, other.second,
                          other.microsecond)
        return result

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.startingMonth]
        return 'BQS' + suffix


class QuarterEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/31/2007, 6/30/2007, ...
    """
    _outputName = 'QuarterEnd'

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = MonthEnd(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = tslib.monthrange(other.year, other.month)

        monthsToGo = 3 - ((other.month - self.startingMonth) % 3)
        if monthsToGo == 3:
            monthsToGo = 0

        if n > 0 and not (other.day >= days_in_month and monthsToGo == 0):
            n = n - 1

        other = other + relativedelta(months=monthsToGo + 3 * n, day=31)

        return other

    def onOffset(self, dt):
        modMonth = (dt.month - self.startingMonth) % 3
        return MonthEnd().onOffset(dt) and modMonth == 0

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.startingMonth]
        return 'Q' + suffix


class QuarterBegin(DateOffset, CacheableOffset):
    _outputName = 'QuarterBegin'

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = MonthBegin(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = tslib.monthrange(other.year, other.month)

        monthsSince = (other.month - self.startingMonth) % 3

        if n <= 0 and monthsSince != 0:
            # make sure you roll forward, so negate
            monthsSince = monthsSince - 3

        if n < 0 and (monthsSince == 0 and other.day > 1):
            # after start, so come back an extra period as if rolled forward
            n = n + 1

        other = other + relativedelta(months=3 * n - monthsSince, day=1)
        return other

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.startingMonth]
        return 'QS' + suffix


class BYearEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business EOM dates"""
    _outputName = 'BusinessYearEnd'

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 12)

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = tslib.monthrange(other.year, self.month)
        lastBDay = (days_in_month -
                    max(((wkday + days_in_month - 1) % 7) - 4, 0))

        years = n
        if n > 0:
            if (other.month < self.month or
                    (other.month == self.month and other.day < lastBDay)):
                years -= 1
        elif n <= 0:
            if (other.month > self.month or
                    (other.month == self.month and other.day > lastBDay)):
                years += 1

        other = other + relativedelta(years=years)

        _, days_in_month = tslib.monthrange(other.year, self.month)
        result = datetime(other.year, self.month, days_in_month,
                          other.hour, other.minute, other.second,
                          other.microsecond)

        if result.weekday() > 4:
            result = result - BDay()

        return result

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.month]
        return 'BA' + suffix


class BYearBegin(DateOffset, CacheableOffset):
    """DateOffset increments between business year begin dates"""
    _outputName = 'BusinessYearBegin'

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 1)

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = tslib.monthrange(other.year, self.month)

        first = _get_firstbday(wkday)

        years = n

        if n > 0:  # roll back first for positive n
            if (other.month < self.month or
                    (other.month == self.month and other.day < first)):
                years -= 1
        elif n <= 0:  # roll forward
            if (other.month > self.month or
                    (other.month == self.month and other.day > first)):
                years += 1

        # set first bday for result
        other = other + relativedelta(years=years)
        wkday, days_in_month = tslib.monthrange(other.year, self.month)
        first = _get_firstbday(wkday)
        return datetime(other.year, self.month, first)

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.month]
        return 'BAS' + suffix


class YearEnd(DateOffset, CacheableOffset):
    """DateOffset increments between calendar year ends"""

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 12)

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        def _increment(date):
            if date.month == self.month:
                _, days_in_month = tslib.monthrange(date.year, self.month)
                if date.day != days_in_month:
                    year = date.year
                else:
                    year = date.year + 1
            elif date.month < self.month:
                year = date.year
            else:
                year = date.year + 1
            _, days_in_month = tslib.monthrange(year, self.month)
            return datetime(year, self.month, days_in_month,
                            date.hour, date.minute, date.second,
                            date.microsecond)

        def _decrement(date):
            year = date.year if date.month > self.month else date.year - 1
            _, days_in_month = tslib.monthrange(year, self.month)
            return datetime(year, self.month, days_in_month,
                            date.hour, date.minute, date.second,
                            date.microsecond)

        def _rollf(date):
            if (date.month != self.month or
                    date.day < tslib.monthrange(date.year, date.month)[1]):
                date = _increment(date)
            return date

        n = self.n
        result = other
        if n > 0:
            while n > 0:
                result = _increment(result)
                n -= 1
        elif n < 0:
            while n < 0:
                result = _decrement(result)
                n += 1
        else:
            # n == 0, roll forward
            result = _rollf(result)

        return result

    def onOffset(self, dt):
        wkday, days_in_month = tslib.monthrange(dt.year, self.month)
        return self.month == dt.month and dt.day == days_in_month

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.month]
        return 'A' + suffix


class YearBegin(DateOffset, CacheableOffset):
    """DateOffset increments between calendar year begin dates"""

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 1)

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        def _increment(date):
            year = date.year
            if date.month >= self.month:
                year += 1
            return datetime(year, self.month, 1, date.hour, date.minute,
                            date.second, date.microsecond)

        def _decrement(date):
            year = date.year
            if date.month < self.month or (date.month == self.month and
                                           date.day == 1):
                year -= 1
            return datetime(year, self.month, 1, date.hour, date.minute,
                            date.second, date.microsecond)

        def _rollf(date):
            if (date.month != self.month) or date.day > 1:
                date = _increment(date)
            return date

        n = self.n
        result = other
        if n > 0:
            while n > 0:
                result = _increment(result)
                n -= 1
        elif n < 0:
            while n < 0:
                result = _decrement(result)
                n += 1
        else:
            # n == 0, roll forward
            result = _rollf(result)

        return result

    def onOffset(self, dt):
        return dt.month == self.month and dt.day == 1

    @property
    def rule_code(self):
        suffix = '-%s' % _month_dict[self.month]
        return 'AS' + suffix


#----------------------------------------------------------------------
# Ticks

import operator


def _tick_comp(op):
    def f(self, other):
        return op(self.delta, other.delta)
    return f


class Tick(DateOffset):
    _inc = timedelta(microseconds=1000)

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
        return self.apply(other)

    def __eq__(self, other):
        if isinstance(other, basestring):
            from pandas.tseries.frequencies import to_offset
            other = to_offset(other)

        if isinstance(other, Tick):
            return self.delta == other.delta
        else:
            return DateOffset.__eq__(self, other)

    # This is identical to DateOffset.__hash__, but has to be redefined here
    # for Python 3, because we've redefined __eq__.
    def __hash__(self):
        return hash(self._params())

    def __ne__(self, other):
        if isinstance(other, basestring):
            from pandas.tseries.frequencies import to_offset
            other = to_offset(other)

        if isinstance(other, Tick):
            return self.delta != other.delta
        else:
            return DateOffset.__ne__(self, other)

    @property
    def delta(self):
        return self.n * self._inc

    @property
    def nanos(self):
        return _delta_to_nanoseconds(self.delta)

    def apply(self, other):
        if type(other) == date:
            other = datetime(other.year, other.month, other.day)

        if isinstance(other, (datetime, timedelta)):
            return other + self.delta
        elif isinstance(other, type(self)):
            return type(self)(self.n + other.n)
        else:  # pragma: no cover
            raise TypeError('Unhandled type: %s' % type(other))

    _rule_base = 'undefined'

    @property
    def rule_code(self):
        return self._rule_base

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
        nanos = _delta_to_nanoseconds(delta)
        if nanos % 1000000 == 0:
            return Milli(nanos // 1000000)
        elif nanos % 1000 == 0:
            return Micro(nanos // 1000)
        else:  # pragma: no cover
            return Nano(nanos)


def _delta_to_nanoseconds(delta):
    if isinstance(delta, Tick):
        delta = delta.delta
    return (delta.days * 24 * 60 * 60 * 1000000
            + delta.seconds * 1000000
            + delta.microseconds) * 1000


class Day(Tick, CacheableOffset):
    _inc = timedelta(1)
    _rule_base = 'D'


class Hour(Tick):
    _inc = timedelta(0, 3600)
    _rule_base = 'H'


class Minute(Tick):
    _inc = timedelta(0, 60)
    _rule_base = 'T'


class Second(Tick):
    _inc = timedelta(0, 1)
    _rule_base = 'S'


class Milli(Tick):
    _rule_base = 'L'


class Micro(Tick):
    _inc = timedelta(microseconds=1)
    _rule_base = 'U'


class Nano(Tick):
    _inc = 1
    _rule_base = 'N'

BDay = BusinessDay
BMonthEnd = BusinessMonthEnd
BMonthBegin = BusinessMonthBegin
CDay = CustomBusinessDay


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

    Note
    ----
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

    if end and not offset.onOffset(end):
        end = offset.rollback(end)

        if periods is None and end < start:
            end = None
            periods = 0

    if end is None:
        end = start + (periods - 1) * offset

    if start is None:
        start = end - (periods - 1) * offset

    cur = start

    next_date = cur
    while cur <= end:
        yield cur

        # faster than cur + offset
        next_date = offset.apply(cur)
        if next_date <= cur:
            raise ValueError('Offset %s did not increment date' % offset)
        cur = next_date
