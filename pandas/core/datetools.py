"""A collection of random tools for dealing with dates in Python"""
from datetime import datetime, timedelta
import sys
import numpy as np
import pandas._tseries as lib

from pandas._tseries import Timestamp

try:
    import dateutil
    from dateutil import parser

    # raise exception if dateutil 2.0 install on 2.x platform
    if (sys.version_info[0] == 2 and
        dateutil.__version__ == '2.0'):  # pragma: no cover
        raise Exception('dateutil 2.0 incompatible with Python 2.x, you must '
                        'install version 1.5!')
except ImportError: # pragma: no cover
    print 'Please install python-dateutil via easy_install or some method!'
    raise # otherwise a 2nd import won't show the message

import calendar #NOTE: replace with _tseries.monthrange

#-------------------------------------------------------------------------------
# Boxing and unboxing

def _dt_box(key, offset=None, tzinfo=None):
    '''
    timestamp-like (int64, python datetime, etc.) => Timestamp
    '''
    return Timestamp(key, offset=offset, tzinfo=None)

def _dt_box_array(arr, offset=None, tzinfo=None):
    if arr is None:
        return arr

    boxfunc = lambda x: _dt_box(x, offset=offset, tzinfo=tzinfo)
    boxer = np.frompyfunc(boxfunc, 1, 1)
    return boxer(arr)

def _dt_unbox(key):
    '''
    Timestamp-like => dt64
    '''
    if type(key) == float:
        raise TypeError("Cannot unbox a float to datetime")
    return np.datetime64(lib.pydt_to_i8(key))

def _dt_unbox_array(arr):
    if arr is None:
        return arr

    unboxer = np.frompyfunc(_dt_unbox, 1, 1)
    return unboxer(arr)

#-------------------------------------------------------------------------------
# Miscellaneous date functions

def format(dt):
    """Returns date in YYYYMMDD format."""
    return dt.strftime('%Y%m%d')

OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)

def ole2datetime(oledt):
    """function for converting excel date to normal date format"""
    val = float(oledt)

    # Excel has a bug where it thinks the date 2/29/1900 exists
    # we just reject any date before 3/1/1900.
    if val < 61:
        raise Exception("Value is outside of acceptable range: %s " % val)

    return OLE_TIME_ZERO + timedelta(days=val)

def to_timestamp(arg, offset=None):
    """ Attempts to convert arg to timestamp """
    if arg is None:
        return arg

    if isinstance(arg, basestring):
        try:
            arg = parser.parse(arg)
        except Exception:
            pass

    return lib.Timestamp(arg, offset=offset)

def to_datetime(arg):
    """Attempts to convert arg to datetime"""
    if arg is None:
        return arg
    elif isinstance(arg, datetime):
        return arg
    try:
        return parser.parse(arg)
    except Exception:
        return arg

def normalize_date(dt):
    if isinstance(dt, np.datetime64):
        dt = _dt_box(dt)
    return datetime(dt.year, dt.month, dt.day)

def _get_firstbday(wkday):
    """
    wkday is the result of calendar.monthrange(year, month)

    If it's a saturday or sunday, increment first business day to reflect this
    """
    firstBDay = 1
    if wkday == 5: # on Saturday
        firstBDay = 3
    elif wkday == 6: # on Sunday
        firstBDay = 2
    return firstBDay

#-------------------------------------------------------------------------------
# DateOffset

class CacheableOffset(object):
    pass

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
      date = rollback(date) # does nothing is date is valid
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
    # For some offsets, want to drop the time information off the
    # first date
    _normalizeFirst = False
    def __init__(self, n=1, **kwds):
        self.n = int(n)
        self.kwds = kwds
        if len(kwds) > 0:
            self._offset = lib.Delta(**kwds)
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

    def _params(self):
        attrs = [(k, v) for k, v in vars(self).iteritems()
                 if k not in ['kwds', '_offset']]
        attrs.extend(self.kwds.items())
        attrs = sorted(set(attrs))

        params = tuple([str(self.__class__)] + attrs)
        return params

    def __repr__(self):
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
        else: # pragma: no cover
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

    def rollback(self, someDate):
        """Roll provided date backward to next offset only if not on offset"""
        if self._normalizeFirst:
            someDate = normalize_date(someDate)

        if not self.onOffset(someDate):
            someDate = someDate - self.__class__(1, **self.kwds)
        return someDate

    def rollforward(self, someDate):
        """Roll provided date forward to next offset only if not on offset"""
        if self._normalizeFirst:
            someDate = normalize_date(someDate)

        if not self.onOffset(someDate):
            someDate = someDate + self.__class__(1, **self.kwds)
        return someDate

    def onOffset(self, someDate):
        if type(self) == DateOffset:
            return True

        # Default (slow) method for determining if some date is a
        # member of the DateRange generated by this offset. Subclasses
        # may have this re-implemented in a nicer way.
        a = someDate
        b = ((someDate + self) - self)
        return a == b


class BDay(DateOffset, CacheableOffset):
    """
    DateOffset subclass representing possibly n business days
    """
    _normalizeFirst = True
    _outputName = 'BusinessDay'
    def __init__(self, n=1, **kwds):
        self.n = int(n)
        self.kwds = kwds
        self.offset = kwds.get('offset', timedelta(0))
        self.normalize = kwds.get('normalize', True)

    def __repr__(self):
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

    def isAnchored(self):
        return (self.n == 1)

    def apply(self, other):
        if isinstance(other, (datetime, Timestamp)):
            n = self.n

            if n == 0 and other.weekday() > 4:
                n = 1

            result = other

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
            raise Exception('Only know how to combine business day with '
                            'datetime or timedelta!')
    @classmethod
    def onOffset(cls, someDate):
        return someDate.weekday() < 5


class MonthEnd(DateOffset, CacheableOffset):
    """DateOffset of one month end"""

    _normalizeFirst = True

    def apply(self, other):
        n = self.n
        _, days_in_month = calendar.monthrange(other.year, other.month)
        if other.day != days_in_month:
            other = other + lib.Delta(months=-1, day=31)
            if n <= 0:
                n = n + 1
        other = other + lib.Delta(months=n, day=31)
        return other

    @classmethod
    def onOffset(cls, someDate):
        __junk, days_in_month = calendar.monthrange(someDate.year,
                                                   someDate.month)
        return someDate.day == days_in_month

class MonthBegin(DateOffset, CacheableOffset):
    """DateOffset of one month at beginning"""

    _normalizeFirst = True

    def apply(self, other):
        n = self.n

        if other.day > 1 and n<=0: #then roll forward if n<=0
            n =+ 1

        other = other + lib.Delta(months=n, day=1)
        return other

    @classmethod
    def onOffset(cls, someDate):
        firstDay, _ = calendar.monthrange(someDate.year, someDate.month)
        return someDate.day == (firstDay + 1)

class BMonthEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business EOM dates"""
    _outputName = 'BusinessMonthEnd'
    _normalizeFirst = True

    def isAnchored(self):
        return (self.n == 1)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = calendar.monthrange(other.year, other.month)
        lastBDay = days_in_month - max(((wkday + days_in_month - 1) % 7) - 4, 0)

        if n > 0 and not other.day >= lastBDay:
            n = n - 1
        elif n <= 0 and other.day > lastBDay:
            n = n + 1
        other = other + lib.Delta(months=n, day=31)

        if other.weekday() > 4:
            other = other - BDay()
        return other

class BMonthBegin(DateOffset, CacheableOffset):
    """DateOffset of one business month at beginning"""

    _normalizeFirst = True

    def apply(self, other):
        n = self.n

        wkday, _ = calendar.monthrange(other.year, other.month)
        firstBDay = _get_firstbday(wkday)

        if other.day > firstBDay and n<=0:
            # as if rolled forward already
            n += 1

        other = other + lib.Delta(months=n)
        wkday, _ = calendar.monthrange(other.year, other.month)
        firstBDay = _get_firstbday(wkday)
        result = datetime(other.year, other.month, firstBDay)
        return result


class Week(DateOffset, CacheableOffset):
    """
    Weekly offset

    Parameters
    ----------
    weekday : int, default None
        Always generate specific day of week. 0 for Monday
    """
    _normalizeFirst = True
    def __init__(self, n=1, **kwds):
        self.n = n
        self.weekday = kwds.get('weekday', None)

        if self.weekday is not None:
            if self.weekday < 0 or self.weekday > 6:
                raise Exception('Day must be 0<=day<=6, got %d' %
                                self.weekday)

        self.inc = timedelta(weeks=1)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.weekday is not None)

    def apply(self, other):
        if self.weekday is None:
            return other + self.n * self.inc

        if self.n > 0:
            k = self.n
            otherDay = other.weekday()
            if otherDay != self.weekday:
                other = other + timedelta((self.weekday - otherDay) % 7)
                k = k - 1
            for i in xrange(k):
                other = other + self.inc
        else:
            k = self.n
            otherDay = other.weekday()
            if otherDay != self.weekday:
                other = other + timedelta((self.weekday - otherDay) % 7)
            for i in xrange(-k):
                other = other - self.inc
        return other

    def onOffset(self, someDate):
        return someDate.weekday() == self.weekday


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
        1: Tuedays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    """
    _normalizeFirst = True
    def __init__(self, n=1, **kwds):
        self.n = n
        self.weekday = kwds['weekday']
        self.week = kwds['week']

        if self.n == 0:
            raise Exception('N cannot be 0')

        if self.weekday < 0 or self.weekday > 6:
            raise Exception('Day must be 0<=day<=6, got %d' %
                            self.weekday)
        if self.week < 0 or self.week > 3:
            raise Exception('Week must be 0<=day<=3, got %d' %
                            self.week)

        self.kwds = kwds

    def apply(self, other):
        offsetOfMonth = self.getOffsetOfMonth(other)

        one_month = lib.Delta(months=1, day=1)

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

        return self.getOffsetOfMonth(other + lib.Delta(months=months, day=1))

    def getOffsetOfMonth(self, someDate):
        w = Week(weekday=self.weekday)
        d = datetime(someDate.year, someDate.month, 1)

        d = w.rollforward(d)

        for i in xrange(self.week):
            d = w.apply(d)

        return d

    def onOffset(self, someDate):
        return someDate == self.getOffsetOfMonth(someDate)

class BQuarterEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...
    """
    _outputName = 'BusinessQuarterEnd'
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = BMonthEnd(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = calendar.monthrange(other.year, other.month)
        lastBDay = days_in_month - max(((wkday + days_in_month - 1) % 7) - 4, 0)

        monthsToGo = 3 - ((other.month - self.startingMonth) % 3)
        if monthsToGo == 3:
            monthsToGo = 0

        if n > 0 and not (other.day >= lastBDay and monthsToGo == 0):
            n = n - 1
        elif n <= 0 and other.day > lastBDay and monthsToGo == 0:
            n = n + 1

        other = other + lib.Delta(months=monthsToGo + 3*n, day=31)

        if other.weekday() > 4:
            other = other - BDay()

        return other

    def onOffset(self, someDate):
        modMonth = (someDate.month - self.startingMonth) % 3
        return BMonthEnd().onOffset(someDate) and modMonth == 0

class BQuarterBegin(DateOffset, CacheableOffset):
    _outputName = "BusinessQuarterBegin"
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = BMonthBegin(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        if self._normalizeFirst:
            other = normalize_date(other)

        wkday, _ = calendar.monthrange(other.year, other.month)

        firstBDay = _get_firstbday(wkday)

        monthsSince = (other.month - self.startingMonth) % 3
        if monthsSince == 3: # on offset
            monthsSince = 0

        if n <= 0 and monthsSince != 0: # make sure to roll forward so negate
            monthsSince = monthsSince - 3

        # roll forward if on same month later than first bday
        if n <= 0 and (monthsSince == 0 and other.day > firstBDay):
            n = n + 1
        # pretend to roll back if on same month but before firstbday
        elif n > 0 and (monthsSince == 0 and other.day < firstBDay):
            n = n - 1

        # get the first bday for result
        other = other + lib.Delta(months=3*n - monthsSince)
        wkday, _ = calendar.monthrange(other.year, other.month)
        firstBDay = _get_firstbday(wkday)
        result = datetime(other.year, other.month, firstBDay)
        return result


class QuarterEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/31/2007, 6/30/2007, ...
    """
    _outputName = 'QuarterEnd'
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = MonthEnd(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = calendar.monthrange(other.year, other.month)

        monthsToGo = 3 - ((other.month - self.startingMonth) % 3)
        if monthsToGo == 3:
            monthsToGo = 0

        if n > 0 and not (other.day >= days_in_month and monthsToGo == 0):
            n = n - 1

        other = other + lib.Delta(months=monthsToGo + 3*n, day=31)

        return other

    def onOffset(self, someDate):
        modMonth = (someDate.month - self.startingMonth) % 3
        return MonthEnd().onOffset(someDate) and modMonth == 0

class QuarterBegin(DateOffset, CacheableOffset):
    _outputName = 'QuarterBegin'
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.n = n
        self.startingMonth = kwds.get('startingMonth', 3)

        self.offset = MonthBegin(3)
        self.kwds = kwds

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    def apply(self, other):
        n = self.n

        wkday, days_in_month = calendar.monthrange(other.year, other.month)

        monthsSince = (other.month - self.startingMonth) % 3

        if monthsSince == 3: # on an offset
            monthsSince = 0

        if n <= 0 and monthsSince != 0:
            # make sure you roll forward, so negate
            monthsSince = monthsSince - 3

        if n < 0 and (monthsSince == 0 and other.day > 1):
            # after start, so come back an extra period as if rolled forward
            n = n + 1

        other = other + lib.Delta(months=3*n - monthsSince, day=1)
        return other


class BYearEnd(DateOffset, CacheableOffset):
    """DateOffset increments between business EOM dates"""
    _outputName = 'BusinessYearEnd'
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 12)

        if self.month < 1 or self.month > 12:
            raise Exception('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        n = self.n

        if self._normalizeFirst:
            other = normalize_date(other)

        wkday, days_in_month = calendar.monthrange(other.year, self.month)
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

        other = other + lib.Delta(years=years)

        _, days_in_month = calendar.monthrange(other.year, self.month)
        result = datetime(other.year, self.month, days_in_month)

        if result.weekday() > 4:
            result = result - BDay()

        return result

class BYearBegin(DateOffset, CacheableOffset):
    """DateOffset increments between business year begin dates"""
    _outputName = 'BusinessYearBegin'
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 1)

        if self.month < 1 or self.month > 12:
            raise Exception('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        n = self.n

        if self._normalizeFirst:
            other = normalize_date(other)

        wkday, days_in_month = calendar.monthrange(other.year, self.month)

        firstBDay = _get_firstbday(wkday)

        years = n


        if n > 0: # roll back first for positive n
            if (other.month < self.month or
                (other.month == self.month and other.day < firstBDay)):
                years -= 1
        elif n <= 0: # roll forward
            if (other.month > self.month or
                (other.month == self.month and other.day > firstBDay)):
                years += 1

        # set first bday for result
        other = other + lib.Delta(years = years)
        wkday, days_in_month = calendar.monthrange(other.year, self.month)
        firstBDay = _get_firstbday(wkday)
        result = datetime(other.year, self.month, firstBDay)
        return result

class YearEnd(DateOffset, CacheableOffset):
    """DateOffset increments between calendar year ends"""
    _normalizeFirst = True

    def __init__(self, n=1, **kwds):
        self.month = kwds.get('month', 12)

        if self.month < 1 or self.month > 12:
            raise Exception('Month must go from 1 to 12')

        DateOffset.__init__(self, n=n, **kwds)

    def apply(self, other):
        n = self.n
        if other.month != 12 or other.day != 31:
            other = datetime(other.year - 1, 12, 31)
            if n <= 0:
                n = n + 1
        other = other + lib.Delta(years=n)
        return other

    @classmethod
    def onOffset(cls, someDate):
        return someDate.month == 12 and someDate.day == 31


class YearBegin(DateOffset, CacheableOffset):
    """DateOffset increments between calendar year begin dates"""
    _normalizeFirst = True

    def apply(self, other):
        n = self.n
        if other.month != 1 or other.day != 1:
            other = datetime(other.year, 1, 1)
            if n <= 0:
                n = n + 1
        other = other + lib.Delta(years = n, day=1)
        return other

    @classmethod
    def onOffset(cls, someDate):
        return someDate.month == 1 and someDate.day == 1

#-------------------------------------------------------------------------------
# Ticks

class Tick(DateOffset):
    _normalizeFirst = False
    _delta = None
    _inc = timedelta(microseconds=1000)

    @property
    def delta(self):
        if self._delta is None:
            self._delta = self.n * self._inc

        return self._delta

    def apply(self, other):
        if isinstance(other, (datetime, timedelta, Timestamp)):
            return other + self.delta
        elif isinstance(other, type(self)):
            return type(self)(self.n + other.n)

class Hour(Tick):
    _inc = timedelta(0, 3600)

class Minute(Tick):
    _inc = timedelta(0, 60)

class Second(Tick):
    _inc = timedelta(0, 1)

day = DateOffset()
bday = BDay(normalize=True)
businessDay = bday
monthEnd = MonthEnd()
yearEnd = YearEnd()
yearBegin = YearBegin()
bmonthEnd = BMonthEnd()
businessMonthEnd = bmonthEnd
bquarterEnd = BQuarterEnd()
quarterEnd = QuarterEnd()
byearEnd = BYearEnd()
week = Week()


# Functions/offsets to roll dates forward
thisMonthEnd = MonthEnd(0)
thisBMonthEnd = BMonthEnd(0)
thisYearEnd = YearEnd(0)
thisYearBegin = YearBegin(0)
thisBQuarterEnd = BQuarterEnd(0)
thisQuarterEnd = QuarterEnd(0)

# Functions to check where a date lies
isBusinessDay = BDay().onOffset
isMonthEnd = MonthEnd().onOffset
isBMonthEnd = BMonthEnd().onOffset

#-------------------------------------------------------------------------------
# Offset names ("time rules") and related functions

_offsetMap = {
    "WEEKDAY"  : BDay(1),
    "EOM"      : BMonthEnd(1),
    "W@MON"    : Week(weekday=0),
    "W@TUE"    : Week(weekday=1),
    "W@WED"    : Week(weekday=2),
    "W@THU"    : Week(weekday=3),
    "W@FRI"    : Week(weekday=4),
    "Q@JAN"    : BQuarterEnd(startingMonth=1),
    "Q@FEB"    : BQuarterEnd(startingMonth=2),
    "Q@MAR"    : BQuarterEnd(startingMonth=3),
    "A@JAN"    : BYearEnd(month=1),
    "A@FEB"    : BYearEnd(month=2),
    "A@MAR"    : BYearEnd(month=3),
    "A@APR"    : BYearEnd(month=4),
    "A@MAY"    : BYearEnd(month=5),
    "A@JUN"    : BYearEnd(month=6),
    "A@JUL"    : BYearEnd(month=7),
    "A@AUG"    : BYearEnd(month=8),
    "A@SEP"    : BYearEnd(month=9),
    "A@OCT"    : BYearEnd(month=10),
    "A@NOV"    : BYearEnd(month=11),
    "A@DEC"    : BYearEnd()
}

_newoffsetMap = {
       # Annual - Calendar
       "A@JAN" : YearEnd(month=1),
       "A@FEB" : YearEnd(month=2),
       "A@MAR" : YearEnd(month=3),
       "A@APR" : YearEnd(month=4),
       "A@MAY" : YearEnd(month=5),
       "A@JUN" : YearEnd(month=6),
       "A@JUL" : YearEnd(month=7),
       "A@AUG" : YearEnd(month=8),
       "A@SEP" : YearEnd(month=9),
       "A@OCT" : YearEnd(month=10),
       "A@NOV" : YearEnd(month=11),
       "A@DEC" : YearEnd(month=12),
       "A"     : YearEnd(month=12),
       # Annual - Calendar (start)
       "AS@JAN" : YearBegin(month=1),
       "AS"     : YearBegin(month=1),
       "AS@FEB" : YearBegin(month=2),
       "AS@MAR" : YearBegin(month=3),
       "AS@APR" : YearBegin(month=4),
       "AS@MAY" : YearBegin(month=5),
       "AS@JUN" : YearBegin(month=6),
       "AS@JUL" : YearBegin(month=7),
       "AS@AUG" : YearBegin(month=8),
       "AS@SEP" : YearBegin(month=9),
       "AS@OCT" : YearBegin(month=10),
       "AS@NOV" : YearBegin(month=11),
       "AS@DEC" : YearBegin(month=12),
       # Annual - Business
       "BA@JAN" : BYearEnd(month=1),
       "BA@FEB" : BYearEnd(month=2),
       "BA@MAR" : BYearEnd(month=3),
       "BA@APR" : BYearEnd(month=4),
       "BA@MAY" : BYearEnd(month=5),
       "BA@JUN" : BYearEnd(month=6),
       "BA@JUL" : BYearEnd(month=7),
       "BA@AUG" : BYearEnd(month=8),
       "BA@SEP" : BYearEnd(month=9),
       "BA@OCT" : BYearEnd(month=10),
       "BA@NOV" : BYearEnd(month=11),
       "BA@DEC" : BYearEnd(month=12),
       "BA"     : BYearEnd(month=12),
       # Annual - Business (Start)
       "BAS@JAN" : BYearBegin(month=1),
       "BAS"     : BYearBegin(month=1),
       "BAS@FEB" : BYearBegin(month=2),
       "BAS@MAR" : BYearBegin(month=3),
       "BAS@APR" : BYearBegin(month=4),
       "BAS@MAY" : BYearBegin(month=5),
       "BAS@JUN" : BYearBegin(month=6),
       "BAS@JUL" : BYearBegin(month=7),
       "BAS@AUG" : BYearBegin(month=8),
       "BAS@SEP" : BYearBegin(month=9),
       "BAS@OCT" : BYearBegin(month=10),
       "BAS@NOV" : BYearBegin(month=11),
       "BAS@DEC" : BYearBegin(month=12),
       # Quarterly - Calendar
       "Q@JAN" : QuarterEnd(startingMonth=1),
       "Q@FEB" : QuarterEnd(startingMonth=2),
       "Q@MAR" : QuarterEnd(startingMonth=3),
       "Q"     : QuarterEnd(startingMonth=3),
       "Q@APR" : QuarterEnd(startingMonth=4),
       "Q@MAY" : QuarterEnd(startingMonth=5),
       "Q@JUN" : QuarterEnd(startingMonth=6),
       "Q@JUL" : QuarterEnd(startingMonth=7),
       "Q@AUG" : QuarterEnd(startingMonth=8),
       "Q@SEP" : QuarterEnd(startingMonth=9),
       "Q@OCT" : QuarterEnd(startingMonth=10),
       "Q@NOV" : QuarterEnd(startingMonth=11),
       "Q@DEC" : QuarterEnd(startingMonth=12),
       # Quarterly - Calendar (Start)
       "QS@JAN" : QuarterBegin(startingMonth=1),
       "QS"     : QuarterBegin(startingMonth=1),
       "QS@FEB" : QuarterBegin(startingMonth=2),
       "QS@MAR" : QuarterBegin(startingMonth=3),
       "QS@APR" : QuarterBegin(startingMonth=4),
       "QS@MAY" : QuarterBegin(startingMonth=5),
       "QS@JUN" : QuarterBegin(startingMonth=6),
       "QS@JUL" : QuarterBegin(startingMonth=7),
       "QS@AUG" : QuarterBegin(startingMonth=8),
       "QS@SEP" : QuarterBegin(startingMonth=9),
       "QS@OCT" : QuarterBegin(startingMonth=10),
       "QS@NOV" : QuarterBegin(startingMonth=11),
       "QS@DEC" : QuarterBegin(startingMonth=12),
       # Quarterly - Business
       "BQ@JAN" : BQuarterEnd(startingMonth=1),
       "BQ@FEB" : BQuarterEnd(startingMonth=2),
       "BQ@MAR" : BQuarterEnd(startingMonth=3),
       "BQ"     : BQuarterEnd(startingMonth=3),
       "BQ@APR" : BQuarterEnd(startingMonth=4),
       "BQ@MAY" : BQuarterEnd(startingMonth=5),
       "BQ@JUN" : BQuarterEnd(startingMonth=6),
       "BQ@JUL" : BQuarterEnd(startingMonth=7),
       "BQ@AUG" : BQuarterEnd(startingMonth=8),
       "BQ@SEP" : BQuarterEnd(startingMonth=9),
       "BQ@OCT" : BQuarterEnd(startingMonth=10),
       "BQ@NOV" : BQuarterEnd(startingMonth=11),
       "BQ@DEC" : BQuarterEnd(startingMonth=12),
       # Quarterly - Business (Start)
       "BQS@JAN" : BQuarterBegin(startingMonth=1),
       "BQS"     : BQuarterBegin(startingMonth=1),
       "BQS@FEB" : BQuarterBegin(startingMonth=2),
       "BQS@MAR" : BQuarterBegin(startingMonth=3),
       "BQS@APR" : BQuarterBegin(startingMonth=4),
       "BQS@MAY" : BQuarterBegin(startingMonth=5),
       "BQS@JUN" : BQuarterBegin(startingMonth=6),
       "BQS@JUL" : BQuarterBegin(startingMonth=7),
       "BQS@AUG" : BQuarterBegin(startingMonth=8),
       "BQS@SEP" : BQuarterBegin(startingMonth=9),
       "BQS@OCT" : BQuarterBegin(startingMonth=10),
       "BQS@NOV" : BQuarterBegin(startingMonth=11),
       "BQS@DEC" : BQuarterBegin(startingMonth=12),
       # Monthly - Calendar
       "M"      : MonthEnd(),
       "EOM"    : MonthEnd(),
       "MS"     : MonthBegin(),
       "SOM"    : MonthBegin(),
       # Monthly - Business
       "BM"     : BMonthEnd(),
       "BEOM"   : BMonthEnd(),
       "BMS"    : BMonthBegin(),
       "BSOM"   : BMonthBegin(),
       # Weekly
       "W@MON" : Week(weekday=0),
       "WS"    : Week(weekday=0),
       "BWS"   : Week(weekday=0),
       "W@TUE" : Week(weekday=1),
       "W@WED" : Week(weekday=2),
       "W@THU" : Week(weekday=3),
       "W@FRI" : Week(weekday=4),
       "BW"    : Week(weekday=4),
       "W@SAT" : Week(weekday=5),
       "W@SUN" : Week(weekday=6),
       "W"     : Week(weekday=6),
       "D"     : DateOffset(),
       "B"     : BDay(),
       "H"     : Hour(),
       "Min"   : Minute(),
       "S"     : Second(),
       "U"     : None,
       None    : None,
        }


for i, weekday in enumerate(['MON', 'TUE', 'WED', 'THU', 'FRI']):
    for iweek in xrange(4):
        _offsetMap['WOM@%d%s' % (iweek + 1, weekday)] = \
            WeekOfMonth(week=iweek, weekday=i)

_offsetNames = dict([(v, k) for k, v in _offsetMap.iteritems()])

def inferTimeRule(index):
    if len(index) < 3:
        raise Exception('Need at least three dates to infer time rule!')

    first, second, third = index[:3]
    for rule, offset in _offsetMap.iteritems():
        if (first + offset) == second and (second + offset) == third:
            return rule

    raise Exception('Could not infer time rule from data!')

def getOffset(name):
    """
    Return DateOffset object associated with rule name

    Example
    -------
    getOffset('EOM') --> BMonthEnd(1)
    """
    offset = _offsetMap.get(name)
    if offset is not None:
        return offset
    else:
        raise Exception('Bad rule name requested: %s!' % name)

def hasOffsetName(offset):
    return offset in _offsetNames

def getOffsetName(offset):
    """
    Return rule name associated with a DateOffset object

    Example
    -------
    getOffsetName(BMonthEnd(1)) --> 'EOM'
    """
    name = _offsetNames.get(offset)
    if name is not None:
        return name
    else:
        raise Exception('Bad offset name requested: %s!' % offset)

def _infer_tzinfo(start, end):
    def _infer(a, b):
        tz = a.tzinfo
        if b and b.tzinfo:
            assert(tz == b.tzinfo)
        return tz
    tz = None
    if start is not None:
        tz = _infer(start, end)
    elif end is not None:
        tz = _infer(end, start)
    return tz

def _will_use_cache(offset):
    return (offset.isAnchored() and isinstance(offset, CacheableOffset))

def _figure_out_timezone(start, end, tzinfo):
    inferred_tz = _infer_tzinfo(start, end)
    tz = inferred_tz
    if inferred_tz is None and tzinfo is not None:
        tz = tzinfo
    elif tzinfo is not None:
        assert(inferred_tz == tzinfo)
        # make tz naive for now

    start = start if start is None else start.replace(tzinfo=None)
    end = end if end is None else end.replace(tzinfo=None)

    return start, end, tz

_CACHE_START = Timestamp(datetime(1950, 1, 1))
_CACHE_END   = Timestamp(datetime(2030, 1, 1))

_daterange_cache = {}

def generate_range(start=_CACHE_START, end=_CACHE_END, periods=None,
                   offset=BDay(), freq=None):
    """
    Generates a sequence of dates corresponding to the specified time
    offset. Similar to dateutil.rrule except uses pandas DateOffset
    objects to represent time increments

    Parameters
    ----------
    start : timestamp-like (default None)
    end : timestamp-like (default None)
    periods : int, optional

    Note
    ----
    * This method is faster for generating weekdays than dateutil.rrule
    * At least two of (start, end, periods) must be specified.
    * If both start and end are specified, the returned dates will
    satisfy start <= date <= end.

    Returns
    -------
    dates : generator object

    See also
    --------
    DateRange, dateutil.rrule
    """

    if freq is not None:
        offset = getOffset(freq)

    if freq is None:
        if offset in _offsetNames:
            freq = _offsetNames[offset]

    start = to_timestamp(start)
    end = to_timestamp(end)

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
    if offset._normalizeFirst:
        cur = normalize_date(cur)

    next_date = cur
    while cur <= end:
        yield cur

        # faster than cur + offset
        next_date = offset.apply(cur)
        if next_date <= cur:
            raise ValueError('Offset %s did not increment date' % offset)
        cur = next_date

def _naive_in_cache_range(start, end):
    if start is None or end is None:
        return False
    else:
        return _in_range(start, end, _CACHE_START, _CACHE_END)

def _in_range(start, end, rng_start, rng_end):
    return start > rng_start and end < rng_end
