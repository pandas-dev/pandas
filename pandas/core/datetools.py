"""A collection of random tools for dealing with dates in Python"""

from datetime import datetime, timedelta
import sys

try:
    import dateutil
    from dateutil import parser
    from dateutil.relativedelta import relativedelta

    # raise exception if dateutil 2.0 install on 2.x platform
    if (sys.version_info[0] == 2 and
        dateutil.__version__ == '2.0'):  # pragma: no cover
        raise Exception('dateutil 2.0 incompatible with Python 2.x, you must '
                        'install version 1.5!')
except ImportError: # pragma: no cover
    print 'Please install python-dateutil via easy_install or some method!'
    raise # otherwise a 2nd import won't show the message

import calendar

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

def to_datetime(arg):
    """Attempts to convert arg to datetime"""
    if arg is None or isinstance(arg, datetime):
        return arg
    try:
        return parser.parse(arg)
    except Exception:
        return arg

def normalize_date(dt):
    return datetime(dt.year, dt.month, dt.day)

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
        return someDate == ((someDate + self) - self)


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
        if isinstance(other, datetime):
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
            other = other + relativedelta(months=-1, day=31)
            if n <= 0:
                n = n + 1
        other = other + relativedelta(months=n, day=31)
        return other

    @classmethod
    def onOffset(cls, someDate):
        __junk, days_in_month = calendar.monthrange(someDate.year,
                                                   someDate.month)
        return someDate.day == days_in_month

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
        other = other + relativedelta(months=n, day=31)

        if other.weekday() > 4:
            other = other - BDay()
        return other


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

        one_month = relativedelta(months=1, day=1)

        if other < offsetOfMonth:
            if self.n > 0:
                months = self.n - 1
            else:
                months = self.n
        elif other == offsetOfMonth:
            months = self.n
        else:
            if self.n > 0:
                months = self.n
            else:
                months = self.n + 1

        return self.getOffsetOfMonth(other + relativedelta(months=months, day=1))

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

        if self.startingMonth < 1 or self.startingMonth > 3:
            raise Exception('Start month must be 1<=day<=3, got %d'
                            % self.startingMonth)

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

        other = other + relativedelta(months=monthsToGo + 3*n, day=31)

        if other.weekday() > 4:
            other = other - BDay()

        return other

    def onOffset(self, someDate):
        modMonth = (someDate.month - self.startingMonth) % 3
        return BMonthEnd().onOffset(someDate) and modMonth == 0

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

        if self.startingMonth < 1 or self.startingMonth > 3:
            raise Exception('Start month must be 1<=day<=3, got %d'
                            % self.startingMonth)

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

        other = other + relativedelta(months=monthsToGo + 3*n, day=31)

        return other

    def onOffset(self, someDate):
        modMonth = (someDate.month - self.startingMonth) % 3
        return MonthEnd().onOffset(someDate) and modMonth == 0

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

        other = other + relativedelta(years=years)

        _, days_in_month = calendar.monthrange(other.year, self.month)
        result = datetime(other.year, self.month, days_in_month)

        if result.weekday() > 4:
            result = result - BDay()

        return result

class YearEnd(DateOffset, CacheableOffset):
    """DateOffset increments between calendar year ends"""
    _normalizeFirst = True

    def apply(self, other):
        n = self.n
        if other.month != 12 or other.day != 31:
            other = datetime(other.year - 1, 12, 31)
            if n <= 0:
                n = n + 1
        other = other + relativedelta(years=n)
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
        other = other + relativedelta(years = n, day=1)
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
        if isinstance(other, (datetime, timedelta)):
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
        if second == (first + offset) and third == (second + offset):
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
