from datetime import date, datetime, timedelta
import operator

import numpy as np

from pandas._libs.tslibs import (
    Timedelta,
    Timestamp,
    ccalendar,
    conversion,
    frequencies as libfrequencies,
    offsets as liboffsets,
)
from pandas._libs.tslibs.offsets import (  # noqa:F401
    FY5253,
    ApplyTypeError,
    BaseOffset,
    BQuarterBegin,
    BQuarterEnd,
    BusinessMixin,
    BusinessMonthBegin,
    BusinessMonthEnd,
    BYearBegin,
    BYearEnd,
    CustomMixin,
    Day,
    Easter,
    FY5253Quarter,
    Hour,
    LastWeekOfMonth,
    Micro,
    Milli,
    Minute,
    MonthBegin,
    MonthEnd,
    Nano,
    QuarterBegin,
    QuarterEnd,
    Second,
    SemiMonthBegin,
    SemiMonthEnd,
    SingleConstructorOffset,
    Tick,
    Week,
    WeekOfMonth,
    YearBegin,
    YearEnd,
    apply_index_wraps,
    apply_wraps,
    is_normalized,
    shift_month,
    to_dt64D,
)
from pandas.util._decorators import cache_readonly, doc

__all__ = [
    "Day",
    "BusinessDay",
    "BDay",
    "CustomBusinessDay",
    "CDay",
    "CBMonthEnd",
    "CBMonthBegin",
    "MonthBegin",
    "BMonthBegin",
    "MonthEnd",
    "BMonthEnd",
    "SemiMonthEnd",
    "SemiMonthBegin",
    "BusinessHour",
    "CustomBusinessHour",
    "YearBegin",
    "BYearBegin",
    "YearEnd",
    "BYearEnd",
    "QuarterBegin",
    "BQuarterBegin",
    "QuarterEnd",
    "BQuarterEnd",
    "LastWeekOfMonth",
    "FY5253Quarter",
    "FY5253",
    "Week",
    "WeekOfMonth",
    "Easter",
    "Hour",
    "Minute",
    "Second",
    "Milli",
    "Micro",
    "Nano",
    "DateOffset",
]


# ---------------------------------------------------------------------
# DateOffset


class OffsetMeta(type):
    """
    Metaclass that allows us to pretend that all BaseOffset subclasses
    inherit from DateOffset (which is needed for backward-compatibility).
    """

    @classmethod
    def __instancecheck__(cls, obj) -> bool:
        return isinstance(obj, BaseOffset)

    @classmethod
    def __subclasscheck__(cls, obj) -> bool:
        return issubclass(obj, BaseOffset)


class DateOffset(liboffsets.RelativeDeltaOffset, metaclass=OffsetMeta):
    """
    Standard kind of date increment used for a date range.

    Works exactly like relativedelta in terms of the keyword args you
    pass in, use of the keyword n is discouraged-- you would be better
    off specifying n in the keywords you use, but regardless it is
    there for you. n is needed for DateOffset subclasses.

    DateOffset work as follows.  Each offset specify a set of dates
    that conform to the DateOffset.  For example, Bday defines this
    set to be the set of dates that are weekdays (M-F).  To test if a
    date is in the set of a DateOffset dateOffset we can use the
    is_on_offset method: dateOffset.is_on_offset(date).

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

    Parameters
    ----------
    n : int, default 1
        The number of time periods the offset represents.
    normalize : bool, default False
        Whether to round the result of a DateOffset addition down to the
        previous midnight.
    **kwds
        Temporal parameter that add to or replace the offset value.

        Parameters that **add** to the offset (like Timedelta):

        - years
        - months
        - weeks
        - days
        - hours
        - minutes
        - seconds
        - microseconds
        - nanoseconds

        Parameters that **replace** the offset value:

        - year
        - month
        - day
        - weekday
        - hour
        - minute
        - second
        - microsecond
        - nanosecond.

    See Also
    --------
    dateutil.relativedelta.relativedelta : The relativedelta type is designed
        to be applied to an existing datetime an can replace specific components of
        that datetime, or represents an interval of time.

    Examples
    --------
    >>> from pandas.tseries.offsets import DateOffset
    >>> ts = pd.Timestamp('2017-01-01 09:10:11')
    >>> ts + DateOffset(months=3)
    Timestamp('2017-04-01 09:10:11')

    >>> ts = pd.Timestamp('2017-01-01 09:10:11')
    >>> ts + DateOffset(months=2)
    Timestamp('2017-03-01 09:10:11')
    """

    pass


class BusinessDay(BusinessMixin):
    """
    DateOffset subclass representing possibly n business days.
    """

    _prefix = "B"
    _attributes = frozenset(["n", "normalize", "offset"])

    def __reduce__(self):
        tup = (self.n, self.normalize, self.offset)
        return type(self), tup

    def _offset_str(self) -> str:
        def get_str(td):
            off_str = ""
            if td.days > 0:
                off_str += str(td.days) + "D"
            if td.seconds > 0:
                s = td.seconds
                hrs = int(s / 3600)
                if hrs != 0:
                    off_str += str(hrs) + "H"
                    s -= hrs * 3600
                mts = int(s / 60)
                if mts != 0:
                    off_str += str(mts) + "Min"
                    s -= mts * 60
                if s != 0:
                    off_str += str(s) + "s"
            if td.microseconds > 0:
                off_str += str(td.microseconds) + "us"
            return off_str

        if isinstance(self.offset, timedelta):
            zero = timedelta(0, 0, 0)
            if self.offset >= zero:
                off_str = "+" + get_str(self.offset)
            else:
                off_str = "-" + get_str(-self.offset)
            return off_str
        else:
            return "+" + repr(self.offset)

    @apply_wraps
    def apply(self, other):
        if isinstance(other, datetime):
            n = self.n
            wday = other.weekday()

            # avoid slowness below by operating on weeks first
            weeks = n // 5
            if n <= 0 and wday > 4:
                # roll forward
                n += 1

            n -= 5 * weeks

            # n is always >= 0 at this point
            if n == 0 and wday > 4:
                # roll back
                days = 4 - wday
            elif wday > 4:
                # roll forward
                days = (7 - wday) + (n - 1)
            elif wday + n <= 4:
                # shift by n days without leaving the current week
                days = n
            else:
                # shift by n days plus 2 to get past the weekend
                days = n + 2

            result = other + timedelta(days=7 * weeks + days)
            if self.offset:
                result = result + self.offset
            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other, normalize=self.normalize)
        else:
            raise ApplyTypeError(
                "Only know how to combine business day with datetime or timedelta."
            )

    @apply_index_wraps
    def apply_index(self, i):
        time = i.to_perioddelta("D")
        # to_period rolls forward to next BDay; track and
        # reduce n where it does when rolling forward
        asper = i.to_period("B")

        if self.n > 0:
            shifted = (i.to_perioddelta("B") - time).asi8 != 0

            # Integer-array addition is deprecated, so we use
            # _time_shift directly
            roll = np.where(shifted, self.n - 1, self.n)
            shifted = asper._addsub_int_array(roll, operator.add)
        else:
            # Integer addition is deprecated, so we use _time_shift directly
            roll = self.n
            shifted = asper._time_shift(roll)

        result = shifted.to_timestamp() + time
        return result

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        return dt.weekday() < 5


class BusinessHour(liboffsets.BusinessHourMixin):
    """
    DateOffset subclass representing possibly n business hours.
    """

    _prefix = "BH"
    _anchor = 0
    _attributes = frozenset(["n", "normalize", "start", "end", "offset"])

    @cache_readonly
    def next_bday(self):
        """
        Used for moving to next business day.
        """
        if self.n >= 0:
            nb_offset = 1
        else:
            nb_offset = -1
        if self._prefix.startswith("C"):
            # CustomBusinessHour
            return CustomBusinessDay(
                n=nb_offset,
                weekmask=self.weekmask,
                holidays=self.holidays,
                calendar=self.calendar,
            )
        else:
            return BusinessDay(n=nb_offset)

    def _next_opening_time(self, other, sign=1):
        """
        If self.n and sign have the same sign, return the earliest opening time
        later than or equal to current time.
        Otherwise the latest opening time earlier than or equal to current
        time.

        Opening time always locates on BusinessDay.
        However, closing time may not if business hour extends over midnight.

        Parameters
        ----------
        other : datetime
            Current time.
        sign : int, default 1.
            Either 1 or -1. Going forward in time if it has the same sign as
            self.n. Going backward in time otherwise.

        Returns
        -------
        result : datetime
            Next opening time.
        """
        earliest_start = self.start[0]
        latest_start = self.start[-1]

        if not self.next_bday.is_on_offset(other):
            # today is not business day
            other = other + sign * self.next_bday
            if self.n * sign >= 0:
                hour, minute = earliest_start.hour, earliest_start.minute
            else:
                hour, minute = latest_start.hour, latest_start.minute
        else:
            if self.n * sign >= 0:
                if latest_start < other.time():
                    # current time is after latest starting time in today
                    other = other + sign * self.next_bday
                    hour, minute = earliest_start.hour, earliest_start.minute
                else:
                    # find earliest starting time no earlier than current time
                    for st in self.start:
                        if other.time() <= st:
                            hour, minute = st.hour, st.minute
                            break
            else:
                if other.time() < earliest_start:
                    # current time is before earliest starting time in today
                    other = other + sign * self.next_bday
                    hour, minute = latest_start.hour, latest_start.minute
                else:
                    # find latest starting time no later than current time
                    for st in reversed(self.start):
                        if other.time() >= st:
                            hour, minute = st.hour, st.minute
                            break

        return datetime(other.year, other.month, other.day, hour, minute)

    def _prev_opening_time(self, other):
        """
        If n is positive, return the latest opening time earlier than or equal
        to current time.
        Otherwise the earliest opening time later than or equal to current
        time.

        Parameters
        ----------
        other : datetime
            Current time.

        Returns
        -------
        result : datetime
            Previous opening time.
        """
        return self._next_opening_time(other, sign=-1)

    @apply_wraps
    def rollback(self, dt):
        """
        Roll provided date backward to next offset only if not on offset.
        """
        if not self.is_on_offset(dt):
            if self.n >= 0:
                dt = self._prev_opening_time(dt)
            else:
                dt = self._next_opening_time(dt)
            return self._get_closing_time(dt)
        return dt

    @apply_wraps
    def rollforward(self, dt):
        """
        Roll provided date forward to next offset only if not on offset.
        """
        if not self.is_on_offset(dt):
            if self.n >= 0:
                return self._next_opening_time(dt)
            else:
                return self._prev_opening_time(dt)
        return dt

    @apply_wraps
    def apply(self, other):
        if isinstance(other, datetime):
            # used for detecting edge condition
            nanosecond = getattr(other, "nanosecond", 0)
            # reset timezone and nanosecond
            # other may be a Timestamp, thus not use replace
            other = datetime(
                other.year,
                other.month,
                other.day,
                other.hour,
                other.minute,
                other.second,
                other.microsecond,
            )
            n = self.n

            # adjust other to reduce number of cases to handle
            if n >= 0:
                if other.time() in self.end or not self._is_on_offset(other):
                    other = self._next_opening_time(other)
            else:
                if other.time() in self.start:
                    # adjustment to move to previous business day
                    other = other - timedelta(seconds=1)
                if not self._is_on_offset(other):
                    other = self._next_opening_time(other)
                    other = self._get_closing_time(other)

            # get total business hours by sec in one business day
            businesshours = sum(
                self._get_business_hours_by_sec(st, en)
                for st, en in zip(self.start, self.end)
            )

            bd, r = divmod(abs(n * 60), businesshours // 60)
            if n < 0:
                bd, r = -bd, -r

            # adjust by business days first
            if bd != 0:
                if isinstance(self, CustomMixin):  # GH 30593
                    skip_bd = CustomBusinessDay(
                        n=bd,
                        weekmask=self.weekmask,
                        holidays=self.holidays,
                        calendar=self.calendar,
                    )
                else:
                    skip_bd = BusinessDay(n=bd)
                # midnight business hour may not on BusinessDay
                if not self.next_bday.is_on_offset(other):
                    prev_open = self._prev_opening_time(other)
                    remain = other - prev_open
                    other = prev_open + skip_bd + remain
                else:
                    other = other + skip_bd

            # remaining business hours to adjust
            bhour_remain = timedelta(minutes=r)

            if n >= 0:
                while bhour_remain != timedelta(0):
                    # business hour left in this business time interval
                    bhour = (
                        self._get_closing_time(self._prev_opening_time(other)) - other
                    )
                    if bhour_remain < bhour:
                        # finish adjusting if possible
                        other += bhour_remain
                        bhour_remain = timedelta(0)
                    else:
                        # go to next business time interval
                        bhour_remain -= bhour
                        other = self._next_opening_time(other + bhour)
            else:
                while bhour_remain != timedelta(0):
                    # business hour left in this business time interval
                    bhour = self._next_opening_time(other) - other
                    if (
                        bhour_remain > bhour
                        or bhour_remain == bhour
                        and nanosecond != 0
                    ):
                        # finish adjusting if possible
                        other += bhour_remain
                        bhour_remain = timedelta(0)
                    else:
                        # go to next business time interval
                        bhour_remain -= bhour
                        other = self._get_closing_time(
                            self._next_opening_time(
                                other + bhour - timedelta(seconds=1)
                            )
                        )

            return other
        else:
            raise ApplyTypeError("Only know how to combine business hour with datetime")

    def is_on_offset(self, dt):
        if self.normalize and not is_normalized(dt):
            return False

        if dt.tzinfo is not None:
            dt = datetime(
                dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, dt.microsecond
            )
        # Valid BH can be on the different BusinessDay during midnight
        # Distinguish by the time spent from previous opening time
        return self._is_on_offset(dt)

    def _is_on_offset(self, dt):
        """
        Slight speedups using calculated values.
        """
        # if self.normalize and not is_normalized(dt):
        #     return False
        # Valid BH can be on the different BusinessDay during midnight
        # Distinguish by the time spent from previous opening time
        if self.n >= 0:
            op = self._prev_opening_time(dt)
        else:
            op = self._next_opening_time(dt)
        span = (dt - op).total_seconds()
        businesshours = 0
        for i, st in enumerate(self.start):
            if op.hour == st.hour and op.minute == st.minute:
                businesshours = self._get_business_hours_by_sec(st, self.end[i])
        if span <= businesshours:
            return True
        else:
            return False


class CustomBusinessDay(CustomMixin, BusinessDay):
    """
    DateOffset subclass representing custom business days excluding holidays.

    Parameters
    ----------
    n : int, default 1
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range.
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        Weekmask of valid business days, passed to ``numpy.busdaycalendar``.
    holidays : list
        List/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``.
    calendar : pd.HolidayCalendar or np.busdaycalendar
    offset : timedelta, default timedelta(0)
    """

    _prefix = "C"
    _attributes = frozenset(
        ["n", "normalize", "weekmask", "holidays", "calendar", "offset"]
    )

    def __reduce__(self):
        # np.holidaycalendar cant be pickled, so pass None there and
        #  it will be re-constructed within __init__
        tup = (self.n, self.normalize, self.weekmask, self.holidays, None, self.offset)
        return type(self), tup

    def __init__(
        self,
        n=1,
        normalize=False,
        weekmask="Mon Tue Wed Thu Fri",
        holidays=None,
        calendar=None,
        offset=timedelta(0),
    ):
        BusinessDay.__init__(self, n, normalize, offset)
        CustomMixin.__init__(self, weekmask, holidays, calendar)

    @apply_wraps
    def apply(self, other):
        if self.n <= 0:
            roll = "forward"
        else:
            roll = "backward"

        if isinstance(other, datetime):
            date_in = other
            np_dt = np.datetime64(date_in.date())

            np_incr_dt = np.busday_offset(
                np_dt, self.n, roll=roll, busdaycal=self.calendar
            )

            dt_date = np_incr_dt.astype(datetime)
            result = datetime.combine(dt_date, date_in.time())

            if self.offset:
                result = result + self.offset
            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other, normalize=self.normalize)
        else:
            raise ApplyTypeError(
                "Only know how to combine trading day with "
                "datetime, datetime64 or timedelta."
            )

    def apply_index(self, i):
        raise NotImplementedError

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        day64 = to_dt64D(dt)
        return np.is_busday(day64, busdaycal=self.calendar)


class CustomBusinessHour(CustomMixin, BusinessHour):
    """
    DateOffset subclass representing possibly n custom business days.
    """

    _prefix = "CBH"
    _anchor = 0
    _attributes = frozenset(
        ["n", "normalize", "weekmask", "holidays", "calendar", "start", "end", "offset"]
    )

    def __init__(
        self,
        n=1,
        normalize=False,
        weekmask="Mon Tue Wed Thu Fri",
        holidays=None,
        calendar=None,
        start="09:00",
        end="17:00",
        offset=timedelta(0),
    ):
        BusinessHour.__init__(self, n, normalize, start=start, end=end, offset=offset)
        CustomMixin.__init__(self, weekmask, holidays, calendar)

    def __reduce__(self):
        # None for self.calendar bc np.busdaycalendar doesnt pickle nicely
        return (
            type(self),
            (
                self.n,
                self.normalize,
                self.weekmask,
                self.holidays,
                None,
                self.start,
                self.end,
                self.offset,
            ),
        )


# ---------------------------------------------------------------------
# Month-Based Offset Classes


@doc(bound="bound")
class _CustomBusinessMonth(CustomMixin, BusinessMixin, liboffsets.MonthOffset):
    """
    DateOffset subclass representing custom business month(s).

    Increments between {bound} of month dates.

    Parameters
    ----------
    n : int, default 1
        The number of months represented.
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range.
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        Weekmask of valid business days, passed to ``numpy.busdaycalendar``.
    holidays : list
        List/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``.
    calendar : pd.HolidayCalendar or np.busdaycalendar
        Calendar to integrate.
    offset : timedelta, default timedelta(0)
        Time offset to apply.
    """

    _attributes = frozenset(
        ["n", "normalize", "weekmask", "holidays", "calendar", "offset"]
    )

    is_on_offset = BaseOffset.is_on_offset  # override MonthOffset method
    apply_index = BaseOffset.apply_index  # override MonthOffset method

    def __init__(
        self,
        n=1,
        normalize=False,
        weekmask="Mon Tue Wed Thu Fri",
        holidays=None,
        calendar=None,
        offset=timedelta(0),
    ):
        BusinessMixin.__init__(self, n, normalize, offset)
        CustomMixin.__init__(self, weekmask, holidays, calendar)

    def __reduce__(self):
        # None for self.calendar bc np.busdaycalendar doesnt pickle nicely
        return (
            type(self),
            (self.n, self.normalize, self.weekmask, self.holidays, None, self.offset),
        )

    @cache_readonly
    def cbday_roll(self):
        """
        Define default roll function to be called in apply method.
        """
        cbday = CustomBusinessDay(n=self.n, normalize=False, **self.kwds)

        if self._prefix.endswith("S"):
            # MonthBegin
            roll_func = cbday.rollforward
        else:
            # MonthEnd
            roll_func = cbday.rollback
        return roll_func

    @cache_readonly
    def m_offset(self):
        if self._prefix.endswith("S"):
            # MonthBegin
            moff = MonthBegin(n=1, normalize=False)
        else:
            # MonthEnd
            moff = MonthEnd(n=1, normalize=False)
        return moff

    @cache_readonly
    def month_roll(self):
        """
        Define default roll function to be called in apply method.
        """
        if self._prefix.endswith("S"):
            # MonthBegin
            roll_func = self.m_offset.rollback
        else:
            # MonthEnd
            roll_func = self.m_offset.rollforward
        return roll_func

    @apply_wraps
    def apply(self, other):
        # First move to month offset
        cur_month_offset_date = self.month_roll(other)

        # Find this custom month offset
        compare_date = self.cbday_roll(cur_month_offset_date)
        n = liboffsets.roll_convention(other.day, self.n, compare_date.day)

        new = cur_month_offset_date + n * self.m_offset
        result = self.cbday_roll(new)
        return result


@doc(_CustomBusinessMonth, bound="end")
class CustomBusinessMonthEnd(_CustomBusinessMonth):
    _prefix = "CBM"


@doc(_CustomBusinessMonth, bound="beginning")
class CustomBusinessMonthBegin(_CustomBusinessMonth):
    _prefix = "CBMS"


# ---------------------------------------------------------------------

BDay = BusinessDay
BMonthEnd = BusinessMonthEnd
BMonthBegin = BusinessMonthBegin
CBMonthEnd = CustomBusinessMonthEnd
CBMonthBegin = CustomBusinessMonthBegin
CDay = CustomBusinessDay

prefix_mapping = {
    offset._prefix: offset
    for offset in [
        YearBegin,  # 'AS'
        YearEnd,  # 'A'
        BYearBegin,  # 'BAS'
        BYearEnd,  # 'BA'
        BusinessDay,  # 'B'
        BusinessMonthBegin,  # 'BMS'
        BusinessMonthEnd,  # 'BM'
        BQuarterEnd,  # 'BQ'
        BQuarterBegin,  # 'BQS'
        BusinessHour,  # 'BH'
        CustomBusinessDay,  # 'C'
        CustomBusinessMonthEnd,  # 'CBM'
        CustomBusinessMonthBegin,  # 'CBMS'
        CustomBusinessHour,  # 'CBH'
        MonthEnd,  # 'M'
        MonthBegin,  # 'MS'
        Nano,  # 'N'
        SemiMonthEnd,  # 'SM'
        SemiMonthBegin,  # 'SMS'
        Week,  # 'W'
        Second,  # 'S'
        Minute,  # 'T'
        Micro,  # 'U'
        QuarterEnd,  # 'Q'
        QuarterBegin,  # 'QS'
        Milli,  # 'L'
        Hour,  # 'H'
        Day,  # 'D'
        WeekOfMonth,  # 'WOM'
        FY5253,
        FY5253Quarter,
    ]
}
