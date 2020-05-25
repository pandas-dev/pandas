from datetime import date, datetime, timedelta
import operator

from dateutil.easter import easter
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
    Hour,
    Micro,
    Milli,
    Minute,
    MonthBegin,
    MonthEnd,
    Nano,
    QuarterBegin,
    QuarterEnd,
    Second,
    SingleConstructorOffset,
    Tick,
    YearBegin,
    YearEnd,
    apply_index_wraps,
    apply_wraps,
    as_datetime,
    is_normalized,
    shift_month,
    to_dt64D,
)
from pandas.errors import AbstractMethodError
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


class DateOffset(BaseOffset, metaclass=OffsetMeta):
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

    _attributes = frozenset(["n", "normalize"] + list(liboffsets.relativedelta_kwds))
    _adjust_dst = False

    def __init__(self, n=1, normalize=False, **kwds):
        BaseOffset.__init__(self, n, normalize)

        off, use_rd = liboffsets._determine_offset(kwds)
        object.__setattr__(self, "_offset", off)
        object.__setattr__(self, "_use_relativedelta", use_rd)
        for key in kwds:
            val = kwds[key]
            object.__setattr__(self, key, val)

    @apply_wraps
    def apply(self, other):
        if self._use_relativedelta:
            other = as_datetime(other)

        if len(self.kwds) > 0:
            tzinfo = getattr(other, "tzinfo", None)
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
                other = conversion.localize_pydatetime(other, tzinfo)

            return Timestamp(other)
        else:
            return other + timedelta(self.n)

    @apply_index_wraps
    def apply_index(self, i):
        """
        Vectorized apply of DateOffset to DatetimeIndex,
        raises NotImplementedError for offsets without a
        vectorized implementation.

        Parameters
        ----------
        i : DatetimeIndex

        Returns
        -------
        y : DatetimeIndex
        """
        kwds = self.kwds
        relativedelta_fast = {
            "years",
            "months",
            "weeks",
            "days",
            "hours",
            "minutes",
            "seconds",
            "microseconds",
        }
        # relativedelta/_offset path only valid for base DateOffset
        if self._use_relativedelta and set(kwds).issubset(relativedelta_fast):

            months = (kwds.get("years", 0) * 12 + kwds.get("months", 0)) * self.n
            if months:
                shifted = liboffsets.shift_months(i.asi8, months)
                i = type(i)(shifted, dtype=i.dtype)

            weeks = (kwds.get("weeks", 0)) * self.n
            if weeks:
                # integer addition on PeriodIndex is deprecated,
                #   so we directly use _time_shift instead
                asper = i.to_period("W")
                shifted = asper._time_shift(weeks)
                i = shifted.to_timestamp() + i.to_perioddelta("W")

            timedelta_kwds = {
                k: v
                for k, v in kwds.items()
                if k in ["days", "hours", "minutes", "seconds", "microseconds"]
            }
            if timedelta_kwds:
                delta = Timedelta(**timedelta_kwds)
                i = i + (self.n * delta)
            return i
        elif not self._use_relativedelta and hasattr(self, "_offset"):
            # timedelta
            return i + (self._offset * self.n)
        else:
            # relativedelta with other keywords
            kwd = set(kwds) - relativedelta_fast
            raise NotImplementedError(
                "DateOffset with relativedelta "
                f"keyword(s) {kwd} not able to be "
                "applied vectorized"
            )

    def is_on_offset(self, dt):
        if self.normalize and not is_normalized(dt):
            return False
        # TODO, see #1395
        return True


class BusinessDay(BusinessMixin, SingleConstructorOffset):
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
# Semi-Month Based Offset Classes


class SemiMonthOffset(SingleConstructorOffset):
    _default_day_of_month = 15
    _min_day_of_month = 2
    _attributes = frozenset(["n", "normalize", "day_of_month"])

    def __init__(self, n=1, normalize=False, day_of_month=None):
        BaseOffset.__init__(self, n, normalize)

        if day_of_month is None:
            day_of_month = self._default_day_of_month

        object.__setattr__(self, "day_of_month", int(day_of_month))
        if not self._min_day_of_month <= self.day_of_month <= 27:
            raise ValueError(
                "day_of_month must be "
                f"{self._min_day_of_month}<=day_of_month<=27, "
                f"got {self.day_of_month}"
            )

    @classmethod
    def _from_name(cls, suffix=None):
        return cls(day_of_month=suffix)

    @property
    def rule_code(self) -> str:
        suffix = f"-{self.day_of_month}"
        return self._prefix + suffix

    @apply_wraps
    def apply(self, other):
        # shift `other` to self.day_of_month, incrementing `n` if necessary
        n = liboffsets.roll_convention(other.day, self.n, self.day_of_month)

        days_in_month = ccalendar.get_days_in_month(other.year, other.month)

        # For SemiMonthBegin on other.day == 1 and
        # SemiMonthEnd on other.day == days_in_month,
        # shifting `other` to `self.day_of_month` _always_ requires
        # incrementing/decrementing `n`, regardless of whether it is
        # initially positive.
        if type(self) is SemiMonthBegin and (self.n <= 0 and other.day == 1):
            n -= 1
        elif type(self) is SemiMonthEnd and (self.n > 0 and other.day == days_in_month):
            n += 1

        return self._apply(n, other)

    def _apply(self, n, other):
        """
        Handle specific apply logic for child classes.
        """
        raise AbstractMethodError(self)

    @apply_index_wraps
    def apply_index(self, i):
        # determine how many days away from the 1st of the month we are
        dti = i
        days_from_start = i.to_perioddelta("M").asi8
        delta = Timedelta(days=self.day_of_month - 1).value

        # get boolean array for each element before the day_of_month
        before_day_of_month = days_from_start < delta

        # get boolean array for each element after the day_of_month
        after_day_of_month = days_from_start > delta

        # determine the correct n for each date in i
        roll = self._get_roll(i, before_day_of_month, after_day_of_month)

        # isolate the time since it will be striped away one the next line
        time = i.to_perioddelta("D")

        # apply the correct number of months

        # integer-array addition on PeriodIndex is deprecated,
        #  so we use _addsub_int_array directly
        asper = i.to_period("M")

        shifted = asper._addsub_int_array(roll // 2, operator.add)
        i = type(dti)(shifted.to_timestamp())

        # apply the correct day
        i = self._apply_index_days(i, roll)

        return i + time

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        """
        Return an array with the correct n for each date in i.

        The roll array is based on the fact that i gets rolled back to
        the first day of the month.
        """
        raise AbstractMethodError(self)

    def _apply_index_days(self, i, roll):
        """
        Apply the correct day for each date in i.
        """
        raise AbstractMethodError(self)


class SemiMonthEnd(SemiMonthOffset):
    """
    Two DateOffset's per month repeating on the last
    day of the month and day_of_month.

    Parameters
    ----------
    n : int
    normalize : bool, default False
    day_of_month : int, {1, 3,...,27}, default 15
    """

    _prefix = "SM"
    _min_day_of_month = 1

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        days_in_month = ccalendar.get_days_in_month(dt.year, dt.month)
        return dt.day in (self.day_of_month, days_in_month)

    def _apply(self, n, other):
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
        """
        Add days portion of offset to DatetimeIndex i.

        Parameters
        ----------
        i : DatetimeIndex
        roll : ndarray[int64_t]

        Returns
        -------
        result : DatetimeIndex
        """
        nanos = (roll % 2) * Timedelta(days=self.day_of_month).value
        i += nanos.astype("timedelta64[ns]")
        return i + Timedelta(days=-1)


class SemiMonthBegin(SemiMonthOffset):
    """
    Two DateOffset's per month repeating on the first
    day of the month and day_of_month.

    Parameters
    ----------
    n : int
    normalize : bool, default False
    day_of_month : int, {2, 3,...,27}, default 15
    """

    _prefix = "SMS"

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        return dt.day in (1, self.day_of_month)

    def _apply(self, n, other):
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
        """
        Add days portion of offset to DatetimeIndex i.

        Parameters
        ----------
        i : DatetimeIndex
        roll : ndarray[int64_t]

        Returns
        -------
        result : DatetimeIndex
        """
        nanos = (roll % 2) * Timedelta(days=self.day_of_month - 1).value
        return i + nanos.astype("timedelta64[ns]")


# ---------------------------------------------------------------------
# Week-Based Offset Classes


class Week(SingleConstructorOffset):
    """
    Weekly offset.

    Parameters
    ----------
    weekday : int, default None
        Always generate specific day of week. 0 for Monday.
    """

    _inc = timedelta(weeks=1)
    _prefix = "W"
    _attributes = frozenset(["n", "normalize", "weekday"])

    def __init__(self, n=1, normalize=False, weekday=None):
        BaseOffset.__init__(self, n, normalize)
        object.__setattr__(self, "weekday", weekday)

        if self.weekday is not None:
            if self.weekday < 0 or self.weekday > 6:
                raise ValueError(f"Day must be 0<=day<=6, got {self.weekday}")

    def is_anchored(self) -> bool:
        return self.n == 1 and self.weekday is not None

    @apply_wraps
    def apply(self, other):
        if self.weekday is None:
            return other + self.n * self._inc

        if not isinstance(other, datetime):
            raise TypeError(
                f"Cannot add {type(other).__name__} to {type(self).__name__}"
            )

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
            # integer addition on PeriodIndex is deprecated,
            #  so we use _time_shift directly
            asper = i.to_period("W")

            shifted = asper._time_shift(self.n)
            return shifted.to_timestamp() + i.to_perioddelta("W")
        else:
            return self._end_apply_index(i)

    def _end_apply_index(self, dtindex):
        """
        Add self to the given DatetimeIndex, specialized for case where
        self.weekday is non-null.

        Parameters
        ----------
        dtindex : DatetimeIndex

        Returns
        -------
        result : DatetimeIndex
        """
        off = dtindex.to_perioddelta("D")

        base, mult = libfrequencies.get_freq_code(self.freqstr)
        base_period = dtindex.to_period(base)

        if self.n > 0:
            # when adding, dates on end roll to next
            normed = dtindex - off + Timedelta(1, "D") - Timedelta(1, "ns")
            roll = np.where(
                base_period.to_timestamp(how="end") == normed, self.n, self.n - 1
            )
            # integer-array addition on PeriodIndex is deprecated,
            #  so we use _addsub_int_array directly
            shifted = base_period._addsub_int_array(roll, operator.add)
            base = shifted.to_timestamp(how="end")
        else:
            # integer addition on PeriodIndex is deprecated,
            #  so we use _time_shift directly
            roll = self.n
            base = base_period._time_shift(roll).to_timestamp(how="end")

        return base + off + Timedelta(1, "ns") - Timedelta(1, "D")

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        elif self.weekday is None:
            return True
        return dt.weekday() == self.weekday

    @property
    def rule_code(self) -> str:
        suffix = ""
        if self.weekday is not None:
            weekday = ccalendar.int_to_weekday[self.weekday]
            suffix = f"-{weekday}"
        return self._prefix + suffix

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            weekday = None
        else:
            weekday = ccalendar.weekday_to_int[suffix]
        return cls(weekday=weekday)


class WeekOfMonth(liboffsets.WeekOfMonthMixin):
    """
    Describes monthly dates like "the Tuesday of the 2nd week of each month".

    Parameters
    ----------
    n : int
    week : int {0, 1, 2, 3, ...}, default 0
        A specific integer for the week of the month.
        e.g. 0 is 1st week of month, 1 is the 2nd week, etc.
    weekday : int {0, 1, ..., 6}, default 0
        A specific integer for the day of the week.

        - 0 is Monday
        - 1 is Tuesday
        - 2 is Wednesday
        - 3 is Thursday
        - 4 is Friday
        - 5 is Saturday
        - 6 is Sunday.
    """

    _prefix = "WOM"
    _attributes = frozenset(["n", "normalize", "week", "weekday"])

    def __init__(self, n=1, normalize=False, week=0, weekday=0):
        liboffsets.WeekOfMonthMixin.__init__(self, n, normalize, weekday)
        object.__setattr__(self, "week", week)

        if self.week < 0 or self.week > 3:
            raise ValueError(f"Week must be 0<=week<=3, got {self.week}")

    def _get_offset_day(self, other: datetime) -> int:
        """
        Find the day in the same month as other that has the same
        weekday as self.weekday and is the self.week'th such day in the month.

        Parameters
        ----------
        other : datetime

        Returns
        -------
        day : int
        """
        mstart = datetime(other.year, other.month, 1)
        wday = mstart.weekday()
        shift_days = (self.weekday - wday) % 7
        return 1 + shift_days + self.week * 7

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            raise ValueError(f"Prefix {repr(cls._prefix)} requires a suffix.")
        # TODO: handle n here...
        # only one digit weeks (1 --> week 0, 2 --> week 1, etc.)
        week = int(suffix[0]) - 1
        weekday = ccalendar.weekday_to_int[suffix[1:]]
        return cls(week=week, weekday=weekday)


class LastWeekOfMonth(liboffsets.WeekOfMonthMixin):
    """
    Describes monthly dates in last week of month like "the last Tuesday of
    each month".

    Parameters
    ----------
    n : int, default 1
    weekday : int {0, 1, ..., 6}, default 0
        A specific integer for the day of the week.

        - 0 is Monday
        - 1 is Tuesday
        - 2 is Wednesday
        - 3 is Thursday
        - 4 is Friday
        - 5 is Saturday
        - 6 is Sunday.
    """

    _prefix = "LWOM"
    _attributes = frozenset(["n", "normalize", "weekday"])

    def __init__(self, n=1, normalize=False, weekday=0):
        liboffsets.WeekOfMonthMixin.__init__(self, n, normalize, weekday)

        if self.n == 0:
            raise ValueError("N cannot be 0")
        object.__setattr__(self, "week", -1)

    def _get_offset_day(self, other: datetime) -> int:
        """
        Find the day in the same month as other that has the same
        weekday as self.weekday and is the last such day in the month.

        Parameters
        ----------
        other: datetime

        Returns
        -------
        day: int
        """
        dim = ccalendar.get_days_in_month(other.year, other.month)
        mend = datetime(other.year, other.month, dim)
        wday = mend.weekday()
        shift_days = (wday - self.weekday) % 7
        return dim - shift_days

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            raise ValueError(f"Prefix {repr(cls._prefix)} requires a suffix.")
        # TODO: handle n here...
        weekday = ccalendar.weekday_to_int[suffix]
        return cls(weekday=weekday)


# ---------------------------------------------------------------------
# Special Offset Classes


class FY5253(liboffsets.FY5253Mixin):
    """
    Describes 52-53 week fiscal year. This is also known as a 4-4-5 calendar.

    It is used by companies that desire that their
    fiscal year always end on the same day of the week.

    It is a method of managing accounting periods.
    It is a common calendar structure for some industries,
    such as retail, manufacturing and parking industry.

    For more information see:
    https://en.wikipedia.org/wiki/4-4-5_calendar

    The year may either:

    - end on the last X day of the Y month.
    - end on the last X day closest to the last day of the Y month.

    X is a specific day of the week.
    Y is a certain month of the year

    Parameters
    ----------
    n : int
    weekday : int {0, 1, ..., 6}, default 0
        A specific integer for the day of the week.

        - 0 is Monday
        - 1 is Tuesday
        - 2 is Wednesday
        - 3 is Thursday
        - 4 is Friday
        - 5 is Saturday
        - 6 is Sunday.

    startingMonth : int {1, 2, ... 12}, default 1
        The month in which the fiscal year ends.

    variation : str, default "nearest"
        Method of employing 4-4-5 calendar.

        There are two options:

        - "nearest" means year end is **weekday** closest to last day of month in year.
        - "last" means year end is final **weekday** of the final month in fiscal year.
    """

    _prefix = "RE"
    _attributes = frozenset(["weekday", "startingMonth", "variation"])

    def __reduce__(self):
        tup = (self.n, self.normalize, self.weekday, self.startingMonth, self.variation)
        return type(self), tup

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        dt = datetime(dt.year, dt.month, dt.day)
        year_end = self.get_year_end(dt)

        if self.variation == "nearest":
            # We have to check the year end of "this" cal year AND the previous
            return year_end == dt or self.get_year_end(shift_month(dt, -1, None)) == dt
        else:
            return year_end == dt

    @apply_wraps
    def apply(self, other):
        norm = Timestamp(other).normalize()

        n = self.n
        prev_year = self.get_year_end(datetime(other.year - 1, self.startingMonth, 1))
        cur_year = self.get_year_end(datetime(other.year, self.startingMonth, 1))
        next_year = self.get_year_end(datetime(other.year + 1, self.startingMonth, 1))

        prev_year = conversion.localize_pydatetime(prev_year, other.tzinfo)
        cur_year = conversion.localize_pydatetime(cur_year, other.tzinfo)
        next_year = conversion.localize_pydatetime(next_year, other.tzinfo)

        # Note: next_year.year == other.year + 1, so we will always
        # have other < next_year
        if norm == prev_year:
            n -= 1
        elif norm == cur_year:
            pass
        elif n > 0:
            if norm < prev_year:
                n -= 2
            elif prev_year < norm < cur_year:
                n -= 1
            elif cur_year < norm < next_year:
                pass
        else:
            if cur_year < norm < next_year:
                n += 1
            elif prev_year < norm < cur_year:
                pass
            elif (
                norm.year == prev_year.year
                and norm < prev_year
                and prev_year - norm <= timedelta(6)
            ):
                # GH#14774, error when next_year.year == cur_year.year
                # e.g. prev_year == datetime(2004, 1, 3),
                # other == datetime(2004, 1, 1)
                n -= 1
            else:
                assert False

        shifted = datetime(other.year + n, self.startingMonth, 1)
        result = self.get_year_end(shifted)
        result = datetime(
            result.year,
            result.month,
            result.day,
            other.hour,
            other.minute,
            other.second,
            other.microsecond,
        )
        return result

    def get_year_end(self, dt):
        assert dt.tzinfo is None

        dim = ccalendar.get_days_in_month(dt.year, self.startingMonth)
        target_date = datetime(dt.year, self.startingMonth, dim)
        wkday_diff = self.weekday - target_date.weekday()
        if wkday_diff == 0:
            # year_end is the same for "last" and "nearest" cases
            return target_date

        if self.variation == "last":
            days_forward = (wkday_diff % 7) - 7

            # days_forward is always negative, so we always end up
            # in the same year as dt
            return target_date + timedelta(days=days_forward)
        else:
            # variation == "nearest":
            days_forward = wkday_diff % 7
            if days_forward <= 3:
                # The upcoming self.weekday is closer than the previous one
                return target_date + timedelta(days_forward)
            else:
                # The previous self.weekday is closer than the upcoming one
                return target_date + timedelta(days_forward - 7)

    @classmethod
    def _parse_suffix(cls, varion_code, startingMonth_code, weekday_code):
        if varion_code == "N":
            variation = "nearest"
        elif varion_code == "L":
            variation = "last"
        else:
            raise ValueError(f"Unable to parse varion_code: {varion_code}")

        startingMonth = ccalendar.MONTH_TO_CAL_NUM[startingMonth_code]
        weekday = ccalendar.weekday_to_int[weekday_code]

        return {
            "weekday": weekday,
            "startingMonth": startingMonth,
            "variation": variation,
        }

    @classmethod
    def _from_name(cls, *args):
        return cls(**cls._parse_suffix(*args))


class FY5253Quarter(liboffsets.FY5253Mixin):
    """
    DateOffset increments between business quarter dates
    for 52-53 week fiscal year (also known as a 4-4-5 calendar).

    It is used by companies that desire that their
    fiscal year always end on the same day of the week.

    It is a method of managing accounting periods.
    It is a common calendar structure for some industries,
    such as retail, manufacturing and parking industry.

    For more information see:
    https://en.wikipedia.org/wiki/4-4-5_calendar

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
    weekday : int {0, 1, ..., 6}, default 0
        A specific integer for the day of the week.

        - 0 is Monday
        - 1 is Tuesday
        - 2 is Wednesday
        - 3 is Thursday
        - 4 is Friday
        - 5 is Saturday
        - 6 is Sunday.

    startingMonth : int {1, 2, ..., 12}, default 1
        The month in which fiscal years end.

    qtr_with_extra_week : int {1, 2, 3, 4}, default 1
        The quarter number that has the leap or 14 week when needed.

    variation : str, default "nearest"
        Method of employing 4-4-5 calendar.

        There are two options:

        - "nearest" means year end is **weekday** closest to last day of month in year.
        - "last" means year end is final **weekday** of the final month in fiscal year.
    """

    _prefix = "REQ"
    _attributes = frozenset(
        ["weekday", "startingMonth", "qtr_with_extra_week", "variation"]
    )

    def __init__(
        self,
        n=1,
        normalize=False,
        weekday=0,
        startingMonth=1,
        qtr_with_extra_week=1,
        variation="nearest",
    ):
        liboffsets.FY5253Mixin.__init__(
            self, n, normalize, weekday, startingMonth, variation
        )
        object.__setattr__(self, "qtr_with_extra_week", qtr_with_extra_week)

    def __reduce__(self):
        tup = (
            self.n,
            self.normalize,
            self.weekday,
            self.startingMonth,
            self.qtr_with_extra_week,
            self.variation,
        )
        return type(self), tup

    @cache_readonly
    def _offset(self):
        return FY5253(
            startingMonth=self.startingMonth,
            weekday=self.weekday,
            variation=self.variation,
        )

    def _rollback_to_year(self, other):
        """
        Roll `other` back to the most recent date that was on a fiscal year
        end.

        Return the date of that year-end, the number of full quarters
        elapsed between that year-end and other, and the remaining Timedelta
        since the most recent quarter-end.

        Parameters
        ----------
        other : datetime or Timestamp

        Returns
        -------
        tuple of
        prev_year_end : Timestamp giving most recent fiscal year end
        num_qtrs : int
        tdelta : Timedelta
        """
        num_qtrs = 0

        norm = Timestamp(other).tz_localize(None)
        start = self._offset.rollback(norm)
        # Note: start <= norm and self._offset.is_on_offset(start)

        if start < norm:
            # roll adjustment
            qtr_lens = self.get_weeks(norm)

            # check that qtr_lens is consistent with self._offset addition
            end = liboffsets.shift_day(start, days=7 * sum(qtr_lens))
            assert self._offset.is_on_offset(end), (start, end, qtr_lens)

            tdelta = norm - start
            for qlen in qtr_lens:
                if qlen * 7 <= tdelta.days:
                    num_qtrs += 1
                    tdelta -= Timedelta(days=qlen * 7)
                else:
                    break
        else:
            tdelta = Timedelta(0)

        # Note: we always have tdelta.value >= 0
        return start, num_qtrs, tdelta

    @apply_wraps
    def apply(self, other):
        # Note: self.n == 0 is not allowed.
        n = self.n

        prev_year_end, num_qtrs, tdelta = self._rollback_to_year(other)
        res = prev_year_end
        n += num_qtrs
        if self.n <= 0 and tdelta.value > 0:
            n += 1

        # Possible speedup by handling years first.
        years = n // 4
        if years:
            res += self._offset * years
            n -= years * 4

        # Add an extra day to make *sure* we are getting the quarter lengths
        # for the upcoming year, not the previous year
        qtr_lens = self.get_weeks(res + Timedelta(days=1))

        # Note: we always have 0 <= n < 4
        weeks = sum(qtr_lens[:n])
        if weeks:
            res = liboffsets.shift_day(res, days=weeks * 7)

        return res

    def get_weeks(self, dt):
        ret = [13] * 4

        year_has_extra_week = self.year_has_extra_week(dt)

        if year_has_extra_week:
            ret[self.qtr_with_extra_week - 1] = 14

        return ret

    def year_has_extra_week(self, dt: datetime) -> bool:
        # Avoid round-down errors --> normalize to get
        # e.g. '370D' instead of '360D23H'
        norm = Timestamp(dt).normalize().tz_localize(None)

        next_year_end = self._offset.rollforward(norm)
        prev_year_end = norm - self._offset
        weeks_in_year = (next_year_end - prev_year_end).days / 7
        assert weeks_in_year in [52, 53], weeks_in_year
        return weeks_in_year == 53

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        if self._offset.is_on_offset(dt):
            return True

        next_year_end = dt - self._offset

        qtr_lens = self.get_weeks(dt)

        current = next_year_end
        for qtr_len in qtr_lens:
            current = liboffsets.shift_day(current, days=qtr_len * 7)
            if dt == current:
                return True
        return False

    @property
    def rule_code(self) -> str:
        suffix = liboffsets.FY5253Mixin.rule_code.__get__(self)
        qtr = self.qtr_with_extra_week
        return f"{suffix}-{qtr}"

    @classmethod
    def _from_name(cls, *args):
        return cls(
            **dict(FY5253._parse_suffix(*args[:-1]), qtr_with_extra_week=int(args[-1]))
        )


class Easter(SingleConstructorOffset):
    """
    DateOffset for the Easter holiday using logic defined in dateutil.

    Right now uses the revised method which is valid in years 1583-4099.
    """

    @apply_wraps
    def apply(self, other):
        current_easter = easter(other.year)
        current_easter = datetime(
            current_easter.year, current_easter.month, current_easter.day
        )
        current_easter = conversion.localize_pydatetime(current_easter, other.tzinfo)

        n = self.n
        if n >= 0 and other < current_easter:
            n -= 1
        elif n < 0 and other > current_easter:
            n += 1
        # TODO: Why does this handle the 0 case the opposite of others?

        # NOTE: easter returns a datetime.date so we have to convert to type of
        # other
        new = easter(other.year + n)
        new = datetime(
            new.year,
            new.month,
            new.day,
            other.hour,
            other.minute,
            other.second,
            other.microsecond,
        )
        return new

    def is_on_offset(self, dt: datetime) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        return date(dt.year, dt.month, dt.day) == easter(dt.year)


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
