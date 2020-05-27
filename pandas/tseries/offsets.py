from datetime import datetime, timedelta

import numpy as np

from pandas._libs.tslibs import offsets as liboffsets
from pandas._libs.tslibs.offsets import (  # noqa:F401
    FY5253,
    ApplyTypeError,
    BaseOffset,
    BQuarterBegin,
    BQuarterEnd,
    BusinessDay,
    BusinessHour,
    BusinessMixin,
    BusinessMonthBegin,
    BusinessMonthEnd,
    BYearBegin,
    BYearEnd,
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


class CustomBusinessDay(BusinessDay):
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
        self._init_custom(weekmask, holidays, calendar)

    def __setstate__(self, state):
        self.holidays = state.pop("holidays")
        self.weekmask = state.pop("weekmask")
        super().__setstate__(state)

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


class CustomBusinessHour(BusinessHour):
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
        self._init_custom(weekmask, holidays, calendar)

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
class _CustomBusinessMonth(BusinessMixin, liboffsets.MonthOffset):
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
        self._init_custom(weekmask, holidays, calendar)

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
