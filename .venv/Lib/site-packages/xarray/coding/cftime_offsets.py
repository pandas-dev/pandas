"""Time offset classes for use with cftime.datetime objects"""

# The offset classes and mechanisms for generating time ranges defined in
# this module were copied/adapted from those defined in pandas.  See in
# particular the objects and methods defined in pandas.tseries.offsets
# and pandas.core.indexes.datetimes.

# For reference, here is a copy of the pandas copyright notice:

# (c) 2011-2012, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2008-2011 AQR Capital Management, LLC
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the copyright holder nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from __future__ import annotations

import re
import warnings
from collections.abc import Mapping
from datetime import datetime, timedelta
from functools import partial
from typing import TYPE_CHECKING, ClassVar, Literal, TypeVar, get_args

import numpy as np
import pandas as pd
from packaging.version import Version

from xarray.coding.cftimeindex import CFTimeIndex
from xarray.coding.times import (
    _is_standard_calendar,
    _parse_iso8601,
    _should_cftime_be_used,
    convert_time_or_go_back,
    format_cftime_datetime,
)
from xarray.compat.pdcompat import (
    count_not_none,
    default_precision_timestamp,
)
from xarray.core.common import _contains_datetime_like_objects, is_np_datetime_like
from xarray.core.types import InclusiveOptions
from xarray.core.utils import attempt_import, emit_user_level_warning

if TYPE_CHECKING:
    from xarray.core.types import (
        PDDatetimeUnitOptions,
        Self,
        TypeAlias,
    )


DayOption: TypeAlias = Literal["start", "end"]
T_FreqStr = TypeVar("T_FreqStr", str, None)


def get_date_type(calendar, use_cftime=True):
    """Return the cftime date type for a given calendar name."""
    if TYPE_CHECKING:
        import cftime
    else:
        cftime = attempt_import("cftime")

    if _is_standard_calendar(calendar) and not use_cftime:
        return default_precision_timestamp

    calendars = {
        "noleap": cftime.DatetimeNoLeap,
        "360_day": cftime.Datetime360Day,
        "365_day": cftime.DatetimeNoLeap,
        "366_day": cftime.DatetimeAllLeap,
        "gregorian": cftime.DatetimeGregorian,
        "proleptic_gregorian": cftime.DatetimeProlepticGregorian,
        "julian": cftime.DatetimeJulian,
        "all_leap": cftime.DatetimeAllLeap,
        "standard": cftime.DatetimeGregorian,
    }
    return calendars[calendar]


class BaseCFTimeOffset:
    _freq: ClassVar[str | None] = None
    _day_option: ClassVar[DayOption | None] = None
    n: int

    def __init__(self, n: int = 1) -> None:
        if not isinstance(n, int):
            raise TypeError(
                "The provided multiple 'n' must be an integer. "
                f"Instead a value of type {type(n)!r} was provided."
            )
        self.n = n

    def rule_code(self) -> str | None:
        return self._freq

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BaseCFTimeOffset):
            return NotImplemented
        return self.n == other.n and self.rule_code() == other.rule_code()

    def __ne__(self, other: object) -> bool:
        return not self == other

    def __add__(self, other):
        return self.__apply__(other)

    def __sub__(self, other):
        if TYPE_CHECKING:
            import cftime
        else:
            cftime = attempt_import("cftime")

        if isinstance(other, cftime.datetime):
            raise TypeError("Cannot subtract a cftime.datetime from a time offset.")
        elif type(other) is type(self):
            return type(self)(self.n - other.n)
        else:
            return NotImplemented

    def __mul__(self, other: int) -> Self:
        if not isinstance(other, int):
            return NotImplemented
        return type(self)(n=other * self.n)

    def __neg__(self) -> Self:
        return self * -1

    def __rmul__(self, other):
        return self.__mul__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        if isinstance(other, BaseCFTimeOffset) and type(self) is not type(other):
            raise TypeError("Cannot subtract cftime offsets of differing types")
        return -self + other

    def __apply__(self, other):
        return NotImplemented

    def onOffset(self, date) -> bool:
        """Check if the given date is in the set of possible dates created
        using a length-one version of this offset class."""
        test_date = (self + date) - self
        return date == test_date

    def rollforward(self, date):
        if self.onOffset(date):
            return date
        else:
            return date + type(self)()

    def rollback(self, date):
        if self.onOffset(date):
            return date
        else:
            return date - type(self)()

    def __str__(self):
        return f"<{type(self).__name__}: n={self.n}>"

    def __repr__(self):
        return str(self)

    def _get_offset_day(self, other):
        # subclass must implement `_day_option`; calling from the base class
        # will raise NotImplementedError.
        return _get_day_of_month(other, self._day_option)


class Tick(BaseCFTimeOffset):
    # analogous https://github.com/pandas-dev/pandas/blob/ccb25ab1d24c4fb9691270706a59c8d319750870/pandas/_libs/tslibs/offsets.pyx#L806

    def _next_higher_resolution(self) -> Tick:
        self_type = type(self)
        if self_type is Day:
            return Hour(self.n * 24)
        if self_type is Hour:
            return Minute(self.n * 60)
        if self_type is Minute:
            return Second(self.n * 60)
        if self_type is Second:
            return Millisecond(self.n * 1000)
        if self_type is Millisecond:
            return Microsecond(self.n * 1000)
        raise ValueError("Could not convert to integer offset at any resolution")

    def __mul__(self, other: int | float) -> Tick:
        if not isinstance(other, int | float):
            return NotImplemented
        if isinstance(other, float):
            n = other * self.n
            # If the new `n` is an integer, we can represent it using the
            #  same BaseCFTimeOffset subclass as self, otherwise we need to move up
            #  to a higher-resolution subclass
            if np.isclose(n % 1, 0):
                return type(self)(int(n))

            new_self = self._next_higher_resolution()
            return new_self * other
        return type(self)(n=other * self.n)

    def as_timedelta(self) -> timedelta:
        """All Tick subclasses must implement an as_timedelta method."""
        raise NotImplementedError


def _get_day_of_month(other, day_option: DayOption) -> int:
    """Find the day in `other`'s month that satisfies a BaseCFTimeOffset's
    onOffset policy, as described by the `day_option` argument.

    Parameters
    ----------
    other : cftime.datetime
    day_option : 'start', 'end'
        'start': returns 1
        'end': returns last day of the month

    Returns
    -------
    day_of_month : int

    """

    if day_option == "start":
        return 1
    elif day_option == "end":
        return other.daysinmonth
    elif day_option is None:
        # Note: unlike `_shift_month`, _get_day_of_month does not
        # allow day_option = None
        raise NotImplementedError()
    raise ValueError(day_option)


def _adjust_n_months(other_day, n, reference_day):
    """Adjust the number of times a monthly offset is applied based
    on the day of a given date, and the reference day provided.
    """
    if n > 0 and other_day < reference_day:
        n = n - 1
    elif n <= 0 and other_day > reference_day:
        n = n + 1
    return n


def _adjust_n_years(other, n, month, reference_day):
    """Adjust the number of times an annual offset is applied based on
    another date, and the reference day provided"""
    if n > 0:
        if other.month < month or (other.month == month and other.day < reference_day):
            n -= 1
    elif other.month > month or (other.month == month and other.day > reference_day):
        n += 1
    return n


def _shift_month(date, months, day_option: DayOption = "start"):
    """Shift the date to a month start or end a given number of months away."""
    _ = attempt_import("cftime")

    has_year_zero = date.has_year_zero
    year = date.year + (date.month + months) // 12
    month = (date.month + months) % 12

    if month == 0:
        month = 12
        year -= 1

    if not has_year_zero:
        if date.year < 0 <= year:
            year += 1
        elif year <= 0 < date.year:
            year -= 1

    # Silence warnings associated with generating dates with years < 1.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="this date/calendar/year zero")

        if day_option == "start":
            day = 1
        elif day_option == "end":
            reference = type(date)(year, month, 1, has_year_zero=has_year_zero)
            day = reference.daysinmonth
        else:
            raise ValueError(day_option)
        return date.replace(year=year, month=month, day=day)


def roll_qtrday(
    other, n: int, month: int, day_option: DayOption, modby: int = 3
) -> int:
    """Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : cftime.datetime
    n : number of periods to increment, before adjusting for rolling
    month : int reference month giving the first month of the year
    day_option : 'start', 'end'
        The convention to use in finding the day in a given month against
        which to compare for rollforward/rollbackward decisions.
    modby : int 3 for quarters, 12 for years

    Returns
    -------
    n : int number of periods to increment

    See Also
    --------
    _get_day_of_month : Find the day in a month provided an offset.
    """

    months_since = other.month % modby - month % modby

    if n > 0:
        if months_since < 0 or (
            months_since == 0 and other.day < _get_day_of_month(other, day_option)
        ):
            # pretend to roll back if on same month but
            # before compare_day
            n -= 1
    elif months_since > 0 or (
        months_since == 0 and other.day > _get_day_of_month(other, day_option)
    ):
        # make sure to roll forward, so negate
        n += 1
    return n


def _validate_month(month: int | None, default_month: int) -> int:
    result_month = default_month if month is None else month
    if not isinstance(result_month, int):
        raise TypeError(
            "'self.month' must be an integer value between 1 "
            "and 12.  Instead, it was set to a value of "
            f"{result_month!r}"
        )
    elif not (1 <= result_month <= 12):
        raise ValueError(
            "'self.month' must be an integer value between 1 "
            "and 12.  Instead, it was set to a value of "
            f"{result_month!r}"
        )
    return result_month


class MonthBegin(BaseCFTimeOffset):
    _freq = "MS"

    def __apply__(self, other):
        n = _adjust_n_months(other.day, self.n, 1)
        return _shift_month(other, n, "start")

    def onOffset(self, date) -> bool:
        """Check if the given date is in the set of possible dates created
        using a length-one version of this offset class."""
        return date.day == 1


class MonthEnd(BaseCFTimeOffset):
    _freq = "ME"

    def __apply__(self, other):
        n = _adjust_n_months(other.day, self.n, other.daysinmonth)
        return _shift_month(other, n, "end")

    def onOffset(self, date) -> bool:
        """Check if the given date is in the set of possible dates created
        using a length-one version of this offset class."""
        return date.day == date.daysinmonth


_MONTH_ABBREVIATIONS = {
    1: "JAN",
    2: "FEB",
    3: "MAR",
    4: "APR",
    5: "MAY",
    6: "JUN",
    7: "JUL",
    8: "AUG",
    9: "SEP",
    10: "OCT",
    11: "NOV",
    12: "DEC",
}


class QuarterOffset(BaseCFTimeOffset):
    """Quarter representation copied off of pandas/tseries/offsets.py"""

    _default_month: ClassVar[int]
    month: int

    def __init__(self, n: int = 1, month: int | None = None) -> None:
        BaseCFTimeOffset.__init__(self, n)
        self.month = _validate_month(month, self._default_month)

    def __apply__(self, other):
        # months_since: find the calendar quarter containing other.month,
        # e.g. if other.month == 8, the calendar quarter is [Jul, Aug, Sep].
        # Then find the month in that quarter containing an onOffset date for
        # self.  `months_since` is the number of months to shift other.month
        # to get to this on-offset month.
        months_since = other.month % 3 - self.month % 3
        qtrs = roll_qtrday(
            other, self.n, self.month, day_option=self._day_option, modby=3
        )
        months = qtrs * 3 - months_since
        return _shift_month(other, months, self._day_option)

    def onOffset(self, date) -> bool:
        """Check if the given date is in the set of possible dates created
        using a length-one version of this offset class."""
        mod_month = (date.month - self.month) % 3
        return mod_month == 0 and date.day == self._get_offset_day(date)

    def __sub__(self, other: Self) -> Self:
        if TYPE_CHECKING:
            import cftime
        else:
            cftime = attempt_import("cftime")

        if isinstance(other, cftime.datetime):
            raise TypeError("Cannot subtract cftime.datetime from offset.")
        if type(other) is type(self) and other.month == self.month:
            return type(self)(self.n - other.n, month=self.month)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, float):
            return NotImplemented
        return type(self)(n=other * self.n, month=self.month)

    def rule_code(self) -> str:
        return f"{self._freq}-{_MONTH_ABBREVIATIONS[self.month]}"

    def __str__(self):
        return f"<{type(self).__name__}: n={self.n}, month={self.month}>"


class QuarterBegin(QuarterOffset):
    # When converting a string to an offset, pandas converts
    # 'QS' to a QuarterBegin offset starting in the month of
    # January.  When creating a QuarterBegin offset directly
    # from the constructor, however, the default month is March.
    # We follow that behavior here.
    _default_month = 3
    _freq = "QS"
    _day_option = "start"

    def rollforward(self, date):
        """Roll date forward to nearest start of quarter"""
        if self.onOffset(date):
            return date
        else:
            return date + QuarterBegin(month=self.month)

    def rollback(self, date):
        """Roll date backward to nearest start of quarter"""
        if self.onOffset(date):
            return date
        else:
            return date - QuarterBegin(month=self.month)


class QuarterEnd(QuarterOffset):
    # When converting a string to an offset, pandas converts
    # 'Q' to a QuarterEnd offset starting in the month of
    # December.  When creating a QuarterEnd offset directly
    # from the constructor, however, the default month is March.
    # We follow that behavior here.
    _default_month = 3
    _freq = "QE"
    _day_option = "end"

    def rollforward(self, date):
        """Roll date forward to nearest end of quarter"""
        if self.onOffset(date):
            return date
        else:
            return date + QuarterEnd(month=self.month)

    def rollback(self, date):
        """Roll date backward to nearest end of quarter"""
        if self.onOffset(date):
            return date
        else:
            return date - QuarterEnd(month=self.month)


class YearOffset(BaseCFTimeOffset):
    _default_month: ClassVar[int]
    month: int

    def __init__(self, n: int = 1, month: int | None = None) -> None:
        BaseCFTimeOffset.__init__(self, n)
        self.month = _validate_month(month, self._default_month)

    def __apply__(self, other):
        reference_day = _get_day_of_month(other, self._day_option)
        years = _adjust_n_years(other, self.n, self.month, reference_day)
        months = years * 12 + (self.month - other.month)
        return _shift_month(other, months, self._day_option)

    def __sub__(self, other):
        if TYPE_CHECKING:
            import cftime
        else:
            cftime = attempt_import("cftime")

        if isinstance(other, cftime.datetime):
            raise TypeError("Cannot subtract cftime.datetime from offset.")
        elif type(other) is type(self) and other.month == self.month:
            return type(self)(self.n - other.n, month=self.month)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, float):
            return NotImplemented
        return type(self)(n=other * self.n, month=self.month)

    def rule_code(self) -> str:
        return f"{self._freq}-{_MONTH_ABBREVIATIONS[self.month]}"

    def __str__(self) -> str:
        return f"<{type(self).__name__}: n={self.n}, month={self.month}>"


class YearBegin(YearOffset):
    _freq = "YS"
    _day_option = "start"
    _default_month = 1

    def onOffset(self, date) -> bool:
        """Check if the given date is in the set of possible dates created
        using a length-one version of this offset class."""
        return date.day == 1 and date.month == self.month

    def rollforward(self, date):
        """Roll date forward to nearest start of year"""
        if self.onOffset(date):
            return date
        else:
            return date + YearBegin(month=self.month)

    def rollback(self, date):
        """Roll date backward to nearest start of year"""
        if self.onOffset(date):
            return date
        else:
            return date - YearBegin(month=self.month)


class YearEnd(YearOffset):
    _freq = "YE"
    _day_option = "end"
    _default_month = 12

    def onOffset(self, date) -> bool:
        """Check if the given date is in the set of possible dates created
        using a length-one version of this offset class."""
        return date.day == date.daysinmonth and date.month == self.month

    def rollforward(self, date):
        """Roll date forward to nearest end of year"""
        if self.onOffset(date):
            return date
        else:
            return date + YearEnd(month=self.month)

    def rollback(self, date):
        """Roll date backward to nearest end of year"""
        if self.onOffset(date):
            return date
        else:
            return date - YearEnd(month=self.month)


class Day(Tick):
    _freq = "D"

    def as_timedelta(self) -> timedelta:
        return timedelta(days=self.n)

    def __apply__(self, other):
        return other + self.as_timedelta()


class Hour(Tick):
    _freq = "h"

    def as_timedelta(self) -> timedelta:
        return timedelta(hours=self.n)

    def __apply__(self, other):
        return other + self.as_timedelta()


class Minute(Tick):
    _freq = "min"

    def as_timedelta(self) -> timedelta:
        return timedelta(minutes=self.n)

    def __apply__(self, other):
        return other + self.as_timedelta()


class Second(Tick):
    _freq = "s"

    def as_timedelta(self) -> timedelta:
        return timedelta(seconds=self.n)

    def __apply__(self, other):
        return other + self.as_timedelta()


class Millisecond(Tick):
    _freq = "ms"

    def as_timedelta(self) -> timedelta:
        return timedelta(milliseconds=self.n)

    def __apply__(self, other):
        return other + self.as_timedelta()


class Microsecond(Tick):
    _freq = "us"

    def as_timedelta(self) -> timedelta:
        return timedelta(microseconds=self.n)

    def __apply__(self, other):
        return other + self.as_timedelta()


def _generate_anchored_offsets(
    base_freq: str, offset: type[YearOffset | QuarterOffset]
) -> dict[str, type[BaseCFTimeOffset]]:
    offsets: dict[str, type[BaseCFTimeOffset]] = {}
    for month, abbreviation in _MONTH_ABBREVIATIONS.items():
        anchored_freq = f"{base_freq}-{abbreviation}"
        offsets[anchored_freq] = partial(offset, month=month)  # type: ignore[assignment]
    return offsets


_FREQUENCIES: Mapping[str, type[BaseCFTimeOffset]] = {
    "A": YearEnd,
    "AS": YearBegin,
    "Y": YearEnd,
    "YE": YearEnd,
    "YS": YearBegin,
    "Q": partial(QuarterEnd, month=12),  # type: ignore[dict-item]
    "QE": partial(QuarterEnd, month=12),  # type: ignore[dict-item]
    "QS": partial(QuarterBegin, month=1),  # type: ignore[dict-item]
    "M": MonthEnd,
    "ME": MonthEnd,
    "MS": MonthBegin,
    "D": Day,
    "H": Hour,
    "h": Hour,
    "T": Minute,
    "min": Minute,
    "S": Second,
    "s": Second,
    "L": Millisecond,
    "ms": Millisecond,
    "U": Microsecond,
    "us": Microsecond,
    **_generate_anchored_offsets("AS", YearBegin),
    **_generate_anchored_offsets("A", YearEnd),
    **_generate_anchored_offsets("YS", YearBegin),
    **_generate_anchored_offsets("Y", YearEnd),
    **_generate_anchored_offsets("YE", YearEnd),
    **_generate_anchored_offsets("QS", QuarterBegin),
    **_generate_anchored_offsets("Q", QuarterEnd),
    **_generate_anchored_offsets("QE", QuarterEnd),
}


_FREQUENCY_CONDITION = "|".join(_FREQUENCIES.keys())
_PATTERN = rf"^((?P<multiple>[+-]?\d+)|())(?P<freq>({_FREQUENCY_CONDITION}))$"


# pandas defines these offsets as "Tick" objects, which for instance have
# distinct behavior from monthly or longer frequencies in resample.
CFTIME_TICKS = (Day, Hour, Minute, Second)


def _generate_anchored_deprecated_frequencies(
    deprecated: str, recommended: str
) -> dict[str, str]:
    pairs = {}
    for abbreviation in _MONTH_ABBREVIATIONS.values():
        anchored_deprecated = f"{deprecated}-{abbreviation}"
        anchored_recommended = f"{recommended}-{abbreviation}"
        pairs[anchored_deprecated] = anchored_recommended
    return pairs


_DEPRECATED_FREQUENCIES: dict[str, str] = {
    "A": "YE",
    "Y": "YE",
    "AS": "YS",
    "Q": "QE",
    "M": "ME",
    "H": "h",
    "T": "min",
    "S": "s",
    "L": "ms",
    "U": "us",
    **_generate_anchored_deprecated_frequencies("A", "YE"),
    **_generate_anchored_deprecated_frequencies("Y", "YE"),
    **_generate_anchored_deprecated_frequencies("AS", "YS"),
    **_generate_anchored_deprecated_frequencies("Q", "QE"),
}


_DEPRECATION_MESSAGE = (
    "{deprecated_freq!r} is deprecated and will be removed in a future "
    "version. Please use {recommended_freq!r} instead of "
    "{deprecated_freq!r}."
)


def _emit_freq_deprecation_warning(deprecated_freq):
    recommended_freq = _DEPRECATED_FREQUENCIES[deprecated_freq]
    message = _DEPRECATION_MESSAGE.format(
        deprecated_freq=deprecated_freq, recommended_freq=recommended_freq
    )
    emit_user_level_warning(message, FutureWarning)


def to_offset(
    freq: BaseCFTimeOffset | str | timedelta | pd.Timedelta | pd.DateOffset,
    warn: bool = True,
) -> BaseCFTimeOffset:
    """Convert a frequency string to the appropriate subclass of
    BaseCFTimeOffset."""
    if isinstance(freq, BaseCFTimeOffset):
        return freq
    if isinstance(freq, timedelta | pd.Timedelta):
        return delta_to_tick(freq)
    if isinstance(freq, pd.DateOffset):
        freq = _legacy_to_new_freq(freq.freqstr)

    match = re.match(_PATTERN, freq)
    if match is None:
        raise ValueError("Invalid frequency string provided")
    freq_data = match.groupdict()

    freq = freq_data["freq"]
    if warn and freq in _DEPRECATED_FREQUENCIES:
        _emit_freq_deprecation_warning(freq)
    multiples = freq_data["multiple"]
    multiples = 1 if multiples is None else int(multiples)
    return _FREQUENCIES[freq](n=multiples)


def delta_to_tick(delta: timedelta | pd.Timedelta) -> Tick:
    """Adapted from pandas.tslib.delta_to_tick"""
    if isinstance(delta, pd.Timedelta) and delta.nanoseconds != 0:
        # pandas.Timedelta has nanoseconds, but these are not supported
        raise ValueError(
            "Unable to convert 'pandas.Timedelta' object with non-zero "
            "nanoseconds to 'CFTimeOffset' object"
        )
    if delta.microseconds == 0:
        if delta.seconds == 0:
            return Day(n=delta.days)
        else:
            seconds = delta.days * 86400 + delta.seconds
            if seconds % 3600 == 0:
                return Hour(n=seconds // 3600)
            elif seconds % 60 == 0:
                return Minute(n=seconds // 60)
            else:
                return Second(n=seconds)
    # Regardless of the days and seconds this will always be a Millisecond
    # or Microsecond object
    elif delta.microseconds % 1_000 == 0:
        return Millisecond(n=delta.microseconds // 1_000)
    else:
        return Microsecond(n=delta.microseconds)


def to_cftime_datetime(date_str_or_date, calendar=None):
    if TYPE_CHECKING:
        import cftime
    else:
        cftime = attempt_import("cftime")

    if isinstance(date_str_or_date, str):
        if calendar is None:
            raise ValueError(
                "If converting a string to a cftime.datetime object, "
                "a calendar type must be provided"
            )
        date, _ = _parse_iso8601(get_date_type(calendar), date_str_or_date)
        return date
    elif isinstance(date_str_or_date, cftime.datetime):
        return date_str_or_date
    elif isinstance(date_str_or_date, datetime | pd.Timestamp):
        return cftime.DatetimeProlepticGregorian(*date_str_or_date.timetuple())
    else:
        raise TypeError(
            "date_str_or_date must be a string or a "
            "subclass of cftime.datetime. Instead got "
            f"{date_str_or_date!r}."
        )


def normalize_date(date):
    """Round datetime down to midnight."""
    return date.replace(hour=0, minute=0, second=0, microsecond=0)


def _get_normalized_cfdate(date, calendar, normalize):
    """convert to cf datetime and round down to midnight if normalize."""
    if date is None:
        return date

    cf_date = to_cftime_datetime(date, calendar)

    return normalize_date(cf_date) if normalize else cf_date


def _generate_linear_date_range(start, end, periods):
    """Generate an equally-spaced sequence of cftime.datetime objects between
    and including two dates (whose length equals the number of periods)."""
    if TYPE_CHECKING:
        import cftime
    else:
        cftime = attempt_import("cftime")

    total_seconds = (end - start).total_seconds()
    values = np.linspace(0.0, total_seconds, periods, endpoint=True)
    units = f"seconds since {format_cftime_datetime(start)}"
    calendar = start.calendar
    return cftime.num2date(
        values, units=units, calendar=calendar, only_use_cftime_datetimes=True
    )


def _generate_linear_date_range_with_freq(start, end, periods, freq):
    """Generate a regular range of cftime.datetime objects with a
    given frequency.

    Adapted from pandas.tseries.offsets.generate_range (now at
    pandas.core.arrays.datetimes._generate_range).

    Parameters
    ----------
    start : cftime.datetime, or None
        Start of range
    end : cftime.datetime, or None
        End of range
    periods : int, or None
        Number of elements in the sequence
    freq: str
        Step size between cftime.datetime objects. Not None.

    Returns
    -------
    A generator object of cftime.datetime objects
    """
    offset = to_offset(freq)

    if start:
        # From pandas GH 56147 / 56832 to account for negative direction and
        # range bounds
        if offset.n >= 0:
            start = offset.rollforward(start)
        else:
            start = offset.rollback(start)

    if periods is None and end < start and offset.n >= 0:
        end = None
        periods = 0

    if end is None:
        end = start + (periods - 1) * offset

    if start is None:
        start = end - (periods - 1) * offset

    current = start
    if offset.n >= 0:
        while current <= end:
            yield current

            next_date = current + offset
            if next_date <= current:
                raise ValueError(f"Offset {offset} did not increment date")
            current = next_date
    else:
        while current >= end:
            yield current

            next_date = current + offset
            if next_date >= current:
                raise ValueError(f"Offset {offset} did not decrement date")
            current = next_date


def cftime_range(
    start=None,
    end=None,
    periods=None,
    freq=None,
    normalize=False,
    name=None,
    inclusive: InclusiveOptions = "both",
    calendar="standard",
) -> CFTimeIndex:
    """Return a fixed frequency CFTimeIndex.

    .. deprecated:: 2025.02.0
        Use :py:func:`~xarray.date_range` with ``use_cftime=True`` instead.

    Parameters
    ----------
    start : str or cftime.datetime, optional
        Left bound for generating dates.
    end : str or cftime.datetime, optional
        Right bound for generating dates.
    periods : int, optional
        Number of periods to generate.
    freq : str or None, default: "D"
        Frequency strings can have multiples, e.g. "5h" and negative values, e.g. "-1D".
    normalize : bool, default: False
        Normalize start/end dates to midnight before generating date range.
    name : str, default: None
        Name of the resulting index
    inclusive : {"both", "neither", "left", "right"}, default "both"
        Include boundaries; whether to set each bound as closed or open.

        .. versionadded:: 2023.02.0
    calendar : str, default: "standard"
        Calendar type for the datetimes.

    Returns
    -------
    CFTimeIndex

    Notes
    -----
    This function is an analog of ``pandas.date_range`` for use in generating
    sequences of ``cftime.datetime`` objects.  It supports most of the
    features of ``pandas.date_range`` (e.g. specifying how the index is
    ``closed`` on either side, or whether or not to ``normalize`` the start and
    end bounds); however, there are some notable exceptions:

    - You cannot specify a ``tz`` (time zone) argument.
    - Start or end dates specified as partial-datetime strings must use the
      `ISO-8601 format <https://en.wikipedia.org/wiki/ISO_8601>`_.
    - It supports many, but not all, frequencies supported by
      ``pandas.date_range``.  For example it does not currently support any of
      the business-related or semi-monthly frequencies.
    - Compound sub-monthly frequencies are not supported, e.g. '1H1min', as
      these can easily be written in terms of the finest common resolution,
      e.g. '61min'.

    Valid simple frequency strings for use with ``cftime``-calendars include
    any multiples of the following.

    +--------+--------------------------+
    | Alias  | Description              |
    +========+==========================+
    | YE     | Year-end frequency       |
    +--------+--------------------------+
    | YS     | Year-start frequency     |
    +--------+--------------------------+
    | QE     | Quarter-end frequency    |
    +--------+--------------------------+
    | QS     | Quarter-start frequency  |
    +--------+--------------------------+
    | ME     | Month-end frequency      |
    +--------+--------------------------+
    | MS     | Month-start frequency    |
    +--------+--------------------------+
    | D      | Day frequency            |
    +--------+--------------------------+
    | h      | Hour frequency           |
    +--------+--------------------------+
    | min    | Minute frequency         |
    +--------+--------------------------+
    | s      | Second frequency         |
    +--------+--------------------------+
    | ms     | Millisecond frequency    |
    +--------+--------------------------+
    | us     | Microsecond frequency    |
    +--------+--------------------------+

    Any multiples of the following anchored offsets are also supported.

    +------------+--------------------------------------------------------------------+
    | Alias      | Description                                                        |
    +============+====================================================================+
    | Y(E,S)-JAN | Annual frequency, anchored at the (end, beginning) of January      |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-FEB | Annual frequency, anchored at the (end, beginning) of February     |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-MAR | Annual frequency, anchored at the (end, beginning) of March        |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-APR | Annual frequency, anchored at the (end, beginning) of April        |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-MAY | Annual frequency, anchored at the (end, beginning) of May          |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-JUN | Annual frequency, anchored at the (end, beginning) of June         |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-JUL | Annual frequency, anchored at the (end, beginning) of July         |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-AUG | Annual frequency, anchored at the (end, beginning) of August       |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-SEP | Annual frequency, anchored at the (end, beginning) of September    |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-OCT | Annual frequency, anchored at the (end, beginning) of October      |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-NOV | Annual frequency, anchored at the (end, beginning) of November     |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-DEC | Annual frequency, anchored at the (end, beginning) of December     |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-JAN | Quarter frequency, anchored at the (end, beginning) of January     |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-FEB | Quarter frequency, anchored at the (end, beginning) of February    |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-MAR | Quarter frequency, anchored at the (end, beginning) of March       |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-APR | Quarter frequency, anchored at the (end, beginning) of April       |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-MAY | Quarter frequency, anchored at the (end, beginning) of May         |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-JUN | Quarter frequency, anchored at the (end, beginning) of June        |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-JUL | Quarter frequency, anchored at the (end, beginning) of July        |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-AUG | Quarter frequency, anchored at the (end, beginning) of August      |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-SEP | Quarter frequency, anchored at the (end, beginning) of September   |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-OCT | Quarter frequency, anchored at the (end, beginning) of October     |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-NOV | Quarter frequency, anchored at the (end, beginning) of November    |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-DEC | Quarter frequency, anchored at the (end, beginning) of December    |
    +------------+--------------------------------------------------------------------+

    Finally, the following calendar aliases are supported.

    +--------------------------------+---------------------------------------+
    | Alias                          | Date type                             |
    +================================+=======================================+
    | standard, gregorian            | ``cftime.DatetimeGregorian``          |
    +--------------------------------+---------------------------------------+
    | proleptic_gregorian            | ``cftime.DatetimeProlepticGregorian`` |
    +--------------------------------+---------------------------------------+
    | noleap, 365_day                | ``cftime.DatetimeNoLeap``             |
    +--------------------------------+---------------------------------------+
    | all_leap, 366_day              | ``cftime.DatetimeAllLeap``            |
    +--------------------------------+---------------------------------------+
    | 360_day                        | ``cftime.Datetime360Day``             |
    +--------------------------------+---------------------------------------+
    | julian                         | ``cftime.DatetimeJulian``             |
    +--------------------------------+---------------------------------------+

    Examples
    --------
    This function returns a ``CFTimeIndex``, populated with ``cftime.datetime``
    objects associated with the specified calendar type, e.g.

    >>> xr.date_range(
    ...     start="2000", periods=6, freq="2MS", calendar="noleap", use_cftime=True
    ... )
    CFTimeIndex([2000-01-01 00:00:00, 2000-03-01 00:00:00, 2000-05-01 00:00:00,
                 2000-07-01 00:00:00, 2000-09-01 00:00:00, 2000-11-01 00:00:00],
                dtype='object', length=6, calendar='noleap', freq='2MS')

    As in the standard pandas function, three of the ``start``, ``end``,
    ``periods``, or ``freq`` arguments must be specified at a given time, with
    the other set to ``None``.  See the `pandas documentation
    <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html>`_
    for more examples of the behavior of ``date_range`` with each of the
    parameters.

    See Also
    --------
    pandas.date_range
    """
    emit_user_level_warning(
        "cftime_range() is deprecated, please use xarray.date_range(..., use_cftime=True) instead.",
        DeprecationWarning,
    )

    return date_range(
        start=start,
        end=end,
        periods=periods,
        freq=freq,
        normalize=normalize,
        name=name,
        inclusive=inclusive,
        calendar=calendar,
        use_cftime=True,
    )


def _cftime_range(
    start=None,
    end=None,
    periods=None,
    freq=None,
    normalize=False,
    name=None,
    inclusive: InclusiveOptions = "both",
    calendar="standard",
) -> CFTimeIndex:
    """Return a fixed frequency CFTimeIndex.

    Parameters
    ----------
    start : str or cftime.datetime, optional
        Left bound for generating dates.
    end : str or cftime.datetime, optional
        Right bound for generating dates.
    periods : int, optional
        Number of periods to generate.
    freq : str or None, default: "D"
        Frequency strings can have multiples, e.g. "5h" and negative values, e.g. "-1D".
    normalize : bool, default: False
        Normalize start/end dates to midnight before generating date range.
    name : str, default: None
        Name of the resulting index
    inclusive : {"both", "neither", "left", "right"}, default "both"
        Include boundaries; whether to set each bound as closed or open.
    calendar : str, default: "standard"
        Calendar type for the datetimes.

    Returns
    -------
    CFTimeIndex

    Notes
    -----
    see cftime_range
    """
    if freq is None and any(arg is None for arg in [periods, start, end]):
        freq = "D"

    # Adapted from pandas.core.indexes.datetimes._generate_range.
    if count_not_none(start, end, periods, freq) != 3:
        raise ValueError(
            "Exactly three of 'start', 'end', 'periods', or 'freq' must be "
            "specified to generate a date range. Note that 'freq' defaults to "
            "'D' in the event that any of 'start', 'end', or 'periods' are "
            "None."
        )

    start = _get_normalized_cfdate(start, calendar, normalize)
    end = _get_normalized_cfdate(end, calendar, normalize)

    if freq is None:
        dates = _generate_linear_date_range(start, end, periods)
    else:
        dates = np.array(
            list(_generate_linear_date_range_with_freq(start, end, periods, freq))
        )

    if not TYPE_CHECKING and inclusive not in get_args(InclusiveOptions):
        raise ValueError(
            f"Argument `inclusive` must be either 'both', 'neither', "
            f"'left', or 'right'.  Got {inclusive}."
        )

    if len(dates) and inclusive != "both":
        if inclusive != "left" and dates[0] == start:
            dates = dates[1:]
        if inclusive != "right" and dates[-1] == end:
            dates = dates[:-1]

    return CFTimeIndex(dates, name=name)


def date_range(
    start=None,
    end=None,
    periods=None,
    freq=None,
    tz=None,
    normalize=False,
    name=None,
    inclusive: InclusiveOptions = "both",
    unit: PDDatetimeUnitOptions = "ns",
    calendar="standard",
    use_cftime=None,
):
    """Return a fixed frequency datetime index.

    The type (:py:class:`xarray.CFTimeIndex` or :py:class:`pandas.DatetimeIndex`)
    of the returned index depends on the requested calendar and on `use_cftime`.

    Parameters
    ----------
    start : str or datetime-like, optional
        Left bound for generating dates.
    end : str or datetime-like, optional
        Right bound for generating dates.
    periods : int, optional
        Number of periods to generate.
    freq : str or None, default: "D"
        Frequency strings can have multiples, e.g. "5h" and negative values, e.g. "-1D".
    tz : str or tzinfo, optional
        Time zone name for returning localized DatetimeIndex, for example
        'Asia/Hong_Kong'. By default, the resulting DatetimeIndex is
        timezone-naive. Only valid with pandas DatetimeIndex.
    normalize : bool, default: False
        Normalize start/end dates to midnight before generating date range.
    name : str, default: None
        Name of the resulting index
    inclusive : {"both", "neither", "left", "right"}, default: "both"
        Include boundaries; whether to set each bound as closed or open.

        .. versionadded:: 2023.02.0
    unit : {"s", "ms", "us", "ns"}, default "ns"
        Specify the desired resolution of the result.

        .. versionadded:: 2024.12.0
    calendar : str, default: "standard"
        Calendar type for the datetimes.
    use_cftime : boolean, optional
        If True, always return a CFTimeIndex.
        If False, return a pd.DatetimeIndex if possible or raise a ValueError.
        If None (default), return a pd.DatetimeIndex if possible,
        otherwise return a CFTimeIndex. Overridden to False if `tz` is not None.

    Returns
    -------
    CFTimeIndex or pd.DatetimeIndex

    Notes
    -----
    When ``use_cftime=True``, or a calendar other than "standard", "gregorian",
    or "proleptic_gregorian" is provided, this function is an analog of ``pandas.date_range``
    for use in generating sequences of ``cftime.datetime`` objects.  It supports most of the
    features of ``pandas.date_range`` (e.g. specifying how the index is
    ``closed`` on either side, or whether or not to ``normalize`` the start and
    end bounds); however, there are some notable exceptions:

    - You cannot specify a ``tz`` (time zone) argument.
    - Start or end dates specified as partial-datetime strings must use the
      `ISO-8601 format <https://en.wikipedia.org/wiki/ISO_8601>`_.
    - It supports many, but not all, frequencies supported by
      ``pandas.date_range``.  For example it does not currently support any of
      the business-related or semi-monthly frequencies.
    - Compound sub-monthly frequencies are not supported, e.g. '1H1min', as
      these can easily be written in terms of the finest common resolution,
      e.g. '61min'.

    Valid simple frequency strings for use with ``cftime``-calendars include
    any multiples of the following.

    +--------+--------------------------+
    | Alias  | Description              |
    +========+==========================+
    | YE     | Year-end frequency       |
    +--------+--------------------------+
    | YS     | Year-start frequency     |
    +--------+--------------------------+
    | QE     | Quarter-end frequency    |
    +--------+--------------------------+
    | QS     | Quarter-start frequency  |
    +--------+--------------------------+
    | ME     | Month-end frequency      |
    +--------+--------------------------+
    | MS     | Month-start frequency    |
    +--------+--------------------------+
    | D      | Day frequency            |
    +--------+--------------------------+
    | h      | Hour frequency           |
    +--------+--------------------------+
    | min    | Minute frequency         |
    +--------+--------------------------+
    | s      | Second frequency         |
    +--------+--------------------------+
    | ms     | Millisecond frequency    |
    +--------+--------------------------+
    | us     | Microsecond frequency    |
    +--------+--------------------------+

    Any multiples of the following anchored offsets are also supported.

    +------------+--------------------------------------------------------------------+
    | Alias      | Description                                                        |
    +============+====================================================================+
    | Y(E,S)-JAN | Annual frequency, anchored at the (end, beginning) of January      |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-FEB | Annual frequency, anchored at the (end, beginning) of February     |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-MAR | Annual frequency, anchored at the (end, beginning) of March        |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-APR | Annual frequency, anchored at the (end, beginning) of April        |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-MAY | Annual frequency, anchored at the (end, beginning) of May          |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-JUN | Annual frequency, anchored at the (end, beginning) of June         |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-JUL | Annual frequency, anchored at the (end, beginning) of July         |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-AUG | Annual frequency, anchored at the (end, beginning) of August       |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-SEP | Annual frequency, anchored at the (end, beginning) of September    |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-OCT | Annual frequency, anchored at the (end, beginning) of October      |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-NOV | Annual frequency, anchored at the (end, beginning) of November     |
    +------------+--------------------------------------------------------------------+
    | Y(E,S)-DEC | Annual frequency, anchored at the (end, beginning) of December     |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-JAN | Quarter frequency, anchored at the (end, beginning) of January     |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-FEB | Quarter frequency, anchored at the (end, beginning) of February    |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-MAR | Quarter frequency, anchored at the (end, beginning) of March       |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-APR | Quarter frequency, anchored at the (end, beginning) of April       |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-MAY | Quarter frequency, anchored at the (end, beginning) of May         |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-JUN | Quarter frequency, anchored at the (end, beginning) of June        |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-JUL | Quarter frequency, anchored at the (end, beginning) of July        |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-AUG | Quarter frequency, anchored at the (end, beginning) of August      |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-SEP | Quarter frequency, anchored at the (end, beginning) of September   |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-OCT | Quarter frequency, anchored at the (end, beginning) of October     |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-NOV | Quarter frequency, anchored at the (end, beginning) of November    |
    +------------+--------------------------------------------------------------------+
    | Q(E,S)-DEC | Quarter frequency, anchored at the (end, beginning) of December    |
    +------------+--------------------------------------------------------------------+

    Finally, the following calendar aliases are supported.

    +--------------------------------+---------------------------------------+----------------------------+
    | Alias                          | Date type                             | Available use_cftime=False |
    +================================+=======================================+============================+
    | standard, gregorian            | ``cftime.DatetimeGregorian``          | True                       |
    +--------------------------------+---------------------------------------+----------------------------+
    | proleptic_gregorian            | ``cftime.DatetimeProlepticGregorian`` | True                       |
    +--------------------------------+---------------------------------------+----------------------------+
    | noleap, 365_day                | ``cftime.DatetimeNoLeap``             | False                      |
    +--------------------------------+---------------------------------------+----------------------------+
    | all_leap, 366_day              | ``cftime.DatetimeAllLeap``            | False                      |
    +--------------------------------+---------------------------------------+----------------------------+
    | 360_day                        | ``cftime.Datetime360Day``             | False                      |
    +--------------------------------+---------------------------------------+----------------------------+
    | julian                         | ``cftime.DatetimeJulian``             | False                      |
    +--------------------------------+---------------------------------------+----------------------------+

    As in the standard pandas function, exactly three of ``start``, ``end``,
    ``periods``, or ``freq`` are required to generate a date range. Note that
    ``freq`` defaults to ``"D"`` in the event that any of ``start``, ``end``,
    or ``periods`` are set to ``None``. See :py:func:`pandas.date_range`.
    for more examples of the behavior of ``date_range`` with each of the
    parameters.

    Examples
    --------
    This function returns a ``CFTimeIndex``, populated with ``cftime.datetime``
    objects associated with the specified calendar type, e.g.

    >>> xr.date_range(
    ...     start="2000", periods=6, freq="2MS", calendar="noleap", use_cftime=True
    ... )
    CFTimeIndex([2000-01-01 00:00:00, 2000-03-01 00:00:00, 2000-05-01 00:00:00,
                 2000-07-01 00:00:00, 2000-09-01 00:00:00, 2000-11-01 00:00:00],
                dtype='object', length=6, calendar='noleap', freq='2MS')

    See also
    --------
    pandas.date_range
    cftime_range
    date_range_like
    """
    if tz is not None:
        use_cftime = False

    if _is_standard_calendar(calendar) and use_cftime is not True:
        try:
            return pd.date_range(
                start=start,
                end=end,
                periods=periods,
                # TODO remove translation once requiring pandas >= 2.2
                freq=_new_to_legacy_freq(freq),
                tz=tz,
                normalize=normalize,
                name=name,
                inclusive=inclusive,
                unit=unit,
            )
        except pd.errors.OutOfBoundsDatetime as err:
            if use_cftime is False:
                raise ValueError(
                    "Date range is invalid for pandas DatetimeIndex, try using `use_cftime=True`."
                ) from err
    elif use_cftime is False:
        raise ValueError(
            f"Invalid calendar {calendar} for pandas DatetimeIndex, try using `use_cftime=True`."
        )

    return _cftime_range(
        start=start,
        end=end,
        periods=periods,
        freq=freq,
        normalize=normalize,
        name=name,
        inclusive=inclusive,
        calendar=calendar,
    )


def _new_to_legacy_freq(freq):
    # xarray will now always return "ME" and "QE" for MonthEnd and QuarterEnd
    # frequencies, but older versions of pandas do not support these as
    # frequency strings.  Until xarray's minimum pandas version is 2.2 or above,
    # we add logic to continue using the deprecated "M" and "Q" frequency
    # strings in these circumstances.

    # NOTE: other conversions ("h" -> "H", ..., "ns" -> "N") not required

    # TODO: remove once requiring pandas >= 2.2
    if not freq or Version(pd.__version__) >= Version("2.2"):
        return freq

    try:
        freq_as_offset = to_offset(freq)
    except ValueError:
        # freq may be valid in pandas but not in xarray
        return freq

    if isinstance(freq_as_offset, MonthEnd) and "ME" in freq:
        freq = freq.replace("ME", "M")
    elif isinstance(freq_as_offset, QuarterEnd) and "QE" in freq:
        freq = freq.replace("QE", "Q")
    elif isinstance(freq_as_offset, YearBegin) and "YS" in freq:
        freq = freq.replace("YS", "AS")
    elif isinstance(freq_as_offset, YearEnd):
        # testing for "Y" is required as this was valid in xarray 2023.11 - 2024.01
        if "Y-" in freq:
            # Check for and replace "Y-" instead of just "Y" to prevent
            # corrupting anchored offsets that contain "Y" in the month
            # abbreviation, e.g. "Y-MAY" -> "A-MAY".
            freq = freq.replace("Y-", "A-")
        elif "YE-" in freq:
            freq = freq.replace("YE-", "A-")
        elif "A-" not in freq and freq.endswith("Y"):
            freq = freq.replace("Y", "A")
        elif freq.endswith("YE"):
            freq = freq.replace("YE", "A")

    return freq


def _legacy_to_new_freq(freq: T_FreqStr) -> T_FreqStr:
    # to avoid internal deprecation warnings when freq is determined using pandas < 2.2

    # TODO: remove once requiring pandas >= 2.2

    if not freq or Version(pd.__version__) >= Version("2.2"):
        return freq

    try:
        freq_as_offset = to_offset(freq, warn=False)
    except ValueError:
        # freq may be valid in pandas but not in xarray
        return freq

    if isinstance(freq_as_offset, MonthEnd) and "ME" not in freq:
        freq = freq.replace("M", "ME")
    elif isinstance(freq_as_offset, QuarterEnd) and "QE" not in freq:
        freq = freq.replace("Q", "QE")
    elif isinstance(freq_as_offset, YearBegin) and "YS" not in freq:
        freq = freq.replace("AS", "YS")
    elif isinstance(freq_as_offset, YearEnd):
        if "A-" in freq:
            # Check for and replace "A-" instead of just "A" to prevent
            # corrupting anchored offsets that contain "Y" in the month
            # abbreviation, e.g. "A-MAY" -> "YE-MAY".
            freq = freq.replace("A-", "YE-")
        elif "Y-" in freq:
            freq = freq.replace("Y-", "YE-")
        elif freq.endswith("A"):
            # the "A-MAY" case is already handled above
            freq = freq.replace("A", "YE")
        elif "YE" not in freq and freq.endswith("Y"):
            # the "Y-MAY" case is already handled above
            freq = freq.replace("Y", "YE")
    elif isinstance(freq_as_offset, Hour):
        freq = freq.replace("H", "h")
    elif isinstance(freq_as_offset, Minute):
        freq = freq.replace("T", "min")
    elif isinstance(freq_as_offset, Second):
        freq = freq.replace("S", "s")
    elif isinstance(freq_as_offset, Millisecond):
        freq = freq.replace("L", "ms")
    elif isinstance(freq_as_offset, Microsecond):
        freq = freq.replace("U", "us")

    return freq


def date_range_like(source, calendar, use_cftime=None):
    """Generate a datetime array with the same frequency, start and end as
    another one, but in a different calendar.

    Parameters
    ----------
    source : DataArray, CFTimeIndex, or pd.DatetimeIndex
        1D datetime array
    calendar : str
        New calendar name.
    use_cftime : bool, optional
        If True, the output uses :py:class:`cftime.datetime` objects.
        If None (default), :py:class:`numpy.datetime64` values are used if possible.
        If False, :py:class:`numpy.datetime64` values are used or an error is raised.

    Returns
    -------
    DataArray
        1D datetime coordinate with the same start, end and frequency as the
        source, but in the new calendar. The start date is assumed to exist in
        the target calendar. If the end date doesn't exist, the code tries 1
        and 2 calendar days before. There is a special case when the source time
        series is daily or coarser and the end of the input range is on the
        last day of the month. Then the output range will also end on the last
        day of the month in the new calendar.
    """
    from xarray.coding.frequencies import infer_freq
    from xarray.core.dataarray import DataArray

    if not isinstance(source, pd.DatetimeIndex | CFTimeIndex) and (
        (isinstance(source, DataArray) and (source.ndim != 1))
        or not _contains_datetime_like_objects(source.variable)
    ):
        raise ValueError(
            "'source' must be a 1D array of datetime objects for inferring its range."
        )

    freq = infer_freq(source)
    if freq is None:
        raise ValueError(
            "`date_range_like` was unable to generate a range as the source frequency was not inferable."
        )

    # TODO remove once requiring pandas >= 2.2
    freq = _legacy_to_new_freq(freq)

    use_cftime = _should_cftime_be_used(source, calendar, use_cftime)

    source_start = source.values.min()
    source_end = source.values.max()

    freq_as_offset = to_offset(freq)
    if freq_as_offset.n < 0:
        source_start, source_end = source_end, source_start

    if is_np_datetime_like(source.dtype):
        # We want to use datetime fields (datetime64 object don't have them)
        source_calendar = "standard"
        source_start = default_precision_timestamp(source_start)
        source_end = default_precision_timestamp(source_end)
    elif isinstance(source, CFTimeIndex):
        source_calendar = source.calendar
    else:  # DataArray
        source_calendar = source.dt.calendar

    if calendar == source_calendar and is_np_datetime_like(source.dtype) ^ use_cftime:
        return source

    date_type = get_date_type(calendar, use_cftime)
    start = convert_time_or_go_back(source_start, date_type)
    end = convert_time_or_go_back(source_end, date_type)

    # For the cases where the source ends on the end of the month, we expect the same in the new calendar.
    if source_end.day == source_end.daysinmonth and isinstance(
        freq_as_offset, YearEnd | QuarterEnd | MonthEnd | Day
    ):
        end = end.replace(day=end.daysinmonth)

    return date_range(
        start=start.isoformat(),
        end=end.isoformat(),
        freq=freq,
        calendar=calendar,
    )
