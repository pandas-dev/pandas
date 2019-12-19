from datetime import timedelta
import re
from typing import Dict, Optional
import warnings

import numpy as np
from pytz import AmbiguousTimeError

from pandas._libs.algos import unique_deltas
from pandas._libs.tslibs import Timedelta, Timestamp
from pandas._libs.tslibs.ccalendar import MONTH_ALIASES, int_to_weekday
from pandas._libs.tslibs.fields import build_field_sarray
import pandas._libs.tslibs.frequencies as libfreqs
from pandas._libs.tslibs.offsets import _offset_to_period_map
import pandas._libs.tslibs.resolution as libresolution
from pandas._libs.tslibs.resolution import Resolution
from pandas._libs.tslibs.timezones import UTC
from pandas._libs.tslibs.tzconversion import tz_convert
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_period_arraylike,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.algorithms import unique

from pandas.tseries.offsets import (
    DateOffset,
    Day,
    Hour,
    Micro,
    Milli,
    Minute,
    Nano,
    Second,
    prefix_mapping,
)

_ONE_MICRO = 1000
_ONE_MILLI = _ONE_MICRO * 1000
_ONE_SECOND = _ONE_MILLI * 1000
_ONE_MINUTE = 60 * _ONE_SECOND
_ONE_HOUR = 60 * _ONE_MINUTE
_ONE_DAY = 24 * _ONE_HOUR

# ---------------------------------------------------------------------
# Offset names ("time rules") and related functions

#: cache of previously seen offsets
_offset_map: Dict[str, DateOffset] = {}


def get_period_alias(offset_str: str) -> Optional[str]:
    """
    Alias to closest period strings BQ->Q etc.
    """
    return _offset_to_period_map.get(offset_str, None)


_name_to_offset_map = {
    "days": Day(1),
    "hours": Hour(1),
    "minutes": Minute(1),
    "seconds": Second(1),
    "milliseconds": Milli(1),
    "microseconds": Micro(1),
    "nanoseconds": Nano(1),
}


def to_offset(freq) -> Optional[DateOffset]:
    """
    Return DateOffset object from string or tuple representation
    or datetime.timedelta object.

    Parameters
    ----------
    freq : str, tuple, datetime.timedelta, DateOffset or None

    Returns
    -------
    DateOffset
        None if freq is None.

    Raises
    ------
    ValueError
        If freq is an invalid frequency

    See Also
    --------
    DateOffset

    Examples
    --------
    >>> to_offset('5min')
    <5 * Minutes>

    >>> to_offset('1D1H')
    <25 * Hours>

    >>> to_offset(('W', 2))
    <2 * Weeks: weekday=6>

    >>> to_offset((2, 'B'))
    <2 * BusinessDays>

    >>> to_offset(datetime.timedelta(days=1))
    <Day>

    >>> to_offset(Hour())
    <Hour>
    """
    if freq is None:
        return None

    if isinstance(freq, DateOffset):
        return freq

    if isinstance(freq, tuple):
        name = freq[0]
        stride = freq[1]
        if isinstance(stride, str):
            name, stride = stride, name
        name, _ = libfreqs._base_and_stride(name)
        delta = _get_offset(name) * stride

    elif isinstance(freq, timedelta):
        delta = None
        freq = Timedelta(freq)
        try:
            for name in freq.components._fields:
                offset = _name_to_offset_map[name]
                stride = getattr(freq.components, name)
                if stride != 0:
                    offset = stride * offset
                    if delta is None:
                        delta = offset
                    else:
                        delta = delta + offset
        except ValueError:
            raise ValueError(libfreqs.INVALID_FREQ_ERR_MSG.format(freq))

    else:
        delta = None
        stride_sign = None
        try:
            splitted = re.split(libfreqs.opattern, freq)
            if splitted[-1] != "" and not splitted[-1].isspace():
                # the last element must be blank
                raise ValueError("last element must be blank")
            for sep, stride, name in zip(
                splitted[0::4], splitted[1::4], splitted[2::4]
            ):
                if sep != "" and not sep.isspace():
                    raise ValueError("separator must be spaces")
                prefix = libfreqs._lite_rule_alias.get(name) or name
                if stride_sign is None:
                    stride_sign = -1 if stride.startswith("-") else 1
                if not stride:
                    stride = 1
                if prefix in Resolution._reso_str_bump_map.keys():
                    stride, name = Resolution.get_stride_from_decimal(
                        float(stride), prefix
                    )
                stride = int(stride)
                offset = _get_offset(name)
                offset = offset * int(np.fabs(stride) * stride_sign)
                if delta is None:
                    delta = offset
                else:
                    delta = delta + offset
        except (ValueError, TypeError):
            raise ValueError(libfreqs.INVALID_FREQ_ERR_MSG.format(freq))

    if delta is None:
        raise ValueError(libfreqs.INVALID_FREQ_ERR_MSG.format(freq))

    return delta


def get_offset(name: str) -> DateOffset:
    """
    Return DateOffset object associated with rule name.

    .. deprecated:: 1.0.0

    Examples
    --------
    get_offset('EOM') --> BMonthEnd(1)
    """
    warnings.warn(
        "get_offset is deprecated and will be removed in a future version, "
        "use to_offset instead",
        FutureWarning,
        stacklevel=2,
    )
    return _get_offset(name)


def _get_offset(name: str) -> DateOffset:
    """
    Return DateOffset object associated with rule name.

    Examples
    --------
    _get_offset('EOM') --> BMonthEnd(1)
    """
    if name not in libfreqs._dont_uppercase:
        name = name.upper()
        name = libfreqs._lite_rule_alias.get(name, name)
        name = libfreqs._lite_rule_alias.get(name.lower(), name)
    else:
        name = libfreqs._lite_rule_alias.get(name, name)

    if name not in _offset_map:
        try:
            split = name.split("-")
            klass = prefix_mapping[split[0]]
            # handles case where there's no suffix (and will TypeError if too
            # many '-')
            offset = klass._from_name(*split[1:])
        except (ValueError, TypeError, KeyError):
            # bad prefix or suffix
            raise ValueError(libfreqs.INVALID_FREQ_ERR_MSG.format(name))
        # cache
        _offset_map[name] = offset

    return _offset_map[name]


# ---------------------------------------------------------------------
# Period codes


def infer_freq(index, warn: bool = True) -> Optional[str]:
    """
    Infer the most likely frequency given the input index. If the frequency is
    uncertain, a warning will be printed.

    Parameters
    ----------
    index : DatetimeIndex or TimedeltaIndex
      if passed a Series will use the values of the series (NOT THE INDEX).
    warn : bool, default True

    Returns
    -------
    str or None
        None if no discernible frequency
        TypeError if the index is not datetime-like
        ValueError if there are less than three values.
    """
    import pandas as pd

    if isinstance(index, ABCSeries):
        values = index._values
        if not (
            is_datetime64_dtype(values)
            or is_timedelta64_dtype(values)
            or values.dtype == object
        ):
            raise TypeError(
                "cannot infer freq from a non-convertible dtype "
                f"on a Series of {index.dtype}"
            )
        index = values

    inferer: _FrequencyInferer
    if is_period_arraylike(index):
        raise TypeError(
            "PeriodIndex given. Check the `freq` attribute "
            "instead of using infer_freq."
        )
    elif is_timedelta64_dtype(index):
        # Allow TimedeltaIndex and TimedeltaArray
        inferer = _TimedeltaFrequencyInferer(index, warn=warn)
        return inferer.get_freq()

    if isinstance(index, pd.Index) and not isinstance(index, pd.DatetimeIndex):
        if isinstance(index, (pd.Int64Index, pd.Float64Index)):
            raise TypeError(
                f"cannot infer freq from a non-convertible index type {type(index)}"
            )
        index = index.values

    if not isinstance(index, pd.DatetimeIndex):
        try:
            index = pd.DatetimeIndex(index)
        except AmbiguousTimeError:
            index = pd.DatetimeIndex(index.asi8)

    inferer = _FrequencyInferer(index, warn=warn)
    return inferer.get_freq()


class _FrequencyInferer:
    """
    Not sure if I can avoid the state machine here
    """

    def __init__(self, index, warn: bool = True):
        self.index = index
        self.values = index.asi8

        # This moves the values, which are implicitly in UTC, to the
        # the timezone so they are in local time
        if hasattr(index, "tz"):
            if index.tz is not None:
                self.values = tz_convert(self.values, UTC, index.tz)

        self.warn = warn

        if len(index) < 3:
            raise ValueError("Need at least 3 dates to infer frequency")

        self.is_monotonic = (
            self.index._is_monotonic_increasing or self.index._is_monotonic_decreasing
        )

    @cache_readonly
    def deltas(self):
        return unique_deltas(self.values)

    @cache_readonly
    def deltas_asi8(self):
        return unique_deltas(self.index.asi8)

    @cache_readonly
    def is_unique(self) -> bool:
        return len(self.deltas) == 1

    @cache_readonly
    def is_unique_asi8(self):
        return len(self.deltas_asi8) == 1

    def get_freq(self) -> Optional[str]:
        """
        Find the appropriate frequency string to describe the inferred
        frequency of self.values

        Returns
        -------
        str or None
        """
        if not self.is_monotonic or not self.index._is_unique:
            return None

        delta = self.deltas[0]
        if _is_multiple(delta, _ONE_DAY):
            return self._infer_daily_rule()

        # Business hourly, maybe. 17: one day / 65: one weekend
        if self.hour_deltas in ([1, 17], [1, 65], [1, 17, 65]):
            return "BH"
        # Possibly intraday frequency.  Here we use the
        # original .asi8 values as the modified values
        # will not work around DST transitions.  See #8772
        elif not self.is_unique_asi8:
            return None

        delta = self.deltas_asi8[0]
        if _is_multiple(delta, _ONE_HOUR):
            # Hours
            return _maybe_add_count("H", delta / _ONE_HOUR)
        elif _is_multiple(delta, _ONE_MINUTE):
            # Minutes
            return _maybe_add_count("T", delta / _ONE_MINUTE)
        elif _is_multiple(delta, _ONE_SECOND):
            # Seconds
            return _maybe_add_count("S", delta / _ONE_SECOND)
        elif _is_multiple(delta, _ONE_MILLI):
            # Milliseconds
            return _maybe_add_count("L", delta / _ONE_MILLI)
        elif _is_multiple(delta, _ONE_MICRO):
            # Microseconds
            return _maybe_add_count("U", delta / _ONE_MICRO)
        else:
            # Nanoseconds
            return _maybe_add_count("N", delta)

    @cache_readonly
    def day_deltas(self):
        return [x / _ONE_DAY for x in self.deltas]

    @cache_readonly
    def hour_deltas(self):
        return [x / _ONE_HOUR for x in self.deltas]

    @cache_readonly
    def fields(self):
        return build_field_sarray(self.values)

    @cache_readonly
    def rep_stamp(self):
        return Timestamp(self.values[0])

    def month_position_check(self):
        return libresolution.month_position_check(self.fields, self.index.dayofweek)

    @cache_readonly
    def mdiffs(self):
        nmonths = self.fields["Y"] * 12 + self.fields["M"]
        return unique_deltas(nmonths.astype("i8"))

    @cache_readonly
    def ydiffs(self):
        return unique_deltas(self.fields["Y"].astype("i8"))

    def _infer_daily_rule(self) -> Optional[str]:
        annual_rule = self._get_annual_rule()
        if annual_rule:
            nyears = self.ydiffs[0]
            month = MONTH_ALIASES[self.rep_stamp.month]
            alias = f"{annual_rule}-{month}"
            return _maybe_add_count(alias, nyears)

        quarterly_rule = self._get_quarterly_rule()
        if quarterly_rule:
            nquarters = self.mdiffs[0] / 3
            mod_dict = {0: 12, 2: 11, 1: 10}
            month = MONTH_ALIASES[mod_dict[self.rep_stamp.month % 3]]
            alias = f"{quarterly_rule}-{month}"
            return _maybe_add_count(alias, nquarters)

        monthly_rule = self._get_monthly_rule()
        if monthly_rule:
            return _maybe_add_count(monthly_rule, self.mdiffs[0])

        if self.is_unique:
            days = self.deltas[0] / _ONE_DAY
            if days % 7 == 0:
                # Weekly
                day = int_to_weekday[self.rep_stamp.weekday()]
                return _maybe_add_count(f"W-{day}", days / 7)
            else:
                return _maybe_add_count("D", days)

        if self._is_business_daily():
            return "B"

        wom_rule = self._get_wom_rule()
        if wom_rule:
            return wom_rule

        return None

    def _get_annual_rule(self) -> Optional[str]:
        if len(self.ydiffs) > 1:
            return None

        if len(unique(self.fields["M"])) > 1:
            return None

        pos_check = self.month_position_check()
        return {"cs": "AS", "bs": "BAS", "ce": "A", "be": "BA"}.get(pos_check)

    def _get_quarterly_rule(self) -> Optional[str]:
        if len(self.mdiffs) > 1:
            return None

        if not self.mdiffs[0] % 3 == 0:
            return None

        pos_check = self.month_position_check()
        return {"cs": "QS", "bs": "BQS", "ce": "Q", "be": "BQ"}.get(pos_check)

    def _get_monthly_rule(self) -> Optional[str]:
        if len(self.mdiffs) > 1:
            return None
        pos_check = self.month_position_check()
        return {"cs": "MS", "bs": "BMS", "ce": "M", "be": "BM"}.get(pos_check)

    def _is_business_daily(self) -> bool:
        # quick check: cannot be business daily
        if self.day_deltas != [1, 3]:
            return False

        # probably business daily, but need to confirm
        first_weekday = self.index[0].weekday()
        shifts = np.diff(self.index.asi8)
        shifts = np.floor_divide(shifts, _ONE_DAY)
        weekdays = np.mod(first_weekday + np.cumsum(shifts), 7)
        return np.all(
            ((weekdays == 0) & (shifts == 3))
            | ((weekdays > 0) & (weekdays <= 4) & (shifts == 1))
        )

    def _get_wom_rule(self) -> Optional[str]:
        #         wdiffs = unique(np.diff(self.index.week))
        # We also need -47, -49, -48 to catch index spanning year boundary
        #     if not lib.ismember(wdiffs, set([4, 5, -47, -49, -48])).all():
        #         return None

        weekdays = unique(self.index.weekday)
        if len(weekdays) > 1:
            return None

        week_of_months = unique((self.index.day - 1) // 7)
        # Only attempt to infer up to WOM-4. See #9425
        week_of_months = week_of_months[week_of_months < 4]
        if len(week_of_months) == 0 or len(week_of_months) > 1:
            return None

        # get which week
        week = week_of_months[0] + 1
        wd = int_to_weekday[weekdays[0]]

        return f"WOM-{week}{wd}"


class _TimedeltaFrequencyInferer(_FrequencyInferer):
    def _infer_daily_rule(self):
        if self.is_unique:
            days = self.deltas[0] / _ONE_DAY
            if days % 7 == 0:
                # Weekly
                wd = int_to_weekday[self.rep_stamp.weekday()]
                alias = f"W-{wd}"
                return _maybe_add_count(alias, days / 7)
            else:
                return _maybe_add_count("D", days)


def _is_multiple(us, mult: int) -> bool:
    return us % mult == 0


def _maybe_add_count(base: str, count: float) -> str:
    if count != 1:
        assert count == int(count)
        count = int(count)
        return f"{count}{base}"
    else:
        return base
