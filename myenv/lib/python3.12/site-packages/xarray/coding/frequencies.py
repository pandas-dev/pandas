"""FrequencyInferer analog for cftime.datetime objects"""

# The infer_freq method and the _CFTimeFrequencyInferer
# subclass defined here were copied and adapted for
# use with cftime.datetime objects based on the source code in
# pandas.tseries.Frequencies._FrequencyInferer

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

import numpy as np
import pandas as pd

from xarray.coding.cftime_offsets import _MONTH_ABBREVIATIONS, _legacy_to_new_freq
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.core.common import _contains_datetime_like_objects

_ONE_MICRO = 1
_ONE_MILLI = _ONE_MICRO * 1000
_ONE_SECOND = _ONE_MILLI * 1000
_ONE_MINUTE = 60 * _ONE_SECOND
_ONE_HOUR = 60 * _ONE_MINUTE
_ONE_DAY = 24 * _ONE_HOUR


def infer_freq(index):
    """
    Infer the most likely frequency given the input index.

    Parameters
    ----------
    index : CFTimeIndex, DataArray, DatetimeIndex, TimedeltaIndex, Series
        If not passed a CFTimeIndex, this simply calls `pandas.infer_freq`.
        If passed a Series or a DataArray will use the values of the series (NOT THE INDEX).

    Returns
    -------
    str or None
        None if no discernible frequency.

    Raises
    ------
    TypeError
        If the index is not datetime-like.
    ValueError
        If there are fewer than three values or the index is not 1D.
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.variable import Variable

    if isinstance(index, (DataArray, pd.Series)):
        if index.ndim != 1:
            raise ValueError("'index' must be 1D")
        elif not _contains_datetime_like_objects(Variable("dim", index)):
            raise ValueError("'index' must contain datetime-like objects")
        dtype = np.asarray(index).dtype
        if dtype == "datetime64[ns]":
            index = pd.DatetimeIndex(index.values)
        elif dtype == "timedelta64[ns]":
            index = pd.TimedeltaIndex(index.values)
        else:
            index = CFTimeIndex(index.values)

    if isinstance(index, CFTimeIndex):
        inferer = _CFTimeFrequencyInferer(index)
        return inferer.get_freq()

    return _legacy_to_new_freq(pd.infer_freq(index))


class _CFTimeFrequencyInferer:  # (pd.tseries.frequencies._FrequencyInferer):
    def __init__(self, index):
        self.index = index
        self.values = index.asi8

        if len(index) < 3:
            raise ValueError("Need at least 3 dates to infer frequency")

        self.is_monotonic = (
            self.index.is_monotonic_decreasing or self.index.is_monotonic_increasing
        )

        self._deltas = None
        self._year_deltas = None
        self._month_deltas = None

    def get_freq(self):
        """Find the appropriate frequency string to describe the inferred frequency of self.index

        Adapted from `pandas.tsseries.frequencies._FrequencyInferer.get_freq` for CFTimeIndexes.

        Returns
        -------
        str or None
        """
        if not self.is_monotonic or not self.index.is_unique:
            return None

        delta = self.deltas[0]  # Smallest delta
        if _is_multiple(delta, _ONE_DAY):
            return self._infer_daily_rule()
        # There is no possible intraday frequency with a non-unique delta
        # Different from pandas: we don't need to manage DST and business offsets in cftime
        elif not len(self.deltas) == 1:
            return None

        if _is_multiple(delta, _ONE_HOUR):
            return _maybe_add_count("h", delta / _ONE_HOUR)
        elif _is_multiple(delta, _ONE_MINUTE):
            return _maybe_add_count("min", delta / _ONE_MINUTE)
        elif _is_multiple(delta, _ONE_SECOND):
            return _maybe_add_count("s", delta / _ONE_SECOND)
        elif _is_multiple(delta, _ONE_MILLI):
            return _maybe_add_count("ms", delta / _ONE_MILLI)
        else:
            return _maybe_add_count("us", delta / _ONE_MICRO)

    def _infer_daily_rule(self):
        annual_rule = self._get_annual_rule()
        if annual_rule:
            nyears = self.year_deltas[0]
            month = _MONTH_ABBREVIATIONS[self.index[0].month]
            alias = f"{annual_rule}-{month}"
            return _maybe_add_count(alias, nyears)

        quartely_rule = self._get_quartely_rule()
        if quartely_rule:
            nquarters = self.month_deltas[0] / 3
            mod_dict = {0: 12, 2: 11, 1: 10}
            month = _MONTH_ABBREVIATIONS[mod_dict[self.index[0].month % 3]]
            alias = f"{quartely_rule}-{month}"
            return _maybe_add_count(alias, nquarters)

        monthly_rule = self._get_monthly_rule()
        if monthly_rule:
            return _maybe_add_count(monthly_rule, self.month_deltas[0])

        if len(self.deltas) == 1:
            # Daily as there is no "Weekly" offsets with CFTime
            days = self.deltas[0] / _ONE_DAY
            return _maybe_add_count("D", days)

        # CFTime has no business freq and no "week of month" (WOM)
        return None

    def _get_annual_rule(self):
        if len(self.year_deltas) > 1:
            return None

        if len(np.unique(self.index.month)) > 1:
            return None

        return {"cs": "YS", "ce": "YE"}.get(month_anchor_check(self.index))

    def _get_quartely_rule(self):
        if len(self.month_deltas) > 1:
            return None

        if self.month_deltas[0] % 3 != 0:
            return None

        return {"cs": "QS", "ce": "QE"}.get(month_anchor_check(self.index))

    def _get_monthly_rule(self):
        if len(self.month_deltas) > 1:
            return None

        return {"cs": "MS", "ce": "ME"}.get(month_anchor_check(self.index))

    @property
    def deltas(self):
        """Sorted unique timedeltas as microseconds."""
        if self._deltas is None:
            self._deltas = _unique_deltas(self.values)
        return self._deltas

    @property
    def year_deltas(self):
        """Sorted unique year deltas."""
        if self._year_deltas is None:
            self._year_deltas = _unique_deltas(self.index.year)
        return self._year_deltas

    @property
    def month_deltas(self):
        """Sorted unique month deltas."""
        if self._month_deltas is None:
            self._month_deltas = _unique_deltas(self.index.year * 12 + self.index.month)
        return self._month_deltas


def _unique_deltas(arr):
    """Sorted unique deltas of numpy array"""
    return np.sort(np.unique(np.diff(arr)))


def _is_multiple(us, mult: int):
    """Whether us is a multiple of mult"""
    return us % mult == 0


def _maybe_add_count(base: str, count: float):
    """If count is greater than 1, add it to the base offset string"""
    if count != 1:
        assert count == int(count)
        count = int(count)
        return f"{count}{base}"
    else:
        return base


def month_anchor_check(dates):
    """Return the monthly offset string.

    Return "cs" if all dates are the first days of the month,
    "ce" if all dates are the last day of the month,
    None otherwise.

    Replicated pandas._libs.tslibs.resolution.month_position_check
    but without business offset handling.
    """
    calendar_end = True
    calendar_start = True

    for date in dates:
        if calendar_start:
            calendar_start &= date.day == 1

        if calendar_end:
            cal = date.day == date.daysinmonth
            calendar_end &= cal
        elif not calendar_start:
            break

    if calendar_end:
        return "ce"
    elif calendar_start:
        return "cs"
    else:
        return None
