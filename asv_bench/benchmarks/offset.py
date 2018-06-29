# -*- coding: utf-8 -*-
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
try:
    import pandas.tseries.holiday  # noqa
except ImportError:
    pass

hcal = pd.tseries.holiday.USFederalHolidayCalendar()
# These offests currently raise a NotImplimentedError with .apply_index()
non_apply = [pd.offsets.Day(),
             pd.offsets.BYearEnd(),
             pd.offsets.BYearBegin(),
             pd.offsets.BQuarterEnd(),
             pd.offsets.BQuarterBegin(),
             pd.offsets.BMonthEnd(),
             pd.offsets.BMonthBegin(),
             pd.offsets.CustomBusinessDay(),
             pd.offsets.CustomBusinessDay(calendar=hcal),
             pd.offsets.CustomBusinessMonthBegin(calendar=hcal),
             pd.offsets.CustomBusinessMonthEnd(calendar=hcal),
             pd.offsets.CustomBusinessMonthEnd(calendar=hcal)]
other_offsets = [pd.offsets.YearEnd(), pd.offsets.YearBegin(),
                 pd.offsets.QuarterEnd(), pd.offsets.QuarterBegin(),
                 pd.offsets.MonthEnd(), pd.offsets.MonthBegin(),
                 pd.offsets.DateOffset(months=2, days=2),
                 pd.offsets.BusinessDay(), pd.offsets.SemiMonthEnd(),
                 pd.offsets.SemiMonthBegin()]
offsets = non_apply + other_offsets


class ApplyIndex(object):

    goal_time = 0.2

    params = other_offsets
    param_names = ['offset']

    def setup(self, offset):
        N = 10000
        self.rng = pd.date_range(start='1/1/2000', periods=N, freq='T')

    def time_apply_index(self, offset):
        offset.apply_index(self.rng)


class OnOffset(object):

    goal_time = 0.2

    params = offsets
    param_names = ['offset']

    def setup(self, offset):
        self.dates = [datetime(2016, m, d)
                      for m in [10, 11, 12]
                      for d in [1, 2, 3, 28, 29, 30, 31]
                      if not (m == 11 and d == 31)]

    def time_on_offset(self, offset):
        for date in self.dates:
            offset.onOffset(date)


class OffsetSeriesArithmetic(object):

    goal_time = 0.2
    params = offsets
    param_names = ['offset']

    def setup(self, offset):
        N = 1000
        rng = pd.date_range(start='1/1/2000', periods=N, freq='T')
        self.data = pd.Series(rng)

    def time_add_offset(self, offset):
        with warnings.catch_warnings(record=True):
            self.data + offset


class OffsetDatetimeIndexArithmetic(object):

    goal_time = 0.2
    params = offsets
    param_names = ['offset']

    def setup(self, offset):
        N = 1000
        self.data = pd.date_range(start='1/1/2000', periods=N, freq='T')

    def time_add_offset(self, offset):
        with warnings.catch_warnings(record=True):
            self.data + offset


class OffestDatetimeArithmetic(object):

    goal_time = 0.2
    params = offsets
    param_names = ['offset']

    def setup(self, offset):
        self.date = datetime(2011, 1, 1)
        self.dt64 = np.datetime64('2011-01-01 09:00Z')

    def time_apply(self, offset):
        offset.apply(self.date)

    def time_apply_np_dt64(self, offset):
        offset.apply(self.dt64)

    def time_add(self, offset):
        self.date + offset

    def time_add_10(self, offset):
        self.date + (10 * offset)

    def time_subtract(self, offset):
        self.date - offset

    def time_subtract_10(self, offset):
        self.date - (10 * offset)
