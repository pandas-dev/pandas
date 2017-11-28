# -*- coding: utf-8 -*-
from datetime import datetime

import numpy as np

import pandas as pd
from pandas import date_range

try:
    import pandas.tseries.holiday
except ImportError:
    pass

hcal = pd.tseries.holiday.USFederalHolidayCalendar()


class ApplyIndex(object):
    goal_time = 0.2

    params = [pd.offsets.YearEnd(), pd.offsets.YearBegin(),
              pd.offsets.BYearEnd(), pd.offsets.BYearBegin(),
              pd.offsets.QuarterEnd(), pd.offsets.QuarterBegin(),
              pd.offsets.BQuarterEnd(), pd.offsets.BQuarterBegin(),
              pd.offsets.MonthEnd(), pd.offsets.MonthBegin(),
              pd.offsets.BMonthEnd(), pd.offsets.BMonthBegin()]

    def setup(self, param):
        self.offset = param

        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')
        self.ser = pd.Series(self.rng)

    def time_apply_index(self, param):
        self.rng + self.offset

    def time_apply_series(self, param):
        self.ser + self.offset


class OnOffset(object):
    goal_time = 0.2

    params = [pd.offsets.QuarterBegin(), pd.offsets.QuarterEnd(),
              pd.offsets.BQuarterBegin(), pd.offsets.BQuarterEnd()]
    param_names = ['offset']

    def setup(self, offset):
        self.offset = offset
        self.dates = [datetime(2016, m, d)
                      for m in [10, 11, 12]
                      for d in [1, 2, 3, 28, 29, 30, 31]
                      if not (m == 11 and d == 31)]

    def time_on_offset(self, offset):
        for date in self.dates:
            self.offset.onOffset(date)


class DatetimeIndexArithmetic(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')
        self.day_offset = pd.offsets.Day()
        self.relativedelta_offset = pd.offsets.DateOffset(months=2, days=2)
        self.busday_offset = pd.offsets.BusinessDay()

    def time_add_offset_delta(self):
        self.rng + self.day_offset

    def time_add_offset_fast(self):
        self.rng + self.relativedelta_offset

    def time_add_offset_slow(self):
        self.rng + self.busday_offset


class SeriesArithmetic(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        rng = date_range(start='20140101', freq='T', periods=self.N)
        self.ser = pd.Series(rng)
        self.day_offset = pd.offsets.Day()
        self.relativedelta_offset = pd.offsets.DateOffset(months=2, days=2)
        self.busday_offset = pd.offsets.BusinessDay()

    def time_add_offset_delta(self):
        self.ser + self.day_offset

    def time_add_offset_fast(self):
        self.ser + self.relativedelta_offset

    def time_add_offset_slow(self):
        self.ser + self.busday_offset


class YearBegin(object):
    goal_time = 0.2

    def setup(self):
        self.date = datetime(2011, 1, 1)
        self.year = pd.offsets.YearBegin()

    def time_timeseries_year_apply(self):
        self.year.apply(self.date)

    def time_timeseries_year_incr(self):
        self.date + self.year


class Day(object):
    goal_time = 0.2

    def setup(self):
        self.date = datetime(2011, 1, 1)
        self.day = pd.offsets.Day()

    def time_timeseries_day_apply(self):
        self.day.apply(self.date)

    def time_timeseries_day_incr(self):
        self.date + self.day


class CBDay(object):
    goal_time = 0.2

    def setup(self):
        self.date = datetime(2011, 1, 1)
        self.dt64 = np.datetime64('2011-01-01 09:00Z')
        self.cday = pd.offsets.CustomBusinessDay()

    def time_custom_bday_decr(self):
        self.date - self.cday

    def time_custom_bday_incr(self):
        self.date + self.cday

    def time_custom_bday_apply(self):
        self.cday.apply(self.date)

    def time_custom_bday_apply_dt64(self):
        self.cday.apply(self.dt64)


class CBDayHolidays(object):
    goal_time = 0.2

    def setup(self):
        self.date = datetime(2011, 1, 1)
        self.cdayh = pd.offsets.CustomBusinessDay(calendar=hcal)

    def time_custom_bday_cal_incr(self):
        self.date + 1 * self.cdayh

    def time_custom_bday_cal_decr(self):
        self.date - 1 * self.cdayh

    def time_custom_bday_cal_incr_n(self):
        self.date + 10 * self.cdayh

    def time_custom_bday_cal_incr_neg_n(self):
        self.date - 10 * self.cdayh


class CBMonthBegin(object):
    goal_time = 0.2

    def setup(self):
        self.date = datetime(2011, 1, 1)
        self.cmb = pd.offsets.CustomBusinessMonthBegin(calendar=hcal)

    def time_custom_bmonthbegin_decr_n(self):
        self.date - (10 * self.cmb)

    def time_custom_bmonthbegin_incr_n(self):
        self.date + (10 * self.cmb)


class CBMonthEnd(object):
    goal_time = 0.2

    def setup(self):
        self.date = datetime(2011, 1, 1)
        self.cme = pd.offsets.CustomBusinessMonthEnd(calendar=hcal)

    def time_custom_bmonthend_incr(self):
        self.date + self.cme

    def time_custom_bmonthend_incr_n(self):
        self.date + (10 * self.cme)

    def time_custom_bmonthend_decr_n(self):
        self.date - (10 * self.cme)


class SemiMonthOffset(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')
        # date is not on an offset which will be slowest case
        self.date = datetime(2011, 1, 2)
        self.semi_month_end = pd.offsets.SemiMonthEnd()
        self.semi_month_begin = pd.offsets.SemiMonthBegin()

    def time_end_apply(self):
        self.semi_month_end.apply(self.date)

    def time_end_incr(self):
        self.date + self.semi_month_end

    def time_end_incr_n(self):
        self.date + 10 * self.semi_month_end

    def time_end_decr(self):
        self.date - self.semi_month_end

    def time_end_decr_n(self):
        self.date - 10 * self.semi_month_end

    def time_end_apply_index(self):
        self.semi_month_end.apply_index(self.rng)

    def time_end_incr_rng(self):
        self.rng + self.semi_month_end

    def time_end_decr_rng(self):
        self.rng - self.semi_month_end

    def time_begin_apply(self):
        self.semi_month_begin.apply(self.date)

    def time_begin_incr(self):
        self.date + self.semi_month_begin

    def time_begin_incr_n(self):
        self.date + 10 * self.semi_month_begin

    def time_begin_decr(self):
        self.date - self.semi_month_begin

    def time_begin_decr_n(self):
        self.date - 10 * self.semi_month_begin

    def time_begin_apply_index(self):
        self.semi_month_begin.apply_index(self.rng)

    def time_begin_incr_rng(self):
        self.rng + self.semi_month_begin

    def time_begin_decr_rng(self):
        self.rng - self.semi_month_begin
