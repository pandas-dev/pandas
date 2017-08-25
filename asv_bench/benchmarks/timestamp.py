from .pandas_vb_common import *
from pandas import to_timedelta, Timestamp


class TimestampProperties(object):
    goal_time = 0.2

    def setup(self):
        self.ts = Timestamp('2017-08-25 08:16:14')

    def time_tz(self):
        self.ts.tz

    def time_offset(self):
        self.ts.offset

    def time_dayofweek(self):
        self.ts.dayofweek

    def time_weekday_name(self):
        self.ts.weekday_name

    def time_dayofyear(self):
        self.ts.dayofyear

    def time_week(self):
        self.ts.week

    def time_quarter(self):
        self.ts.quarter

    def time_days_in_month(self):
        self.ts.days_in_month

    def time_freqstr(self):
        self.ts.freqstr

    def time_is_month_start(self):
        self.ts.is_month_start

    def time_is_month_end(self):
        self.ts.is_month_end

    def time_is_quarter_start(self):
        self.ts.is_quarter_start

    def time_is_quarter_end(self):
        self.ts.is_quarter_end

    def time_is_year_start(self):
        self.ts.is_quarter_end

    def time_is_year_end(self):
        self.ts.is_quarter_end

    def time_is_leap_year(self):
        self.ts.is_quarter_end

    def time_microsecond(self):
        self.ts.microsecond
