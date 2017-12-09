import numpy as np
import pandas as pd

from pandas import to_timedelta, Timestamp, Timedelta


class ToTimedelta(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randint(0, 1000, size=10000)
        self.arr2 = ['{0} days'.format(i) for i in self.arr]

        self.arr3 = np.random.randint(0, 60, size=10000)
        self.arr3 = ['00:00:{0:02d}'.format(i) for i in self.arr3]

        self.arr4 = list(self.arr2)
        self.arr4[-1] = 'apple'

    def time_convert_int(self):
        to_timedelta(self.arr, unit='s')

    def time_convert_string(self):
        to_timedelta(self.arr2)

    def time_convert_string_seconds(self):
        to_timedelta(self.arr3)

    def time_convert_coerce(self):
        to_timedelta(self.arr4, errors='coerce')

    def time_convert_ignore(self):
        to_timedelta(self.arr4, errors='ignore')


class TimedeltaOps(object):
    goal_time = 0.2

    def setup(self):
        self.td = to_timedelta(np.arange(1000000))
        self.ts = Timestamp('2000')

    def time_add_td_ts(self):
        self.td + self.ts


class TimedeltaProperties(object):
    goal_time = 0.2

    def setup(self):
        self.td = Timedelta(days=365, minutes=35, seconds=25, milliseconds=35)

    def time_timedelta_days(self):
        self.td.days

    def time_timedelta_seconds(self):
        self.td.seconds

    def time_timedelta_microseconds(self):
        self.td.microseconds

    def time_timedelta_nanoseconds(self):
        self.td.nanoseconds


class DatetimeAccessor(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.series = pd.Series(
            pd.timedelta_range('1 days', periods=self.N, freq='h'))

    def time_dt_accessor(self):
        self.series.dt

    def time_timedelta_dt_accessor_days(self):
        self.series.dt.days

    def time_timedelta_dt_accessor_seconds(self):
        self.series.dt.seconds

    def time_timedelta_dt_accessor_microseconds(self):
        self.series.dt.microseconds

    def time_timedelta_dt_accessor_nanoseconds(self):
        self.series.dt.nanoseconds
