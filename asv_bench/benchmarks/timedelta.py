from .pandas_vb_common import *
from pandas import to_timedelta, Timestamp


class timedelta_convert_int(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randint(0, 1000, size=10000)

    def time_timedelta_convert_int(self):
        to_timedelta(self.arr, unit='s')


class timedelta_convert_string(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randint(0, 1000, size=10000)
        self.arr = ['{0} days'.format(i) for i in self.arr]

    def time_timedelta_convert_string(self):
        to_timedelta(self.arr)


class timedelta_convert_string_seconds(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randint(0, 60, size=10000)
        self.arr = ['00:00:{0:02d}'.format(i) for i in self.arr]

    def time_timedelta_convert_string_seconds(self):
        to_timedelta(self.arr)


class timedelta_convert_bad_parse(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randint(0, 1000, size=10000)
        self.arr = ['{0} days'.format(i) for i in self.arr]
        self.arr[-1] = 'apple'

    def time_timedelta_convert_coerce(self):
        to_timedelta(self.arr, errors='coerce')

    def time_timedelta_convert_ignore(self):
        to_timedelta(self.arr, errors='ignore')


class timedelta_add_overflow(object):
    goal_time = 0.2

    def setup(self):
        self.td = to_timedelta(np.arange(1000000))
        self.ts = Timestamp('2000')

    def test_add_td_ts(self):
        self.td + self.ts
