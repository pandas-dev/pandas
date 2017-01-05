from .pandas_vb_common import *
from pandas import to_timedelta, Timestamp


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


class Ops(object):
    goal_time = 0.2

    def setup(self):
        self.td = to_timedelta(np.arange(1000000))
        self.ts = Timestamp('2000')

    def test_add_td_ts(self):
        self.td + self.ts
