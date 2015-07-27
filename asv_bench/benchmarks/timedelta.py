from pandas_vb_common import *
from pandas import to_timedelta


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