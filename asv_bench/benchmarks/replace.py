from pandas_vb_common import *
from pandas.compat import range
from datetime import timedelta


class replace_fillna(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        try:
            self.rng = date_range('1/1/2000', periods=self.N, freq='min')
        except NameError:
            self.rng = DatetimeIndex('1/1/2000', periods=self.N, offset=datetools.Minute())
            self.date_range = DateRange
        self.ts = Series(np.random.randn(self.N), index=self.rng)

    def time_replace_fillna(self):
        self.ts.fillna(0.0, inplace=True)


class replace_large_dict(object):
    goal_time = 0.2

    def setup(self):
        self.n = (10 ** 6)
        self.start_value = (10 ** 5)
        self.to_rep = dict(((i, (self.start_value + i)) for i in range(self.n)))
        self.s = Series(np.random.randint(self.n, size=(10 ** 3)))

    def time_replace_large_dict(self):
        self.s.replace(self.to_rep, inplace=True)


class replace_replacena(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        try:
            self.rng = date_range('1/1/2000', periods=self.N, freq='min')
        except NameError:
            self.rng = DatetimeIndex('1/1/2000', periods=self.N, offset=datetools.Minute())
            self.date_range = DateRange
        self.ts = Series(np.random.randn(self.N), index=self.rng)

    def time_replace_replacena(self):
        self.ts.replace(np.nan, 0.0, inplace=True)