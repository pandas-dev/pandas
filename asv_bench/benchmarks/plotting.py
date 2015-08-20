from pandas_vb_common import *
try:
    from pandas import date_range
except ImportError:

    def date_range(start=None, end=None, periods=None, freq=None):
        return DatetimeIndex(start, end, periods=periods, offset=freq)


class plot_timeseries_period(object):
    goal_time = 0.2

    def setup(self):
        self.N = 2000
        self.M = 5
        self.df = DataFrame(np.random.randn(self.N, self.M), index=date_range('1/1/1975', periods=self.N))

    def time_plot_timeseries_period(self):
        self.df.plot()