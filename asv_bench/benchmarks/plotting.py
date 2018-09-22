import numpy as np
from pandas import DataFrame, Series, DatetimeIndex, date_range
try:
    from pandas.plotting import andrews_curves
except ImportError:
    from pandas.tools.plotting import andrews_curves
import matplotlib
matplotlib.use('Agg')

from .pandas_vb_common import setup  # noqa


class Plotting(object):

    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randn(1000000))
        self.df = DataFrame({'col': self.s})

    def time_series_plot(self):
        self.s.plot()

    def time_frame_plot(self):
        self.df.plot()


class TimeseriesPlotting(object):

    goal_time = 0.2

    def setup(self):
        N = 2000
        M = 5
        idx = date_range('1/1/1975', periods=N)
        self.df = DataFrame(np.random.randn(N, M), index=idx)

        idx_irregular = DatetimeIndex(np.concatenate((idx.values[0:10],
                                                      idx.values[12:])))
        self.df2 = DataFrame(np.random.randn(len(idx_irregular), M),
                             index=idx_irregular)

    def time_plot_regular(self):
        self.df.plot()

    def time_plot_regular_compat(self):
        self.df.plot(x_compat=True)

    def time_plot_irregular(self):
        self.df2.plot()


class Misc(object):

    goal_time = 0.6

    def setup(self):
        N = 500
        M = 10
        self.df = DataFrame(np.random.randn(N, M))
        self.df['Name'] = ["A"] * N

    def time_plot_andrews_curves(self):
        andrews_curves(self.df, "Name")
