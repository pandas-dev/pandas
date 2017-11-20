from .pandas_vb_common import *
try:
    from pandas import date_range
except ImportError:
    def date_range(start=None, end=None, periods=None, freq=None):
        return DatetimeIndex(start, end, periods=periods, offset=freq)
try:
    from pandas.plotting import andrews_curves
except ImportError:
    from pandas.tools.plotting import andrews_curves


class Plotting(object):
    goal_time = 0.2

    def setup(self):
        import matplotlib
        matplotlib.use('Agg')
        self.s = Series(np.random.randn(1000000))
        self.df = DataFrame({'col': self.s})

    def time_series_plot(self):
        self.s.plot()

    def time_frame_plot(self):
        self.df.plot()


class TimeseriesPlotting(object):
    goal_time = 0.2

    def setup(self):
        import matplotlib
        matplotlib.use('Agg')
        N = 2000
        M = 5
        idx = date_range('1/1/1975', periods=N)
        self.df = DataFrame(np.random.randn(N, M), index=idx)

        idx_irregular = pd.DatetimeIndex(np.concatenate((idx.values[0:10],
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
        import matplotlib
        matplotlib.use('Agg')
        self.N = 500
        self.M = 10
        data_dict = {x: np.random.randn(self.N) for x in range(self.M)}
        data_dict["Name"] = ["A"] * self.N
        self.df = DataFrame(data_dict)

    def time_plot_andrews_curves(self):
        andrews_curves(self.df, "Name")
