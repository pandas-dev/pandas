import numpy as np
import pandas as pd

from .pandas_vb_common import setup  # noqa


class FrameOps(object):

    goal_time = 0.2
    param_names = ['op', 'use_bottleneck', 'dtype', 'axis']
    params = [['mean', 'sum', 'median', 'std'],
              [True, False],
              ['float', 'int'],
              [0, 1]]

    def setup(self, op, use_bottleneck, dtype, axis):
        df = pd.DataFrame(np.random.randn(100000, 4)).astype(dtype)
        try:
            pd.options.compute.use_bottleneck = use_bottleneck
        except:
            from pandas.core import nanops
            nanops._USE_BOTTLENECK = use_bottleneck
        self.df_func = getattr(df, op)

    def time_op(self, op, use_bottleneck, dtype, axis):
        self.df_func(axis=axis)


class FrameMultiIndexOps(object):

    goal_time = 0.2
    params = ([0, 1, [0, 1]], ['mean', 'sum', 'median'])
    param_names = ['level', 'op']

    def setup(self, level, op):
        levels = [np.arange(10), np.arange(100), np.arange(100)]
        labels = [np.arange(10).repeat(10000),
                  np.tile(np.arange(100).repeat(100), 10),
                  np.tile(np.tile(np.arange(100), 100), 10)]
        index = pd.MultiIndex(levels=levels, labels=labels)
        df = pd.DataFrame(np.random.randn(len(index), 4), index=index)
        self.df_func = getattr(df, op)

    def time_op(self, level, op):
        self.df_func(level=level)


class SeriesOps(object):

    goal_time = 0.2
    param_names = ['op', 'use_bottleneck', 'dtype']
    params = [['mean', 'sum', 'median', 'std'],
              [True, False],
              ['float', 'int']]

    def setup(self, op, use_bottleneck, dtype):
        s = pd.Series(np.random.randn(100000)).astype(dtype)
        try:
            pd.options.compute.use_bottleneck = use_bottleneck
        except:
            from pandas.core import nanops
            nanops._USE_BOTTLENECK = use_bottleneck
        self.s_func = getattr(s, op)

    def time_op(self, op, use_bottleneck, dtype):
        self.s_func()


class SeriesMultiIndexOps(object):

    goal_time = 0.2
    params = ([0, 1, [0, 1]], ['mean', 'sum', 'median'])
    param_names = ['level', 'op']

    def setup(self, level, op):
        levels = [np.arange(10), np.arange(100), np.arange(100)]
        labels = [np.arange(10).repeat(10000),
                  np.tile(np.arange(100).repeat(100), 10),
                  np.tile(np.tile(np.arange(100), 100), 10)]
        index = pd.MultiIndex(levels=levels, labels=labels)
        s = pd.Series(np.random.randn(len(index)), index=index)
        self.s_func = getattr(s, op)

    def time_op(self, level, op):
        self.s_func(level=level)


class Rank(object):

    goal_time = 0.2
    params = [['DataFrame', 'Series'], [True, False]]
    param_names = ['constructor', 'pct']

    def setup(self, constructor, pct):
        values = np.random.randn(10**5)
        self.data = getattr(pd, constructor)(values)

    def time_rank(self, constructor, pct):
        self.data.rank(pct=pct)

    def time_average_old(self, constructor, pct):
        self.data.rank(pct=pct) / len(self.data)


class Correlation(object):

    goal_time = 0.2
    params = ['spearman', 'kendall', 'pearson']
    param_names = ['method']

    def setup(self, method):
        self.df = pd.DataFrame(np.random.randn(1000, 30))

    def time_corr(self, method):
        self.df.corr(method=method)
