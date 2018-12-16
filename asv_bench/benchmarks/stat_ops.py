import numpy as np
import pandas as pd


ops = ['mean', 'sum', 'median', 'std', 'skew', 'kurt', 'mad', 'prod', 'sem',
       'var']


class FrameOps(object):

    params = [ops, ['float', 'int'], [0, 1]]
    param_names = ['op', 'dtype', 'axis']

    def setup(self, op, dtype, axis):
        df = pd.DataFrame(np.random.randn(100000, 4)).astype(dtype)
        self.df_func = getattr(df, op)

    def time_op(self, op, dtype, axis):
        self.df_func(axis=axis)


class FrameMultiIndexOps(object):

    params = ([0, 1, [0, 1]], ops)
    param_names = ['level', 'op']

    def setup(self, level, op):
        levels = [np.arange(10), np.arange(100), np.arange(100)]
        codes = [np.arange(10).repeat(10000),
                 np.tile(np.arange(100).repeat(100), 10),
                 np.tile(np.tile(np.arange(100), 100), 10)]
        index = pd.MultiIndex(levels=levels, codes=codes)
        df = pd.DataFrame(np.random.randn(len(index), 4), index=index)
        self.df_func = getattr(df, op)

    def time_op(self, level, op):
        self.df_func(level=level)


class SeriesOps(object):

    params = [ops, ['float', 'int']]
    param_names = ['op', 'dtype']

    def setup(self, op, dtype):
        s = pd.Series(np.random.randn(100000)).astype(dtype)
        self.s_func = getattr(s, op)

    def time_op(self, op, dtype):
        self.s_func()


class SeriesMultiIndexOps(object):

    params = ([0, 1, [0, 1]], ops)
    param_names = ['level', 'op']

    def setup(self, level, op):
        levels = [np.arange(10), np.arange(100), np.arange(100)]
        codes = [np.arange(10).repeat(10000),
                 np.tile(np.arange(100).repeat(100), 10),
                 np.tile(np.tile(np.arange(100), 100), 10)]
        index = pd.MultiIndex(levels=levels, codes=codes)
        s = pd.Series(np.random.randn(len(index)), index=index)
        self.s_func = getattr(s, op)

    def time_op(self, level, op):
        self.s_func(level=level)


class Rank(object):

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

    params = [['spearman', 'kendall', 'pearson']]
    param_names = ['method']

    def setup(self, method):
        self.df = pd.DataFrame(np.random.randn(1000, 30))
        self.s = pd.Series(np.random.randn(1000))
        self.s2 = pd.Series(np.random.randn(1000))

    def time_corr(self, method):
        self.df.corr(method=method)

    def time_corr_series(self, method):
        self.s.corr(self.s2, method=method)


class Covariance(object):

    def setup(self):
        self.s = pd.Series(np.random.randn(100000))
        self.s2 = pd.Series(np.random.randn(100000))

    def time_cov_series(self):
        self.s.cov(self.s2)


from .pandas_vb_common import setup  # noqa: F401
