import numpy as np
import pandas as pd


ops = ['mean', 'sum', 'median', 'std', 'skew', 'kurt', 'mad', 'prod', 'sem',
       'var']


class FrameOps:

    params = [ops, ['float', 'int'], [0, 1], [True, False]]
    param_names = ['op', 'dtype', 'axis', 'use_bottleneck']

    def setup(self, op, dtype, axis, use_bottleneck):
        df = pd.DataFrame(np.random.randn(100000, 4)).astype(dtype)
        try:
            pd.options.compute.use_bottleneck = use_bottleneck
        except TypeError:
            from pandas.core import nanops
            nanops._USE_BOTTLENECK = use_bottleneck
        self.df_func = getattr(df, op)

    def time_op(self, op, dtype, axis, use_bottleneck):
        self.df_func(axis=axis)


class FrameMultiIndexOps:

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


class SeriesOps:

    params = [ops, ['float', 'int'], [True, False]]
    param_names = ['op', 'dtype', 'use_bottleneck']

    def setup(self, op, dtype, use_bottleneck):
        s = pd.Series(np.random.randn(100000)).astype(dtype)
        try:
            pd.options.compute.use_bottleneck = use_bottleneck
        except TypeError:
            from pandas.core import nanops
            nanops._USE_BOTTLENECK = use_bottleneck
        self.s_func = getattr(s, op)

    def time_op(self, op, dtype, use_bottleneck):
        self.s_func()


class SeriesMultiIndexOps:

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


class Rank:

    params = [['DataFrame', 'Series'], [True, False]]
    param_names = ['constructor', 'pct']

    def setup(self, constructor, pct):
        values = np.random.randn(10**5)
        self.data = getattr(pd, constructor)(values)

    def time_rank(self, constructor, pct):
        self.data.rank(pct=pct)

    def time_average_old(self, constructor, pct):
        self.data.rank(pct=pct) / len(self.data)


class Correlation:

    params = [['spearman', 'kendall', 'pearson'], [True, False]]
    param_names = ['method', 'use_bottleneck']

    def setup(self, method, use_bottleneck):
        try:
            pd.options.compute.use_bottleneck = use_bottleneck
        except TypeError:
            from pandas.core import nanops
            nanops._USE_BOTTLENECK = use_bottleneck
        self.df = pd.DataFrame(np.random.randn(1000, 30))
        self.df2 = pd.DataFrame(np.random.randn(1000, 30))
        self.s = pd.Series(np.random.randn(1000))
        self.s2 = pd.Series(np.random.randn(1000))

    def time_corr(self, method, use_bottleneck):
        self.df.corr(method=method)

    def time_corr_series(self, method, use_bottleneck):
        self.s.corr(self.s2, method=method)

    def time_corrwith_cols(self, method, use_bottleneck):
        self.df.corrwith(self.df2, method=method)

    def time_corrwith_rows(self, method, use_bottleneck):
        self.df.corrwith(self.df2, axis=1, method=method)


class Covariance:

    params = [[True, False]]
    param_names = ['use_bottleneck']

    def setup(self, use_bottleneck):
        try:
            pd.options.compute.use_bottleneck = use_bottleneck
        except TypeError:
            from pandas.core import nanops
            nanops._USE_BOTTLENECK = use_bottleneck
        self.s = pd.Series(np.random.randn(100000))
        self.s2 = pd.Series(np.random.randn(100000))

    def time_cov_series(self, use_bottleneck):
        self.s.cov(self.s2)


from .pandas_vb_common import setup  # noqa: F401
