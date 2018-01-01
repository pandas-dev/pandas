import pandas as pd
import numpy as np

from .pandas_vb_common import setup  # noqa


class Methods(object):

    sample_time = 0.2
    params = (['DataFrame', 'Series'],
              [10, 1000],
              ['int', 'float'],
              ['median', 'mean', 'max', 'min', 'std', 'count', 'skew', 'kurt',
               'sum', 'corr', 'cov'])
    param_names = ['contructor', 'window', 'dtype', 'method']

    def setup(self, contructor, window, dtype, method):
        N = 10**5
        arr = np.random.random(N).astype(dtype)
        self.roll = getattr(pd, contructor)(arr).rolling(window)

    def time_rolling(self, contructor, window, dtype, method):
        getattr(self.roll, method)()


class Quantile(object):

    sample_time = 0.2
    params = (['DataFrame', 'Series'],
              [10, 1000],
              ['int', 'float'],
              [0, 0.5, 1])
    param_names = ['contructor', 'window', 'dtype', 'percentile']

    def setup(self, contructor, window, dtype, percentile):
        N = 10**5
        arr = np.random.random(N).astype(dtype)
        self.roll = getattr(pd, contructor)(arr).rolling(window)

    def time_quantile(self, contructor, window, dtype, percentile):
        self.roll.quantile(percentile)


class DepreciatedRolling(object):

    sample_time = 0.2
    params = ['rolling_median', 'rolling_mean', 'rolling_min', 'rolling_max',
              'rolling_var', 'rolling_skew', 'rolling_kurt', 'rolling_std']
    param_names = ['method']

    def setup(self, method):
        self.arr = np.random.randn(100000)
        self.win = 100

    def time_method(self, method):
        getattr(pd, method)(self.arr, self.win)
