import pandas as pd
import numpy as np

from .pandas_vb_common import setup


class Methods(object):

    sample_time = 0.2
    params = (['DataFrame', 'Series'],
              [10, 1000],
              [10**4, 10**5],
              ['int', 'float'],
              ['median', 'mean', 'max', 'min', 'std', 'count', 'skew', 'kurt',
               'sum', 'corr', 'cov'])
    param_names = ['contructor', 'window', 'num_data', 'dtype', 'method']

    def setup(self, contructor, window, num_data, dtype, method):
        arr = np.random.random(num_data).astype(dtype)
        self.data = getattr(pd, contructor)(arr)

    def time_rolling(self, contructor, window, num_data, dtype, method):
        getattr(self.data.rolling(window), method)()


class Quantile(object):

    sample_time = 0.2
    params = (['DataFrame', 'Series'],
              [10, 1000],
              [10**4, 10**5],
              [0, 0.5, 1])
    param_names = ['contructor', 'window', 'num_data', 'dtype', 'percentile']

    def setup(self, contructor, window, num_data, dtype, percentile):
        arr = np.random.random(num_data).astype(dtype)
        self.data = getattr(pd, contructor)(arr)

    def time_quantile(self, contructor, window, num_data, dtype, percentile):
        self.data.rolling(window).quantile(0.5)
