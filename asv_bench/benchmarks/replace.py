import numpy as np
import pandas as pd


class FillNa(object):

    params = [True, False]
    param_names = ['inplace']

    def setup(self, inplace):
        N = 10**6
        rng = pd.date_range('1/1/2000', periods=N, freq='min')
        data = np.random.randn(N)
        data[::2] = np.nan
        self.ts = pd.Series(data, index=rng)

    def time_fillna(self, inplace):
        self.ts.fillna(0.0, inplace=inplace)

    def time_replace(self, inplace):
        self.ts.replace(np.nan, 0.0, inplace=inplace)


class ReplaceDict(object):

    params = [True, False]
    param_names = ['inplace']

    def setup(self, inplace):
        N = 10**5
        start_value = 10**5
        self.to_rep = dict(enumerate(np.arange(N) + start_value))
        self.s = pd.Series(np.random.randint(N, size=10**3))

    def time_replace_series(self, inplace):
        self.s.replace(self.to_rep, inplace=inplace)


class Convert(object):

    params = (['DataFrame', 'Series'], ['Timestamp', 'Timedelta'])
    param_names = ['constructor', 'replace_data']

    def setup(self, constructor, replace_data):
        N = 10**3
        data = {'Series': pd.Series(np.random.randint(N, size=N)),
                'DataFrame': pd.DataFrame({'A': np.random.randint(N, size=N),
                                           'B': np.random.randint(N, size=N)})}
        self.to_replace = {i: getattr(pd, replace_data) for i in range(N)}
        self.data = data[constructor]

    def time_replace(self, constructor, replace_data):
        self.data.replace(self.to_replace)


from .pandas_vb_common import setup  # noqa: F401
