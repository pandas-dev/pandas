from datetime import datetime

import numpy as np
import pandas.util.testing as tm
from pandas import Series, date_range, NaT

from .pandas_vb_common import setup  # noqa


class SeriesConstructor(object):

    goal_time = 0.2
    params = [None, 'dict']
    param_names = ['data']

    def setup(self, data):
        self.idx = date_range(start=datetime(2015, 10, 26),
                              end=datetime(2016, 1, 1),
                              freq='50s')
        dict_data = dict(zip(self.idx, range(len(self.idx))))
        self.data = None if data is None else dict_data

    def time_constructor(self, data):
        Series(data=self.data, index=self.idx)


class IsIn(object):

    goal_time = 0.2
    params = ['int64', 'object']
    param_names = ['dtype']

    def setup(self, dtype):
        self.s = Series(np.random.randint(1, 10, 100000)).astype(dtype)
        self.values = [1, 2]

    def time_isin(self, dtypes):
        self.s.isin(self.values)


class IsInForObjects(object):

    def setup(self):
        self.s_nans = Series(np.full(10**4, np.nan)).astype(np.object)
        self.vals_nans = np.full(10**4, np.nan).astype(np.object)
        self.s_short = Series(np.arange(2)).astype(np.object)
        self.s_long = Series(np.arange(10**5)).astype(np.object)
        self.vals_short = np.arange(2).astype(np.object)
        self.vals_long = np.arange(10**5).astype(np.object)
        # because of nans floats are special:
        self.s_long_floats = Series(np.arange(10**5, 
                                       dtype=np.float)).astype(np.object)
        self.vals_long_floats = np.arange(10**5, 
                                       dtype=np.float).astype(np.object)

    def time_isin_nans(self):
        # if nan-objects are different objects,
        # this has the potential to trigger O(n^2) running time
        self.s_nans.isin(self.vals_nans)

    def time_isin_short_series_long_values(self):
        # running time dominated by the preprocessing
        self.s_short.isin(self.vals_long)

    def time_isin_long_series_short_values(self):
        # running time dominated by look-up
        self.s_long.isin(self.vals_short)

    def time_isin_long_series_long_values(self):
        # no dominating part
        self.s_long.isin(self.vals_long)

    def time_isin_long_series_long_values_floats(self):
        # no dominating part
        self.s_long_floats.isin(self.vals_long_floats)


class NSort(object):

    goal_time = 0.2
    params = ['first', 'last', 'all']
    param_names = ['keep']

    def setup(self, keep):
        self.s = Series(np.random.randint(1, 10, 100000))

    def time_nlargest(self, keep):
        self.s.nlargest(3, keep=keep)

    def time_nsmallest(self, keep):
        self.s.nsmallest(3, keep=keep)


class Dropna(object):

    goal_time = 0.2
    params = ['int', 'datetime']
    param_names = ['dtype']

    def setup(self, dtype):
        N = 10**6
        data = {'int': np.random.randint(1, 10, N),
                'datetime': date_range('2000-01-01', freq='S', periods=N)}
        self.s = Series(data[dtype])
        if dtype == 'datetime':
            self.s[np.random.randint(1, N, 100)] = NaT

    def time_dropna(self, dtype):
        self.s.dropna()


class Map(object):

    goal_time = 0.2
    params = ['dict', 'Series']
    param_names = 'mapper'

    def setup(self, mapper):
        map_size = 1000
        map_data = Series(map_size - np.arange(map_size))
        self.map_data = map_data if mapper == 'Series' else map_data.to_dict()
        self.s = Series(np.random.randint(0, map_size, 10000))

    def time_map(self, mapper):
        self.s.map(self.map_data)


class Clip(object):

    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randn(50))

    def time_clip(self):
        self.s.clip(0, 1)


class ValueCounts(object):

    goal_time = 0.2
    params = ['int', 'float', 'object']
    param_names = ['dtype']

    def setup(self, dtype):
        self.s = Series(np.random.randint(0, 1000, size=100000)).astype(dtype)

    def time_value_counts(self, dtype):
        self.s.value_counts()


class Dir(object):

    goal_time = 0.2

    def setup(self):
        self.s = Series(index=tm.makeStringIndex(10000))

    def time_dir_strings(self):
        dir(self.s)


class SeriesGetattr(object):
    # https://github.com/pandas-dev/pandas/issues/19764
    goal_time = 0.2

    def setup(self):
        self.s = Series(1,
                        index=date_range("2012-01-01", freq='s',
                                         periods=int(1e6)))

    def time_series_datetimeindex_repr(self):
        getattr(self.s, 'a', None)
