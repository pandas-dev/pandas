import warnings
from importlib import import_module

import numpy as np
import pandas as pd
from pandas.util import testing as tm

for imp in ['pandas.util', 'pandas.tools.hashing']:
    try:
        hashing = import_module(imp)
        break
    except (ImportError, TypeError, ValueError):
        pass

from .pandas_vb_common import setup # noqa


class Factorize(object):

    goal_time = 0.2

    params = [True, False]
    param_names = ['sort']

    def setup(self, sort):
        N = 10**5
        self.int_idx = pd.Int64Index(np.arange(N).repeat(5))
        self.float_idx = pd.Float64Index(np.random.randn(N).repeat(5))
        self.string_idx = tm.makeStringIndex(N)

    def time_factorize_int(self, sort):
        self.int_idx.factorize(sort=sort)

    def time_factorize_float(self, sort):
        self.float_idx.factorize(sort=sort)

    def time_factorize_string(self, sort):
        self.string_idx.factorize(sort=sort)


class Duplicated(object):

    goal_time = 0.2

    params = ['first', 'last', False]
    param_names = ['keep']

    def setup(self, keep):
        N = 10**5
        self.int_idx = pd.Int64Index(np.arange(N).repeat(5))
        self.float_idx = pd.Float64Index(np.random.randn(N).repeat(5))
        self.string_idx = tm.makeStringIndex(N)

    def time_duplicated_int(self, keep):
        self.int_idx.duplicated(keep=keep)

    def time_duplicated_float(self, keep):
        self.float_idx.duplicated(keep=keep)

    def time_duplicated_string(self, keep):
        self.string_idx.duplicated(keep=keep)


class DuplicatedUniqueIndex(object):

    goal_time = 0.2

    def setup(self):
        N = 10**5
        self.idx_int_dup = pd.Int64Index(np.arange(N * 5))
        # cache is_unique
        self.idx_int_dup.is_unique

    def time_duplicated_unique_int(self):
        self.idx_int_dup.duplicated()


class Match(object):

    goal_time = 0.2

    def setup(self):
        self.uniques = tm.makeStringIndex(1000).values
        self.all = self.uniques.repeat(10)

    def time_match_string(self):
        with warnings.catch_warnings(record=True):
            pd.match(self.all, self.uniques)


class Hashing(object):

    goal_time = 0.2

    def setup_cache(self):
        N = 10**5

        df = pd.DataFrame(
            {'strings': pd.Series(tm.makeStringIndex(10000).take(
                np.random.randint(0, 10000, size=N))),
             'floats': np.random.randn(N),
             'ints': np.arange(N),
             'dates': pd.date_range('20110101', freq='s', periods=N),
             'timedeltas': pd.timedelta_range('1 day', freq='s', periods=N)})
        df['categories'] = df['strings'].astype('category')
        df.iloc[10:20] = np.nan
        return df

    def time_frame(self, df):
        hashing.hash_pandas_object(df)

    def time_series_int(self, df):
        hashing.hash_pandas_object(df['ints'])

    def time_series_string(self, df):
        hashing.hash_pandas_object(df['strings'])

    def time_series_float(self, df):
        hashing.hash_pandas_object(df['floats'])

    def time_series_categorical(self, df):
        hashing.hash_pandas_object(df['categories'])

    def time_series_timedeltas(self, df):
        hashing.hash_pandas_object(df['timedeltas'])

    def time_series_dates(self, df):
        hashing.hash_pandas_object(df['dates'])
