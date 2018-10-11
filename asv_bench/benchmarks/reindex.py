import numpy as np
import pandas.util.testing as tm
from pandas import (DataFrame, Series, DatetimeIndex, MultiIndex, Index,
                    date_range)
from .pandas_vb_common import lib


class Reindex(object):

    goal_time = 0.2

    def setup(self):
        rng = DatetimeIndex(start='1/1/1970', periods=10000, freq='1min')
        self.df = DataFrame(np.random.rand(10000, 10), index=rng,
                            columns=range(10))
        self.df['foo'] = 'bar'
        self.rng_subset = Index(rng[::2])
        self.df2 = DataFrame(index=range(10000),
                             data=np.random.rand(10000, 30), columns=range(30))
        N = 5000
        K = 200
        level1 = tm.makeStringIndex(N).values.repeat(K)
        level2 = np.tile(tm.makeStringIndex(K).values, N)
        index = MultiIndex.from_arrays([level1, level2])
        self.s = Series(np.random.randn(N * K), index=index)
        self.s_subset = self.s[::2]

    def time_reindex_dates(self):
        self.df.reindex(self.rng_subset)

    def time_reindex_columns(self):
        self.df2.reindex(columns=self.df.columns[1:5])

    def time_reindex_multiindex(self):
        self.s.reindex(self.s_subset.index)


class ReindexMethod(object):

    goal_time = 0.2
    params = ['pad', 'backfill']
    param_names = ['method']

    def setup(self, method):
        N = 100000
        self.idx = date_range('1/1/2000', periods=N, freq='1min')
        self.ts = Series(np.random.randn(N), index=self.idx)[::2]

    def time_reindex_method(self, method):
        self.ts.reindex(self.idx, method=method)


class Fillna(object):

    goal_time = 0.2
    params = ['pad', 'backfill']
    param_names = ['method']

    def setup(self, method):
        N = 100000
        self.idx = date_range('1/1/2000', periods=N, freq='1min')
        ts = Series(np.random.randn(N), index=self.idx)[::2]
        self.ts_reindexed = ts.reindex(self.idx)
        self.ts_float32 = self.ts_reindexed.astype('float32')

    def time_reindexed(self, method):
        self.ts_reindexed.fillna(method=method)

    def time_float_32(self, method):
        self.ts_float32.fillna(method=method)


class LevelAlign(object):

    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(
            levels=[np.arange(10), np.arange(100), np.arange(100)],
            labels=[np.arange(10).repeat(10000),
                    np.tile(np.arange(100).repeat(100), 10),
                    np.tile(np.tile(np.arange(100), 100), 10)])
        self.df = DataFrame(np.random.randn(len(self.index), 4),
                            index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4),
                                  index=self.index.levels[1])

    def time_align_level(self):
        self.df.align(self.df_level, level=1, copy=False)

    def time_reindex_level(self):
        self.df_level.reindex(self.index, level=1)


class DropDuplicates(object):

    goal_time = 0.2
    params = [True, False]
    param_names = ['inplace']

    def setup(self, inplace):
        N = 10000
        K = 10
        key1 = tm.makeStringIndex(N).values.repeat(K)
        key2 = tm.makeStringIndex(N).values.repeat(K)
        self.df = DataFrame({'key1': key1, 'key2': key2,
                             'value': np.random.randn(N * K)})
        self.df_nan = self.df.copy()
        self.df_nan.iloc[:10000, :] = np.nan

        self.s = Series(np.random.randint(0, 1000, size=10000))
        self.s_str = Series(np.tile(tm.makeStringIndex(1000).values, 10))

        N = 1000000
        K = 10000
        key1 = np.random.randint(0, K, size=N)
        self.df_int = DataFrame({'key1': key1})
        self.df_bool = DataFrame(np.random.randint(0, 2, size=(K, 10),
                                                   dtype=bool))

    def time_frame_drop_dups(self, inplace):
        self.df.drop_duplicates(['key1', 'key2'], inplace=inplace)

    def time_frame_drop_dups_na(self, inplace):
        self.df_nan.drop_duplicates(['key1', 'key2'], inplace=inplace)

    def time_series_drop_dups_int(self, inplace):
        self.s.drop_duplicates(inplace=inplace)

    def time_series_drop_dups_string(self, inplace):
        self.s_str.drop_duplicates(inplace=inplace)

    def time_frame_drop_dups_int(self, inplace):
        self.df_int.drop_duplicates(inplace=inplace)

    def time_frame_drop_dups_bool(self, inplace):
        self.df_bool.drop_duplicates(inplace=inplace)


class Align(object):
    # blog "pandas escaped the zoo"
    goal_time = 0.2

    def setup(self):
        n = 50000
        indices = tm.makeStringIndex(n)
        subsample_size = 40000
        self.x = Series(np.random.randn(n), indices)
        self.y = Series(np.random.randn(subsample_size),
                        index=np.random.choice(indices, subsample_size,
                                               replace=False))

    def time_align_series_irregular_string(self):
        self.x + self.y


class LibFastZip(object):

    goal_time = 0.2

    def setup(self):
        N = 10000
        K = 10
        key1 = tm.makeStringIndex(N).values.repeat(K)
        key2 = tm.makeStringIndex(N).values.repeat(K)
        col_array = np.vstack([key1, key2, np.random.randn(N * K)])
        col_array2 = col_array.copy()
        col_array2[:, :10000] = np.nan
        self.col_array_list = list(col_array)

    def time_lib_fast_zip(self):
        lib.fast_zip(self.col_array_list)


from .pandas_vb_common import setup  # noqa: F401
