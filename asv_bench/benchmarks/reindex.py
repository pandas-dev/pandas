from .pandas_vb_common import *
from random import shuffle


class Reindexing(object):
    goal_time = 0.2

    def setup(self):
        self.rng = DatetimeIndex(start='1/1/1970', periods=10000, freq='1min')
        self.df = DataFrame(np.random.rand(10000, 10), index=self.rng,
                            columns=range(10))
        self.df['foo'] = 'bar'
        self.rng2 = Index(self.rng[::2])

        self.df2 = DataFrame(index=range(10000),
                             data=np.random.rand(10000, 30), columns=range(30))

        # multi-index
        N = 5000
        K = 200
        level1 = tm.makeStringIndex(N).values.repeat(K)
        level2 = np.tile(tm.makeStringIndex(K).values, N)
        index = MultiIndex.from_arrays([level1, level2])
        self.s1 = Series(np.random.randn((N * K)), index=index)
        self.s2 = self.s1[::2]

    def time_reindex_dates(self):
        self.df.reindex(self.rng2)

    def time_reindex_columns(self):
        self.df2.reindex(columns=self.df.columns[1:5])

    def time_reindex_multiindex(self):
        self.s1.reindex(self.s2.index)


#----------------------------------------------------------------------
# Pad / backfill


class FillMethod(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq='1min')
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

    def pad(self, source_series, target_index):
        try:
            source_series.reindex(target_index, method='pad')
        except:
            source_series.reindex(target_index, fillMethod='pad')

    def backfill(self, source_series, target_index):
        try:
            source_series.reindex(target_index, method='backfill')
        except:
            source_series.reindex(target_index, fillMethod='backfill')

    def time_backfill_dates(self):
        self.backfill(self.ts2, self.ts.index)

    def time_pad_daterange(self):
        self.pad(self.ts2, self.ts.index)

    def time_backfill(self):
        self.ts3.fillna(method='backfill')

    def time_backfill_float32(self):
        self.ts4.fillna(method='backfill')

    def time_pad(self):
        self.ts3.fillna(method='pad')

    def time_pad_float32(self):
        self.ts4.fillna(method='pad')


#----------------------------------------------------------------------
# align on level


class LevelAlign(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(
            levels=[np.arange(10), np.arange(100), np.arange(100)],
            labels=[np.arange(10).repeat(10000),
                    np.tile(np.arange(100).repeat(100), 10),
                    np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4),
                            index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4),
                                  index=self.index.levels[1])

    def time_align_level(self):
        self.df.align(self.df_level, level=1, copy=False)

    def time_reindex_level(self):
        self.df_level.reindex(self.df.index, level=1)


#----------------------------------------------------------------------
# drop_duplicates


class Duplicates(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2,
                             'value': np.random.randn((self.N * self.K)),})
        self.col_array_list = list(self.df.values.T)

        self.df2 = self.df.copy()
        self.df2.ix[:10000, :] = np.nan

        self.s = Series(np.random.randint(0, 1000, size=10000))
        self.s2 = Series(np.tile(tm.makeStringIndex(1000).values, 10))

        np.random.seed(1234)
        self.N = 1000000
        self.K = 10000
        self.key1 = np.random.randint(0, self.K, size=self.N)
        self.df_int = DataFrame({'key1': self.key1})

    def time_frame_drop_dups(self):
        self.df.drop_duplicates(['key1', 'key2'])

    def time_frame_drop_dups_inplace(self):
        self.df.drop_duplicates(['key1', 'key2'], inplace=True)

    def time_frame_drop_dups_na(self):
        self.df2.drop_duplicates(['key1', 'key2'])

    def time_frame_drop_dups_na_inplace(self):
        self.df2.drop_duplicates(['key1', 'key2'], inplace=True)

    def time_series_drop_dups_int(self):
        self.s.drop_duplicates()

    def time_series_drop_dups_string(self):
        self.s2.drop_duplicates()

    def time_frame_drop_dups_int(self):
        self.df_int.drop_duplicates()


#----------------------------------------------------------------------
# blog "pandas escaped the zoo"


class Align(object):
    goal_time = 0.2

    def setup(self):
        n = 50000
        indices = tm.makeStringIndex(n)
        subsample_size = 40000

        def sample(values, k):
            sampler = np.arange(len(values))
            shuffle(sampler)
            return values.take(sampler[:k])

        self.x = Series(np.random.randn(50000), indices)
        self.y = Series(np.random.randn(subsample_size),
                        index=sample(indices, subsample_size))

    def time_align_series_irregular_string(self):
        (self.x + self.y)


class LibFastZip(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)

        self.df2 = self.df.copy()
        self.df2.ix[:10000, :] = np.nan
        self.col_array_list2 = list(self.df2.values.T)

    def time_lib_fast_zip(self):
        lib.fast_zip(self.col_array_list)

    def time_lib_fast_zip_fillna(self):
        lib.fast_zip_fillna(self.col_array_list2)
