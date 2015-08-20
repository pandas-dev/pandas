from pandas_vb_common import *
from random import shuffle


class dataframe_reindex(object):
    goal_time = 0.2

    def setup(self):
        self.rng = DatetimeIndex(start='1/1/1970', periods=10000, freq=datetools.Minute())
        self.df = DataFrame(np.random.rand(10000, 10), index=self.rng, columns=range(10))
        self.df['foo'] = 'bar'
        self.rng2 = Index(self.rng[::2])

    def time_dataframe_reindex(self):
        self.df.reindex(self.rng2)


class frame_drop_dup_inplace(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)

    def time_frame_drop_dup_inplace(self):
        self.df.drop_duplicates(['key1', 'key2'], inplace=True)


class frame_drop_dup_na_inplace(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)
        self.df.ix[:10000, :] = np.nan

    def time_frame_drop_dup_na_inplace(self):
        self.df.drop_duplicates(['key1', 'key2'], inplace=True)


class frame_drop_duplicates(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)

    def time_frame_drop_duplicates(self):
        self.df.drop_duplicates(['key1', 'key2'])


class frame_drop_duplicates_na(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)
        self.df.ix[:10000, :] = np.nan

    def time_frame_drop_duplicates_na(self):
        self.df.drop_duplicates(['key1', 'key2'])


class frame_fillna_many_columns_pad(object):
    goal_time = 0.2

    def setup(self):
        self.values = np.random.randn(1000, 1000)
        self.values[::2] = np.nan
        self.df = DataFrame(self.values)

    def time_frame_fillna_many_columns_pad(self):
        self.df.fillna(method='pad')


class frame_reindex_columns(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(index=range(10000), data=np.random.rand(10000, 30), columns=range(30))

    def time_frame_reindex_columns(self):
        self.df.reindex(columns=self.df.columns[1:5])


class frame_sort_index_by_columns(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)

    def time_frame_sort_index_by_columns(self):
        self.df.sort_index(by=['key1', 'key2'])


class lib_fast_zip(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)

    def time_lib_fast_zip(self):
        lib.fast_zip(self.col_array_list)


class lib_fast_zip_fillna(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)
        self.df.ix[:10000, :] = np.nan

    def time_lib_fast_zip_fillna(self):
        lib.fast_zip_fillna(self.col_array_list)


class reindex_daterange_backfill(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

        def pad(source_series, target_index):
            try:
                source_series.reindex(target_index, method='pad')
            except:
                source_series.reindex(target_index, fillMethod='pad')

        def backfill(source_series, target_index):
            try:
                source_series.reindex(target_index, method='backfill')
            except:
                source_series.reindex(target_index, fillMethod='backfill')

    def time_reindex_daterange_backfill(self):
        backfill(self.ts2, self.ts.index)


class reindex_daterange_pad(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

        def pad(source_series, target_index):
            try:
                source_series.reindex(target_index, method='pad')
            except:
                source_series.reindex(target_index, fillMethod='pad')

        def backfill(source_series, target_index):
            try:
                source_series.reindex(target_index, method='backfill')
            except:
                source_series.reindex(target_index, fillMethod='backfill')

    def time_reindex_daterange_pad(self):
        pad(self.ts2, self.ts.index)


class reindex_fillna_backfill(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

        def pad(source_series, target_index):
            try:
                source_series.reindex(target_index, method='pad')
            except:
                source_series.reindex(target_index, fillMethod='pad')

        def backfill(source_series, target_index):
            try:
                source_series.reindex(target_index, method='backfill')
            except:
                source_series.reindex(target_index, fillMethod='backfill')

    def time_reindex_fillna_backfill(self):
        self.ts3.fillna(method='backfill')


class reindex_fillna_backfill_float32(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

        def pad(source_series, target_index):
            try:
                source_series.reindex(target_index, method='pad')
            except:
                source_series.reindex(target_index, fillMethod='pad')

        def backfill(source_series, target_index):
            try:
                source_series.reindex(target_index, method='backfill')
            except:
                source_series.reindex(target_index, fillMethod='backfill')

    def time_reindex_fillna_backfill_float32(self):
        self.ts4.fillna(method='backfill')


class reindex_fillna_pad(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

        def pad(source_series, target_index):
            try:
                source_series.reindex(target_index, method='pad')
            except:
                source_series.reindex(target_index, fillMethod='pad')

        def backfill(source_series, target_index):
            try:
                source_series.reindex(target_index, method='backfill')
            except:
                source_series.reindex(target_index, fillMethod='backfill')

    def time_reindex_fillna_pad(self):
        self.ts3.fillna(method='pad')


class reindex_fillna_pad_float32(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)
        self.ts2 = self.ts[::2]
        self.ts3 = self.ts2.reindex(self.ts.index)
        self.ts4 = self.ts3.astype('float32')

        def pad(source_series, target_index):
            try:
                source_series.reindex(target_index, method='pad')
            except:
                source_series.reindex(target_index, fillMethod='pad')

        def backfill(source_series, target_index):
            try:
                source_series.reindex(target_index, method='backfill')
            except:
                source_series.reindex(target_index, fillMethod='backfill')

    def time_reindex_fillna_pad_float32(self):
        self.ts4.fillna(method='pad')


class reindex_frame_level_align(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4), index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4), index=self.index.levels[1])

    def time_reindex_frame_level_align(self):
        self.df.align(self.df_level, level=1, copy=False)


class reindex_frame_level_reindex(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4), index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4), index=self.index.levels[1])

    def time_reindex_frame_level_reindex(self):
        self.df_level.reindex(self.df.index, level=1)


class reindex_multiindex(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000
        self.K = 20
        self.level1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.level2 = np.tile(tm.makeStringIndex(self.K).values, self.N)
        self.index = MultiIndex.from_arrays([self.level1, self.level2])
        self.s1 = Series(np.random.randn((self.N * self.K)), index=self.index)
        self.s2 = self.s1[::2]

    def time_reindex_multiindex(self):
        self.s1.reindex(self.s2.index)


class series_align_irregular_string(object):
    goal_time = 0.2

    def setup(self):
        self.n = 50000
        self.indices = tm.makeStringIndex(self.n)

        def sample(values, k):
            self.sampler = np.arange(len(values))
            shuffle(self.sampler)
            return values.take(self.sampler[:k])
        self.subsample_size = 40000
        self.x = Series(np.random.randn(50000), self.indices)
        self.y = Series(np.random.randn(self.subsample_size), index=sample(self.indices, self.subsample_size))

    def time_series_align_irregular_string(self):
        (self.x + self.y)


class series_drop_duplicates_int(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randint(0, 1000, size=10000))
        self.s2 = Series(np.tile(tm.makeStringIndex(1000).values, 10))

    def time_series_drop_duplicates_int(self):
        self.s.drop_duplicates()


class series_drop_duplicates_string(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randint(0, 1000, size=10000))
        self.s2 = Series(np.tile(tm.makeStringIndex(1000).values, 10))

    def time_series_drop_duplicates_string(self):
        self.s2.drop_duplicates()