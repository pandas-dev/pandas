from pandas_vb_common import *


class append_frame_single_homogenous(object):
    goal_time = 0.2

    def setup(self):
        self.df1 = pd.DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
        self.df2 = self.df1.copy()
        self.df2.index = np.arange(10000, 20000)
        self.mdf1 = self.df1.copy()
        self.mdf1['obj1'] = 'bar'
        self.mdf1['obj2'] = 'bar'
        self.mdf1['int1'] = 5
        try:
            self.mdf1.consolidate(inplace=True)
        except:
            pass
        self.mdf2 = self.mdf1.copy()
        self.mdf2.index = self.df2.index

    def time_append_frame_single_homogenous(self):
        self.df1.append(self.df2)


class append_frame_single_mixed(object):
    goal_time = 0.2

    def setup(self):
        self.df1 = pd.DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
        self.df2 = self.df1.copy()
        self.df2.index = np.arange(10000, 20000)
        self.mdf1 = self.df1.copy()
        self.mdf1['obj1'] = 'bar'
        self.mdf1['obj2'] = 'bar'
        self.mdf1['int1'] = 5
        try:
            self.mdf1.consolidate(inplace=True)
        except:
            pass
        self.mdf2 = self.mdf1.copy()
        self.mdf2.index = self.df2.index

    def time_append_frame_single_mixed(self):
        self.mdf1.append(self.mdf2)


class concat_empty_frames1(object):
    goal_time = 0.2

    def setup(self):
        self.df = pd.DataFrame(dict(A=range(10000)), index=date_range('20130101', periods=10000, freq='s'))
        self.empty = pd.DataFrame()

    def time_concat_empty_frames1(self):
        concat([self.df, self.empty])


class concat_empty_frames2(object):
    goal_time = 0.2

    def setup(self):
        self.df = pd.DataFrame(dict(A=range(10000)), index=date_range('20130101', periods=10000, freq='s'))
        self.empty = pd.DataFrame()

    def time_concat_empty_frames2(self):
        concat([self.empty, self.df])


class concat_series_axis1(object):
    goal_time = 0.2

    def setup(self):
        self.n = 1000
        self.indices = tm.makeStringIndex(1000)
        self.s = Series(self.n, index=self.indices)
        self.pieces = [self.s[i:(- i)] for i in range(1, 10)]
        self.pieces = (self.pieces * 50)

    def time_concat_series_axis1(self):
        concat(self.pieces, axis=1)


class concat_small_frames(object):
    goal_time = 0.2

    def setup(self):
        self.df = pd.DataFrame(randn(5, 4))

    def time_concat_small_frames(self):
        concat(([self.df] * 1000))


class i8merge(object):
    goal_time = 0.2

    def setup(self):
        (low, high, n) = (((-1) << 10), (1 << 10), (1 << 20))
        self.left = pd.DataFrame(np.random.randint(low, high, (n, 7)), columns=list('ABCDEFG'))
        self.left['left'] = self.left.sum(axis=1)
        self.i = np.random.permutation(len(self.left))
        self.right = self.left.iloc[self.i].copy()
        self.right.columns = (self.right.columns[:(-1)].tolist() + ['right'])
        self.right.index = np.arange(len(self.right))
        self.right['right'] *= (-1)

    def time_i8merge(self):
        merge(self.left, self.right, how='outer')


class join_dataframe_index_multi(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = tm.makeStringIndex(10).values
        self.level2 = tm.makeStringIndex(1000).values
        self.label1 = np.arange(10).repeat(1000)
        self.label2 = np.tile(np.arange(1000), 10)
        self.key1 = np.tile(self.level1.take(self.label1), 10)
        self.key2 = np.tile(self.level2.take(self.label2), 10)
        self.shuf = np.arange(100000)
        random.shuffle(self.shuf)
        try:
            self.index2 = MultiIndex(levels=[self.level1, self.level2], labels=[self.label1, self.label2])
            self.index3 = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
            self.df_multi = DataFrame(np.random.randn(len(self.index2), 4), index=self.index2, columns=['A', 'B', 'C', 'D'])
        except:
            pass
        try:
            self.DataFrame = DataMatrix
        except:
            pass
        self.df = pd.DataFrame({'data1': np.random.randn(100000), 'data2': np.random.randn(100000), 'key1': self.key1, 'key2': self.key2, })
        self.df_key1 = pd.DataFrame(np.random.randn(len(self.level1), 4), index=self.level1, columns=['A', 'B', 'C', 'D'])
        self.df_key2 = pd.DataFrame(np.random.randn(len(self.level2), 4), index=self.level2, columns=['A', 'B', 'C', 'D'])
        self.df_shuf = self.df.reindex(self.df.index[self.shuf])

    def time_join_dataframe_index_multi(self):
        self.df.join(self.df_multi, on=['key1', 'key2'])


class join_dataframe_index_single_key_bigger(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = tm.makeStringIndex(10).values
        self.level2 = tm.makeStringIndex(1000).values
        self.label1 = np.arange(10).repeat(1000)
        self.label2 = np.tile(np.arange(1000), 10)
        self.key1 = np.tile(self.level1.take(self.label1), 10)
        self.key2 = np.tile(self.level2.take(self.label2), 10)
        self.shuf = np.arange(100000)
        random.shuffle(self.shuf)
        try:
            self.index2 = MultiIndex(levels=[self.level1, self.level2], labels=[self.label1, self.label2])
            self.index3 = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
            self.df_multi = DataFrame(np.random.randn(len(self.index2), 4), index=self.index2, columns=['A', 'B', 'C', 'D'])
        except:
            pass
        try:
            self.DataFrame = DataMatrix
        except:
            pass
        self.df = pd.DataFrame({'data1': np.random.randn(100000), 'data2': np.random.randn(100000), 'key1': self.key1, 'key2': self.key2, })
        self.df_key1 = pd.DataFrame(np.random.randn(len(self.level1), 4), index=self.level1, columns=['A', 'B', 'C', 'D'])
        self.df_key2 = pd.DataFrame(np.random.randn(len(self.level2), 4), index=self.level2, columns=['A', 'B', 'C', 'D'])
        self.df_shuf = self.df.reindex(self.df.index[self.shuf])

    def time_join_dataframe_index_single_key_bigger(self):
        self.df.join(self.df_key2, on='key2')


class join_dataframe_index_single_key_bigger_sort(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = tm.makeStringIndex(10).values
        self.level2 = tm.makeStringIndex(1000).values
        self.label1 = np.arange(10).repeat(1000)
        self.label2 = np.tile(np.arange(1000), 10)
        self.key1 = np.tile(self.level1.take(self.label1), 10)
        self.key2 = np.tile(self.level2.take(self.label2), 10)
        self.shuf = np.arange(100000)
        random.shuffle(self.shuf)
        try:
            self.index2 = MultiIndex(levels=[self.level1, self.level2], labels=[self.label1, self.label2])
            self.index3 = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
            self.df_multi = DataFrame(np.random.randn(len(self.index2), 4), index=self.index2, columns=['A', 'B', 'C', 'D'])
        except:
            pass
        try:
            self.DataFrame = DataMatrix
        except:
            pass
        self.df = pd.DataFrame({'data1': np.random.randn(100000), 'data2': np.random.randn(100000), 'key1': self.key1, 'key2': self.key2, })
        self.df_key1 = pd.DataFrame(np.random.randn(len(self.level1), 4), index=self.level1, columns=['A', 'B', 'C', 'D'])
        self.df_key2 = pd.DataFrame(np.random.randn(len(self.level2), 4), index=self.level2, columns=['A', 'B', 'C', 'D'])
        self.df_shuf = self.df.reindex(self.df.index[self.shuf])

    def time_join_dataframe_index_single_key_bigger_sort(self):
        self.df_shuf.join(self.df_key2, on='key2', sort=True)


class join_dataframe_index_single_key_small(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = tm.makeStringIndex(10).values
        self.level2 = tm.makeStringIndex(1000).values
        self.label1 = np.arange(10).repeat(1000)
        self.label2 = np.tile(np.arange(1000), 10)
        self.key1 = np.tile(self.level1.take(self.label1), 10)
        self.key2 = np.tile(self.level2.take(self.label2), 10)
        self.shuf = np.arange(100000)
        random.shuffle(self.shuf)
        try:
            self.index2 = MultiIndex(levels=[self.level1, self.level2], labels=[self.label1, self.label2])
            self.index3 = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
            self.df_multi = DataFrame(np.random.randn(len(self.index2), 4), index=self.index2, columns=['A', 'B', 'C', 'D'])
        except:
            pass
        try:
            self.DataFrame = DataMatrix
        except:
            pass
        self.df = pd.DataFrame({'data1': np.random.randn(100000), 'data2': np.random.randn(100000), 'key1': self.key1, 'key2': self.key2, })
        self.df_key1 = pd.DataFrame(np.random.randn(len(self.level1), 4), index=self.level1, columns=['A', 'B', 'C', 'D'])
        self.df_key2 = pd.DataFrame(np.random.randn(len(self.level2), 4), index=self.level2, columns=['A', 'B', 'C', 'D'])
        self.df_shuf = self.df.reindex(self.df.index[self.shuf])

    def time_join_dataframe_index_single_key_small(self):
        self.df.join(self.df_key1, on='key1')


class join_dataframe_integer_2key(object):
    goal_time = 0.2

    def setup(self):
        self.df = pd.DataFrame({'key1': np.tile(np.arange(500).repeat(10), 2), 'key2': np.tile(np.arange(250).repeat(10), 4), 'value': np.random.randn(10000), })
        self.df2 = pd.DataFrame({'key1': np.arange(500), 'value2': randn(500), })
        self.df3 = self.df[:5000]

    def time_join_dataframe_integer_2key(self):
        merge(self.df, self.df3)


class join_dataframe_integer_key(object):
    goal_time = 0.2

    def setup(self):
        self.df = pd.DataFrame({'key1': np.tile(np.arange(500).repeat(10), 2), 'key2': np.tile(np.arange(250).repeat(10), 4), 'value': np.random.randn(10000), })
        self.df2 = pd.DataFrame({'key1': np.arange(500), 'value2': randn(500), })
        self.df3 = self.df[:5000]

    def time_join_dataframe_integer_key(self):
        merge(self.df, self.df2, on='key1')


class join_non_unique_equal(object):
    goal_time = 0.2

    def setup(self):
        self.date_index = date_range('01-Jan-2013', '23-Jan-2013', freq='T')
        self.daily_dates = self.date_index.to_period('D').to_timestamp('S', 'S')
        self.fracofday = (self.date_index.view(np.ndarray) - self.daily_dates.view(np.ndarray))
        self.fracofday = (self.fracofday.astype('timedelta64[ns]').astype(np.float64) / 86400000000000.0)
        self.fracofday = TimeSeries(self.fracofday, self.daily_dates)
        self.index = date_range(self.date_index.min().to_period('A').to_timestamp('D', 'S'), self.date_index.max().to_period('A').to_timestamp('D', 'E'), freq='D')
        self.temp = TimeSeries(1.0, self.index)

    def time_join_non_unique_equal(self):
        (self.fracofday * self.temp[self.fracofday.index])


class left_outer_join_index(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(2718281)
        self.n = 50000
        self.left = pd.DataFrame(np.random.randint(1, (self.n / 500), (self.n, 2)), columns=['jim', 'joe'])
        self.right = pd.DataFrame(np.random.randint(1, (self.n / 500), (self.n, 2)), columns=['jolie', 'jolia']).set_index('jolie')

    def time_left_outer_join_index(self):
        self.left.join(self.right, on='jim')


class merge_2intkey_nosort(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.indices = tm.makeStringIndex(self.N).values
        self.indices2 = tm.makeStringIndex(self.N).values
        self.key = np.tile(self.indices[:8000], 10)
        self.key2 = np.tile(self.indices2[:8000], 10)
        self.left = pd.DataFrame({'key': self.key, 'key2': self.key2, 'value': np.random.randn(80000), })
        self.right = pd.DataFrame({'key': self.indices[2000:], 'key2': self.indices2[2000:], 'value2': np.random.randn(8000), })

    def time_merge_2intkey_nosort(self):
        merge(self.left, self.right, sort=False)


class merge_2intkey_sort(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.indices = tm.makeStringIndex(self.N).values
        self.indices2 = tm.makeStringIndex(self.N).values
        self.key = np.tile(self.indices[:8000], 10)
        self.key2 = np.tile(self.indices2[:8000], 10)
        self.left = pd.DataFrame({'key': self.key, 'key2': self.key2, 'value': np.random.randn(80000), })
        self.right = pd.DataFrame({'key': self.indices[2000:], 'key2': self.indices2[2000:], 'value2': np.random.randn(8000), })

    def time_merge_2intkey_sort(self):
        merge(self.left, self.right, sort=True)


class series_align_int64_index(object):
    goal_time = 0.2

    def setup(self):
        self.n = 1000000

        def sample(values, k):
            self.sampler = np.random.permutation(len(values))
            return values.take(self.sampler[:k])
        self.sz = 500000
        self.rng = np.arange(0, 10000000000000, 10000000)
        self.stamps = (np.datetime64(datetime.now()).view('i8') + self.rng)
        self.idx1 = np.sort(sample(self.stamps, self.sz))
        self.idx2 = np.sort(sample(self.stamps, self.sz))
        self.ts1 = Series(np.random.randn(self.sz), self.idx1)
        self.ts2 = Series(np.random.randn(self.sz), self.idx2)

    def time_series_align_int64_index(self):
        (self.ts1 + self.ts2)


class series_align_left_monotonic(object):
    goal_time = 0.2

    def setup(self):
        self.n = 1000000

        def sample(values, k):
            self.sampler = np.random.permutation(len(values))
            return values.take(self.sampler[:k])
        self.sz = 500000
        self.rng = np.arange(0, 10000000000000, 10000000)
        self.stamps = (np.datetime64(datetime.now()).view('i8') + self.rng)
        self.idx1 = np.sort(sample(self.stamps, self.sz))
        self.idx2 = np.sort(sample(self.stamps, self.sz))
        self.ts1 = Series(np.random.randn(self.sz), self.idx1)
        self.ts2 = Series(np.random.randn(self.sz), self.idx2)

    def time_series_align_left_monotonic(self):
        self.ts1.align(self.ts2, join='left')