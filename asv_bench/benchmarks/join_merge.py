from .pandas_vb_common import *

try:
    from pandas import merge_ordered
except ImportError:
    from pandas import ordered_merge as merge_ordered


#----------------------------------------------------------------------
# Append

class Append(object):
    goal_time = 0.2

    def setup(self):
        self.df1 = pd.DataFrame(np.random.randn(10000, 4),
                                columns=['A', 'B', 'C', 'D'])
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

    def time_append_homogenous(self):
        self.df1.append(self.df2)

    def time_append_mixed(self):
        self.mdf1.append(self.mdf2)


#----------------------------------------------------------------------
# Concat

class Concat(object):
    goal_time = 0.2

    def setup(self):
        self.n = 1000
        self.indices = tm.makeStringIndex(1000)
        self.s = Series(self.n, index=self.indices)
        self.pieces = [self.s[i:(- i)] for i in range(1, 10)]
        self.pieces = (self.pieces * 50)

        self.df_small = pd.DataFrame(randn(5, 4))

        # empty
        self.df = pd.DataFrame(dict(A=range(10000)), index=date_range('20130101', periods=10000, freq='s'))
        self.empty = pd.DataFrame()

    def time_concat_series_axis1(self):
        concat(self.pieces, axis=1)

    def time_concat_small_frames(self):
        concat(([self.df_small] * 1000))

    def time_concat_empty_frames1(self):
        concat([self.df, self.empty])

    def time_concat_empty_frames2(self):
        concat([self.empty, self.df])


class ConcatPanels(object):
    goal_time = 0.2

    def setup(self):
        dataset = np.zeros((10000, 200, 2), dtype=np.float32)
        self.panels_f = [pd.Panel(np.copy(dataset, order='F'))
                         for i in range(20)]
        self.panels_c = [pd.Panel(np.copy(dataset, order='C'))
                         for i in range(20)]

    def time_c_ordered_axis0(self):
        concat(self.panels_c, axis=0, ignore_index=True)

    def time_f_ordered_axis0(self):
        concat(self.panels_f, axis=0, ignore_index=True)

    def time_c_ordered_axis1(self):
        concat(self.panels_c, axis=1, ignore_index=True)

    def time_f_ordered_axis1(self):
        concat(self.panels_f, axis=1, ignore_index=True)

    def time_c_ordered_axis2(self):
        concat(self.panels_c, axis=2, ignore_index=True)

    def time_f_ordered_axis2(self):
        concat(self.panels_f, axis=2, ignore_index=True)


class ConcatFrames(object):
    goal_time = 0.2

    def setup(self):
        dataset = np.zeros((10000, 200), dtype=np.float32)

        self.frames_f = [pd.DataFrame(np.copy(dataset, order='F'))
                         for i in range(20)]
        self.frames_c = [pd.DataFrame(np.copy(dataset, order='C'))
                         for i in range(20)]

    def time_c_ordered_axis0(self):
        concat(self.frames_c, axis=0, ignore_index=True)

    def time_f_ordered_axis0(self):
        concat(self.frames_f, axis=0, ignore_index=True)

    def time_c_ordered_axis1(self):
        concat(self.frames_c, axis=1, ignore_index=True)

    def time_f_ordered_axis1(self):
        concat(self.frames_f, axis=1, ignore_index=True)


#----------------------------------------------------------------------
# Joins

class Join(object):
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
            self.index2 = MultiIndex(levels=[self.level1, self.level2],
                                     labels=[self.label1, self.label2])
            self.index3 = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)],
                                     labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
            self.df_multi = DataFrame(np.random.randn(len(self.index2), 4),
                                      index=self.index2,
                                      columns=['A', 'B', 'C', 'D'])
        except:
            pass
        self.df = pd.DataFrame({'data1': np.random.randn(100000),
                                'data2': np.random.randn(100000),
                                'key1': self.key1,
                                'key2': self.key2})
        self.df_key1 = pd.DataFrame(np.random.randn(len(self.level1), 4),
                                    index=self.level1,
                                    columns=['A', 'B', 'C', 'D'])
        self.df_key2 = pd.DataFrame(np.random.randn(len(self.level2), 4),
                                    index=self.level2,
                                    columns=['A', 'B', 'C', 'D'])
        self.df_shuf = self.df.reindex(self.df.index[self.shuf])

    def time_join_dataframe_index_multi(self):
        self.df.join(self.df_multi, on=['key1', 'key2'])

    def time_join_dataframe_index_single_key_bigger(self):
        self.df.join(self.df_key2, on='key2')

    def time_join_dataframe_index_single_key_bigger_sort(self):
        self.df_shuf.join(self.df_key2, on='key2', sort=True)

    def time_join_dataframe_index_single_key_small(self):
        self.df.join(self.df_key1, on='key1')


class JoinIndex(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(2718281)
        self.n = 50000
        self.left = pd.DataFrame(np.random.randint(1, (self.n / 500), (self.n, 2)), columns=['jim', 'joe'])
        self.right = pd.DataFrame(np.random.randint(1, (self.n / 500), (self.n, 2)), columns=['jolie', 'jolia']).set_index('jolie')

    def time_left_outer_join_index(self):
        self.left.join(self.right, on='jim')


class join_non_unique_equal(object):
    # outer join of non-unique
    # GH 6329

    goal_time = 0.2

    def setup(self):
        self.date_index = date_range('01-Jan-2013', '23-Jan-2013', freq='T')
        self.daily_dates = self.date_index.to_period('D').to_timestamp('S', 'S')
        self.fracofday = (self.date_index.view(np.ndarray) - self.daily_dates.view(np.ndarray))
        self.fracofday = (self.fracofday.astype('timedelta64[ns]').astype(np.float64) / 86400000000000.0)
        self.fracofday = Series(self.fracofday, self.daily_dates)
        self.index = date_range(self.date_index.min().to_period('A').to_timestamp('D', 'S'), self.date_index.max().to_period('A').to_timestamp('D', 'E'), freq='D')
        self.temp = Series(1.0, self.index)

    def time_join_non_unique_equal(self):
        (self.fracofday * self.temp[self.fracofday.index])


#----------------------------------------------------------------------
# Merges

class Merge(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.indices = tm.makeStringIndex(self.N).values
        self.indices2 = tm.makeStringIndex(self.N).values
        self.key = np.tile(self.indices[:8000], 10)
        self.key2 = np.tile(self.indices2[:8000], 10)
        self.left = pd.DataFrame({'key': self.key, 'key2': self.key2,
                                  'value': np.random.randn(80000)})
        self.right = pd.DataFrame({'key': self.indices[2000:],
                                   'key2': self.indices2[2000:],
                                   'value2': np.random.randn(8000)})

        self.df = pd.DataFrame({'key1': np.tile(np.arange(500).repeat(10), 2),
                                'key2': np.tile(np.arange(250).repeat(10), 4),
                                'value': np.random.randn(10000)})
        self.df2 = pd.DataFrame({'key1': np.arange(500), 'value2': randn(500)})
        self.df3 = self.df[:5000]

    def time_merge_2intkey_nosort(self):
        merge(self.left, self.right, sort=False)

    def time_merge_2intkey_sort(self):
        merge(self.left, self.right, sort=True)

    def time_merge_dataframe_integer_2key(self):
        merge(self.df, self.df3)

    def time_merge_dataframe_integer_key(self):
        merge(self.df, self.df2, on='key1')


class i8merge(object):
    goal_time = 0.2

    def setup(self):
        (low, high, n) = (((-1) << 10), (1 << 10), (1 << 20))
        self.left = pd.DataFrame(np.random.randint(low, high, (n, 7)),
                                 columns=list('ABCDEFG'))
        self.left['left'] = self.left.sum(axis=1)
        self.i = np.random.permutation(len(self.left))
        self.right = self.left.iloc[self.i].copy()
        self.right.columns = (self.right.columns[:(-1)].tolist() + ['right'])
        self.right.index = np.arange(len(self.right))
        self.right['right'] *= (-1)

    def time_i8merge(self):
        merge(self.left, self.right, how='outer')


#----------------------------------------------------------------------
# Ordered merge

class MergeOrdered(object):

    def setup(self):

        groups = tm.makeStringIndex(10).values

        self.left = pd.DataFrame({'group': groups.repeat(5000),
                                  'key' : np.tile(np.arange(0, 10000, 2), 10),
                                  'lvalue': np.random.randn(50000)})

        self.right = pd.DataFrame({'key' : np.arange(10000),
                                   'rvalue' : np.random.randn(10000)})

    def time_merge_ordered(self):
        merge_ordered(self.left, self.right, on='key', left_by='group')


# ----------------------------------------------------------------------
# asof merge

class MergeAsof(object):

    def setup(self):
        import string
        np.random.seed(0)
        one_count = 200000
        two_count = 1000000

        self.df1 = pd.DataFrame(
            {'time': np.random.randint(0, one_count / 20, one_count),
             'key': np.random.choice(list(string.uppercase), one_count),
             'key2': np.random.randint(0, 25, one_count),
             'value1': np.random.randn(one_count)})
        self.df2 = pd.DataFrame(
            {'time': np.random.randint(0, two_count / 20, two_count),
             'key': np.random.choice(list(string.uppercase), two_count),
             'key2': np.random.randint(0, 25, two_count),
             'value2': np.random.randn(two_count)})

        self.df1 = self.df1.sort_values('time')
        self.df2 = self.df2.sort_values('time')

        self.df1['time32'] = np.int32(self.df1.time)
        self.df2['time32'] = np.int32(self.df2.time)

        self.df1a = self.df1[['time', 'value1']]
        self.df2a = self.df2[['time', 'value2']]
        self.df1b = self.df1[['time', 'key', 'value1']]
        self.df2b = self.df2[['time', 'key', 'value2']]
        self.df1c = self.df1[['time', 'key2', 'value1']]
        self.df2c = self.df2[['time', 'key2', 'value2']]
        self.df1d = self.df1[['time32', 'value1']]
        self.df2d = self.df2[['time32', 'value2']]
        self.df1e = self.df1[['time', 'key', 'key2', 'value1']]
        self.df2e = self.df2[['time', 'key', 'key2', 'value2']]

    def time_noby(self):
        merge_asof(self.df1a, self.df2a, on='time')

    def time_by_object(self):
        merge_asof(self.df1b, self.df2b, on='time', by='key')

    def time_by_int(self):
        merge_asof(self.df1c, self.df2c, on='time', by='key2')

    def time_on_int32(self):
        merge_asof(self.df1d, self.df2d, on='time32')

    def time_multiby(self):
        merge_asof(self.df1e, self.df2e, on='time', by=['key', 'key2'])


#----------------------------------------------------------------------
# data alignment

class Align(object):
    goal_time = 0.2

    def setup(self):
        self.n = 1000000
        self.sz = 500000
        self.rng = np.arange(0, 10000000000000, 10000000)
        self.stamps = (np.datetime64(datetime.now()).view('i8') + self.rng)
        self.idx1 = np.sort(self.sample(self.stamps, self.sz))
        self.idx2 = np.sort(self.sample(self.stamps, self.sz))
        self.ts1 = Series(np.random.randn(self.sz), self.idx1)
        self.ts2 = Series(np.random.randn(self.sz), self.idx2)

    def sample(self, values, k):
        self.sampler = np.random.permutation(len(values))
        return values.take(self.sampler[:k])

    def time_series_align_int64_index(self):
        (self.ts1 + self.ts2)

    def time_series_align_left_monotonic(self):
        self.ts1.align(self.ts2, join='left')
