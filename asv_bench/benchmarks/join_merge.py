import warnings
import string

import numpy as np
import pandas.util.testing as tm
from pandas import (DataFrame, Series, Panel, MultiIndex,
                    date_range, concat, merge, merge_asof)

try:
    from pandas import merge_ordered
except ImportError:
    from pandas import ordered_merge as merge_ordered


class Append(object):

    def setup(self):
        self.df1 = DataFrame(np.random.randn(10000, 4),
                             columns=['A', 'B', 'C', 'D'])
        self.df2 = self.df1.copy()
        self.df2.index = np.arange(10000, 20000)
        self.mdf1 = self.df1.copy()
        self.mdf1['obj1'] = 'bar'
        self.mdf1['obj2'] = 'bar'
        self.mdf1['int1'] = 5
        try:
            with warnings.catch_warnings(record=True):
                self.mdf1.consolidate(inplace=True)
        except (AttributeError, TypeError):
            pass
        self.mdf2 = self.mdf1.copy()
        self.mdf2.index = self.df2.index

    def time_append_homogenous(self):
        self.df1.append(self.df2)

    def time_append_mixed(self):
        self.mdf1.append(self.mdf2)


class Concat(object):

    params = [0, 1]
    param_names = ['axis']

    def setup(self, axis):
        N = 1000
        s = Series(N, index=tm.makeStringIndex(N))
        self.series = [s[i:- i] for i in range(1, 10)] * 50
        self.small_frames = [DataFrame(np.random.randn(5, 4))] * 1000
        df = DataFrame({'A': range(N)},
                       index=date_range('20130101', periods=N, freq='s'))
        self.empty_left = [DataFrame(), df]
        self.empty_right = [df, DataFrame()]

    def time_concat_series(self, axis):
        concat(self.series, axis=axis)

    def time_concat_small_frames(self, axis):
        concat(self.small_frames, axis=axis)

    def time_concat_empty_right(self, axis):
        concat(self.empty_right, axis=axis)

    def time_concat_empty_left(self, axis):
        concat(self.empty_left, axis=axis)


class ConcatPanels(object):

    params = ([0, 1, 2], [True, False])
    param_names = ['axis', 'ignore_index']

    def setup(self, axis, ignore_index):
        with warnings.catch_warnings(record=True):
            panel_c = Panel(np.zeros((10000, 200, 2),
                                     dtype=np.float32,
                                     order='C'))
            self.panels_c = [panel_c] * 20
            panel_f = Panel(np.zeros((10000, 200, 2),
                            dtype=np.float32,
                            order='F'))
            self.panels_f = [panel_f] * 20

    def time_c_ordered(self, axis, ignore_index):
        with warnings.catch_warnings(record=True):
            concat(self.panels_c, axis=axis, ignore_index=ignore_index)

    def time_f_ordered(self, axis, ignore_index):
        with warnings.catch_warnings(record=True):
            concat(self.panels_f, axis=axis, ignore_index=ignore_index)


class ConcatDataFrames(object):

    params = ([0, 1], [True, False])
    param_names = ['axis', 'ignore_index']

    def setup(self, axis, ignore_index):
        frame_c = DataFrame(np.zeros((10000, 200),
                            dtype=np.float32, order='C'))
        self.frame_c = [frame_c] * 20
        frame_f = DataFrame(np.zeros((10000, 200),
                            dtype=np.float32, order='F'))
        self.frame_f = [frame_f] * 20

    def time_c_ordered(self, axis, ignore_index):
        concat(self.frame_c, axis=axis, ignore_index=ignore_index)

    def time_f_ordered(self, axis, ignore_index):
        concat(self.frame_f, axis=axis, ignore_index=ignore_index)


class Join(object):

    params = [True, False]
    param_names = ['sort']

    def setup(self, sort):
        level1 = tm.makeStringIndex(10).values
        level2 = tm.makeStringIndex(1000).values
        label1 = np.arange(10).repeat(1000)
        label2 = np.tile(np.arange(1000), 10)
        index2 = MultiIndex(levels=[level1, level2],
                            labels=[label1, label2])
        self.df_multi = DataFrame(np.random.randn(len(index2), 4),
                                  index=index2,
                                  columns=['A', 'B', 'C', 'D'])

        self.key1 = np.tile(level1.take(label1), 10)
        self.key2 = np.tile(level2.take(label2), 10)
        self.df = DataFrame({'data1': np.random.randn(100000),
                             'data2': np.random.randn(100000),
                             'key1': self.key1,
                             'key2': self.key2})

        self.df_key1 = DataFrame(np.random.randn(len(level1), 4),
                                 index=level1,
                                 columns=['A', 'B', 'C', 'D'])
        self.df_key2 = DataFrame(np.random.randn(len(level2), 4),
                                 index=level2,
                                 columns=['A', 'B', 'C', 'D'])

        shuf = np.arange(100000)
        np.random.shuffle(shuf)
        self.df_shuf = self.df.reindex(self.df.index[shuf])

    def time_join_dataframe_index_multi(self, sort):
        self.df.join(self.df_multi, on=['key1', 'key2'], sort=sort)

    def time_join_dataframe_index_single_key_bigger(self, sort):
        self.df.join(self.df_key2, on='key2', sort=sort)

    def time_join_dataframe_index_single_key_small(self, sort):
        self.df.join(self.df_key1, on='key1', sort=sort)

    def time_join_dataframe_index_shuffle_key_bigger_sort(self, sort):
        self.df_shuf.join(self.df_key2, on='key2', sort=sort)


class JoinIndex(object):

    def setup(self):
        N = 50000
        self.left = DataFrame(np.random.randint(1, N / 500, (N, 2)),
                              columns=['jim', 'joe'])
        self.right = DataFrame(np.random.randint(1, N / 500, (N, 2)),
                               columns=['jolie', 'jolia']).set_index('jolie')

    def time_left_outer_join_index(self):
        self.left.join(self.right, on='jim')


class JoinNonUnique(object):
    # outer join of non-unique
    # GH 6329
    def setup(self):
        date_index = date_range('01-Jan-2013', '23-Jan-2013', freq='T')
        daily_dates = date_index.to_period('D').to_timestamp('S', 'S')
        self.fracofday = date_index.values - daily_dates.values
        self.fracofday = self.fracofday.astype('timedelta64[ns]')
        self.fracofday = self.fracofday.astype(np.float64) / 86400000000000.0
        self.fracofday = Series(self.fracofday, daily_dates)
        index = date_range(date_index.min(), date_index.max(), freq='D')
        self.temp = Series(1.0, index)[self.fracofday.index]

    def time_join_non_unique_equal(self):
        self.fracofday * self.temp


class Merge(object):

    params = [True, False]
    param_names = ['sort']

    def setup(self, sort):
        N = 10000
        indices = tm.makeStringIndex(N).values
        indices2 = tm.makeStringIndex(N).values
        key = np.tile(indices[:8000], 10)
        key2 = np.tile(indices2[:8000], 10)
        self.left = DataFrame({'key': key, 'key2': key2,
                               'value': np.random.randn(80000)})
        self.right = DataFrame({'key': indices[2000:],
                                'key2': indices2[2000:],
                                'value2': np.random.randn(8000)})

        self.df = DataFrame({'key1': np.tile(np.arange(500).repeat(10), 2),
                             'key2': np.tile(np.arange(250).repeat(10), 4),
                             'value': np.random.randn(10000)})
        self.df2 = DataFrame({'key1': np.arange(500),
                              'value2': np.random.randn(500)})
        self.df3 = self.df[:5000]

    def time_merge_2intkey(self, sort):
        merge(self.left, self.right, sort=sort)

    def time_merge_dataframe_integer_2key(self, sort):
        merge(self.df, self.df3, sort=sort)

    def time_merge_dataframe_integer_key(self, sort):
        merge(self.df, self.df2, on='key1', sort=sort)


class I8Merge(object):

    params = ['inner', 'outer', 'left', 'right']
    param_names = ['how']

    def setup(self, how):
        low, high, n = -1000, 1000, 10**6
        self.left = DataFrame(np.random.randint(low, high, (n, 7)),
                              columns=list('ABCDEFG'))
        self.left['left'] = self.left.sum(axis=1)
        self.right = self.left.sample(frac=1).rename({'left': 'right'}, axis=1)
        self.right = self.right.reset_index(drop=True)
        self.right['right'] *= -1

    def time_i8merge(self, how):
        merge(self.left, self.right, how=how)


class MergeCategoricals(object):

    def setup(self):
        self.left_object = DataFrame(
            {'X': np.random.choice(range(0, 10), size=(10000,)),
             'Y': np.random.choice(['one', 'two', 'three'], size=(10000,))})

        self.right_object = DataFrame(
            {'X': np.random.choice(range(0, 10), size=(10000,)),
             'Z': np.random.choice(['jjj', 'kkk', 'sss'], size=(10000,))})

        self.left_cat = self.left_object.assign(
            Y=self.left_object['Y'].astype('category'))
        self.right_cat = self.right_object.assign(
            Z=self.right_object['Z'].astype('category'))

    def time_merge_object(self):
        merge(self.left_object, self.right_object, on='X')

    def time_merge_cat(self):
        merge(self.left_cat, self.right_cat, on='X')


class MergeOrdered(object):

    def setup(self):
        groups = tm.makeStringIndex(10).values
        self.left = DataFrame({'group': groups.repeat(5000),
                               'key': np.tile(np.arange(0, 10000, 2), 10),
                               'lvalue': np.random.randn(50000)})
        self.right = DataFrame({'key': np.arange(10000),
                                'rvalue': np.random.randn(10000)})

    def time_merge_ordered(self):
        merge_ordered(self.left, self.right, on='key', left_by='group')


class MergeAsof(object):

    def setup(self):
        one_count = 200000
        two_count = 1000000

        df1 = DataFrame(
            {'time': np.random.randint(0, one_count / 20, one_count),
             'key': np.random.choice(list(string.ascii_uppercase), one_count),
             'key2': np.random.randint(0, 25, one_count),
             'value1': np.random.randn(one_count)})
        df2 = DataFrame(
            {'time': np.random.randint(0, two_count / 20, two_count),
             'key': np.random.choice(list(string.ascii_uppercase), two_count),
             'key2': np.random.randint(0, 25, two_count),
             'value2': np.random.randn(two_count)})

        df1 = df1.sort_values('time')
        df2 = df2.sort_values('time')

        df1['time32'] = np.int32(df1.time)
        df2['time32'] = np.int32(df2.time)

        self.df1a = df1[['time', 'value1']]
        self.df2a = df2[['time', 'value2']]
        self.df1b = df1[['time', 'key', 'value1']]
        self.df2b = df2[['time', 'key', 'value2']]
        self.df1c = df1[['time', 'key2', 'value1']]
        self.df2c = df2[['time', 'key2', 'value2']]
        self.df1d = df1[['time32', 'value1']]
        self.df2d = df2[['time32', 'value2']]
        self.df1e = df1[['time', 'key', 'key2', 'value1']]
        self.df2e = df2[['time', 'key', 'key2', 'value2']]

    def time_on_int(self):
        merge_asof(self.df1a, self.df2a, on='time')

    def time_on_int32(self):
        merge_asof(self.df1d, self.df2d, on='time32')

    def time_by_object(self):
        merge_asof(self.df1b, self.df2b, on='time', by='key')

    def time_by_int(self):
        merge_asof(self.df1c, self.df2c, on='time', by='key2')

    def time_multiby(self):
        merge_asof(self.df1e, self.df2e, on='time', by=['key', 'key2'])


class Align(object):

    def setup(self):
        size = 5 * 10**5
        rng = np.arange(0, 10**13, 10**7)
        stamps = np.datetime64('now').view('i8') + rng
        idx1 = np.sort(np.random.choice(stamps, size, replace=False))
        idx2 = np.sort(np.random.choice(stamps, size, replace=False))
        self.ts1 = Series(np.random.randn(size), idx1)
        self.ts2 = Series(np.random.randn(size), idx2)

    def time_series_align_int64_index(self):
        self.ts1 + self.ts2

    def time_series_align_left_monotonic(self):
        self.ts1.align(self.ts2, join='left')


from .pandas_vb_common import setup  # noqa: F401
