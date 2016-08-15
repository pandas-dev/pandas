import nose

import numpy as np
from numpy.random import randn

from datetime import datetime
from pandas.compat import StringIO
import pandas as pd
from pandas import (DataFrame, concat,
                    read_csv, isnull, Series, date_range,
                    Index, Panel, MultiIndex, Timestamp,
                    DatetimeIndex, Categorical)
from pandas.types.concat import union_categoricals
from pandas.util import testing as tm
from pandas.util.testing import (assert_frame_equal,
                                 makeCustomDataframe as mkdf,
                                 assert_almost_equal)


class ConcatenateBase(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.frame = DataFrame(tm.getSeriesData())
        self.mixed_frame = self.frame.copy()
        self.mixed_frame['foo'] = 'bar'


class TestAppend(ConcatenateBase):

    def test_append(self):
        begin_index = self.frame.index[:5]
        end_index = self.frame.index[5:]

        begin_frame = self.frame.reindex(begin_index)
        end_frame = self.frame.reindex(end_index)

        appended = begin_frame.append(end_frame)
        assert_almost_equal(appended['A'], self.frame['A'])

        del end_frame['A']
        partial_appended = begin_frame.append(end_frame)
        self.assertIn('A', partial_appended)

        partial_appended = end_frame.append(begin_frame)
        self.assertIn('A', partial_appended)

        # mixed type handling
        appended = self.mixed_frame[:5].append(self.mixed_frame[5:])
        assert_frame_equal(appended, self.mixed_frame)

        # what to test here
        mixed_appended = self.mixed_frame[:5].append(self.frame[5:])
        mixed_appended2 = self.frame[:5].append(self.mixed_frame[5:])

        # all equal except 'foo' column
        assert_frame_equal(
            mixed_appended.reindex(columns=['A', 'B', 'C', 'D']),
            mixed_appended2.reindex(columns=['A', 'B', 'C', 'D']))

        # append empty
        empty = DataFrame({})

        appended = self.frame.append(empty)
        assert_frame_equal(self.frame, appended)
        self.assertIsNot(appended, self.frame)

        appended = empty.append(self.frame)
        assert_frame_equal(self.frame, appended)
        self.assertIsNot(appended, self.frame)

        # overlap
        self.assertRaises(ValueError, self.frame.append, self.frame,
                          verify_integrity=True)

        # new columns
        # GH 6129
        df = DataFrame({'a': {'x': 1, 'y': 2}, 'b': {'x': 3, 'y': 4}})
        row = Series([5, 6, 7], index=['a', 'b', 'c'], name='z')
        expected = DataFrame({'a': {'x': 1, 'y': 2, 'z': 5}, 'b': {
                             'x': 3, 'y': 4, 'z': 6}, 'c': {'z': 7}})
        result = df.append(row)
        assert_frame_equal(result, expected)

    def test_append_length0_frame(self):
        df = DataFrame(columns=['A', 'B', 'C'])
        df3 = DataFrame(index=[0, 1], columns=['A', 'B'])
        df5 = df.append(df3)

        expected = DataFrame(index=[0, 1], columns=['A', 'B', 'C'])
        assert_frame_equal(df5, expected)

    def test_append_records(self):
        arr1 = np.zeros((2,), dtype=('i4,f4,a10'))
        arr1[:] = [(1, 2., 'Hello'), (2, 3., "World")]

        arr2 = np.zeros((3,), dtype=('i4,f4,a10'))
        arr2[:] = [(3, 4., 'foo'),
                   (5, 6., "bar"),
                   (7., 8., 'baz')]

        df1 = DataFrame(arr1)
        df2 = DataFrame(arr2)

        result = df1.append(df2, ignore_index=True)
        expected = DataFrame(np.concatenate((arr1, arr2)))
        assert_frame_equal(result, expected)

    def test_append_different_columns(self):
        df = DataFrame({'bools': np.random.randn(10) > 0,
                        'ints': np.random.randint(0, 10, 10),
                        'floats': np.random.randn(10),
                        'strings': ['foo', 'bar'] * 5})

        a = df[:5].ix[:, ['bools', 'ints', 'floats']]
        b = df[5:].ix[:, ['strings', 'ints', 'floats']]

        appended = a.append(b)
        self.assertTrue(isnull(appended['strings'][0:4]).all())
        self.assertTrue(isnull(appended['bools'][5:]).all())

    def test_append_many(self):
        chunks = [self.frame[:5], self.frame[5:10],
                  self.frame[10:15], self.frame[15:]]

        result = chunks[0].append(chunks[1:])
        tm.assert_frame_equal(result, self.frame)

        chunks[-1] = chunks[-1].copy()
        chunks[-1]['foo'] = 'bar'
        result = chunks[0].append(chunks[1:])
        tm.assert_frame_equal(result.ix[:, self.frame.columns], self.frame)
        self.assertTrue((result['foo'][15:] == 'bar').all())
        self.assertTrue(result['foo'][:15].isnull().all())

    def test_append_preserve_index_name(self):
        # #980
        df1 = DataFrame(data=None, columns=['A', 'B', 'C'])
        df1 = df1.set_index(['A'])
        df2 = DataFrame(data=[[1, 4, 7], [2, 5, 8], [3, 6, 9]],
                        columns=['A', 'B', 'C'])
        df2 = df2.set_index(['A'])

        result = df1.append(df2)
        self.assertEqual(result.index.name, 'A')

    def test_append_dtype_coerce(self):

        # GH 4993
        # appending with datetime will incorrectly convert datetime64
        import datetime as dt
        from pandas import NaT

        df1 = DataFrame(index=[1, 2], data=[dt.datetime(2013, 1, 1, 0, 0),
                                            dt.datetime(2013, 1, 2, 0, 0)],
                        columns=['start_time'])
        df2 = DataFrame(index=[4, 5], data=[[dt.datetime(2013, 1, 3, 0, 0),
                                             dt.datetime(2013, 1, 3, 6, 10)],
                                            [dt.datetime(2013, 1, 4, 0, 0),
                                             dt.datetime(2013, 1, 4, 7, 10)]],
                        columns=['start_time', 'end_time'])

        expected = concat([Series([NaT, NaT, dt.datetime(2013, 1, 3, 6, 10),
                                   dt.datetime(2013, 1, 4, 7, 10)],
                                  name='end_time'),
                           Series([dt.datetime(2013, 1, 1, 0, 0),
                                   dt.datetime(2013, 1, 2, 0, 0),
                                   dt.datetime(2013, 1, 3, 0, 0),
                                   dt.datetime(2013, 1, 4, 0, 0)],
                                  name='start_time')], axis=1)
        result = df1.append(df2, ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_missing_column_proper_upcast(self):
        df1 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='i8')})
        df2 = DataFrame({'B': np.array([True, False, True, False],
                                       dtype=bool)})

        appended = df1.append(df2, ignore_index=True)
        self.assertEqual(appended['A'].dtype, 'f8')
        self.assertEqual(appended['B'].dtype, 'O')


class TestConcatenate(ConcatenateBase):

    def test_concat_copy(self):

        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randint(0, 10, size=4).reshape(4, 1))
        df3 = DataFrame({5: 'foo'}, index=range(4))

        # these are actual copies
        result = concat([df, df2, df3], axis=1, copy=True)
        for b in result._data.blocks:
            self.assertIsNone(b.values.base)

        # these are the same
        result = concat([df, df2, df3], axis=1, copy=False)
        for b in result._data.blocks:
            if b.is_float:
                self.assertTrue(
                    b.values.base is df._data.blocks[0].values.base)
            elif b.is_integer:
                self.assertTrue(
                    b.values.base is df2._data.blocks[0].values.base)
            elif b.is_object:
                self.assertIsNotNone(b.values.base)

        # float block was consolidated
        df4 = DataFrame(np.random.randn(4, 1))
        result = concat([df, df2, df3, df4], axis=1, copy=False)
        for b in result._data.blocks:
            if b.is_float:
                self.assertIsNone(b.values.base)
            elif b.is_integer:
                self.assertTrue(
                    b.values.base is df2._data.blocks[0].values.base)
            elif b.is_object:
                self.assertIsNotNone(b.values.base)

    def test_concat_with_group_keys(self):
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randn(4, 4))

        # axis=0
        df = DataFrame(np.random.randn(3, 4))
        df2 = DataFrame(np.random.randn(4, 4))

        result = concat([df, df2], keys=[0, 1])
        exp_index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1, 1],
                                            [0, 1, 2, 0, 1, 2, 3]])
        expected = DataFrame(np.r_[df.values, df2.values],
                             index=exp_index)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df], keys=[0, 1])
        exp_index2 = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1],
                                             [0, 1, 2, 0, 1, 2]])
        expected = DataFrame(np.r_[df.values, df.values],
                             index=exp_index2)
        tm.assert_frame_equal(result, expected)

        # axis=1
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randn(4, 4))

        result = concat([df, df2], keys=[0, 1], axis=1)
        expected = DataFrame(np.c_[df.values, df2.values],
                             columns=exp_index)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df], keys=[0, 1], axis=1)
        expected = DataFrame(np.c_[df.values, df.values],
                             columns=exp_index2)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_specific_levels(self):
        df = DataFrame(np.random.randn(10, 4))
        pieces = [df.ix[:, [0, 1]], df.ix[:, [2]], df.ix[:, [3]]]
        level = ['three', 'two', 'one', 'zero']
        result = concat(pieces, axis=1, keys=['one', 'two', 'three'],
                        levels=[level],
                        names=['group_key'])

        self.assert_index_equal(result.columns.levels[0],
                                Index(level, name='group_key'))
        self.assertEqual(result.columns.names[0], 'group_key')

    def test_concat_dataframe_keys_bug(self):
        t1 = DataFrame({
            'value': Series([1, 2, 3], index=Index(['a', 'b', 'c'],
                                                   name='id'))})
        t2 = DataFrame({
            'value': Series([7, 8], index=Index(['a', 'b'], name='id'))})

        # it works
        result = concat([t1, t2], axis=1, keys=['t1', 't2'])
        self.assertEqual(list(result.columns), [('t1', 'value'),
                                                ('t2', 'value')])

    def test_concat_series_partial_columns_names(self):
        # GH10698
        foo = Series([1, 2], name='foo')
        bar = Series([1, 2])
        baz = Series([4, 5])

        result = concat([foo, bar, baz], axis=1)
        expected = DataFrame({'foo': [1, 2], 0: [1, 2], 1: [
                             4, 5]}, columns=['foo', 0, 1])
        tm.assert_frame_equal(result, expected)

        result = concat([foo, bar, baz], axis=1, keys=[
                        'red', 'blue', 'yellow'])
        expected = DataFrame({'red': [1, 2], 'blue': [1, 2], 'yellow': [
                             4, 5]}, columns=['red', 'blue', 'yellow'])
        tm.assert_frame_equal(result, expected)

        result = concat([foo, bar, baz], axis=1, ignore_index=True)
        expected = DataFrame({0: [1, 2], 1: [1, 2], 2: [4, 5]})
        tm.assert_frame_equal(result, expected)

    def test_concat_dict(self):
        frames = {'foo': DataFrame(np.random.randn(4, 3)),
                  'bar': DataFrame(np.random.randn(4, 3)),
                  'baz': DataFrame(np.random.randn(4, 3)),
                  'qux': DataFrame(np.random.randn(4, 3))}

        sorted_keys = sorted(frames)

        result = concat(frames)
        expected = concat([frames[k] for k in sorted_keys], keys=sorted_keys)
        tm.assert_frame_equal(result, expected)

        result = concat(frames, axis=1)
        expected = concat([frames[k] for k in sorted_keys], keys=sorted_keys,
                          axis=1)
        tm.assert_frame_equal(result, expected)

        keys = ['baz', 'foo', 'bar']
        result = concat(frames, keys=keys)
        expected = concat([frames[k] for k in keys], keys=keys)
        tm.assert_frame_equal(result, expected)

    def test_concat_ignore_index(self):
        frame1 = DataFrame({"test1": ["a", "b", "c"],
                            "test2": [1, 2, 3],
                            "test3": [4.5, 3.2, 1.2]})
        frame2 = DataFrame({"test3": [5.2, 2.2, 4.3]})
        frame1.index = Index(["x", "y", "z"])
        frame2.index = Index(["x", "y", "q"])

        v1 = concat([frame1, frame2], axis=1, ignore_index=True)

        nan = np.nan
        expected = DataFrame([[nan, nan, nan, 4.3],
                              ['a', 1, 4.5, 5.2],
                              ['b', 2, 3.2, 2.2],
                              ['c', 3, 1.2, nan]],
                             index=Index(["q", "x", "y", "z"]))

        tm.assert_frame_equal(v1, expected)

    def test_concat_multiindex_with_keys(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        frame = DataFrame(np.random.randn(10, 3), index=index,
                          columns=Index(['A', 'B', 'C'], name='exp'))
        result = concat([frame, frame], keys=[0, 1], names=['iteration'])

        self.assertEqual(result.index.names, ('iteration',) + index.names)
        tm.assert_frame_equal(result.ix[0], frame)
        tm.assert_frame_equal(result.ix[1], frame)
        self.assertEqual(result.index.nlevels, 3)

    def test_concat_multiindex_with_tz(self):
        # GH 6606
        df = DataFrame({'dt': [datetime(2014, 1, 1),
                               datetime(2014, 1, 2),
                               datetime(2014, 1, 3)],
                        'b': ['A', 'B', 'C'],
                        'c': [1, 2, 3], 'd': [4, 5, 6]})
        df['dt'] = df['dt'].apply(lambda d: Timestamp(d, tz='US/Pacific'))
        df = df.set_index(['dt', 'b'])

        exp_idx1 = DatetimeIndex(['2014-01-01', '2014-01-02',
                                  '2014-01-03'] * 2,
                                 tz='US/Pacific', name='dt')
        exp_idx2 = Index(['A', 'B', 'C'] * 2, name='b')
        exp_idx = MultiIndex.from_arrays([exp_idx1, exp_idx2])
        expected = DataFrame({'c': [1, 2, 3] * 2, 'd': [4, 5, 6] * 2},
                             index=exp_idx, columns=['c', 'd'])

        result = concat([df, df])
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_and_levels(self):
        df = DataFrame(np.random.randn(1, 3))
        df2 = DataFrame(np.random.randn(1, 4))

        levels = [['foo', 'baz'], ['one', 'two']]
        names = ['first', 'second']
        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        levels=levels,
                        names=names)
        expected = concat([df, df2, df, df2])
        exp_index = MultiIndex(levels=levels + [[0]],
                               labels=[[0, 0, 1, 1], [0, 1, 0, 1],
                                       [0, 0, 0, 0]],
                               names=names + [None])
        expected.index = exp_index

        assert_frame_equal(result, expected)

        # no names

        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        levels=levels)
        self.assertEqual(result.index.names, (None,) * 3)

        # no levels
        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        names=['first', 'second'])
        self.assertEqual(result.index.names, ('first', 'second') + (None,))
        self.assert_index_equal(result.index.levels[0],
                                Index(['baz', 'foo'], name='first'))

    def test_concat_keys_levels_no_overlap(self):
        # GH #1406
        df = DataFrame(np.random.randn(1, 3), index=['a'])
        df2 = DataFrame(np.random.randn(1, 4), index=['b'])

        self.assertRaises(ValueError, concat, [df, df],
                          keys=['one', 'two'], levels=[['foo', 'bar', 'baz']])

        self.assertRaises(ValueError, concat, [df, df2],
                          keys=['one', 'two'], levels=[['foo', 'bar', 'baz']])

    def test_concat_rename_index(self):
        a = DataFrame(np.random.rand(3, 3),
                      columns=list('ABC'),
                      index=Index(list('abc'), name='index_a'))
        b = DataFrame(np.random.rand(3, 3),
                      columns=list('ABC'),
                      index=Index(list('abc'), name='index_b'))

        result = concat([a, b], keys=['key0', 'key1'],
                        names=['lvl0', 'lvl1'])

        exp = concat([a, b], keys=['key0', 'key1'], names=['lvl0'])
        names = list(exp.index.names)
        names[1] = 'lvl1'
        exp.index.set_names(names, inplace=True)

        tm.assert_frame_equal(result, exp)
        self.assertEqual(result.index.names, exp.index.names)

    def test_crossed_dtypes_weird_corner(self):
        columns = ['A', 'B', 'C', 'D']
        df1 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='f8'),
                         'B': np.array([1, 2, 3, 4], dtype='i8'),
                         'C': np.array([1, 2, 3, 4], dtype='f8'),
                         'D': np.array([1, 2, 3, 4], dtype='i8')},
                        columns=columns)

        df2 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='i8'),
                         'B': np.array([1, 2, 3, 4], dtype='f8'),
                         'C': np.array([1, 2, 3, 4], dtype='i8'),
                         'D': np.array([1, 2, 3, 4], dtype='f8')},
                        columns=columns)

        appended = df1.append(df2, ignore_index=True)
        expected = DataFrame(np.concatenate([df1.values, df2.values], axis=0),
                             columns=columns)
        tm.assert_frame_equal(appended, expected)

        df = DataFrame(np.random.randn(1, 3), index=['a'])
        df2 = DataFrame(np.random.randn(1, 4), index=['b'])
        result = concat(
            [df, df2], keys=['one', 'two'], names=['first', 'second'])
        self.assertEqual(result.index.names, ('first', 'second'))

    def test_dups_index(self):
        # GH 4771

        # single dtypes
        df = DataFrame(np.random.randint(0, 10, size=40).reshape(
            10, 4), columns=['A', 'A', 'C', 'C'])

        result = concat([df, df], axis=1)
        assert_frame_equal(result.iloc[:, :4], df)
        assert_frame_equal(result.iloc[:, 4:], df)

        result = concat([df, df], axis=0)
        assert_frame_equal(result.iloc[:10], df)
        assert_frame_equal(result.iloc[10:], df)

        # multi dtypes
        df = concat([DataFrame(np.random.randn(10, 4),
                               columns=['A', 'A', 'B', 'B']),
                     DataFrame(np.random.randint(0, 10, size=20)
                               .reshape(10, 2),
                               columns=['A', 'C'])],
                    axis=1)

        result = concat([df, df], axis=1)
        assert_frame_equal(result.iloc[:, :6], df)
        assert_frame_equal(result.iloc[:, 6:], df)

        result = concat([df, df], axis=0)
        assert_frame_equal(result.iloc[:10], df)
        assert_frame_equal(result.iloc[10:], df)

        # append
        result = df.iloc[0:8, :].append(df.iloc[8:])
        assert_frame_equal(result, df)

        result = df.iloc[0:8, :].append(df.iloc[8:9]).append(df.iloc[9:10])
        assert_frame_equal(result, df)

        expected = concat([df, df], axis=0)
        result = df.append(df)
        assert_frame_equal(result, expected)

    def test_with_mixed_tuples(self):
        # 10697
        # columns have mixed tuples, so handle properly
        df1 = DataFrame({u'A': 'foo', (u'B', 1): 'bar'}, index=range(2))
        df2 = DataFrame({u'B': 'foo', (u'B', 1): 'bar'}, index=range(2))

        # it works
        concat([df1, df2])

    def test_handle_empty_objects(self):
        df = DataFrame(np.random.randn(10, 4), columns=list('abcd'))

        baz = df[:5].copy()
        baz['foo'] = 'bar'
        empty = df[5:5]

        frames = [baz, empty, empty, df[5:]]
        concatted = concat(frames, axis=0)

        expected = df.ix[:, ['a', 'b', 'c', 'd', 'foo']]
        expected['foo'] = expected['foo'].astype('O')
        expected.loc[0:4, 'foo'] = 'bar'

        tm.assert_frame_equal(concatted, expected)

        # empty as first element with time series
        # GH3259
        df = DataFrame(dict(A=range(10000)), index=date_range(
            '20130101', periods=10000, freq='s'))
        empty = DataFrame()
        result = concat([df, empty], axis=1)
        assert_frame_equal(result, df)
        result = concat([empty, df], axis=1)
        assert_frame_equal(result, df)

        result = concat([df, empty])
        assert_frame_equal(result, df)
        result = concat([empty, df])
        assert_frame_equal(result, df)

    def test_concat_mixed_objs(self):

        # concat mixed series/frames
        # G2385

        # axis 1
        index = date_range('01-Jan-2013', periods=10, freq='H')
        arr = np.arange(10, dtype='int64')
        s1 = Series(arr, index=index)
        s2 = Series(arr, index=index)
        df = DataFrame(arr.reshape(-1, 1), index=index)

        expected = DataFrame(np.repeat(arr, 2).reshape(-1, 2),
                             index=index, columns=[0, 0])
        result = concat([df, df], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr, 2).reshape(-1, 2),
                             index=index, columns=[0, 1])
        result = concat([s1, s2], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=[0, 1, 2])
        result = concat([s1, s2, s1], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr, 5).reshape(-1, 5),
                             index=index, columns=[0, 0, 1, 2, 3])
        result = concat([s1, df, s2, s2, s1], axis=1)
        assert_frame_equal(result, expected)

        # with names
        s1.name = 'foo'
        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=['foo', 0, 0])
        result = concat([s1, df, s2], axis=1)
        assert_frame_equal(result, expected)

        s2.name = 'bar'
        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=['foo', 0, 'bar'])
        result = concat([s1, df, s2], axis=1)
        assert_frame_equal(result, expected)

        # ignore index
        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=[0, 1, 2])
        result = concat([s1, df, s2], axis=1, ignore_index=True)
        assert_frame_equal(result, expected)

        # axis 0
        expected = DataFrame(np.tile(arr, 3).reshape(-1, 1),
                             index=index.tolist() * 3, columns=[0])
        result = concat([s1, df, s2])
        assert_frame_equal(result, expected)

        expected = DataFrame(np.tile(arr, 3).reshape(-1, 1), columns=[0])
        result = concat([s1, df, s2], ignore_index=True)
        assert_frame_equal(result, expected)

        # invalid concatente of mixed dims
        panel = tm.makePanel()
        self.assertRaises(ValueError, lambda: concat([panel, s1], axis=1))

    def test_empty_dtype_coerce(self):

        # xref to #12411
        # xref to #12045
        # xref to #11594
        # see below

        # 10571
        df1 = DataFrame(data=[[1, None], [2, None]], columns=['a', 'b'])
        df2 = DataFrame(data=[[3, None], [4, None]], columns=['a', 'b'])
        result = concat([df1, df2])
        expected = df1.dtypes
        tm.assert_series_equal(result.dtypes, expected)

    def test_dtype_coerceion(self):

        # 12411
        df = DataFrame({'date': [pd.Timestamp('20130101').tz_localize('UTC'),
                                 pd.NaT]})

        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 12045
        import datetime
        df = DataFrame({'date': [datetime.datetime(2012, 1, 1),
                                 datetime.datetime(1012, 1, 2)]})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 11594
        df = DataFrame({'text': ['some words'] + [None] * 9})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

    def test_panel_concat_other_axes(self):
        panel = tm.makePanel()

        p1 = panel.ix[:, :5, :]
        p2 = panel.ix[:, 5:, :]

        result = concat([p1, p2], axis=1)
        tm.assert_panel_equal(result, panel)

        p1 = panel.ix[:, :, :2]
        p2 = panel.ix[:, :, 2:]

        result = concat([p1, p2], axis=2)
        tm.assert_panel_equal(result, panel)

        # if things are a bit misbehaved
        p1 = panel.ix[:2, :, :2]
        p2 = panel.ix[:, :, 2:]
        p1['ItemC'] = 'baz'

        result = concat([p1, p2], axis=2)

        expected = panel.copy()
        expected['ItemC'] = expected['ItemC'].astype('O')
        expected.ix['ItemC', :, :2] = 'baz'
        tm.assert_panel_equal(result, expected)

    def test_panel_concat_buglet(self):
        # #2257
        def make_panel():
            index = 5
            cols = 3

            def df():
                return DataFrame(np.random.randn(index, cols),
                                 index=["I%s" % i for i in range(index)],
                                 columns=["C%s" % i for i in range(cols)])
            return Panel(dict([("Item%s" % x, df()) for x in ['A', 'B', 'C']]))

        panel1 = make_panel()
        panel2 = make_panel()

        panel2 = panel2.rename_axis(dict([(x, "%s_1" % x)
                                          for x in panel2.major_axis]),
                                    axis=1)

        panel3 = panel2.rename_axis(lambda x: '%s_1' % x, axis=1)
        panel3 = panel3.rename_axis(lambda x: '%s_1' % x, axis=2)

        # it works!
        concat([panel1, panel3], axis=1, verify_integrity=True)

    def test_panel4d_concat(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            p4d = tm.makePanel4D()

            p1 = p4d.ix[:, :, :5, :]
            p2 = p4d.ix[:, :, 5:, :]

            result = concat([p1, p2], axis=2)
            tm.assert_panel4d_equal(result, p4d)

            p1 = p4d.ix[:, :, :, :2]
            p2 = p4d.ix[:, :, :, 2:]

            result = concat([p1, p2], axis=3)
            tm.assert_panel4d_equal(result, p4d)

    def test_panel4d_concat_mixed_type(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            p4d = tm.makePanel4D()

            # if things are a bit misbehaved
            p1 = p4d.ix[:, :2, :, :2]
            p2 = p4d.ix[:, :, :, 2:]
            p1['L5'] = 'baz'

            result = concat([p1, p2], axis=3)

            p2['L5'] = np.nan
            expected = concat([p1, p2], axis=3)
            expected = expected.ix[result.labels]

            tm.assert_panel4d_equal(result, expected)

    def test_concat_series(self):

        ts = tm.makeTimeSeries()
        ts.name = 'foo'

        pieces = [ts[:5], ts[5:15], ts[15:]]

        result = concat(pieces)
        tm.assert_series_equal(result, ts)
        self.assertEqual(result.name, ts.name)

        result = concat(pieces, keys=[0, 1, 2])
        expected = ts.copy()

        ts.index = DatetimeIndex(np.array(ts.index.values, dtype='M8[ns]'))

        exp_labels = [np.repeat([0, 1, 2], [len(x) for x in pieces]),
                      np.arange(len(ts))]
        exp_index = MultiIndex(levels=[[0, 1, 2], ts.index],
                               labels=exp_labels)
        expected.index = exp_index
        tm.assert_series_equal(result, expected)

    def test_concat_series_axis1(self):
        ts = tm.makeTimeSeries()

        pieces = [ts[:-2], ts[2:], ts[2:-2]]

        result = concat(pieces, axis=1)
        expected = DataFrame(pieces).T
        assert_frame_equal(result, expected)

        result = concat(pieces, keys=['A', 'B', 'C'], axis=1)
        expected = DataFrame(pieces, index=['A', 'B', 'C']).T
        assert_frame_equal(result, expected)

        # preserve series names, #2489
        s = Series(randn(5), name='A')
        s2 = Series(randn(5), name='B')

        result = concat([s, s2], axis=1)
        expected = DataFrame({'A': s, 'B': s2})
        assert_frame_equal(result, expected)

        s2.name = None
        result = concat([s, s2], axis=1)
        self.assertTrue(np.array_equal(
            result.columns, Index(['A', 0], dtype='object')))

        # must reindex, #2603
        s = Series(randn(3), index=['c', 'a', 'b'], name='A')
        s2 = Series(randn(4), index=['d', 'a', 'b', 'c'], name='B')
        result = concat([s, s2], axis=1)
        expected = DataFrame({'A': s, 'B': s2})
        assert_frame_equal(result, expected)

    def test_concat_single_with_key(self):
        df = DataFrame(np.random.randn(10, 4))

        result = concat([df], keys=['foo'])
        expected = concat([df, df], keys=['foo', 'bar'])
        tm.assert_frame_equal(result, expected[:10])

    def test_concat_exclude_none(self):
        df = DataFrame(np.random.randn(10, 4))

        pieces = [df[:5], None, None, df[5:]]
        result = concat(pieces)
        tm.assert_frame_equal(result, df)
        self.assertRaises(ValueError, concat, [None, None])

    def test_concat_datetime64_block(self):
        from pandas.tseries.index import date_range

        rng = date_range('1/1/2000', periods=10)

        df = DataFrame({'time': rng})

        result = concat([df, df])
        self.assertTrue((result.iloc[:10]['time'] == rng).all())
        self.assertTrue((result.iloc[10:]['time'] == rng).all())

    def test_concat_timedelta64_block(self):
        from pandas import to_timedelta

        rng = to_timedelta(np.arange(10), unit='s')

        df = DataFrame({'time': rng})

        result = concat([df, df])
        self.assertTrue((result.iloc[:10]['time'] == rng).all())
        self.assertTrue((result.iloc[10:]['time'] == rng).all())

    def test_concat_keys_with_none(self):
        # #1649
        df0 = DataFrame([[10, 20, 30], [10, 20, 30], [10, 20, 30]])

        result = concat(dict(a=None, b=df0, c=df0[:2], d=df0[:1], e=df0))
        expected = concat(dict(b=df0, c=df0[:2], d=df0[:1], e=df0))
        tm.assert_frame_equal(result, expected)

        result = concat([None, df0, df0[:2], df0[:1], df0],
                        keys=['a', 'b', 'c', 'd', 'e'])
        expected = concat([df0, df0[:2], df0[:1], df0],
                          keys=['b', 'c', 'd', 'e'])
        tm.assert_frame_equal(result, expected)

    def test_union_categorical(self):
        # GH 13361
        data = [
            (list('abc'), list('abd'), list('abcabd')),
            ([0, 1, 2], [2, 3, 4], [0, 1, 2, 2, 3, 4]),
            ([0, 1.2, 2], [2, 3.4, 4], [0, 1.2, 2, 2, 3.4, 4]),

            (['b', 'b', np.nan, 'a'], ['a', np.nan, 'c'],
             ['b', 'b', np.nan, 'a', 'a', np.nan, 'c']),

            (pd.date_range('2014-01-01', '2014-01-05'),
             pd.date_range('2014-01-06', '2014-01-07'),
             pd.date_range('2014-01-01', '2014-01-07')),

            (pd.date_range('2014-01-01', '2014-01-05', tz='US/Central'),
             pd.date_range('2014-01-06', '2014-01-07', tz='US/Central'),
             pd.date_range('2014-01-01', '2014-01-07', tz='US/Central')),

            (pd.period_range('2014-01-01', '2014-01-05'),
             pd.period_range('2014-01-06', '2014-01-07'),
             pd.period_range('2014-01-01', '2014-01-07')),
        ]

        for a, b, combined in data:
            result = union_categoricals([Categorical(a), Categorical(b)])
            expected = Categorical(combined)
            tm.assert_categorical_equal(result, expected,
                                        check_category_order=True)

        # new categories ordered by appearance
        s = Categorical(['x', 'y', 'z'])
        s2 = Categorical(['a', 'b', 'c'])
        result = union_categoricals([s, s2])
        expected = Categorical(['x', 'y', 'z', 'a', 'b', 'c'],
                               categories=['x', 'y', 'z', 'a', 'b', 'c'])
        tm.assert_categorical_equal(result, expected)

        s = Categorical([0, 1.2, 2], ordered=True)
        s2 = Categorical([0, 1.2, 2], ordered=True)
        result = union_categoricals([s, s2])
        expected = Categorical([0, 1.2, 2, 0, 1.2, 2], ordered=True)
        tm.assert_categorical_equal(result, expected)

        # must exactly match types
        s = Categorical([0, 1.2, 2])
        s2 = Categorical([2, 3, 4])
        msg = 'dtype of categories must be the same'
        with tm.assertRaisesRegexp(TypeError, msg):
            union_categoricals([s, s2])

        msg = 'No Categoricals to union'
        with tm.assertRaisesRegexp(ValueError, msg):
            union_categoricals([])

    def test_union_categoricals_nan(self):
        # GH 13759
        res = union_categoricals([pd.Categorical([1, 2, np.nan]),
                                  pd.Categorical([3, 2, np.nan])])
        exp = Categorical([1, 2, np.nan, 3, 2, np.nan])
        tm.assert_categorical_equal(res, exp)

        res = union_categoricals([pd.Categorical(['A', 'B']),
                                  pd.Categorical(['B', 'B', np.nan])])
        exp = Categorical(['A', 'B', 'B', 'B', np.nan])
        tm.assert_categorical_equal(res, exp)

        val1 = [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-03-01'),
                pd.NaT]
        val2 = [pd.NaT, pd.Timestamp('2011-01-01'),
                pd.Timestamp('2011-02-01')]

        res = union_categoricals([pd.Categorical(val1), pd.Categorical(val2)])
        exp = Categorical(val1 + val2,
                          categories=[pd.Timestamp('2011-01-01'),
                                      pd.Timestamp('2011-03-01'),
                                      pd.Timestamp('2011-02-01')])
        tm.assert_categorical_equal(res, exp)

        # all NaN
        res = union_categoricals([pd.Categorical([np.nan, np.nan]),
                                  pd.Categorical(['X'])])
        exp = Categorical([np.nan, np.nan, 'X'])
        tm.assert_categorical_equal(res, exp)

        res = union_categoricals([pd.Categorical([np.nan, np.nan]),
                                  pd.Categorical([np.nan, np.nan])])
        exp = Categorical([np.nan, np.nan, np.nan, np.nan])
        tm.assert_categorical_equal(res, exp)

    def test_union_categoricals_empty(self):
        # GH 13759
        res = union_categoricals([pd.Categorical([]),
                                  pd.Categorical([])])
        exp = Categorical([])
        tm.assert_categorical_equal(res, exp)

        res = union_categoricals([pd.Categorical([]),
                                  pd.Categorical([1.0])])
        exp = Categorical([1.0])
        tm.assert_categorical_equal(res, exp)

        # to make dtype equal
        nanc = pd.Categorical(np.array([np.nan], dtype=np.float64))
        res = union_categoricals([nanc,
                                  pd.Categorical([])])
        tm.assert_categorical_equal(res, nanc)

    def test_union_categorical_same_category(self):
        # check fastpath
        c1 = Categorical([1, 2, 3, 4], categories=[1, 2, 3, 4])
        c2 = Categorical([3, 2, 1, np.nan], categories=[1, 2, 3, 4])
        res = union_categoricals([c1, c2])
        exp = Categorical([1, 2, 3, 4, 3, 2, 1, np.nan],
                          categories=[1, 2, 3, 4])
        tm.assert_categorical_equal(res, exp)

        c1 = Categorical(['z', 'z', 'z'], categories=['x', 'y', 'z'])
        c2 = Categorical(['x', 'x', 'x'], categories=['x', 'y', 'z'])
        res = union_categoricals([c1, c2])
        exp = Categorical(['z', 'z', 'z', 'x', 'x', 'x'],
                          categories=['x', 'y', 'z'])
        tm.assert_categorical_equal(res, exp)

    def test_union_categoricals_ordered(self):
        c1 = Categorical([1, 2, 3], ordered=True)
        c2 = Categorical([1, 2, 3], ordered=False)

        msg = 'Categorical.ordered must be the same'
        with tm.assertRaisesRegexp(TypeError, msg):
            union_categoricals([c1, c2])

        res = union_categoricals([c1, c1])
        exp = Categorical([1, 2, 3, 1, 2, 3], ordered=True)
        tm.assert_categorical_equal(res, exp)

        c1 = Categorical([1, 2, 3, np.nan], ordered=True)
        c2 = Categorical([3, 2], categories=[1, 2, 3], ordered=True)

        res = union_categoricals([c1, c2])
        exp = Categorical([1, 2, 3, np.nan, 3, 2], ordered=True)
        tm.assert_categorical_equal(res, exp)

        c1 = Categorical([1, 2, 3], ordered=True)
        c2 = Categorical([1, 2, 3], categories=[3, 2, 1], ordered=True)

        msg = "to union ordered Categoricals, all categories must be the same"
        with tm.assertRaisesRegexp(TypeError, msg):
            union_categoricals([c1, c2])

    def test_union_categoricals_sort(self):
        # GH 13846
        c1 = Categorical(['x', 'y', 'z'])
        c2 = Categorical(['a', 'b', 'c'])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical(['x', 'y', 'z', 'a', 'b', 'c'],
                               categories=['a', 'b', 'c', 'x', 'y', 'z'])
        tm.assert_categorical_equal(result, expected)

        # fastpath
        c1 = Categorical(['a', 'b'], categories=['b', 'a', 'c'])
        c2 = Categorical(['b', 'c'], categories=['b', 'a', 'c'])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical(['a', 'b', 'b', 'c'],
                               categories=['a', 'b', 'c'])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical(['a', 'b'], categories=['c', 'a', 'b'])
        c2 = Categorical(['b', 'c'], categories=['c', 'a', 'b'])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical(['a', 'b', 'b', 'c'],
                               categories=['a', 'b', 'c'])
        tm.assert_categorical_equal(result, expected)

        # fastpath - skip resort
        c1 = Categorical(['a', 'b'], categories=['a', 'b', 'c'])
        c2 = Categorical(['b', 'c'], categories=['a', 'b', 'c'])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical(['a', 'b', 'b', 'c'],
                               categories=['a', 'b', 'c'])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical(['x', np.nan])
        c2 = Categorical([np.nan, 'b'])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical(['x', np.nan, np.nan, 'b'],
                               categories=['b', 'x'])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical([np.nan])
        c2 = Categorical([np.nan])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical([np.nan, np.nan], categories=[])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical([])
        c2 = Categorical([])
        result = union_categoricals([c1, c2], sort_categories=True)
        expected = Categorical([])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical(['b', 'a'], categories=['b', 'a', 'c'], ordered=True)
        c2 = Categorical(['a', 'c'], categories=['b', 'a', 'c'], ordered=True)
        with tm.assertRaises(TypeError):
            union_categoricals([c1, c2], sort_categories=True)

    def test_union_categoricals_sort_false(self):
        # GH 13846
        c1 = Categorical(['x', 'y', 'z'])
        c2 = Categorical(['a', 'b', 'c'])
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical(['x', 'y', 'z', 'a', 'b', 'c'],
                               categories=['x', 'y', 'z', 'a', 'b', 'c'])
        tm.assert_categorical_equal(result, expected)

        # fastpath
        c1 = Categorical(['a', 'b'], categories=['b', 'a', 'c'])
        c2 = Categorical(['b', 'c'], categories=['b', 'a', 'c'])
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical(['a', 'b', 'b', 'c'],
                               categories=['b', 'a', 'c'])
        tm.assert_categorical_equal(result, expected)

        # fastpath - skip resort
        c1 = Categorical(['a', 'b'], categories=['a', 'b', 'c'])
        c2 = Categorical(['b', 'c'], categories=['a', 'b', 'c'])
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical(['a', 'b', 'b', 'c'],
                               categories=['a', 'b', 'c'])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical(['x', np.nan])
        c2 = Categorical([np.nan, 'b'])
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical(['x', np.nan, np.nan, 'b'],
                               categories=['x', 'b'])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical([np.nan])
        c2 = Categorical([np.nan])
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical([np.nan, np.nan], categories=[])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical([])
        c2 = Categorical([])
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical([])
        tm.assert_categorical_equal(result, expected)

        c1 = Categorical(['b', 'a'], categories=['b', 'a', 'c'], ordered=True)
        c2 = Categorical(['a', 'c'], categories=['b', 'a', 'c'], ordered=True)
        result = union_categoricals([c1, c2], sort_categories=False)
        expected = Categorical(['b', 'a', 'a', 'c'],
                               categories=['b', 'a', 'c'], ordered=True)
        tm.assert_categorical_equal(result, expected)

    def test_concat_bug_1719(self):
        ts1 = tm.makeTimeSeries()
        ts2 = tm.makeTimeSeries()[::2]

        # to join with union
        # these two are of different length!
        left = concat([ts1, ts2], join='outer', axis=1)
        right = concat([ts2, ts1], join='outer', axis=1)

        self.assertEqual(len(left), len(right))

    def test_concat_bug_2972(self):
        ts0 = Series(np.zeros(5))
        ts1 = Series(np.ones(5))
        ts0.name = ts1.name = 'same name'
        result = concat([ts0, ts1], axis=1)

        expected = DataFrame({0: ts0, 1: ts1})
        expected.columns = ['same name', 'same name']
        assert_frame_equal(result, expected)

    def test_concat_bug_3602(self):

        # GH 3602, duplicate columns
        df1 = DataFrame({'firmNo': [0, 0, 0, 0], 'stringvar': [
                        'rrr', 'rrr', 'rrr', 'rrr'], 'prc': [6, 6, 6, 6]})
        df2 = DataFrame({'misc': [1, 2, 3, 4], 'prc': [
                        6, 6, 6, 6], 'C': [9, 10, 11, 12]})
        expected = DataFrame([[0, 6, 'rrr', 9, 1, 6],
                              [0, 6, 'rrr', 10, 2, 6],
                              [0, 6, 'rrr', 11, 3, 6],
                              [0, 6, 'rrr', 12, 4, 6]])
        expected.columns = ['firmNo', 'prc', 'stringvar', 'C', 'misc', 'prc']

        result = concat([df1, df2], axis=1)
        assert_frame_equal(result, expected)

    def test_concat_series_axis1_same_names_ignore_index(self):
        dates = date_range('01-Jan-2013', '01-Jan-2014', freq='MS')[0:-1]
        s1 = Series(randn(len(dates)), index=dates, name='value')
        s2 = Series(randn(len(dates)), index=dates, name='value')

        result = concat([s1, s2], axis=1, ignore_index=True)
        self.assertTrue(np.array_equal(result.columns, [0, 1]))

    def test_concat_iterables(self):
        from collections import deque, Iterable

        # GH8645 check concat works with tuples, list, generators, and weird
        # stuff like deque and custom iterables
        df1 = DataFrame([1, 2, 3])
        df2 = DataFrame([4, 5, 6])
        expected = DataFrame([1, 2, 3, 4, 5, 6])
        assert_frame_equal(concat((df1, df2), ignore_index=True), expected)
        assert_frame_equal(concat([df1, df2], ignore_index=True), expected)
        assert_frame_equal(concat((df for df in (df1, df2)),
                                  ignore_index=True), expected)
        assert_frame_equal(
            concat(deque((df1, df2)), ignore_index=True), expected)

        class CustomIterator1(object):

            def __len__(self):
                return 2

            def __getitem__(self, index):
                try:
                    return {0: df1, 1: df2}[index]
                except KeyError:
                    raise IndexError
        assert_frame_equal(pd.concat(CustomIterator1(),
                                     ignore_index=True), expected)

        class CustomIterator2(Iterable):

            def __iter__(self):
                yield df1
                yield df2
        assert_frame_equal(pd.concat(CustomIterator2(),
                                     ignore_index=True), expected)

    def test_concat_invalid(self):

        # trying to concat a ndframe with a non-ndframe
        df1 = mkdf(10, 2)
        for obj in [1, dict(), [1, 2], (1, 2)]:
            self.assertRaises(TypeError, lambda x: concat([df1, obj]))

    def test_concat_invalid_first_argument(self):
        df1 = mkdf(10, 2)
        df2 = mkdf(10, 2)
        self.assertRaises(TypeError, concat, df1, df2)

        # generator ok though
        concat(DataFrame(np.random.rand(5, 5)) for _ in range(3))

        # text reader ok
        # GH6583
        data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

        reader = read_csv(StringIO(data), chunksize=1)
        result = concat(reader, ignore_index=True)
        expected = read_csv(StringIO(data))
        assert_frame_equal(result, expected)

    def test_concat_NaT_series(self):
        # GH 11693
        # test for merging NaT series with datetime series.
        x = Series(date_range('20151124 08:00', '20151124 09:00',
                              freq='1h', tz='US/Eastern'))
        y = Series(pd.NaT, index=[0, 1], dtype='datetime64[ns, US/Eastern]')
        expected = Series([x[0], x[1], pd.NaT, pd.NaT])

        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # all NaT with tz
        expected = Series(pd.NaT, index=range(4),
                          dtype='datetime64[ns, US/Eastern]')
        result = pd.concat([y, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # without tz
        x = pd.Series(pd.date_range('20151124 08:00',
                                    '20151124 09:00', freq='1h'))
        y = pd.Series(pd.date_range('20151124 10:00',
                                    '20151124 11:00', freq='1h'))
        y[:] = pd.NaT
        expected = pd.Series([x[0], x[1], pd.NaT, pd.NaT])
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # all NaT without tz
        x[:] = pd.NaT
        expected = pd.Series(pd.NaT, index=range(4),
                             dtype='datetime64[ns]')
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    def test_concat_tz_frame(self):
        df2 = DataFrame(dict(A=pd.Timestamp('20130102', tz='US/Eastern'),
                             B=pd.Timestamp('20130603', tz='CET')),
                        index=range(5))

        # concat
        df3 = pd.concat([df2.A.to_frame(), df2.B.to_frame()], axis=1)
        assert_frame_equal(df2, df3)

    def test_concat_tz_series(self):
        # GH 11755
        # tz and no tz
        x = Series(date_range('20151124 08:00',
                              '20151124 09:00',
                              freq='1h', tz='UTC'))
        y = Series(date_range('2012-01-01', '2012-01-02'))
        expected = Series([x[0], x[1], y[0], y[1]],
                          dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # GH 11887
        # concat tz and object
        x = Series(date_range('20151124 08:00',
                              '20151124 09:00',
                              freq='1h', tz='UTC'))
        y = Series(['a', 'b'])
        expected = Series([x[0], x[1], y[0], y[1]],
                          dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # 12217
        # 12306 fixed I think

        # Concat'ing two UTC times
        first = pd.DataFrame([[datetime(2016, 1, 1)]])
        first[0] = first[0].dt.tz_localize('UTC')

        second = pd.DataFrame([[datetime(2016, 1, 2)]])
        second[0] = second[0].dt.tz_localize('UTC')

        result = pd.concat([first, second])
        self.assertEqual(result[0].dtype, 'datetime64[ns, UTC]')

        # Concat'ing two London times
        first = pd.DataFrame([[datetime(2016, 1, 1)]])
        first[0] = first[0].dt.tz_localize('Europe/London')

        second = pd.DataFrame([[datetime(2016, 1, 2)]])
        second[0] = second[0].dt.tz_localize('Europe/London')

        result = pd.concat([first, second])
        self.assertEqual(result[0].dtype, 'datetime64[ns, Europe/London]')

        # Concat'ing 2+1 London times
        first = pd.DataFrame([[datetime(2016, 1, 1)], [datetime(2016, 1, 2)]])
        first[0] = first[0].dt.tz_localize('Europe/London')

        second = pd.DataFrame([[datetime(2016, 1, 3)]])
        second[0] = second[0].dt.tz_localize('Europe/London')

        result = pd.concat([first, second])
        self.assertEqual(result[0].dtype, 'datetime64[ns, Europe/London]')

        # Concat'ing 1+2 London times
        first = pd.DataFrame([[datetime(2016, 1, 1)]])
        first[0] = first[0].dt.tz_localize('Europe/London')

        second = pd.DataFrame([[datetime(2016, 1, 2)], [datetime(2016, 1, 3)]])
        second[0] = second[0].dt.tz_localize('Europe/London')

        result = pd.concat([first, second])
        self.assertEqual(result[0].dtype, 'datetime64[ns, Europe/London]')

    def test_concat_tz_series_with_datetimelike(self):
        # GH 12620
        # tz and timedelta
        x = [pd.Timestamp('2011-01-01', tz='US/Eastern'),
             pd.Timestamp('2011-02-01', tz='US/Eastern')]
        y = [pd.Timedelta('1 day'), pd.Timedelta('2 day')]
        result = concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y, dtype='object'))

        # tz and period
        y = [pd.Period('2011-03', freq='M'), pd.Period('2011-04', freq='M')]
        result = concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y, dtype='object'))

    def test_concat_tz_series_tzlocal(self):
        # GH 13583
        tm._skip_if_no_dateutil()
        import dateutil
        x = [pd.Timestamp('2011-01-01', tz=dateutil.tz.tzlocal()),
             pd.Timestamp('2011-02-01', tz=dateutil.tz.tzlocal())]
        y = [pd.Timestamp('2012-01-01', tz=dateutil.tz.tzlocal()),
             pd.Timestamp('2012-02-01', tz=dateutil.tz.tzlocal())]
        result = concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y))
        self.assertEqual(result.dtype, 'datetime64[ns, tzlocal()]')

    def test_concat_period_series(self):
        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.PeriodIndex(['2015-10-01', '2016-01-01'], freq='D'))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        self.assertEqual(result.dtype, 'object')

        # different freq
        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.PeriodIndex(['2015-10-01', '2016-01-01'], freq='M'))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        self.assertEqual(result.dtype, 'object')

        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='M'))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        self.assertEqual(result.dtype, 'object')

        # non-period
        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.DatetimeIndex(['2015-11-01', '2015-12-01']))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        self.assertEqual(result.dtype, 'object')

        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(['A', 'B'])
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        self.assertEqual(result.dtype, 'object')

    def test_concat_empty_series(self):
        # GH 11082
        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series(name='y')
        res = pd.concat([s1, s2], axis=1)
        exp = pd.DataFrame({'x': [1, 2, 3], 'y': [np.nan, np.nan, np.nan]})
        tm.assert_frame_equal(res, exp)

        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series(name='y')
        res = pd.concat([s1, s2], axis=0)
        # name will be reset
        exp = pd.Series([1, 2, 3])
        tm.assert_series_equal(res, exp)

        # empty Series with no name
        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series(name=None)
        res = pd.concat([s1, s2], axis=1)
        exp = pd.DataFrame({'x': [1, 2, 3], 0: [np.nan, np.nan, np.nan]},
                           columns=['x', 0])
        tm.assert_frame_equal(res, exp)

    def test_default_index(self):
        # is_series and ignore_index
        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series([4, 5, 6], name='y')
        res = pd.concat([s1, s2], axis=1, ignore_index=True)
        self.assertIsInstance(res.columns, pd.RangeIndex)
        exp = pd.DataFrame([[1, 4], [2, 5], [3, 6]])
        # use check_index_type=True to check the result have
        # RangeIndex (default index)
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

        # is_series and all inputs have no names
        s1 = pd.Series([1, 2, 3])
        s2 = pd.Series([4, 5, 6])
        res = pd.concat([s1, s2], axis=1, ignore_index=False)
        self.assertIsInstance(res.columns, pd.RangeIndex)
        exp = pd.DataFrame([[1, 4], [2, 5], [3, 6]])
        exp.columns = pd.RangeIndex(2)
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

        # is_dataframe and ignore_index
        df1 = pd.DataFrame({'A': [1, 2], 'B': [5, 6]})
        df2 = pd.DataFrame({'A': [3, 4], 'B': [7, 8]})

        res = pd.concat([df1, df2], axis=0, ignore_index=True)
        exp = pd.DataFrame([[1, 5], [2, 6], [3, 7], [4, 8]],
                           columns=['A', 'B'])
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

        res = pd.concat([df1, df2], axis=1, ignore_index=True)
        exp = pd.DataFrame([[1, 5, 3, 7], [2, 6, 4, 8]])
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

    def test_concat_multiindex_rangeindex(self):
        # GH13542
        # when multi-index levels are RangeIndex objects
        # there is a bug in concat with objects of len 1

        df = DataFrame(np.random.randn(9, 2))
        df.index = MultiIndex(levels=[pd.RangeIndex(3), pd.RangeIndex(3)],
                              labels=[np.repeat(np.arange(3), 3),
                                      np.tile(np.arange(3), 3)])

        res = concat([df.iloc[[2, 3, 4], :], df.iloc[[5], :]])
        exp = df.iloc[[2, 3, 4, 5], :]
        tm.assert_frame_equal(res, exp)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
