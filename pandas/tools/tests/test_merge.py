# pylint: disable=E1103

import nose

from datetime import datetime
from numpy.random import randn
from numpy import nan
import numpy as np
import random

import pandas as pd
from pandas.compat import range, lrange, lzip, zip, StringIO
from pandas import compat
from pandas.tseries.index import DatetimeIndex
from pandas.tools.merge import merge, concat, ordered_merge, MergeError
from pandas.util.testing import (assert_frame_equal, assert_series_equal,
                                 assert_almost_equal,
                                 makeCustomDataframe as mkdf,
                                 assertRaisesRegexp)
from pandas import isnull, DataFrame, Index, MultiIndex, Panel, Series, date_range, read_table, read_csv
import pandas.algos as algos
import pandas.util.testing as tm

a_ = np.array

N = 50
NGROUPS = 8
JOIN_TYPES = ['inner', 'outer', 'left', 'right']


def get_test_data(ngroups=NGROUPS, n=N):
    unique_groups = lrange(ngroups)
    arr = np.asarray(np.tile(unique_groups, n // ngroups))

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)])

    random.shuffle(arr)
    return arr


class TestMerge(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        # aggregate multiple columns
        self.df = DataFrame({'key1': get_test_data(),
                             'key2': get_test_data(),
                             'data1': np.random.randn(N),
                             'data2': np.random.randn(N)})

        # exclude a couple keys for fun
        self.df = self.df[self.df['key2'] > 1]

        self.df2 = DataFrame({'key1': get_test_data(n=N // 5),
                              'key2': get_test_data(ngroups=NGROUPS // 2,
                                                    n=N // 5),
                              'value': np.random.randn(N // 5)})

        index, data = tm.getMixedTypeDict()
        self.target = DataFrame(data, index=index)

        # Join on string value
        self.source = DataFrame({'MergedA': data['A'], 'MergedD': data['D']},
                                index=data['C'])

        self.left = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'e', 'a'],
                               'v1': np.random.randn(7)})
        self.right = DataFrame({'v2': np.random.randn(4)},
                               index=['d', 'b', 'c', 'a'])

    def test_cython_left_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype=np.int64)
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype=np.int64)
        max_group = 5

        ls, rs = algos.left_outer_join(left, right, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        exp_li = a_([0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                     6, 6, 7, 7, 8, 8, 9, 10])
        exp_ri = a_([0, 0, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                     4, 5, 4, 5, 4, 5, -1, -1])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        self.assert_numpy_array_equal(ls, exp_ls)
        self.assert_numpy_array_equal(rs, exp_rs)

    def test_cython_right_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype=np.int64)
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype=np.int64)
        max_group = 5

        rs, ls = algos.left_outer_join(right, left, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        #            0        1        1        1
        exp_li = a_([0, 1, 2, 3, 4, 5, 3, 4, 5, 3, 4, 5,
                     #            2        2        4
                     6, 7, 8, 6, 7, 8, -1])
        exp_ri = a_([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
                     4, 4, 4, 5, 5, 5, 6])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        self.assert_numpy_array_equal(ls, exp_ls)
        self.assert_numpy_array_equal(rs, exp_rs)

    def test_cython_inner_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype=np.int64)
        right = a_([1, 1, 0, 4, 2, 2, 1, 4], dtype=np.int64)
        max_group = 5

        ls, rs = algos.inner_join(left, right, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        exp_li = a_([0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                     6, 6, 7, 7, 8, 8])
        exp_ri = a_([0, 0, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                     4, 5, 4, 5, 4, 5])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        self.assert_numpy_array_equal(ls, exp_ls)
        self.assert_numpy_array_equal(rs, exp_rs)

    def test_left_outer_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='left')

        joined_both = merge(self.df, self.df2)
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='left')

    def test_right_outer_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2', how='right')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='right')

        joined_both = merge(self.df, self.df2, how='right')
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='right')

    def test_full_outer_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2', how='outer')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='outer')

        joined_both = merge(self.df, self.df2, how='outer')
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='outer')

    def test_inner_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2', how='inner')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='inner')

        joined_both = merge(self.df, self.df2, how='inner')
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='inner')

    def test_handle_overlap(self):
        joined = merge(self.df, self.df2, on='key2',
                       suffixes=['.foo', '.bar'])

        self.assertIn('key1.foo', joined)
        self.assertIn('key1.bar', joined)

    def test_handle_overlap_arbitrary_key(self):
        joined = merge(self.df, self.df2,
                       left_on='key2', right_on='key1',
                       suffixes=['.foo', '.bar'])
        self.assertIn('key1.foo', joined)
        self.assertIn('key2.bar', joined)

    def test_merge_common(self):
        joined = merge(self.df, self.df2)
        exp = merge(self.df, self.df2, on=['key1', 'key2'])
        tm.assert_frame_equal(joined, exp)

    def test_join_on(self):
        target = self.target
        source = self.source

        merged = target.join(source, on='C')
        self.assert_numpy_array_equal(merged['MergedA'], target['A'])
        self.assert_numpy_array_equal(merged['MergedD'], target['D'])

        # join with duplicates (fix regression from DataFrame/Matrix merge)
        df = DataFrame({'key': ['a', 'a', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1, 2]}, index=['a', 'b', 'c'])
        joined = df.join(df2, on='key')
        expected = DataFrame({'key': ['a', 'a', 'b', 'b', 'c'],
                              'value': [0, 0, 1, 1, 2]})
        assert_frame_equal(joined, expected)

        # Test when some are missing
        df_a = DataFrame([[1], [2], [3]], index=['a', 'b', 'c'],
                         columns=['one'])
        df_b = DataFrame([['foo'], ['bar']], index=[1, 2],
                         columns=['two'])
        df_c = DataFrame([[1], [2]], index=[1, 2],
                         columns=['three'])
        joined = df_a.join(df_b, on='one')
        joined = joined.join(df_c, on='one')
        self.assertTrue(np.isnan(joined['two']['c']))
        self.assertTrue(np.isnan(joined['three']['c']))

        # merge column not p resent
        self.assertRaises(KeyError, target.join, source, on='E')

        # overlap
        source_copy = source.copy()
        source_copy['A'] = 0
        self.assertRaises(ValueError, target.join, source_copy, on='A')

    def test_join_on_fails_with_different_right_index(self):
        with tm.assertRaises(ValueError):
            df = DataFrame({'a': tm.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)})
            df2 = DataFrame({'a': tm.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)},
                            index=tm.makeCustomIndex(10, 2))
            merge(df, df2, left_on='a', right_index=True)

    def test_join_on_fails_with_different_left_index(self):
        with tm.assertRaises(ValueError):
            df = DataFrame({'a': tm.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)},
                           index=tm.makeCustomIndex(10, 2))
            df2 = DataFrame({'a': tm.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)})
            merge(df, df2, right_on='b', left_index=True)

    def test_join_on_fails_with_different_column_counts(self):
        with tm.assertRaises(ValueError):
            df = DataFrame({'a': tm.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)})
            df2 = DataFrame({'a': tm.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)},
                            index=tm.makeCustomIndex(10, 2))
            merge(df, df2, right_on='a', left_on=['a', 'b'])

    def test_join_on_pass_vector(self):
        expected = self.target.join(self.source, on='C')
        del expected['C']

        join_col = self.target.pop('C')
        result = self.target.join(self.source, on=join_col)
        assert_frame_equal(result, expected)

    def test_join_with_len0(self):
        # nothing to merge
        merged = self.target.join(self.source.reindex([]), on='C')
        for col in self.source:
            self.assertIn(col, merged)
            self.assertTrue(merged[col].isnull().all())

        merged2 = self.target.join(self.source.reindex([]), on='C',
                                   how='inner')
        self.assertTrue(merged2.columns.equals(merged.columns))
        self.assertEqual(len(merged2), 0)

    def test_join_on_inner(self):
        df = DataFrame({'key': ['a', 'a', 'd', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1]}, index=['a', 'b'])

        joined = df.join(df2, on='key', how='inner')

        expected = df.join(df2, on='key')
        expected = expected[expected['value'].notnull()]
        self.assert_numpy_array_equal(joined['key'], expected['key'])
        self.assert_numpy_array_equal(joined['value'], expected['value'])
        self.assertTrue(joined.index.equals(expected.index))

    def test_join_on_singlekey_list(self):
        df = DataFrame({'key': ['a', 'a', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1, 2]}, index=['a', 'b', 'c'])

        # corner cases
        joined = df.join(df2, on=['key'])
        expected = df.join(df2, on='key')

        assert_frame_equal(joined, expected)

    def test_join_on_series(self):
        result = self.target.join(self.source['MergedA'], on='C')
        expected = self.target.join(self.source[['MergedA']], on='C')
        assert_frame_equal(result, expected)

    def test_join_on_series_buglet(self):
        # GH #638
        df = DataFrame({'a': [1, 1]})
        ds = Series([2], index=[1], name='b')
        result = df.join(ds, on='a')
        expected = DataFrame({'a': [1, 1],
                              'b': [2, 2]}, index=df.index)
        tm.assert_frame_equal(result, expected)

    def test_join_index_mixed(self):
        df1 = DataFrame({'A': 1., 'B': 2, 'C': 'foo', 'D': True},
                        index=np.arange(10),
                        columns=['A', 'B', 'C', 'D'])
        self.assertEqual(df1['B'].dtype, np.int64)
        self.assertEqual(df1['D'].dtype, np.bool_)

        df2 = DataFrame({'A': 1., 'B': 2, 'C': 'foo', 'D': True},
                        index=np.arange(0, 10, 2),
                        columns=['A', 'B', 'C', 'D'])

        # overlap
        joined = df1.join(df2, lsuffix='_one', rsuffix='_two')
        expected_columns = ['A_one', 'B_one', 'C_one', 'D_one',
                            'A_two', 'B_two', 'C_two', 'D_two']
        df1.columns = expected_columns[:4]
        df2.columns = expected_columns[4:]
        expected = _join_by_hand(df1, df2)
        assert_frame_equal(joined, expected)

        # no overlapping blocks
        df1 = DataFrame(index=np.arange(10))
        df1['bool'] = True
        df1['string'] = 'foo'

        df2 = DataFrame(index=np.arange(5, 15))
        df2['int'] = 1
        df2['float'] = 1.

        for kind in JOIN_TYPES:

            joined = df1.join(df2, how=kind)
            expected = _join_by_hand(df1, df2, how=kind)
            assert_frame_equal(joined, expected)

            joined = df2.join(df1, how=kind)
            expected = _join_by_hand(df2, df1, how=kind)
            assert_frame_equal(joined, expected)

    def test_join_empty_bug(self):
        # generated an exception in 0.4.3
        x = DataFrame()
        x.join(DataFrame([3], index=[0], columns=['A']), how='outer')

    def test_join_unconsolidated(self):
        # GH #331
        a = DataFrame(randn(30, 2), columns=['a', 'b'])
        c = Series(randn(30))
        a['c'] = c
        d = DataFrame(randn(30, 1), columns=['q'])

        # it works!
        a.join(d)
        d.join(a)

    def test_join_multiindex(self):
        index1 = MultiIndex.from_arrays([['a', 'a', 'a', 'b', 'b', 'b'],
                                         [1, 2, 3, 1, 2, 3]],
                                        names=['first', 'second'])

        index2 = MultiIndex.from_arrays([['b', 'b', 'b', 'c', 'c', 'c'],
                                         [1, 2, 3, 1, 2, 3]],
                                        names=['first', 'second'])

        df1 = DataFrame(data=np.random.randn(6), index=index1,
                        columns=['var X'])
        df2 = DataFrame(data=np.random.randn(6), index=index2,
                        columns=['var Y'])

        df1 = df1.sortlevel(0)
        df2 = df2.sortlevel(0)

        joined = df1.join(df2, how='outer')
        ex_index = index1._tuple_index.union(index2._tuple_index)
        expected = df1.reindex(ex_index).join(df2.reindex(ex_index))
        expected.index.names = index1.names
        assert_frame_equal(joined, expected)
        self.assertEqual(joined.index.names, index1.names)

        df1 = df1.sortlevel(1)
        df2 = df2.sortlevel(1)

        joined = df1.join(df2, how='outer').sortlevel(0)
        ex_index = index1._tuple_index.union(index2._tuple_index)
        expected = df1.reindex(ex_index).join(df2.reindex(ex_index))
        expected.index.names = index1.names

        assert_frame_equal(joined, expected)
        self.assertEqual(joined.index.names, index1.names)

    def test_join_inner_multiindex(self):
        key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
                'qux', 'snap']
        key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
                'three', 'one']

        data = np.random.randn(len(key1))
        data = DataFrame({'key1': key1, 'key2': key2,
                         'data': data})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        to_join = DataFrame(np.random.randn(10, 3), index=index,
                            columns=['j_one', 'j_two', 'j_three'])

        joined = data.join(to_join, on=['key1', 'key2'], how='inner')
        expected = merge(data, to_join.reset_index(),
                         left_on=['key1', 'key2'],
                         right_on=['first', 'second'], how='inner',
                         sort=False)

        expected2 = merge(to_join, data,
                          right_on=['key1', 'key2'], left_index=True,
                          how='inner', sort=False)
        assert_frame_equal(joined, expected2.reindex_like(joined))

        expected2 = merge(to_join, data, right_on=['key1', 'key2'],
                          left_index=True, how='inner', sort=False)

        expected = expected.drop(['first', 'second'], axis=1)
        expected.index = joined.index

        self.assertTrue(joined.index.is_monotonic)
        assert_frame_equal(joined, expected)

        # _assert_same_contents(expected, expected2.ix[:, expected.columns])

    def test_join_hierarchical_mixed(self):
        df = DataFrame([(1, 2, 3), (4, 5, 6)], columns=['a', 'b', 'c'])
        new_df = df.groupby(['a']).agg({'b': [np.mean, np.sum]})
        other_df = DataFrame(
            [(1, 2, 3), (7, 10, 6)], columns=['a', 'b', 'd'])
        other_df.set_index('a', inplace=True)

        result = merge(new_df, other_df, left_index=True, right_index=True)
        self.assertTrue(('b', 'mean') in result)
        self.assertTrue('b' in result)

    def test_join_float64_float32(self):

        a = DataFrame(randn(10, 2), columns=['a', 'b'], dtype = np.float64)
        b = DataFrame(randn(10, 1), columns=['c'], dtype = np.float32)
        joined = a.join(b)
        self.assertEqual(joined.dtypes['a'], 'float64')
        self.assertEqual(joined.dtypes['b'], 'float64')
        self.assertEqual(joined.dtypes['c'], 'float32')

        a = np.random.randint(0, 5, 100).astype('int64')
        b = np.random.random(100).astype('float64')
        c = np.random.random(100).astype('float32')
        df = DataFrame({'a': a, 'b': b, 'c': c})
        xpdf = DataFrame({'a': a, 'b': b, 'c': c })
        s = DataFrame(np.random.random(5).astype('float32'), columns=['md'])
        rs = df.merge(s, left_on='a', right_index=True)
        self.assertEqual(rs.dtypes['a'], 'int64')
        self.assertEqual(rs.dtypes['b'], 'float64')
        self.assertEqual(rs.dtypes['c'], 'float32')
        self.assertEqual(rs.dtypes['md'], 'float32')

        xp = xpdf.merge(s, left_on='a', right_index=True)
        assert_frame_equal(rs, xp)

    def test_join_many_non_unique_index(self):
        df1 = DataFrame({"a": [1, 1], "b": [1, 1], "c": [10, 20]})
        df2 = DataFrame({"a": [1, 1], "b": [1, 2], "d": [100, 200]})
        df3 = DataFrame({"a": [1, 1], "b": [1, 2], "e": [1000, 2000]})
        idf1 = df1.set_index(["a", "b"])
        idf2 = df2.set_index(["a", "b"])
        idf3 = df3.set_index(["a", "b"])

        result = idf1.join([idf2, idf3], how='outer')

        df_partially_merged = merge(df1, df2, on=['a', 'b'], how='outer')
        expected = merge(df_partially_merged, df3, on=['a', 'b'], how='outer')

        result = result.reset_index()

        result['a'] = result['a'].astype(np.float64)
        result['b'] = result['b'].astype(np.float64)

        assert_frame_equal(result, expected.ix[:, result.columns])

        df1 = DataFrame({"a": [1, 1, 1], "b": [1, 1, 1], "c": [10, 20, 30]})
        df2 = DataFrame({"a": [1, 1, 1], "b": [1, 1, 2], "d": [100, 200, 300]})
        df3 = DataFrame(
            {"a": [1, 1, 1], "b": [1, 1, 2], "e": [1000, 2000, 3000]})
        idf1 = df1.set_index(["a", "b"])
        idf2 = df2.set_index(["a", "b"])
        idf3 = df3.set_index(["a", "b"])
        result = idf1.join([idf2, idf3], how='inner')

        df_partially_merged = merge(df1, df2, on=['a', 'b'], how='inner')
        expected = merge(df_partially_merged, df3, on=['a', 'b'], how='inner')

        result = result.reset_index()

        assert_frame_equal(result, expected.ix[:, result.columns])

    def test_merge_index_singlekey_right_vs_left(self):
        left = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'e', 'a'],
                          'v1': np.random.randn(7)})
        right = DataFrame({'v2': np.random.randn(4)},
                          index=['d', 'b', 'c', 'a'])

        merged1 = merge(left, right, left_on='key',
                        right_index=True, how='left', sort=False)
        merged2 = merge(right, left, right_on='key',
                        left_index=True, how='right', sort=False)
        assert_frame_equal(merged1, merged2.ix[:, merged1.columns])

        merged1 = merge(left, right, left_on='key',
                        right_index=True, how='left', sort=True)
        merged2 = merge(right, left, right_on='key',
                        left_index=True, how='right', sort=True)
        assert_frame_equal(merged1, merged2.ix[:, merged1.columns])

    def test_merge_index_singlekey_inner(self):
        left = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'e', 'a'],
                          'v1': np.random.randn(7)})
        right = DataFrame({'v2': np.random.randn(4)},
                          index=['d', 'b', 'c', 'a'])

        # inner join
        result = merge(left, right, left_on='key', right_index=True,
                       how='inner')
        expected = left.join(right, on='key').ix[result.index]
        assert_frame_equal(result, expected)

        result = merge(right, left, right_on='key', left_index=True,
                       how='inner')
        expected = left.join(right, on='key').ix[result.index]
        assert_frame_equal(result, expected.ix[:, result.columns])

    def test_merge_misspecified(self):
        self.assertRaises(ValueError, merge, self.left, self.right,
                          left_index=True)
        self.assertRaises(ValueError, merge, self.left, self.right,
                          right_index=True)

        self.assertRaises(ValueError, merge, self.left, self.left,
                          left_on='key', on='key')

        self.assertRaises(ValueError, merge, self.df, self.df2,
                          left_on=['key1'], right_on=['key1', 'key2'])

    def test_merge_overlap(self):
        merged = merge(self.left, self.left, on='key')
        exp_len = (self.left['key'].value_counts() ** 2).sum()
        self.assertEqual(len(merged), exp_len)
        self.assertIn('v1_x', merged)
        self.assertIn('v1_y', merged)

    def test_merge_different_column_key_names(self):
        left = DataFrame({'lkey': ['foo', 'bar', 'baz', 'foo'],
                          'value': [1, 2, 3, 4]})
        right = DataFrame({'rkey': ['foo', 'bar', 'qux', 'foo'],
                           'value': [5, 6, 7, 8]})

        merged = left.merge(right, left_on='lkey', right_on='rkey',
                            how='outer', sort=True)

        assert_almost_equal(merged['lkey'],
                            ['bar', 'baz', 'foo', 'foo', 'foo', 'foo', np.nan])
        assert_almost_equal(merged['rkey'],
                            ['bar', np.nan, 'foo', 'foo', 'foo', 'foo', 'qux'])
        assert_almost_equal(merged['value_x'], [2, 3, 1, 1, 4, 4, np.nan])
        assert_almost_equal(merged['value_y'], [6, np.nan, 5, 8, 5, 8, 7])

    def test_merge_copy(self):
        left = DataFrame({'a': 0, 'b': 1}, index=lrange(10))
        right = DataFrame({'c': 'foo', 'd': 'bar'}, index=lrange(10))

        merged = merge(left, right, left_index=True,
                       right_index=True, copy=True)

        merged['a'] = 6
        self.assertTrue((left['a'] == 0).all())

        merged['d'] = 'peekaboo'
        self.assertTrue((right['d'] == 'bar').all())

    def test_merge_nocopy(self):
        left = DataFrame({'a': 0, 'b': 1}, index=lrange(10))
        right = DataFrame({'c': 'foo', 'd': 'bar'}, index=lrange(10))

        merged = merge(left, right, left_index=True,
                       right_index=True, copy=False)

        merged['a'] = 6
        self.assertTrue((left['a'] == 6).all())

        merged['d'] = 'peekaboo'
        self.assertTrue((right['d'] == 'peekaboo').all())

    def test_join_sort(self):
        left = DataFrame({'key': ['foo', 'bar', 'baz', 'foo'],
                          'value': [1, 2, 3, 4]})
        right = DataFrame({'value2': ['a', 'b', 'c']},
                          index=['bar', 'baz', 'foo'])

        joined = left.join(right, on='key', sort=True)
        expected = DataFrame({'key': ['bar', 'baz', 'foo', 'foo'],
                              'value': [2, 3, 1, 4],
                              'value2': ['a', 'b', 'c', 'c']},
                             index=[1, 2, 0, 3])
        assert_frame_equal(joined, expected)

        # smoke test
        joined = left.join(right, on='key', sort=False)
        self.assert_numpy_array_equal(joined.index, lrange(4))

    def test_intelligently_handle_join_key(self):
        # #733, be a bit more 1337 about not returning unconsolidated DataFrame

        left = DataFrame({'key': [1, 1, 2, 2, 3],
                          'value': lrange(5)}, columns=['value', 'key'])
        right = DataFrame({'key': [1, 1, 2, 3, 4, 5],
                           'rvalue': lrange(6)})

        joined = merge(left, right, on='key', how='outer')
        expected = DataFrame({'key': [1, 1, 1, 1, 2, 2, 3, 4, 5.],
                              'value': np.array([0, 0, 1, 1, 2, 3, 4,
                                                 np.nan, np.nan]),
                              'rvalue': np.array([0, 1, 0, 1, 2, 2, 3, 4, 5])},
                             columns=['value', 'key', 'rvalue'])
        assert_frame_equal(joined, expected, check_dtype=False)

        self.assertTrue(joined._data.is_consolidated())

    def test_handle_join_key_pass_array(self):
        left = DataFrame({'key': [1, 1, 2, 2, 3],
                          'value': lrange(5)}, columns=['value', 'key'])
        right = DataFrame({'rvalue': lrange(6)})
        key = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on='key', right_on=key, how='outer')
        merged2 = merge(right, left, left_on=key, right_on='key', how='outer')

        assert_series_equal(merged['key'], merged2['key'])
        self.assertTrue(merged['key'].notnull().all())
        self.assertTrue(merged2['key'].notnull().all())

        left = DataFrame({'value': lrange(5)}, columns=['value'])
        right = DataFrame({'rvalue': lrange(6)})
        lkey = np.array([1, 1, 2, 2, 3])
        rkey = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on=lkey, right_on=rkey, how='outer')
        self.assert_numpy_array_equal(merged['key_0'],
                                      np.array([1, 1, 1, 1, 2, 2, 3, 4, 5]))

        left = DataFrame({'value': lrange(3)})
        right = DataFrame({'rvalue': lrange(6)})

        key = np.array([0, 1, 1, 2, 2, 3])
        merged = merge(left, right, left_index=True, right_on=key, how='outer')
        self.assert_numpy_array_equal(merged['key_0'], key)

    def test_mixed_type_join_with_suffix(self):
        # GH #916
        df = DataFrame(np.random.randn(20, 6),
                       columns=['a', 'b', 'c', 'd', 'e', 'f'])
        df.insert(0, 'id', 0)
        df.insert(5, 'dt', 'foo')

        grouped = df.groupby('id')
        mn = grouped.mean()
        cn = grouped.count()

        # it works!
        mn.join(cn, rsuffix='_right')

    def test_no_overlap_more_informative_error(self):
        dt = datetime.now()
        df1 = DataFrame({'x': ['a']}, index=[dt])

        df2 = DataFrame({'y': ['b', 'c']}, index=[dt, dt])
        self.assertRaises(MergeError, merge, df1, df2)

    def test_merge_non_unique_indexes(self):

        dt = datetime(2012, 5, 1)
        dt2 = datetime(2012, 5, 2)
        dt3 = datetime(2012, 5, 3)
        dt4 = datetime(2012, 5, 4)

        df1 = DataFrame({'x': ['a']}, index=[dt])
        df2 = DataFrame({'y': ['b', 'c']}, index=[dt, dt])
        _check_merge(df1, df2)

        # Not monotonic
        df1 = DataFrame({'x': ['a', 'b', 'q']}, index=[dt2, dt, dt4])
        df2 = DataFrame({'y': ['c', 'd', 'e', 'f', 'g', 'h']},
                        index=[dt3, dt3, dt2, dt2, dt, dt])
        _check_merge(df1, df2)

        df1 = DataFrame({'x': ['a', 'b']}, index=[dt, dt])
        df2 = DataFrame({'y': ['c', 'd']}, index=[dt, dt])
        _check_merge(df1, df2)

    def test_merge_non_unique_index_many_to_many(self):
        dt = datetime(2012, 5, 1)
        dt2 = datetime(2012, 5, 2)
        dt3 = datetime(2012, 5, 3)
        df1 = DataFrame({'x': ['a', 'b', 'c', 'd']},
                        index=[dt2, dt2, dt, dt])
        df2 = DataFrame({'y': ['e', 'f', 'g', ' h', 'i']},
                        index=[dt2, dt2, dt3, dt, dt])
        _check_merge(df1, df2)

    def test_left_merge_empty_dataframe(self):
        left = DataFrame({'key': [1], 'value': [2]})
        right = DataFrame({'key': []})

        result = merge(left, right, on='key', how='left')
        assert_frame_equal(result, left)

        result = merge(right, left, on='key', how='right')
        assert_frame_equal(result, left)

    def test_merge_left_empty_right_empty(self):
        # GH 10824
        left = pd.DataFrame([], columns=['a', 'b', 'c'])
        right = pd.DataFrame([], columns=['x', 'y', 'z'])

        exp_in = pd.DataFrame([], columns=['a', 'b', 'c', 'x', 'y', 'z'],
                              dtype=object)

        for kwarg in [dict(left_index=True, right_index=True),
                      dict(left_index=True, right_on='x'),
                      dict(left_on='a', right_index=True),
                      dict(left_on='a', right_on='x')]:

            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp_in)
            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp_in)
            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp_in)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp_in)

    def test_merge_left_empty_right_notempty(self):
        # GH 10824
        left = pd.DataFrame([], columns=['a', 'b', 'c'])
        right = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                             columns=['x', 'y', 'z'])

        exp_out = pd.DataFrame({'a': np.array([np.nan]*3, dtype=object),
                                'b': np.array([np.nan]*3, dtype=object),
                                'c': np.array([np.nan]*3, dtype=object),
                                'x': [1, 4, 7],
                                'y': [2, 5, 8],
                                'z': [3, 6, 9]},
                               columns=['a', 'b', 'c', 'x', 'y', 'z'])
        exp_in = exp_out[0:0] # make empty DataFrame keeping dtype

        for kwarg in [dict(left_index=True, right_index=True),
                      dict(left_index=True, right_on='x'),
                      dict(left_on='a', right_index=True),
                      dict(left_on='a', right_on='x')]:

            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp_in)
            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp_in)

            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp_out)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp_out)

    def test_merge_left_notempty_right_empty(self):
        # GH 10824
        left = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            columns=['a', 'b', 'c'])
        right = pd.DataFrame([], columns=['x', 'y', 'z'])

        exp_out = pd.DataFrame({'a': [1, 4, 7],
                                'b': [2, 5, 8],
                                'c': [3, 6, 9],
                                'x': np.array([np.nan]*3, dtype=object),
                                'y': np.array([np.nan]*3, dtype=object),
                                'z': np.array([np.nan]*3, dtype=object)},
                               columns=['a', 'b', 'c', 'x', 'y', 'z'])
        exp_in = exp_out[0:0] # make empty DataFrame keeping dtype

        for kwarg in [dict(left_index=True, right_index=True),
                      dict(left_index=True, right_on='x'),
                      dict(left_on='a', right_index=True),
                      dict(left_on='a', right_on='x')]:

            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp_in)
            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp_in)

            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp_out)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp_out)

    def test_merge_nosort(self):
        # #2098, anything to do?

        from datetime import datetime

        d = {"var1": np.random.randint(0, 10, size=10),
             "var2": np.random.randint(0, 10, size=10),
             "var3": [datetime(2012, 1, 12), datetime(2011, 2, 4),
                      datetime(
                      2010, 2, 3), datetime(2012, 1, 12),
                      datetime(
                      2011, 2, 4), datetime(2012, 4, 3),
                      datetime(
                      2012, 3, 4), datetime(2008, 5, 1),
                      datetime(2010, 2, 3), datetime(2012, 2, 3)]}
        df = DataFrame.from_dict(d)
        var3 = df.var3.unique()
        var3.sort()
        new = DataFrame.from_dict({"var3": var3,
                                   "var8": np.random.random(7)})

        result = df.merge(new, on="var3", sort=False)
        exp = merge(df, new, on='var3', sort=False)
        assert_frame_equal(result, exp)

        self.assertTrue((df.var3.unique() == result.var3.unique()).all())

    def test_merge_nan_right(self):
        df1 = DataFrame({"i1" : [0, 1], "i2" : [0, 1]})
        df2 = DataFrame({"i1" : [0], "i3" : [0]})
        result = df1.join(df2, on="i1", rsuffix="_")
        expected = DataFrame({'i1': {0: 0.0, 1: 1}, 'i2': {0: 0, 1: 1},
                              'i1_': {0: 0, 1: np.nan}, 'i3': {0: 0.0, 1: np.nan},
                               None: {0: 0, 1: 0}}).set_index(None).reset_index()[['i1', 'i2', 'i1_', 'i3']]
        assert_frame_equal(result, expected, check_dtype=False)

        df1 = DataFrame({"i1" : [0, 1], "i2" : [0.5, 1.5]})
        df2 = DataFrame({"i1" : [0], "i3" : [0.7]})
        result = df1.join(df2, rsuffix="_", on='i1')
        expected = DataFrame({'i1': {0: 0, 1: 1}, 'i1_': {0: 0.0, 1: nan},
                              'i2': {0: 0.5, 1: 1.5}, 'i3': {0: 0.69999999999999996,
                              1: nan}})[['i1', 'i2', 'i1_', 'i3']]
        assert_frame_equal(result, expected)

    def test_merge_type(self):
        class NotADataFrame(DataFrame):
            @property
            def _constructor(self):
                return NotADataFrame

        nad = NotADataFrame(self.df)
        result = nad.merge(self.df2, on='key1')

        tm.assertIsInstance(result, NotADataFrame)

    def test_append_dtype_coerce(self):

        # GH 4993
        # appending with datetime will incorrectly convert datetime64
        import datetime as dt
        from pandas import NaT

        df1 = DataFrame(index=[1,2], data=[dt.datetime(2013,1,1,0,0),
                                           dt.datetime(2013,1,2,0,0)],
                        columns=['start_time'])
        df2 = DataFrame(index=[4,5], data=[[dt.datetime(2013,1,3,0,0),
                                            dt.datetime(2013,1,3,6,10)],
                                           [dt.datetime(2013,1,4,0,0),
                                            dt.datetime(2013,1,4,7,10)]],
                        columns=['start_time','end_time'])

        expected = concat([
            Series([NaT,NaT,dt.datetime(2013,1,3,6,10),dt.datetime(2013,1,4,7,10)],name='end_time'),
            Series([dt.datetime(2013,1,1,0,0),dt.datetime(2013,1,2,0,0),dt.datetime(2013,1,3,0,0),dt.datetime(2013,1,4,0,0)],name='start_time'),
            ],axis=1)
        result = df1.append(df2,ignore_index=True)
        assert_frame_equal(result, expected)

    def test_join_append_timedeltas(self):

        import datetime as dt
        from pandas import NaT

        # timedelta64 issues with join/merge
        # GH 5695

        d = {'d': dt.datetime(2013, 11, 5, 5, 56), 't': dt.timedelta(0, 22500)}
        df = DataFrame(columns=list('dt'))
        df = df.append(d, ignore_index=True)
        result = df.append(d, ignore_index=True)
        expected = DataFrame({'d': [dt.datetime(2013, 11, 5, 5, 56),
                                    dt.datetime(2013, 11, 5, 5, 56) ],
                              't': [ dt.timedelta(0, 22500),
                                     dt.timedelta(0, 22500) ]})
        assert_frame_equal(result, expected)

        td = np.timedelta64(300000000)
        lhs = DataFrame(Series([td,td],index=["A","B"]))
        rhs = DataFrame(Series([td],index=["A"]))

        from pandas import NaT
        result = lhs.join(rhs,rsuffix='r', how="left")
        expected = DataFrame({ '0' : Series([td,td],index=list('AB')), '0r' : Series([td,NaT],index=list('AB')) })
        assert_frame_equal(result, expected)

    def test_overlapping_columns_error_message(self):
        df = DataFrame({'key': [1, 2, 3],
                        'v1': [4, 5, 6],
                        'v2': [7, 8, 9]})
        df2 = DataFrame({'key': [1, 2, 3],
                         'v1': [4, 5, 6],
                         'v2': [7, 8, 9]})

        df.columns = ['key', 'foo', 'foo']
        df2.columns = ['key', 'bar', 'bar']
        expected = DataFrame({'key': [1, 2, 3],
                         'v1': [4, 5, 6],
                         'v2': [7, 8, 9],
                         'v3': [4, 5, 6],
                         'v4': [7, 8, 9]})
        expected.columns = ['key', 'foo', 'foo', 'bar', 'bar']
        assert_frame_equal(merge(df, df2), expected)

        # #2649, #10639
        df2.columns = ['key1', 'foo', 'foo']
        self.assertRaises(ValueError, merge, df, df2)

    def test_indicator(self):
        # PR #10054. xref #7412 and closes #8790.
        df1 = pd.DataFrame({'col1':[0,1], 'col_left':['a','b'], 'col_conflict':[1,2]})
        df1_copy = df1.copy()

        df2 = pd.DataFrame({'col1':[1,2,3,4,5],'col_right':[2,2,2,2,2],
                            'col_conflict':[1,2,3,4,5]})
        df2_copy = df2.copy()

        df_result = pd.DataFrame({'col1':[0,1,2,3,4,5],
                'col_conflict_x':[1,2,np.nan,np.nan,np.nan,np.nan],
                'col_left':['a','b', np.nan,np.nan,np.nan,np.nan],
                'col_conflict_y':[np.nan,1,2,3,4,5],
                'col_right':[np.nan, 2,2,2,2,2]},
                dtype='float64')
        df_result['_merge'] = pd.Categorical(['left_only','both','right_only',
            'right_only','right_only','right_only']
            , categories=['left_only', 'right_only', 'both'])

        df_result = df_result[['col1', 'col_conflict_x', 'col_left',
                               'col_conflict_y', 'col_right', '_merge' ]]

        test = pd.merge(df1, df2, on='col1', how='outer', indicator=True)
        assert_frame_equal(test, df_result)
        test = df1.merge(df2, on='col1', how='outer', indicator=True)
        assert_frame_equal(test, df_result)

        # No side effects
        assert_frame_equal(df1, df1_copy)
        assert_frame_equal(df2, df2_copy)

        # Check with custom name
        df_result_custom_name = df_result
        df_result_custom_name = df_result_custom_name.rename(columns={'_merge':'custom_name'})

        test_custom_name = pd.merge(df1, df2, on='col1', how='outer', indicator='custom_name')
        assert_frame_equal(test_custom_name, df_result_custom_name)
        test_custom_name = df1.merge(df2, on='col1', how='outer', indicator='custom_name')
        assert_frame_equal(test_custom_name, df_result_custom_name)

        # Check only accepts strings and booleans
        with tm.assertRaises(ValueError):
            pd.merge(df1, df2, on='col1', how='outer', indicator=5)
        with tm.assertRaises(ValueError):
            df1.merge(df2, on='col1', how='outer', indicator=5)

        # Check result integrity

        test2 = pd.merge(df1, df2, on='col1', how='left', indicator=True)
        self.assertTrue((test2._merge != 'right_only').all())
        test2 = df1.merge(df2, on='col1', how='left', indicator=True)
        self.assertTrue((test2._merge != 'right_only').all())

        test3 = pd.merge(df1, df2, on='col1', how='right', indicator=True)
        self.assertTrue((test3._merge != 'left_only').all())
        test3 = df1.merge(df2, on='col1', how='right', indicator=True)
        self.assertTrue((test3._merge != 'left_only').all())

        test4 = pd.merge(df1, df2, on='col1', how='inner', indicator=True)
        self.assertTrue((test4._merge == 'both').all())
        test4 = df1.merge(df2, on='col1', how='inner', indicator=True)
        self.assertTrue((test4._merge == 'both').all())

        # Check if working name in df
        for i in ['_right_indicator', '_left_indicator', '_merge']:
            df_badcolumn = pd.DataFrame({'col1':[1,2], i:[2,2]})

            with tm.assertRaises(ValueError):
                pd.merge(df1, df_badcolumn, on='col1', how='outer', indicator=True)
            with tm.assertRaises(ValueError):
                df1.merge(df_badcolumn, on='col1', how='outer', indicator=True)

        # Check for name conflict with custom name
        df_badcolumn = pd.DataFrame({'col1':[1,2], 'custom_column_name':[2,2]})

        with tm.assertRaises(ValueError):
            pd.merge(df1, df_badcolumn, on='col1', how='outer', indicator='custom_column_name')
        with tm.assertRaises(ValueError):
            df1.merge(df_badcolumn, on='col1', how='outer', indicator='custom_column_name')

        # Merge on multiple columns
        df3 = pd.DataFrame({'col1':[0,1], 'col2':['a','b']})

        df4 = pd.DataFrame({'col1':[1,1,3], 'col2':['b','x','y']})

        hand_coded_result = pd.DataFrame({'col1':[0,1,1,3.0],
                                         'col2':['a','b','x','y']})
        hand_coded_result['_merge'] = pd.Categorical(
            ['left_only','both','right_only','right_only']
            , categories=['left_only', 'right_only', 'both'])

        test5 = pd.merge(df3, df4, on=['col1', 'col2'], how='outer', indicator=True)
        assert_frame_equal(test5, hand_coded_result)
        test5 = df3.merge(df4, on=['col1', 'col2'], how='outer', indicator=True)
        assert_frame_equal(test5, hand_coded_result)


def _check_merge(x, y):
    for how in ['inner', 'left', 'outer']:
        result = x.join(y, how=how)

        expected = merge(x.reset_index(), y.reset_index(), how=how,
                         sort=True)
        expected = expected.set_index('index')

        assert_frame_equal(result, expected, check_names=False)  # TODO check_names on merge?


class TestMergeMulti(tm.TestCase):

    def setUp(self):
        self.index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                        ['one', 'two', 'three']],
                                labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                        [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                                names=['first', 'second'])
        self.to_join = DataFrame(np.random.randn(10, 3), index=self.index,
                                 columns=['j_one', 'j_two', 'j_three'])

        # a little relevant example with NAs
        key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
                'qux', 'snap']
        key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
                'three', 'one']

        data = np.random.randn(len(key1))
        self.data = DataFrame({'key1': key1, 'key2': key2,
                               'data': data})

    def test_merge_on_multikey(self):
        joined = self.data.join(self.to_join, on=['key1', 'key2'])

        join_key = Index(lzip(self.data['key1'], self.data['key2']))
        indexer = self.to_join.index.get_indexer(join_key)
        ex_values = self.to_join.values.take(indexer, axis=0)
        ex_values[indexer == -1] = np.nan
        expected = self.data.join(DataFrame(ex_values,
                                            columns=self.to_join.columns))

        # TODO: columns aren't in the same order yet
        assert_frame_equal(joined, expected.ix[:, joined.columns])

        left = self.data.join(self.to_join, on=['key1', 'key2'], sort=True)
        right = expected.ix[:, joined.columns].sort_values(['key1', 'key2'],
                                                           kind='mergesort')
        assert_frame_equal(left, right)

    def test_left_join_multi_index(self):
        icols = ['1st', '2nd', '3rd']

        def bind_cols(df):
            iord = lambda a: 0 if a != a else ord(a)
            f = lambda ts: ts.map(iord) - ord('a')
            return f(df['1st']) + f(df['3rd'])* 1e2 + df['2nd'].fillna(0) * 1e4

        def run_asserts(left, right):
            for sort in [False, True]:
                res = left.join(right, on=icols, how='left', sort=sort)

                self.assertTrue(len(left) < len(res) + 1)
                self.assertFalse(res['4th'].isnull().any())
                self.assertFalse(res['5th'].isnull().any())

                tm.assert_series_equal(res['4th'], - res['5th'], check_names=False)
                result = bind_cols(res.iloc[:, :-2])
                tm.assert_series_equal(res['4th'], result, check_names=False)
                self.assertTrue(result.name is None)

                if sort:
                    tm.assert_frame_equal(res,
                                          res.sort_values(icols, kind='mergesort'))

                out = merge(left, right.reset_index(), on=icols,
                            sort=sort, how='left')

                res.index = np.arange(len(res))
                tm.assert_frame_equal(out, res)

        lc = list(map(chr, np.arange(ord('a'), ord('z') + 1)))
        left = DataFrame(np.random.choice(lc, (5000, 2)),
                         columns=['1st', '3rd'])
        left.insert(1, '2nd', np.random.randint(0, 1000, len(left)))

        i = np.random.permutation(len(left))
        right = left.iloc[i].copy()

        left['4th'] = bind_cols(left)
        right['5th'] = - bind_cols(right)
        right.set_index(icols, inplace=True)

        run_asserts(left, right)

        # inject some nulls
        left.loc[1::23, '1st'] = np.nan
        left.loc[2::37, '2nd'] = np.nan
        left.loc[3::43, '3rd'] = np.nan
        left['4th'] = bind_cols(left)

        i = np.random.permutation(len(left))
        right = left.iloc[i, :-1]
        right['5th'] = - bind_cols(right)
        right.set_index(icols, inplace=True)

        run_asserts(left, right)

    def test_merge_right_vs_left(self):
        # compare left vs right merge with multikey
        for sort in [False, True]:
            merged1 = self.data.merge(self.to_join, left_on=['key1', 'key2'],
                    right_index=True, how='left', sort=sort)

            merged2 = self.to_join.merge(self.data, right_on=['key1', 'key2'],
                    left_index=True, how='right', sort=sort)

            merged2 = merged2.ix[:, merged1.columns]
            assert_frame_equal(merged1, merged2)

    def test_compress_group_combinations(self):

        # ~ 40000000 possible unique groups
        key1 = tm.rands_array(10, 10000)
        key1 = np.tile(key1, 2)
        key2 = key1[::-1]

        df = DataFrame({'key1': key1, 'key2': key2,
                        'value1': np.random.randn(20000)})

        df2 = DataFrame({'key1': key1[::2], 'key2': key2[::2],
                         'value2': np.random.randn(10000)})

        # just to hit the label compression code path
        merged = merge(df, df2, how='outer')

    def test_left_join_index_preserve_order(self):

        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'v': np.array(np.arange(24),dtype=np.int64) })

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'),'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'),'v2'] = 7

        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result.sort_values(['k1', 'k2'], kind='mergesort'),
                              left.join(right, on=['k1', 'k2'], sort=True))

        # test join with multi dtypes blocks
        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'k3' : np.array([0, 1, 2]*8, dtype=np.float32),
                          'v': np.array(np.arange(24),dtype=np.int32) })

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'),'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'),'v2'] = 7

        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result.sort_values(['k1', 'k2'], kind='mergesort'),
                              left.join(right, on=['k1', 'k2'], sort=True))

        # do a right join for an extra test
        joined = merge(right, left, left_index=True,
                       right_on=['k1', 'k2'], how='right')
        tm.assert_frame_equal(joined.ix[:, expected.columns], expected)

    def test_left_join_index_multi_match_multiindex(self):
        left = DataFrame([
            ['X', 'Y', 'C', 'a'],
            ['W', 'Y', 'C', 'e'],
            ['V', 'Q', 'A', 'h'],
            ['V', 'R', 'D', 'i'],
            ['X', 'Y', 'D', 'b'],
            ['X', 'Y', 'A', 'c'],
            ['W', 'Q', 'B', 'f'],
            ['W', 'R', 'C', 'g'],
            ['V', 'Y', 'C', 'j'],
            ['X', 'Y', 'B', 'd']],
            columns=['cola', 'colb', 'colc', 'tag'],
            index=[3, 2, 0, 1, 7, 6, 4, 5, 9, 8])

        right = DataFrame([
            ['W', 'R', 'C',  0],
            ['W', 'Q', 'B',  3],
            ['W', 'Q', 'B',  8],
            ['X', 'Y', 'A',  1],
            ['X', 'Y', 'A',  4],
            ['X', 'Y', 'B',  5],
            ['X', 'Y', 'C',  6],
            ['X', 'Y', 'C',  9],
            ['X', 'Q', 'C', -6],
            ['X', 'R', 'C', -9],
            ['V', 'Y', 'C',  7],
            ['V', 'R', 'D',  2],
            ['V', 'R', 'D', -1],
            ['V', 'Q', 'A', -3]],
            columns=['col1', 'col2', 'col3', 'val'])

        right.set_index(['col1', 'col2', 'col3'], inplace=True)
        result = left.join(right, on=['cola', 'colb', 'colc'], how='left')

        expected = DataFrame([
            ['X', 'Y', 'C', 'a',   6],
            ['X', 'Y', 'C', 'a',   9],
            ['W', 'Y', 'C', 'e', nan],
            ['V', 'Q', 'A', 'h',  -3],
            ['V', 'R', 'D', 'i',   2],
            ['V', 'R', 'D', 'i',  -1],
            ['X', 'Y', 'D', 'b', nan],
            ['X', 'Y', 'A', 'c',   1],
            ['X', 'Y', 'A', 'c',   4],
            ['W', 'Q', 'B', 'f',   3],
            ['W', 'Q', 'B', 'f',   8],
            ['W', 'R', 'C', 'g',   0],
            ['V', 'Y', 'C', 'j',   7],
            ['X', 'Y', 'B', 'd',   5]],
            columns=['cola', 'colb', 'colc', 'tag', 'val'],
            index=[3, 3, 2, 0, 1, 1, 7, 6, 6, 4, 4, 5, 9, 8])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on=['cola', 'colb', 'colc'],
                           how='left', sort=True)

        tm.assert_frame_equal(result,
                expected.sort_values(['cola', 'colb', 'colc'], kind='mergesort'))

        # GH7331 - maintain left frame order in left merge
        right.reset_index(inplace=True)
        right.columns = left.columns[:3].tolist() + right.columns[-1:].tolist()
        result = merge(left, right, how='left', on=left.columns[:-1].tolist())
        expected.index = np.arange(len(expected))
        tm.assert_frame_equal(result, expected)

    def test_left_join_index_multi_match(self):
        left = DataFrame([
            ['c', 0],
            ['b', 1],
            ['a', 2],
            ['b', 3]],
            columns=['tag', 'val'],
            index=[2, 0, 1, 3])

        right = DataFrame([
            ['a', 'v'],
            ['c', 'w'],
            ['c', 'x'],
            ['d', 'y'],
            ['a', 'z'],
            ['c', 'r'],
            ['e', 'q'],
            ['c', 's']],
            columns=['tag', 'char'])

        right.set_index('tag', inplace=True)
        result = left.join(right, on='tag', how='left')

        expected = DataFrame([
            ['c', 0, 'w'],
            ['c', 0, 'x'],
            ['c', 0, 'r'],
            ['c', 0, 's'],
            ['b', 1, nan],
            ['a', 2, 'v'],
            ['a', 2, 'z'],
            ['b', 3, nan]],
            columns=['tag', 'val', 'char'],
            index=[2, 2, 2, 2, 0, 1, 1, 3])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on='tag', how='left', sort=True)
        tm.assert_frame_equal(result, expected.sort_values('tag', kind='mergesort'))

        # GH7331 - maintain left frame order in left merge
        result = merge(left, right.reset_index(), how='left', on='tag')
        expected.index = np.arange(len(expected))
        tm.assert_frame_equal(result, expected)

    def test_join_multi_dtypes(self):

        # test with multi dtypes in the join index
        def _test(dtype1,dtype2):
            left = DataFrame({'k1': np.array([0, 1, 2] * 8, dtype=dtype1),
                              'k2': ['foo', 'bar'] * 12,
                              'v': np.array(np.arange(24),dtype=np.int64) })

            index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
            right = DataFrame({'v2': np.array([5, 7], dtype=dtype2)}, index=index)

            result = left.join(right, on=['k1', 'k2'])

            expected = left.copy()

            if dtype2.kind == 'i':
                dtype2 = np.dtype('float64')
            expected['v2'] = np.array(np.nan,dtype=dtype2)
            expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'),'v2'] = 5
            expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'),'v2'] = 7

            tm.assert_frame_equal(result, expected)

            result = left.join(right, on=['k1', 'k2'], sort=True)
            expected.sort_values(['k1', 'k2'], kind='mergesort', inplace=True)
            tm.assert_frame_equal(result, expected)

        for d1 in [np.int64,np.int32,np.int16,np.int8,np.uint8]:
            for d2 in [np.int64,np.float64,np.float32,np.float16]:
                _test(np.dtype(d1),np.dtype(d2))

    def test_left_merge_na_buglet(self):
        left = DataFrame({'id': list('abcde'), 'v1': randn(5),
                          'v2': randn(5), 'dummy': list('abcde'),
                          'v3': randn(5)},
                         columns=['id', 'v1', 'v2', 'dummy', 'v3'])
        right = DataFrame({'id': ['a', 'b', np.nan, np.nan, np.nan],
                           'sv3': [1.234, 5.678, np.nan, np.nan, np.nan]})

        merged = merge(left, right, on='id', how='left')

        rdf = right.drop(['id'], axis=1)
        expected = left.join(rdf)
        tm.assert_frame_equal(merged, expected)

    def test_merge_na_keys(self):
        data = [[1950, "A", 1.5],
                [1950, "B", 1.5],
                [1955, "B", 1.5],
                [1960, "B", np.nan],
                [1970, "B", 4.],
                [1950, "C", 4.],
                [1960, "C", np.nan],
                [1965, "C", 3.],
                [1970, "C", 4.]]

        frame = DataFrame(data, columns=["year", "panel", "data"])

        other_data = [[1960, 'A', np.nan],
                      [1970, 'A', np.nan],
                      [1955, 'A', np.nan],
                      [1965, 'A', np.nan],
                      [1965, 'B', np.nan],
                      [1955, 'C', np.nan]]
        other = DataFrame(other_data, columns=['year', 'panel', 'data'])

        result = frame.merge(other, how='outer')

        expected = frame.fillna(-999).merge(other.fillna(-999), how='outer')
        expected = expected.replace(-999, np.nan)

        tm.assert_frame_equal(result, expected)

    def test_int64_overflow_issues(self):
        from itertools import product
        from collections import defaultdict
        from pandas.core.groupby import _int64_overflow_possible

        # #2690, combinatorial explosion
        df1 = DataFrame(np.random.randn(1000, 7),
                        columns=list('ABCDEF') + ['G1'])
        df2 = DataFrame(np.random.randn(1000, 7),
                        columns=list('ABCDEF') + ['G2'])

        # it works!
        result = merge(df1, df2, how='outer')
        self.assertTrue(len(result) == 2000)

        low, high, n = -1 << 10, 1 << 10, 1 << 20
        left = DataFrame(np.random.randint(low, high, (n, 7)),
                         columns=list('ABCDEFG'))
        left['left'] = left.sum(axis=1)

        # one-2-one match
        i = np.random.permutation(len(left))
        right = left.iloc[i].copy()
        right.columns = right.columns[:-1].tolist() + ['right']
        right.index = np.arange(len(right))
        right['right'] *= -1

        out = merge(left, right, how='outer')
        self.assertEqual(len(out), len(left))
        assert_series_equal(out['left'], - out['right'], check_names=False)
        result = out.iloc[:, :-2].sum(axis=1)
        assert_series_equal(out['left'], result, check_names=False)
        self.assertTrue(result.name is None)

        out.sort_values(out.columns.tolist(), inplace=True)
        out.index = np.arange(len(out))
        for how in ['left', 'right', 'outer', 'inner']:
            assert_frame_equal(out, merge(left, right, how=how, sort=True))

        # check that left merge w/ sort=False maintains left frame order
        out = merge(left, right, how='left', sort=False)
        assert_frame_equal(left, out[left.columns.tolist()])

        out = merge(right, left, how='left', sort=False)
        assert_frame_equal(right, out[right.columns.tolist()])

        # one-2-many/none match
        n = 1 << 11
        left = DataFrame(np.random.randint(low, high, (n, 7)).astype('int64'),
                         columns=list('ABCDEFG'))

        # confirm that this is checking what it is supposed to check
        shape = left.apply(pd.Series.nunique).values
        self.assertTrue(_int64_overflow_possible(shape))

        # add duplicates to left frame
        left = pd.concat([left, left], ignore_index=True)

        right = DataFrame(np.random.randint(low, high, (n // 2, 7)).astype('int64'),
                          columns=list('ABCDEFG'))

        # add duplicates & overlap with left to the right frame
        i = np.random.choice(len(left), n)
        right = pd.concat([right, right, left.iloc[i]], ignore_index=True)

        left['left'] = np.random.randn(len(left))
        right['right'] = np.random.randn(len(right))

        # shuffle left & right frames
        i = np.random.permutation(len(left))
        left = left.iloc[i].copy()
        left.index = np.arange(len(left))

        i = np.random.permutation(len(right))
        right = right.iloc[i].copy()
        right.index = np.arange(len(right))

        # manually compute outer merge
        ldict, rdict = defaultdict(list), defaultdict(list)

        for idx, row in left.set_index(list('ABCDEFG')).iterrows():
            ldict[idx].append(row['left'])

        for idx, row in right.set_index(list('ABCDEFG')).iterrows():
            rdict[idx].append(row['right'])

        vals = []
        for k, lval in ldict.items():
            rval = rdict.get(k, [np.nan])
            for lv, rv in product(lval, rval):
                vals.append(k + tuple([lv, rv]))

        for k, rval in rdict.items():
            if k not in ldict:
                for rv in rval:
                    vals.append(k + tuple([np.nan, rv]))

        def align(df):
            df = df.sort_values(df.columns.tolist())
            df.index = np.arange(len(df))
            return df

        def verify_order(df):
            kcols = list('ABCDEFG')
            assert_frame_equal(df[kcols].copy(),
                               df[kcols].sort_values(kcols, kind='mergesort'))

        out = DataFrame(vals, columns=list('ABCDEFG') + ['left', 'right'])
        out = align(out)

        jmask = {'left': out['left'].notnull(),
                 'right': out['right'].notnull(),
                 'inner': out['left'].notnull() & out['right'].notnull(),
                 'outer': np.ones(len(out), dtype='bool')}

        for how in 'left', 'right', 'outer', 'inner':
            mask = jmask[how]
            frame = align(out[mask].copy())
            self.assertTrue(mask.all() ^ mask.any() or how == 'outer')

            for sort in [False, True]:
                res = merge(left, right, how=how, sort=sort)
                if sort:
                    verify_order(res)

                # as in GH9092 dtypes break with outer/right join
                assert_frame_equal(frame, align(res),
                                   check_dtype=how not in ('right', 'outer'))


    def test_join_multi_levels(self):

        # GH 3662
        # merge multi-levels

        household = DataFrame(dict(household_id = [1,2,3],
                                   male = [0,1,0],
                                   wealth = [196087.3,316478.7,294750]),
                              columns = ['household_id','male','wealth']).set_index('household_id')
        portfolio = DataFrame(dict(household_id = [1,2,2,3,3,3,4],
                                   asset_id = ["nl0000301109","nl0000289783","gb00b03mlx29","gb00b03mlx29","lu0197800237","nl0000289965",np.nan],
                                   name = ["ABN Amro","Robeco","Royal Dutch Shell","Royal Dutch Shell","AAB Eastern Europe Equity Fund","Postbank BioTech Fonds",np.nan],
                                   share = [1.0,0.4,0.6,0.15,0.6,0.25,1.0]),
                              columns = ['household_id','asset_id','name','share']).set_index(['household_id','asset_id'])
        result = household.join(portfolio, how='inner')
        expected = DataFrame(dict(male = [0,1,1,0,0,0],
                                  wealth = [ 196087.3, 316478.7, 316478.7, 294750.0, 294750.0, 294750.0 ],
                                  name = ['ABN Amro','Robeco','Royal Dutch Shell','Royal Dutch Shell','AAB Eastern Europe Equity Fund','Postbank BioTech Fonds'],
                                  share = [1.00,0.40,0.60,0.15,0.60,0.25],
                                  household_id = [1,2,2,3,3,3],
                                  asset_id = ['nl0000301109','nl0000289783','gb00b03mlx29','gb00b03mlx29','lu0197800237','nl0000289965']),
                             ).set_index(['household_id','asset_id']).reindex(columns=['male','wealth','name','share'])
        assert_frame_equal(result,expected)

        assert_frame_equal(result,expected)

        # equivalency
        result2 = merge(household.reset_index(),portfolio.reset_index(),on=['household_id'],how='inner').set_index(['household_id','asset_id'])
        assert_frame_equal(result2,expected)

        result = household.join(portfolio, how='outer')
        expected = concat([expected,DataFrame(dict(share = [1.00]),
                                              index=MultiIndex.from_tuples([(4,np.nan)],
                                                                           names=['household_id','asset_id']))],
                          axis=0).reindex(columns=expected.columns)
        assert_frame_equal(result,expected)

        # invalid cases
        household.index.name = 'foo'
        def f():
            household.join(portfolio, how='inner')
        self.assertRaises(ValueError, f)

        portfolio2 = portfolio.copy()
        portfolio2.index.set_names(['household_id','foo'])
        def f():
            portfolio2.join(portfolio, how='inner')
        self.assertRaises(ValueError, f)

    def test_join_multi_levels2(self):

        # some more advanced merges
        # GH6360
        household = DataFrame(dict(household_id = [1,2,2,3,3,3,4],
                                   asset_id = ["nl0000301109","nl0000301109","gb00b03mlx29","gb00b03mlx29","lu0197800237","nl0000289965",np.nan],
                                   share = [1.0,0.4,0.6,0.15,0.6,0.25,1.0]),
                              columns = ['household_id','asset_id','share']).set_index(['household_id','asset_id'])

        log_return = DataFrame(dict(
            asset_id = ["gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "lu0197800237", "lu0197800237"],
            t = [233, 234, 235, 180, 181],
            log_return = [.09604978, -.06524096, .03532373, .03025441, .036997]
                )).set_index(["asset_id","t"])

        expected = DataFrame(dict(
            household_id = [2, 2, 2, 3, 3, 3, 3, 3],
            asset_id = ["gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "lu0197800237", "lu0197800237"],
            t = [233, 234, 235, 233, 234, 235, 180, 181],
            share = [0.6, 0.6, 0.6, 0.15, 0.15, 0.15, 0.6, 0.6],
            log_return = [.09604978, -.06524096, .03532373, .09604978, -.06524096, .03532373, .03025441, .036997]
            )).set_index(["household_id", "asset_id", "t"]).reindex(columns=['share','log_return'])

        def f():
            household.join(log_return, how='inner')
        self.assertRaises(NotImplementedError, f)

        # this is the equivalency
        result = merge(household.reset_index(),log_return.reset_index(),on=['asset_id'],how='inner').set_index(['household_id','asset_id','t'])
        assert_frame_equal(result,expected)

        expected = DataFrame(dict(
            household_id = [1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4],
            asset_id = ["nl0000301109", "nl0000289783", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29", "lu0197800237", "lu0197800237", "nl0000289965", None],
            t = [None, None, 233, 234, 235, 233, 234, 235, 180, 181, None, None],
            share = [1.0, 0.4, 0.6, 0.6, 0.6, 0.15, 0.15, 0.15, 0.6, 0.6, 0.25, 1.0],
            log_return = [None, None, .09604978, -.06524096, .03532373, .09604978, -.06524096, .03532373, .03025441, .036997, None, None]
            )).set_index(["household_id", "asset_id", "t"])

        def f():
            household.join(log_return, how='outer')
        self.assertRaises(NotImplementedError, f)

def _check_join(left, right, result, join_col, how='left',
                lsuffix='_x', rsuffix='_y'):

    # some smoke tests
    for c in join_col:
        assert(result[c].notnull().all())

    left_grouped = left.groupby(join_col)
    right_grouped = right.groupby(join_col)

    for group_key, group in result.groupby(join_col):
        l_joined = _restrict_to_columns(group, left.columns, lsuffix)
        r_joined = _restrict_to_columns(group, right.columns, rsuffix)

        try:
            lgroup = left_grouped.get_group(group_key)
        except KeyError:
            if how in ('left', 'inner'):
                raise AssertionError('key %s should not have been in the join'
                                     % str(group_key))

            _assert_all_na(l_joined, left.columns, join_col)
        else:
            _assert_same_contents(l_joined, lgroup)

        try:
            rgroup = right_grouped.get_group(group_key)
        except KeyError:
            if how in ('right', 'inner'):
                raise AssertionError('key %s should not have been in the join'
                                     % str(group_key))

            _assert_all_na(r_joined, right.columns, join_col)
        else:
            _assert_same_contents(r_joined, rgroup)


def _restrict_to_columns(group, columns, suffix):
    found = [c for c in group.columns
             if c in columns or c.replace(suffix, '') in columns]

     # filter
    group = group.ix[:, found]

    # get rid of suffixes, if any
    group = group.rename(columns=lambda x: x.replace(suffix, ''))

    # put in the right order...
    group = group.ix[:, columns]

    return group


def _assert_same_contents(join_chunk, source):
    NA_SENTINEL = -1234567  # drop_duplicates not so NA-friendly...

    jvalues = join_chunk.fillna(NA_SENTINEL).drop_duplicates().values
    svalues = source.fillna(NA_SENTINEL).drop_duplicates().values

    rows = set(tuple(row) for row in jvalues)
    assert(len(rows) == len(source))
    assert(all(tuple(row) in rows for row in svalues))


def _assert_all_na(join_chunk, source_columns, join_col):
    for c in source_columns:
        if c in join_col:
            continue
        assert(join_chunk[c].isnull().all())


def _join_by_hand(a, b, how='left'):
    join_index = a.index.join(b.index, how=how)

    a_re = a.reindex(join_index)
    b_re = b.reindex(join_index)

    result_columns = a.columns.append(b.columns)

    for col, s in compat.iteritems(b_re):
        a_re[col] = s
    return a_re.reindex(columns=result_columns)


class TestConcatenate(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.frame = DataFrame(tm.getSeriesData())
        self.mixed_frame = self.frame.copy()
        self.mixed_frame['foo'] = 'bar'

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
        expected = DataFrame({'a': {'x': 1, 'y': 2, 'z': 5}, 'b': {'x': 3, 'y': 4, 'z': 6}, 'c' : {'z' : 7}})
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

    def test_join_many(self):
        df = DataFrame(np.random.randn(10, 6), columns=list('abcdef'))
        df_list = [df[['a', 'b']], df[['c', 'd']], df[['e', 'f']]]

        joined = df_list[0].join(df_list[1:])
        tm.assert_frame_equal(joined, df)

        df_list = [df[['a', 'b']][:-2],
                   df[['c', 'd']][2:], df[['e', 'f']][1:9]]

        def _check_diff_index(df_list, result, exp_index):
            reindexed = [x.reindex(exp_index) for x in df_list]
            expected = reindexed[0].join(reindexed[1:])
            tm.assert_frame_equal(result, expected)

        # different join types
        joined = df_list[0].join(df_list[1:], how='outer')
        _check_diff_index(df_list, joined, df.index)

        joined = df_list[0].join(df_list[1:])
        _check_diff_index(df_list, joined, df_list[0].index)

        joined = df_list[0].join(df_list[1:], how='inner')
        _check_diff_index(df_list, joined, df.index[2:8])

        self.assertRaises(ValueError, df_list[0].join, df_list[1:], on='a')

    def test_join_many_mixed(self):
        df = DataFrame(np.random.randn(8, 4), columns=['A', 'B', 'C', 'D'])
        df['key'] = ['foo', 'bar'] * 4
        df1 = df.ix[:, ['A', 'B']]
        df2 = df.ix[:, ['C', 'D']]
        df3 = df.ix[:, ['key']]

        result = df1.join([df2, df3])
        assert_frame_equal(result, df)

    def test_append_missing_column_proper_upcast(self):
        df1 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='i8')})
        df2 = DataFrame({'B': np.array([True, False, True, False],
                                       dtype=bool)})

        appended = df1.append(df2, ignore_index=True)
        self.assertEqual(appended['A'].dtype, 'f8')
        self.assertEqual(appended['B'].dtype, 'O')

    def test_concat_copy(self):

        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randint(0,10,size=4).reshape(4,1))
        df3 = DataFrame({5 : 'foo'},index=range(4))

        # these are actual copies
        result = concat([df,df2,df3],axis=1,copy=True)
        for b in result._data.blocks:
            self.assertIsNone(b.values.base)

        # these are the same
        result = concat([df,df2,df3],axis=1,copy=False)
        for b in result._data.blocks:
            if b.is_float:
                self.assertTrue(b.values.base is df._data.blocks[0].values.base)
            elif b.is_integer:
                self.assertTrue(b.values.base is df2._data.blocks[0].values.base)
            elif b.is_object:
                self.assertIsNotNone(b.values.base)

        # float block was consolidated
        df4 = DataFrame(np.random.randn(4,1))
        result = concat([df,df2,df3,df4],axis=1,copy=False)
        for b in result._data.blocks:
            if b.is_float:
                self.assertIsNone(b.values.base)
            elif b.is_integer:
                self.assertTrue(b.values.base is df2._data.blocks[0].values.base)
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

        self.assert_numpy_array_equal(result.columns.levels[0], level)
        self.assertEqual(result.columns.names[0], 'group_key')

    def test_concat_dataframe_keys_bug(self):
        t1 = DataFrame({'value': Series([1, 2, 3],
                       index=Index(['a', 'b', 'c'], name='id'))})
        t2 = DataFrame({'value': Series([7, 8],
                       index=Index(['a', 'b'], name='id'))})

        # it works
        result = concat([t1, t2], axis=1, keys=['t1', 't2'])
        self.assertEqual(list(result.columns), [('t1', 'value'),
                                                ('t2', 'value')])

    def test_concat_series_partial_columns_names(self):
        # GH10698
        foo = pd.Series([1,2], name='foo')
        bar = pd.Series([1,2])
        baz = pd.Series([4,5])

        result = pd.concat([foo, bar, baz], axis=1)
        expected = DataFrame({'foo' : [1,2], 0 : [1,2], 1 : [4,5]}, columns=['foo',0,1])
        tm.assert_frame_equal(result, expected)

        result = pd.concat([foo, bar, baz], axis=1, keys=['red','blue','yellow'])
        expected = DataFrame({'red' : [1,2], 'blue' : [1,2], 'yellow' : [4,5]}, columns=['red','blue','yellow'])
        tm.assert_frame_equal(result, expected)

        result = pd.concat([foo, bar, baz], axis=1, ignore_index=True)
        expected = DataFrame({0 : [1,2], 1 : [1,2], 2 : [4,5]})
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
        df['dt'] = df['dt'].apply(lambda d: pd.Timestamp(d, tz='US/Pacific'))
        df = df.set_index(['dt', 'b'])

        exp_idx1 = pd.DatetimeIndex(['2014-01-01', '2014-01-02', '2014-01-03'] * 2,
                                    tz='US/Pacific', name='dt')
        exp_idx2 = Index(['A', 'B', 'C'] * 2, name='b')
        exp_idx = pd.MultiIndex.from_arrays([exp_idx1, exp_idx2])
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
        self.assert_numpy_array_equal(result.index.levels[0], ['baz', 'foo'])

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
        df = DataFrame(np.random.randint(0,10,size=40).reshape(10,4),columns=['A','A','C','C'])

        result = concat([df,df],axis=1)
        assert_frame_equal(result.iloc[:,:4],df)
        assert_frame_equal(result.iloc[:,4:],df)

        result = concat([df,df],axis=0)
        assert_frame_equal(result.iloc[:10],df)
        assert_frame_equal(result.iloc[10:],df)

        # multi dtypes
        df = concat([DataFrame(np.random.randn(10,4),columns=['A','A','B','B']),
                     DataFrame(np.random.randint(0,10,size=20).reshape(10,2),columns=['A','C'])],
                    axis=1)

        result = concat([df,df],axis=1)
        assert_frame_equal(result.iloc[:,:6],df)
        assert_frame_equal(result.iloc[:,6:],df)

        result = concat([df,df],axis=0)
        assert_frame_equal(result.iloc[:10],df)
        assert_frame_equal(result.iloc[10:],df)

        # append
        result = df.iloc[0:8,:].append(df.iloc[8:])
        assert_frame_equal(result, df)

        result = df.iloc[0:8,:].append(df.iloc[8:9]).append(df.iloc[9:10])
        assert_frame_equal(result, df)

        expected = concat([df,df],axis=0)
        result = df.append(df)
        assert_frame_equal(result, expected)

    def test_with_mixed_tuples(self):
        # 10697
        # columns have mixed tuples, so handle properly
        df1 = DataFrame({ u'A' : 'foo', (u'B',1) : 'bar' },index=range(2))
        df2 = DataFrame({ u'B' : 'foo', (u'B',1) : 'bar' },index=range(2))
        result = concat([df1,df2])

    def test_join_dups(self):

        # joining dups
        df = concat([DataFrame(np.random.randn(10,4),columns=['A','A','B','B']),
                     DataFrame(np.random.randint(0,10,size=20).reshape(10,2),columns=['A','C'])],
                    axis=1)

        expected = concat([df,df],axis=1)
        result = df.join(df,rsuffix='_2')
        result.columns = expected.columns
        assert_frame_equal(result, expected)

        # GH 4975, invalid join on dups
        w = DataFrame(np.random.randn(4,2), columns=["x", "y"])
        x = DataFrame(np.random.randn(4,2), columns=["x", "y"])
        y = DataFrame(np.random.randn(4,2), columns=["x", "y"])
        z = DataFrame(np.random.randn(4,2), columns=["x", "y"])

        dta = x.merge(y, left_index=True, right_index=True).merge(z, left_index=True, right_index=True, how="outer")
        dta = dta.merge(w, left_index=True, right_index=True)
        expected = concat([x,y,z,w],axis=1)
        expected.columns=['x_x','y_x','x_y','y_y','x_x','y_x','x_y','y_y']
        assert_frame_equal(dta,expected)

    def test_handle_empty_objects(self):
        df = DataFrame(np.random.randn(10, 4), columns=list('abcd'))

        baz = df[:5].copy()
        baz['foo'] = 'bar'
        empty = df[5:5]

        frames = [baz, empty, empty, df[5:]]
        concatted = concat(frames, axis=0)

        expected = df.ix[:, ['a', 'b', 'c', 'd', 'foo']]
        expected['foo'] = expected['foo'].astype('O')
        expected.loc[0:4,'foo'] = 'bar'

        tm.assert_frame_equal(concatted, expected)

        # empty as first element with time series
        # GH3259
        df = DataFrame(dict(A = range(10000)),index=date_range('20130101',periods=10000,freq='s'))
        empty = DataFrame()
        result = concat([df,empty],axis=1)
        assert_frame_equal(result, df)
        result = concat([empty,df],axis=1)
        assert_frame_equal(result, df)

        result = concat([df,empty])
        assert_frame_equal(result, df)
        result = concat([empty,df])
        assert_frame_equal(result, df)

    def test_concat_mixed_objs(self):

        # concat mixed series/frames
        # G2385

        # axis 1
        index=date_range('01-Jan-2013', periods=10, freq='H')
        arr = np.arange(10, dtype='int64')
        s1 = Series(arr, index=index)
        s2 = Series(arr, index=index)
        df = DataFrame(arr.reshape(-1,1), index=index)

        expected = DataFrame(np.repeat(arr,2).reshape(-1,2), index=index, columns = [0, 0])
        result = concat([df,df], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr,2).reshape(-1,2), index=index, columns = [0, 1])
        result = concat([s1,s2], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr,3).reshape(-1,3), index=index, columns = [0, 1, 2])
        result = concat([s1,s2,s1], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr,5).reshape(-1,5), index=index, columns = [0, 0, 1, 2, 3])
        result = concat([s1,df,s2,s2,s1], axis=1)
        assert_frame_equal(result, expected)

        # with names
        s1.name = 'foo'
        expected = DataFrame(np.repeat(arr,3).reshape(-1,3), index=index, columns = ['foo', 0, 0])
        result = concat([s1,df,s2], axis=1)
        assert_frame_equal(result, expected)

        s2.name = 'bar'
        expected = DataFrame(np.repeat(arr,3).reshape(-1,3), index=index, columns = ['foo', 0, 'bar'])
        result = concat([s1,df,s2], axis=1)
        assert_frame_equal(result, expected)

        # ignore index
        expected = DataFrame(np.repeat(arr,3).reshape(-1,3), index=index, columns = [0, 1, 2])
        result = concat([s1,df,s2], axis=1, ignore_index=True)
        assert_frame_equal(result, expected)

        # axis 0
        expected = DataFrame(np.tile(arr,3).reshape(-1,1), index=index.tolist() * 3, columns = [0])
        result = concat([s1,df,s2])
        assert_frame_equal(result, expected)

        expected = DataFrame(np.tile(arr,3).reshape(-1,1), columns = [0])
        result = concat([s1,df,s2], ignore_index=True)
        assert_frame_equal(result, expected)

        # invalid concatente of mixed dims
        panel = tm.makePanel()
        self.assertRaises(ValueError, lambda : concat([panel,s1],axis=1))

    def test_panel_join(self):
        panel = tm.makePanel()
        tm.add_nans(panel)

        p1 = panel.ix[:2, :10, :3]
        p2 = panel.ix[2:, 5:, 2:]

        # left join
        result = p1.join(p2)
        expected = p1.copy()
        expected['ItemC'] = p2['ItemC']
        tm.assert_panel_equal(result, expected)

        # right join
        result = p1.join(p2, how='right')
        expected = p2.copy()
        expected['ItemA'] = p1['ItemA']
        expected['ItemB'] = p1['ItemB']
        expected = expected.reindex(items=['ItemA', 'ItemB', 'ItemC'])
        tm.assert_panel_equal(result, expected)

        # inner join
        result = p1.join(p2, how='inner')
        expected = panel.ix[:, 5:10, 2:3]
        tm.assert_panel_equal(result, expected)

        # outer join
        result = p1.join(p2, how='outer')
        expected = p1.reindex(major=panel.major_axis,
                              minor=panel.minor_axis)
        expected = expected.join(p2.reindex(major=panel.major_axis,
                                            minor=panel.minor_axis))
        tm.assert_panel_equal(result, expected)

    def test_panel_join_overlap(self):
        panel = tm.makePanel()
        tm.add_nans(panel)

        p1 = panel.ix[['ItemA', 'ItemB', 'ItemC']]
        p2 = panel.ix[['ItemB', 'ItemC']]

        # Expected index is
        #
        # ItemA, ItemB_p1, ItemC_p1, ItemB_p2, ItemC_p2
        joined = p1.join(p2, lsuffix='_p1', rsuffix='_p2')
        p1_suf = p1.ix[['ItemB', 'ItemC']].add_suffix('_p1')
        p2_suf = p2.ix[['ItemB', 'ItemC']].add_suffix('_p2')
        no_overlap = panel.ix[['ItemA']]
        expected = no_overlap.join(p1_suf.join(p2_suf))
        tm.assert_panel_equal(joined, expected)

    def test_panel_join_many(self):
        tm.K = 10
        panel = tm.makePanel()
        tm.K = 4

        panels = [panel.ix[:2], panel.ix[2:6], panel.ix[6:]]

        joined = panels[0].join(panels[1:])
        tm.assert_panel_equal(joined, panel)

        panels = [panel.ix[:2, :-5], panel.ix[2:6, 2:], panel.ix[6:, 5:-7]]

        data_dict = {}
        for p in panels:
            data_dict.update(compat.iteritems(p))

        joined = panels[0].join(panels[1:], how='inner')
        expected = Panel.from_dict(data_dict, intersect=True)
        tm.assert_panel_equal(joined, expected)

        joined = panels[0].join(panels[1:], how='outer')
        expected = Panel.from_dict(data_dict, intersect=False)
        tm.assert_panel_equal(joined, expected)

        # edge cases
        self.assertRaises(ValueError, panels[0].join, panels[1:],
                          how='outer', lsuffix='foo', rsuffix='bar')
        self.assertRaises(ValueError, panels[0].join, panels[1:],
                          how='right')

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
        self.assertTrue(np.array_equal(result.columns, Index(['A', 0], dtype='object')))

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

        rng = to_timedelta(np.arange(10),unit='s')

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

    def test_concat_bug_1719(self):
        ts1 = tm.makeTimeSeries()
        ts2 = tm.makeTimeSeries()[::2]

        ## to join with union
        ## these two are of different length!
        left = concat([ts1, ts2], join='outer', axis=1)
        right = concat([ts2, ts1], join='outer', axis=1)

        self.assertEqual(len(left), len(right))

    def test_concat_bug_2972(self):
        ts0 = Series(np.zeros(5))
        ts1 = Series(np.ones(5))
        ts0.name = ts1.name = 'same name'
        result = concat([ts0, ts1], axis=1)

        expected = DataFrame({0: ts0, 1: ts1})
        expected.columns=['same name', 'same name']
        assert_frame_equal(result, expected)

    def test_concat_bug_3602(self):

        # GH 3602, duplicate columns
        df1 = DataFrame({'firmNo' : [0,0,0,0], 'stringvar' : ['rrr', 'rrr', 'rrr', 'rrr'], 'prc' : [6,6,6,6] })
        df2 = DataFrame({'misc' : [1,2,3,4], 'prc' : [6,6,6,6], 'C' : [9,10,11,12]})
        expected = DataFrame([[0,6,'rrr',9,1,6],
                              [0,6,'rrr',10,2,6],
                              [0,6,'rrr',11,3,6],
                              [0,6,'rrr',12,4,6]])
        expected.columns = ['firmNo','prc','stringvar','C','misc','prc']

        result = concat([df1,df2],axis=1)
        assert_frame_equal(result,expected)

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
        assert_frame_equal(pd.concat((df1, df2), ignore_index=True), expected)
        assert_frame_equal(pd.concat([df1, df2], ignore_index=True), expected)
        assert_frame_equal(pd.concat((df for df in (df1, df2)), ignore_index=True), expected)
        assert_frame_equal(pd.concat(deque((df1, df2)), ignore_index=True), expected)
        class CustomIterator1(object):
            def __len__(self):
                return 2
            def __getitem__(self, index):
                try:
                    return {0: df1, 1: df2}[index]
                except KeyError:
                    raise IndexError
        assert_frame_equal(pd.concat(CustomIterator1(), ignore_index=True), expected)
        class CustomIterator2(Iterable):
            def __iter__(self):
                yield df1
                yield df2
        assert_frame_equal(pd.concat(CustomIterator2(), ignore_index=True), expected)

    def test_concat_invalid(self):

        # trying to concat a ndframe with a non-ndframe
        df1 = mkdf(10, 2)
        for obj in [1, dict(), [1, 2], (1, 2) ]:
            self.assertRaises(TypeError, lambda x: concat([ df1, obj ]))

    def test_concat_invalid_first_argument(self):
        df1 = mkdf(10, 2)
        df2 = mkdf(10, 2)
        self.assertRaises(TypeError, concat, df1, df2)

        # generator ok though
        concat(DataFrame(np.random.rand(5,5)) for _ in range(3))

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
        assert_frame_equal(result,expected)

class TestOrderedMerge(tm.TestCase):

    def setUp(self):
        self.left = DataFrame({'key': ['a', 'c', 'e'],
                               'lvalue': [1, 2., 3]})

        self.right = DataFrame({'key': ['b', 'c', 'd', 'f'],
                                'rvalue': [1, 2, 3., 4]})

    # GH #813

    def test_basic(self):
        result = ordered_merge(self.left, self.right, on='key')
        expected = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'f'],
                              'lvalue': [1, nan, 2, nan, 3, nan],
                              'rvalue': [nan, 1, 2, 3, nan, 4]})

        assert_frame_equal(result, expected)

    def test_ffill(self):
        result = ordered_merge(
            self.left, self.right, on='key', fill_method='ffill')
        expected = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'f'],
                              'lvalue': [1., 1, 2, 2, 3, 3.],
                              'rvalue': [nan, 1, 2, 3, 3, 4]})
        assert_frame_equal(result, expected)

    def test_multigroup(self):
        left = concat([self.left, self.left], ignore_index=True)
        # right = concat([self.right, self.right], ignore_index=True)

        left['group'] = ['a'] * 3 + ['b'] * 3
        # right['group'] = ['a'] * 4 + ['b'] * 4

        result = ordered_merge(left, self.right, on='key', left_by='group',
                               fill_method='ffill')
        expected = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'f'] * 2,
                              'lvalue': [1., 1, 2, 2, 3, 3.] * 2,
                              'rvalue': [nan, 1, 2, 3, 3, 4] * 2})
        expected['group'] = ['a'] * 6 + ['b'] * 6

        assert_frame_equal(result, expected.ix[:, result.columns])

        result2 = ordered_merge(self.right, left, on='key', right_by='group',
                                fill_method='ffill')
        assert_frame_equal(result, result2.ix[:, result.columns])

        result = ordered_merge(left, self.right, on='key', left_by='group')
        self.assertTrue(result['group'].notnull().all())

    def test_merge_type(self):
        class NotADataFrame(DataFrame):
            @property
            def _constructor(self):
                return NotADataFrame

        nad = NotADataFrame(self.left)
        result = nad.merge(self.right, on='key')

        tm.assertIsInstance(result, NotADataFrame)

    def test_empty_sequence_concat(self):
        # GH 9157
        empty_pat = "[Nn]o objects"
        none_pat = "objects.*None"
        test_cases = [
            ((), empty_pat),
            ([], empty_pat),
            ({}, empty_pat),
            ([None], none_pat),
            ([None, None], none_pat)
        ]
        for df_seq, pattern in test_cases:
            assertRaisesRegexp(ValueError, pattern, pd.concat, df_seq)

        pd.concat([pd.DataFrame()])
        pd.concat([None, pd.DataFrame()])
        pd.concat([pd.DataFrame(), None])

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
