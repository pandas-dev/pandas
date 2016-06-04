# pylint: disable=E1103

import nose

from datetime import datetime
from numpy.random import randn
from numpy import nan
import numpy as np
import random

import pandas as pd
from pandas.compat import range, lrange, lzip
from pandas.tools.merge import merge, concat, MergeError
from pandas.util.testing import (assert_frame_equal,
                                 assert_series_equal,
                                 slow)
from pandas import (DataFrame, Index, MultiIndex,
                    Series, date_range, Categorical,
                    compat)
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
        self.assert_series_equal(merged['MergedA'], target['A'],
                                 check_names=False)
        self.assert_series_equal(merged['MergedD'], target['D'],
                                 check_names=False)

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
            df = DataFrame({'a': np.random.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)})
            df2 = DataFrame({'a': np.random.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)},
                            index=tm.makeCustomIndex(10, 2))
            merge(df, df2, left_on='a', right_index=True)

    def test_join_on_fails_with_different_left_index(self):
        with tm.assertRaises(ValueError):
            df = DataFrame({'a': np.random.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)},
                           index=tm.makeCustomIndex(10, 2))
            df2 = DataFrame({'a': np.random.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)})
            merge(df, df2, right_on='b', left_index=True)

    def test_join_on_fails_with_different_column_counts(self):
        with tm.assertRaises(ValueError):
            df = DataFrame({'a': np.random.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)})
            df2 = DataFrame({'a': np.random.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)},
                            index=tm.makeCustomIndex(10, 2))
            merge(df, df2, right_on='a', left_on=['a', 'b'])

    def test_join_on_fails_with_wrong_object_type(self):
        # GH12081
        wrongly_typed = [Series([0, 1]), 2, 'str', None, np.array([0, 1])]
        df = DataFrame({'a': [1, 1]})

        for obj in wrongly_typed:
            with tm.assertRaisesRegexp(ValueError, str(type(obj))):
                merge(obj, df, left_on='a', right_on='a')
            with tm.assertRaisesRegexp(ValueError, str(type(obj))):
                merge(df, obj, left_on='a', right_on='a')

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
        self.assert_index_equal(merged2.columns, merged.columns)
        self.assertEqual(len(merged2), 0)

    def test_join_on_inner(self):
        df = DataFrame({'key': ['a', 'a', 'd', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1]}, index=['a', 'b'])

        joined = df.join(df2, on='key', how='inner')

        expected = df.join(df2, on='key')
        expected = expected[expected['value'].notnull()]
        self.assert_series_equal(joined['key'], expected['key'],
                                 check_dtype=False)
        self.assert_series_equal(joined['value'], expected['value'],
                                 check_dtype=False)
        self.assert_index_equal(joined.index, expected.index)

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
        # GH 2024
        df = DataFrame([(1, 2, 3), (4, 5, 6)], columns=['a', 'b', 'c'])
        new_df = df.groupby(['a']).agg({'b': [np.mean, np.sum]})
        other_df = DataFrame(
            [(1, 2, 3), (7, 10, 6)], columns=['a', 'b', 'd'])
        other_df.set_index('a', inplace=True)
        # GH 9455, 12219
        with tm.assert_produces_warning(UserWarning):
            result = merge(new_df, other_df, left_index=True, right_index=True)
        self.assertTrue(('b', 'mean') in result)
        self.assertTrue('b' in result)

    def test_join_float64_float32(self):

        a = DataFrame(randn(10, 2), columns=['a', 'b'], dtype=np.float64)
        b = DataFrame(randn(10, 1), columns=['c'], dtype=np.float32)
        joined = a.join(b)
        self.assertEqual(joined.dtypes['a'], 'float64')
        self.assertEqual(joined.dtypes['b'], 'float64')
        self.assertEqual(joined.dtypes['c'], 'float32')

        a = np.random.randint(0, 5, 100).astype('int64')
        b = np.random.random(100).astype('float64')
        c = np.random.random(100).astype('float32')
        df = DataFrame({'a': a, 'b': b, 'c': c})
        xpdf = DataFrame({'a': a, 'b': b, 'c': c})
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
        expected = expected[result.columns]
        expected['a'] = expected.a.astype('int64')
        expected['b'] = expected.b.astype('int64')
        assert_frame_equal(result, expected)

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

        # GH 11519
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'three',
                              'two', 'two', 'one', 'three'],
                        'C': np.random.randn(8),
                        'D': np.random.randn(8)})
        s = Series(np.repeat(np.arange(8), 2),
                   index=np.repeat(np.arange(8), 2), name='TEST')
        inner = df.join(s, how='inner')
        outer = df.join(s, how='outer')
        left = df.join(s, how='left')
        right = df.join(s, how='right')
        assert_frame_equal(inner, outer)
        assert_frame_equal(inner, left)
        assert_frame_equal(inner, right)

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

        exp = pd.Series(['bar', 'baz', 'foo', 'foo', 'foo', 'foo', np.nan],
                        name='lkey')
        tm.assert_series_equal(merged['lkey'], exp)

        exp = pd.Series(['bar', np.nan, 'foo', 'foo', 'foo', 'foo', 'qux'],
                        name='rkey')
        tm.assert_series_equal(merged['rkey'], exp)

        exp = pd.Series([2, 3, 1, 1, 4, 4, np.nan], name='value_x')
        tm.assert_series_equal(merged['value_x'], exp)

        exp = pd.Series([6, np.nan, 5, 8, 5, 8, 7], name='value_y')
        tm.assert_series_equal(merged['value_y'], exp)

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
        self.assert_index_equal(joined.index, pd.Index(lrange(4)))

    def test_intelligently_handle_join_key(self):
        # #733, be a bit more 1337 about not returning unconsolidated DataFrame

        left = DataFrame({'key': [1, 1, 2, 2, 3],
                          'value': lrange(5)}, columns=['value', 'key'])
        right = DataFrame({'key': [1, 1, 2, 3, 4, 5],
                           'rvalue': lrange(6)})

        joined = merge(left, right, on='key', how='outer')
        expected = DataFrame({'key': [1, 1, 1, 1, 2, 2, 3, 4, 5],
                              'value': np.array([0, 0, 1, 1, 2, 3, 4,
                                                 np.nan, np.nan]),
                              'rvalue': [0, 1, 0, 1, 2, 2, 3, 4, 5]},
                             columns=['value', 'key', 'rvalue'])
        assert_frame_equal(joined, expected)

    def test_merge_join_key_dtype_cast(self):
        # #8596

        df1 = DataFrame({'key': [1], 'v1': [10]})
        df2 = DataFrame({'key': [2], 'v1': [20]})
        df = merge(df1, df2, how='outer')
        self.assertEqual(df['key'].dtype, 'int64')

        df1 = DataFrame({'key': [True], 'v1': [1]})
        df2 = DataFrame({'key': [False], 'v1': [0]})
        df = merge(df1, df2, how='outer')

        # GH13169
        # this really should be bool
        self.assertEqual(df['key'].dtype, 'object')

        df1 = DataFrame({'val': [1]})
        df2 = DataFrame({'val': [2]})
        lkey = np.array([1])
        rkey = np.array([2])
        df = merge(df1, df2, left_on=lkey, right_on=rkey, how='outer')
        self.assertEqual(df['key_0'].dtype, 'int64')

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
        self.assert_series_equal(merged['key_0'],
                                 Series([1, 1, 1, 1, 2, 2, 3, 4, 5],
                                        name='key_0'))

        left = DataFrame({'value': lrange(3)})
        right = DataFrame({'rvalue': lrange(6)})

        key = np.array([0, 1, 1, 2, 2, 3], dtype=np.int64)
        merged = merge(left, right, left_index=True, right_on=key, how='outer')
        self.assert_series_equal(merged['key_0'], Series(key, name='key_0'))

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
                              index=pd.Index([], dtype=object),
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

        exp_out = pd.DataFrame({'a': np.array([np.nan] * 3, dtype=object),
                                'b': np.array([np.nan] * 3, dtype=object),
                                'c': np.array([np.nan] * 3, dtype=object),
                                'x': [1, 4, 7],
                                'y': [2, 5, 8],
                                'z': [3, 6, 9]},
                               columns=['a', 'b', 'c', 'x', 'y', 'z'])
        exp_in = exp_out[0:0]  # make empty DataFrame keeping dtype
        # result will have object dtype
        exp_in.index = exp_in.index.astype(object)

        def check1(exp, kwarg):
            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp)

        def check2(exp, kwarg):
            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp)

        for kwarg in [dict(left_index=True, right_index=True),
                      dict(left_index=True, right_on='x')]:
            check1(exp_in, kwarg)
            check2(exp_out, kwarg)

        kwarg = dict(left_on='a', right_index=True)
        check1(exp_in, kwarg)
        exp_out['a'] = [0, 1, 2]
        check2(exp_out, kwarg)

        kwarg = dict(left_on='a', right_on='x')
        check1(exp_in, kwarg)
        exp_out['a'] = np.array([np.nan] * 3, dtype=object)
        check2(exp_out, kwarg)

    def test_merge_left_notempty_right_empty(self):
        # GH 10824
        left = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            columns=['a', 'b', 'c'])
        right = pd.DataFrame([], columns=['x', 'y', 'z'])

        exp_out = pd.DataFrame({'a': [1, 4, 7],
                                'b': [2, 5, 8],
                                'c': [3, 6, 9],
                                'x': np.array([np.nan] * 3, dtype=object),
                                'y': np.array([np.nan] * 3, dtype=object),
                                'z': np.array([np.nan] * 3, dtype=object)},
                               columns=['a', 'b', 'c', 'x', 'y', 'z'])
        exp_in = exp_out[0:0]  # make empty DataFrame keeping dtype
        # result will have object dtype
        exp_in.index = exp_in.index.astype(object)

        def check1(exp, kwarg):
            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp)

        def check2(exp, kwarg):
            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp)

            for kwarg in [dict(left_index=True, right_index=True),
                          dict(left_index=True, right_on='x'),
                          dict(left_on='a', right_index=True),
                          dict(left_on='a', right_on='x')]:
                check1(exp_in, kwarg)
                check2(exp_out, kwarg)

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
        df1 = DataFrame({"i1": [0, 1], "i2": [0, 1]})
        df2 = DataFrame({"i1": [0], "i3": [0]})
        result = df1.join(df2, on="i1", rsuffix="_")
        expected = (DataFrame({'i1': {0: 0.0, 1: 1}, 'i2': {0: 0, 1: 1},
                               'i1_': {0: 0, 1: np.nan},
                               'i3': {0: 0.0, 1: np.nan},
                               None: {0: 0, 1: 0}})
                    .set_index(None)
                    .reset_index()[['i1', 'i2', 'i1_', 'i3']])
        assert_frame_equal(result, expected, check_dtype=False)

        df1 = DataFrame({"i1": [0, 1], "i2": [0.5, 1.5]})
        df2 = DataFrame({"i1": [0], "i3": [0.7]})
        result = df1.join(df2, rsuffix="_", on='i1')
        expected = (DataFrame({'i1': {0: 0, 1: 1}, 'i1_': {0: 0.0, 1: nan},
                               'i2': {0: 0.5, 1: 1.5},
                               'i3': {0: 0.69999999999999996,
                                      1: nan}})
                    [['i1', 'i2', 'i1_', 'i3']])
        assert_frame_equal(result, expected)

    def test_merge_type(self):
        class NotADataFrame(DataFrame):

            @property
            def _constructor(self):
                return NotADataFrame

        nad = NotADataFrame(self.df)
        result = nad.merge(self.df2, on='key1')

        tm.assertIsInstance(result, NotADataFrame)

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
        assert_series_equal(result.dtypes, expected)

    def test_dtype_coerceion(self):

        # 12411
        df = DataFrame({'date': [pd.Timestamp('20130101').tz_localize('UTC'),
                                 pd.NaT]})

        result = concat([df.iloc[[0]], df.iloc[[1]]])
        assert_series_equal(result.dtypes, df.dtypes)

        # 12045
        import datetime
        df = DataFrame({'date': [datetime.datetime(2012, 1, 1),
                                 datetime.datetime(1012, 1, 2)]})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        assert_series_equal(result.dtypes, df.dtypes)

        # 11594
        df = DataFrame({'text': ['some words'] + [None] * 9})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        assert_series_equal(result.dtypes, df.dtypes)

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
                                    dt.datetime(2013, 11, 5, 5, 56)],
                              't': [dt.timedelta(0, 22500),
                                    dt.timedelta(0, 22500)]})
        assert_frame_equal(result, expected)

        td = np.timedelta64(300000000)
        lhs = DataFrame(Series([td, td], index=["A", "B"]))
        rhs = DataFrame(Series([td], index=["A"]))

        result = lhs.join(rhs, rsuffix='r', how="left")
        expected = DataFrame({'0': Series([td, td], index=list('AB')),
                              '0r': Series([td, NaT], index=list('AB'))})
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

    def test_merge_on_datetime64tz(self):

        # GH11405
        left = pd.DataFrame({'key': pd.date_range('20151010', periods=2,
                                                  tz='US/Eastern'),
                             'value': [1, 2]})
        right = pd.DataFrame({'key': pd.date_range('20151011', periods=3,
                                                   tz='US/Eastern'),
                              'value': [1, 2, 3]})

        expected = DataFrame({'key': pd.date_range('20151010', periods=4,
                                                   tz='US/Eastern'),
                              'value_x': [1, 2, np.nan, np.nan],
                              'value_y': [np.nan, 1, 2, 3]})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)

        left = pd.DataFrame({'value': pd.date_range('20151010', periods=2,
                                                    tz='US/Eastern'),
                             'key': [1, 2]})
        right = pd.DataFrame({'value': pd.date_range('20151011', periods=2,
                                                     tz='US/Eastern'),
                              'key': [2, 3]})
        expected = DataFrame({
            'value_x': list(pd.date_range('20151010', periods=2,
                                          tz='US/Eastern')) + [pd.NaT],
            'value_y': [pd.NaT] + list(pd.date_range('20151011', periods=2,
                                                     tz='US/Eastern')),
            'key': [1, 2, 3]})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)
        self.assertEqual(result['value_x'].dtype, 'datetime64[ns, US/Eastern]')
        self.assertEqual(result['value_y'].dtype, 'datetime64[ns, US/Eastern]')

    def test_merge_on_periods(self):
        left = pd.DataFrame({'key': pd.period_range('20151010', periods=2,
                                                    freq='D'),
                             'value': [1, 2]})
        right = pd.DataFrame({'key': pd.period_range('20151011', periods=3,
                                                     freq='D'),
                              'value': [1, 2, 3]})

        expected = DataFrame({'key': pd.period_range('20151010', periods=4,
                                                     freq='D'),
                              'value_x': [1, 2, np.nan, np.nan],
                              'value_y': [np.nan, 1, 2, 3]})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)

        left = pd.DataFrame({'value': pd.period_range('20151010', periods=2,
                                                      freq='D'),
                             'key': [1, 2]})
        right = pd.DataFrame({'value': pd.period_range('20151011', periods=2,
                                                       freq='D'),
                              'key': [2, 3]})

        exp_x = pd.period_range('20151010', periods=2, freq='D')
        exp_y = pd.period_range('20151011', periods=2, freq='D')
        expected = DataFrame({'value_x': list(exp_x) + [pd.NaT],
                              'value_y': [pd.NaT] + list(exp_y),
                              'key': [1, 2, 3]})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)
        self.assertEqual(result['value_x'].dtype, 'object')
        self.assertEqual(result['value_y'].dtype, 'object')

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

    def test_indicator(self):
        # PR #10054. xref #7412 and closes #8790.
        df1 = DataFrame({'col1': [0, 1], 'col_left': [
                        'a', 'b'], 'col_conflict': [1, 2]})
        df1_copy = df1.copy()

        df2 = DataFrame({'col1': [1, 2, 3, 4, 5], 'col_right': [2, 2, 2, 2, 2],
                         'col_conflict': [1, 2, 3, 4, 5]})
        df2_copy = df2.copy()

        df_result = DataFrame({
            'col1': [0, 1, 2, 3, 4, 5],
            'col_conflict_x': [1, 2, np.nan, np.nan, np.nan, np.nan],
            'col_left': ['a', 'b', np.nan, np.nan, np.nan, np.nan],
            'col_conflict_y': [np.nan, 1, 2, 3, 4, 5],
            'col_right': [np.nan, 2, 2, 2, 2, 2]})
        df_result['_merge'] = Categorical(
            ['left_only', 'both', 'right_only',
             'right_only', 'right_only', 'right_only'],
            categories=['left_only', 'right_only', 'both'])

        df_result = df_result[['col1', 'col_conflict_x', 'col_left',
                               'col_conflict_y', 'col_right', '_merge']]

        test = merge(df1, df2, on='col1', how='outer', indicator=True)
        assert_frame_equal(test, df_result)
        test = df1.merge(df2, on='col1', how='outer', indicator=True)
        assert_frame_equal(test, df_result)

        # No side effects
        assert_frame_equal(df1, df1_copy)
        assert_frame_equal(df2, df2_copy)

        # Check with custom name
        df_result_custom_name = df_result
        df_result_custom_name = df_result_custom_name.rename(
            columns={'_merge': 'custom_name'})

        test_custom_name = merge(
            df1, df2, on='col1', how='outer', indicator='custom_name')
        assert_frame_equal(test_custom_name, df_result_custom_name)
        test_custom_name = df1.merge(
            df2, on='col1', how='outer', indicator='custom_name')
        assert_frame_equal(test_custom_name, df_result_custom_name)

        # Check only accepts strings and booleans
        with tm.assertRaises(ValueError):
            merge(df1, df2, on='col1', how='outer', indicator=5)
        with tm.assertRaises(ValueError):
            df1.merge(df2, on='col1', how='outer', indicator=5)

        # Check result integrity

        test2 = merge(df1, df2, on='col1', how='left', indicator=True)
        self.assertTrue((test2._merge != 'right_only').all())
        test2 = df1.merge(df2, on='col1', how='left', indicator=True)
        self.assertTrue((test2._merge != 'right_only').all())

        test3 = merge(df1, df2, on='col1', how='right', indicator=True)
        self.assertTrue((test3._merge != 'left_only').all())
        test3 = df1.merge(df2, on='col1', how='right', indicator=True)
        self.assertTrue((test3._merge != 'left_only').all())

        test4 = merge(df1, df2, on='col1', how='inner', indicator=True)
        self.assertTrue((test4._merge == 'both').all())
        test4 = df1.merge(df2, on='col1', how='inner', indicator=True)
        self.assertTrue((test4._merge == 'both').all())

        # Check if working name in df
        for i in ['_right_indicator', '_left_indicator', '_merge']:
            df_badcolumn = DataFrame({'col1': [1, 2], i: [2, 2]})

            with tm.assertRaises(ValueError):
                merge(df1, df_badcolumn, on='col1',
                      how='outer', indicator=True)
            with tm.assertRaises(ValueError):
                df1.merge(df_badcolumn, on='col1', how='outer', indicator=True)

        # Check for name conflict with custom name
        df_badcolumn = DataFrame(
            {'col1': [1, 2], 'custom_column_name': [2, 2]})

        with tm.assertRaises(ValueError):
            merge(df1, df_badcolumn, on='col1', how='outer',
                  indicator='custom_column_name')
        with tm.assertRaises(ValueError):
            df1.merge(df_badcolumn, on='col1', how='outer',
                      indicator='custom_column_name')

        # Merge on multiple columns
        df3 = DataFrame({'col1': [0, 1], 'col2': ['a', 'b']})

        df4 = DataFrame({'col1': [1, 1, 3], 'col2': ['b', 'x', 'y']})

        hand_coded_result = DataFrame({'col1': [0, 1, 1, 3],
                                       'col2': ['a', 'b', 'x', 'y']})
        hand_coded_result['_merge'] = Categorical(
            ['left_only', 'both', 'right_only', 'right_only'],
            categories=['left_only', 'right_only', 'both'])

        test5 = merge(df3, df4, on=['col1', 'col2'],
                      how='outer', indicator=True)
        assert_frame_equal(test5, hand_coded_result)
        test5 = df3.merge(df4, on=['col1', 'col2'],
                          how='outer', indicator=True)
        assert_frame_equal(test5, hand_coded_result)


def _check_merge(x, y):
    for how in ['inner', 'left', 'outer']:
        result = x.join(y, how=how)

        expected = merge(x.reset_index(), y.reset_index(), how=how,
                         sort=True)
        expected = expected.set_index('index')

        # TODO check_names on merge?
        assert_frame_equal(result, expected, check_names=False)


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
            return (f(df['1st']) + f(df['3rd']) * 1e2 +
                    df['2nd'].fillna(0) * 1e4)

        def run_asserts(left, right):
            for sort in [False, True]:
                res = left.join(right, on=icols, how='left', sort=sort)

                self.assertTrue(len(left) < len(res) + 1)
                self.assertFalse(res['4th'].isnull().any())
                self.assertFalse(res['5th'].isnull().any())

                tm.assert_series_equal(
                    res['4th'], - res['5th'], check_names=False)
                result = bind_cols(res.iloc[:, :-2])
                tm.assert_series_equal(res['4th'], result, check_names=False)
                self.assertTrue(result.name is None)

                if sort:
                    tm.assert_frame_equal(
                        res, res.sort_values(icols, kind='mergesort'))

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
                                         left_index=True, how='right',
                                         sort=sort)

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
        merge(df, df2, how='outer')

    def test_left_join_index_preserve_order(self):

        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'v': np.array(np.arange(24), dtype=np.int64)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(
            result.sort_values(['k1', 'k2'], kind='mergesort'),
            left.join(right, on=['k1', 'k2'], sort=True))

        # test join with multi dtypes blocks
        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'k3': np.array([0, 1, 2] * 8, dtype=np.float32),
                          'v': np.array(np.arange(24), dtype=np.int32)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(
            result.sort_values(['k1', 'k2'], kind='mergesort'),
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
            ['W', 'R', 'C', 0],
            ['W', 'Q', 'B', 3],
            ['W', 'Q', 'B', 8],
            ['X', 'Y', 'A', 1],
            ['X', 'Y', 'A', 4],
            ['X', 'Y', 'B', 5],
            ['X', 'Y', 'C', 6],
            ['X', 'Y', 'C', 9],
            ['X', 'Q', 'C', -6],
            ['X', 'R', 'C', -9],
            ['V', 'Y', 'C', 7],
            ['V', 'R', 'D', 2],
            ['V', 'R', 'D', -1],
            ['V', 'Q', 'A', -3]],
            columns=['col1', 'col2', 'col3', 'val'])

        right.set_index(['col1', 'col2', 'col3'], inplace=True)
        result = left.join(right, on=['cola', 'colb', 'colc'], how='left')

        expected = DataFrame([
            ['X', 'Y', 'C', 'a', 6],
            ['X', 'Y', 'C', 'a', 9],
            ['W', 'Y', 'C', 'e', nan],
            ['V', 'Q', 'A', 'h', -3],
            ['V', 'R', 'D', 'i', 2],
            ['V', 'R', 'D', 'i', -1],
            ['X', 'Y', 'D', 'b', nan],
            ['X', 'Y', 'A', 'c', 1],
            ['X', 'Y', 'A', 'c', 4],
            ['W', 'Q', 'B', 'f', 3],
            ['W', 'Q', 'B', 'f', 8],
            ['W', 'R', 'C', 'g', 0],
            ['V', 'Y', 'C', 'j', 7],
            ['X', 'Y', 'B', 'd', 5]],
            columns=['cola', 'colb', 'colc', 'tag', 'val'],
            index=[3, 3, 2, 0, 1, 1, 7, 6, 6, 4, 4, 5, 9, 8])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on=['cola', 'colb', 'colc'],
                           how='left', sort=True)

        tm.assert_frame_equal(
            result,
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
        tm.assert_frame_equal(
            result, expected.sort_values('tag', kind='mergesort'))

        # GH7331 - maintain left frame order in left merge
        result = merge(left, right.reset_index(), how='left', on='tag')
        expected.index = np.arange(len(expected))
        tm.assert_frame_equal(result, expected)

    def test_join_multi_dtypes(self):

        # test with multi dtypes in the join index
        def _test(dtype1, dtype2):
            left = DataFrame({'k1': np.array([0, 1, 2] * 8, dtype=dtype1),
                              'k2': ['foo', 'bar'] * 12,
                              'v': np.array(np.arange(24), dtype=np.int64)})

            index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
            right = DataFrame(
                {'v2': np.array([5, 7], dtype=dtype2)}, index=index)

            result = left.join(right, on=['k1', 'k2'])

            expected = left.copy()

            if dtype2.kind == 'i':
                dtype2 = np.dtype('float64')
            expected['v2'] = np.array(np.nan, dtype=dtype2)
            expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
            expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

            tm.assert_frame_equal(result, expected)

            result = left.join(right, on=['k1', 'k2'], sort=True)
            expected.sort_values(['k1', 'k2'], kind='mergesort', inplace=True)
            tm.assert_frame_equal(result, expected)

        for d1 in [np.int64, np.int32, np.int16, np.int8, np.uint8]:
            for d2 in [np.int64, np.float64, np.float32, np.float16]:
                _test(np.dtype(d1), np.dtype(d2))

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

    @slow
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
        shape = left.apply(Series.nunique).values
        self.assertTrue(_int64_overflow_possible(shape))

        # add duplicates to left frame
        left = concat([left, left], ignore_index=True)

        right = DataFrame(np.random.randint(low, high, (n // 2, 7))
                          .astype('int64'),
                          columns=list('ABCDEFG'))

        # add duplicates & overlap with left to the right frame
        i = np.random.choice(len(left), n)
        right = concat([right, right, left.iloc[i]], ignore_index=True)

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
        household = (
            DataFrame(
                dict(household_id=[1, 2, 3],
                     male=[0, 1, 0],
                     wealth=[196087.3, 316478.7, 294750]),
                columns=['household_id', 'male', 'wealth'])
            .set_index('household_id'))
        portfolio = (
            DataFrame(
                dict(household_id=[1, 2, 2, 3, 3, 3, 4],
                     asset_id=["nl0000301109", "nl0000289783", "gb00b03mlx29",
                               "gb00b03mlx29", "lu0197800237", "nl0000289965",
                               np.nan],
                     name=["ABN Amro", "Robeco", "Royal Dutch Shell",
                           "Royal Dutch Shell",
                           "AAB Eastern Europe Equity Fund",
                           "Postbank BioTech Fonds", np.nan],
                     share=[1.0, 0.4, 0.6, 0.15, 0.6, 0.25, 1.0]),
                columns=['household_id', 'asset_id', 'name', 'share'])
            .set_index(['household_id', 'asset_id']))
        result = household.join(portfolio, how='inner')
        expected = (
            DataFrame(
                dict(male=[0, 1, 1, 0, 0, 0],
                     wealth=[196087.3, 316478.7, 316478.7,
                             294750.0, 294750.0, 294750.0],
                     name=['ABN Amro', 'Robeco', 'Royal Dutch Shell',
                           'Royal Dutch Shell',
                           'AAB Eastern Europe Equity Fund',
                           'Postbank BioTech Fonds'],
                     share=[1.00, 0.40, 0.60, 0.15, 0.60, 0.25],
                     household_id=[1, 2, 2, 3, 3, 3],
                     asset_id=['nl0000301109', 'nl0000289783', 'gb00b03mlx29',
                               'gb00b03mlx29', 'lu0197800237',
                               'nl0000289965']))
            .set_index(['household_id', 'asset_id'])
            .reindex(columns=['male', 'wealth', 'name', 'share']))
        assert_frame_equal(result, expected)

        assert_frame_equal(result, expected)

        # equivalency
        result2 = (merge(household.reset_index(), portfolio.reset_index(),
                         on=['household_id'], how='inner')
                   .set_index(['household_id', 'asset_id']))
        assert_frame_equal(result2, expected)

        result = household.join(portfolio, how='outer')
        expected = (concat([
            expected,
            (DataFrame(
                dict(share=[1.00]),
                index=MultiIndex.from_tuples(
                    [(4, np.nan)],
                    names=['household_id', 'asset_id'])))
        ], axis=0).reindex(columns=expected.columns))
        assert_frame_equal(result, expected)

        # invalid cases
        household.index.name = 'foo'

        def f():
            household.join(portfolio, how='inner')
        self.assertRaises(ValueError, f)

        portfolio2 = portfolio.copy()
        portfolio2.index.set_names(['household_id', 'foo'])

        def f():
            portfolio2.join(portfolio, how='inner')
        self.assertRaises(ValueError, f)

    def test_join_multi_levels2(self):

        # some more advanced merges
        # GH6360
        household = (
            DataFrame(
                dict(household_id=[1, 2, 2, 3, 3, 3, 4],
                     asset_id=["nl0000301109", "nl0000301109", "gb00b03mlx29",
                               "gb00b03mlx29", "lu0197800237", "nl0000289965",
                               np.nan],
                     share=[1.0, 0.4, 0.6, 0.15, 0.6, 0.25, 1.0]),
                columns=['household_id', 'asset_id', 'share'])
            .set_index(['household_id', 'asset_id']))

        log_return = DataFrame(dict(
            asset_id=["gb00b03mlx29", "gb00b03mlx29",
                      "gb00b03mlx29", "lu0197800237", "lu0197800237"],
            t=[233, 234, 235, 180, 181],
            log_return=[.09604978, -.06524096, .03532373, .03025441, .036997]
        )).set_index(["asset_id", "t"])

        expected = (
            DataFrame(dict(
                household_id=[2, 2, 2, 3, 3, 3, 3, 3],
                asset_id=["gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "lu0197800237", "lu0197800237"],
                t=[233, 234, 235, 233, 234, 235, 180, 181],
                share=[0.6, 0.6, 0.6, 0.15, 0.15, 0.15, 0.6, 0.6],
                log_return=[.09604978, -.06524096, .03532373,
                            .09604978, -.06524096, .03532373,
                            .03025441, .036997]
            ))
            .set_index(["household_id", "asset_id", "t"])
            .reindex(columns=['share', 'log_return']))

        def f():
            household.join(log_return, how='inner')
        self.assertRaises(NotImplementedError, f)

        # this is the equivalency
        result = (merge(household.reset_index(), log_return.reset_index(),
                        on=['asset_id'], how='inner')
                  .set_index(['household_id', 'asset_id', 't']))
        assert_frame_equal(result, expected)

        expected = (
            DataFrame(dict(
                household_id=[1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4],
                asset_id=["nl0000301109", "nl0000289783", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29",
                          "lu0197800237", "lu0197800237",
                          "nl0000289965", None],
                t=[None, None, 233, 234, 235, 233, 234,
                   235, 180, 181, None, None],
                share=[1.0, 0.4, 0.6, 0.6, 0.6, 0.15,
                       0.15, 0.15, 0.6, 0.6, 0.25, 1.0],
                log_return=[None, None, .09604978, -.06524096, .03532373,
                            .09604978, -.06524096, .03532373,
                            .03025441, .036997, None, None]
            ))
            .set_index(["household_id", "asset_id", "t"]))

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


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
