# pylint: disable=E1103

import nose
import unittest

from numpy.random import randn
import numpy as np
import random

from pandas import *
from pandas.tools.merge import merge, concat
from pandas.util.testing import (assert_frame_equal, assert_series_equal,
                                 assert_almost_equal, rands)
import pandas._tseries as lib
import pandas.util.testing as tm

a_ = np.array

N = 50
NGROUPS = 8
JOIN_TYPES = ['inner', 'outer', 'left', 'right']


def get_test_data(ngroups=NGROUPS, n=N):
    unique_groups = range(ngroups)
    arr = np.asarray(np.tile(unique_groups, n // ngroups))

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)])

    random.shuffle(arr)
    return arr

class TestMerge(unittest.TestCase):

    def setUp(self):
        # aggregate multiple columns
        self.df = DataFrame({'key1': get_test_data(),
                             'key2': get_test_data(),
                             'data1': np.random.randn(N),
                             'data2': np.random.randn(N)})

        # exclude a couple keys for fun
        self.df = self.df[self.df['key2'] > 1]

        self.df2 = DataFrame({'key1'  : get_test_data(n=N//5),
                              'key2'  : get_test_data(ngroups=NGROUPS//2,
                                                      n=N//5),
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
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype='i4')
        max_group = 5

        ls, rs = lib.left_outer_join(left, right, max_group)

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

        self.assert_(np.array_equal(ls, exp_ls))
        self.assert_(np.array_equal(rs, exp_rs))

    def test_cython_right_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype='i4')
        max_group = 5

        rs, ls  = lib.left_outer_join(right, left, max_group)

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

        self.assert_(np.array_equal(ls, exp_ls))
        self.assert_(np.array_equal(rs, exp_rs))

    def test_cython_inner_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1, 4], dtype='i4')
        max_group = 5

        ls, rs = lib.inner_join(left, right, max_group)

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

        self.assert_(np.array_equal(ls, exp_ls))
        self.assert_(np.array_equal(rs, exp_rs))

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

        self.assert_('key1.foo' in joined)
        self.assert_('key1.bar' in joined)

    def test_handle_overlap_arbitrary_key(self):
        joined = merge(self.df, self.df2,
                       left_on='key2', right_on='key1',
                       suffixes=['.foo', '.bar'])
        self.assert_('key1.foo' in joined)
        self.assert_('key2.bar' in joined)

    def test_merge_common(self):
        joined = merge(self.df, self.df2)
        exp = merge(self.df, self.df2, on=['key1', 'key2'])
        tm.assert_frame_equal(joined, exp)

    def test_join_on(self):
        target = self.target
        source = self.source

        merged = target.join(source, on='C')
        self.assert_(np.array_equal(merged['MergedA'], target['A']))
        self.assert_(np.array_equal(merged['MergedD'], target['D']))

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
        self.assert_(np.isnan(joined['two']['c']))
        self.assert_(np.isnan(joined['three']['c']))

        # merge column not p resent
        self.assertRaises(Exception, target.join, source, on='E')

        # overlap
        source_copy = source.copy()
        source_copy['A'] = 0
        self.assertRaises(Exception, target.join, source_copy, on='A')

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
            self.assert_(col in merged)
            self.assert_(merged[col].isnull().all())

        merged2 = self.target.join(self.source.reindex([]), on='C',
                                   how='inner')
        self.assert_(merged2.columns.equals(merged.columns))
        self.assertEqual(len(merged2), 0)

    def test_join_on_inner(self):
        df = DataFrame({'key': ['a', 'a', 'd', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1]}, index=['a', 'b'])

        joined = df.join(df2, on='key', how='inner')

        expected = df.join(df2, on='key')
        expected = expected[expected['value'].notnull()]
        self.assert_(np.array_equal(joined['key'], expected['key']))
        self.assert_(np.array_equal(joined['value'], expected['value']))
        self.assert_(joined.index.equals(expected.index))

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
        expected = DataFrame({'a' : [1, 1],
                              'b': [2, 2]}, index=df.index)
        tm.assert_frame_equal(result, expected)

    def test_join_index_mixed(self):

        df1 = DataFrame({'A': 1., 'B': 2, 'C': 'foo', 'D': True},
                        index=np.arange(10),
                        columns=['A', 'B', 'C', 'D'])
        self.assert_(df1['B'].dtype == np.int64)
        self.assert_(df1['D'].dtype == np.bool_)

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
        a = DataFrame(randn(30,2), columns=['a','b'])
        c = Series(randn(30))
        a['c'] = c
        d = DataFrame(randn(30,1), columns=['q'])

        # it works!
        a.join(d)
        d.join(a)

    def test_join_multiindex(self):
        index1 = MultiIndex.from_arrays([['a','a','a','b','b','b'],
                                         [1,2,3,1,2,3]],
                                        names=['first', 'second'])

        index2 = MultiIndex.from_arrays([['b','b','b','c','c','c'],
                                         [1,2,3,1,2,3]],
                                        names=['first', 'second'])

        df1 = DataFrame(data=np.random.randn(6), index=index1,
                        columns=['var X'])
        df2 = DataFrame(data=np.random.randn(6), index=index2,
                        columns=['var Y'])

        df1 = df1.sortlevel(0)
        df2 = df2.sortlevel(0)

        joined = df1.join(df2, how='outer')
        ex_index = index1.get_tuple_index() + index2.get_tuple_index()
        expected = df1.reindex(ex_index).join(df2.reindex(ex_index))
        assert_frame_equal(joined, expected)
        self.assertEqual(joined.index.names, index1.names)

        df1 = df1.sortlevel(1)
        df2 = df2.sortlevel(1)

        joined = df1.join(df2, how='outer').sortlevel(0)
        ex_index = index1.get_tuple_index() + index2.get_tuple_index()
        expected = df1.reindex(ex_index).join(df2.reindex(ex_index))

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

        self.assert_(joined.index.is_monotonic)
        assert_frame_equal(joined, expected)

        # _assert_same_contents(expected, expected2.ix[:, expected.columns])

    def test_join_float64_float32(self):
        a = DataFrame(randn(10,2), columns=['a','b'])
        b = DataFrame(randn(10,1), columns=['c']).astype(np.float32)
        joined = a.join(b)
        expected = a.join(b.astype('f8'))
        assert_frame_equal(joined, expected)

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
        self.assertRaises(Exception, merge, self.left, self.right,
                          left_index=True)
        self.assertRaises(Exception, merge, self.left, self.right,
                          right_index=True)

        self.assertRaises(Exception, merge, self.left, self.left,
                          left_on='key', on='key')

        self.assertRaises(Exception, merge, self.df, self.df2,
                          left_on=['key1'], right_on=['key1', 'key2'])

    def test_merge_overlap(self):
        merged = merge(self.left, self.left, on='key')
        exp_len = (self.left['key'].value_counts() ** 2).sum()
        self.assertEqual(len(merged), exp_len)
        self.assert_('v1.x' in merged)
        self.assert_('v1.y' in merged)

    def test_merge_different_column_key_names(self):
        left = DataFrame({'lkey': ['foo', 'bar', 'baz', 'foo'],
                          'value': [1, 2, 3, 4]})
        right = DataFrame({'rkey': ['foo', 'bar', 'qux', 'foo'],
                           'value' : [5, 6, 7, 8]})

        merged = left.merge(right, left_on='lkey', right_on='rkey',
                            how='outer')

        assert_almost_equal(merged['lkey'],
                            ['bar', 'baz', 'foo', 'foo', 'foo', 'foo', np.nan])
        assert_almost_equal(merged['rkey'],
                            ['bar', np.nan, 'foo', 'foo', 'foo', 'foo', 'qux'])
        assert_almost_equal(merged['value.x'], [2, 3, 1, 1, 4, 4, np.nan])
        assert_almost_equal(merged['value.y'], [6, np.nan, 5, 8, 5, 8, 7])

    def test_merge_nocopy(self):
        left = DataFrame({'a' : 0, 'b' : 1}, index=range(10))
        right = DataFrame({'c' : 'foo', 'd' : 'bar'}, index=range(10))

        merged = merge(left, right, left_index=True,
                       right_index=True, copy=False)

        merged['a'] = 6
        self.assert_((left['a'] == 6).all())

        merged['d'] = 'peekaboo'
        self.assert_((right['d'] == 'peekaboo').all())

    def test_join_sort(self):
        left = DataFrame({'key' : ['foo', 'bar', 'baz', 'foo'],
                          'value' : [1, 2, 3, 4]})
        right = DataFrame({'value2' : ['a', 'b', 'c']},
                          index=['bar', 'baz', 'foo'])

        joined = left.join(right, on='key', sort=True)
        expected = DataFrame({'key' : ['bar', 'baz', 'foo', 'foo'],
                              'value' : [2, 3, 1, 4],
                              'value2' : ['a', 'b', 'c', 'c']},
                             index=[1, 2, 0, 3])
        assert_frame_equal(joined, expected)

        # smoke test
        joined = left.join(right, on='key', sort=False)
        self.assert_(np.array_equal(joined.index, range(4)))

    def test_intelligently_handle_join_key(self):
        # #733, be a bit more 1337 about not returning unconsolidated DataFrame

        left = DataFrame({'key' : [1, 1, 2, 2, 3],
                          'value' : range(5)}, columns=['value', 'key'])
        right = DataFrame({'key' : [1, 1, 2, 3, 4, 5],
                           'rvalue' : range(6)})

        joined = merge(left, right, on='key', how='outer')
        expected = DataFrame({'key' : [1, 1, 1, 1, 2, 2, 3, 4, 5.],
                              'value' : np.array([0, 0, 1, 1, 2, 3, 4,
                                                  np.nan, np.nan]),
                              'rvalue' : np.array([0, 1, 0, 1, 2, 2, 3, 4, 5])},
                             columns=['value', 'key', 'rvalue'])
        assert_frame_equal(joined, expected)

        self.assert_(joined._data.is_consolidated())

    def test_handle_join_key_pass_array(self):
        left = DataFrame({'key' : [1, 1, 2, 2, 3],
                          'value' : range(5)}, columns=['value', 'key'])
        right = DataFrame({'rvalue' : range(6)})
        key = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on='key', right_on=key, how='outer')
        merged2 = merge(right, left, left_on=key, right_on='key', how='outer')

        assert_series_equal(merged['key'], merged2['key'])
        self.assert_(merged['key'].notnull().all())
        self.assert_(merged2['key'].notnull().all())

        left = DataFrame({'value' : range(5)}, columns=['value'])
        right = DataFrame({'rvalue' : range(6)})
        lkey = np.array([1, 1, 2, 2, 3])
        rkey = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on=lkey, right_on=rkey, how='outer')
        self.assert_(np.array_equal(merged['key_0'],
                                    np.array([1, 1, 1, 1, 2, 2, 3, 4, 5])))

        left = DataFrame({'value': range(3)})
        right = DataFrame({'rvalue' : range(6)})

        key = np.array([0, 1, 1, 2, 2, 3])
        merged = merge(left, right, left_index=True, right_on=key, how='outer')
        self.assert_(np.array_equal(merged['key_0'], key))

class TestMergeMulti(unittest.TestCase):

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
        self.data = DataFrame({'key1' : key1, 'key2' : key2,
                               'data' : data})

    def test_merge_on_multikey(self):
        joined = self.data.join(self.to_join, on=['key1', 'key2'])

        join_key = Index(zip(self.data['key1'], self.data['key2']))
        indexer = self.to_join.index.get_indexer(join_key)
        ex_values = self.to_join.values.take(indexer, axis=0)
        ex_values[indexer == -1] = np.nan
        expected = self.data.join(DataFrame(ex_values,
                                            columns=self.to_join.columns))

        # TODO: columns aren't in the same order yet
        assert_frame_equal(joined, expected.ix[:, joined.columns])

    def test_merge_right_vs_left(self):
        # compare left vs right merge with multikey
        merged1 = self.data.merge(self.to_join, left_on=['key1', 'key2'],
                                  right_index=True, how='left')
        merged2 = self.to_join.merge(self.data, right_on=['key1', 'key2'],
                                     left_index=True, how='right')
        merged2 = merged2.ix[:, merged1.columns]
        assert_frame_equal(merged1, merged2)

    def test_compress_group_combinations(self):

        # ~ 40000000 possible unique groups
        key1 = np.array([rands(10) for _ in xrange(10000)], dtype='O')
        key1 = np.tile(key1, 2)
        key2 = key1[::-1]

        df = DataFrame({'key1' : key1, 'key2' : key2,
                        'value1' : np.random.randn(20000)})

        df2 = DataFrame({'key1' : key1[::2], 'key2' : key2[::2],
                         'value2' : np.random.randn(10000)})

        # just to hit the label compression code path
        merged = merge(df, df2, how='outer')

    def test_left_join_index_preserve_order(self):

        left = DataFrame({'k1' : [0, 1, 2] * 8,
                          'k2' : ['foo', 'bar'] * 12,
                          'v' : np.arange(24)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2' : [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected['v2'][(expected.k1 == 2) & (expected.k2 == 'bar')] = 5
        expected['v2'][(expected.k1 == 1) & (expected.k2 == 'foo')] = 7

        tm.assert_frame_equal(result, expected)

        # do a right join for an extra test
        joined = merge(right, left, left_index=True,
                       right_on=['k1', 'k2'], how='right')
        tm.assert_frame_equal(joined.ix[:, expected.columns], expected)

    def test_left_merge_na_buglet(self):
        left = DataFrame({'id': list('abcde'), 'v1': randn(5),
                          'v2': randn(5), 'dummy' : list('abcde'),
                          'v3' : randn(5)},
                         columns=['id', 'v1', 'v2', 'dummy', 'v3'])
        right = DataFrame({'id' : ['a', 'b', np.nan, np.nan, np.nan],
                           'sv3' : [1.234, 5.678, np.nan, np.nan, np.nan]})

        merged = merge(left, right, on='id', how='left')

        rdf = right.drop(['id'], axis=1)
        expected = left.join(rdf)
        tm.assert_frame_equal(merged, expected)

def _check_join(left, right, result, join_col, how='left',
                lsuffix='.x', rsuffix='.y'):

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
    NA_SENTINEL = -1234567 # drop_duplicates not so NA-friendly...

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

    for col, s in b_re.iteritems():
        a_re[col] = s
    return a_re.reindex(columns=result_columns)

class TestConcatenate(unittest.TestCase):

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
        self.assert_('A' in partial_appended)

        partial_appended = end_frame.append(begin_frame)
        self.assert_('A' in partial_appended)

        # mixed type handling
        appended = self.mixed_frame[:5].append(self.mixed_frame[5:])
        assert_frame_equal(appended, self.mixed_frame)

        # what to test here
        mixed_appended = self.mixed_frame[:5].append(self.frame[5:])
        mixed_appended2 = self.frame[:5].append(self.mixed_frame[5:])

        # all equal except 'foo' column
        assert_frame_equal(mixed_appended.reindex(columns=['A', 'B', 'C', 'D']),
                           mixed_appended2.reindex(columns=['A', 'B', 'C', 'D']))

        # append empty
        empty = DataFrame({})

        appended = self.frame.append(empty)
        assert_frame_equal(self.frame, appended)
        self.assert_(appended is not self.frame)

        appended = empty.append(self.frame)
        assert_frame_equal(self.frame, appended)
        self.assert_(appended is not self.frame)

        # overlap
        self.assertRaises(Exception, self.frame.append, self.frame)

    def test_append_records(self):
        arr1 = np.zeros((2,),dtype=('i4,f4,a10'))
        arr1[:] = [(1,2.,'Hello'),(2,3.,"World")]

        arr2 = np.zeros((3,),dtype=('i4,f4,a10'))
        arr2[:] = [(3, 4.,'foo'),
                   (5, 6.,"bar"),
                   (7., 8., 'baz')]

        df1 = DataFrame(arr1)
        df2 = DataFrame(arr2)

        result = df1.append(df2, ignore_index=True)
        expected = DataFrame(np.concatenate((arr1, arr2)))
        assert_frame_equal(result, expected)

    def test_append_different_columns(self):
        df = DataFrame({'bools' : np.random.randn(10) > 0,
                        'ints' : np.random.randint(0, 10, 10),
                        'floats' : np.random.randn(10),
                        'strings' : ['foo', 'bar'] * 5})

        a = df[:5].ix[:, ['bools', 'ints', 'floats']]
        b = df[5:].ix[:, ['strings', 'ints', 'floats']]

        appended = a.append(b)
        self.assert_(isnull(appended['strings'][0:4]).all())
        self.assert_(isnull(appended['bools'][5:]).all())

    def test_append_many(self):
        chunks = [self.frame[:5], self.frame[5:10],
                  self.frame[10:15], self.frame[15:]]

        result = chunks[0].append(chunks[1:])
        tm.assert_frame_equal(result, self.frame)

        chunks[-1]['foo'] = 'bar'
        result = chunks[0].append(chunks[1:])
        tm.assert_frame_equal(result.ix[:, self.frame.columns], self.frame)
        self.assert_((result['foo'][15:] == 'bar').all())
        self.assert_(result['foo'][:15].isnull().all())

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
        df = DataFrame(np.random.randn(8, 4), columns=['A','B','C','D'])
        df['key'] = ['foo', 'bar'] * 4
        df1 = df.ix[:, ['A', 'B']]
        df2 = df.ix[:, ['C', 'D']]
        df3 = df.ix[:, ['key']]

        result = df1.join([df2, df3])
        assert_frame_equal(result, df)

    def test_append_missing_column_proper_upcast(self):
        df1 = DataFrame({'A' : np.array([1,2, 3, 4], dtype='i8')})
        df2 = DataFrame({'B' : np.array([True,False, True, False],
                                        dtype=bool)})

        appended = df1.append(df2, ignore_index=True)
        self.assert_(appended['A'].dtype == 'f8')
        self.assert_(appended['B'].dtype == 'O')

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

        self.assert_(np.array_equal(result.columns.levels[0], level))
        self.assertEqual(result.columns.names[0], 'group_key')

    def test_concat_dataframe_keys_bug(self):
        t1 = DataFrame({'value': Series([1,2,3],
                       index=Index(['a', 'b', 'c'], name='id'))})
        t2 = DataFrame({'value': Series([7, 8],
                       index=Index(['a', 'b'], name = 'id'))})

        # it works
        result = concat([t1, t2], axis=1, keys=['t1', 't2'])
        self.assertEqual(list(result.columns), [('t1', 'value'),
                                                ('t2', 'value')])

    def test_concat_dict(self):
        frames = {'foo' : DataFrame(np.random.randn(4, 3)),
                  'bar' : DataFrame(np.random.randn(4, 3)),
                  'baz' : DataFrame(np.random.randn(4, 3)),
                  'qux' : DataFrame(np.random.randn(4, 3))}

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

    def test_concat_multiindex_with_keys(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        frame = DataFrame(np.random.randn(10, 3), index=index,
                          columns=Index(['A', 'B', 'C'], name='exp'))
        result = concat([frame, frame], keys=[0, 1], names=['iteration'])

        self.assertEqual(result.index.names, ['iteration'] + index.names)
        tm.assert_frame_equal(result.ix[0], frame)
        tm.assert_frame_equal(result.ix[1], frame)
        self.assertEqual(result.index.nlevels, 3)

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
        self.assertEqual(result.index.names, [None] * 3)

        # no levels
        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        names=['first', 'second'])
        self.assertEqual(result.index.names, ['first', 'second'] + [None])
        self.assert_(np.array_equal(result.index.levels[0], ['baz', 'foo']))

    def test_crossed_dtypes_weird_corner(self):
        columns = ['A', 'B', 'C', 'D']
        df1 = DataFrame({'A' : np.array([1, 2, 3, 4], dtype='f8'),
                         'B' : np.array([1, 2, 3, 4], dtype='i8'),
                         'C' : np.array([1, 2, 3, 4], dtype='f8'),
                         'D' : np.array([1, 2, 3, 4], dtype='i8')},
                        columns=columns)

        df2 = DataFrame({'A' : np.array([1, 2, 3, 4], dtype='i8'),
                         'B' : np.array([1, 2, 3, 4], dtype='f8'),
                         'C' : np.array([1, 2, 3, 4], dtype='i8'),
                         'D' : np.array([1, 2, 3, 4], dtype='f8')},
                        columns=columns)

        appended = df1.append(df2, ignore_index=True)
        expected = DataFrame(np.concatenate([df1.values, df2.values], axis=0),
                             columns=columns)
        tm.assert_frame_equal(appended, expected)

    def test_handle_empty_objects(self):
        df = DataFrame(np.random.randn(10, 4), columns=list('abcd'))

        baz = df[:5]
        baz['foo'] = 'bar'
        empty = df[5:5]

        frames = [baz, empty, empty, df[5:]]
        concatted = concat(frames, axis=0)

        expected = df.ix[:, ['a', 'b', 'c', 'd', 'foo']]
        expected['foo'] = expected['foo'].astype('O')
        expected['foo'][:5] = 'bar'

        tm.assert_frame_equal(concatted, expected)

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

        joined = p1.join(p2, lsuffix='_p1', rsuffix='_p2')
        p1_suf = p1.ix[['ItemB', 'ItemC']].add_suffix('_p1')
        p2_suf = p2.ix[['ItemB', 'ItemC']].add_suffix('_p2')
        no_overlap = panel.ix[['ItemA']]
        expected = p1_suf.join(p2_suf).join(no_overlap)
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
            data_dict.update(p.iterkv())

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

    def test_concat_series(self):
        ts = tm.makeTimeSeries()
        ts.name = 'foo'

        pieces = [ts[:5], ts[5:15], ts[15:]]

        result = concat(pieces)
        tm.assert_series_equal(result, ts)
        self.assertEqual(result.name, ts.name)

        result = concat(pieces, keys=[0, 1, 2])
        expected = ts.copy()

        exp_labels = [np.repeat([0, 1, 2], [len(x) for x in pieces]),
                      np.arange(len(ts))]
        exp_index = MultiIndex(levels=[[0, 1, 2], ts.index],
                               labels=exp_labels)
        expected.index = exp_index
        tm.assert_series_equal(result, expected)

        self.assertRaises(Exception, concat, pieces, axis=1)

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
        self.assertRaises(Exception, concat, [None, None])

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)


