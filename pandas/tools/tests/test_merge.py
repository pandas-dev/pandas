import nose
import unittest

from numpy.random import randn
import numpy as np
import random

from pandas import *
from pandas.tools.merge import merge
from pandas.util.testing import assert_frame_equal, assert_series_equal
import pandas._tseries as lib
import pandas.util.testing as tm

a_ = np.array

N = 50
NGROUPS = 8
JOIN_TYPES = ['inner', 'outer', 'left', 'right']


def get_test_data(ngroups=NGROUPS, n=N):
    unique_groups = range(ngroups)
    arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
                         dtype=object)

    random.shuffle(arr)
    return arr

class TestMerge(unittest.TestCase):

    def setUp(self):
        # aggregate multiple columns
        self.df = DataFrame({'key1' : get_test_data(),
                             'key2' : get_test_data(),
                             'data1' : np.random.randn(N),
                             'data2' : np.random.randn(N)})

        # exclude a couple keys for fun
        self.df = self.df[self.df['key2'] > 1]

        self.df2 = DataFrame({'key1'  : get_test_data(n=N//5),
                              'key2'  : get_test_data(ngroups=NGROUPS//2,
                                                      n=N//5),
                              'value' : np.random.randn(N // 5)})

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

    def test_cython_full_outer_join(self):
        pass

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

    def test_merge_common(self):
        joined = merge(self.df, self.df2)
        exp = merge(self.df, self.df2, on=['key1', 'key2'])
        tm.assert_frame_equal(joined, exp)

    def test_merge_index(self):
        pass

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


class TestDataFrameJoin(unittest.TestCase):

    def setUp(self):
        index, data = tm.getMixedTypeDict()
        self.target = DataFrame(data, index=index)

        # Join on string value
        self.source = DataFrame({'MergedA' : data['A'], 'MergedD' : data['D']},
                                index=data['C'])

    def test_join_on(self):
        target = self.target
        source = self.source

        merged = target.join(source, on='C')
        self.assert_(np.array_equal(merged['MergedA'], target['A']))
        self.assert_(np.array_equal(merged['MergedD'], target['D']))

        # join with duplicates (fix regression from DataFrame/Matrix merge)
        df = DataFrame({'key' : ['a', 'a', 'b', 'b', 'c']})
        df2 = DataFrame({'value' : [0, 1, 2]}, index=['a', 'b', 'c'])
        joined = df.join(df2, on='key')
        expected = DataFrame({'key' : ['a', 'a', 'b', 'b', 'c'],
                              'value' : [0, 0, 1, 1, 2]})
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
        df = DataFrame({'key' : ['a', 'a', 'd', 'b', 'b', 'c']})
        df2 = DataFrame({'value' : [0, 1]}, index=['a', 'b'])

        joined = df.join(df2, on='key', how='inner')

        expected = df.join(df2, on='key')
        expected = expected[expected['value'].notnull()]
        self.assert_(np.array_equal(joined['key'], expected['key']))
        self.assert_(np.array_equal(joined['value'], expected['value']))
        self.assert_(joined.index.equals(expected.index))

    def test_join_on_singlekey_list(self):
        df = DataFrame({'key' : ['a', 'a', 'b', 'b', 'c']})
        df2 = DataFrame({'value' : [0, 1, 2]}, index=['a', 'b', 'c'])

        # corner cases
        joined = df.join(df2, on=['key'])
        expected = df.join(df2, on='key')

        assert_frame_equal(joined, expected)

    def test_join_on_multikey(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        to_join = DataFrame(np.random.randn(10, 3), index=index,
                            columns=['j_one', 'j_two', 'j_three'])

        # a little relevant example with NAs
        key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
                'qux', 'snap']
        key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
                'three', 'one']

        data = np.random.randn(len(key1))
        data = DataFrame({'key1' : key1, 'key2' : key2,
                          'data' : data})

        joined = data.join(to_join, on=['key1', 'key2'])

        join_key = Index(zip(key1, key2))
        indexer = to_join.index.get_indexer(join_key)
        ex_values = to_join.values.take(indexer, axis=0)
        ex_values[indexer == -1] = np.nan
        expected = data.join(DataFrame(ex_values, columns=to_join.columns))

        # TODO: columns aren't in the same order yet
        assert_frame_equal(joined, expected.ix[:, joined.columns])

    def test_join_on_series(self):
        result = self.target.join(self.source['MergedA'], on='C')
        expected = self.target.join(self.source[['MergedA']], on='C')
        assert_frame_equal(result, expected)

    def test_join_index_mixed(self):

        df1 = DataFrame({'A' : 1., 'B' : 2, 'C' : 'foo', 'D' : True},
                        index=np.arange(10),
                        columns=['A', 'B', 'C', 'D'])
        self.assert_(df1['B'].dtype == np.int64)
        self.assert_(df1['D'].dtype == np.bool_)

        df2 = DataFrame({'A' : 1., 'B' : 2, 'C' : 'foo', 'D' : True},
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
        d = DataFrame(randn(30,1), columns=['d'])

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

    def test_join_float64_float32(self):
        a = DataFrame(randn(10,2), columns=['a','b'])
        b = DataFrame(randn(10,1), columns=['c']).astype(np.float32)
        joined = a.join(b)
        expected = a.join(b.astype('f8'))
        assert_frame_equal(joined, expected)

def _join_by_hand(a, b, how='left'):
    join_index = a.index.join(b.index, how=how)

    a_re = a.reindex(join_index)
    b_re = b.reindex(join_index)

    result_columns = a.columns.append(b.columns)

    for col, s in b_re.iteritems():
        a_re[col] = s
    return a_re.reindex(columns=result_columns)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)



