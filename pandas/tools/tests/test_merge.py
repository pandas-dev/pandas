import nose
import unittest

import numpy as np
import random

from pandas import *
from pandas.tools.merge import merge
import pandas._sandbox as sbx

a_ = np.array

N = 100
NGROUPS = 8

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

        ls, rs = sbx.left_outer_join(left, right, max_group)

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

        rs, ls  = sbx.left_outer_join(right, left, max_group)

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
        raise nose.SkipTest

        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1, 4], dtype='i4')
        max_group = 5

        ls, rs = sbx.inner_join(left, right, max_group)

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

    # def test_full_outer_join(self):
    #     joined_key2 = merge(self.df, self.df2, on='key2', how='outer')
    #     _check_join(self.df, self.df2, joined_key2, ['key2'], how='outer')

    #     joined_both = merge(self.df, self.df2, how='outer')
    #     _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
    #                 how='outer')

    def test_handle_overlap(self):
        pass

    def test_merge_common(self):
        pass

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
            if how == 'left':
                raise AssertionError('key %s should not have been in the join'
                                     % str(group_key))

            _assert_all_na(l_joined, left.columns, join_col)
        else:
            _assert_same_contents(l_joined, lgroup)

        try:
            rgroup = right_grouped.get_group(group_key)
        except KeyError:
            if how == 'right':
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

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)



