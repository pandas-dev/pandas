# -*- coding: utf-8 -*-

from __future__ import print_function

import pandas as pd

from pandas.tests.frame.common import TestData

import pandas.util.testing as tm


class TestDataFrameJoin(TestData):

    def test_join(self):
        df1 = pd.DataFrame({'a': [20, 10, 0]}, index=[2, 1, 0])
        df2 = pd.DataFrame({'b': [100, 200, 300]}, index=[1, 2, 3])

        result = df1.join(df2)
        expected = pd.DataFrame({'a': [20, 10, 0], 'b': [200, 100, None]},
                                index=[2, 1, 0])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='left')
        expected = pd.DataFrame({'a': [20, 10, 0], 'b': [200, 100, None]},
                                index=[2, 1, 0])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='right')
        expected = pd.DataFrame({'a': [10, 20, None], 'b': [100, 200, 300]},
                                index=[1, 2, 3])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='inner')
        expected = pd.DataFrame({'a': [20, 10], 'b': [200, 100]},
                                index=[2, 1])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='outer')
        expected = pd.DataFrame({'a': [0, 10, 20, None],
                                 'b': [None, 100, 200, 300]},
                                index=[0, 1, 2, 3])
        tm.assert_frame_equal(result, expected)

    def test_join_sort(self):
        df1 = pd.DataFrame({'a': [20, 10, 0]}, index=[2, 1, 0])
        df2 = pd.DataFrame({'b': [100, 200, 300]}, index=[1, 2, 3])

        result = df1.join(df2, sort=True)
        expected = pd.DataFrame({'a': [0, 10, 20], 'b': [None, 100, 200]},
                                index=[0, 1, 2])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='left', sort=True)
        expected = pd.DataFrame({'a': [0, 10, 20], 'b': [None, 100, 200]},
                                index=[0, 1, 2])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='right', sort=True)
        expected = pd.DataFrame({'a': [10, 20, None], 'b': [100, 200, 300]},
                                index=[1, 2, 3])
        tm.assert_frame_equal(result, expected)

        result = df2.join(df1, how='right', sort=True)
        expected = pd.DataFrame([[None, 0], [100, 10], [200, 20]],
                                columns=['b', 'a'], index=[0, 1, 2])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='inner', sort=True)
        expected = pd.DataFrame({'a': [10, 20], 'b': [100, 200]},
                                index=[1, 2])
        tm.assert_frame_equal(result, expected)

        result = df1.join(df2, how='outer', sort=True)
        expected = pd.DataFrame({'a': [0, 10, 20, None],
                                 'b': [None, 100, 200, 300]},
                                index=[0, 1, 2, 3])
        tm.assert_frame_equal(result, expected)

    def test_join_index(self):
        # left / right

        f = self.frame.loc[self.frame.index[:10], ['A', 'B']]
        f2 = self.frame.loc[self.frame.index[5:], ['C', 'D']].iloc[::-1]

        joined = f.join(f2)
        tm.assert_index_equal(f.index, joined.index)
        expected_columns = pd.Index(['A', 'B', 'C', 'D'])
        tm.assert_index_equal(joined.columns, expected_columns)

        joined = f.join(f2, how='left')
        tm.assert_index_equal(joined.index, f.index)
        tm.assert_index_equal(joined.columns, expected_columns)

        joined = f.join(f2, how='right')
        tm.assert_index_equal(joined.index, f2.index)
        tm.assert_index_equal(joined.columns, expected_columns)

        # inner

        joined = f.join(f2, how='inner')
        tm.assert_index_equal(joined.index, f.index[5:10])
        tm.assert_index_equal(joined.columns, expected_columns)

        # outer

        joined = f.join(f2, how='outer')
        tm.assert_index_equal(joined.index, self.frame.index.sort_values())
        tm.assert_index_equal(joined.columns, expected_columns)

        tm.assertRaisesRegexp(ValueError, 'join method', f.join, f2, how='foo')

        # corner case - overlapping columns
        for how in ('outer', 'left', 'inner'):
            with tm.assertRaisesRegexp(ValueError, 'columns overlap but '
                                       'no suffix'):
                self.frame.join(self.frame, how=how)

    def test_join_index_more(self):
        af = self.frame.loc[:, ['A', 'B']]
        bf = self.frame.loc[::2, ['C', 'D']]

        expected = af.copy()
        expected['C'] = self.frame['C'][::2]
        expected['D'] = self.frame['D'][::2]

        result = af.join(bf)
        tm.assert_frame_equal(result, expected)

        result = af.join(bf, how='right')
        tm.assert_frame_equal(result, expected[::2])

        result = bf.join(af, how='right')
        tm.assert_frame_equal(result, expected.loc[:, result.columns])

    def test_join_index_series(self):
        df = self.frame.copy()
        s = df.pop(self.frame.columns[-1])
        joined = df.join(s)

        # TODO should this check_names ?
        tm.assert_frame_equal(joined, self.frame, check_names=False)

        s.name = None
        tm.assertRaisesRegexp(ValueError, 'must have a name', df.join, s)

    def test_join_overlap(self):
        df1 = self.frame.loc[:, ['A', 'B', 'C']]
        df2 = self.frame.loc[:, ['B', 'C', 'D']]

        joined = df1.join(df2, lsuffix='_df1', rsuffix='_df2')
        df1_suf = df1.loc[:, ['B', 'C']].add_suffix('_df1')
        df2_suf = df2.loc[:, ['B', 'C']].add_suffix('_df2')

        no_overlap = self.frame.loc[:, ['A', 'D']]
        expected = df1_suf.join(df2_suf).join(no_overlap)

        # column order not necessarily sorted
        tm.assert_frame_equal(joined, expected.loc[:, joined.columns])
