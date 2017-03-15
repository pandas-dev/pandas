# -*- coding: utf-8 -*-

from __future__ import print_function

import pandas as pd

from pandas.tests.frame.common import TestData

import pandas.util.testing as tm
from pandas.util.testing import (assertRaisesRegexp,
                                 assert_frame_equal)


class TestDataFrameJoin(tm.TestCase, TestData):

    def test_join(self):
        df1 = pd.DataFrame({'a': [20, 10, 0]}, index=[2, 1, 0])
        df2 = pd.DataFrame({'b': [100, 200, 300]}, index=[1, 2, 3])

        result = df1.join(df2, how='inner')
        expected = pd.DataFrame({'a': [20, 10], 'b': [200, 100]}, index=[2, 1])
        self.assert_frame_equal(result, expected)

    def test_join_index(self):
        # left / right

        f = self.frame.loc[self.frame.index[:10], ['A', 'B']]
        f2 = self.frame.loc[self.frame.index[5:], ['C', 'D']].iloc[::-1]

        joined = f.join(f2)
        self.assert_index_equal(f.index, joined.index)
        expected_columns = pd.Index(['A', 'B', 'C', 'D'])
        self.assert_index_equal(joined.columns, expected_columns)

        joined = f.join(f2, how='left')
        self.assert_index_equal(joined.index, f.index)
        self.assert_index_equal(joined.columns, expected_columns)

        joined = f.join(f2, how='right')
        self.assert_index_equal(joined.index, f2.index)
        self.assert_index_equal(joined.columns, expected_columns)

        # inner

        joined = f.join(f2, how='inner')
        self.assert_index_equal(joined.index, f.index[5:10])
        self.assert_index_equal(joined.columns, expected_columns)

        # outer

        joined = f.join(f2, how='outer')
        self.assert_index_equal(joined.index, self.frame.index.sort_values())
        self.assert_index_equal(joined.columns, expected_columns)

        assertRaisesRegexp(ValueError, 'join method', f.join, f2, how='foo')

        # corner case - overlapping columns
        for how in ('outer', 'left', 'inner'):
            with assertRaisesRegexp(ValueError, 'columns overlap but '
                                    'no suffix'):
                self.frame.join(self.frame, how=how)

    def test_join_index_more(self):
        af = self.frame.loc[:, ['A', 'B']]
        bf = self.frame.loc[::2, ['C', 'D']]

        expected = af.copy()
        expected['C'] = self.frame['C'][::2]
        expected['D'] = self.frame['D'][::2]

        result = af.join(bf)
        assert_frame_equal(result, expected)

        result = af.join(bf, how='right')
        assert_frame_equal(result, expected[::2])

        result = bf.join(af, how='right')
        assert_frame_equal(result, expected.loc[:, result.columns])

    def test_join_index_series(self):
        df = self.frame.copy()
        s = df.pop(self.frame.columns[-1])
        joined = df.join(s)

        # TODO should this check_names ?
        assert_frame_equal(joined, self.frame, check_names=False)

        s.name = None
        assertRaisesRegexp(ValueError, 'must have a name', df.join, s)

    def test_join_overlap(self):
        df1 = self.frame.loc[:, ['A', 'B', 'C']]
        df2 = self.frame.loc[:, ['B', 'C', 'D']]

        joined = df1.join(df2, lsuffix='_df1', rsuffix='_df2')
        df1_suf = df1.loc[:, ['B', 'C']].add_suffix('_df1')
        df2_suf = df2.loc[:, ['B', 'C']].add_suffix('_df2')

        no_overlap = self.frame.loc[:, ['A', 'D']]
        expected = df1_suf.join(df2_suf).join(no_overlap)

        # column order not necessarily sorted
        assert_frame_equal(joined, expected.loc[:, joined.columns])