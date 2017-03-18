# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import pandas as pd
import pandas.util.testing as tm


class TestDataFrameMerge(object):

    def test_merge_on_indexes(self):
        df1 = pd.DataFrame({'a': [20, 10, 0]}, index=[2, 1, 0])
        df2 = pd.DataFrame({'b': [100, 200, 300]}, index=[1, 2, 3])

        # default how='inner'
        result = df1.merge(df2, left_index=True, right_index=True)
        expected = pd.DataFrame({'a': [20, 10], 'b': [200, 100]},
                                index=[2, 1])
        tm.assert_frame_equal(result, expected)

        # how='left'
        result = df1.merge(df2, left_index=True, right_index=True, how='left')
        expected = pd.DataFrame({'a': [20, 10, 0], 'b': [200, 100, np.nan]},
                                index=[2, 1, 0])
        tm.assert_frame_equal(result, expected)

        # how='right'
        result = df1.merge(df2, left_index=True, right_index=True, how='right')
        expected = pd.DataFrame({'a': [10, 20, np.nan], 'b': [100, 200, 300]},
                                index=[1, 2, 3])
        tm.assert_frame_equal(result, expected)

        # how='inner'
        result = df1.merge(df2, left_index=True, right_index=True, how='inner')
        expected = pd.DataFrame({'a': [20, 10], 'b': [200, 100]},
                                index=[2, 1])
        tm.assert_frame_equal(result, expected)

        # how='outer'
        result = df1.merge(df2, left_index=True, right_index=True, how='outer')
        expected = pd.DataFrame({'a': [0, 10, 20, np.nan],
                                 'b': [np.nan, 100, 200, 300]},
                                index=[0, 1, 2, 3])
        tm.assert_frame_equal(result, expected)

    def test_merge_on_indexes_sort(self):
        df1 = pd.DataFrame({'a': [20, 10, 0]}, index=[2, 1, 0])
        df2 = pd.DataFrame({'b': [100, 200, 300]}, index=[1, 2, 3])

        # default how='inner'
        result = df1.merge(df2, left_index=True, right_index=True, sort=True)
        expected = pd.DataFrame({'a': [10, 20], 'b': [100, 200]},
                                index=[1, 2])
        tm.assert_frame_equal(result, expected)

        # how='left'
        result = df1.merge(df2, left_index=True, right_index=True, how='left',
                           sort=True)
        expected = pd.DataFrame({'a': [0, 10, 20], 'b': [np.nan, 100, 200]},
                                index=[0, 1, 2])
        tm.assert_frame_equal(result, expected)

        # how='right' (already sorted)
        result = df1.merge(df2, left_index=True, right_index=True,
                           how='right', sort=True)
        expected = pd.DataFrame({'a': [10, 20, np.nan], 'b': [100, 200, 300]},
                                index=[1, 2, 3])
        tm.assert_frame_equal(result, expected)

        # how='right'
        result = df2.merge(df1, left_index=True, right_index=True,
                           how='right', sort=True)
        expected = pd.DataFrame([[np.nan, 0], [100, 10], [200, 20]],
                                columns=['b', 'a'], index=[0, 1, 2])
        tm.assert_frame_equal(result, expected)

        # how='inner'
        result = df1.merge(df2, left_index=True, right_index=True,
                           how='inner', sort=True)
        expected = pd.DataFrame({'a': [10, 20], 'b': [100, 200]},
                                index=[1, 2])
        tm.assert_frame_equal(result, expected)

        # how='outer'
        result = df1.merge(df2, left_index=True, right_index=True,
                           how='outer', sort=True)
        expected = pd.DataFrame({'a': [0, 10, 20, np.nan],
                                 'b': [np.nan, 100, 200, 300]},
                                index=[0, 1, 2, 3])
        tm.assert_frame_equal(result, expected)
