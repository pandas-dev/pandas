# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pandas.util.testing as tm


class TestSparseGroupBy(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.dense = pd.DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                                         'foo', 'bar', 'foo', 'foo'],
                                   'B': ['one', 'one', 'two', 'three',
                                         'two', 'two', 'one', 'three'],
                                   'C': np.random.randn(8),
                                   'D': np.random.randn(8),
                                   'E': [np.nan, np.nan, 1, 2,
                                         np.nan, 1, np.nan, np.nan]})
        self.sparse = self.dense.to_sparse()

    def test_first_last_nth(self):
        # tests for first / last / nth
        sparse_grouped = self.sparse.groupby('A')
        dense_grouped = self.dense.groupby('A')

        tm.assert_frame_equal(sparse_grouped.first(),
                              dense_grouped.first())
        tm.assert_frame_equal(sparse_grouped.last(),
                              dense_grouped.last())
        tm.assert_frame_equal(sparse_grouped.nth(1),
                              dense_grouped.nth(1))

    def test_aggfuncs(self):
        sparse_grouped = self.sparse.groupby('A')
        dense_grouped = self.dense.groupby('A')

        tm.assert_frame_equal(sparse_grouped.mean(),
                              dense_grouped.mean())

        # ToDo: sparse sum includes str column
        # tm.assert_frame_equal(sparse_grouped.sum(),
        #                       dense_grouped.sum())

        tm.assert_frame_equal(sparse_grouped.count(),
                              dense_grouped.count())
