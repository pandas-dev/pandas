# -*- coding: utf-8 -*-
import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm


class TestSparseGroupBy(object):

    def setup_method(self, method):
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

        # TODO: shouldn't these all be spares or not?
        tm.assert_frame_equal(sparse_grouped.first(),
                              dense_grouped.first())
        tm.assert_frame_equal(sparse_grouped.last(),
                              dense_grouped.last())
        tm.assert_frame_equal(sparse_grouped.nth(1),
                              dense_grouped.nth(1).to_sparse())

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


@pytest.mark.parametrize("fill_value", [0, np.nan])
def test_groupby_includes_fill_value(fill_value):
    # https://github.com/pandas-dev/pandas/issues/5078
    df = pd.DataFrame({'a': [fill_value, 1, fill_value, fill_value],
                       'b': [fill_value, 1, fill_value, fill_value]})
    sdf = df.to_sparse(fill_value=fill_value)
    result = sdf.groupby('a').sum()
    expected = df.groupby('a').sum()
    tm.assert_frame_equal(result, expected,
                          check_index_type=False)
