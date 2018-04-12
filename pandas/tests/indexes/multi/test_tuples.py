# -*- coding: utf-8 -*-

from pandas import MultiIndex

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base

class TestTuples(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_from_tuples(self):
        tm.assert_raises_regex(TypeError, 'Cannot infer number of levels '
                                          'from empty list',
                               MultiIndex.from_tuples, [])

        expected = MultiIndex(levels=[[1, 3], [2, 4]],
                              labels=[[0, 1], [0, 1]],
                              names=['a', 'b'])

        # input tuples
        result = MultiIndex.from_tuples(((1, 2), (3, 4)), names=['a', 'b'])
        tm.assert_index_equal(result, expected)

    def test_from_tuples_iterator(self):
        # GH 18434
        # input iterator for tuples
        expected = MultiIndex(levels=[[1, 3], [2, 4]],
                              labels=[[0, 1], [0, 1]],
                              names=['a', 'b'])

        result = MultiIndex.from_tuples(zip([1, 3], [2, 4]), names=['a', 'b'])
        tm.assert_index_equal(result, expected)

        # input non-iterables
        with tm.assert_raises_regex(
                TypeError, 'Input must be a list / sequence of tuple-likes.'):
            MultiIndex.from_tuples(0)

    def test_from_tuples_empty(self):
        # GH 16777
        result = MultiIndex.from_tuples([], names=['a', 'b'])
        expected = MultiIndex.from_arrays(arrays=[[], []],
                                          names=['a', 'b'])
        tm.assert_index_equal(result, expected)
