# -*- coding: utf-8 -*-

import pytest


import pandas as pd

from pandas import (MultiIndex, date_range)
from pandas.compat import lrange
from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike

import pandas.util.testing as tm

from .common import Base


class TestFromProduct(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_from_product(self):

        first = ['foo', 'bar', 'buz']
        second = ['a', 'b', 'c']
        names = ['first', 'second']
        result = MultiIndex.from_product([first, second], names=names)

        tuples = [('foo', 'a'), ('foo', 'b'), ('foo', 'c'), ('bar', 'a'),
                  ('bar', 'b'), ('bar', 'c'), ('buz', 'a'), ('buz', 'b'),
                  ('buz', 'c')]
        expected = MultiIndex.from_tuples(tuples, names=names)

        tm.assert_index_equal(result, expected)

    def test_from_product_iterator(self):
        # GH 18434
        first = ['foo', 'bar', 'buz']
        second = ['a', 'b', 'c']
        names = ['first', 'second']
        tuples = [('foo', 'a'), ('foo', 'b'), ('foo', 'c'), ('bar', 'a'),
                  ('bar', 'b'), ('bar', 'c'), ('buz', 'a'), ('buz', 'b'),
                  ('buz', 'c')]
        expected = MultiIndex.from_tuples(tuples, names=names)

        # iterator as input
        result = MultiIndex.from_product(iter([first, second]), names=names)
        tm.assert_index_equal(result, expected)

        # Invalid non-iterable input
        with tm.assert_raises_regex(
                TypeError, "Input must be a list / sequence of iterables."):
            MultiIndex.from_product(0)

    def test_from_product_empty(self):
        # 0 levels
        with tm.assert_raises_regex(
                ValueError, "Must pass non-zero number of levels/labels"):
            MultiIndex.from_product([])

        # 1 level
        result = MultiIndex.from_product([[]], names=['A'])
        expected = pd.Index([], name='A')
        tm.assert_index_equal(result.levels[0], expected)

        # 2 levels
        l1 = [[], ['foo', 'bar', 'baz'], []]
        l2 = [[], [], ['a', 'b', 'c']]
        names = ['A', 'B']
        for first, second in zip(l1, l2):
            result = MultiIndex.from_product([first, second], names=names)
            expected = MultiIndex(levels=[first, second],
                                  labels=[[], []], names=names)
            tm.assert_index_equal(result, expected)

        # GH12258
        names = ['A', 'B', 'C']
        for N in range(4):
            lvl2 = lrange(N)
            result = MultiIndex.from_product([[], lvl2, []], names=names)
            expected = MultiIndex(levels=[[], lvl2, []],
                                  labels=[[], [], []], names=names)
            tm.assert_index_equal(result, expected)

    def test_from_product_invalid_input(self):
        invalid_inputs = [1, [1], [1, 2], [[1], 2],
                          'a', ['a'], ['a', 'b'], [['a'], 'b']]
        for i in invalid_inputs:
            pytest.raises(TypeError, MultiIndex.from_product, iterables=i)

    def test_from_product_datetimeindex(self):
        dt_index = date_range('2000-01-01', periods=2)
        mi = pd.MultiIndex.from_product([[1, 2], dt_index])
        etalon = construct_1d_object_array_from_listlike([(1, pd.Timestamp(
            '2000-01-01')), (1, pd.Timestamp('2000-01-02')), (2, pd.Timestamp(
                '2000-01-01')), (2, pd.Timestamp('2000-01-02'))])
        tm.assert_numpy_array_equal(mi.values, etalon)

    def test_from_product_index_series_categorical(self):
        # GH13743
        first = ['foo', 'bar']
        for ordered in [False, True]:
            idx = pd.CategoricalIndex(list("abcaab"), categories=list("bac"),
                                      ordered=ordered)
            expected = pd.CategoricalIndex(list("abcaab") + list("abcaab"),
                                           categories=list("bac"),
                                           ordered=ordered)

            for arr in [idx, pd.Series(idx), idx.values]:
                result = pd.MultiIndex.from_product([first, arr])
                tm.assert_index_equal(result.get_level_values(1), expected)
