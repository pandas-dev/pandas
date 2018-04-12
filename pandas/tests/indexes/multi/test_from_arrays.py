# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas as pd

from pandas import Index, MultiIndex
from pandas._libs.tslib import Timestamp

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base

class TestFromArrays(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_from_arrays(self):
        arrays = []
        for lev, lab in zip(self.index.levels, self.index.labels):
            arrays.append(np.asarray(lev).take(lab))

        # list of arrays as input
        result = MultiIndex.from_arrays(arrays, names=self.index.names)
        tm.assert_index_equal(result, self.index)

        # infer correctly
        result = MultiIndex.from_arrays([[pd.NaT, Timestamp('20130101')],
                                         ['a', 'b']])
        assert result.levels[0].equals(Index([Timestamp('20130101')]))
        assert result.levels[1].equals(Index(['a', 'b']))

    def test_from_arrays_iterator(self):
        # GH 18434
        arrays = []
        for lev, lab in zip(self.index.levels, self.index.labels):
            arrays.append(np.asarray(lev).take(lab))

        # iterator as input
        result = MultiIndex.from_arrays(iter(arrays), names=self.index.names)
        tm.assert_index_equal(result, self.index)

        # invalid iterator input
        with tm.assert_raises_regex(
                TypeError, "Input must be a list / sequence of array-likes."):
            MultiIndex.from_arrays(0)

    def test_from_arrays_index_series_datetimetz(self):
        idx1 = pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                             tz='US/Eastern')
        idx2 = pd.date_range('2015-01-01 10:00', freq='H', periods=3,
                             tz='Asia/Tokyo')
        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_series_timedelta(self):
        idx1 = pd.timedelta_range('1 days', freq='D', periods=3)
        idx2 = pd.timedelta_range('2 hours', freq='H', periods=3)
        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_series_period(self):
        idx1 = pd.period_range('2011-01-01', freq='D', periods=3)
        idx2 = pd.period_range('2015-01-01', freq='H', periods=3)
        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_datetimelike_mixed(self):
        idx1 = pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                             tz='US/Eastern')
        idx2 = pd.date_range('2015-01-01 10:00', freq='H', periods=3)
        idx3 = pd.timedelta_range('1 days', freq='D', periods=3)
        idx4 = pd.period_range('2011-01-01', freq='D', periods=3)

        result = pd.MultiIndex.from_arrays([idx1, idx2, idx3, idx4])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)
        tm.assert_index_equal(result.get_level_values(2), idx3)
        tm.assert_index_equal(result.get_level_values(3), idx4)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1),
                                             pd.Series(idx2),
                                             pd.Series(idx3),
                                             pd.Series(idx4)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)
        tm.assert_index_equal(result2.get_level_values(2), idx3)
        tm.assert_index_equal(result2.get_level_values(3), idx4)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_series_categorical(self):
        # GH13743
        idx1 = pd.CategoricalIndex(list("abcaab"), categories=list("bac"),
                                   ordered=False)
        idx2 = pd.CategoricalIndex(list("abcaab"), categories=list("bac"),
                                   ordered=True)

        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        result3 = pd.MultiIndex.from_arrays([idx1.values, idx2.values])
        tm.assert_index_equal(result3.get_level_values(0), idx1)
        tm.assert_index_equal(result3.get_level_values(1), idx2)

    def test_from_arrays_empty(self):
        # 0 levels
        with tm.assert_raises_regex(
                ValueError, "Must pass non-zero number of levels/labels"):
            MultiIndex.from_arrays(arrays=[])

        # 1 level
        result = MultiIndex.from_arrays(arrays=[[]], names=['A'])
        assert isinstance(result, MultiIndex)
        expected = Index([], name='A')
        tm.assert_index_equal(result.levels[0], expected)

        # N levels
        for N in [2, 3]:
            arrays = [[]] * N
            names = list('ABC')[:N]
            result = MultiIndex.from_arrays(arrays=arrays, names=names)
            expected = MultiIndex(levels=[[]] * N, labels=[[]] * N,
                                  names=names)
            tm.assert_index_equal(result, expected)

    def test_from_arrays_invalid_input(self):
        invalid_inputs = [1, [1], [1, 2], [[1], 2],
                          'a', ['a'], ['a', 'b'], [['a'], 'b']]
        for i in invalid_inputs:
            pytest.raises(TypeError, MultiIndex.from_arrays, arrays=i)

    def test_from_arrays_different_lengths(self):
        # see gh-13599
        idx1 = [1, 2, 3]
        idx2 = ['a', 'b']
        tm.assert_raises_regex(ValueError, '^all arrays must '
                                           'be same length$',
                               MultiIndex.from_arrays, [idx1, idx2])

        idx1 = []
        idx2 = ['a', 'b']
        tm.assert_raises_regex(ValueError, '^all arrays must '
                                           'be same length$',
                               MultiIndex.from_arrays, [idx1, idx2])

        idx1 = [1, 2, 3]
        idx2 = []
        tm.assert_raises_regex(ValueError, '^all arrays must '
                                           'be same length$',
                               MultiIndex.from_arrays, [idx1, idx2])
