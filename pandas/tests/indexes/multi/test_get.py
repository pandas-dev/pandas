# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas as pd

from pandas import (CategoricalIndex, Index, MultiIndex)
from pandas.compat import lrange
from pandas.core.indexes.base import InvalidIndexError

import pandas.util.testing as tm

from pandas.util.testing import assert_almost_equal

from pandas.tests.indexes.common import Base


class TestGet(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_getitem(self):
        # scalar
        assert self.index[2] == ('bar', 'one')

        # slice
        result = self.index[2:5]
        expected = self.index[[2, 3, 4]]
        assert result.equals(expected)

        # boolean
        result = self.index[[True, False, True, False, True, True]]
        result2 = self.index[np.array([True, False, True, False, True, True])]
        expected = self.index[[0, 2, 4, 5]]
        assert result.equals(expected)
        assert result2.equals(expected)

    def test_getitem_group_select(self):
        sorted_idx, _ = self.index.sortlevel(0)
        assert sorted_idx.get_loc('baz') == slice(3, 4)
        assert sorted_idx.get_loc('foo') == slice(0, 2)

    def test_get_loc(self):
        assert self.index.get_loc(('foo', 'two')) == 1
        assert self.index.get_loc(('baz', 'two')) == 3
        pytest.raises(KeyError, self.index.get_loc, ('bar', 'two'))
        pytest.raises(KeyError, self.index.get_loc, 'quux')

        pytest.raises(NotImplementedError, self.index.get_loc, 'foo',
                      method='nearest')

        # 3 levels
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])
        pytest.raises(KeyError, index.get_loc, (1, 1))
        assert index.get_loc((2, 0)) == slice(3, 5)

    def test_get_loc_duplicates(self):
        index = Index([2, 2, 2, 2])
        result = index.get_loc(2)
        expected = slice(0, 4)
        assert result == expected
        # pytest.raises(Exception, index.get_loc, 2)

        index = Index(['c', 'a', 'a', 'b', 'b'])
        rs = index.get_loc('c')
        xp = 0
        assert rs == xp

    def test_get_value_duplicates(self):
        index = MultiIndex(levels=[['D', 'B', 'C'],
                                   [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                                   [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])

        assert index.get_loc('D') == slice(0, 3)
        with pytest.raises(KeyError):
            index._engine.get_value(np.array([]), 'D')

    def test_get_loc_level(self):
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        loc, new_index = index.get_loc_level((0, 1))
        expected = slice(1, 2)
        exp_index = index[expected].droplevel(0).droplevel(0)
        assert loc == expected
        assert new_index.equals(exp_index)

        loc, new_index = index.get_loc_level((0, 1, 0))
        expected = 1
        assert loc == expected
        assert new_index is None

        pytest.raises(KeyError, index.get_loc_level, (2, 2))

        index = MultiIndex(levels=[[2000], lrange(4)], labels=[np.array(
            [0, 0, 0, 0]), np.array([0, 1, 2, 3])])
        result, new_index = index.get_loc_level((2000, slice(None, None)))
        expected = slice(None, None)
        assert result == expected
        assert new_index.equals(index.droplevel(0))

    @pytest.mark.parametrize('level', [0, 1])
    @pytest.mark.parametrize('null_val', [np.nan, pd.NaT, None])
    def test_get_loc_nan(self, level, null_val):
        # GH 18485 : NaN in MultiIndex
        levels = [['a', 'b'], ['c', 'd']]
        key = ['b', 'd']
        levels[level] = np.array([0, null_val], dtype=type(null_val))
        key[level] = null_val
        idx = MultiIndex.from_product(levels)
        assert idx.get_loc(tuple(key)) == 3

    def test_get_loc_missing_nan(self):
        # GH 8569
        idx = MultiIndex.from_arrays([[1.0, 2.0], [3.0, 4.0]])
        assert isinstance(idx.get_loc(1), slice)
        pytest.raises(KeyError, idx.get_loc, 3)
        pytest.raises(KeyError, idx.get_loc, np.nan)
        pytest.raises(KeyError, idx.get_loc, [np.nan])

    @pytest.mark.parametrize('dtype1', [int, float, bool, str])
    @pytest.mark.parametrize('dtype2', [int, float, bool, str])
    def test_get_loc_multiple_dtypes(self, dtype1, dtype2):
        # GH 18520
        levels = [np.array([0, 1]).astype(dtype1),
                  np.array([0, 1]).astype(dtype2)]
        idx = pd.MultiIndex.from_product(levels)
        assert idx.get_loc(idx[2]) == 2

    @pytest.mark.parametrize('level', [0, 1])
    @pytest.mark.parametrize('dtypes', [[int, float], [float, int]])
    def test_get_loc_implicit_cast(self, level, dtypes):
        # GH 18818, GH 15994 : as flat index, cast int to float and vice-versa
        levels = [['a', 'b'], ['c', 'd']]
        key = ['b', 'd']
        lev_dtype, key_dtype = dtypes
        levels[level] = np.array([0, 1], dtype=lev_dtype)
        key[level] = key_dtype(1)
        idx = MultiIndex.from_product(levels)
        assert idx.get_loc(tuple(key)) == 3

    def test_get_loc_cast_bool(self):
        # GH 19086 : int is casted to bool, but not vice-versa
        levels = [[False, True], np.arange(2, dtype='int64')]
        idx = MultiIndex.from_product(levels)

        assert idx.get_loc((0, 1)) == 1
        assert idx.get_loc((1, 0)) == 2

        pytest.raises(KeyError, idx.get_loc, (False, True))
        pytest.raises(KeyError, idx.get_loc, (True, False))

    def test_get_indexer(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3, 3], dtype=np.intp)
        minor_labels = np.array([0, 1, 0, 0, 1, 0, 1], dtype=np.intp)

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        idx1 = index[:5]
        idx2 = index[[1, 3, 5]]

        r1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, np.array([1, 3, -1], dtype=np.intp))

        r1 = idx2.get_indexer(idx1, method='pad')
        e1 = np.array([-1, 0, 0, 1, 1], dtype=np.intp)
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='pad')
        assert_almost_equal(r2, e1[::-1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        e1 = np.array([0, 0, 1, 1, 2], dtype=np.intp)
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='backfill')
        assert_almost_equal(r2, e1[::-1])

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

        # pass non-MultiIndex
        r1 = idx1.get_indexer(idx2.values)
        rexp1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, rexp1)

        r1 = idx1.get_indexer([1, 2, 3])
        assert (r1 == [-1, -1, -1]).all()

        # create index with duplicates
        idx1 = Index(lrange(10) + lrange(10))
        idx2 = Index(lrange(20))

        msg = "Reindexing only valid with uniquely valued Index objects"
        with tm.assert_raises_regex(InvalidIndexError, msg):
            idx1.get_indexer(idx2)

    def test_get_indexer_nearest(self):
        midx = MultiIndex.from_tuples([('a', 1), ('b', 2)])
        with pytest.raises(NotImplementedError):
            midx.get_indexer(['a'], method='nearest')
        with pytest.raises(NotImplementedError):
            midx.get_indexer(['a'], method='pad', tolerance=2)

    def test_get_level_values(self):
        result = self.index.get_level_values(0)
        expected = Index(['foo', 'foo', 'bar', 'baz', 'qux', 'qux'],
                         name='first')
        tm.assert_index_equal(result, expected)
        assert result.name == 'first'

        result = self.index.get_level_values('first')
        expected = self.index.get_level_values(0)
        tm.assert_index_equal(result, expected)

        # GH 10460
        index = MultiIndex(
            levels=[CategoricalIndex(['A', 'B']),
                    CategoricalIndex([1, 2, 3])],
            labels=[np.array([0, 0, 0, 1, 1, 1]),
                    np.array([0, 1, 2, 0, 1, 2])])

        exp = CategoricalIndex(['A', 'A', 'A', 'B', 'B', 'B'])
        tm.assert_index_equal(index.get_level_values(0), exp)
        exp = CategoricalIndex([1, 2, 3, 1, 2, 3])
        tm.assert_index_equal(index.get_level_values(1), exp)

    def test_get_level_values_int_with_na(self):
        # GH 17924
        arrays = [['a', 'b', 'b'], [1, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        result = index.get_level_values(1)
        expected = Index([1, np.nan, 2])
        tm.assert_index_equal(result, expected)

        arrays = [['a', 'b', 'b'], [np.nan, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        result = index.get_level_values(1)
        expected = Index([np.nan, np.nan, 2])
        tm.assert_index_equal(result, expected)

    def test_get_level_values_na(self):
        arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
        index = pd.MultiIndex.from_arrays(arrays)
        result = index.get_level_values(0)
        expected = pd.Index([np.nan, np.nan, np.nan])
        tm.assert_index_equal(result, expected)

        result = index.get_level_values(1)
        expected = pd.Index(['a', np.nan, 1])
        tm.assert_index_equal(result, expected)

        arrays = [['a', 'b', 'b'], pd.DatetimeIndex([0, 1, pd.NaT])]
        index = pd.MultiIndex.from_arrays(arrays)
        result = index.get_level_values(1)
        expected = pd.DatetimeIndex([0, 1, pd.NaT])
        tm.assert_index_equal(result, expected)

        arrays = [[], []]
        index = pd.MultiIndex.from_arrays(arrays)
        result = index.get_level_values(0)
        expected = pd.Index([], dtype=object)
        tm.assert_index_equal(result, expected)

    def test_get_level_values_all_na(self):
        # GH 17924 when level entirely consists of nan
        arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
        index = pd.MultiIndex.from_arrays(arrays)
        result = index.get_level_values(0)
        expected = pd.Index([np.nan, np.nan, np.nan], dtype=np.float64)
        tm.assert_index_equal(result, expected)

        result = index.get_level_values(1)
        expected = pd.Index(['a', np.nan, 1], dtype=object)
        tm.assert_index_equal(result, expected)

    def test_get_unique_index(self):
        idx = self.index[[0, 1, 0, 1, 1, 0, 0]]
        expected = self.index._shallow_copy(idx[[0, 1]])

        for dropna in [False, True]:
            result = idx._get_unique_index(dropna=dropna)
            assert result.unique
            tm.assert_index_equal(result, expected)

    def test_get_level_number_integer(self):
        self.index.names = [1, 0]
        assert self.index._get_level_number(1) == 0
        assert self.index._get_level_number(0) == 1
        pytest.raises(IndexError, self.index._get_level_number, 2)
        tm.assert_raises_regex(KeyError, 'Level fourth not found',
                               self.index._get_level_number, 'fourth')
