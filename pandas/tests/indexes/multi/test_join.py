# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestJoin(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def setup_method(self, method):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])
        self.index_names = ['first', 'second']
        self.indices = dict(index=MultiIndex(levels=[major_axis, minor_axis],
                                             labels=[major_labels, minor_labels
                                                     ], names=self.index_names,
                                             verify_integrity=False))
        self.setup_indices()

    def create_index(self):
        return self.index

    @pytest.mark.parametrize('other',
                             [Index(['three', 'one', 'two']),
                              Index(['one']),
                              Index(['one', 'three'])])
    def test_join_level(self, other, join_type):
        join_index, lidx, ridx = other.join(self.index, how=join_type,
                                            level='second',
                                            return_indexers=True)

        exp_level = other.join(self.index.levels[1], how=join_type)
        assert join_index.levels[0].equals(self.index.levels[0])
        assert join_index.levels[1].equals(exp_level)

        # pare down levels
        mask = np.array(
            [x[1] in exp_level for x in self.index], dtype=bool)
        exp_values = self.index.values[mask]
        tm.assert_numpy_array_equal(join_index.values, exp_values)

        if join_type in ('outer', 'inner'):
            join_index2, ridx2, lidx2 = \
                self.index.join(other, how=join_type, level='second',
                                return_indexers=True)

            assert join_index.equals(join_index2)
            tm.assert_numpy_array_equal(lidx, lidx2)
            tm.assert_numpy_array_equal(ridx, ridx2)
            tm.assert_numpy_array_equal(join_index2.values, exp_values)

    def test_join_level_corner_case(self):
        # some corner cases
        idx = Index(['three', 'one', 'two'])
        result = idx.join(self.index, level='second')
        assert isinstance(result, MultiIndex)

        tm.assert_raises_regex(TypeError, "Join.*MultiIndex.*ambiguous",
                               self.index.join, self.index, level=1)

    def test_join_self(self, join_type):
        res = self.index
        joined = res.join(res, how=join_type)
        assert res is joined

    def test_join_multi(self):
        # GH 10665
        midx = pd.MultiIndex.from_product(
            [np.arange(4), np.arange(4)], names=['a', 'b'])
        idx = pd.Index([1, 2, 5], name='b')

        # inner
        jidx, lidx, ridx = midx.join(idx, how='inner', return_indexers=True)
        exp_idx = pd.MultiIndex.from_product(
            [np.arange(4), [1, 2]], names=['a', 'b'])
        exp_lidx = np.array([1, 2, 5, 6, 9, 10, 13, 14], dtype=np.intp)
        exp_ridx = np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=np.intp)
        tm.assert_index_equal(jidx, exp_idx)
        tm.assert_numpy_array_equal(lidx, exp_lidx)
        tm.assert_numpy_array_equal(ridx, exp_ridx)
        # flip
        jidx, ridx, lidx = idx.join(midx, how='inner', return_indexers=True)
        tm.assert_index_equal(jidx, exp_idx)
        tm.assert_numpy_array_equal(lidx, exp_lidx)
        tm.assert_numpy_array_equal(ridx, exp_ridx)

        # keep MultiIndex
        jidx, lidx, ridx = midx.join(idx, how='left', return_indexers=True)
        exp_ridx = np.array([-1, 0, 1, -1, -1, 0, 1, -1, -1, 0, 1, -1, -1, 0,
                             1, -1], dtype=np.intp)
        tm.assert_index_equal(jidx, midx)
        assert lidx is None
        tm.assert_numpy_array_equal(ridx, exp_ridx)
        # flip
        jidx, ridx, lidx = idx.join(midx, how='right', return_indexers=True)
        tm.assert_index_equal(jidx, midx)
        assert lidx is None
        tm.assert_numpy_array_equal(ridx, exp_ridx)
