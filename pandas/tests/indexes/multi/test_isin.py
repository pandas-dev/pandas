# -*- coding: utf-8 -*-


import pytest

import numpy as np


from pandas import (Index, MultiIndex)
from pandas.compat import PYPY
import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestIsIn(Base):
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

    def test_isin(self):
        values = [('foo', 2), ('bar', 3), ('quux', 4)]

        idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
            4)])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(result, expected)

        # empty, return dtype bool
        idx = MultiIndex.from_arrays([[], []])
        result = idx.isin(values)
        assert len(result) == 0
        assert result.dtype == np.bool_

    @pytest.mark.skipif(PYPY, reason="tuples cmp recursively on PyPy")
    def test_isin_nan_not_pypy(self):
        idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
        tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                    np.array([False, False]))
        tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                    np.array([False, False]))

    @pytest.mark.skipif(not PYPY, reason="tuples cmp recursively on PyPy")
    def test_isin_nan_pypy(self):
        idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
        tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                    np.array([False, True]))
        tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                    np.array([False, True]))

    def test_isin_level_kwarg(self):
        idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
            4)])

        vals_0 = ['foo', 'bar', 'quux']
        vals_1 = [2, 3, 10]

        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=0))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=-2))

        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=1))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=-1))

        pytest.raises(IndexError, idx.isin, vals_0, level=5)
        pytest.raises(IndexError, idx.isin, vals_0, level=-5)

        pytest.raises(KeyError, idx.isin, vals_0, level=1.0)
        pytest.raises(KeyError, idx.isin, vals_1, level=-1.0)
        pytest.raises(KeyError, idx.isin, vals_1, level='A')

        idx.names = ['A', 'B']
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level='A'))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level='B'))

        pytest.raises(KeyError, idx.isin, vals_1, level='C')
