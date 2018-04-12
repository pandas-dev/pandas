# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)

from pandas.tests.indexes.common import Base


class TestMonotonic(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_is_monotonic_increasing(self):
        i = MultiIndex.from_product([np.arange(10),
                                     np.arange(10)], names=['one', 'two'])
        assert i.is_monotonic
        assert i._is_strictly_monotonic_increasing
        assert Index(i.values).is_monotonic
        assert i._is_strictly_monotonic_increasing

        i = MultiIndex.from_product([np.arange(10, 0, -1),
                                     np.arange(10)], names=['one', 'two'])
        assert not i.is_monotonic
        assert not i._is_strictly_monotonic_increasing
        assert not Index(i.values).is_monotonic
        assert not Index(i.values)._is_strictly_monotonic_increasing

        i = MultiIndex.from_product([np.arange(10),
                                     np.arange(10, 0, -1)],
                                    names=['one', 'two'])
        assert not i.is_monotonic
        assert not i._is_strictly_monotonic_increasing
        assert not Index(i.values).is_monotonic
        assert not Index(i.values)._is_strictly_monotonic_increasing

        i = MultiIndex.from_product([[1.0, np.nan, 2.0], ['a', 'b', 'c']])
        assert not i.is_monotonic
        assert not i._is_strictly_monotonic_increasing
        assert not Index(i.values).is_monotonic
        assert not Index(i.values)._is_strictly_monotonic_increasing

        # string ordering
        i = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                               ['one', 'two', 'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
        assert not i.is_monotonic
        assert not Index(i.values).is_monotonic
        assert not i._is_strictly_monotonic_increasing
        assert not Index(i.values)._is_strictly_monotonic_increasing

        i = MultiIndex(levels=[['bar', 'baz', 'foo', 'qux'],
                               ['mom', 'next', 'zenith']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
        assert i.is_monotonic
        assert Index(i.values).is_monotonic
        assert i._is_strictly_monotonic_increasing
        assert Index(i.values)._is_strictly_monotonic_increasing

        # mixed levels, hits the TypeError
        i = MultiIndex(
            levels=[[1, 2, 3, 4], ['gb00b03mlx29', 'lu0197800237',
                                   'nl0000289783',
                                   'nl0000289965', 'nl0000301109']],
            labels=[[0, 1, 1, 2, 2, 2, 3], [4, 2, 0, 0, 1, 3, -1]],
            names=['household_id', 'asset_id'])

        assert not i.is_monotonic
        assert not i._is_strictly_monotonic_increasing

        # empty
        i = MultiIndex.from_arrays([[], []])
        assert i.is_monotonic
        assert Index(i.values).is_monotonic
        assert i._is_strictly_monotonic_increasing
        assert Index(i.values)._is_strictly_monotonic_increasing

    def test_is_monotonic_decreasing(self):
        i = MultiIndex.from_product([np.arange(9, -1, -1),
                                     np.arange(9, -1, -1)],
                                    names=['one', 'two'])
        assert i.is_monotonic_decreasing
        assert i._is_strictly_monotonic_decreasing
        assert Index(i.values).is_monotonic_decreasing
        assert i._is_strictly_monotonic_decreasing

        i = MultiIndex.from_product([np.arange(10),
                                     np.arange(10, 0, -1)],
                                    names=['one', 'two'])
        assert not i.is_monotonic_decreasing
        assert not i._is_strictly_monotonic_decreasing
        assert not Index(i.values).is_monotonic_decreasing
        assert not Index(i.values)._is_strictly_monotonic_decreasing

        i = MultiIndex.from_product([np.arange(10, 0, -1),
                                     np.arange(10)], names=['one', 'two'])
        assert not i.is_monotonic_decreasing
        assert not i._is_strictly_monotonic_decreasing
        assert not Index(i.values).is_monotonic_decreasing
        assert not Index(i.values)._is_strictly_monotonic_decreasing

        i = MultiIndex.from_product([[2.0, np.nan, 1.0], ['c', 'b', 'a']])
        assert not i.is_monotonic_decreasing
        assert not i._is_strictly_monotonic_decreasing
        assert not Index(i.values).is_monotonic_decreasing
        assert not Index(i.values)._is_strictly_monotonic_decreasing

        # string ordering
        i = MultiIndex(levels=[['qux', 'foo', 'baz', 'bar'],
                               ['three', 'two', 'one']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
        assert not i.is_monotonic_decreasing
        assert not Index(i.values).is_monotonic_decreasing
        assert not i._is_strictly_monotonic_decreasing
        assert not Index(i.values)._is_strictly_monotonic_decreasing

        i = MultiIndex(levels=[['qux', 'foo', 'baz', 'bar'],
                               ['zenith', 'next', 'mom']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
        assert i.is_monotonic_decreasing
        assert Index(i.values).is_monotonic_decreasing
        assert i._is_strictly_monotonic_decreasing
        assert Index(i.values)._is_strictly_monotonic_decreasing

        # mixed levels, hits the TypeError
        i = MultiIndex(
            levels=[[4, 3, 2, 1], ['nl0000301109', 'nl0000289965',
                                   'nl0000289783', 'lu0197800237',
                                   'gb00b03mlx29']],
            labels=[[0, 1, 1, 2, 2, 2, 3], [4, 2, 0, 0, 1, 3, -1]],
            names=['household_id', 'asset_id'])

        assert not i.is_monotonic_decreasing
        assert not i._is_strictly_monotonic_decreasing

        # empty
        i = MultiIndex.from_arrays([[], []])
        assert i.is_monotonic_decreasing
        assert Index(i.values).is_monotonic_decreasing
        assert i._is_strictly_monotonic_decreasing
        assert Index(i.values)._is_strictly_monotonic_decreasing

    def test_is_strictly_monotonic_increasing(self):
        idx = pd.MultiIndex(levels=[['bar', 'baz'], ['mom', 'next']],
                            labels=[[0, 0, 1, 1], [0, 0, 0, 1]])
        assert idx.is_monotonic_increasing
        assert not idx._is_strictly_monotonic_increasing

    def test_is_strictly_monotonic_decreasing(self):
        idx = pd.MultiIndex(levels=[['baz', 'bar'], ['next', 'mom']],
                            labels=[[0, 0, 1, 1], [0, 0, 0, 1]])
        assert idx.is_monotonic_decreasing
        assert not idx._is_strictly_monotonic_decreasing
