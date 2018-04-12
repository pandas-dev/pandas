# -*- coding: utf-8 -*-

from datetime import timedelta

import numpy as np

from pandas import (Index, MultiIndex)
from pandas.compat import lrange

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestSlice(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_slice_keep_name(self):
        x = MultiIndex.from_tuples([('a', 'b'), (1, 2), ('c', 'd')],
                                   names=['x', 'y'])
        assert x[1:].names == x.names

    def test_slice_locs(self):
        df = tm.makeTimeDataFrame()
        stacked = df.stack()
        idx = stacked.index

        slob = slice(*idx.slice_locs(df.index[5], df.index[15]))
        sliced = stacked[slob]
        expected = df[5:16].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

        slob = slice(*idx.slice_locs(df.index[5] + timedelta(seconds=30),
                                     df.index[15] - timedelta(seconds=30)))
        sliced = stacked[slob]
        expected = df[6:15].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

    def test_slice_locs_with_type_mismatch(self):
        df = tm.makeTimeDataFrame()
        stacked = df.stack()
        idx = stacked.index
        tm.assert_raises_regex(TypeError, '^Level type mismatch',
                               idx.slice_locs, (1, 3))
        tm.assert_raises_regex(TypeError, '^Level type mismatch',
                               idx.slice_locs,
                               df.index[5] + timedelta(
                                   seconds=30), (5, 2))
        df = tm.makeCustomDataframe(5, 5)
        stacked = df.stack()
        idx = stacked.index
        with tm.assert_raises_regex(TypeError, '^Level type mismatch'):
            idx.slice_locs(timedelta(seconds=30))
        # TODO: Try creating a UnicodeDecodeError in exception message
        with tm.assert_raises_regex(TypeError, '^Level type mismatch'):
            idx.slice_locs(df.index[1], (16, "a"))

    def test_slice_locs_not_sorted(self):
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        tm.assert_raises_regex(KeyError, "[Kk]ey length.*greater than "
                                         "MultiIndex lexsort depth",
                               index.slice_locs, (1, 0, 1), (2, 1, 0))

        # works
        sorted_index, _ = index.sortlevel(0)
        # should there be a test case here???
        sorted_index.slice_locs((1, 0, 1), (2, 1, 0))

    def test_slice_locs_partial(self):
        sorted_idx, _ = self.index.sortlevel(0)

        result = sorted_idx.slice_locs(('foo', 'two'), ('qux', 'one'))
        assert result == (1, 5)

        result = sorted_idx.slice_locs(None, ('qux', 'one'))
        assert result == (0, 5)

        result = sorted_idx.slice_locs(('foo', 'two'), None)
        assert result == (1, len(sorted_idx))

        result = sorted_idx.slice_locs('bar', 'baz')
        assert result == (2, 4)

    def test_slice_locs_not_contained(self):
        # some searchsorted action

        index = MultiIndex(levels=[[0, 2, 4, 6], [0, 2, 4]],
                           labels=[[0, 0, 0, 1, 1, 2, 3, 3, 3],
                                   [0, 1, 2, 1, 2, 2, 0, 1, 2]], sortorder=0)

        result = index.slice_locs((1, 0), (5, 2))
        assert result == (3, 6)

        result = index.slice_locs(1, 5)
        assert result == (3, 6)

        result = index.slice_locs((2, 2), (5, 2))
        assert result == (3, 6)

        result = index.slice_locs(2, 5)
        assert result == (3, 6)

        result = index.slice_locs((1, 0), (6, 3))
        assert result == (3, 8)

        result = index.slice_locs(-1, 10)
        assert result == (0, len(index))
