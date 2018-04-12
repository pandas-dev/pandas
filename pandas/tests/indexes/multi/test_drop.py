# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)
from pandas.compat import lrange
from pandas.errors import PerformanceWarning

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestDrop(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_drop(self):
        dropped = self.index.drop([('foo', 'two'), ('qux', 'one')])

        index = MultiIndex.from_tuples([('foo', 'two'), ('qux', 'one')])
        dropped2 = self.index.drop(index)

        expected = self.index[[0, 2, 3, 5]]
        tm.assert_index_equal(dropped, expected)
        tm.assert_index_equal(dropped2, expected)

        dropped = self.index.drop(['bar'])
        expected = self.index[[0, 1, 3, 4, 5]]
        tm.assert_index_equal(dropped, expected)

        dropped = self.index.drop('foo')
        expected = self.index[[2, 3, 4, 5]]
        tm.assert_index_equal(dropped, expected)

        index = MultiIndex.from_tuples([('bar', 'two')])
        pytest.raises(KeyError, self.index.drop, [('bar', 'two')])
        pytest.raises(KeyError, self.index.drop, index)
        pytest.raises(KeyError, self.index.drop, ['foo', 'two'])

        # partially correct argument
        mixed_index = MultiIndex.from_tuples([('qux', 'one'), ('bar', 'two')])
        pytest.raises(KeyError, self.index.drop, mixed_index)

        # error='ignore'
        dropped = self.index.drop(index, errors='ignore')
        expected = self.index[[0, 1, 2, 3, 4, 5]]
        tm.assert_index_equal(dropped, expected)

        dropped = self.index.drop(mixed_index, errors='ignore')
        expected = self.index[[0, 1, 2, 3, 5]]
        tm.assert_index_equal(dropped, expected)

        dropped = self.index.drop(['foo', 'two'], errors='ignore')
        expected = self.index[[2, 3, 4, 5]]
        tm.assert_index_equal(dropped, expected)

        # mixed partial / full drop
        dropped = self.index.drop(['foo', ('qux', 'one')])
        expected = self.index[[2, 3, 5]]
        tm.assert_index_equal(dropped, expected)

        # mixed partial / full drop / error='ignore'
        mixed_index = ['foo', ('qux', 'one'), 'two']
        pytest.raises(KeyError, self.index.drop, mixed_index)
        dropped = self.index.drop(mixed_index, errors='ignore')
        expected = self.index[[2, 3, 5]]
        tm.assert_index_equal(dropped, expected)

    def test_droplevel_with_names(self):
        index = self.index[self.index.get_loc('foo')]
        dropped = index.droplevel(0)
        assert dropped.name == 'second'

        index = MultiIndex(
            levels=[Index(lrange(4)), Index(lrange(4)), Index(lrange(4))],
            labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])],
            names=['one', 'two', 'three'])
        dropped = index.droplevel(0)
        assert dropped.names == ('two', 'three')

        dropped = index.droplevel('two')
        expected = index.droplevel(1)
        assert dropped.equals(expected)

    def test_droplevel_multiple(self):
        index = MultiIndex(
            levels=[Index(lrange(4)), Index(lrange(4)), Index(lrange(4))],
            labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])],
            names=['one', 'two', 'three'])

        dropped = index[:2].droplevel(['three', 'one'])
        expected = index[:2].droplevel(2).droplevel(0)
        assert dropped.equals(expected)

    def test_drop_not_lexsorted(self):
        # GH 12078

        # define the lexsorted version of the multi-index
        tuples = [('a', ''), ('b1', 'c1'), ('b2', 'c2')]
        lexsorted_mi = MultiIndex.from_tuples(tuples, names=['b', 'c'])
        assert lexsorted_mi.is_lexsorted()

        # and the not-lexsorted version
        df = pd.DataFrame(columns=['a', 'b', 'c', 'd'],
                          data=[[1, 'b1', 'c1', 3], [1, 'b2', 'c2', 4]])
        df = df.pivot_table(index='a', columns=['b', 'c'], values='d')
        df = df.reset_index()
        not_lexsorted_mi = df.columns
        assert not not_lexsorted_mi.is_lexsorted()

        # compare the results
        tm.assert_index_equal(lexsorted_mi, not_lexsorted_mi)
        with tm.assert_produces_warning(PerformanceWarning):
            tm.assert_index_equal(lexsorted_mi.drop('a'),
                                  not_lexsorted_mi.drop('a'))
