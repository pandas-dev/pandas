# -*- coding: utf-8 -*-

import pytest

import pandas as pd

from pandas import MultiIndex

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestUnique(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    @pytest.mark.parametrize('names', [None, ['first', 'second']])
    def test_unique(self, names):
        mi = pd.MultiIndex.from_arrays([[1, 2, 1, 2], [1, 1, 1, 2]],
                                       names=names)

        res = mi.unique()
        exp = pd.MultiIndex.from_arrays([[1, 2, 2], [1, 1, 2]], names=mi.names)
        tm.assert_index_equal(res, exp)

        mi = pd.MultiIndex.from_arrays([list('aaaa'), list('abab')],
                                       names=names)
        res = mi.unique()
        exp = pd.MultiIndex.from_arrays([list('aa'), list('ab')],
                                        names=mi.names)
        tm.assert_index_equal(res, exp)

        mi = pd.MultiIndex.from_arrays([list('aaaa'), list('aaaa')],
                                       names=names)
        res = mi.unique()
        exp = pd.MultiIndex.from_arrays([['a'], ['a']], names=mi.names)
        tm.assert_index_equal(res, exp)

        # GH #20568 - empty MI
        mi = pd.MultiIndex.from_arrays([[], []], names=names)
        res = mi.unique()
        tm.assert_index_equal(mi, res)

    @pytest.mark.parametrize('level', [0, 'first', 1, 'second'])
    def test_unique_level(self, level):
        # GH #17896 - with level= argument
        result = self.index.unique(level=level)
        expected = self.index.get_level_values(level).unique()
        tm.assert_index_equal(result, expected)

        # With already unique level
        mi = pd.MultiIndex.from_arrays([[1, 3, 2, 4], [1, 3, 2, 5]],
                                       names=['first', 'second'])
        result = mi.unique(level=level)
        expected = mi.get_level_values(level)
        tm.assert_index_equal(result, expected)

        # With empty MI
        mi = pd.MultiIndex.from_arrays([[], []], names=['first', 'second'])
        result = mi.unique(level=level)
        expected = mi.get_level_values(level)

    def test_unique_datetimelike(self):
        idx1 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-01',
                                 '2015-01-01', 'NaT', 'NaT'])
        idx2 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-02',
                                 '2015-01-02', 'NaT', '2015-01-01'],
                                tz='Asia/Tokyo')
        result = pd.MultiIndex.from_arrays([idx1, idx2]).unique()

        eidx1 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', 'NaT', 'NaT'])
        eidx2 = pd.DatetimeIndex(['2015-01-01', '2015-01-02',
                                  'NaT', '2015-01-01'],
                                 tz='Asia/Tokyo')
        exp = pd.MultiIndex.from_arrays([eidx1, eidx2])
        tm.assert_index_equal(result, exp)
