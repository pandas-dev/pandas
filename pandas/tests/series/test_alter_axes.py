# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime

import numpy as np
import pytest

from pandas.compat import lrange, range, zip

from pandas import DataFrame, Index, MultiIndex, RangeIndex, Series
import pandas.util.testing as tm


class TestSeriesAlterAxes(object):

    def test_setindex(self, string_series):
        # wrong type
        msg = (r"Index\(\.\.\.\) must be called with a collection of some"
               r" kind, None was passed")
        with pytest.raises(TypeError, match=msg):
            string_series.index = None

        # wrong length
        msg = ("Length mismatch: Expected axis has 30 elements, new"
               " values have 29 elements")
        with pytest.raises(ValueError, match=msg):
            string_series.index = np.arange(len(string_series) - 1)

        # works
        string_series.index = np.arange(len(string_series))
        assert isinstance(string_series.index, Index)

    # Renaming

    def test_rename(self, datetime_series):
        ts = datetime_series
        renamer = lambda x: x.strftime('%Y%m%d')
        renamed = ts.rename(renamer)
        assert renamed.index[0] == renamer(ts.index[0])

        # dict
        rename_dict = dict(zip(ts.index, renamed.index))
        renamed2 = ts.rename(rename_dict)
        tm.assert_series_equal(renamed, renamed2)

        # partial dict
        s = Series(np.arange(4), index=['a', 'b', 'c', 'd'], dtype='int64')
        renamed = s.rename({'b': 'foo', 'd': 'bar'})
        tm.assert_index_equal(renamed.index, Index(['a', 'foo', 'c', 'bar']))

        # index with name
        renamer = Series(np.arange(4),
                         index=Index(['a', 'b', 'c', 'd'], name='name'),
                         dtype='int64')
        renamed = renamer.rename({})
        assert renamed.index.name == renamer.index.name

    def test_rename_by_series(self):
        s = Series(range(5), name='foo')
        renamer = Series({1: 10, 2: 20})
        result = s.rename(renamer)
        expected = Series(range(5), index=[0, 10, 20, 3, 4], name='foo')
        tm.assert_series_equal(result, expected)

    def test_rename_set_name(self):
        s = Series(range(4), index=list('abcd'))
        for name in ['foo', 123, 123., datetime(2001, 11, 11), ('foo',)]:
            result = s.rename(name)
            assert result.name == name
            tm.assert_numpy_array_equal(result.index.values, s.index.values)
            assert s.name is None

    def test_rename_set_name_inplace(self):
        s = Series(range(3), index=list('abc'))
        for name in ['foo', 123, 123., datetime(2001, 11, 11), ('foo',)]:
            s.rename(name, inplace=True)
            assert s.name == name

            exp = np.array(['a', 'b', 'c'], dtype=np.object_)
            tm.assert_numpy_array_equal(s.index.values, exp)

    def test_rename_axis_supported(self):
        # Supporting axis for compatibility, detailed in GH-18589
        s = Series(range(5))
        s.rename({}, axis=0)
        s.rename({}, axis='index')
        with pytest.raises(ValueError, match='No axis named 5'):
            s.rename({}, axis=5)

    def test_set_name_attribute(self):
        s = Series([1, 2, 3])
        s2 = Series([1, 2, 3], name='bar')
        for name in [7, 7., 'name', datetime(2001, 1, 1), (1,), u"\u05D0"]:
            s.name = name
            assert s.name == name
            s2.name = name
            assert s2.name == name

    def test_set_name(self):
        s = Series([1, 2, 3])
        s2 = s._set_name('foo')
        assert s2.name == 'foo'
        assert s.name is None
        assert s is not s2

    def test_rename_inplace(self, datetime_series):
        renamer = lambda x: x.strftime('%Y%m%d')
        expected = renamer(datetime_series.index[0])

        datetime_series.rename(renamer, inplace=True)
        assert datetime_series.index[0] == expected

    def test_set_index_makes_timeseries(self):
        idx = tm.makeDateIndex(10)

        s = Series(lrange(10))
        s.index = idx
        assert s.index.is_all_dates

    def test_reset_index(self):
        df = tm.makeDataFrame()[:5]
        ser = df.stack()
        ser.index.names = ['hash', 'category']

        ser.name = 'value'
        df = ser.reset_index()
        assert 'value' in df

        df = ser.reset_index(name='value2')
        assert 'value2' in df

        # check inplace
        s = ser.reset_index(drop=True)
        s2 = ser
        s2.reset_index(drop=True, inplace=True)
        tm.assert_series_equal(s, s2)

        # level
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           codes=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2],
                                  [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        rs = s.reset_index(level=1)
        assert len(rs.columns) == 2

        rs = s.reset_index(level=[0, 2], drop=True)
        tm.assert_index_equal(rs.index, Index(index.get_level_values(1)))
        assert isinstance(rs, Series)

    def test_reset_index_name(self):
        s = Series([1, 2, 3], index=Index(range(3), name='x'))
        assert s.reset_index().index.name is None
        assert s.reset_index(drop=True).index.name is None

    def test_reset_index_level(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]],
                       columns=['A', 'B', 'C'])

        for levels in ['A', 'B'], [0, 1]:
            # With MultiIndex
            s = df.set_index(['A', 'B'])['C']

            result = s.reset_index(level=levels[0])
            tm.assert_frame_equal(result, df.set_index('B'))

            result = s.reset_index(level=levels[:1])
            tm.assert_frame_equal(result, df.set_index('B'))

            result = s.reset_index(level=levels)
            tm.assert_frame_equal(result, df)

            result = df.set_index(['A', 'B']).reset_index(level=levels,
                                                          drop=True)
            tm.assert_frame_equal(result, df[['C']])

            with pytest.raises(KeyError, match='Level E '):
                s.reset_index(level=['A', 'E'])

            # With single-level Index
            s = df.set_index('A')['B']

            result = s.reset_index(level=levels[0])
            tm.assert_frame_equal(result, df[['A', 'B']])

            result = s.reset_index(level=levels[:1])
            tm.assert_frame_equal(result, df[['A', 'B']])

            result = s.reset_index(level=levels[0], drop=True)
            tm.assert_series_equal(result, df['B'])

            with pytest.raises(IndexError, match='Too many levels'):
                s.reset_index(level=[0, 1, 2])

        # Check that .reset_index([],drop=True) doesn't fail
        result = Series(range(4)).reset_index([], drop=True)
        expected = Series(range(4))
        tm.assert_series_equal(result, expected)

    def test_reset_index_range(self):
        # GH 12071
        s = Series(range(2), name='A', dtype='int64')
        series_result = s.reset_index()
        assert isinstance(series_result.index, RangeIndex)
        series_expected = DataFrame([[0, 0], [1, 1]],
                                    columns=['index', 'A'],
                                    index=RangeIndex(stop=2))
        tm.assert_frame_equal(series_result, series_expected)

    def test_reorder_levels(self):
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           codes=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2],
                                  [0, 1, 0, 1, 0, 1]],
                           names=['L0', 'L1', 'L2'])
        s = Series(np.arange(6), index=index)

        # no change, position
        result = s.reorder_levels([0, 1, 2])
        tm.assert_series_equal(s, result)

        # no change, labels
        result = s.reorder_levels(['L0', 'L1', 'L2'])
        tm.assert_series_equal(s, result)

        # rotate, position
        result = s.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(levels=[['one', 'two', 'three'], [0, 1], ['bar']],
                           codes=[[0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1],
                                  [0, 0, 0, 0, 0, 0]],
                           names=['L1', 'L2', 'L0'])
        expected = Series(np.arange(6), index=e_idx)
        tm.assert_series_equal(result, expected)

    def test_rename_axis_mapper(self):
        # GH 19978
        mi = MultiIndex.from_product([['a', 'b', 'c'], [1, 2]],
                                     names=['ll', 'nn'])
        s = Series([i for i in range(len(mi))], index=mi)

        result = s.rename_axis(index={'ll': 'foo'})
        assert result.index.names == ['foo', 'nn']

        result = s.rename_axis(index=str.upper, axis=0)
        assert result.index.names == ['LL', 'NN']

        result = s.rename_axis(index=['foo', 'goo'])
        assert result.index.names == ['foo', 'goo']

        with pytest.raises(TypeError, match='unexpected'):
            s.rename_axis(columns='wrong')

    def test_rename_axis_inplace(self, datetime_series):
        # GH 15704
        expected = datetime_series.rename_axis('foo')
        result = datetime_series
        no_return = result.rename_axis('foo', inplace=True)

        assert no_return is None
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('kwargs', [{'mapper': None}, {'index': None}, {}])
    def test_rename_axis_none(self, kwargs):
        # GH 25034
        index = Index(list('abc'), name='foo')
        df = Series([1, 2, 3], index=index)

        result = df.rename_axis(**kwargs)
        expected_index = index.rename(None) if kwargs else index
        expected = Series([1, 2, 3], index=expected_index)
        tm.assert_series_equal(result, expected)

    def test_set_axis_inplace_axes(self, axis_series):
        # GH14636
        ser = Series(np.arange(4), index=[1, 3, 5, 7], dtype='int64')

        expected = ser.copy()
        expected.index = list('abcd')

        # inplace=True
        # The FutureWarning comes from the fact that we would like to have
        # inplace default to False some day
        for inplace, warn in [(None, FutureWarning), (True, None)]:
            result = ser.copy()
            kwargs = {'inplace': inplace}
            with tm.assert_produces_warning(warn):
                result.set_axis(list('abcd'), axis=axis_series, **kwargs)
            tm.assert_series_equal(result, expected)

    def test_set_axis_inplace(self):
        # GH14636

        s = Series(np.arange(4), index=[1, 3, 5, 7], dtype='int64')

        expected = s.copy()
        expected.index = list('abcd')

        # inplace=False
        result = s.set_axis(list('abcd'), axis=0, inplace=False)
        tm.assert_series_equal(expected, result)

        # omitting the "axis" parameter
        with tm.assert_produces_warning(None):
            result = s.set_axis(list('abcd'), inplace=False)
        tm.assert_series_equal(result, expected)

        # wrong values for the "axis" parameter
        for axis in [2, 'foo']:
            with pytest.raises(ValueError, match='No axis named'):
                s.set_axis(list('abcd'), axis=axis, inplace=False)

    def test_set_axis_prior_to_deprecation_signature(self):
        s = Series(np.arange(4), index=[1, 3, 5, 7], dtype='int64')

        expected = s.copy()
        expected.index = list('abcd')

        for axis in [0, 'index']:
            with tm.assert_produces_warning(FutureWarning):
                result = s.set_axis(0, list('abcd'), inplace=False)
            tm.assert_series_equal(result, expected)

    def test_reset_index_drop_errors(self):
        #  GH 20925

        # KeyError raised for series index when passed level name is missing
        s = Series(range(4))
        with pytest.raises(KeyError, match='must be same as name'):
            s.reset_index('wrong', drop=True)
        with pytest.raises(KeyError, match='must be same as name'):
            s.reset_index('wrong')

        # KeyError raised for series when level to be dropped is missing
        s = Series(range(4), index=MultiIndex.from_product([[1, 2]] * 2))
        with pytest.raises(KeyError, match='not found'):
            s.reset_index('wrong', drop=True)

    def test_droplevel(self):
        # GH20342
        ser = Series([1, 2, 3, 4])
        ser.index = MultiIndex.from_arrays([(1, 2, 3, 4), (5, 6, 7, 8)],
                                           names=['a', 'b'])
        expected = ser.reset_index('b', drop=True)
        result = ser.droplevel('b', axis='index')
        tm.assert_series_equal(result, expected)
        # test that droplevel raises ValueError on axis != 0
        with pytest.raises(ValueError):
            ser.droplevel(1, axis='columns')
