# -*- coding: utf-8 -*-

from __future__ import print_function

from datetime import datetime, timedelta
import itertools

import numpy as np
import pytest

from pandas.compat import StringIO

import pandas as pd
from pandas import (
    Categorical, DataFrame, Series, Timestamp, compat, date_range,
    option_context)
from pandas.core.arrays import IntervalArray, integer_array
from pandas.core.internals.blocks import IntBlock
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_almost_equal, assert_frame_equal, assert_series_equal)

# Segregated collection of methods that require the BlockManager internal data
# structure


class TestDataFrameBlockInternals():
    def test_setitem_invalidates_datetime_index_freq(self):
        # GH#24096 altering a datetime64tz column inplace invalidates the
        #  `freq` attribute on the underlying DatetimeIndex

        dti = date_range('20130101', periods=3, tz='US/Eastern')
        ts = dti[1]

        df = DataFrame({'B': dti})
        assert df['B']._values.freq == 'D'

        df.iloc[1, 0] = pd.NaT
        assert df['B']._values.freq is None

        # check that the DatetimeIndex was not altered in place
        assert dti.freq == 'D'
        assert dti[1] == ts

    def test_cast_internals(self, float_frame):
        casted = DataFrame(float_frame._data, dtype=int)
        expected = DataFrame(float_frame._series, dtype=int)
        assert_frame_equal(casted, expected)

        casted = DataFrame(float_frame._data, dtype=np.int32)
        expected = DataFrame(float_frame._series, dtype=np.int32)
        assert_frame_equal(casted, expected)

    def test_consolidate(self, float_frame):
        float_frame['E'] = 7.
        consolidated = float_frame._consolidate()
        assert len(consolidated._data.blocks) == 1

        # Ensure copy, do I want this?
        recons = consolidated._consolidate()
        assert recons is not consolidated
        tm.assert_frame_equal(recons, consolidated)

        float_frame['F'] = 8.
        assert len(float_frame._data.blocks) == 3

        float_frame._consolidate(inplace=True)
        assert len(float_frame._data.blocks) == 1

    def test_consolidate_inplace(self, float_frame):
        frame = float_frame.copy()  # noqa

        # triggers in-place consolidation
        for letter in range(ord('A'), ord('Z')):
            float_frame[chr(letter)] = chr(letter)

    def test_values_consolidate(self, float_frame):
        float_frame['E'] = 7.
        assert not float_frame._data.is_consolidated()
        _ = float_frame.values  # noqa
        assert float_frame._data.is_consolidated()

    def test_modify_values(self, float_frame):
        float_frame.values[5] = 5
        assert (float_frame.values[5] == 5).all()

        # unconsolidated
        float_frame['E'] = 7.
        float_frame.values[6] = 6
        assert (float_frame.values[6] == 6).all()

    def test_boolean_set_uncons(self, float_frame):
        float_frame['E'] = 7.

        expected = float_frame.values.copy()
        expected[expected > 1] = 2

        float_frame[float_frame > 1] = 2
        assert_almost_equal(expected, float_frame.values)

    def test_values_numeric_cols(self, float_frame):
        float_frame['foo'] = 'bar'

        values = float_frame[['A', 'B', 'C', 'D']].values
        assert values.dtype == np.float64

    def test_values_lcd(self, mixed_float_frame, mixed_int_frame):

        # mixed lcd
        values = mixed_float_frame[['A', 'B', 'C', 'D']].values
        assert values.dtype == np.float64

        values = mixed_float_frame[['A', 'B', 'C']].values
        assert values.dtype == np.float32

        values = mixed_float_frame[['C']].values
        assert values.dtype == np.float16

        # GH 10364
        # B uint64 forces float because there are other signed int types
        values = mixed_int_frame[['A', 'B', 'C', 'D']].values
        assert values.dtype == np.float64

        values = mixed_int_frame[['A', 'D']].values
        assert values.dtype == np.int64

        # B uint64 forces float because there are other signed int types
        values = mixed_int_frame[['A', 'B', 'C']].values
        assert values.dtype == np.float64

        # as B and C are both unsigned, no forcing to float is needed
        values = mixed_int_frame[['B', 'C']].values
        assert values.dtype == np.uint64

        values = mixed_int_frame[['A', 'C']].values
        assert values.dtype == np.int32

        values = mixed_int_frame[['C', 'D']].values
        assert values.dtype == np.int64

        values = mixed_int_frame[['A']].values
        assert values.dtype == np.int32

        values = mixed_int_frame[['C']].values
        assert values.dtype == np.uint8

    def test_constructor_with_convert(self):
        # this is actually mostly a test of lib.maybe_convert_objects
        # #2845
        df = DataFrame({'A': [2 ** 63 - 1]})
        result = df['A']
        expected = Series(np.asarray([2 ** 63 - 1], np.int64), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [2 ** 63]})
        result = df['A']
        expected = Series(np.asarray([2 ** 63], np.uint64), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [datetime(2005, 1, 1), True]})
        result = df['A']
        expected = Series(np.asarray([datetime(2005, 1, 1), True], np.object_),
                          name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [None, 1]})
        result = df['A']
        expected = Series(np.asarray([np.nan, 1], np.float_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [1.0, 2]})
        result = df['A']
        expected = Series(np.asarray([1.0, 2], np.float_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [1.0 + 2.0j, 3]})
        result = df['A']
        expected = Series(np.asarray([1.0 + 2.0j, 3], np.complex_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [1.0 + 2.0j, 3.0]})
        result = df['A']
        expected = Series(np.asarray([1.0 + 2.0j, 3.0], np.complex_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [1.0 + 2.0j, True]})
        result = df['A']
        expected = Series(np.asarray([1.0 + 2.0j, True], np.object_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [1.0, None]})
        result = df['A']
        expected = Series(np.asarray([1.0, np.nan], np.float_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [1.0 + 2.0j, None]})
        result = df['A']
        expected = Series(np.asarray(
            [1.0 + 2.0j, np.nan], np.complex_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [2.0, 1, True, None]})
        result = df['A']
        expected = Series(np.asarray(
            [2.0, 1, True, None], np.object_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A': [2.0, 1, datetime(2006, 1, 1), None]})
        result = df['A']
        expected = Series(np.asarray([2.0, 1, datetime(2006, 1, 1),
                                      None], np.object_), name='A')
        assert_series_equal(result, expected)

    def test_construction_with_mixed(self, float_string_frame):
        # test construction edge cases with mixed types

        # f7u12, this does not work without extensive workaround
        data = [[datetime(2001, 1, 5), np.nan, datetime(2001, 1, 2)],
                [datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 1)]]
        df = DataFrame(data)

        # check dtypes
        result = df.get_dtype_counts().sort_values()
        expected = Series({'datetime64[ns]': 3})

        # mixed-type frames
        float_string_frame['datetime'] = datetime.now()
        float_string_frame['timedelta'] = timedelta(days=1, seconds=1)
        assert float_string_frame['datetime'].dtype == 'M8[ns]'
        assert float_string_frame['timedelta'].dtype == 'm8[ns]'
        result = float_string_frame.get_dtype_counts().sort_values()
        expected = Series({'float64': 4,
                           'object': 1,
                           'datetime64[ns]': 1,
                           'timedelta64[ns]': 1}).sort_values()
        assert_series_equal(result, expected)

    def test_construction_with_conversions(self):

        # convert from a numpy array of non-ns timedelta64
        arr = np.array([1, 2, 3], dtype='timedelta64[s]')
        df = DataFrame(index=range(3))
        df['A'] = arr
        expected = DataFrame({'A': pd.timedelta_range('00:00:01', periods=3,
                                                      freq='s')},
                             index=range(3))
        assert_frame_equal(df, expected)

        expected = DataFrame({
            'dt1': Timestamp('20130101'),
            'dt2': date_range('20130101', periods=3),
            # 'dt3' : date_range('20130101 00:00:01',periods=3,freq='s'),
        }, index=range(3))

        df = DataFrame(index=range(3))
        df['dt1'] = np.datetime64('2013-01-01')
        df['dt2'] = np.array(['2013-01-01', '2013-01-02', '2013-01-03'],
                             dtype='datetime64[D]')

        # df['dt3'] = np.array(['2013-01-01 00:00:01','2013-01-01
        # 00:00:02','2013-01-01 00:00:03'],dtype='datetime64[s]')

        assert_frame_equal(df, expected)

    def test_constructor_compound_dtypes(self):
        # GH 5191
        # compound dtypes should raise not-implementederror

        def f(dtype):
            data = list(itertools.repeat((datetime(2001, 1, 1),
                                          "aa", 20), 9))
            return DataFrame(data=data,
                             columns=["A", "B", "C"],
                             dtype=dtype)

        pytest.raises(NotImplementedError, f,
                      [("A", "datetime64[h]"),
                       ("B", "str"),
                       ("C", "int32")])

        # these work (though results may be unexpected)
        f('int64')
        f('float64')

        # 10822
        # invalid error message on dt inference
        if not compat.is_platform_windows():
            f('M8[ns]')

    def test_equals_different_blocks(self):
        # GH 9330
        df0 = pd.DataFrame({"A": ["x", "y"], "B": [1, 2],
                            "C": ["w", "z"]})
        df1 = df0.reset_index()[["A", "B", "C"]]
        # this assert verifies that the above operations have
        # induced a block rearrangement
        assert (df0._data.blocks[0].dtype != df1._data.blocks[0].dtype)

        # do the real tests
        assert_frame_equal(df0, df1)
        assert df0.equals(df1)
        assert df1.equals(df0)

    def test_copy_blocks(self, float_frame):
        # API/ENH 9607
        df = DataFrame(float_frame, copy=True)
        column = df.columns[0]

        # use the default copy=True, change a column

        # deprecated 0.21.0
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            blocks = df.as_blocks()
        for dtype, _df in blocks.items():
            if column in _df:
                _df.loc[:, column] = _df[column] + 1

        # make sure we did not change the original DataFrame
        assert not _df[column].equals(df[column])

    def test_no_copy_blocks(self, float_frame):
        # API/ENH 9607
        df = DataFrame(float_frame, copy=True)
        column = df.columns[0]

        # use the copy=False, change a column

        # deprecated 0.21.0
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            blocks = df.as_blocks(copy=False)
        for dtype, _df in blocks.items():
            if column in _df:
                _df.loc[:, column] = _df[column] + 1

        # make sure we did change the original DataFrame
        assert _df[column].equals(df[column])

    def test_copy(self, float_frame, float_string_frame):
        cop = float_frame.copy()
        cop['E'] = cop['A']
        assert 'E' not in float_frame

        # copy objects
        copy = float_string_frame.copy()
        assert copy._data is not float_string_frame._data

    def test_pickle(self, float_string_frame, empty_frame, timezone_frame):
        unpickled = tm.round_trip_pickle(float_string_frame)
        assert_frame_equal(float_string_frame, unpickled)

        # buglet
        float_string_frame._data.ndim

        # empty
        unpickled = tm.round_trip_pickle(empty_frame)
        repr(unpickled)

        # tz frame
        unpickled = tm.round_trip_pickle(timezone_frame)
        assert_frame_equal(timezone_frame, unpickled)

    def test_consolidate_datetime64(self):
        # numpy vstack bug

        data = """\
starting,ending,measure
2012-06-21 00:00,2012-06-23 07:00,77
2012-06-23 07:00,2012-06-23 16:30,65
2012-06-23 16:30,2012-06-25 08:00,77
2012-06-25 08:00,2012-06-26 12:00,0
2012-06-26 12:00,2012-06-27 08:00,77
"""
        df = pd.read_csv(StringIO(data), parse_dates=[0, 1])

        ser_starting = df.starting
        ser_starting.index = ser_starting.values
        ser_starting = ser_starting.tz_localize('US/Eastern')
        ser_starting = ser_starting.tz_convert('UTC')
        ser_starting.index.name = 'starting'

        ser_ending = df.ending
        ser_ending.index = ser_ending.values
        ser_ending = ser_ending.tz_localize('US/Eastern')
        ser_ending = ser_ending.tz_convert('UTC')
        ser_ending.index.name = 'ending'

        df.starting = ser_starting.index
        df.ending = ser_ending.index

        tm.assert_index_equal(pd.DatetimeIndex(
            df.starting), ser_starting.index)
        tm.assert_index_equal(pd.DatetimeIndex(df.ending), ser_ending.index)

    def test_is_mixed_type(self, float_frame, float_string_frame):
        assert not float_frame._is_mixed_type
        assert float_string_frame._is_mixed_type

    def test_get_numeric_data(self):
        # TODO(wesm): unused?
        intname = np.dtype(np.int_).name  # noqa
        floatname = np.dtype(np.float_).name  # noqa

        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                        'f': Timestamp('20010102')},
                       index=np.arange(10))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, 'float64': 1,
                           datetime64name: 1, objectname: 1})
        result = result.sort_index()
        expected = expected.sort_index()
        assert_series_equal(result, expected)

        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                        'd': np.array([1.] * 10, dtype='float32'),
                        'e': np.array([1] * 10, dtype='int32'),
                        'f': np.array([1] * 10, dtype='int16'),
                        'g': Timestamp('20010102')},
                       index=np.arange(10))

        result = df._get_numeric_data()
        expected = df.loc[:, ['a', 'b', 'd', 'e', 'f']]
        assert_frame_equal(result, expected)

        only_obj = df.loc[:, ['c', 'g']]
        result = only_obj._get_numeric_data()
        expected = df.loc[:, []]
        assert_frame_equal(result, expected)

        df = DataFrame.from_dict(
            {'a': [1, 2], 'b': ['foo', 'bar'], 'c': [np.pi, np.e]})
        result = df._get_numeric_data()
        expected = DataFrame.from_dict({'a': [1, 2], 'c': [np.pi, np.e]})
        assert_frame_equal(result, expected)

        df = result.copy()
        result = df._get_numeric_data()
        expected = df
        assert_frame_equal(result, expected)

    def test_get_numeric_data_extension_dtype(self):
        # GH 22290
        df = DataFrame({
            'A': integer_array([-10, np.nan, 0, 10, 20, 30], dtype='Int64'),
            'B': Categorical(list('abcabc')),
            'C': integer_array([0, 1, 2, 3, np.nan, 5], dtype='UInt8'),
            'D': IntervalArray.from_breaks(range(7))})
        result = df._get_numeric_data()
        expected = df.loc[:, ['A', 'C']]
        assert_frame_equal(result, expected)

    def test_convert_objects(self, float_string_frame):

        oops = float_string_frame.T.T
        converted = oops._convert(datetime=True)
        assert_frame_equal(converted, float_string_frame)
        assert converted['A'].dtype == np.float64

        # force numeric conversion
        float_string_frame['H'] = '1.'
        float_string_frame['I'] = '1'

        # add in some items that will be nan
        length = len(float_string_frame)
        float_string_frame['J'] = '1.'
        float_string_frame['K'] = '1'
        float_string_frame.loc[0:5, ['J', 'K']] = 'garbled'
        converted = float_string_frame._convert(datetime=True, numeric=True)
        assert converted['H'].dtype == 'float64'
        assert converted['I'].dtype == 'int64'
        assert converted['J'].dtype == 'float64'
        assert converted['K'].dtype == 'float64'
        assert len(converted['J'].dropna()) == length - 5
        assert len(converted['K'].dropna()) == length - 5

        # via astype
        converted = float_string_frame.copy()
        converted['H'] = converted['H'].astype('float64')
        converted['I'] = converted['I'].astype('int64')
        assert converted['H'].dtype == 'float64'
        assert converted['I'].dtype == 'int64'

        # via astype, but errors
        converted = float_string_frame.copy()
        with pytest.raises(ValueError, match='invalid literal'):
            converted['H'].astype('int32')

        # mixed in a single column
        df = DataFrame(dict(s=Series([1, 'na', 3, 4])))
        result = df._convert(datetime=True, numeric=True)
        expected = DataFrame(dict(s=Series([1, np.nan, 3, 4])))
        assert_frame_equal(result, expected)

    def test_convert_objects_no_conversion(self):
        mixed1 = DataFrame(
            {'a': [1, 2, 3], 'b': [4.0, 5, 6], 'c': ['x', 'y', 'z']})
        mixed2 = mixed1._convert(datetime=True)
        assert_frame_equal(mixed1, mixed2)

    def test_infer_objects(self):
        # GH 11221
        df = DataFrame({'a': ['a', 1, 2, 3],
                        'b': ['b', 2.0, 3.0, 4.1],
                        'c': ['c', datetime(2016, 1, 1),
                              datetime(2016, 1, 2),
                              datetime(2016, 1, 3)],
                        'd': [1, 2, 3, 'd']},
                       columns=['a', 'b', 'c', 'd'])
        df = df.iloc[1:].infer_objects()

        assert df['a'].dtype == 'int64'
        assert df['b'].dtype == 'float64'
        assert df['c'].dtype == 'M8[ns]'
        assert df['d'].dtype == 'object'

        expected = DataFrame({'a': [1, 2, 3],
                              'b': [2.0, 3.0, 4.1],
                              'c': [datetime(2016, 1, 1),
                                    datetime(2016, 1, 2),
                                    datetime(2016, 1, 3)],
                              'd': [2, 3, 'd']},
                             columns=['a', 'b', 'c', 'd'])
        # reconstruct frame to verify inference is same
        tm.assert_frame_equal(df.reset_index(drop=True), expected)

    def test_stale_cached_series_bug_473(self):

        # this is chained, but ok
        with option_context('chained_assignment', None):
            Y = DataFrame(np.random.random((4, 4)), index=('a', 'b', 'c', 'd'),
                          columns=('e', 'f', 'g', 'h'))
            repr(Y)
            Y['e'] = Y['e'].astype('object')
            Y['g']['c'] = np.NaN
            repr(Y)
            result = Y.sum()  # noqa
            exp = Y['g'].sum()  # noqa
            assert pd.isna(Y['g']['c'])

    def test_get_X_columns(self):
        # numeric and object columns

        df = DataFrame({'a': [1, 2, 3],
                        'b': [True, False, True],
                        'c': ['foo', 'bar', 'baz'],
                        'd': [None, None, None],
                        'e': [3.14, 0.577, 2.773]})

        tm.assert_index_equal(df._get_numeric_data().columns,
                              pd.Index(['a', 'b', 'e']))

    def test_strange_column_corruption_issue(self):
        # (wesm) Unclear how exactly this is related to internal matters
        df = DataFrame(index=[0, 1])
        df[0] = np.nan
        wasCol = {}
        # uncommenting these makes the results match
        # for col in xrange(100, 200):
        #    wasCol[col] = 1
        #    df[col] = np.nan

        for i, dt in enumerate(df.index):
            for col in range(100, 200):
                if col not in wasCol:
                    wasCol[col] = 1
                    df[col] = np.nan
                df[col][dt] = i

        myid = 100

        first = len(df.loc[pd.isna(df[myid]), [myid]])
        second = len(df.loc[pd.isna(df[myid]), [myid]])
        assert first == second == 0

    def test_constructor_no_pandas_array(self):
        # Ensure that PandasArray isn't allowed inside Series
        # See https://github.com/pandas-dev/pandas/issues/23995 for more.
        arr = pd.Series([1, 2, 3]).array
        result = pd.DataFrame({"A": arr})
        expected = pd.DataFrame({"A": [1, 2, 3]})
        tm.assert_frame_equal(result, expected)
        assert isinstance(result._data.blocks[0], IntBlock)
