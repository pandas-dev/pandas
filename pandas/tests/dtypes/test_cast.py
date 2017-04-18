# -*- coding: utf-8 -*-

"""
These test the private routines in types/cast.py

"""

import pytest
import decimal
from datetime import datetime, timedelta, date
import numpy as np

import pandas as pd
from pandas import (Timedelta, Timestamp, DatetimeIndex,
                    to_numeric, _np_version_under1p9)

from pandas.core.dtypes.cast import (
    maybe_downcast_to_dtype,
    maybe_convert_objects,
    infer_dtype_from_scalar,
    infer_dtype_from_array,
    maybe_convert_string_to_object,
    maybe_convert_scalar,
    find_common_type)
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    DatetimeTZDtype,
    PeriodDtype)
from pandas.util import testing as tm

from numpy import iinfo


class TestMaybeDowncast(tm.TestCase):

    def test_downcast_conv(self):
        # test downcasting

        arr = np.array([8.5, 8.6, 8.7, 8.8, 8.9999999999995])
        result = maybe_downcast_to_dtype(arr, 'infer')
        assert (np.array_equal(result, arr))

        arr = np.array([8., 8., 8., 8., 8.9999999999995])
        result = maybe_downcast_to_dtype(arr, 'infer')
        expected = np.array([8, 8, 8, 8, 9])
        assert (np.array_equal(result, expected))

        arr = np.array([8., 8., 8., 8., 9.0000000000005])
        result = maybe_downcast_to_dtype(arr, 'infer')
        expected = np.array([8, 8, 8, 8, 9])
        assert (np.array_equal(result, expected))

        # conversions

        expected = np.array([1, 2])
        for dtype in [np.float64, object, np.int64]:
            arr = np.array([1.0, 2.0], dtype=dtype)
            result = maybe_downcast_to_dtype(arr, 'infer')
            tm.assert_almost_equal(result, expected, check_dtype=False)

        for dtype in [np.float64, object]:
            expected = np.array([1.0, 2.0, np.nan], dtype=dtype)
            arr = np.array([1.0, 2.0, np.nan], dtype=dtype)
            result = maybe_downcast_to_dtype(arr, 'infer')
            tm.assert_almost_equal(result, expected)

        # empties
        for dtype in [np.int32, np.float64, np.float32, np.bool_,
                      np.int64, object]:
            arr = np.array([], dtype=dtype)
            result = maybe_downcast_to_dtype(arr, 'int64')
            tm.assert_almost_equal(result, np.array([], dtype=np.int64))
            assert result.dtype == np.int64

    def test_datetimelikes_nan(self):
        arr = np.array([1, 2, np.nan])
        exp = np.array([1, 2, np.datetime64('NaT')], dtype='datetime64[ns]')
        res = maybe_downcast_to_dtype(arr, 'datetime64[ns]')
        tm.assert_numpy_array_equal(res, exp)

        exp = np.array([1, 2, np.timedelta64('NaT')], dtype='timedelta64[ns]')
        res = maybe_downcast_to_dtype(arr, 'timedelta64[ns]')
        tm.assert_numpy_array_equal(res, exp)

    def test_datetime_with_timezone(self):
        # GH 15426
        ts = Timestamp("2016-01-01 12:00:00", tz='US/Pacific')
        exp = DatetimeIndex([ts, ts])
        res = maybe_downcast_to_dtype(exp, exp.dtype)
        tm.assert_index_equal(res, exp)

        res = maybe_downcast_to_dtype(exp.asi8, exp.dtype)
        tm.assert_index_equal(res, exp)


class TestInferDtype(object):

    def test_infer_dtype_from_scalar(self):
        # Test that _infer_dtype_from_scalar is returning correct dtype for int
        # and float.

        for dtypec in [np.uint8, np.int8, np.uint16, np.int16, np.uint32,
                       np.int32, np.uint64, np.int64]:
            data = dtypec(12)
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == type(data)

        data = 12
        dtype, val = infer_dtype_from_scalar(data)
        assert dtype == np.int64

        for dtypec in [np.float16, np.float32, np.float64]:
            data = dtypec(12)
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == dtypec

        data = np.float(12)
        dtype, val = infer_dtype_from_scalar(data)
        assert dtype == np.float64

        for data in [True, False]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == np.bool_

        for data in [np.complex64(1), np.complex128(1)]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == np.complex_

        for data in [np.datetime64(1, 'ns'), Timestamp(1),
                     datetime(2000, 1, 1, 0, 0)]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == 'M8[ns]'

        for data in [np.timedelta64(1, 'ns'), Timedelta(1),
                     timedelta(1)]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == 'm8[ns]'

        for data in [date(2000, 1, 1),
                     Timestamp(1, tz='US/Eastern'), 'foo']:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == np.object_

    @pytest.mark.parametrize(
        "arr, expected",
        [('foo', np.object_),
         (b'foo', np.object_),
         (1, np.int_),
         (1.5, np.float_),
         ([1], np.int_),
         (np.array([1]), np.int_),
         ([np.nan, 1, ''], np.object_),
         (np.array([[1.0, 2.0]]), np.float_),
         (Timestamp('20160101'), np.object_),
         (np.datetime64('2016-01-01'), np.dtype('<M8[D]')),
         ])
    def test_infer_dtype_from_array(self, arr, expected):

        # these infer specifically to numpy dtypes
        dtype, _ = infer_dtype_from_array(arr)
        assert dtype == expected


class TestMaybe(tm.TestCase):

    def test_maybe_convert_string_to_array(self):
        result = maybe_convert_string_to_object('x')
        tm.assert_numpy_array_equal(result, np.array(['x'], dtype=object))
        self.assertTrue(result.dtype == object)

        result = maybe_convert_string_to_object(1)
        self.assertEqual(result, 1)

        arr = np.array(['x', 'y'], dtype=str)
        result = maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 'y'], dtype=object))
        self.assertTrue(result.dtype == object)

        # unicode
        arr = np.array(['x', 'y']).astype('U')
        result = maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 'y'], dtype=object))
        self.assertTrue(result.dtype == object)

        # object
        arr = np.array(['x', 2], dtype=object)
        result = maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 2], dtype=object))
        self.assertTrue(result.dtype == object)

    def test_maybe_convert_scalar(self):

        # pass thru
        result = maybe_convert_scalar('x')
        self.assertEqual(result, 'x')
        result = maybe_convert_scalar(np.array([1]))
        self.assertEqual(result, np.array([1]))

        # leave scalar dtype
        result = maybe_convert_scalar(np.int64(1))
        self.assertEqual(result, np.int64(1))
        result = maybe_convert_scalar(np.int32(1))
        self.assertEqual(result, np.int32(1))
        result = maybe_convert_scalar(np.float32(1))
        self.assertEqual(result, np.float32(1))
        result = maybe_convert_scalar(np.int64(1))
        self.assertEqual(result, np.float64(1))

        # coerce
        result = maybe_convert_scalar(1)
        self.assertEqual(result, np.int64(1))
        result = maybe_convert_scalar(1.0)
        self.assertEqual(result, np.float64(1))
        result = maybe_convert_scalar(Timestamp('20130101'))
        self.assertEqual(result, Timestamp('20130101').value)
        result = maybe_convert_scalar(datetime(2013, 1, 1))
        self.assertEqual(result, Timestamp('20130101').value)
        result = maybe_convert_scalar(Timedelta('1 day 1 min'))
        self.assertEqual(result, Timedelta('1 day 1 min').value)


class TestConvert(tm.TestCase):

    def test_maybe_convert_objects_copy(self):
        values = np.array([1, 2])

        out = maybe_convert_objects(values, copy=False)
        self.assertTrue(values is out)

        out = maybe_convert_objects(values, copy=True)
        self.assertTrue(values is not out)

        values = np.array(['apply', 'banana'])
        out = maybe_convert_objects(values, copy=False)
        self.assertTrue(values is out)

        out = maybe_convert_objects(values, copy=True)
        self.assertTrue(values is not out)


class TestCommonTypes(tm.TestCase):

    def test_numpy_dtypes(self):
        # (source_types, destination_type)
        testcases = (
            # identity
            ((np.int64,), np.int64),
            ((np.uint64,), np.uint64),
            ((np.float32,), np.float32),
            ((np.object,), np.object),

            # into ints
            ((np.int16, np.int64), np.int64),
            ((np.int32, np.uint32), np.int64),
            ((np.uint16, np.uint64), np.uint64),

            # into floats
            ((np.float16, np.float32), np.float32),
            ((np.float16, np.int16), np.float32),
            ((np.float32, np.int16), np.float32),
            ((np.uint64, np.int64), np.float64),
            ((np.int16, np.float64), np.float64),
            ((np.float16, np.int64), np.float64),

            # into others
            ((np.complex128, np.int32), np.complex128),
            ((np.object, np.float32), np.object),
            ((np.object, np.int16), np.object),

            # bool with int
            ((np.dtype('bool'), np.int64), np.object),
            ((np.dtype('bool'), np.int32), np.object),
            ((np.dtype('bool'), np.int16), np.object),
            ((np.dtype('bool'), np.int8), np.object),
            ((np.dtype('bool'), np.uint64), np.object),
            ((np.dtype('bool'), np.uint32), np.object),
            ((np.dtype('bool'), np.uint16), np.object),
            ((np.dtype('bool'), np.uint8), np.object),

            # bool with float
            ((np.dtype('bool'), np.float64), np.object),
            ((np.dtype('bool'), np.float32), np.object),

            ((np.dtype('datetime64[ns]'), np.dtype('datetime64[ns]')),
             np.dtype('datetime64[ns]')),
            ((np.dtype('timedelta64[ns]'), np.dtype('timedelta64[ns]')),
             np.dtype('timedelta64[ns]')),

            ((np.dtype('datetime64[ns]'), np.dtype('datetime64[ms]')),
             np.dtype('datetime64[ns]')),
            ((np.dtype('timedelta64[ms]'), np.dtype('timedelta64[ns]')),
             np.dtype('timedelta64[ns]')),

            ((np.dtype('datetime64[ns]'), np.dtype('timedelta64[ns]')),
             np.object),
            ((np.dtype('datetime64[ns]'), np.int64), np.object)
        )
        for src, common in testcases:
            self.assertEqual(find_common_type(src), common)

        with tm.assertRaises(ValueError):
            # empty
            find_common_type([])

    def test_categorical_dtype(self):
        dtype = CategoricalDtype()
        self.assertEqual(find_common_type([dtype]), 'category')
        self.assertEqual(find_common_type([dtype, dtype]), 'category')
        self.assertEqual(find_common_type([np.object, dtype]), np.object)

    def test_datetimetz_dtype(self):
        dtype = DatetimeTZDtype(unit='ns', tz='US/Eastern')
        self.assertEqual(find_common_type([dtype, dtype]),
                         'datetime64[ns, US/Eastern]')

        for dtype2 in [DatetimeTZDtype(unit='ns', tz='Asia/Tokyo'),
                       np.dtype('datetime64[ns]'), np.object, np.int64]:
            self.assertEqual(find_common_type([dtype, dtype2]), np.object)
            self.assertEqual(find_common_type([dtype2, dtype]), np.object)

    def test_period_dtype(self):
        dtype = PeriodDtype(freq='D')
        self.assertEqual(find_common_type([dtype, dtype]), 'period[D]')

        for dtype2 in [DatetimeTZDtype(unit='ns', tz='Asia/Tokyo'),
                       PeriodDtype(freq='2D'), PeriodDtype(freq='H'),
                       np.dtype('datetime64[ns]'), np.object, np.int64]:
            self.assertEqual(find_common_type([dtype, dtype2]), np.object)
            self.assertEqual(find_common_type([dtype2, dtype]), np.object)


class TestToNumeric(tm.TestCase):

    def test_series(self):
        s = pd.Series(['1', '-3.14', '7'])
        res = to_numeric(s)
        expected = pd.Series([1, -3.14, 7])
        tm.assert_series_equal(res, expected)

        s = pd.Series(['1', '-3.14', 7])
        res = to_numeric(s)
        tm.assert_series_equal(res, expected)

    def test_series_numeric(self):
        s = pd.Series([1, 3, 4, 5], index=list('ABCD'), name='XXX')
        res = to_numeric(s)
        tm.assert_series_equal(res, s)

        s = pd.Series([1., 3., 4., 5.], index=list('ABCD'), name='XXX')
        res = to_numeric(s)
        tm.assert_series_equal(res, s)

        # bool is regarded as numeric
        s = pd.Series([True, False, True, True],
                      index=list('ABCD'), name='XXX')
        res = to_numeric(s)
        tm.assert_series_equal(res, s)

    def test_error(self):
        s = pd.Series([1, -3.14, 'apple'])
        msg = 'Unable to parse string "apple" at position 2'
        with tm.assertRaisesRegexp(ValueError, msg):
            to_numeric(s, errors='raise')

        res = to_numeric(s, errors='ignore')
        expected = pd.Series([1, -3.14, 'apple'])
        tm.assert_series_equal(res, expected)

        res = to_numeric(s, errors='coerce')
        expected = pd.Series([1, -3.14, np.nan])
        tm.assert_series_equal(res, expected)

        s = pd.Series(['orange', 1, -3.14, 'apple'])
        msg = 'Unable to parse string "orange" at position 0'
        with tm.assertRaisesRegexp(ValueError, msg):
            to_numeric(s, errors='raise')

    def test_error_seen_bool(self):
        s = pd.Series([True, False, 'apple'])
        msg = 'Unable to parse string "apple" at position 2'
        with tm.assertRaisesRegexp(ValueError, msg):
            to_numeric(s, errors='raise')

        res = to_numeric(s, errors='ignore')
        expected = pd.Series([True, False, 'apple'])
        tm.assert_series_equal(res, expected)

        # coerces to float
        res = to_numeric(s, errors='coerce')
        expected = pd.Series([1., 0., np.nan])
        tm.assert_series_equal(res, expected)

    def test_list(self):
        s = ['1', '-3.14', '7']
        res = to_numeric(s)
        expected = np.array([1, -3.14, 7])
        tm.assert_numpy_array_equal(res, expected)

    def test_list_numeric(self):
        s = [1, 3, 4, 5]
        res = to_numeric(s)
        tm.assert_numpy_array_equal(res, np.array(s, dtype=np.int64))

        s = [1., 3., 4., 5.]
        res = to_numeric(s)
        tm.assert_numpy_array_equal(res, np.array(s))

        # bool is regarded as numeric
        s = [True, False, True, True]
        res = to_numeric(s)
        tm.assert_numpy_array_equal(res, np.array(s))

    def test_numeric(self):
        s = pd.Series([1, -3.14, 7], dtype='O')
        res = to_numeric(s)
        expected = pd.Series([1, -3.14, 7])
        tm.assert_series_equal(res, expected)

        s = pd.Series([1, -3.14, 7])
        res = to_numeric(s)
        tm.assert_series_equal(res, expected)

        # GH 14827
        df = pd.DataFrame(dict(
            a=[1.2, decimal.Decimal(3.14), decimal.Decimal("infinity"), '0.1'],
            b=[1.0, 2.0, 3.0, 4.0],
        ))
        expected = pd.DataFrame(dict(
            a=[1.2, 3.14, np.inf, 0.1],
            b=[1.0, 2.0, 3.0, 4.0],
        ))

        # Test to_numeric over one column
        df_copy = df.copy()
        df_copy['a'] = df_copy['a'].apply(to_numeric)
        tm.assert_frame_equal(df_copy, expected)

        # Test to_numeric over multiple columns
        df_copy = df.copy()
        df_copy[['a', 'b']] = df_copy[['a', 'b']].apply(to_numeric)
        tm.assert_frame_equal(df_copy, expected)

    def test_numeric_lists_and_arrays(self):
        # Test to_numeric with embedded lists and arrays
        df = pd.DataFrame(dict(
            a=[[decimal.Decimal(3.14), 1.0], decimal.Decimal(1.6), 0.1]
        ))
        df['a'] = df['a'].apply(to_numeric)
        expected = pd.DataFrame(dict(
            a=[[3.14, 1.0], 1.6, 0.1],
        ))
        tm.assert_frame_equal(df, expected)

        df = pd.DataFrame(dict(
            a=[np.array([decimal.Decimal(3.14), 1.0]), 0.1]
        ))
        df['a'] = df['a'].apply(to_numeric)
        expected = pd.DataFrame(dict(
            a=[[3.14, 1.0], 0.1],
        ))
        tm.assert_frame_equal(df, expected)

    def test_all_nan(self):
        s = pd.Series(['a', 'b', 'c'])
        res = to_numeric(s, errors='coerce')
        expected = pd.Series([np.nan, np.nan, np.nan])
        tm.assert_series_equal(res, expected)

    def test_type_check(self):
        # GH 11776
        df = pd.DataFrame({'a': [1, -3.14, 7], 'b': ['4', '5', '6']})
        with tm.assertRaisesRegexp(TypeError, "1-d array"):
            to_numeric(df)
        for errors in ['ignore', 'raise', 'coerce']:
            with tm.assertRaisesRegexp(TypeError, "1-d array"):
                to_numeric(df, errors=errors)

    def test_scalar(self):
        self.assertEqual(pd.to_numeric(1), 1)
        self.assertEqual(pd.to_numeric(1.1), 1.1)

        self.assertEqual(pd.to_numeric('1'), 1)
        self.assertEqual(pd.to_numeric('1.1'), 1.1)

        with tm.assertRaises(ValueError):
            to_numeric('XX', errors='raise')

        self.assertEqual(to_numeric('XX', errors='ignore'), 'XX')
        self.assertTrue(np.isnan(to_numeric('XX', errors='coerce')))

    def test_numeric_dtypes(self):
        idx = pd.Index([1, 2, 3], name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, idx)

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(idx, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, idx.values)

        idx = pd.Index([1., np.nan, 3., np.nan], name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, idx)

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(idx, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, idx.values)

    def test_str(self):
        idx = pd.Index(['1', '2', '3'], name='xxx')
        exp = np.array([1, 2, 3], dtype='int64')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(exp, name='xxx'))

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(exp, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, exp)

        idx = pd.Index(['1.5', '2.7', '3.4'], name='xxx')
        exp = np.array([1.5, 2.7, 3.4])
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(exp, name='xxx'))

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(exp, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, exp)

    def test_datetimelike(self):
        for tz in [None, 'US/Eastern', 'Asia/Tokyo']:
            idx = pd.date_range('20130101', periods=3, tz=tz, name='xxx')
            res = pd.to_numeric(idx)
            tm.assert_index_equal(res, pd.Index(idx.asi8, name='xxx'))

            res = pd.to_numeric(pd.Series(idx, name='xxx'))
            tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

            res = pd.to_numeric(idx.values)
            tm.assert_numpy_array_equal(res, idx.asi8)

    def test_timedelta(self):
        idx = pd.timedelta_range('1 days', periods=3, freq='D', name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(idx.asi8, name='xxx'))

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, idx.asi8)

    def test_period(self):
        idx = pd.period_range('2011-01', periods=3, freq='M', name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(idx.asi8, name='xxx'))

        # ToDo: enable when we can support native PeriodDtype
        # res = pd.to_numeric(pd.Series(idx, name='xxx'))
        # tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

    def test_non_hashable(self):
        # Test for Bug #13324
        s = pd.Series([[10.0, 2], 1.0, 'apple'])
        res = pd.to_numeric(s, errors='coerce')
        tm.assert_series_equal(res, pd.Series([np.nan, 1.0, np.nan]))

        res = pd.to_numeric(s, errors='ignore')
        tm.assert_series_equal(res, pd.Series([[10.0, 2], 1.0, 'apple']))

        with self.assertRaisesRegexp(TypeError, "Invalid object type"):
            pd.to_numeric(s)

    def test_downcast(self):
        # see gh-13352
        mixed_data = ['1', 2, 3]
        int_data = [1, 2, 3]
        date_data = np.array(['1970-01-02', '1970-01-03',
                              '1970-01-04'], dtype='datetime64[D]')

        invalid_downcast = 'unsigned-integer'
        msg = 'invalid downcasting method provided'

        smallest_int_dtype = np.dtype(np.typecodes['Integer'][0])
        smallest_uint_dtype = np.dtype(np.typecodes['UnsignedInteger'][0])

        # support below np.float32 is rare and far between
        float_32_char = np.dtype(np.float32).char
        smallest_float_dtype = float_32_char

        for data in (mixed_data, int_data, date_data):
            with self.assertRaisesRegexp(ValueError, msg):
                pd.to_numeric(data, downcast=invalid_downcast)

            expected = np.array([1, 2, 3], dtype=np.int64)

            res = pd.to_numeric(data)
            tm.assert_numpy_array_equal(res, expected)

            res = pd.to_numeric(data, downcast=None)
            tm.assert_numpy_array_equal(res, expected)

            expected = np.array([1, 2, 3], dtype=smallest_int_dtype)

            for signed_downcast in ('integer', 'signed'):
                res = pd.to_numeric(data, downcast=signed_downcast)
                tm.assert_numpy_array_equal(res, expected)

            expected = np.array([1, 2, 3], dtype=smallest_uint_dtype)
            res = pd.to_numeric(data, downcast='unsigned')
            tm.assert_numpy_array_equal(res, expected)

            expected = np.array([1, 2, 3], dtype=smallest_float_dtype)
            res = pd.to_numeric(data, downcast='float')
            tm.assert_numpy_array_equal(res, expected)

        # if we can't successfully cast the given
        # data to a numeric dtype, do not bother
        # with the downcast parameter
        data = ['foo', 2, 3]
        expected = np.array(data, dtype=object)
        res = pd.to_numeric(data, errors='ignore',
                            downcast='unsigned')
        tm.assert_numpy_array_equal(res, expected)

        # cannot cast to an unsigned integer because
        # we have a negative number
        data = ['-1', 2, 3]
        expected = np.array([-1, 2, 3], dtype=np.int64)
        res = pd.to_numeric(data, downcast='unsigned')
        tm.assert_numpy_array_equal(res, expected)

        # cannot cast to an integer (signed or unsigned)
        # because we have a float number
        data = (['1.1', 2, 3],
                [10000.0, 20000, 3000, 40000.36, 50000, 50000.00])
        expected = (np.array([1.1, 2, 3], dtype=np.float64),
                    np.array([10000.0, 20000, 3000,
                              40000.36, 50000, 50000.00], dtype=np.float64))

        for _data, _expected in zip(data, expected):
            for downcast in ('integer', 'signed', 'unsigned'):
                res = pd.to_numeric(_data, downcast=downcast)
                tm.assert_numpy_array_equal(res, _expected)

        # the smallest integer dtype need not be np.(u)int8
        data = ['256', 257, 258]

        for downcast, expected_dtype in zip(
                ['integer', 'signed', 'unsigned'],
                [np.int16, np.int16, np.uint16]):
            expected = np.array([256, 257, 258], dtype=expected_dtype)
            res = pd.to_numeric(data, downcast=downcast)
            tm.assert_numpy_array_equal(res, expected)

    def test_downcast_limits(self):
        # Test the limits of each downcast. Bug: #14401.
        # Check to make sure numpy is new enough to run this test.
        if _np_version_under1p9:
            pytest.skip("Numpy version is under 1.9")

        i = 'integer'
        u = 'unsigned'
        dtype_downcast_min_max = [
            ('int8', i, [iinfo(np.int8).min, iinfo(np.int8).max]),
            ('int16', i, [iinfo(np.int16).min, iinfo(np.int16).max]),
            ('int32', i, [iinfo(np.int32).min, iinfo(np.int32).max]),
            ('int64', i, [iinfo(np.int64).min, iinfo(np.int64).max]),
            ('uint8', u, [iinfo(np.uint8).min, iinfo(np.uint8).max]),
            ('uint16', u, [iinfo(np.uint16).min, iinfo(np.uint16).max]),
            ('uint32', u, [iinfo(np.uint32).min, iinfo(np.uint32).max]),
            ('uint64', u, [iinfo(np.uint64).min, iinfo(np.uint64).max]),
            ('int16', i, [iinfo(np.int8).min, iinfo(np.int8).max + 1]),
            ('int32', i, [iinfo(np.int16).min, iinfo(np.int16).max + 1]),
            ('int64', i, [iinfo(np.int32).min, iinfo(np.int32).max + 1]),
            ('int16', i, [iinfo(np.int8).min - 1, iinfo(np.int16).max]),
            ('int32', i, [iinfo(np.int16).min - 1, iinfo(np.int32).max]),
            ('int64', i, [iinfo(np.int32).min - 1, iinfo(np.int64).max]),
            ('uint16', u, [iinfo(np.uint8).min, iinfo(np.uint8).max + 1]),
            ('uint32', u, [iinfo(np.uint16).min, iinfo(np.uint16).max + 1]),
            ('uint64', u, [iinfo(np.uint32).min, iinfo(np.uint32).max + 1])
        ]

        for dtype, downcast, min_max in dtype_downcast_min_max:
            series = pd.to_numeric(pd.Series(min_max), downcast=downcast)
            assert series.dtype == dtype
