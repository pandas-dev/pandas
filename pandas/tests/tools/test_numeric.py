import pytest
import decimal

import numpy as np
import pandas as pd
from pandas import to_numeric

from pandas.util import testing as tm
from numpy import iinfo


class TestToNumeric(object):

    def test_empty(self):
        # see gh-16302
        s = pd.Series([], dtype=object)

        res = to_numeric(s)
        expected = pd.Series([], dtype=np.int64)

        tm.assert_series_equal(res, expected)

        # Original issue example
        res = to_numeric(s, errors='coerce', downcast='integer')
        expected = pd.Series([], dtype=np.int8)

        tm.assert_series_equal(res, expected)

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
        with tm.assert_raises_regex(ValueError, msg):
            to_numeric(s, errors='raise')

        res = to_numeric(s, errors='ignore')
        expected = pd.Series([1, -3.14, 'apple'])
        tm.assert_series_equal(res, expected)

        res = to_numeric(s, errors='coerce')
        expected = pd.Series([1, -3.14, np.nan])
        tm.assert_series_equal(res, expected)

        s = pd.Series(['orange', 1, -3.14, 'apple'])
        msg = 'Unable to parse string "orange" at position 0'
        with tm.assert_raises_regex(ValueError, msg):
            to_numeric(s, errors='raise')

    def test_error_seen_bool(self):
        s = pd.Series([True, False, 'apple'])
        msg = 'Unable to parse string "apple" at position 2'
        with tm.assert_raises_regex(ValueError, msg):
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

    @pytest.mark.parametrize("errors", [None, "ignore", "raise", "coerce"])
    def test_type_check(self, errors):
        # see gh-11776
        df = pd.DataFrame({"a": [1, -3.14, 7], "b": ["4", "5", "6"]})
        kwargs = dict(errors=errors) if errors is not None else dict()
        error_ctx = tm.assert_raises_regex(TypeError, "1-d array")

        with error_ctx:
            to_numeric(df, **kwargs)

    def test_scalar(self):
        assert pd.to_numeric(1) == 1
        assert pd.to_numeric(1.1) == 1.1

        assert pd.to_numeric('1') == 1
        assert pd.to_numeric('1.1') == 1.1

        with pytest.raises(ValueError):
            to_numeric('XX', errors='raise')

        assert to_numeric('XX', errors='ignore') == 'XX'
        assert np.isnan(to_numeric('XX', errors='coerce'))

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

    def test_datetime_like(self, tz_naive_fixture):
        idx = pd.date_range("20130101", periods=3,
                            tz=tz_naive_fixture, name="xxx")
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(idx.asi8, name="xxx"))

        res = pd.to_numeric(pd.Series(idx, name="xxx"))
        tm.assert_series_equal(res, pd.Series(idx.asi8, name="xxx"))

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

        # TODO: enable when we can support native PeriodDtype
        # res = pd.to_numeric(pd.Series(idx, name='xxx'))
        # tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

    def test_non_hashable(self):
        # Test for Bug #13324
        s = pd.Series([[10.0, 2], 1.0, 'apple'])
        res = pd.to_numeric(s, errors='coerce')
        tm.assert_series_equal(res, pd.Series([np.nan, 1.0, np.nan]))

        res = pd.to_numeric(s, errors='ignore')
        tm.assert_series_equal(res, pd.Series([[10.0, 2], 1.0, 'apple']))

        with tm.assert_raises_regex(TypeError, "Invalid object type"):
            pd.to_numeric(s)

    @pytest.mark.parametrize("data", [
        ["1", 2, 3],
        [1, 2, 3],
        np.array(["1970-01-02", "1970-01-03",
                  "1970-01-04"], dtype="datetime64[D]")
    ])
    def test_downcast_basic(self, data):
        # see gh-13352
        invalid_downcast = "unsigned-integer"
        msg = "invalid downcasting method provided"

        with tm.assert_raises_regex(ValueError, msg):
            pd.to_numeric(data, downcast=invalid_downcast)

        expected = np.array([1, 2, 3], dtype=np.int64)

        # Basic function tests.
        res = pd.to_numeric(data)
        tm.assert_numpy_array_equal(res, expected)

        res = pd.to_numeric(data, downcast=None)
        tm.assert_numpy_array_equal(res, expected)

        # Basic dtype support.
        smallest_uint_dtype = np.dtype(np.typecodes["UnsignedInteger"][0])

        # Support below np.float32 is rare and far between.
        float_32_char = np.dtype(np.float32).char
        smallest_float_dtype = float_32_char

        expected = np.array([1, 2, 3], dtype=smallest_uint_dtype)
        res = pd.to_numeric(data, downcast="unsigned")
        tm.assert_numpy_array_equal(res, expected)

        expected = np.array([1, 2, 3], dtype=smallest_float_dtype)
        res = pd.to_numeric(data, downcast="float")
        tm.assert_numpy_array_equal(res, expected)

    @pytest.mark.parametrize("signed_downcast", ["integer", "signed"])
    @pytest.mark.parametrize("data", [
        ["1", 2, 3],
        [1, 2, 3],
        np.array(["1970-01-02", "1970-01-03",
                  "1970-01-04"], dtype="datetime64[D]")
    ])
    def test_signed_downcast(self, data, signed_downcast):
        # see gh-13352
        smallest_int_dtype = np.dtype(np.typecodes["Integer"][0])
        expected = np.array([1, 2, 3], dtype=smallest_int_dtype)

        res = pd.to_numeric(data, downcast=signed_downcast)
        tm.assert_numpy_array_equal(res, expected)

    def test_ignore_downcast_invalid_data(self):
        # If we can't successfully cast the given
        # data to a numeric dtype, do not bother
        # with the downcast parameter.
        data = ["foo", 2, 3]
        expected = np.array(data, dtype=object)

        res = pd.to_numeric(data, errors="ignore",
                            downcast="unsigned")
        tm.assert_numpy_array_equal(res, expected)

    def test_ignore_downcast_neg_to_unsigned(self):
        # Cannot cast to an unsigned integer
        # because we have a negative number.
        data = ["-1", 2, 3]
        expected = np.array([-1, 2, 3], dtype=np.int64)

        res = pd.to_numeric(data, downcast="unsigned")
        tm.assert_numpy_array_equal(res, expected)

    @pytest.mark.parametrize("downcast", ["integer", "signed", "unsigned"])
    @pytest.mark.parametrize("data,expected", [
        (["1.1", 2, 3],
         np.array([1.1, 2, 3], dtype=np.float64)),
        ([10000.0, 20000, 3000, 40000.36, 50000, 50000.00],
         np.array([10000.0, 20000, 3000,
                   40000.36, 50000, 50000.00], dtype=np.float64))
    ])
    def test_ignore_downcast_cannot_convert_float(
            self, data, expected, downcast):
        # Cannot cast to an integer (signed or unsigned)
        # because we have a float number.
        res = pd.to_numeric(data, downcast=downcast)
        tm.assert_numpy_array_equal(res, expected)

    @pytest.mark.parametrize("downcast,expected_dtype", [
        ("integer", np.int16),
        ("signed", np.int16),
        ("unsigned", np.uint16)
    ])
    def test_downcast_not8bit(self, downcast, expected_dtype):
        # the smallest integer dtype need not be np.(u)int8
        data = ["256", 257, 258]

        expected = np.array([256, 257, 258], dtype=expected_dtype)
        res = pd.to_numeric(data, downcast=downcast)
        tm.assert_numpy_array_equal(res, expected)

    @pytest.mark.parametrize("dtype,downcast,min_max", [
        ("int8", "integer", [iinfo(np.int8).min,
                             iinfo(np.int8).max]),
        ("int16", "integer", [iinfo(np.int16).min,
                              iinfo(np.int16).max]),
        ('int32', "integer", [iinfo(np.int32).min,
                              iinfo(np.int32).max]),
        ('int64', "integer", [iinfo(np.int64).min,
                              iinfo(np.int64).max]),
        ('uint8', "unsigned", [iinfo(np.uint8).min,
                               iinfo(np.uint8).max]),
        ('uint16', "unsigned", [iinfo(np.uint16).min,
                                iinfo(np.uint16).max]),
        ('uint32', "unsigned", [iinfo(np.uint32).min,
                                iinfo(np.uint32).max]),
        ('uint64', "unsigned", [iinfo(np.uint64).min,
                                iinfo(np.uint64).max]),
        ('int16', "integer", [iinfo(np.int8).min,
                              iinfo(np.int8).max + 1]),
        ('int32', "integer", [iinfo(np.int16).min,
                              iinfo(np.int16).max + 1]),
        ('int64', "integer", [iinfo(np.int32).min,
                              iinfo(np.int32).max + 1]),
        ('int16', "integer", [iinfo(np.int8).min - 1,
                              iinfo(np.int16).max]),
        ('int32', "integer", [iinfo(np.int16).min - 1,
                              iinfo(np.int32).max]),
        ('int64', "integer", [iinfo(np.int32).min - 1,
                              iinfo(np.int64).max]),
        ('uint16', "unsigned", [iinfo(np.uint8).min,
                                iinfo(np.uint8).max + 1]),
        ('uint32', "unsigned", [iinfo(np.uint16).min,
                                iinfo(np.uint16).max + 1]),
        ('uint64', "unsigned", [iinfo(np.uint32).min,
                                iinfo(np.uint32).max + 1])
    ])
    def test_downcast_limits(self, dtype, downcast, min_max):
        # see gh-14404: test the limits of each downcast.
        series = pd.to_numeric(pd.Series(min_max), downcast=downcast)
        assert series.dtype == dtype

    def test_coerce_uint64_conflict(self):
        # see gh-17007 and gh-17125
        #
        # Still returns float despite the uint64-nan conflict,
        # which would normally force the casting to object.
        df = pd.DataFrame({"a": [200, 300, "", "NaN", 30000000000000000000]})
        expected = pd.Series([200, 300, np.nan, np.nan,
                              30000000000000000000], dtype=float, name="a")
        result = to_numeric(df["a"], errors="coerce")
        tm.assert_series_equal(result, expected)

        s = pd.Series(["12345678901234567890", "1234567890", "ITEM"])
        expected = pd.Series([12345678901234567890,
                              1234567890, np.nan], dtype=float)
        result = to_numeric(s, errors="coerce")
        tm.assert_series_equal(result, expected)

        # For completeness, check against "ignore" and "raise"
        result = to_numeric(s, errors="ignore")
        tm.assert_series_equal(result, s)

        msg = "Unable to parse string"
        with tm.assert_raises_regex(ValueError, msg):
            to_numeric(s, errors="raise")
