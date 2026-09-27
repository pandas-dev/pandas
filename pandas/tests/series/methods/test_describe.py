import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas.core.dtypes.common import (
    is_complex_dtype,
    is_extension_array_dtype,
)

import pandas as pd
import pandas._testing as tm


class TestSeriesDescribe:
    def test_describe_ints(self):
        ser = pd.Series([0, 1, 2, 3, 4], name="int_data")
        result = ser.describe()
        expected = pd.Series(
            [5, 2, ser.std(), 0, 1, 2, 3, 4],
            name="int_data",
            index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

    def test_describe_bools(self):
        ser = pd.Series([True, True, False, False, False], name="bool_data")
        result = ser.describe()
        expected = pd.Series(
            [5, 2, False, 3], name="bool_data", index=["count", "unique", "top", "freq"]
        )
        tm.assert_series_equal(result, expected)

    def test_describe_strs(self):
        ser = pd.Series(["a", "a", "b", "c", "d"], name="str_data")
        result = ser.describe()
        expected = pd.Series(
            [5, 4, "a", 2], name="str_data", index=["count", "unique", "top", "freq"]
        )
        tm.assert_series_equal(result, expected)

    def test_describe_timedelta64(self):
        ser = pd.Series(
            [
                pd.Timedelta("1 days"),
                pd.Timedelta("2 days"),
                pd.Timedelta("3 days"),
                pd.Timedelta("4 days"),
                pd.Timedelta("5 days"),
            ],
            name="timedelta_data",
        )
        result = ser.describe()
        expected = pd.Series(
            [5, ser[2], ser.std(), ser[0], ser[1], ser[2], ser[3], ser[4]],
            name="timedelta_data",
            index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

    def test_describe_period(self):
        ser = pd.Series(
            [
                pd.Period("2020-01", "M"),
                pd.Period("2020-01", "M"),
                pd.Period("2019-12", "M"),
            ],
            name="period_data",
        )
        result = ser.describe()
        expected = pd.Series(
            [3, 2, ser[0], 2],
            name="period_data",
            index=["count", "unique", "top", "freq"],
        )
        tm.assert_series_equal(result, expected)

    def test_describe_empty_object(self):
        # https://github.com/pandas-dev/pandas/issues/27183
        s = pd.Series([None, None], dtype=object)
        result = s.describe()
        expected = pd.Series(
            [0, 0, np.nan, np.nan],
            dtype=object,
            index=["count", "unique", "top", "freq"],
        )
        tm.assert_series_equal(result, expected)

        result = s[:0].describe()
        tm.assert_series_equal(result, expected)
        # ensure NaN, not None
        assert np.isnan(result.iloc[2])
        assert np.isnan(result.iloc[3])

    def test_describe_with_tz(self, tz_naive_fixture):
        # GH 21332
        tz = tz_naive_fixture
        name = str(tz_naive_fixture)
        start = pd.Timestamp(2018, 1, 1)
        end = pd.Timestamp(2018, 1, 5)
        s = pd.Series(pd.date_range(start, end, tz=tz), name=name)
        result = s.describe()
        expected = pd.Series(
            [
                5,
                pd.Timestamp(2018, 1, 3).tz_localize(tz),
                start.tz_localize(tz),
                s[1],
                s[2],
                s[3],
                end.tz_localize(tz),
            ],
            name=name,
            index=["count", "mean", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

    def test_describe_with_tz_numeric(self):
        name = tz = "CET"
        start = pd.Timestamp(2018, 1, 1)
        end = pd.Timestamp(2018, 1, 5)
        s = pd.Series(pd.date_range(start, end, tz=tz), name=name)

        result = s.describe()

        expected = pd.Series(
            [
                5,
                pd.Timestamp("2018-01-03 00:00:00", tz=tz),
                pd.Timestamp("2018-01-01 00:00:00", tz=tz),
                pd.Timestamp("2018-01-02 00:00:00", tz=tz),
                pd.Timestamp("2018-01-03 00:00:00", tz=tz),
                pd.Timestamp("2018-01-04 00:00:00", tz=tz),
                pd.Timestamp("2018-01-05 00:00:00", tz=tz),
            ],
            name=name,
            index=["count", "mean", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

    def test_datetime_is_numeric_includes_datetime(self):
        s = pd.Series(pd.date_range("2012", periods=3))
        result = s.describe()
        expected = pd.Series(
            [
                3,
                pd.Timestamp("2012-01-02"),
                pd.Timestamp("2012-01-01"),
                pd.Timestamp("2012-01-01T12:00:00"),
                pd.Timestamp("2012-01-02"),
                pd.Timestamp("2012-01-02T12:00:00"),
                pd.Timestamp("2012-01-03"),
            ],
            index=["count", "mean", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

    def test_numeric_result_dtype(self, any_numeric_dtype):
        # GH#48340 - describe should always return float on non-complex numeric input
        if is_extension_array_dtype(any_numeric_dtype):
            dtype = "Float64"
        else:
            dtype = "complex128" if is_complex_dtype(any_numeric_dtype) else None

        ser = pd.Series([0, 1], dtype=any_numeric_dtype)
        if dtype == "complex128":
            with pytest.raises(
                TypeError, match=r"^a must be an array of real numbers$"
            ):
                ser.describe()
            return
        result = ser.describe()
        expected = pd.Series(
            [
                2.0,
                0.5,
                ser.std(),
                0,
                0.25,
                0.5,
                0.75,
                1.0,
            ],
            index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
            dtype=dtype,
        )
        tm.assert_series_equal(result, expected)

    def test_describe_one_element_ea(self):
        # GH#52515
        ser = pd.Series([0.0], dtype="Float64")
        with tm.assert_produces_warning(None):
            result = ser.describe()
        expected = pd.Series(
            [1, 0, pd.NA, 0, 0, 0, 0, 0],
            dtype="Float64",
            index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"include": [np.number]},
            {"exclude": [object]},
            {"include": "all"},
            {"include": [np.number], "exclude": [object]},
        ],
    )
    def test_describe_include_exclude_deprecated(self, kwargs):
        # GH#54193 - include/exclude have no effect on Series; deprecate them.
        ser = pd.Series([1, 2, 3], name="int_data")
        msg = "'include' and 'exclude' arguments are deprecated for Series.describe"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = ser.describe(**kwargs)
        expected = ser.describe()
        tm.assert_series_equal(result, expected)
