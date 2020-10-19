import numpy as np

from pandas import Period, Series, Timedelta, Timestamp, date_range
import pandas._testing as tm


class TestSeriesDescribe:
    def test_describe(self):
        s = Series([0, 1, 2, 3, 4], name="int_data")
        result = s.describe()
        expected = Series(
            [5, 2, s.std(), 0, 1, 2, 3, 4],
            name="int_data",
            index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

        s = Series([True, True, False, False, False], name="bool_data")
        result = s.describe()
        expected = Series(
            [5, 2, False, 3], name="bool_data", index=["count", "unique", "top", "freq"]
        )
        tm.assert_series_equal(result, expected)

        s = Series(["a", "a", "b", "c", "d"], name="str_data")
        result = s.describe()
        expected = Series(
            [5, 4, "a", 2], name="str_data", index=["count", "unique", "top", "freq"]
        )
        tm.assert_series_equal(result, expected)

        s = Series(
            [
                Timedelta("1 days"),
                Timedelta("2 days"),
                Timedelta("3 days"),
                Timedelta("4 days"),
                Timedelta("5 days"),
            ],
            name="timedelta_data",
        )
        result = s.describe()
        expected = Series(
            [5, s[2], s.std(), s[0], s[1], s[2], s[3], s[4]],
            name="timedelta_data",
            index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)

        s = Series(
            [Period("2020-01", "M"), Period("2020-01", "M"), Period("2019-12", "M")],
            name="period_data",
        )
        result = s.describe()
        expected = Series(
            [3, 2, s[0], 2],
            name="period_data",
            index=["count", "unique", "top", "freq"],
        )
        tm.assert_series_equal(result, expected)

    def test_describe_empty_object(self):
        # https://github.com/pandas-dev/pandas/issues/27183
        s = Series([None, None], dtype=object)
        result = s.describe()
        expected = Series(
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
        start = Timestamp(2018, 1, 1)
        end = Timestamp(2018, 1, 5)
        s = Series(date_range(start, end, tz=tz), name=name)
        result = s.describe(datetime_is_numeric=True)
        expected = Series(
            [
                5,
                Timestamp(2018, 1, 3).tz_localize(tz),
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

    def test_describe_with_tz_warns(self):
        name = tz = "CET"
        start = Timestamp(2018, 1, 1)
        end = Timestamp(2018, 1, 5)
        s = Series(date_range(start, end, tz=tz), name=name)

        with tm.assert_produces_warning(FutureWarning):
            result = s.describe()

        expected = Series(
            [
                5,
                5,
                s.value_counts().index[0],
                1,
                start.tz_localize(tz),
                end.tz_localize(tz),
            ],
            name=name,
            index=["count", "unique", "top", "freq", "first", "last"],
        )
        tm.assert_series_equal(result, expected)

    def test_datetime_is_numeric_includes_datetime(self):
        s = Series(date_range("2012", periods=3))
        result = s.describe(datetime_is_numeric=True)
        expected = Series(
            [
                3,
                Timestamp("2012-01-02"),
                Timestamp("2012-01-01"),
                Timestamp("2012-01-01T12:00:00"),
                Timestamp("2012-01-02"),
                Timestamp("2012-01-02T12:00:00"),
                Timestamp("2012-01-03"),
            ],
            index=["count", "mean", "min", "25%", "50%", "75%", "max"],
        )
        tm.assert_series_equal(result, expected)
