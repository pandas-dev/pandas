import numpy as np

from pandas import Series, Timestamp, date_range
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
