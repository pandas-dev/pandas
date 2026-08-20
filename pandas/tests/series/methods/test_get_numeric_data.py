from pandas import (
    Series,
    date_range,
)
import pandas._testing as tm


class TestGetNumericData:
    def test_get_numeric_data_preserve_dtype(self):
        # get the numeric data
        obj = Series([1, 2, 3])
        result = obj._get_numeric_data()
        tm.assert_series_equal(result, obj)

        # returned object is a shallow copy
        result.iloc[0] = 0
        assert obj.iloc[0] == 1

        obj = Series([1, "2", 3.0])
        result = obj._get_numeric_data()
        expected = Series([], dtype=object)
        tm.assert_series_equal(result, expected)

        obj = Series([True, False, True])
        result = obj._get_numeric_data()
        tm.assert_series_equal(result, obj)

        obj = Series(date_range("20130101", periods=3, unit="ns"))
        result = obj._get_numeric_data()
        expected = Series([], dtype="M8[ns]")
        tm.assert_series_equal(result, expected)
