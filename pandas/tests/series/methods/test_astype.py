from pandas import Series, date_range
import pandas._testing as tm


class TestAstype:
    def test_astype_dt64_to_str(self):
        # GH#10442 : testing astype(str) is correct for Series/DatetimeIndex
        dti = date_range("2012-01-01", periods=3)
        result = Series(dti).astype(str)
        expected = Series(["2012-01-01", "2012-01-02", "2012-01-03"], dtype=object)
        tm.assert_series_equal(result, expected)

    def test_astype_dt64tz_to_str(self):
        # GH#10442 : testing astype(str) is correct for Series/DatetimeIndex
        dti_tz = date_range("2012-01-01", periods=3, tz="US/Eastern")
        result = Series(dti_tz).astype(str)
        expected = Series(
            [
                "2012-01-01 00:00:00-05:00",
                "2012-01-02 00:00:00-05:00",
                "2012-01-03 00:00:00-05:00",
            ],
            dtype=object,
        )
        tm.assert_series_equal(result, expected)
