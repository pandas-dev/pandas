from pandas import Series, date_range
import pandas._testing as tm


class TestView:
    def test_view_tz(self):
        # GH#24024
        ser = Series(date_range("2000", periods=4, tz="US/Central"))
        result = ser.view("i8")
        expected = Series(
            [
                946706400000000000,
                946792800000000000,
                946879200000000000,
                946965600000000000,
            ]
        )
        tm.assert_series_equal(result, expected)
