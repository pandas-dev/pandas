import pytest

from pandas import MultiIndex, Series
import pandas._testing as tm


class TestDropLevel:
    def test_droplevel(self):
        # GH#20342
        ser = Series([1, 2, 3, 4])
        ser.index = MultiIndex.from_arrays(
            [(1, 2, 3, 4), (5, 6, 7, 8)], names=["a", "b"]
        )
        expected = ser.reset_index("b", drop=True)
        result = ser.droplevel("b", axis="index")
        tm.assert_series_equal(result, expected)
        # test that droplevel raises ValueError on axis != 0
        with pytest.raises(ValueError):
            ser.droplevel(1, axis="columns")
