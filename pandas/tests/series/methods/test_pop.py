import numpy as np

from pandas import Series
import pandas._testing as tm


def test_pop():
    # GH#6600
    ser = Series([0, 4, 0], index=["A", "B", "C"], name=4)

    result = ser.pop("B")
    assert result == 4

    expected = Series([0, 0], index=["A", "C"], name=4)
    tm.assert_series_equal(ser, expected)


def test_pop_object_dtype_preserves_numpy_scalars():
    # GH#64266
    left, right = np.int8(1), np.int8(3)
    ser = Series([left, right], dtype=object)
    result = ser.pop(0)
    assert result is left

    expected = Series([right], index=[1], dtype=object)
    tm.assert_series_equal(ser, expected)
