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
    # GH#64266 future.python_scalars should not unbox numpy scalars
    #  stored in object dtype
    ser = Series([np.int8(1), np.int8(3)], dtype=object)
    result = ser.pop(0)
    assert isinstance(result, np.int8)
