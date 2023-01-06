import numpy as np
import pytest

from pandas import (
    NA,
    Series,
)
import pandas._testing as tm


@pytest.mark.parametrize("dtype", ["int64", "float64"])
def test_to_numpy_na_value(dtype):
    # GH#48951
    ser = Series([1, 2, NA, 4])
    result = ser.to_numpy(dtype=dtype, na_value=0)
    expected = np.array([1, 2, 0, 4], dtype=dtype)
    tm.assert_numpy_array_equal(result, expected)
