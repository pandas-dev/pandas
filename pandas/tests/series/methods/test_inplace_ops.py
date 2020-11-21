import pytest

from pandas import Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "dtype1, dtype2, dtype_expected, dtype_mul",
    (
        ["Int64"] * 4,
        ["float"] * 4,
        ["Int64"] + ["float"] * 3,
        pytest.param(
            "Int64",
            "Float64",
            "Float64",
            "Float64",
            marks=pytest.mark.xfail(reason="Not implemented yet"),
        ),
    ),
)
def test_series_inplace_ops(dtype1, dtype2, dtype_expected, dtype_mul):

    ser1 = Series([1], dtype=dtype1)
    ser2 = Series([2], dtype=dtype2)
    ser1 += ser2
    expected = Series([3], dtype=dtype_expected)
    tm.assert_series_equal(ser1, expected)

    ser1 -= ser2
    expected = Series([1], dtype=dtype_expected)
    tm.assert_series_equal(ser1, expected)

    ser1 *= ser2
    expected = Series([2], dtype=dtype_mul)
    tm.assert_series_equal(ser1, expected)
