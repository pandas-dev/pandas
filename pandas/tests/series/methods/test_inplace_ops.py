import pytest

from pandas import Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "ser1, ser2, expected_add, expected_sub, expected_mul",
    (
        [
            Series([1], dtype="Int64"),
            Series([2], dtype="Int64"),
            Series([3], dtype="Int64"),
            Series([-1], dtype="Int64"),
            Series([2], dtype="Int64"),
        ],
        [
            Series([1], dtype="float"),
            Series([2.0], dtype="float"),
            Series([3.0], dtype="float"),
            Series([-1.0], dtype="float"),
            Series([2.0], dtype="float"),
        ],
        [
            Series([1], dtype="Int64"),
            Series([2.0], dtype="float"),
            Series([3.0], dtype="float"),
            Series([-1.0], dtype="float"),
            Series([2], dtype="float"),
        ],
        [
            Series([1.0], dtype="float"),
            Series([2], dtype="Int64"),
            Series([3.0], dtype="float"),
            Series([-1.0], dtype="float"),
            Series([2], dtype="float"),
        ],
        pytest.param(
            Series([1], dtype="Int64"),
            Series([2.0], dtype="Float64"),
            Series([3.0], dtype="Float64"),
            Series([-1.0], dtype="Float64"),
            Series([2], dtype="Float64"),
            marks=pytest.mark.xfail(reason="Not implemented yet"),
        ),
    ),
)
def test_series_inplace_ops(ser1, ser2, expected_add, expected_sub, expected_mul):

    res = ser1.copy()
    res += ser2
    tm.assert_series_equal(res, expected_add)

    res = ser1.copy()
    res -= ser2
    tm.assert_series_equal(res, expected_sub)

    res = ser1.copy()
    res *= ser2
    tm.assert_series_equal(res, expected_mul)
