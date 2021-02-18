import pytest

from pandas import Series
import pandas._testing as tm


class TestSeriesUnaryOps:
    # __neg__, __pos__, __inv__

    def test_neg(self):
        ser = tm.makeStringSeries()
        ser.name = "series"
        tm.assert_series_equal(-ser, -1 * ser)

    def test_invert(self):
        ser = tm.makeStringSeries()
        ser.name = "series"
        tm.assert_series_equal(-(ser < 0), ~(ser < 0))

    @pytest.mark.parametrize(
        "source, target",
        [
            ([1, 2, 3], [-1, -2, -3]),
            ([1, 2, None], [-1, -2, None]),
            ([-1, 0, 1], [1, 0, -1]),
        ],
    )
    def test_unary_minus_nullable_int(
        self, any_signed_nullable_int_dtype, source, target
    ):
        dtype = any_signed_nullable_int_dtype
        ser = Series(source, dtype=dtype)
        result = -ser
        expected = Series(target, dtype=dtype)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("source", [[1, 2, 3], [1, 2, None], [-1, 0, 1]])
    def test_unary_plus_nullable_int(self, any_signed_nullable_int_dtype, source):
        dtype = any_signed_nullable_int_dtype
        expected = Series(source, dtype=dtype)
        result = +expected
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "source, target",
        [
            ([1, 2, 3], [1, 2, 3]),
            ([1, -2, None], [1, 2, None]),
            ([-1, 0, 1], [1, 0, 1]),
        ],
    )
    def test_abs_nullable_int(self, any_signed_nullable_int_dtype, source, target):
        dtype = any_signed_nullable_int_dtype
        ser = Series(source, dtype=dtype)
        result = abs(ser)
        expected = Series(target, dtype=dtype)
        tm.assert_series_equal(result, expected)
