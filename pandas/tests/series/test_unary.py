import pytest

from pandas.errors import Pandas4Warning

from pandas import Series
import pandas._testing as tm


class TestSeriesUnaryOps:
    # __neg__, __pos__, __invert__

    def test_neg(self):
        ser = Series(range(5), dtype="float64", name="series")
        tm.assert_series_equal(-ser, -1 * ser)

    def test_invert(self):
        ser = Series(range(5), dtype="float64", name="series")
        tm.assert_series_equal(-(ser < 0), ~(ser < 0))

    def test_invert_object_deprecated(self):
        # GH#51567 ~ on object dtype dispatches to the objects' own __invert__
        ser = Series([1, 2], dtype=object)

        msg = "__invert__ .* on object dtype is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = ~ser
        tm.assert_series_equal(result, Series([-2, -3], dtype=object))

    def test_invert_object_bool_mask_deprecated(self):
        # GH#16873, GH#31035 a boolean mask cast to object silently becomes
        #  integers instead of being negated
        ser = Series([True, False], dtype=object)

        msg = "__invert__ .* on object dtype is deprecated"
        # Python itself deprecated ~ on bool in 3.12, so on new enough versions
        #  there is a second warning here that we don't care about.
        with tm.assert_produces_warning(
            Pandas4Warning, match=msg, raise_on_extra_warnings=False
        ):
            result = ~ser
        tm.assert_series_equal(result, Series([-2, -1], dtype=object))

    def test_invert_object_raising_still_warns(self):
        # GH#51567 elements that have no __invert__ warn before raising
        ser = Series([1.0, 2.0], dtype=object)

        msg = "__invert__ .* on object dtype is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            with pytest.raises(TypeError, match="bad operand type for unary ~"):
                ~ser

    @pytest.mark.parametrize(
        "source, neg_target, abs_target",
        [
            ([1, 2, 3], [-1, -2, -3], [1, 2, 3]),
            ([1, 2, None], [-1, -2, None], [1, 2, None]),
        ],
    )
    def test_all_numeric_unary_operators(
        self, any_numeric_ea_dtype, source, neg_target, abs_target
    ):
        # GH38794
        dtype = any_numeric_ea_dtype
        ser = Series(source, dtype=dtype)
        neg_result, pos_result, abs_result = -ser, +ser, abs(ser)
        if dtype.startswith("U"):
            neg_target = -Series(source, dtype=dtype)
        else:
            neg_target = Series(neg_target, dtype=dtype)

        abs_target = Series(abs_target, dtype=dtype)

        tm.assert_series_equal(neg_result, neg_target)
        tm.assert_series_equal(pos_result, ser)
        tm.assert_series_equal(abs_result, abs_target)

    @pytest.mark.parametrize("op", ["__neg__", "__abs__"])
    def test_unary_float_op_mask(self, float_ea_dtype, op):
        dtype = float_ea_dtype
        ser = Series([1.1, 2.2, 3.3], dtype=dtype)
        result = getattr(ser, op)()
        target = result.copy(deep=True)
        ser[0] = None
        tm.assert_series_equal(result, target)
