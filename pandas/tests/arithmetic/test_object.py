# Arithmetic tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for object dtype
from decimal import Decimal
import operator

import numpy as np
import pytest

import pandas as pd
from pandas import Series, Timestamp
import pandas._testing as tm
from pandas.core import ops

# ------------------------------------------------------------------
# Comparisons


class TestObjectComparisons:
    def test_comparison_object_numeric_nas(self):
        ser = Series(np.random.randn(10), dtype=object)
        shifted = ser.shift(2)

        ops = ["lt", "le", "gt", "ge", "eq", "ne"]
        for op in ops:
            func = getattr(operator, op)

            result = func(ser, shifted)
            expected = func(ser.astype(float), shifted.astype(float))
            tm.assert_series_equal(result, expected)

    def test_object_comparisons(self):
        ser = Series(["a", "b", np.nan, "c", "a"])

        result = ser == "a"
        expected = Series([True, False, False, False, True])
        tm.assert_series_equal(result, expected)

        result = ser < "a"
        expected = Series([False, False, False, False, False])
        tm.assert_series_equal(result, expected)

        result = ser != "a"
        expected = -(ser == "a")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [None, object])
    def test_more_na_comparisons(self, dtype):
        left = Series(["a", np.nan, "c"], dtype=dtype)
        right = Series(["a", np.nan, "d"], dtype=dtype)

        result = left == right
        expected = Series([True, False, False])
        tm.assert_series_equal(result, expected)

        result = left != right
        expected = Series([False, True, True])
        tm.assert_series_equal(result, expected)

        result = left == np.nan
        expected = Series([False, False, False])
        tm.assert_series_equal(result, expected)

        result = left != np.nan
        expected = Series([True, True, True])
        tm.assert_series_equal(result, expected)


# ------------------------------------------------------------------
# Arithmetic


class TestArithmetic:

    # TODO: parametrize
    def test_pow_ops_object(self):
        # GH#22922
        # pow is weird with masking & 1, so testing here
        a = Series([1, np.nan, 1, np.nan], dtype=object)
        b = Series([1, np.nan, np.nan, 1], dtype=object)
        result = a ** b
        expected = Series(a.values ** b.values, dtype=object)
        tm.assert_series_equal(result, expected)

        result = b ** a
        expected = Series(b.values ** a.values, dtype=object)

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    @pytest.mark.parametrize("other", ["category", "Int64"])
    def test_add_extension_scalar(self, other, box_with_array, op):
        # GH#22378
        # Check that scalars satisfying is_extension_array_dtype(obj)
        # do not incorrectly try to dispatch to an ExtensionArray operation

        arr = pd.Series(["a", "b", "c"])
        expected = pd.Series([op(x, other) for x in arr])

        arr = tm.box_expected(arr, box_with_array)
        expected = tm.box_expected(expected, box_with_array)

        result = op(arr, other)
        tm.assert_equal(result, expected)

    def test_objarr_add_str(self, box):
        ser = pd.Series(["x", np.nan, "x"])
        expected = pd.Series(["xa", np.nan, "xa"])

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = ser + "a"
        tm.assert_equal(result, expected)

    def test_objarr_radd_str(self, box):
        ser = pd.Series(["x", np.nan, "x"])
        expected = pd.Series(["ax", np.nan, "ax"])

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = "a" + ser
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize(
        "data",
        [
            [1, 2, 3],
            [1.1, 2.2, 3.3],
            [Timestamp("2011-01-01"), Timestamp("2011-01-02"), pd.NaT],
            ["x", "y", 1],
        ],
    )
    @pytest.mark.parametrize("dtype", [None, object])
    def test_objarr_radd_str_invalid(self, dtype, data, box_with_array):
        ser = Series(data, dtype=dtype)

        ser = tm.box_expected(ser, box_with_array)
        msg = (
            "can only concatenate str|"
            "did not contain a loop with signature matching types|"
            "unsupported operand type|"
            "must be str"
        )
        with pytest.raises(TypeError, match=msg):
            "foo_" + ser

    @pytest.mark.parametrize("op", [operator.add, ops.radd, operator.sub, ops.rsub])
    def test_objarr_add_invalid(self, op, box_with_array):
        # invalid ops
        box = box_with_array

        obj_ser = tm.makeObjectSeries()
        obj_ser.name = "objects"

        obj_ser = tm.box_expected(obj_ser, box)
        with pytest.raises(Exception):
            op(obj_ser, 1)
        with pytest.raises(Exception):
            op(obj_ser, np.array(1, dtype=np.int64))

    # TODO: Moved from tests.series.test_operators; needs cleanup
    def test_operators_na_handling(self):
        ser = Series(["foo", "bar", "baz", np.nan])
        result = "prefix_" + ser
        expected = pd.Series(["prefix_foo", "prefix_bar", "prefix_baz", np.nan])
        tm.assert_series_equal(result, expected)

        result = ser + "_suffix"
        expected = pd.Series(["foo_suffix", "bar_suffix", "baz_suffix", np.nan])
        tm.assert_series_equal(result, expected)

    # TODO: parametrize over box
    @pytest.mark.parametrize("dtype", [None, object])
    def test_series_with_dtype_radd_timedelta(self, dtype):
        # note this test is _not_ aimed at timedelta64-dtyped Series
        ser = pd.Series(
            [pd.Timedelta("1 days"), pd.Timedelta("2 days"), pd.Timedelta("3 days")],
            dtype=dtype,
        )
        expected = pd.Series(
            [pd.Timedelta("4 days"), pd.Timedelta("5 days"), pd.Timedelta("6 days")]
        )

        result = pd.Timedelta("3 days") + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.Timedelta("3 days")
        tm.assert_series_equal(result, expected)

    # TODO: cleanup & parametrize over box
    def test_mixed_timezone_series_ops_object(self):
        # GH#13043
        ser = pd.Series(
            [
                pd.Timestamp("2015-01-01", tz="US/Eastern"),
                pd.Timestamp("2015-01-01", tz="Asia/Tokyo"),
            ],
            name="xxx",
        )
        assert ser.dtype == object

        exp = pd.Series(
            [
                pd.Timestamp("2015-01-02", tz="US/Eastern"),
                pd.Timestamp("2015-01-02", tz="Asia/Tokyo"),
            ],
            name="xxx",
        )
        tm.assert_series_equal(ser + pd.Timedelta("1 days"), exp)
        tm.assert_series_equal(pd.Timedelta("1 days") + ser, exp)

        # object series & object series
        ser2 = pd.Series(
            [
                pd.Timestamp("2015-01-03", tz="US/Eastern"),
                pd.Timestamp("2015-01-05", tz="Asia/Tokyo"),
            ],
            name="xxx",
        )
        assert ser2.dtype == object
        exp = pd.Series([pd.Timedelta("2 days"), pd.Timedelta("4 days")], name="xxx")
        tm.assert_series_equal(ser2 - ser, exp)
        tm.assert_series_equal(ser - ser2, -exp)

        ser = pd.Series(
            [pd.Timedelta("01:00:00"), pd.Timedelta("02:00:00")],
            name="xxx",
            dtype=object,
        )
        assert ser.dtype == object

        exp = pd.Series(
            [pd.Timedelta("01:30:00"), pd.Timedelta("02:30:00")], name="xxx"
        )
        tm.assert_series_equal(ser + pd.Timedelta("00:30:00"), exp)
        tm.assert_series_equal(pd.Timedelta("00:30:00") + ser, exp)

    # TODO: cleanup & parametrize over box
    def test_iadd_preserves_name(self):
        # GH#17067, GH#19723 __iadd__ and __isub__ should preserve index name
        ser = pd.Series([1, 2, 3])
        ser.index.name = "foo"

        ser.index += 1
        assert ser.index.name == "foo"

        ser.index -= 1
        assert ser.index.name == "foo"

    def test_add_string(self):
        # from bug report
        index = pd.Index(["a", "b", "c"])
        index2 = index + "foo"

        assert "a" not in index2
        assert "afoo" in index2

    def test_iadd_string(self):
        index = pd.Index(["a", "b", "c"])
        # doesn't fail test unless there is a check before `+=`
        assert "a" in index

        index += "_x"
        assert "a_x" in index

    def test_add(self):
        index = tm.makeStringIndex(100)
        expected = pd.Index(index.values * 2)
        tm.assert_index_equal(index + index, expected)
        tm.assert_index_equal(index + index.tolist(), expected)
        tm.assert_index_equal(index.tolist() + index, expected)

        # test add and radd
        index = pd.Index(list("abc"))
        expected = pd.Index(["a1", "b1", "c1"])
        tm.assert_index_equal(index + "1", expected)
        expected = pd.Index(["1a", "1b", "1c"])
        tm.assert_index_equal("1" + index, expected)

    def test_sub_fail(self):
        index = tm.makeStringIndex(100)
        with pytest.raises(TypeError):
            index - "a"
        with pytest.raises(TypeError):
            index - index
        with pytest.raises(TypeError):
            index - index.tolist()
        with pytest.raises(TypeError):
            index.tolist() - index

    def test_sub_object(self):
        # GH#19369
        index = pd.Index([Decimal(1), Decimal(2)])
        expected = pd.Index([Decimal(0), Decimal(1)])

        result = index - Decimal(1)
        tm.assert_index_equal(result, expected)

        result = index - pd.Index([Decimal(1), Decimal(1)])
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            index - "foo"

        with pytest.raises(TypeError):
            index - np.array([2, "foo"])

    def test_rsub_object(self):
        # GH#19369
        index = pd.Index([Decimal(1), Decimal(2)])
        expected = pd.Index([Decimal(1), Decimal(0)])

        result = Decimal(2) - index
        tm.assert_index_equal(result, expected)

        result = np.array([Decimal(2), Decimal(2)]) - index
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            "foo" - index

        with pytest.raises(TypeError):
            np.array([True, pd.Timestamp.now()]) - index
