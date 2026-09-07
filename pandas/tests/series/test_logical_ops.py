from datetime import datetime
import operator

import numpy as np
import pytest

from pandas.compat import HAS_PYARROW

import pandas as pd
import pandas._testing as tm
from pandas.core import ops


class TestSeriesLogicalOps:
    @pytest.mark.parametrize("bool_op", [operator.and_, operator.or_, operator.xor])
    def test_bool_operators_with_nas(self, bool_op):
        # boolean &, |, ^ should work with object arrays and propagate NAs
        ser = pd.Series(pd.bdate_range("1/1/2000", periods=10), dtype=object)
        ser[::2] = np.nan

        mask = ser.isna()
        filled = ser.fillna(ser[0])

        result = bool_op(ser < ser[9], ser > ser[3])

        expected = bool_op(filled < filled[9], filled > filled[3])
        expected[mask] = False
        tm.assert_series_equal(result, expected)

    def test_logical_operators_bool_dtype_with_empty(self):
        # GH#9016: support bitwise op for integer types
        index = list("bca")

        s_tft = pd.Series([True, False, True], index=index)
        s_fff = pd.Series([False, False, False], index=index)
        s_empty = pd.Series([], dtype=object)

        res = s_tft & s_empty
        expected = s_fff.sort_index()
        tm.assert_series_equal(res, expected)

        res = s_tft | s_empty
        expected = s_tft.sort_index()
        tm.assert_series_equal(res, expected)

    def test_logical_operators_int_dtype_with_int_dtype(self):
        # GH#9016: support bitwise op for integer types

        s_0123 = pd.Series(range(4), dtype="int64")
        s_3333 = pd.Series([3] * 4)
        s_4444 = pd.Series([4] * 4)

        res = s_0123 & s_3333
        expected = pd.Series(range(4), dtype="int64")
        tm.assert_series_equal(res, expected)

        res = s_0123 | s_4444
        expected = pd.Series(range(4, 8), dtype="int64")
        tm.assert_series_equal(res, expected)

        s_1111 = pd.Series([1] * 4, dtype="int8")
        res = s_0123 & s_1111
        expected = pd.Series([0, 1, 0, 1], dtype="int64")
        tm.assert_series_equal(res, expected)

        res = s_0123.astype(np.int16) | s_1111.astype(np.int32)
        expected = pd.Series([1, 1, 3, 3], dtype="int32")
        tm.assert_series_equal(res, expected)

    def test_logical_operators_int_dtype_with_int_scalar(self):
        # GH#9016: support bitwise op for integer types
        s_0123 = pd.Series(range(4), dtype="int64")

        res = s_0123 & 0
        expected = pd.Series([0] * 4)
        tm.assert_series_equal(res, expected)

        res = s_0123 & 1
        expected = pd.Series([0, 1, 0, 1])
        tm.assert_series_equal(res, expected)

    def test_logical_operators_int_dtype_with_float(self):
        # GH#9016: support bitwise op for integer types
        s_0123 = pd.Series(range(4), dtype="int64")

        err_msg = (
            r"Logical ops \(and, or, xor\) between Pandas objects and "
            "dtype-less sequences"
        )

        msg = "Cannot perform.+with a dtyped.+array and scalar of type"
        with pytest.raises(TypeError, match=msg):
            s_0123 & np.nan
        with pytest.raises(TypeError, match=msg):
            s_0123 & 3.14
        msg = "unsupported operand type.+for &:"
        with pytest.raises(TypeError, match=err_msg):
            s_0123 & [0.1, 4, 3.14, 2]
        with pytest.raises(TypeError, match=msg):
            s_0123 & np.array([0.1, 4, 3.14, 2])
        with pytest.raises(TypeError, match=msg):
            s_0123 & pd.Series([0.1, 4, -3.14, 2])

    def test_logical_operators_int_dtype_with_str(self):
        s_1111 = pd.Series([1] * 4, dtype="int8")

        err_msg = (
            r"Logical ops \(and, or, xor\) between Pandas objects and "
            "dtype-less sequences"
        )

        msg = "Cannot perform 'and_' with a dtyped.+array and scalar of type"
        with pytest.raises(TypeError, match=msg):
            s_1111 & "a"
        with pytest.raises(TypeError, match=err_msg):
            s_1111 & ["a", "b", "c", "d"]

    def test_logical_operators_int_dtype_with_bool(self):
        # GH#9016: support bitwise op for integer types
        s_0123 = pd.Series(range(4), dtype="int64")

        expected = pd.Series([False] * 4)

        result = s_0123 & False
        tm.assert_series_equal(result, expected)

        msg = (
            r"Logical ops \(and, or, xor\) between Pandas objects and "
            "dtype-less sequences"
        )
        with pytest.raises(TypeError, match=msg):
            s_0123 & [False]

        with pytest.raises(TypeError, match=msg):
            s_0123 & (False,)

        result = s_0123 ^ False
        expected = pd.Series([False, True, True, True])
        tm.assert_series_equal(result, expected)

    def test_logical_operators_int_dtype_with_object(self, using_infer_string):
        # GH#9016: support bitwise op for integer types
        s_0123 = pd.Series(range(4), dtype="int64")

        result = s_0123 & pd.Series([False, np.nan, False, False])
        expected = pd.Series([False] * 4)
        tm.assert_series_equal(result, expected)

        s_abNd = pd.Series(["a", "b", np.nan, "d"])
        # pyarrow-backed str routes through the pandas op; object dtype and the
        # python-backed str fallback hit the Python operator instead
        if using_infer_string and HAS_PYARROW:
            msg = "'rand_' not supported"
        else:
            msg = r"unsupported operand type\(s\) for &: 'int' and 'str'"
        with pytest.raises(TypeError, match=msg):
            s_0123 & s_abNd

    def test_logical_operators_bool_dtype_with_int(self):
        index = list("bca")

        s_tft = pd.Series([True, False, True], index=index)
        s_fff = pd.Series([False, False, False], index=index)

        res = s_tft & 0
        expected = s_fff
        tm.assert_series_equal(res, expected)

        res = s_tft & 1
        expected = s_tft
        tm.assert_series_equal(res, expected)

    def test_logical_ops_bool_dtype_with_ndarray(self):
        # make sure we operate on ndarray the same as Series
        left = pd.Series([True, True, True, False, True])
        right = [True, False, None, True, np.nan]

        msg = (
            r"Logical ops \(and, or, xor\) between Pandas objects and "
            "dtype-less sequences"
        )

        expected = pd.Series([True, False, False, False, False])
        with pytest.raises(TypeError, match=msg):
            left & right
        result = left & np.array(right)
        tm.assert_series_equal(result, expected)
        result = left & pd.Index(right)
        tm.assert_series_equal(result, expected)
        result = left & pd.Series(right)
        tm.assert_series_equal(result, expected)

        expected = pd.Series([True, True, True, True, True])
        with pytest.raises(TypeError, match=msg):
            left | right
        result = left | np.array(right)
        tm.assert_series_equal(result, expected)
        result = left | pd.Index(right)
        tm.assert_series_equal(result, expected)
        result = left | pd.Series(right)
        tm.assert_series_equal(result, expected)

        expected = pd.Series([False, True, True, True, True])
        with pytest.raises(TypeError, match=msg):
            left ^ right
        result = left ^ np.array(right)
        tm.assert_series_equal(result, expected)
        result = left ^ pd.Index(right)
        tm.assert_series_equal(result, expected)
        result = left ^ pd.Series(right)
        tm.assert_series_equal(result, expected)

    def test_logical_operators_int_dtype_with_bool_dtype_and_reindex(self):
        # GH#9016: support bitwise op for integer types

        index = list("bca")

        s_tft = pd.Series([True, False, True], index=index)
        s_tft = pd.Series([True, False, True], index=index)
        s_tff = pd.Series([True, False, False], index=index)

        s_0123 = pd.Series(range(4), dtype="int64")

        # s_0123 will be all false now because of reindexing like s_tft
        expected = pd.Series([False] * 7, index=[0, 1, 2, 3, "a", "b", "c"])
        result = s_tft & s_0123
        tm.assert_series_equal(result, expected)

        # GH#52538: no longer to object type when reindex is needed;
        # matches DataFrame behavior
        msg = r"unsupported operand type\(s\) for &: 'float' and 'bool'"
        with pytest.raises(TypeError, match=msg):
            s_0123 & s_tft

        s_a0b1c0 = pd.Series([1], list("b"))

        res = s_tft & s_a0b1c0
        expected = s_tff.reindex(list("abc"))
        tm.assert_series_equal(res, expected)

        res = s_tft | s_a0b1c0
        expected = s_tft.reindex(list("abc"))
        tm.assert_series_equal(res, expected)

    def test_scalar_na_logical_ops_corners(self):
        s = pd.Series([2, 3, 4, 5, 6, 7, 8, 9, 10])

        msg = "Cannot perform.+with a dtyped.+array and scalar of type"
        with pytest.raises(TypeError, match=msg):
            s & datetime(2005, 1, 1)

        s = pd.Series([2, 3, 4, 5, 6, 7, 8, 9, datetime(2005, 1, 1)])
        s[::2] = np.nan

        expected = pd.Series(True, index=s.index)
        expected[::2] = False

        msg = (
            r"Logical ops \(and, or, xor\) between Pandas objects and "
            "dtype-less sequences"
        )
        with pytest.raises(TypeError, match=msg):
            s & list(s)

    def test_scalar_na_logical_ops_corners_aligns(self):
        s = pd.Series([2, 3, 4, 5, 6, 7, 8, 9, datetime(2005, 1, 1)])
        s[::2] = np.nan
        d = pd.DataFrame({"A": s})

        expected = pd.DataFrame(False, index=range(9), columns=["A", *list(range(9))])

        result = s & d
        tm.assert_frame_equal(result, expected)

        result = d & s
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("op", [operator.and_, operator.or_, operator.xor])
    def test_logical_ops_with_index(self, op):
        # GH#22092, GH#19792
        ser = pd.Series([True, True, False, False])
        idx1 = pd.Index([True, False, True, False])
        idx2 = pd.Index([1, 0, 1, 0])

        expected = pd.Series([op(ser[n], idx1[n]) for n in range(len(ser))])

        result = op(ser, idx1)
        tm.assert_series_equal(result, expected)

        expected = pd.Series([op(ser[n], idx2[n]) for n in range(len(ser))], dtype=bool)

        result = op(ser, idx2)
        tm.assert_series_equal(result, expected)

    def test_reversed_xor_with_index_returns_series(self):
        # GH#22092, GH#19792 pre-2.0 these were aliased to setops
        ser = pd.Series([True, True, False, False])
        idx1 = pd.Index([True, False, True, False], dtype=bool)
        idx2 = pd.Index([1, 0, 1, 0])

        expected = pd.Series([False, True, True, False])
        result = idx1 ^ ser
        tm.assert_series_equal(result, expected)

        result = idx2 ^ ser
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "op",
        [
            ops.rand_,
            ops.ror_,
        ],
    )
    def test_reversed_logical_op_with_index_returns_series(self, op):
        # GH#22092, GH#19792
        ser = pd.Series([True, True, False, False])
        idx1 = pd.Index([True, False, True, False])
        idx2 = pd.Index([1, 0, 1, 0])

        expected = pd.Series(op(idx1.values, ser.values))
        result = op(ser, idx1)
        tm.assert_series_equal(result, expected)

        expected = op(ser, pd.Series(idx2))
        result = op(ser, idx2)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "op, expected",
        [
            (ops.rand_, [False, False]),
            (ops.ror_, [True, True]),
            (ops.rxor, [True, True]),
        ],
    )
    def test_reverse_ops_with_index(self, op, expected):
        # https://github.com/pandas-dev/pandas/pull/23628
        # multi-set Index ops are buggy, so let's avoid duplicates...
        # GH#49503
        ser = pd.Series([True, False])
        idx = pd.Index([False, True])

        result = op(ser, idx)
        expected = pd.Series(expected)
        tm.assert_series_equal(result, expected)

    def test_logical_ops_label_based(self, using_infer_string):
        # GH#4947
        # logical ops should be label based

        a = pd.Series([True, False, True], list("bca"))
        b = pd.Series([False, True, False], list("abc"))

        expected = pd.Series([False, True, False], list("abc"))
        result = a & b
        tm.assert_series_equal(result, expected)

        expected = pd.Series([True, True, False], list("abc"))
        result = a | b
        tm.assert_series_equal(result, expected)

        expected = pd.Series([True, False, False], list("abc"))
        result = a ^ b
        tm.assert_series_equal(result, expected)

        # rhs is bigger
        a = pd.Series([True, False, True], list("bca"))
        b = pd.Series([False, True, False, True], list("abcd"))

        expected = pd.Series([False, True, False, False], list("abcd"))
        result = a & b
        tm.assert_series_equal(result, expected)

        expected = pd.Series([True, True, False, False], list("abcd"))
        result = a | b
        tm.assert_series_equal(result, expected)

        # filling

        # vs empty
        empty = pd.Series([], dtype=object)

        result = a & empty
        expected = pd.Series([False, False, False], list("abc"))
        tm.assert_series_equal(result, expected)

        result = a | empty
        expected = pd.Series([True, True, False], list("abc"))
        tm.assert_series_equal(result, expected)

        # vs non-matching
        result = a & pd.Series([1], ["z"])
        expected = pd.Series([False, False, False, False], list("abcz"))
        tm.assert_series_equal(result, expected)

        result = a | pd.Series([1], ["z"])
        expected = pd.Series([True, True, False, False], list("abcz"))
        tm.assert_series_equal(result, expected)

        # identity
        # we would like s[s|e] == s to hold for any e, whether empty or not
        for e in [
            empty.copy(),
            pd.Series([1], ["z"]),
            pd.Series(np.nan, b.index),
            pd.Series(np.nan, a.index),
        ]:
            result = a[a | e]
            tm.assert_series_equal(result, a[a])

        for e in [pd.Series(["z"])]:
            if using_infer_string:
                # TODO(infer_string) should this behave differently?
                # -> https://github.com/pandas-dev/pandas/issues/60234
                with pytest.raises(
                    TypeError,
                    match="|".join(
                        ["not supported for dtype", "unsupported operand type"]
                    ),
                ):
                    result = a[a | e]
            else:
                result = a[a | e]
            tm.assert_series_equal(result, a[a])

        # vs scalars
        index = list("bca")
        t = pd.Series([True, False, True])

        for v in [True, 1, 2]:
            result = pd.Series([True, False, True], index=index) | v
            expected = pd.Series([True, True, True], index=index)
            tm.assert_series_equal(result, expected)

        msg = "Cannot perform.+with a dtyped.+array and scalar of type"
        for v in [np.nan, "foo"]:
            with pytest.raises(TypeError, match=msg):
                t | v

        for v in [False, 0]:
            result = pd.Series([True, False, True], index=index) | v
            expected = pd.Series([True, False, True], index=index)
            tm.assert_series_equal(result, expected)

        for v in [True, 1]:
            result = pd.Series([True, False, True], index=index) & v
            expected = pd.Series([True, False, True], index=index)
            tm.assert_series_equal(result, expected)

        for v in [False, 0]:
            result = pd.Series([True, False, True], index=index) & v
            expected = pd.Series([False, False, False], index=index)
            tm.assert_series_equal(result, expected)
        msg = "Cannot perform.+with a dtyped.+array and scalar of type"
        for v in [np.nan]:
            with pytest.raises(TypeError, match=msg):
                t & v

    def test_logical_ops_df_compat(self):
        # GH#1134
        s1 = pd.Series([True, False, True], index=list("ABC"), name="x")
        s2 = pd.Series([True, True, False], index=list("ABD"), name="x")

        exp = pd.Series([True, False, False, False], index=list("ABCD"), name="x")
        tm.assert_series_equal(s1 & s2, exp)
        tm.assert_series_equal(s2 & s1, exp)

        # True | np.nan => True
        exp_or1 = pd.Series([True, True, True, False], index=list("ABCD"), name="x")
        tm.assert_series_equal(s1 | s2, exp_or1)
        # np.nan | True => np.nan, filled with False
        exp_or = pd.Series([True, True, False, False], index=list("ABCD"), name="x")
        tm.assert_series_equal(s2 | s1, exp_or)

        # DataFrame doesn't fill nan with False
        tm.assert_frame_equal(s1.to_frame() & s2.to_frame(), exp.to_frame())
        tm.assert_frame_equal(s2.to_frame() & s1.to_frame(), exp.to_frame())

        exp = pd.DataFrame({"x": [True, True, np.nan, np.nan]}, index=list("ABCD"))
        tm.assert_frame_equal(s1.to_frame() | s2.to_frame(), exp_or1.to_frame())
        tm.assert_frame_equal(s2.to_frame() | s1.to_frame(), exp_or.to_frame())

        # different length
        s3 = pd.Series([True, False, True], index=list("ABC"), name="x")
        s4 = pd.Series([True, True, True, True], index=list("ABCD"), name="x")

        exp = pd.Series([True, False, True, False], index=list("ABCD"), name="x")
        tm.assert_series_equal(s3 & s4, exp)
        tm.assert_series_equal(s4 & s3, exp)

        # np.nan | True => np.nan, filled with False
        exp_or1 = pd.Series([True, True, True, False], index=list("ABCD"), name="x")
        tm.assert_series_equal(s3 | s4, exp_or1)
        # True | np.nan => True
        exp_or = pd.Series([True, True, True, True], index=list("ABCD"), name="x")
        tm.assert_series_equal(s4 | s3, exp_or)

        tm.assert_frame_equal(s3.to_frame() & s4.to_frame(), exp.to_frame())
        tm.assert_frame_equal(s4.to_frame() & s3.to_frame(), exp.to_frame())

        tm.assert_frame_equal(s3.to_frame() | s4.to_frame(), exp_or1.to_frame())
        tm.assert_frame_equal(s4.to_frame() | s3.to_frame(), exp_or.to_frame())

    def test_int_dtype_different_index_not_bool(self):
        # GH 52500
        ser1 = pd.Series([1, 2, 3], index=[10, 11, 23], name="a")
        ser2 = pd.Series([10, 20, 30], index=[11, 10, 23], name="a")
        result = np.bitwise_xor(ser1, ser2)
        expected = pd.Series([21, 8, 29], index=[10, 11, 23], name="a")
        tm.assert_series_equal(result, expected)

        result = ser1 ^ ser2
        tm.assert_series_equal(result, expected)
