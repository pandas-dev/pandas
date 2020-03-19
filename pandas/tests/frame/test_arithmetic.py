from collections import deque
from datetime import datetime
import operator
import re

import numpy as np
import pytest
import pytz

import pandas as pd
import pandas._testing as tm
from pandas.tests.frame.common import _check_mixed_float, _check_mixed_int

# -------------------------------------------------------------------
# Comparisons


class TestFrameComparisons:
    # Specifically _not_ flex-comparisons

    def test_frame_in_list(self):
        # GH#12689 this should raise at the DataFrame level, not blocks
        df = pd.DataFrame(np.random.randn(6, 4), columns=list("ABCD"))
        msg = "The truth value of a DataFrame is ambiguous"
        with pytest.raises(ValueError, match=msg):
            df in [None]

    def test_comparison_invalid(self):
        def check(df, df2):

            for (x, y) in [(df, df2), (df2, df)]:
                # we expect the result to match Series comparisons for
                # == and !=, inequalities should raise
                result = x == y
                expected = pd.DataFrame(
                    {col: x[col] == y[col] for col in x.columns},
                    index=x.index,
                    columns=x.columns,
                )
                tm.assert_frame_equal(result, expected)

                result = x != y
                expected = pd.DataFrame(
                    {col: x[col] != y[col] for col in x.columns},
                    index=x.index,
                    columns=x.columns,
                )
                tm.assert_frame_equal(result, expected)

                msg = re.escape(
                    "Invalid comparison between dtype=datetime64[ns] and ndarray"
                )
                with pytest.raises(TypeError, match=msg):
                    x >= y
                with pytest.raises(TypeError, match=msg):
                    x > y
                with pytest.raises(TypeError, match=msg):
                    x < y
                with pytest.raises(TypeError, match=msg):
                    x <= y

        # GH4968
        # invalid date/int comparisons
        df = pd.DataFrame(np.random.randint(10, size=(10, 1)), columns=["a"])
        df["dates"] = pd.date_range("20010101", periods=len(df))

        df2 = df.copy()
        df2["dates"] = df["a"]
        check(df, df2)

        df = pd.DataFrame(np.random.randint(10, size=(10, 2)), columns=["a", "b"])
        df2 = pd.DataFrame(
            {
                "a": pd.date_range("20010101", periods=len(df)),
                "b": pd.date_range("20100101", periods=len(df)),
            }
        )
        check(df, df2)

    def test_timestamp_compare(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH#4982
        df = pd.DataFrame(
            {
                "dates1": pd.date_range("20010101", periods=10),
                "dates2": pd.date_range("20010102", periods=10),
                "intcol": np.random.randint(1000000000, size=10),
                "floatcol": np.random.randn(10),
                "stringcol": list(tm.rands(10)),
            }
        )
        df.loc[np.random.rand(len(df)) > 0.5, "dates2"] = pd.NaT
        ops = {"gt": "lt", "lt": "gt", "ge": "le", "le": "ge", "eq": "eq", "ne": "ne"}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            # no nats
            if left in ["eq", "ne"]:
                expected = left_f(df, pd.Timestamp("20010109"))
                result = right_f(pd.Timestamp("20010109"), df)
                tm.assert_frame_equal(result, expected)
            else:
                msg = (
                    "'(<|>)=?' not supported between "
                    "instances of 'Timestamp' and 'float'"
                )
                with pytest.raises(TypeError, match=msg):
                    left_f(df, pd.Timestamp("20010109"))
                with pytest.raises(TypeError, match=msg):
                    right_f(pd.Timestamp("20010109"), df)
            # nats
            expected = left_f(df, pd.Timestamp("nat"))
            result = right_f(pd.Timestamp("nat"), df)
            tm.assert_frame_equal(result, expected)

    def test_mixed_comparison(self):
        # GH#13128, GH#22163 != datetime64 vs non-dt64 should be False,
        # not raise TypeError
        # (this appears to be fixed before GH#22163, not sure when)
        df = pd.DataFrame([["1989-08-01", 1], ["1989-08-01", 2]])
        other = pd.DataFrame([["a", "b"], ["c", "d"]])

        result = df == other
        assert not result.any().any()

        result = df != other
        assert result.all().all()

    def test_df_boolean_comparison_error(self):
        # GH#4576, GH#22880
        # comparing DataFrame against list/tuple with len(obj) matching
        #  len(df.columns) is supported as of GH#22800
        df = pd.DataFrame(np.arange(6).reshape((3, 2)))

        expected = pd.DataFrame([[False, False], [True, False], [False, False]])

        result = df == (2, 2)
        tm.assert_frame_equal(result, expected)

        result = df == [2, 2]
        tm.assert_frame_equal(result, expected)

    def test_df_float_none_comparison(self):
        df = pd.DataFrame(
            np.random.randn(8, 3), index=range(8), columns=["A", "B", "C"]
        )

        result = df.__eq__(None)
        assert not result.any().any()

    def test_df_string_comparison(self):
        df = pd.DataFrame([{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}])
        mask_a = df.a > 1
        tm.assert_frame_equal(df[mask_a], df.loc[1:1, :])
        tm.assert_frame_equal(df[-mask_a], df.loc[0:0, :])

        mask_b = df.b == "foo"
        tm.assert_frame_equal(df[mask_b], df.loc[0:0, :])
        tm.assert_frame_equal(df[-mask_b], df.loc[1:1, :])


class TestFrameFlexComparisons:
    # TODO: test_bool_flex_frame needs a better name
    def test_bool_flex_frame(self):
        data = np.random.randn(5, 3)
        other_data = np.random.randn(5, 3)
        df = pd.DataFrame(data)
        other = pd.DataFrame(other_data)
        ndim_5 = np.ones(df.shape + (1, 3))

        # Unaligned
        def _check_unaligned_frame(meth, op, df, other):
            part_o = other.loc[3:, 1:].copy()
            rs = meth(part_o)
            xp = op(df, part_o.reindex(index=df.index, columns=df.columns))
            tm.assert_frame_equal(rs, xp)

        # DataFrame
        assert df.eq(df).values.all()
        assert not df.ne(df).values.any()
        for op in ["eq", "ne", "gt", "lt", "ge", "le"]:
            f = getattr(df, op)
            o = getattr(operator, op)
            # No NAs
            tm.assert_frame_equal(f(other), o(df, other))
            _check_unaligned_frame(f, o, df, other)
            # ndarray
            tm.assert_frame_equal(f(other.values), o(df, other.values))
            # scalar
            tm.assert_frame_equal(f(0), o(df, 0))
            # NAs
            msg = "Unable to coerce to Series/DataFrame"
            tm.assert_frame_equal(f(np.nan), o(df, np.nan))
            with pytest.raises(ValueError, match=msg):
                f(ndim_5)

        # Series
        def _test_seq(df, idx_ser, col_ser):
            idx_eq = df.eq(idx_ser, axis=0)
            col_eq = df.eq(col_ser)
            idx_ne = df.ne(idx_ser, axis=0)
            col_ne = df.ne(col_ser)
            tm.assert_frame_equal(col_eq, df == pd.Series(col_ser))
            tm.assert_frame_equal(col_eq, -col_ne)
            tm.assert_frame_equal(idx_eq, -idx_ne)
            tm.assert_frame_equal(idx_eq, df.T.eq(idx_ser).T)
            tm.assert_frame_equal(col_eq, df.eq(list(col_ser)))
            tm.assert_frame_equal(idx_eq, df.eq(pd.Series(idx_ser), axis=0))
            tm.assert_frame_equal(idx_eq, df.eq(list(idx_ser), axis=0))

            idx_gt = df.gt(idx_ser, axis=0)
            col_gt = df.gt(col_ser)
            idx_le = df.le(idx_ser, axis=0)
            col_le = df.le(col_ser)

            tm.assert_frame_equal(col_gt, df > pd.Series(col_ser))
            tm.assert_frame_equal(col_gt, -col_le)
            tm.assert_frame_equal(idx_gt, -idx_le)
            tm.assert_frame_equal(idx_gt, df.T.gt(idx_ser).T)

            idx_ge = df.ge(idx_ser, axis=0)
            col_ge = df.ge(col_ser)
            idx_lt = df.lt(idx_ser, axis=0)
            col_lt = df.lt(col_ser)
            tm.assert_frame_equal(col_ge, df >= pd.Series(col_ser))
            tm.assert_frame_equal(col_ge, -col_lt)
            tm.assert_frame_equal(idx_ge, -idx_lt)
            tm.assert_frame_equal(idx_ge, df.T.ge(idx_ser).T)

        idx_ser = pd.Series(np.random.randn(5))
        col_ser = pd.Series(np.random.randn(3))
        _test_seq(df, idx_ser, col_ser)

        # list/tuple
        _test_seq(df, idx_ser.values, col_ser.values)

        # NA
        df.loc[0, 0] = np.nan
        rs = df.eq(df)
        assert not rs.loc[0, 0]
        rs = df.ne(df)
        assert rs.loc[0, 0]
        rs = df.gt(df)
        assert not rs.loc[0, 0]
        rs = df.lt(df)
        assert not rs.loc[0, 0]
        rs = df.ge(df)
        assert not rs.loc[0, 0]
        rs = df.le(df)
        assert not rs.loc[0, 0]

    def test_bool_flex_frame_complex_dtype(self):
        # complex
        arr = np.array([np.nan, 1, 6, np.nan])
        arr2 = np.array([2j, np.nan, 7, None])
        df = pd.DataFrame({"a": arr})
        df2 = pd.DataFrame({"a": arr2})

        msg = "|".join(
            [
                "'>' not supported between instances of '.*' and 'complex'",
                r"unorderable types: .*complex\(\)",  # PY35
            ]
        )
        with pytest.raises(TypeError, match=msg):
            # inequalities are not well-defined for complex numbers
            df.gt(df2)
        with pytest.raises(TypeError, match=msg):
            # regression test that we get the same behavior for Series
            df["a"].gt(df2["a"])
        with pytest.raises(TypeError, match=msg):
            # Check that we match numpy behavior here
            df.values > df2.values

        rs = df.ne(df2)
        assert rs.values.all()

        arr3 = np.array([2j, np.nan, None])
        df3 = pd.DataFrame({"a": arr3})

        with pytest.raises(TypeError, match=msg):
            # inequalities are not well-defined for complex numbers
            df3.gt(2j)
        with pytest.raises(TypeError, match=msg):
            # regression test that we get the same behavior for Series
            df3["a"].gt(2j)
        with pytest.raises(TypeError, match=msg):
            # Check that we match numpy behavior here
            df3.values > 2j

    def test_bool_flex_frame_object_dtype(self):
        # corner, dtype=object
        df1 = pd.DataFrame({"col": ["foo", np.nan, "bar"]})
        df2 = pd.DataFrame({"col": ["foo", datetime.now(), "bar"]})
        result = df1.ne(df2)
        exp = pd.DataFrame({"col": [False, True, False]})
        tm.assert_frame_equal(result, exp)

    def test_flex_comparison_nat(self):
        # GH 15697, GH 22163 df.eq(pd.NaT) should behave like df == pd.NaT,
        # and _definitely_ not be NaN
        df = pd.DataFrame([pd.NaT])

        result = df == pd.NaT
        # result.iloc[0, 0] is a np.bool_ object
        assert result.iloc[0, 0].item() is False

        result = df.eq(pd.NaT)
        assert result.iloc[0, 0].item() is False

        result = df != pd.NaT
        assert result.iloc[0, 0].item() is True

        result = df.ne(pd.NaT)
        assert result.iloc[0, 0].item() is True

    @pytest.mark.parametrize("opname", ["eq", "ne", "gt", "lt", "ge", "le"])
    def test_df_flex_cmp_constant_return_types(self, opname):
        # GH 15077, non-empty DataFrame
        df = pd.DataFrame({"x": [1, 2, 3], "y": [1.0, 2.0, 3.0]})
        const = 2

        result = getattr(df, opname)(const).dtypes.value_counts()
        tm.assert_series_equal(result, pd.Series([2], index=[np.dtype(bool)]))

    @pytest.mark.parametrize("opname", ["eq", "ne", "gt", "lt", "ge", "le"])
    def test_df_flex_cmp_constant_return_types_empty(self, opname):
        # GH 15077 empty DataFrame
        df = pd.DataFrame({"x": [1, 2, 3], "y": [1.0, 2.0, 3.0]})
        const = 2

        empty = df.iloc[:0]
        result = getattr(empty, opname)(const).dtypes.value_counts()
        tm.assert_series_equal(result, pd.Series([2], index=[np.dtype(bool)]))


# -------------------------------------------------------------------
# Arithmetic


class TestFrameFlexArithmetic:
    def test_floordiv_axis0(self):
        # make sure we df.floordiv(ser, axis=0) matches column-wise result
        arr = np.arange(3)
        ser = pd.Series(arr)
        df = pd.DataFrame({"A": ser, "B": ser})

        result = df.floordiv(ser, axis=0)

        expected = pd.DataFrame({col: df[col] // ser for col in df.columns})

        tm.assert_frame_equal(result, expected)

        result2 = df.floordiv(ser.values, axis=0)
        tm.assert_frame_equal(result2, expected)

    @pytest.mark.slow
    @pytest.mark.parametrize("opname", ["floordiv", "pow"])
    def test_floordiv_axis0_numexpr_path(self, opname):
        # case that goes through numexpr and has to fall back to masked_arith_op
        op = getattr(operator, opname)

        arr = np.arange(10 ** 6).reshape(100, -1)
        df = pd.DataFrame(arr)
        df["C"] = 1.0

        ser = df[0]
        result = getattr(df, opname)(ser, axis=0)

        expected = pd.DataFrame({col: op(df[col], ser) for col in df.columns})
        tm.assert_frame_equal(result, expected)

        result2 = getattr(df, opname)(ser.values, axis=0)
        tm.assert_frame_equal(result2, expected)

    def test_df_add_td64_columnwise(self):
        # GH 22534 Check that column-wise addition broadcasts correctly
        dti = pd.date_range("2016-01-01", periods=10)
        tdi = pd.timedelta_range("1", periods=10)
        tser = pd.Series(tdi)
        df = pd.DataFrame({0: dti, 1: tdi})

        result = df.add(tser, axis=0)
        expected = pd.DataFrame({0: dti + tdi, 1: tdi + tdi})
        tm.assert_frame_equal(result, expected)

    def test_df_add_flex_filled_mixed_dtypes(self):
        # GH 19611
        dti = pd.date_range("2016-01-01", periods=3)
        ser = pd.Series(["1 Day", "NaT", "2 Days"], dtype="timedelta64[ns]")
        df = pd.DataFrame({"A": dti, "B": ser})
        other = pd.DataFrame({"A": ser, "B": ser})
        fill = pd.Timedelta(days=1).to_timedelta64()
        result = df.add(other, fill_value=fill)

        expected = pd.DataFrame(
            {
                "A": pd.Series(
                    ["2016-01-02", "2016-01-03", "2016-01-05"], dtype="datetime64[ns]"
                ),
                "B": ser * 2,
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_frame(
        self, all_arithmetic_operators, float_frame, mixed_float_frame
    ):
        # one instance of parametrized fixture
        op = all_arithmetic_operators

        def f(x, y):
            # r-versions not in operator-stdlib; get op without "r" and invert
            if op.startswith("__r"):
                return getattr(operator, op.replace("__r", "__"))(y, x)
            return getattr(operator, op)(x, y)

        result = getattr(float_frame, op)(2 * float_frame)
        expected = f(float_frame, 2 * float_frame)
        tm.assert_frame_equal(result, expected)

        # vs mix float
        result = getattr(mixed_float_frame, op)(2 * mixed_float_frame)
        expected = f(mixed_float_frame, 2 * mixed_float_frame)
        tm.assert_frame_equal(result, expected)
        _check_mixed_float(result, dtype=dict(C=None))

    @pytest.mark.parametrize("op", ["__add__", "__sub__", "__mul__"])
    def test_arith_flex_frame_mixed(
        self, op, int_frame, mixed_int_frame, mixed_float_frame
    ):
        f = getattr(operator, op)

        # vs mix int
        result = getattr(mixed_int_frame, op)(2 + mixed_int_frame)
        expected = f(mixed_int_frame, 2 + mixed_int_frame)

        # no overflow in the uint
        dtype = None
        if op in ["__sub__"]:
            dtype = dict(B="uint64", C=None)
        elif op in ["__add__", "__mul__"]:
            dtype = dict(C=None)
        tm.assert_frame_equal(result, expected)
        _check_mixed_int(result, dtype=dtype)

        # vs mix float
        result = getattr(mixed_float_frame, op)(2 * mixed_float_frame)
        expected = f(mixed_float_frame, 2 * mixed_float_frame)
        tm.assert_frame_equal(result, expected)
        _check_mixed_float(result, dtype=dict(C=None))

        # vs plain int
        result = getattr(int_frame, op)(2 * int_frame)
        expected = f(int_frame, 2 * int_frame)
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_frame_raise(self, all_arithmetic_operators, float_frame):
        # one instance of parametrized fixture
        op = all_arithmetic_operators

        # Check that arrays with dim >= 3 raise
        for dim in range(3, 6):
            arr = np.ones((1,) * dim)
            msg = "Unable to coerce to Series/DataFrame"
            with pytest.raises(ValueError, match=msg):
                getattr(float_frame, op)(arr)

    def test_arith_flex_frame_corner(self, float_frame):

        const_add = float_frame.add(1)
        tm.assert_frame_equal(const_add, float_frame + 1)

        # corner cases
        result = float_frame.add(float_frame[:0])
        tm.assert_frame_equal(result, float_frame * np.nan)

        result = float_frame[:0].add(float_frame)
        tm.assert_frame_equal(result, float_frame * np.nan)

        with pytest.raises(NotImplementedError, match="fill_value"):
            float_frame.add(float_frame.iloc[0], fill_value=3)

        with pytest.raises(NotImplementedError, match="fill_value"):
            float_frame.add(float_frame.iloc[0], axis="index", fill_value=3)

    def test_arith_flex_series(self, simple_frame):
        df = simple_frame

        row = df.xs("a")
        col = df["two"]
        # after arithmetic refactor, add truediv here
        ops = ["add", "sub", "mul", "mod"]
        for op in ops:
            f = getattr(df, op)
            op = getattr(operator, op)
            tm.assert_frame_equal(f(row), op(df, row))
            tm.assert_frame_equal(f(col, axis=0), op(df.T, col).T)

        # special case for some reason
        tm.assert_frame_equal(df.add(row, axis=None), df + row)

        # cases which will be refactored after big arithmetic refactor
        tm.assert_frame_equal(df.div(row), df / row)
        tm.assert_frame_equal(df.div(col, axis=0), (df.T / col).T)

        # broadcasting issue in GH 7325
        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype="int64")
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0], axis="index")
        tm.assert_frame_equal(result, expected)

        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype="float64")
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0], axis="index")
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_zero_len_raises(self):
        # GH 19522 passing fill_value to frame flex arith methods should
        # raise even in the zero-length special cases
        ser_len0 = pd.Series([], dtype=object)
        df_len0 = pd.DataFrame(columns=["A", "B"])
        df = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"])

        with pytest.raises(NotImplementedError, match="fill_value"):
            df.add(ser_len0, fill_value="E")

        with pytest.raises(NotImplementedError, match="fill_value"):
            df_len0.sub(df["A"], axis=None, fill_value=3)


class TestFrameArithmetic:
    def test_td64_op_nat_casting(self):
        # Make sure we don't accidentally treat timedelta64(NaT) as datetime64
        #  when calling dispatch_to_series in DataFrame arithmetic
        ser = pd.Series(["NaT", "NaT"], dtype="timedelta64[ns]")
        df = pd.DataFrame([[1, 2], [3, 4]])

        result = df * ser
        expected = pd.DataFrame({0: ser, 1: ser})
        tm.assert_frame_equal(result, expected)

    def test_df_add_2d_array_rowlike_broadcasts(self):
        # GH#23000
        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=["A", "B", "C"])

        rowlike = arr[[1], :]  # shape --> (1, ncols)
        assert rowlike.shape == (1, df.shape[1])

        expected = pd.DataFrame(
            [[2, 4], [4, 6], [6, 8]],
            columns=df.columns,
            index=df.index,
            # specify dtype explicitly to avoid failing
            # on 32bit builds
            dtype=arr.dtype,
        )
        result = df + rowlike
        tm.assert_frame_equal(result, expected)
        result = rowlike + df
        tm.assert_frame_equal(result, expected)

    def test_df_add_2d_array_collike_broadcasts(self):
        # GH#23000
        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=["A", "B", "C"])

        collike = arr[:, [1]]  # shape --> (nrows, 1)
        assert collike.shape == (df.shape[0], 1)

        expected = pd.DataFrame(
            [[1, 2], [5, 6], [9, 10]],
            columns=df.columns,
            index=df.index,
            # specify dtype explicitly to avoid failing
            # on 32bit builds
            dtype=arr.dtype,
        )
        result = df + collike
        tm.assert_frame_equal(result, expected)
        result = collike + df
        tm.assert_frame_equal(result, expected)

    def test_df_arith_2d_array_rowlike_broadcasts(self, all_arithmetic_operators):
        # GH#23000
        opname = all_arithmetic_operators

        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=["A", "B", "C"])

        rowlike = arr[[1], :]  # shape --> (1, ncols)
        assert rowlike.shape == (1, df.shape[1])

        exvals = [
            getattr(df.loc["A"], opname)(rowlike.squeeze()),
            getattr(df.loc["B"], opname)(rowlike.squeeze()),
            getattr(df.loc["C"], opname)(rowlike.squeeze()),
        ]

        expected = pd.DataFrame(exvals, columns=df.columns, index=df.index)

        if opname in ["__rmod__", "__rfloordiv__"]:
            # exvals will have dtypes [f8, i8, i8] so expected will be
            #   all-f8, but the DataFrame operation will return mixed dtypes
            # use exvals[-1].dtype instead of "i8" for compat with 32-bit
            # systems/pythons
            expected[False] = expected[False].astype(exvals[-1].dtype)

        result = getattr(df, opname)(rowlike)
        tm.assert_frame_equal(result, expected)

    def test_df_arith_2d_array_collike_broadcasts(self, all_arithmetic_operators):
        # GH#23000
        opname = all_arithmetic_operators

        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=["A", "B", "C"])

        collike = arr[:, [1]]  # shape --> (nrows, 1)
        assert collike.shape == (df.shape[0], 1)

        exvals = {
            True: getattr(df[True], opname)(collike.squeeze()),
            False: getattr(df[False], opname)(collike.squeeze()),
        }

        dtype = None
        if opname in ["__rmod__", "__rfloordiv__"]:
            # Series ops may return mixed int/float dtypes in cases where
            #   DataFrame op will return all-float.  So we upcast `expected`
            dtype = np.common_type(*[x.values for x in exvals.values()])

        expected = pd.DataFrame(exvals, columns=df.columns, index=df.index, dtype=dtype)

        result = getattr(df, opname)(collike)
        tm.assert_frame_equal(result, expected)

    def test_df_bool_mul_int(self):
        # GH 22047, GH 22163 multiplication by 1 should result in int dtype,
        # not object dtype
        df = pd.DataFrame([[False, True], [False, False]])
        result = df * 1

        # On appveyor this comes back as np.int32 instead of np.int64,
        # so we check dtype.kind instead of just dtype
        kinds = result.dtypes.apply(lambda x: x.kind)
        assert (kinds == "i").all()

        result = 1 * df
        kinds = result.dtypes.apply(lambda x: x.kind)
        assert (kinds == "i").all()

    def test_arith_mixed(self):

        left = pd.DataFrame({"A": ["a", "b", "c"], "B": [1, 2, 3]})

        result = left + left
        expected = pd.DataFrame({"A": ["aa", "bb", "cc"], "B": [2, 4, 6]})
        tm.assert_frame_equal(result, expected)

    def test_arith_getitem_commute(self):
        df = pd.DataFrame({"A": [1.1, 3.3], "B": [2.5, -3.9]})

        def _test_op(df, op):
            result = op(df, 1)

            if not df.columns.is_unique:
                raise ValueError("Only unique columns supported by this test")

            for col in result.columns:
                tm.assert_series_equal(result[col], op(df[col], 1))

        _test_op(df, operator.add)
        _test_op(df, operator.sub)
        _test_op(df, operator.mul)
        _test_op(df, operator.truediv)
        _test_op(df, operator.floordiv)
        _test_op(df, operator.pow)

        _test_op(df, lambda x, y: y + x)
        _test_op(df, lambda x, y: y - x)
        _test_op(df, lambda x, y: y * x)
        _test_op(df, lambda x, y: y / x)
        _test_op(df, lambda x, y: y ** x)

        _test_op(df, lambda x, y: x + y)
        _test_op(df, lambda x, y: x - y)
        _test_op(df, lambda x, y: x * y)
        _test_op(df, lambda x, y: x / y)
        _test_op(df, lambda x, y: x ** y)

    @pytest.mark.parametrize(
        "values", [[1, 2], (1, 2), np.array([1, 2]), range(1, 3), deque([1, 2])]
    )
    def test_arith_alignment_non_pandas_object(self, values):
        # GH#17901
        df = pd.DataFrame({"A": [1, 1], "B": [1, 1]})
        expected = pd.DataFrame({"A": [2, 2], "B": [3, 3]})
        result = df + values
        tm.assert_frame_equal(result, expected)

    def test_arith_non_pandas_object(self):
        df = pd.DataFrame(
            np.arange(1, 10, dtype="f8").reshape(3, 3),
            columns=["one", "two", "three"],
            index=["a", "b", "c"],
        )

        val1 = df.xs("a").values
        added = pd.DataFrame(df.values + val1, index=df.index, columns=df.columns)
        tm.assert_frame_equal(df + val1, added)

        added = pd.DataFrame((df.values.T + val1).T, index=df.index, columns=df.columns)
        tm.assert_frame_equal(df.add(val1, axis=0), added)

        val2 = list(df["two"])

        added = pd.DataFrame(df.values + val2, index=df.index, columns=df.columns)
        tm.assert_frame_equal(df + val2, added)

        added = pd.DataFrame((df.values.T + val2).T, index=df.index, columns=df.columns)
        tm.assert_frame_equal(df.add(val2, axis="index"), added)

        val3 = np.random.rand(*df.shape)
        added = pd.DataFrame(df.values + val3, index=df.index, columns=df.columns)
        tm.assert_frame_equal(df.add(val3), added)

    def test_operations_with_interval_categories_index(self, all_arithmetic_operators):
        # GH#27415
        op = all_arithmetic_operators
        ind = pd.CategoricalIndex(pd.interval_range(start=0.0, end=2.0))
        data = [1, 2]
        df = pd.DataFrame([data], columns=ind)
        num = 10
        result = getattr(df, op)(num)
        expected = pd.DataFrame([[getattr(n, op)(num) for n in data]], columns=ind)
        tm.assert_frame_equal(result, expected)

    def test_frame_with_frame_reindex(self):
        # GH#31623
        df = pd.DataFrame(
            {
                "foo": [pd.Timestamp("2019"), pd.Timestamp("2020")],
                "bar": [pd.Timestamp("2018"), pd.Timestamp("2021")],
            },
            columns=["foo", "bar"],
        )
        df2 = df[["foo"]]

        result = df - df2

        expected = pd.DataFrame(
            {"foo": [pd.Timedelta(0), pd.Timedelta(0)], "bar": [np.nan, np.nan]},
            columns=["bar", "foo"],
        )
        tm.assert_frame_equal(result, expected)


def test_frame_with_zero_len_series_corner_cases():
    # GH#28600
    # easy all-float case
    df = pd.DataFrame(np.random.randn(6).reshape(3, 2), columns=["A", "B"])
    ser = pd.Series(dtype=np.float64)

    result = df + ser
    expected = pd.DataFrame(df.values * np.nan, columns=df.columns)
    tm.assert_frame_equal(result, expected)

    result = df == ser
    expected = pd.DataFrame(False, index=df.index, columns=df.columns)
    tm.assert_frame_equal(result, expected)

    # non-float case should not raise on comparison
    df2 = pd.DataFrame(df.values.view("M8[ns]"), columns=df.columns)
    result = df2 == ser
    expected = pd.DataFrame(False, index=df.index, columns=df.columns)
    tm.assert_frame_equal(result, expected)


def test_zero_len_frame_with_series_corner_cases():
    # GH#28600
    df = pd.DataFrame(columns=["A", "B"], dtype=np.float64)
    ser = pd.Series([1, 2], index=["A", "B"])

    result = df + ser
    expected = df
    tm.assert_frame_equal(result, expected)


def test_frame_single_columns_object_sum_axis_1():
    # GH 13758
    data = {
        "One": pd.Series(["A", 1.2, np.nan]),
    }
    df = pd.DataFrame(data)
    result = df.sum(axis=1)
    expected = pd.Series(["A", 1.2, 0])
    tm.assert_series_equal(result, expected)


# -------------------------------------------------------------------
# Unsorted
#  These arithmetic tests were previously in other files, eventually
#  should be parametrized and put into tests.arithmetic


class TestFrameArithmeticUnsorted:
    def test_frame_add_tz_mismatch_converts_to_utc(self):
        rng = pd.date_range("1/1/2011", periods=10, freq="H", tz="US/Eastern")
        df = pd.DataFrame(np.random.randn(len(rng)), index=rng, columns=["a"])

        df_moscow = df.tz_convert("Europe/Moscow")
        result = df + df_moscow
        assert result.index.tz is pytz.utc

        result = df_moscow + df
        assert result.index.tz is pytz.utc

    def test_align_frame(self):
        rng = pd.period_range("1/1/2000", "1/1/2010", freq="A")
        ts = pd.DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts + ts[::2]
        expected = ts + ts
        expected.values[1::2] = np.nan
        tm.assert_frame_equal(result, expected)

        half = ts[::2]
        result = ts + half.take(np.random.permutation(len(half)))
        tm.assert_frame_equal(result, expected)


def test_pow_with_realignment():
    # GH#32685 pow has special semantics for operating with null values
    left = pd.DataFrame({"A": [0, 1, 2]})
    right = pd.DataFrame(index=[0, 1, 2])

    result = left ** right
    expected = pd.DataFrame({"A": [np.nan, 1.0, np.nan]})
    tm.assert_frame_equal(result, expected)


# TODO: move to tests.arithmetic and parametrize
def test_pow_nan_with_zero():
    left = pd.DataFrame({"A": [np.nan, np.nan, np.nan]})
    right = pd.DataFrame({"A": [0, 0, 0]})

    expected = pd.DataFrame({"A": [1.0, 1.0, 1.0]})

    result = left ** right
    tm.assert_frame_equal(result, expected)

    result = left["A"] ** right["A"]
    tm.assert_series_equal(result, expected["A"])
