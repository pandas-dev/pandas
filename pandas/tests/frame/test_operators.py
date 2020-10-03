from decimal import Decimal
import operator
import re

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm


class TestDataFrameUnaryOperators:
    # __pos__, __neg__, __inv__

    @pytest.mark.parametrize(
        "df,expected",
        [
            (pd.DataFrame({"a": [-1, 1]}), pd.DataFrame({"a": [1, -1]})),
            (pd.DataFrame({"a": [False, True]}), pd.DataFrame({"a": [True, False]})),
            (
                pd.DataFrame({"a": pd.Series(pd.to_timedelta([-1, 1]))}),
                pd.DataFrame({"a": pd.Series(pd.to_timedelta([1, -1]))}),
            ),
        ],
    )
    def test_neg_numeric(self, df, expected):
        tm.assert_frame_equal(-df, expected)
        tm.assert_series_equal(-df["a"], expected["a"])

    @pytest.mark.parametrize(
        "df, expected",
        [
            (np.array([1, 2], dtype=object), np.array([-1, -2], dtype=object)),
            ([Decimal("1.0"), Decimal("2.0")], [Decimal("-1.0"), Decimal("-2.0")]),
        ],
    )
    def test_neg_object(self, df, expected):
        # GH#21380
        df = pd.DataFrame({"a": df})
        expected = pd.DataFrame({"a": expected})
        tm.assert_frame_equal(-df, expected)
        tm.assert_series_equal(-df["a"], expected["a"])

    @pytest.mark.parametrize(
        "df",
        [
            pd.DataFrame({"a": ["a", "b"]}),
            pd.DataFrame({"a": pd.to_datetime(["2017-01-22", "1970-01-01"])}),
        ],
    )
    def test_neg_raises(self, df):
        msg = (
            "bad operand type for unary -: 'str'|"
            r"Unary negative expects numeric dtype, not datetime64\[ns\]"
        )
        with pytest.raises(TypeError, match=msg):
            (-df)
        with pytest.raises(TypeError, match=msg):
            (-df["a"])

    def test_invert(self, float_frame):
        df = float_frame

        tm.assert_frame_equal(-(df < 0), ~(df < 0))

    def test_invert_mixed(self):
        shape = (10, 5)
        df = pd.concat(
            [
                pd.DataFrame(np.zeros(shape, dtype="bool")),
                pd.DataFrame(np.zeros(shape, dtype=int)),
            ],
            axis=1,
            ignore_index=True,
        )
        result = ~df
        expected = pd.concat(
            [
                pd.DataFrame(np.ones(shape, dtype="bool")),
                pd.DataFrame(-np.ones(shape, dtype=int)),
            ],
            axis=1,
            ignore_index=True,
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "df",
        [
            pd.DataFrame({"a": [-1, 1]}),
            pd.DataFrame({"a": [False, True]}),
            pd.DataFrame({"a": pd.Series(pd.to_timedelta([-1, 1]))}),
        ],
    )
    def test_pos_numeric(self, df):
        # GH#16073
        tm.assert_frame_equal(+df, df)
        tm.assert_series_equal(+df["a"], df["a"])

    @pytest.mark.parametrize(
        "df",
        [
            # numpy changing behavior in the future
            pytest.param(
                pd.DataFrame({"a": ["a", "b"]}),
                marks=[pytest.mark.filterwarnings("ignore")],
            ),
            pd.DataFrame({"a": np.array([-1, 2], dtype=object)}),
            pd.DataFrame({"a": [Decimal("-1.0"), Decimal("2.0")]}),
        ],
    )
    def test_pos_object(self, df):
        # GH#21380
        tm.assert_frame_equal(+df, df)
        tm.assert_series_equal(+df["a"], df["a"])

    @pytest.mark.parametrize(
        "df", [pd.DataFrame({"a": pd.to_datetime(["2017-01-22", "1970-01-01"])})]
    )
    def test_pos_raises(self, df):
        msg = "Unary plus expects .* dtype, not datetime64\\[ns\\]"
        with pytest.raises(TypeError, match=msg):
            (+df)
        with pytest.raises(TypeError, match=msg):
            (+df["a"])


class TestDataFrameLogicalOperators:
    # &, |, ^

    def test_logical_ops_empty_frame(self):
        # GH#5808
        # empty frames, non-mixed dtype
        df = DataFrame(index=[1])

        result = df & df
        tm.assert_frame_equal(result, df)

        result = df | df
        tm.assert_frame_equal(result, df)

        df2 = DataFrame(index=[1, 2])
        result = df & df2
        tm.assert_frame_equal(result, df2)

        dfa = DataFrame(index=[1], columns=["A"])

        result = dfa & dfa
        expected = DataFrame(False, index=[1], columns=["A"])
        tm.assert_frame_equal(result, expected)

    def test_logical_ops_bool_frame(self):
        # GH#5808
        df1a_bool = DataFrame(True, index=[1], columns=["A"])

        result = df1a_bool & df1a_bool
        tm.assert_frame_equal(result, df1a_bool)

        result = df1a_bool | df1a_bool
        tm.assert_frame_equal(result, df1a_bool)

    def test_logical_ops_int_frame(self):
        # GH#5808
        df1a_int = DataFrame(1, index=[1], columns=["A"])
        df1a_bool = DataFrame(True, index=[1], columns=["A"])

        result = df1a_int | df1a_bool
        tm.assert_frame_equal(result, df1a_bool)

        # Check that this matches Series behavior
        res_ser = df1a_int["A"] | df1a_bool["A"]
        tm.assert_series_equal(res_ser, df1a_bool["A"])

    def test_logical_ops_invalid(self):
        # GH#5808

        df1 = DataFrame(1.0, index=[1], columns=["A"])
        df2 = DataFrame(True, index=[1], columns=["A"])
        msg = re.escape("unsupported operand type(s) for |: 'float' and 'bool'")
        with pytest.raises(TypeError, match=msg):
            df1 | df2

        df1 = DataFrame("foo", index=[1], columns=["A"])
        df2 = DataFrame(True, index=[1], columns=["A"])
        msg = re.escape("unsupported operand type(s) for |: 'str' and 'bool'")
        with pytest.raises(TypeError, match=msg):
            df1 | df2

    def test_logical_operators(self):
        def _check_bin_op(op):
            result = op(df1, df2)
            expected = DataFrame(
                op(df1.values, df2.values), index=df1.index, columns=df1.columns
            )
            assert result.values.dtype == np.bool_
            tm.assert_frame_equal(result, expected)

        def _check_unary_op(op):
            result = op(df1)
            expected = DataFrame(op(df1.values), index=df1.index, columns=df1.columns)
            assert result.values.dtype == np.bool_
            tm.assert_frame_equal(result, expected)

        df1 = {
            "a": {"a": True, "b": False, "c": False, "d": True, "e": True},
            "b": {"a": False, "b": True, "c": False, "d": False, "e": False},
            "c": {"a": False, "b": False, "c": True, "d": False, "e": False},
            "d": {"a": True, "b": False, "c": False, "d": True, "e": True},
            "e": {"a": True, "b": False, "c": False, "d": True, "e": True},
        }

        df2 = {
            "a": {"a": True, "b": False, "c": True, "d": False, "e": False},
            "b": {"a": False, "b": True, "c": False, "d": False, "e": False},
            "c": {"a": True, "b": False, "c": True, "d": False, "e": False},
            "d": {"a": False, "b": False, "c": False, "d": True, "e": False},
            "e": {"a": False, "b": False, "c": False, "d": False, "e": True},
        }

        df1 = DataFrame(df1)
        df2 = DataFrame(df2)

        _check_bin_op(operator.and_)
        _check_bin_op(operator.or_)
        _check_bin_op(operator.xor)

        _check_unary_op(operator.inv)  # TODO: belongs elsewhere

    def test_logical_with_nas(self):
        d = DataFrame({"a": [np.nan, False], "b": [True, True]})

        # GH4947
        # bool comparisons should return bool
        result = d["a"] | d["b"]
        expected = Series([False, True])
        tm.assert_series_equal(result, expected)

        # GH4604, automatic casting here
        result = d["a"].fillna(False) | d["b"]
        expected = Series([True, True])
        tm.assert_series_equal(result, expected)

        result = d["a"].fillna(False, downcast=False) | d["b"]
        expected = Series([True, True])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "left, right, op, expected",
        [
            (
                [True, False, np.nan],
                [True, False, True],
                operator.and_,
                [True, False, False],
            ),
            (
                [True, False, True],
                [True, False, np.nan],
                operator.and_,
                [True, False, False],
            ),
            (
                [True, False, np.nan],
                [True, False, True],
                operator.or_,
                [True, False, False],
            ),
            (
                [True, False, True],
                [True, False, np.nan],
                operator.or_,
                [True, False, True],
            ),
        ],
    )
    def test_logical_operators_nans(self, left, right, op, expected):
        # GH 13896
        result = op(DataFrame(left), DataFrame(right))
        expected = DataFrame(expected)

        tm.assert_frame_equal(result, expected)
