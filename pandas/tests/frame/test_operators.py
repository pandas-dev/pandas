from decimal import Decimal
import operator

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, MultiIndex, Series
import pandas.core.common as com
from pandas.tests.frame.common import _check_mixed_float
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_frame_equal,
    assert_numpy_array_equal,
    assert_series_equal,
)


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
        assert_frame_equal(-df, expected)
        assert_series_equal(-df["a"], expected["a"])

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
        assert_frame_equal(-df, expected)
        assert_series_equal(-df["a"], expected["a"])

    @pytest.mark.parametrize(
        "df",
        [
            pd.DataFrame({"a": ["a", "b"]}),
            pd.DataFrame({"a": pd.to_datetime(["2017-01-22", "1970-01-01"])}),
        ],
    )
    def test_neg_raises(self, df):
        with pytest.raises(TypeError):
            (-df)
        with pytest.raises(TypeError):
            (-df["a"])

    def test_invert(self, float_frame):
        df = float_frame

        assert_frame_equal(-(df < 0), ~(df < 0))

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
        assert_frame_equal(+df, df)
        assert_series_equal(+df["a"], df["a"])

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
        assert_frame_equal(+df, df)
        assert_series_equal(+df["a"], df["a"])

    @pytest.mark.parametrize(
        "df", [pd.DataFrame({"a": pd.to_datetime(["2017-01-22", "1970-01-01"])})]
    )
    def test_pos_raises(self, df):
        with pytest.raises(TypeError):
            (+df)
        with pytest.raises(TypeError):
            (+df["a"])


class TestDataFrameLogicalOperators:
    # &, |, ^

    def test_logical_ops_empty_frame(self):
        # GH#5808
        # empty frames, non-mixed dtype
        df = DataFrame(index=[1])

        result = df & df
        assert_frame_equal(result, df)

        result = df | df
        assert_frame_equal(result, df)

        df2 = DataFrame(index=[1, 2])
        result = df & df2
        assert_frame_equal(result, df2)

        dfa = DataFrame(index=[1], columns=["A"])

        result = dfa & dfa
        assert_frame_equal(result, dfa)

    def test_logical_ops_bool_frame(self):
        # GH#5808
        df1a_bool = DataFrame(True, index=[1], columns=["A"])

        result = df1a_bool & df1a_bool
        assert_frame_equal(result, df1a_bool)

        result = df1a_bool | df1a_bool
        assert_frame_equal(result, df1a_bool)

    def test_logical_ops_int_frame(self):
        # GH#5808
        df1a_int = DataFrame(1, index=[1], columns=["A"])
        df1a_bool = DataFrame(True, index=[1], columns=["A"])

        result = df1a_int | df1a_bool
        assert_frame_equal(result, df1a_int)

    def test_logical_ops_invalid(self):
        # GH#5808

        df1 = DataFrame(1.0, index=[1], columns=["A"])
        df2 = DataFrame(True, index=[1], columns=["A"])
        with pytest.raises(TypeError):
            df1 | df2

        df1 = DataFrame("foo", index=[1], columns=["A"])
        df2 = DataFrame(True, index=[1], columns=["A"])
        with pytest.raises(TypeError):
            df1 | df2

    def test_logical_operators(self):
        def _check_bin_op(op):
            result = op(df1, df2)
            expected = DataFrame(
                op(df1.values, df2.values), index=df1.index, columns=df1.columns
            )
            assert result.values.dtype == np.bool_
            assert_frame_equal(result, expected)

        def _check_unary_op(op):
            result = op(df1)
            expected = DataFrame(op(df1.values), index=df1.index, columns=df1.columns)
            assert result.values.dtype == np.bool_
            assert_frame_equal(result, expected)

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
        assert_series_equal(result, expected)

        # GH4604, automatic casting here
        result = d["a"].fillna(False) | d["b"]
        expected = Series([True, True])
        assert_series_equal(result, expected)

        result = d["a"].fillna(False, downcast=False) | d["b"]
        expected = Series([True, True])
        assert_series_equal(result, expected)


class TestDataFrameOperators:
    @pytest.mark.parametrize(
        "op", [operator.add, operator.sub, operator.mul, operator.truediv]
    )
    def test_operators_none_as_na(self, op):
        df = DataFrame(
            {"col1": [2, 5.0, 123, None], "col2": [1, 2, 3, 4]}, dtype=object
        )

        # since filling converts dtypes from object, changed expected to be
        # object
        filled = df.fillna(np.nan)
        result = op(df, 3)
        expected = op(filled, 3).astype(object)
        expected[com.isna(expected)] = None
        assert_frame_equal(result, expected)

        result = op(df, df)
        expected = op(filled, filled).astype(object)
        expected[com.isna(expected)] = None
        assert_frame_equal(result, expected)

        result = op(df, df.fillna(7))
        assert_frame_equal(result, expected)

        result = op(df.fillna(7), df)
        assert_frame_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize("op,res", [("__eq__", False), ("__ne__", True)])
    # TODO: not sure what's correct here.
    @pytest.mark.filterwarnings("ignore:elementwise:FutureWarning")
    def test_logical_typeerror_with_non_valid(self, op, res, float_frame):
        # we are comparing floats vs a string
        result = getattr(float_frame, op)("foo")
        assert bool(result.all().all()) is res

    def test_binary_ops_align(self):

        # test aligning binary ops

        # GH 6681
        index = MultiIndex.from_product(
            [list("abc"), ["one", "two", "three"], [1, 2, 3]],
            names=["first", "second", "third"],
        )

        df = DataFrame(
            np.arange(27 * 3).reshape(27, 3),
            index=index,
            columns=["value1", "value2", "value3"],
        ).sort_index()

        idx = pd.IndexSlice
        for op in ["add", "sub", "mul", "div", "truediv"]:
            opa = getattr(operator, op, None)
            if opa is None:
                continue

            x = Series([1.0, 10.0, 100.0], [1, 2, 3])
            result = getattr(df, op)(x, level="third", axis=0)

            expected = pd.concat(
                [opa(df.loc[idx[:, :, i], :], v) for i, v in x.items()]
            ).sort_index()
            assert_frame_equal(result, expected)

            x = Series([1.0, 10.0], ["two", "three"])
            result = getattr(df, op)(x, level="second", axis=0)

            expected = (
                pd.concat([opa(df.loc[idx[:, i], :], v) for i, v in x.items()])
                .reindex_like(df)
                .sort_index()
            )
            assert_frame_equal(result, expected)

        # GH9463 (alignment level of dataframe with series)

        midx = MultiIndex.from_product([["A", "B"], ["a", "b"]])
        df = DataFrame(np.ones((2, 4), dtype="int64"), columns=midx)
        s = pd.Series({"a": 1, "b": 2})

        df2 = df.copy()
        df2.columns.names = ["lvl0", "lvl1"]
        s2 = s.copy()
        s2.index.name = "lvl1"

        # different cases of integer/string level names:
        res1 = df.mul(s, axis=1, level=1)
        res2 = df.mul(s2, axis=1, level=1)
        res3 = df2.mul(s, axis=1, level=1)
        res4 = df2.mul(s2, axis=1, level=1)
        res5 = df2.mul(s, axis=1, level="lvl1")
        res6 = df2.mul(s2, axis=1, level="lvl1")

        exp = DataFrame(
            np.array([[1, 2, 1, 2], [1, 2, 1, 2]], dtype="int64"), columns=midx
        )

        for res in [res1, res2]:
            assert_frame_equal(res, exp)

        exp.columns.names = ["lvl0", "lvl1"]
        for res in [res3, res4, res5, res6]:
            assert_frame_equal(res, exp)

    def test_dti_tz_convert_to_utc(self):
        base = pd.DatetimeIndex(["2011-01-01", "2011-01-02", "2011-01-03"], tz="UTC")
        idx1 = base.tz_convert("Asia/Tokyo")[:2]
        idx2 = base.tz_convert("US/Eastern")[1:]

        df1 = DataFrame({"A": [1, 2]}, index=idx1)
        df2 = DataFrame({"A": [1, 1]}, index=idx2)
        exp = DataFrame({"A": [np.nan, 3, np.nan]}, index=base)
        assert_frame_equal(df1 + df2, exp)

    def test_combineFrame(self, float_frame, mixed_float_frame, mixed_int_frame):
        frame_copy = float_frame.reindex(float_frame.index[::2])

        del frame_copy["D"]
        frame_copy["C"][:5] = np.nan

        added = float_frame + frame_copy

        indexer = added["A"].dropna().index
        exp = (float_frame["A"] * 2).copy()

        tm.assert_series_equal(added["A"].dropna(), exp.loc[indexer])

        exp.loc[~exp.index.isin(indexer)] = np.nan
        tm.assert_series_equal(added["A"], exp.loc[added["A"].index])

        assert np.isnan(added["C"].reindex(frame_copy.index)[:5]).all()

        # assert(False)

        assert np.isnan(added["D"]).all()

        self_added = float_frame + float_frame
        tm.assert_index_equal(self_added.index, float_frame.index)

        added_rev = frame_copy + float_frame
        assert np.isnan(added["D"]).all()
        assert np.isnan(added_rev["D"]).all()

        # corner cases

        # empty
        plus_empty = float_frame + DataFrame()
        assert np.isnan(plus_empty.values).all()

        empty_plus = DataFrame() + float_frame
        assert np.isnan(empty_plus.values).all()

        empty_empty = DataFrame() + DataFrame()
        assert empty_empty.empty

        # out of order
        reverse = float_frame.reindex(columns=float_frame.columns[::-1])

        assert_frame_equal(reverse + float_frame, float_frame * 2)

        # mix vs float64, upcast
        added = float_frame + mixed_float_frame
        _check_mixed_float(added, dtype="float64")
        added = mixed_float_frame + float_frame
        _check_mixed_float(added, dtype="float64")

        # mix vs mix
        added = mixed_float_frame + mixed_float_frame
        _check_mixed_float(added, dtype=dict(C=None))

        # with int
        added = float_frame + mixed_int_frame
        _check_mixed_float(added, dtype="float64")

    def test_combineSeries(
        self, float_frame, mixed_float_frame, mixed_int_frame, datetime_frame
    ):

        # Series
        series = float_frame.xs(float_frame.index[0])

        added = float_frame + series

        for key, s in added.items():
            assert_series_equal(s, float_frame[key] + series[key])

        larger_series = series.to_dict()
        larger_series["E"] = 1
        larger_series = Series(larger_series)
        larger_added = float_frame + larger_series

        for key, s in float_frame.items():
            assert_series_equal(larger_added[key], s + series[key])
        assert "E" in larger_added
        assert np.isnan(larger_added["E"]).all()

        # no upcast needed
        added = mixed_float_frame + series
        _check_mixed_float(added)

        # vs mix (upcast) as needed
        added = mixed_float_frame + series.astype("float32")
        _check_mixed_float(added, dtype=dict(C=None))
        added = mixed_float_frame + series.astype("float16")
        _check_mixed_float(added, dtype=dict(C=None))

        # these raise with numexpr.....as we are adding an int64 to an
        # uint64....weird vs int

        # added = mixed_int_frame + (100*series).astype('int64')
        # _check_mixed_int(added, dtype = dict(A = 'int64', B = 'float64', C =
        # 'int64', D = 'int64'))
        # added = mixed_int_frame + (100*series).astype('int32')
        # _check_mixed_int(added, dtype = dict(A = 'int32', B = 'float64', C =
        # 'int32', D = 'int64'))

        # TimeSeries
        ts = datetime_frame["A"]

        # 10890
        # we no longer allow auto timeseries broadcasting
        # and require explicit broadcasting
        added = datetime_frame.add(ts, axis="index")

        for key, col in datetime_frame.items():
            result = col + ts
            assert_series_equal(added[key], result, check_names=False)
            assert added[key].name == key
            if col.name == ts.name:
                assert result.name == "A"
            else:
                assert result.name is None

        smaller_frame = datetime_frame[:-5]
        smaller_added = smaller_frame.add(ts, axis="index")

        tm.assert_index_equal(smaller_added.index, datetime_frame.index)

        smaller_ts = ts[:-5]
        smaller_added2 = datetime_frame.add(smaller_ts, axis="index")
        assert_frame_equal(smaller_added, smaller_added2)

        # length 0, result is all-nan
        result = datetime_frame.add(ts[:0], axis="index")
        expected = DataFrame(
            np.nan, index=datetime_frame.index, columns=datetime_frame.columns
        )
        assert_frame_equal(result, expected)

        # Frame is all-nan
        result = datetime_frame[:0].add(ts, axis="index")
        expected = DataFrame(
            np.nan, index=datetime_frame.index, columns=datetime_frame.columns
        )
        assert_frame_equal(result, expected)

        # empty but with non-empty index
        frame = datetime_frame[:1].reindex(columns=[])
        result = frame.mul(ts, axis="index")
        assert len(result) == len(ts)

    def test_combineFunc(self, float_frame, mixed_float_frame):
        result = float_frame * 2
        tm.assert_numpy_array_equal(result.values, float_frame.values * 2)

        # vs mix
        result = mixed_float_frame * 2
        for c, s in result.items():
            tm.assert_numpy_array_equal(s.values, mixed_float_frame[c].values * 2)
        _check_mixed_float(result, dtype=dict(C=None))

        result = DataFrame() * 2
        assert result.index.equals(DataFrame().index)
        assert len(result.columns) == 0

    def test_comparisons(self, simple_frame, float_frame):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame()

        row = simple_frame.xs("a")
        ndim_5 = np.ones(df1.shape + (1, 1, 1))

        def test_comp(func):
            result = func(df1, df2)
            tm.assert_numpy_array_equal(result.values, func(df1.values, df2.values))

            with pytest.raises(ValueError, match="dim must be <= 2"):
                func(df1, ndim_5)

            result2 = func(simple_frame, row)
            tm.assert_numpy_array_equal(
                result2.values, func(simple_frame.values, row.values)
            )

            result3 = func(float_frame, 0)
            tm.assert_numpy_array_equal(result3.values, func(float_frame.values, 0))

            msg = "Can only compare identically-labeled DataFrame"
            with pytest.raises(ValueError, match=msg):
                func(simple_frame, simple_frame[:2])

        test_comp(operator.eq)
        test_comp(operator.ne)
        test_comp(operator.lt)
        test_comp(operator.gt)
        test_comp(operator.ge)
        test_comp(operator.le)

    def test_comparison_protected_from_errstate(self):
        missing_df = tm.makeDataFrame()
        missing_df.iloc[0]["A"] = np.nan
        with np.errstate(invalid="ignore"):
            expected = missing_df.values < 0
        with np.errstate(invalid="raise"):
            result = (missing_df < 0).values
        tm.assert_numpy_array_equal(result, expected)

    def test_boolean_comparison(self):

        # GH 4576
        # boolean comparisons with a tuple/list give unexpected results
        df = DataFrame(np.arange(6).reshape((3, 2)))
        b = np.array([2, 2])
        b_r = np.atleast_2d([2, 2])
        b_c = b_r.T
        lst = [2, 2, 2]
        tup = tuple(lst)

        # gt
        expected = DataFrame([[False, False], [False, True], [True, True]])
        result = df > b
        assert_frame_equal(result, expected)

        result = df.values > b
        assert_numpy_array_equal(result, expected.values)

        msg1d = "Unable to coerce to Series, length must be 2: given 3"
        msg2d = "Unable to coerce to DataFrame, shape must be"
        msg2db = "operands could not be broadcast together with shapes"
        with pytest.raises(ValueError, match=msg1d):
            # wrong shape
            df > lst

        with pytest.raises(ValueError, match=msg1d):
            # wrong shape
            result = df > tup

        # broadcasts like ndarray (GH#23000)
        result = df > b_r
        assert_frame_equal(result, expected)

        result = df.values > b_r
        assert_numpy_array_equal(result, expected.values)

        with pytest.raises(ValueError, match=msg2d):
            df > b_c

        with pytest.raises(ValueError, match=msg2db):
            df.values > b_c

        # ==
        expected = DataFrame([[False, False], [True, False], [False, False]])
        result = df == b
        assert_frame_equal(result, expected)

        with pytest.raises(ValueError, match=msg1d):
            result = df == lst

        with pytest.raises(ValueError, match=msg1d):
            result = df == tup

        # broadcasts like ndarray (GH#23000)
        result = df == b_r
        assert_frame_equal(result, expected)

        result = df.values == b_r
        assert_numpy_array_equal(result, expected.values)

        with pytest.raises(ValueError, match=msg2d):
            df == b_c

        assert df.values.shape != b_c.shape

        # with alignment
        df = DataFrame(
            np.arange(6).reshape((3, 2)), columns=list("AB"), index=list("abc")
        )
        expected.index = df.index
        expected.columns = df.columns

        with pytest.raises(ValueError, match=msg1d):
            result = df == lst

        with pytest.raises(ValueError, match=msg1d):
            result = df == tup

    def test_combine_generic(self, float_frame):
        df1 = float_frame
        df2 = float_frame.loc[float_frame.index[:-5], ["A", "B", "C"]]

        combined = df1.combine(df2, np.add)
        combined2 = df2.combine(df1, np.add)
        assert combined["D"].isna().all()
        assert combined2["D"].isna().all()

        chunk = combined.loc[combined.index[:-5], ["A", "B", "C"]]
        chunk2 = combined2.loc[combined2.index[:-5], ["A", "B", "C"]]

        exp = (
            float_frame.loc[float_frame.index[:-5], ["A", "B", "C"]].reindex_like(chunk)
            * 2
        )
        assert_frame_equal(chunk, exp)
        assert_frame_equal(chunk2, exp)

    def test_inplace_ops_alignment(self):

        # inplace ops / ops alignment
        # GH 8511

        columns = list("abcdefg")
        X_orig = DataFrame(
            np.arange(10 * len(columns)).reshape(-1, len(columns)),
            columns=columns,
            index=range(10),
        )
        Z = 100 * X_orig.iloc[:, 1:-1].copy()
        block1 = list("bedcf")
        subs = list("bcdef")

        # add
        X = X_orig.copy()
        result1 = (X[block1] + Z).reindex(columns=subs)

        X[block1] += Z
        result2 = X.reindex(columns=subs)

        X = X_orig.copy()
        result3 = (X[block1] + Z[block1]).reindex(columns=subs)

        X[block1] += Z[block1]
        result4 = X.reindex(columns=subs)

        assert_frame_equal(result1, result2)
        assert_frame_equal(result1, result3)
        assert_frame_equal(result1, result4)

        # sub
        X = X_orig.copy()
        result1 = (X[block1] - Z).reindex(columns=subs)

        X[block1] -= Z
        result2 = X.reindex(columns=subs)

        X = X_orig.copy()
        result3 = (X[block1] - Z[block1]).reindex(columns=subs)

        X[block1] -= Z[block1]
        result4 = X.reindex(columns=subs)

        assert_frame_equal(result1, result2)
        assert_frame_equal(result1, result3)
        assert_frame_equal(result1, result4)

    def test_inplace_ops_identity(self):

        # GH 5104
        # make sure that we are actually changing the object
        s_orig = Series([1, 2, 3])
        df_orig = DataFrame(np.random.randint(0, 5, size=10).reshape(-1, 5))

        # no dtype change
        s = s_orig.copy()
        s2 = s
        s += 1
        assert_series_equal(s, s2)
        assert_series_equal(s_orig + 1, s)
        assert s is s2
        assert s._data is s2._data

        df = df_orig.copy()
        df2 = df
        df += 1
        assert_frame_equal(df, df2)
        assert_frame_equal(df_orig + 1, df)
        assert df is df2
        assert df._data is df2._data

        # dtype change
        s = s_orig.copy()
        s2 = s
        s += 1.5
        assert_series_equal(s, s2)
        assert_series_equal(s_orig + 1.5, s)

        df = df_orig.copy()
        df2 = df
        df += 1.5
        assert_frame_equal(df, df2)
        assert_frame_equal(df_orig + 1.5, df)
        assert df is df2
        assert df._data is df2._data

        # mixed dtype
        arr = np.random.randint(0, 10, size=5)
        df_orig = DataFrame({"A": arr.copy(), "B": "foo"})
        df = df_orig.copy()
        df2 = df
        df["A"] += 1
        expected = DataFrame({"A": arr.copy() + 1, "B": "foo"})
        assert_frame_equal(df, expected)
        assert_frame_equal(df2, expected)
        assert df._data is df2._data

        df = df_orig.copy()
        df2 = df
        df["A"] += 1.5
        expected = DataFrame({"A": arr.copy() + 1.5, "B": "foo"})
        assert_frame_equal(df, expected)
        assert_frame_equal(df2, expected)
        assert df._data is df2._data

    @pytest.mark.parametrize(
        "op",
        [
            "add",
            "and",
            "div",
            "floordiv",
            "mod",
            "mul",
            "or",
            "pow",
            "sub",
            "truediv",
            "xor",
        ],
    )
    def test_inplace_ops_identity2(self, op):

        if op == "div":
            return

        df = DataFrame({"a": [1.0, 2.0, 3.0], "b": [1, 2, 3]})

        operand = 2
        if op in ("and", "or", "xor"):
            # cannot use floats for boolean ops
            df["a"] = [True, False, True]

        df_copy = df.copy()
        iop = "__i{}__".format(op)
        op = "__{}__".format(op)

        # no id change and value is correct
        getattr(df, iop)(operand)
        expected = getattr(df_copy, op)(operand)
        assert_frame_equal(df, expected)
        expected = id(df)
        assert id(df) == expected

    def test_alignment_non_pandas(self):
        index = ["A", "B", "C"]
        columns = ["X", "Y", "Z"]
        df = pd.DataFrame(np.random.randn(3, 3), index=index, columns=columns)

        align = pd.core.ops._align_method_FRAME
        for val in [
            [1, 2, 3],
            (1, 2, 3),
            np.array([1, 2, 3], dtype=np.int64),
            range(1, 4),
        ]:

            tm.assert_series_equal(
                align(df, val, "index"), Series([1, 2, 3], index=df.index)
            )
            tm.assert_series_equal(
                align(df, val, "columns"), Series([1, 2, 3], index=df.columns)
            )

        # length mismatch
        msg = "Unable to coerce to Series, length must be 3: given 2"
        for val in [[1, 2], (1, 2), np.array([1, 2]), range(1, 3)]:

            with pytest.raises(ValueError, match=msg):
                align(df, val, "index")

            with pytest.raises(ValueError, match=msg):
                align(df, val, "columns")

        val = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        tm.assert_frame_equal(
            align(df, val, "index"), DataFrame(val, index=df.index, columns=df.columns)
        )
        tm.assert_frame_equal(
            align(df, val, "columns"),
            DataFrame(val, index=df.index, columns=df.columns),
        )

        # shape mismatch
        msg = "Unable to coerce to DataFrame, shape must be"
        val = np.array([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError, match=msg):
            align(df, val, "index")

        with pytest.raises(ValueError, match=msg):
            align(df, val, "columns")

        val = np.zeros((3, 3, 3))
        with pytest.raises(ValueError):
            align(df, val, "index")
        with pytest.raises(ValueError):
            align(df, val, "columns")

    def test_no_warning(self, all_arithmetic_operators):
        df = pd.DataFrame({"A": [0.0, 0.0], "B": [0.0, None]})
        b = df["B"]
        with tm.assert_produces_warning(None):
            getattr(df, all_arithmetic_operators)(b, 0)


class TestTranspose:
    def test_transpose_tzaware_1col_single_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")

        df = pd.DataFrame(dti)
        assert (df.dtypes == dti.dtype).all()
        res = df.T
        assert (res.dtypes == dti.dtype).all()

    def test_transpose_tzaware_2col_single_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")

        df3 = pd.DataFrame({"A": dti, "B": dti})
        assert (df3.dtypes == dti.dtype).all()
        res3 = df3.T
        assert (res3.dtypes == dti.dtype).all()

    def test_transpose_tzaware_2col_mixed_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")
        dti2 = dti.tz_convert("US/Pacific")

        df4 = pd.DataFrame({"A": dti, "B": dti2})
        assert (df4.dtypes == [dti.dtype, dti2.dtype]).all()
        assert (df4.T.dtypes == object).all()
        tm.assert_frame_equal(df4.T.T, df4)

    def test_transpose_object_to_tzaware_mixed_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")
        dti2 = dti.tz_convert("US/Pacific")

        # mixed all-tzaware dtypes
        df2 = pd.DataFrame([dti, dti2])
        assert (df2.dtypes == object).all()
        res2 = df2.T
        assert (res2.dtypes == [dti.dtype, dti2.dtype]).all()
