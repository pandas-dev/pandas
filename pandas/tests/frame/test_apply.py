from collections import OrderedDict
from datetime import datetime
from itertools import chain
import operator
import warnings

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import DataFrame, MultiIndex, Series, Timestamp, date_range, notna
import pandas._testing as tm
from pandas.conftest import _get_cython_table_params
from pandas.core.apply import frame_apply
from pandas.core.base import SpecificationError


@pytest.fixture
def int_frame_const_col():
    """
    Fixture for DataFrame of ints which are constant per column

    Columns are ['A', 'B', 'C'], with values (per column): [1, 2, 3]
    """
    df = DataFrame(
        np.tile(np.arange(3, dtype="int64"), 6).reshape(6, -1) + 1,
        columns=["A", "B", "C"],
    )
    return df


class TestDataFrameApply:
    def test_apply(self, float_frame):
        with np.errstate(all="ignore"):
            # ufunc
            applied = float_frame.apply(np.sqrt)
            tm.assert_series_equal(np.sqrt(float_frame["A"]), applied["A"])

            # aggregator
            applied = float_frame.apply(np.mean)
            assert applied["A"] == np.mean(float_frame["A"])

            d = float_frame.index[0]
            applied = float_frame.apply(np.mean, axis=1)
            assert applied[d] == np.mean(float_frame.xs(d))
            assert applied.index is float_frame.index  # want this

        # invalid axis
        df = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=["a", "a", "c"])
        with pytest.raises(ValueError):
            df.apply(lambda x: x, 2)

        # GH 9573
        df = DataFrame({"c0": ["A", "A", "B", "B"], "c1": ["C", "C", "D", "D"]})
        df = df.apply(lambda ts: ts.astype("category"))

        assert df.shape == (4, 2)
        assert isinstance(df["c0"].dtype, CategoricalDtype)
        assert isinstance(df["c1"].dtype, CategoricalDtype)

    def test_apply_mixed_datetimelike(self):
        # mixed datetimelike
        # GH 7778
        df = DataFrame(
            {
                "A": date_range("20130101", periods=3),
                "B": pd.to_timedelta(np.arange(3), unit="s"),
            }
        )
        result = df.apply(lambda x: x, axis=1)
        tm.assert_frame_equal(result, df)

    def test_apply_empty(self, float_frame):
        # empty
        empty_frame = DataFrame()

        applied = empty_frame.apply(np.sqrt)
        assert applied.empty

        applied = empty_frame.apply(np.mean)
        assert applied.empty

        no_rows = float_frame[:0]
        result = no_rows.apply(lambda x: x.mean())
        expected = Series(np.nan, index=float_frame.columns)
        tm.assert_series_equal(result, expected)

        no_cols = float_frame.loc[:, []]
        result = no_cols.apply(lambda x: x.mean(), axis=1)
        expected = Series(np.nan, index=float_frame.index)
        tm.assert_series_equal(result, expected)

        # GH 2476
        expected = DataFrame(index=["a"])
        result = expected.apply(lambda x: x["a"], axis=1)
        tm.assert_frame_equal(expected, result)

    def test_apply_with_reduce_empty(self):
        # reduce with an empty DataFrame
        empty_frame = DataFrame()

        x = []
        result = empty_frame.apply(x.append, axis=1, result_type="expand")
        tm.assert_frame_equal(result, empty_frame)
        result = empty_frame.apply(x.append, axis=1, result_type="reduce")
        expected = Series([], index=pd.Index([], dtype=object), dtype=np.float64)
        tm.assert_series_equal(result, expected)

        empty_with_cols = DataFrame(columns=["a", "b", "c"])
        result = empty_with_cols.apply(x.append, axis=1, result_type="expand")
        tm.assert_frame_equal(result, empty_with_cols)
        result = empty_with_cols.apply(x.append, axis=1, result_type="reduce")
        expected = Series([], index=pd.Index([], dtype=object), dtype=np.float64)
        tm.assert_series_equal(result, expected)

        # Ensure that x.append hasn't been called
        assert x == []

    @pytest.mark.parametrize("func", ["sum", "prod", "any", "all"])
    def test_apply_funcs_over_empty(self, func):
        # GH 28213
        df = DataFrame(columns=["a", "b", "c"])

        result = df.apply(getattr(np, func))
        expected = getattr(df, func)()
        tm.assert_series_equal(result, expected)

    def test_nunique_empty(self):
        # GH 28213
        df = DataFrame(columns=["a", "b", "c"])

        result = df.nunique()
        expected = Series(0, index=df.columns)
        tm.assert_series_equal(result, expected)

        result = df.T.nunique()
        expected = Series([], index=pd.Index([]), dtype=np.float64)
        tm.assert_series_equal(result, expected)

    def test_apply_standard_nonunique(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=["a", "a", "c"])

        result = df.apply(lambda s: s[0], axis=1)
        expected = Series([1, 4, 7], ["a", "a", "c"])
        tm.assert_series_equal(result, expected)

        result = df.T.apply(lambda s: s[0], axis=0)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("func", ["sum", "mean", "min", "max", "std"])
    @pytest.mark.parametrize(
        "args,kwds",
        [
            pytest.param([], {}, id="no_args_or_kwds"),
            pytest.param([1], {}, id="axis_from_args"),
            pytest.param([], {"axis": 1}, id="axis_from_kwds"),
            pytest.param([], {"numeric_only": True}, id="optional_kwds"),
            pytest.param([1, None], {"numeric_only": True}, id="args_and_kwds"),
        ],
    )
    def test_apply_with_string_funcs(self, float_frame, func, args, kwds):
        result = float_frame.apply(func, *args, **kwds)
        expected = getattr(float_frame, func)(*args, **kwds)
        tm.assert_series_equal(result, expected)

    def test_apply_broadcast(self, float_frame, int_frame_const_col):

        # scalars
        result = float_frame.apply(np.mean, result_type="broadcast")
        expected = DataFrame([float_frame.mean()], index=float_frame.index)
        tm.assert_frame_equal(result, expected)

        result = float_frame.apply(np.mean, axis=1, result_type="broadcast")
        m = float_frame.mean(axis=1)
        expected = DataFrame({c: m for c in float_frame.columns})
        tm.assert_frame_equal(result, expected)

        # lists
        result = float_frame.apply(
            lambda x: list(range(len(float_frame.columns))),
            axis=1,
            result_type="broadcast",
        )
        m = list(range(len(float_frame.columns)))
        expected = DataFrame(
            [m] * len(float_frame.index),
            dtype="float64",
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(result, expected)

        result = float_frame.apply(
            lambda x: list(range(len(float_frame.index))), result_type="broadcast"
        )
        m = list(range(len(float_frame.index)))
        expected = DataFrame(
            {c: m for c in float_frame.columns},
            dtype="float64",
            index=float_frame.index,
        )
        tm.assert_frame_equal(result, expected)

        # preserve columns
        df = int_frame_const_col
        result = df.apply(lambda x: [1, 2, 3], axis=1, result_type="broadcast")
        tm.assert_frame_equal(result, df)

        df = int_frame_const_col
        result = df.apply(
            lambda x: Series([1, 2, 3], index=list("abc")),
            axis=1,
            result_type="broadcast",
        )
        expected = df.copy()
        tm.assert_frame_equal(result, expected)

    def test_apply_broadcast_error(self, int_frame_const_col):
        df = int_frame_const_col

        # > 1 ndim
        with pytest.raises(ValueError):
            df.apply(
                lambda x: np.array([1, 2]).reshape(-1, 2),
                axis=1,
                result_type="broadcast",
            )

        # cannot broadcast
        with pytest.raises(ValueError):
            df.apply(lambda x: [1, 2], axis=1, result_type="broadcast")

        with pytest.raises(ValueError):
            df.apply(lambda x: Series([1, 2]), axis=1, result_type="broadcast")

    def test_apply_raw(self, float_frame):
        result0 = float_frame.apply(np.mean, raw=True)
        result1 = float_frame.apply(np.mean, axis=1, raw=True)

        expected0 = float_frame.apply(lambda x: x.values.mean())
        expected1 = float_frame.apply(lambda x: x.values.mean(), axis=1)

        tm.assert_series_equal(result0, expected0)
        tm.assert_series_equal(result1, expected1)

        # no reduction
        result = float_frame.apply(lambda x: x * 2, raw=True)
        expected = float_frame * 2
        tm.assert_frame_equal(result, expected)

    def test_apply_axis1(self, float_frame):
        d = float_frame.index[0]
        tapplied = float_frame.apply(np.mean, axis=1)
        assert tapplied[d] == np.mean(float_frame.xs(d))

    def test_apply_ignore_failures(self, float_string_frame):
        result = frame_apply(
            float_string_frame, np.mean, 0, ignore_failures=True
        ).apply_standard()
        expected = float_string_frame._get_numeric_data().apply(np.mean)
        tm.assert_series_equal(result, expected)

    def test_apply_mixed_dtype_corner(self):
        df = DataFrame({"A": ["foo"], "B": [1.0]})
        result = df[:0].apply(np.mean, axis=1)
        # the result here is actually kind of ambiguous, should it be a Series
        # or a DataFrame?
        expected = Series(np.nan, index=pd.Index([], dtype="int64"))
        tm.assert_series_equal(result, expected)

        df = DataFrame({"A": ["foo"], "B": [1.0]})
        result = df.apply(lambda x: x["A"], axis=1)
        expected = Series(["foo"], index=[0])
        tm.assert_series_equal(result, expected)

        result = df.apply(lambda x: x["B"], axis=1)
        expected = Series([1.0], index=[0])
        tm.assert_series_equal(result, expected)

    def test_apply_empty_infer_type(self):
        no_cols = DataFrame(index=["a", "b", "c"])
        no_index = DataFrame(columns=["a", "b", "c"])

        def _check(df, f):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("ignore", RuntimeWarning)
                test_res = f(np.array([], dtype="f8"))
            is_reduction = not isinstance(test_res, np.ndarray)

            def _checkit(axis=0, raw=False):
                result = df.apply(f, axis=axis, raw=raw)
                if is_reduction:
                    agg_axis = df._get_agg_axis(axis)
                    assert isinstance(result, Series)
                    assert result.index is agg_axis
                else:
                    assert isinstance(result, DataFrame)

            _checkit()
            _checkit(axis=1)
            _checkit(raw=True)
            _checkit(axis=0, raw=True)

        with np.errstate(all="ignore"):
            _check(no_cols, lambda x: x)
            _check(no_cols, lambda x: x.mean())
            _check(no_index, lambda x: x)
            _check(no_index, lambda x: x.mean())

        result = no_cols.apply(lambda x: x.mean(), result_type="broadcast")
        assert isinstance(result, DataFrame)

    def test_apply_with_args_kwds(self, float_frame):
        def add_some(x, howmuch=0):
            return x + howmuch

        def agg_and_add(x, howmuch=0):
            return x.mean() + howmuch

        def subtract_and_divide(x, sub, divide=1):
            return (x - sub) / divide

        result = float_frame.apply(add_some, howmuch=2)
        expected = float_frame.apply(lambda x: x + 2)
        tm.assert_frame_equal(result, expected)

        result = float_frame.apply(agg_and_add, howmuch=2)
        expected = float_frame.apply(lambda x: x.mean() + 2)
        tm.assert_series_equal(result, expected)

        result = float_frame.apply(subtract_and_divide, args=(2,), divide=2)
        expected = float_frame.apply(lambda x: (x - 2.0) / 2.0)
        tm.assert_frame_equal(result, expected)

    def test_apply_yield_list(self, float_frame):
        result = float_frame.apply(list)
        tm.assert_frame_equal(result, float_frame)

    def test_apply_reduce_Series(self, float_frame):
        float_frame.loc[::2, "A"] = np.nan
        expected = float_frame.mean(1)
        result = float_frame.apply(np.mean, axis=1)
        tm.assert_series_equal(result, expected)

    def test_apply_reduce_rows_to_dict(self):
        # GH 25196
        data = pd.DataFrame([[1, 2], [3, 4]])
        expected = pd.Series([{0: 1, 1: 3}, {0: 2, 1: 4}])
        result = data.apply(dict)
        tm.assert_series_equal(result, expected)

    def test_apply_differently_indexed(self):
        df = DataFrame(np.random.randn(20, 10))

        result0 = df.apply(Series.describe, axis=0)
        expected0 = DataFrame(
            {i: v.describe() for i, v in df.items()}, columns=df.columns
        )
        tm.assert_frame_equal(result0, expected0)

        result1 = df.apply(Series.describe, axis=1)
        expected1 = DataFrame(
            {i: v.describe() for i, v in df.T.items()}, columns=df.index
        ).T
        tm.assert_frame_equal(result1, expected1)

    def test_apply_modify_traceback(self):
        data = DataFrame(
            {
                "A": [
                    "foo",
                    "foo",
                    "foo",
                    "foo",
                    "bar",
                    "bar",
                    "bar",
                    "bar",
                    "foo",
                    "foo",
                    "foo",
                ],
                "B": [
                    "one",
                    "one",
                    "one",
                    "two",
                    "one",
                    "one",
                    "one",
                    "two",
                    "two",
                    "two",
                    "one",
                ],
                "C": [
                    "dull",
                    "dull",
                    "shiny",
                    "dull",
                    "dull",
                    "shiny",
                    "shiny",
                    "dull",
                    "shiny",
                    "shiny",
                    "shiny",
                ],
                "D": np.random.randn(11),
                "E": np.random.randn(11),
                "F": np.random.randn(11),
            }
        )

        data.loc[4, "C"] = np.nan

        def transform(row):
            if row["C"].startswith("shin") and row["A"] == "foo":
                row["D"] = 7
            return row

        def transform2(row):
            if notna(row["C"]) and row["C"].startswith("shin") and row["A"] == "foo":
                row["D"] = 7
            return row

        msg = "'float' object has no attribute 'startswith'"
        with pytest.raises(AttributeError, match=msg):
            data.apply(transform, axis=1)

    def test_apply_bug(self):

        # GH 6125
        positions = pd.DataFrame(
            [
                [1, "ABC0", 50],
                [1, "YUM0", 20],
                [1, "DEF0", 20],
                [2, "ABC1", 50],
                [2, "YUM1", 20],
                [2, "DEF1", 20],
            ],
            columns=["a", "market", "position"],
        )

        def f(r):
            return r["market"]

        expected = positions.apply(f, axis=1)

        positions = DataFrame(
            [
                [datetime(2013, 1, 1), "ABC0", 50],
                [datetime(2013, 1, 2), "YUM0", 20],
                [datetime(2013, 1, 3), "DEF0", 20],
                [datetime(2013, 1, 4), "ABC1", 50],
                [datetime(2013, 1, 5), "YUM1", 20],
                [datetime(2013, 1, 6), "DEF1", 20],
            ],
            columns=["a", "market", "position"],
        )
        result = positions.apply(f, axis=1)
        tm.assert_series_equal(result, expected)

    def test_apply_convert_objects(self):
        data = DataFrame(
            {
                "A": [
                    "foo",
                    "foo",
                    "foo",
                    "foo",
                    "bar",
                    "bar",
                    "bar",
                    "bar",
                    "foo",
                    "foo",
                    "foo",
                ],
                "B": [
                    "one",
                    "one",
                    "one",
                    "two",
                    "one",
                    "one",
                    "one",
                    "two",
                    "two",
                    "two",
                    "one",
                ],
                "C": [
                    "dull",
                    "dull",
                    "shiny",
                    "dull",
                    "dull",
                    "shiny",
                    "shiny",
                    "dull",
                    "shiny",
                    "shiny",
                    "shiny",
                ],
                "D": np.random.randn(11),
                "E": np.random.randn(11),
                "F": np.random.randn(11),
            }
        )

        result = data.apply(lambda x: x, axis=1)
        tm.assert_frame_equal(result._convert(datetime=True), data)

    def test_apply_attach_name(self, float_frame):
        result = float_frame.apply(lambda x: x.name)
        expected = Series(float_frame.columns, index=float_frame.columns)
        tm.assert_series_equal(result, expected)

        result = float_frame.apply(lambda x: x.name, axis=1)
        expected = Series(float_frame.index, index=float_frame.index)
        tm.assert_series_equal(result, expected)

        # non-reductions
        result = float_frame.apply(lambda x: np.repeat(x.name, len(x)))
        expected = DataFrame(
            np.tile(float_frame.columns, (len(float_frame.index), 1)),
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(result, expected)

        result = float_frame.apply(lambda x: np.repeat(x.name, len(x)), axis=1)
        expected = Series(
            np.repeat(t[0], len(float_frame.columns)) for t in float_frame.itertuples()
        )
        expected.index = float_frame.index
        tm.assert_series_equal(result, expected)

    def test_apply_multi_index(self, float_frame):
        index = MultiIndex.from_arrays([["a", "a", "b"], ["c", "d", "d"]])
        s = DataFrame([[1, 2], [3, 4], [5, 6]], index=index, columns=["col1", "col2"])
        result = s.apply(lambda x: Series({"min": min(x), "max": max(x)}), 1)
        expected = DataFrame(
            [[1, 2], [3, 4], [5, 6]], index=index, columns=["min", "max"]
        )
        tm.assert_frame_equal(result, expected, check_like=True)

    def test_apply_dict(self):

        # GH 8735
        A = DataFrame([["foo", "bar"], ["spam", "eggs"]])
        A_dicts = Series(
            [dict([(0, "foo"), (1, "spam")]), dict([(0, "bar"), (1, "eggs")])]
        )
        B = DataFrame([[0, 1], [2, 3]])
        B_dicts = Series([dict([(0, 0), (1, 2)]), dict([(0, 1), (1, 3)])])
        fn = lambda x: x.to_dict()

        for df, dicts in [(A, A_dicts), (B, B_dicts)]:
            reduce_true = df.apply(fn, result_type="reduce")
            reduce_false = df.apply(fn, result_type="expand")
            reduce_none = df.apply(fn)

            tm.assert_series_equal(reduce_true, dicts)
            tm.assert_frame_equal(reduce_false, df)
            tm.assert_series_equal(reduce_none, dicts)

    def test_applymap(self, float_frame):
        applied = float_frame.applymap(lambda x: x * 2)
        tm.assert_frame_equal(applied, float_frame * 2)
        float_frame.applymap(type)

        # GH 465: function returning tuples
        result = float_frame.applymap(lambda x: (x, x))
        assert isinstance(result["A"][0], tuple)

        # GH 2909: object conversion to float in constructor?
        df = DataFrame(data=[1, "a"])
        result = df.applymap(lambda x: x)
        assert result.dtypes[0] == object

        df = DataFrame(data=[1.0, "a"])
        result = df.applymap(lambda x: x)
        assert result.dtypes[0] == object

        # GH 2786
        df = DataFrame(np.random.random((3, 4)))
        df2 = df.copy()
        cols = ["a", "a", "a", "a"]
        df.columns = cols

        expected = df2.applymap(str)
        expected.columns = cols
        result = df.applymap(str)
        tm.assert_frame_equal(result, expected)

        # datetime/timedelta
        df["datetime"] = Timestamp("20130101")
        df["timedelta"] = pd.Timedelta("1 min")
        result = df.applymap(str)
        for f in ["datetime", "timedelta"]:
            assert result.loc[0, f] == str(df.loc[0, f])

        # GH 8222
        empty_frames = [
            pd.DataFrame(),
            pd.DataFrame(columns=list("ABC")),
            pd.DataFrame(index=list("ABC")),
            pd.DataFrame({"A": [], "B": [], "C": []}),
        ]
        for frame in empty_frames:
            for func in [round, lambda x: x]:
                result = frame.applymap(func)
                tm.assert_frame_equal(result, frame)

    def test_applymap_box_timestamps(self):
        # GH 2689, GH 2627
        ser = pd.Series(date_range("1/1/2000", periods=10))

        def func(x):
            return (x.hour, x.day, x.month)

        # it works!
        pd.DataFrame(ser).applymap(func)

    def test_applymap_box(self):
        # ufunc will not be boxed. Same test cases as the test_map_box
        df = pd.DataFrame(
            {
                "a": [pd.Timestamp("2011-01-01"), pd.Timestamp("2011-01-02")],
                "b": [
                    pd.Timestamp("2011-01-01", tz="US/Eastern"),
                    pd.Timestamp("2011-01-02", tz="US/Eastern"),
                ],
                "c": [pd.Timedelta("1 days"), pd.Timedelta("2 days")],
                "d": [
                    pd.Period("2011-01-01", freq="M"),
                    pd.Period("2011-01-02", freq="M"),
                ],
            }
        )

        result = df.applymap(lambda x: type(x).__name__)
        expected = pd.DataFrame(
            {
                "a": ["Timestamp", "Timestamp"],
                "b": ["Timestamp", "Timestamp"],
                "c": ["Timedelta", "Timedelta"],
                "d": ["Period", "Period"],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_frame_apply_dont_convert_datetime64(self):
        from pandas.tseries.offsets import BDay

        df = DataFrame({"x1": [datetime(1996, 1, 1)]})

        df = df.applymap(lambda x: x + BDay())
        df = df.applymap(lambda x: x + BDay())

        assert df.x1.dtype == "M8[ns]"

    def test_apply_non_numpy_dtype(self):
        # GH 12244
        df = DataFrame(
            {"dt": pd.date_range("2015-01-01", periods=3, tz="Europe/Brussels")}
        )
        result = df.apply(lambda x: x)
        tm.assert_frame_equal(result, df)

        result = df.apply(lambda x: x + pd.Timedelta("1day"))
        expected = DataFrame(
            {"dt": pd.date_range("2015-01-02", periods=3, tz="Europe/Brussels")}
        )
        tm.assert_frame_equal(result, expected)

        df = DataFrame({"dt": ["a", "b", "c", "a"]}, dtype="category")
        result = df.apply(lambda x: x)
        tm.assert_frame_equal(result, df)

    def test_apply_dup_names_multi_agg(self):
        # GH 21063
        df = pd.DataFrame([[0, 1], [2, 3]], columns=["a", "a"])
        expected = pd.DataFrame([[0, 1]], columns=["a", "a"], index=["min"])
        result = df.agg(["min"])

        tm.assert_frame_equal(result, expected)

    def test_apply_nested_result_axis_1(self):
        # GH 13820
        def apply_list(row):
            return [2 * row["A"], 2 * row["C"], 2 * row["B"]]

        df = pd.DataFrame(np.zeros((4, 4)), columns=list("ABCD"))
        result = df.apply(apply_list, axis=1)
        expected = Series(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )
        tm.assert_series_equal(result, expected)

    def test_apply_noreduction_tzaware_object(self):
        # https://github.com/pandas-dev/pandas/issues/31505
        df = pd.DataFrame({"foo": [pd.Timestamp("2020", tz="UTC")]}, dtype="object")
        result = df.apply(lambda x: x)
        tm.assert_frame_equal(result, df)
        result = df.apply(lambda x: x.copy())
        tm.assert_frame_equal(result, df)


class TestInferOutputShape:
    # the user has supplied an opaque UDF where
    # they are transforming the input that requires
    # us to infer the output

    def test_infer_row_shape(self):
        # GH 17437
        # if row shape is changing, infer it
        df = pd.DataFrame(np.random.rand(10, 2))
        result = df.apply(np.fft.fft, axis=0)
        assert result.shape == (10, 2)

        result = df.apply(np.fft.rfft, axis=0)
        assert result.shape == (6, 2)

    def test_with_dictlike_columns(self):
        # GH 17602
        df = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
        result = df.apply(lambda x: {"s": x["a"] + x["b"]}, axis=1)
        expected = Series([{"s": 3} for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

        df["tm"] = [
            pd.Timestamp("2017-05-01 00:00:00"),
            pd.Timestamp("2017-05-02 00:00:00"),
        ]
        result = df.apply(lambda x: {"s": x["a"] + x["b"]}, axis=1)
        tm.assert_series_equal(result, expected)

        # compose a series
        result = (df["a"] + df["b"]).apply(lambda x: {"s": x})
        expected = Series([{"s": 3}, {"s": 3}])
        tm.assert_series_equal(result, expected)

        # GH 18775
        df = DataFrame()
        df["author"] = ["X", "Y", "Z"]
        df["publisher"] = ["BBC", "NBC", "N24"]
        df["date"] = pd.to_datetime(
            ["17-10-2010 07:15:30", "13-05-2011 08:20:35", "15-01-2013 09:09:09"]
        )
        result = df.apply(lambda x: {}, axis=1)
        expected = Series([{}, {}, {}])
        tm.assert_series_equal(result, expected)

    def test_with_dictlike_columns_with_infer(self):
        # GH 17602
        df = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
        result = df.apply(
            lambda x: {"s": x["a"] + x["b"]}, axis=1, result_type="expand"
        )
        expected = DataFrame({"s": [3, 3]})
        tm.assert_frame_equal(result, expected)

        df["tm"] = [
            pd.Timestamp("2017-05-01 00:00:00"),
            pd.Timestamp("2017-05-02 00:00:00"),
        ]
        result = df.apply(
            lambda x: {"s": x["a"] + x["b"]}, axis=1, result_type="expand"
        )
        tm.assert_frame_equal(result, expected)

    def test_with_listlike_columns(self):
        # GH 17348
        df = DataFrame(
            {
                "a": Series(np.random.randn(4)),
                "b": ["a", "list", "of", "words"],
                "ts": date_range("2016-10-01", periods=4, freq="H"),
            }
        )

        result = df[["a", "b"]].apply(tuple, axis=1)
        expected = Series([t[1:] for t in df[["a", "b"]].itertuples()])
        tm.assert_series_equal(result, expected)

        result = df[["a", "ts"]].apply(tuple, axis=1)
        expected = Series([t[1:] for t in df[["a", "ts"]].itertuples()])
        tm.assert_series_equal(result, expected)

        # GH 18919
        df = DataFrame(
            {"x": Series([["a", "b"], ["q"]]), "y": Series([["z"], ["q", "t"]])}
        )
        df.index = MultiIndex.from_tuples([("i0", "j0"), ("i1", "j1")])

        result = df.apply(lambda row: [el for el in row["x"] if el in row["y"]], axis=1)
        expected = Series([[], ["q"]], index=df.index)
        tm.assert_series_equal(result, expected)

    def test_infer_output_shape_columns(self):
        # GH 18573

        df = DataFrame(
            {
                "number": [1.0, 2.0],
                "string": ["foo", "bar"],
                "datetime": [
                    pd.Timestamp("2017-11-29 03:30:00"),
                    pd.Timestamp("2017-11-29 03:45:00"),
                ],
            }
        )
        result = df.apply(lambda row: (row.number, row.string), axis=1)
        expected = Series([(t.number, t.string) for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

    def test_infer_output_shape_listlike_columns(self):
        # GH 16353

        df = DataFrame(np.random.randn(6, 3), columns=["A", "B", "C"])

        result = df.apply(lambda x: [1, 2, 3], axis=1)
        expected = Series([[1, 2, 3] for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

        result = df.apply(lambda x: [1, 2], axis=1)
        expected = Series([[1, 2] for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

        # GH 17970
        df = DataFrame({"a": [1, 2, 3]}, index=list("abc"))

        result = df.apply(lambda row: np.ones(1), axis=1)
        expected = Series([np.ones(1) for t in df.itertuples()], index=df.index)
        tm.assert_series_equal(result, expected)

        result = df.apply(lambda row: np.ones(2), axis=1)
        expected = Series([np.ones(2) for t in df.itertuples()], index=df.index)
        tm.assert_series_equal(result, expected)

        # GH 17892
        df = pd.DataFrame(
            {
                "a": [
                    pd.Timestamp("2010-02-01"),
                    pd.Timestamp("2010-02-04"),
                    pd.Timestamp("2010-02-05"),
                    pd.Timestamp("2010-02-06"),
                ],
                "b": [9, 5, 4, 3],
                "c": [5, 3, 4, 2],
                "d": [1, 2, 3, 4],
            }
        )

        def fun(x):
            return (1, 2)

        result = df.apply(fun, axis=1)
        expected = Series([(1, 2) for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

    def test_consistent_coerce_for_shapes(self):
        # we want column names to NOT be propagated
        # just because the shape matches the input shape
        df = DataFrame(np.random.randn(4, 3), columns=["A", "B", "C"])

        result = df.apply(lambda x: [1, 2, 3], axis=1)
        expected = Series([[1, 2, 3] for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

        result = df.apply(lambda x: [1, 2], axis=1)
        expected = Series([[1, 2] for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

    def test_consistent_names(self, int_frame_const_col):
        # if a Series is returned, we should use the resulting index names
        df = int_frame_const_col

        result = df.apply(
            lambda x: Series([1, 2, 3], index=["test", "other", "cols"]), axis=1
        )
        expected = int_frame_const_col.rename(
            columns={"A": "test", "B": "other", "C": "cols"}
        )
        tm.assert_frame_equal(result, expected)

        result = df.apply(lambda x: Series([1, 2], index=["test", "other"]), axis=1)
        expected = expected[["test", "other"]]
        tm.assert_frame_equal(result, expected)

    def test_result_type(self, int_frame_const_col):
        # result_type should be consistent no matter which
        # path we take in the code
        df = int_frame_const_col

        result = df.apply(lambda x: [1, 2, 3], axis=1, result_type="expand")
        expected = df.copy()
        expected.columns = [0, 1, 2]
        tm.assert_frame_equal(result, expected)

        result = df.apply(lambda x: [1, 2], axis=1, result_type="expand")
        expected = df[["A", "B"]].copy()
        expected.columns = [0, 1]
        tm.assert_frame_equal(result, expected)

        # broadcast result
        result = df.apply(lambda x: [1, 2, 3], axis=1, result_type="broadcast")
        expected = df.copy()
        tm.assert_frame_equal(result, expected)

        columns = ["other", "col", "names"]
        result = df.apply(
            lambda x: Series([1, 2, 3], index=columns), axis=1, result_type="broadcast"
        )
        expected = df.copy()
        tm.assert_frame_equal(result, expected)

        # series result
        result = df.apply(lambda x: Series([1, 2, 3], index=x.index), axis=1)
        expected = df.copy()
        tm.assert_frame_equal(result, expected)

        # series result with other index
        columns = ["other", "col", "names"]
        result = df.apply(lambda x: Series([1, 2, 3], index=columns), axis=1)
        expected = df.copy()
        expected.columns = columns
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("result_type", ["foo", 1])
    def test_result_type_error(self, result_type, int_frame_const_col):
        # allowed result_type
        df = int_frame_const_col

        with pytest.raises(ValueError):
            df.apply(lambda x: [1, 2, 3], axis=1, result_type=result_type)

    @pytest.mark.parametrize(
        "box",
        [lambda x: list(x), lambda x: tuple(x), lambda x: np.array(x, dtype="int64")],
        ids=["list", "tuple", "array"],
    )
    def test_consistency_for_boxed(self, box, int_frame_const_col):
        # passing an array or list should not affect the output shape
        df = int_frame_const_col

        result = df.apply(lambda x: box([1, 2]), axis=1)
        expected = Series([box([1, 2]) for t in df.itertuples()])
        tm.assert_series_equal(result, expected)

        result = df.apply(lambda x: box([1, 2]), axis=1, result_type="expand")
        expected = int_frame_const_col[["A", "B"]].rename(columns={"A": 0, "B": 1})
        tm.assert_frame_equal(result, expected)


def zip_frames(frames, axis=1):
    """
    take a list of frames, zip them together under the
    assumption that these all have the first frames' index/columns.

    Returns
    -------
    new_frame : DataFrame
    """
    if axis == 1:
        columns = frames[0].columns
        zipped = [f.loc[:, c] for c in columns for f in frames]
        return pd.concat(zipped, axis=1)
    else:
        index = frames[0].index
        zipped = [f.loc[i, :] for i in index for f in frames]
        return pd.DataFrame(zipped)


class TestDataFrameAggregate:
    def test_agg_transform(self, axis, float_frame):
        other_axis = 1 if axis in {0, "index"} else 0

        with np.errstate(all="ignore"):

            f_abs = np.abs(float_frame)
            f_sqrt = np.sqrt(float_frame)

            # ufunc
            result = float_frame.transform(np.sqrt, axis=axis)
            expected = f_sqrt.copy()
            tm.assert_frame_equal(result, expected)

            result = float_frame.apply(np.sqrt, axis=axis)
            tm.assert_frame_equal(result, expected)

            result = float_frame.transform(np.sqrt, axis=axis)
            tm.assert_frame_equal(result, expected)

            # list-like
            result = float_frame.apply([np.sqrt], axis=axis)
            expected = f_sqrt.copy()
            if axis in {0, "index"}:
                expected.columns = pd.MultiIndex.from_product(
                    [float_frame.columns, ["sqrt"]]
                )
            else:
                expected.index = pd.MultiIndex.from_product(
                    [float_frame.index, ["sqrt"]]
                )
            tm.assert_frame_equal(result, expected)

            result = float_frame.transform([np.sqrt], axis=axis)
            tm.assert_frame_equal(result, expected)

            # multiple items in list
            # these are in the order as if we are applying both
            # functions per series and then concatting
            result = float_frame.apply([np.abs, np.sqrt], axis=axis)
            expected = zip_frames([f_abs, f_sqrt], axis=other_axis)
            if axis in {0, "index"}:
                expected.columns = pd.MultiIndex.from_product(
                    [float_frame.columns, ["absolute", "sqrt"]]
                )
            else:
                expected.index = pd.MultiIndex.from_product(
                    [float_frame.index, ["absolute", "sqrt"]]
                )
            tm.assert_frame_equal(result, expected)

            result = float_frame.transform([np.abs, "sqrt"], axis=axis)
            tm.assert_frame_equal(result, expected)

    def test_transform_and_agg_err(self, axis, float_frame):
        # cannot both transform and agg
        with pytest.raises(ValueError):
            float_frame.transform(["max", "min"], axis=axis)

        with pytest.raises(ValueError):
            with np.errstate(all="ignore"):
                float_frame.agg(["max", "sqrt"], axis=axis)

        with pytest.raises(ValueError):
            with np.errstate(all="ignore"):
                float_frame.transform(["max", "sqrt"], axis=axis)

        df = pd.DataFrame({"A": range(5), "B": 5})

        def f():
            with np.errstate(all="ignore"):
                df.agg({"A": ["abs", "sum"], "B": ["mean", "max"]}, axis=axis)

    @pytest.mark.parametrize("method", ["abs", "shift", "pct_change", "cumsum", "rank"])
    def test_transform_method_name(self, method):
        # GH 19760
        df = pd.DataFrame({"A": [-1, 2]})
        result = df.transform(method)
        expected = operator.methodcaller(method)(df)
        tm.assert_frame_equal(result, expected)

    def test_demo(self):
        # demonstration tests
        df = pd.DataFrame({"A": range(5), "B": 5})

        result = df.agg(["min", "max"])
        expected = DataFrame(
            {"A": [0, 4], "B": [5, 5]}, columns=["A", "B"], index=["min", "max"]
        )
        tm.assert_frame_equal(result, expected)

        result = df.agg({"A": ["min", "max"], "B": ["sum", "max"]})
        expected = DataFrame(
            {"A": [4.0, 0.0, np.nan], "B": [5.0, np.nan, 25.0]},
            columns=["A", "B"],
            index=["max", "min", "sum"],
        )
        tm.assert_frame_equal(result.reindex_like(expected), expected)

    def test_agg_multiple_mixed_no_warning(self):
        # GH 20909
        mdf = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": [1.0, 2.0, 3.0],
                "C": ["foo", "bar", "baz"],
                "D": pd.date_range("20130101", periods=3),
            }
        )
        expected = pd.DataFrame(
            {
                "A": [1, 6],
                "B": [1.0, 6.0],
                "C": ["bar", "foobarbaz"],
                "D": [pd.Timestamp("2013-01-01"), pd.NaT],
            },
            index=["min", "sum"],
        )
        # sorted index
        with tm.assert_produces_warning(None):
            result = mdf.agg(["min", "sum"])

        tm.assert_frame_equal(result, expected)

        with tm.assert_produces_warning(None):
            result = mdf[["D", "C", "B", "A"]].agg(["sum", "min"])

        # For backwards compatibility, the result's index is
        # still sorted by function name, so it's ['min', 'sum']
        # not ['sum', 'min'].
        expected = expected[["D", "C", "B", "A"]]
        tm.assert_frame_equal(result, expected)

    def test_agg_dict_nested_renaming_depr(self):

        df = pd.DataFrame({"A": range(5), "B": 5})

        # nested renaming
        msg = r"nested renamer is not supported"
        with pytest.raises(SpecificationError, match=msg):
            df.agg({"A": {"foo": "min"}, "B": {"bar": "max"}})

    def test_agg_reduce(self, axis, float_frame):
        other_axis = 1 if axis in {0, "index"} else 0
        name1, name2 = float_frame.axes[other_axis].unique()[:2].sort_values()

        # all reducers
        expected = pd.concat(
            [
                float_frame.mean(axis=axis),
                float_frame.max(axis=axis),
                float_frame.sum(axis=axis),
            ],
            axis=1,
        )
        expected.columns = ["mean", "max", "sum"]
        expected = expected.T if axis in {0, "index"} else expected

        result = float_frame.agg(["mean", "max", "sum"], axis=axis)
        tm.assert_frame_equal(result, expected)

        # dict input with scalars
        func = OrderedDict([(name1, "mean"), (name2, "sum")])
        result = float_frame.agg(func, axis=axis)
        expected = Series(
            [
                float_frame.loc(other_axis)[name1].mean(),
                float_frame.loc(other_axis)[name2].sum(),
            ],
            index=[name1, name2],
        )
        tm.assert_series_equal(result, expected)

        # dict input with lists
        func = OrderedDict([(name1, ["mean"]), (name2, ["sum"])])
        result = float_frame.agg(func, axis=axis)
        expected = DataFrame(
            {
                name1: Series(
                    [float_frame.loc(other_axis)[name1].mean()], index=["mean"]
                ),
                name2: Series(
                    [float_frame.loc(other_axis)[name2].sum()], index=["sum"]
                ),
            }
        )
        expected = expected.T if axis in {1, "columns"} else expected
        tm.assert_frame_equal(result, expected)

        # dict input with lists with multiple
        func = OrderedDict([(name1, ["mean", "sum"]), (name2, ["sum", "max"])])
        result = float_frame.agg(func, axis=axis)
        expected = DataFrame(
            OrderedDict(
                [
                    (
                        name1,
                        Series(
                            [
                                float_frame.loc(other_axis)[name1].mean(),
                                float_frame.loc(other_axis)[name1].sum(),
                            ],
                            index=["mean", "sum"],
                        ),
                    ),
                    (
                        name2,
                        Series(
                            [
                                float_frame.loc(other_axis)[name2].sum(),
                                float_frame.loc(other_axis)[name2].max(),
                            ],
                            index=["sum", "max"],
                        ),
                    ),
                ]
            )
        )
        expected = expected.T if axis in {1, "columns"} else expected
        tm.assert_frame_equal(result, expected)

    def test_nuiscance_columns(self):

        # GH 15015
        df = DataFrame(
            {
                "A": [1, 2, 3],
                "B": [1.0, 2.0, 3.0],
                "C": ["foo", "bar", "baz"],
                "D": pd.date_range("20130101", periods=3),
            }
        )

        result = df.agg("min")
        expected = Series([1, 1.0, "bar", pd.Timestamp("20130101")], index=df.columns)
        tm.assert_series_equal(result, expected)

        result = df.agg(["min"])
        expected = DataFrame(
            [[1, 1.0, "bar", pd.Timestamp("20130101")]],
            index=["min"],
            columns=df.columns,
        )
        tm.assert_frame_equal(result, expected)

        result = df.agg("sum")
        expected = Series([6, 6.0, "foobarbaz"], index=["A", "B", "C"])
        tm.assert_series_equal(result, expected)

        result = df.agg(["sum"])
        expected = DataFrame(
            [[6, 6.0, "foobarbaz"]], index=["sum"], columns=["A", "B", "C"]
        )
        tm.assert_frame_equal(result, expected)

    def test_non_callable_aggregates(self):

        # GH 16405
        # 'size' is a property of frame/series
        # validate that this is working
        df = DataFrame(
            {"A": [None, 2, 3], "B": [1.0, np.nan, 3.0], "C": ["foo", None, "bar"]}
        )

        # Function aggregate
        result = df.agg({"A": "count"})
        expected = Series({"A": 2})

        tm.assert_series_equal(result, expected)

        # Non-function aggregate
        result = df.agg({"A": "size"})
        expected = Series({"A": 3})

        tm.assert_series_equal(result, expected)

        # Mix function and non-function aggs
        result1 = df.agg(["count", "size"])
        result2 = df.agg(
            {"A": ["count", "size"], "B": ["count", "size"], "C": ["count", "size"]}
        )
        expected = pd.DataFrame(
            {
                "A": {"count": 2, "size": 3},
                "B": {"count": 2, "size": 3},
                "C": {"count": 2, "size": 3},
            }
        )

        tm.assert_frame_equal(result1, result2, check_like=True)
        tm.assert_frame_equal(result2, expected, check_like=True)

        # Just functional string arg is same as calling df.arg()
        result = df.agg("count")
        expected = df.count()

        tm.assert_series_equal(result, expected)

        # Just a string attribute arg same as calling df.arg
        result = df.agg("size")
        expected = df.size

        assert result == expected

    def test_agg_listlike_result(self):
        # GH-29587 user defined function returning list-likes
        df = DataFrame(
            {"A": [2, 2, 3], "B": [1.5, np.nan, 1.5], "C": ["foo", None, "bar"]}
        )

        def func(group_col):
            return list(group_col.dropna().unique())

        result = df.agg(func)
        expected = pd.Series([[2, 3], [1.5], ["foo", "bar"]], index=["A", "B", "C"])
        tm.assert_series_equal(result, expected)

        result = df.agg([func])
        expected = expected.to_frame("func").T
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "df, func, expected",
        chain(
            _get_cython_table_params(
                DataFrame(),
                [
                    ("sum", Series(dtype="float64")),
                    ("max", Series(dtype="float64")),
                    ("min", Series(dtype="float64")),
                    ("all", Series(dtype=bool)),
                    ("any", Series(dtype=bool)),
                    ("mean", Series(dtype="float64")),
                    ("prod", Series(dtype="float64")),
                    ("std", Series(dtype="float64")),
                    ("var", Series(dtype="float64")),
                    ("median", Series(dtype="float64")),
                ],
            ),
            _get_cython_table_params(
                DataFrame([[np.nan, 1], [1, 2]]),
                [
                    ("sum", Series([1.0, 3])),
                    ("max", Series([1.0, 2])),
                    ("min", Series([1.0, 1])),
                    ("all", Series([True, True])),
                    ("any", Series([True, True])),
                    ("mean", Series([1, 1.5])),
                    ("prod", Series([1.0, 2])),
                    ("std", Series([np.nan, 0.707107])),
                    ("var", Series([np.nan, 0.5])),
                    ("median", Series([1, 1.5])),
                ],
            ),
        ),
    )
    def test_agg_cython_table(self, df, func, expected, axis):
        # GH 21224
        # test reducing functions in
        # pandas.core.base.SelectionMixin._cython_table
        result = df.agg(func, axis=axis)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "df, func, expected",
        chain(
            _get_cython_table_params(
                DataFrame(), [("cumprod", DataFrame()), ("cumsum", DataFrame())]
            ),
            _get_cython_table_params(
                DataFrame([[np.nan, 1], [1, 2]]),
                [
                    ("cumprod", DataFrame([[np.nan, 1], [1, 2]])),
                    ("cumsum", DataFrame([[np.nan, 1], [1, 3]])),
                ],
            ),
        ),
    )
    def test_agg_cython_table_transform(self, df, func, expected, axis):
        # GH 21224
        # test transforming functions in
        # pandas.core.base.SelectionMixin._cython_table (cumprod, cumsum)
        if axis == "columns" or axis == 1:
            # operating blockwise doesn't let us preserve dtypes
            expected = expected.astype("float64")

        result = df.agg(func, axis=axis)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "df, func, expected",
        _get_cython_table_params(
            DataFrame([["a", "b"], ["b", "a"]]), [["cumprod", TypeError]]
        ),
    )
    def test_agg_cython_table_raises(self, df, func, expected, axis):
        # GH 21224
        with pytest.raises(expected):
            df.agg(func, axis=axis)

    @pytest.mark.parametrize("num_cols", [2, 3, 5])
    def test_frequency_is_original(self, num_cols):
        # GH 22150
        index = pd.DatetimeIndex(["1950-06-30", "1952-10-24", "1953-05-29"])
        original = index.copy()
        df = DataFrame(1, index=index, columns=range(num_cols))
        df.apply(lambda x: x)
        assert index.freq == original.freq

    def test_apply_datetime_tz_issue(self):
        # GH 29052

        timestamps = [
            pd.Timestamp("2019-03-15 12:34:31.909000+0000", tz="UTC"),
            pd.Timestamp("2019-03-15 12:34:34.359000+0000", tz="UTC"),
            pd.Timestamp("2019-03-15 12:34:34.660000+0000", tz="UTC"),
        ]
        df = DataFrame(data=[0, 1, 2], index=timestamps)
        result = df.apply(lambda x: x.name, axis=1)
        expected = pd.Series(index=timestamps, data=timestamps)

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("df", [pd.DataFrame({"A": ["a", None], "B": ["c", "d"]})])
    @pytest.mark.parametrize("method", ["min", "max", "sum"])
    def test_consistency_of_aggregates_of_columns_with_missing_values(self, df, method):
        # GH 16832
        none_in_first_column_result = getattr(df[["A", "B"]], method)()
        none_in_second_column_result = getattr(df[["B", "A"]], method)()

        tm.assert_series_equal(
            none_in_first_column_result, none_in_second_column_result
        )
