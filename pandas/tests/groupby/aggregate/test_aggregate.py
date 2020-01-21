"""
test .agg behavior / note that .apply is tested generally in test_groupby.py
"""
import functools

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series, concat
import pandas._testing as tm
from pandas.core.base import SpecificationError
from pandas.core.groupby.grouper import Grouping


def test_agg_regression1(tsframe):
    grouped = tsframe.groupby([lambda x: x.year, lambda x: x.month])
    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)


def test_agg_must_agg(df):
    grouped = df.groupby("A")["C"]

    msg = "Must produce aggregated value"
    with pytest.raises(Exception, match=msg):
        grouped.agg(lambda x: x.describe())
    with pytest.raises(Exception, match=msg):
        grouped.agg(lambda x: x.index[:2])


def test_agg_ser_multi_key(df):
    # TODO(wesm): unused
    ser = df.C  # noqa

    f = lambda x: x.sum()
    results = df.C.groupby([df.A, df.B]).aggregate(f)
    expected = df.groupby(["A", "B"]).sum()["C"]
    tm.assert_series_equal(results, expected)


def test_groupby_aggregation_mixed_dtype():

    # GH 6212
    expected = DataFrame(
        {
            "v1": [5, 5, 7, np.nan, 3, 3, 4, 1],
            "v2": [55, 55, 77, np.nan, 33, 33, 44, 11],
        },
        index=MultiIndex.from_tuples(
            [
                (1, 95),
                (1, 99),
                (2, 95),
                (2, 99),
                ("big", "damp"),
                ("blue", "dry"),
                ("red", "red"),
                ("red", "wet"),
            ],
            names=["by1", "by2"],
        ),
    )

    df = DataFrame(
        {
            "v1": [1, 3, 5, 7, 8, 3, 5, np.nan, 4, 5, 7, 9],
            "v2": [11, 33, 55, 77, 88, 33, 55, np.nan, 44, 55, 77, 99],
            "by1": ["red", "blue", 1, 2, np.nan, "big", 1, 2, "red", 1, np.nan, 12],
            "by2": [
                "wet",
                "dry",
                99,
                95,
                np.nan,
                "damp",
                95,
                99,
                "red",
                99,
                np.nan,
                np.nan,
            ],
        }
    )

    g = df.groupby(["by1", "by2"])
    result = g[["v1", "v2"]].mean()
    tm.assert_frame_equal(result, expected)


def test_groupby_aggregation_multi_level_column():
    # GH 29772
    lst = [
        [True, True, True, False],
        [True, False, np.nan, False],
        [True, True, np.nan, False],
        [True, True, np.nan, False],
    ]
    df = pd.DataFrame(
        data=lst,
        columns=pd.MultiIndex.from_tuples([("A", 0), ("A", 1), ("B", 0), ("B", 1)]),
    )

    result = df.groupby(level=1, axis=1).sum()
    expected = pd.DataFrame({0: [2.0, 1, 1, 1], 1: [1, 0, 1, 1]})

    tm.assert_frame_equal(result, expected)


def test_agg_apply_corner(ts, tsframe):
    # nothing to group, all NA
    grouped = ts.groupby(ts * np.nan)
    assert ts.dtype == np.float64

    # groupby float64 values results in Float64Index
    exp = Series([], dtype=np.float64, index=pd.Index([], dtype=np.float64))
    tm.assert_series_equal(grouped.sum(), exp)
    tm.assert_series_equal(grouped.agg(np.sum), exp)
    tm.assert_series_equal(grouped.apply(np.sum), exp, check_index_type=False)

    # DataFrame
    grouped = tsframe.groupby(tsframe["A"] * np.nan)
    exp_df = DataFrame(
        columns=tsframe.columns, dtype=float, index=pd.Index([], dtype=np.float64)
    )
    tm.assert_frame_equal(grouped.sum(), exp_df, check_names=False)
    tm.assert_frame_equal(grouped.agg(np.sum), exp_df, check_names=False)
    tm.assert_frame_equal(grouped.apply(np.sum), exp_df.iloc[:, :0], check_names=False)


def test_agg_grouping_is_list_tuple(ts):
    df = tm.makeTimeDataFrame()

    grouped = df.groupby(lambda x: x.year)
    grouper = grouped.grouper.groupings[0].grouper
    grouped.grouper.groupings[0] = Grouping(ts.index, list(grouper))

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)

    grouped.grouper.groupings[0] = Grouping(ts.index, tuple(grouper))

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)


def test_agg_python_multiindex(mframe):
    grouped = mframe.groupby(["A", "B"])

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "groupbyfunc", [lambda x: x.weekday(), [lambda x: x.month, lambda x: x.weekday()]]
)
def test_aggregate_str_func(tsframe, groupbyfunc):
    grouped = tsframe.groupby(groupbyfunc)

    # single series
    result = grouped["A"].agg("std")
    expected = grouped["A"].std()
    tm.assert_series_equal(result, expected)

    # group frame by function name
    result = grouped.aggregate("var")
    expected = grouped.var()
    tm.assert_frame_equal(result, expected)

    # group frame by function dict
    result = grouped.agg({"A": "var", "B": "std", "C": "mean", "D": "sem"})
    expected = DataFrame(
        {
            "A": grouped["A"].var(),
            "B": grouped["B"].std(),
            "C": grouped["C"].mean(),
            "D": grouped["D"].sem(),
        }
    )
    tm.assert_frame_equal(result, expected)


def test_aggregate_item_by_item(df):
    grouped = df.groupby("A")

    aggfun = lambda ser: ser.size
    result = grouped.agg(aggfun)
    foo = (df.A == "foo").sum()
    bar = (df.A == "bar").sum()
    K = len(result.columns)

    # GH5782
    # odd comparisons can result here, so cast to make easy
    exp = pd.Series(
        np.array([foo] * K), index=list("BCD"), dtype=np.float64, name="foo"
    )
    tm.assert_series_equal(result.xs("foo"), exp)

    exp = pd.Series(
        np.array([bar] * K), index=list("BCD"), dtype=np.float64, name="bar"
    )
    tm.assert_almost_equal(result.xs("bar"), exp)

    def aggfun(ser):
        return ser.size

    result = DataFrame().groupby(df.A).agg(aggfun)
    assert isinstance(result, DataFrame)
    assert len(result) == 0


def test_wrap_agg_out(three_group):
    grouped = three_group.groupby(["A", "B"])

    def func(ser):
        if ser.dtype == np.object:
            raise TypeError
        else:
            return ser.sum()

    result = grouped.aggregate(func)
    exp_grouped = three_group.loc[:, three_group.columns != "C"]
    expected = exp_grouped.groupby(["A", "B"]).aggregate(func)
    tm.assert_frame_equal(result, expected)


def test_agg_multiple_functions_maintain_order(df):
    # GH #610
    funcs = [("mean", np.mean), ("max", np.max), ("min", np.min)]
    result = df.groupby("A")["C"].agg(funcs)
    exp_cols = Index(["mean", "max", "min"])

    tm.assert_index_equal(result.columns, exp_cols)


def test_multiple_functions_tuples_and_non_tuples(df):
    # #1359
    funcs = [("foo", "mean"), "std"]
    ex_funcs = [("foo", "mean"), ("std", "std")]

    result = df.groupby("A")["C"].agg(funcs)
    expected = df.groupby("A")["C"].agg(ex_funcs)
    tm.assert_frame_equal(result, expected)

    result = df.groupby("A").agg(funcs)
    expected = df.groupby("A").agg(ex_funcs)
    tm.assert_frame_equal(result, expected)


def test_more_flexible_frame_multi_function(df):
    grouped = df.groupby("A")

    exmean = grouped.agg({"C": np.mean, "D": np.mean})
    exstd = grouped.agg({"C": np.std, "D": np.std})

    expected = concat([exmean, exstd], keys=["mean", "std"], axis=1)
    expected = expected.swaplevel(0, 1, axis=1).sort_index(level=0, axis=1)

    d = {"C": [np.mean, np.std], "D": [np.mean, np.std]}
    result = grouped.aggregate(d)

    tm.assert_frame_equal(result, expected)

    # be careful
    result = grouped.aggregate({"C": np.mean, "D": [np.mean, np.std]})
    expected = grouped.aggregate({"C": np.mean, "D": [np.mean, np.std]})
    tm.assert_frame_equal(result, expected)

    def foo(x):
        return np.mean(x)

    def bar(x):
        return np.std(x, ddof=1)

    # this uses column selection & renaming
    msg = r"nested renamer is not supported"
    with pytest.raises(SpecificationError, match=msg):
        d = dict([["C", np.mean], ["D", dict([["foo", np.mean], ["bar", np.std]])]])
        grouped.aggregate(d)

    # But without renaming, these functions are OK
    d = {"C": [np.mean], "D": [foo, bar]}
    grouped.aggregate(d)


def test_multi_function_flexible_mix(df):
    # GH #1268
    grouped = df.groupby("A")

    # Expected
    d = {"C": {"foo": "mean", "bar": "std"}, "D": {"sum": "sum"}}
    # this uses column selection & renaming
    msg = r"nested renamer is not supported"
    with pytest.raises(SpecificationError, match=msg):
        grouped.aggregate(d)

    # Test 1
    d = {"C": {"foo": "mean", "bar": "std"}, "D": "sum"}
    # this uses column selection & renaming
    with pytest.raises(SpecificationError, match=msg):
        grouped.aggregate(d)

    # Test 2
    d = {"C": {"foo": "mean", "bar": "std"}, "D": "sum"}
    # this uses column selection & renaming
    with pytest.raises(SpecificationError, match=msg):
        grouped.aggregate(d)


def test_groupby_agg_coercing_bools():
    # issue 14873
    dat = pd.DataFrame({"a": [1, 1, 2, 2], "b": [0, 1, 2, 3], "c": [None, None, 1, 1]})
    gp = dat.groupby("a")

    index = Index([1, 2], name="a")

    result = gp["b"].aggregate(lambda x: (x != 0).all())
    expected = Series([False, True], index=index, name="b")
    tm.assert_series_equal(result, expected)

    result = gp["c"].aggregate(lambda x: x.isnull().all())
    expected = Series([True, False], index=index, name="c")
    tm.assert_series_equal(result, expected)


def test_order_aggregate_multiple_funcs():
    # GH 25692
    df = pd.DataFrame({"A": [1, 1, 2, 2], "B": [1, 2, 3, 4]})

    res = df.groupby("A").agg(["sum", "max", "mean", "ohlc", "min"])
    result = res.columns.levels[1]

    expected = pd.Index(["sum", "max", "mean", "ohlc", "min"])

    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
@pytest.mark.parametrize("how", ["first", "last", "min", "max", "mean", "median"])
def test_uint64_type_handling(dtype, how):
    # GH 26310
    df = pd.DataFrame({"x": 6903052872240755750, "y": [1, 2]})
    expected = df.groupby("y").agg({"x": how})
    df.x = df.x.astype(dtype)
    result = df.groupby("y").agg({"x": how})
    result.x = result.x.astype(np.int64)
    tm.assert_frame_equal(result, expected, check_exact=True)


def test_func_duplicates_raises():
    # GH28426
    msg = "Function names"
    df = pd.DataFrame({"A": [0, 0, 1, 1], "B": [1, 2, 3, 4]})
    with pytest.raises(SpecificationError, match=msg):
        df.groupby("A").agg(["min", "min"])


class TestNamedAggregationSeries:
    def test_series_named_agg(self):
        df = pd.Series([1, 2, 3, 4])
        gr = df.groupby([0, 0, 1, 1])
        result = gr.agg(a="sum", b="min")
        expected = pd.DataFrame(
            {"a": [3, 7], "b": [1, 3]}, columns=["a", "b"], index=[0, 1]
        )
        tm.assert_frame_equal(result, expected)

        result = gr.agg(b="min", a="sum")
        expected = expected[["b", "a"]]
        tm.assert_frame_equal(result, expected)

    def test_no_args_raises(self):
        gr = pd.Series([1, 2]).groupby([0, 1])
        with pytest.raises(TypeError, match="Must provide"):
            gr.agg()

        # but we do allow this
        result = gr.agg([])
        expected = pd.DataFrame()
        tm.assert_frame_equal(result, expected)

    def test_series_named_agg_duplicates_no_raises(self):
        # GH28426
        gr = pd.Series([1, 2, 3]).groupby([0, 0, 1])
        grouped = gr.agg(a="sum", b="sum")
        expected = pd.DataFrame({"a": [3, 3], "b": [3, 3]})
        tm.assert_frame_equal(expected, grouped)

    def test_mangled(self):
        gr = pd.Series([1, 2, 3]).groupby([0, 0, 1])
        result = gr.agg(a=lambda x: 0, b=lambda x: 1)
        expected = pd.DataFrame({"a": [0, 0], "b": [1, 1]})
        tm.assert_frame_equal(result, expected)


class TestNamedAggregationDataFrame:
    def test_agg_relabel(self):
        df = pd.DataFrame(
            {"group": ["a", "a", "b", "b"], "A": [0, 1, 2, 3], "B": [5, 6, 7, 8]}
        )
        result = df.groupby("group").agg(a_max=("A", "max"), b_max=("B", "max"))
        expected = pd.DataFrame(
            {"a_max": [1, 3], "b_max": [6, 8]},
            index=pd.Index(["a", "b"], name="group"),
            columns=["a_max", "b_max"],
        )
        tm.assert_frame_equal(result, expected)

        # order invariance
        p98 = functools.partial(np.percentile, q=98)
        result = df.groupby("group").agg(
            b_min=("B", "min"),
            a_min=("A", min),
            a_mean=("A", np.mean),
            a_max=("A", "max"),
            b_max=("B", "max"),
            a_98=("A", p98),
        )
        expected = pd.DataFrame(
            {
                "b_min": [5, 7],
                "a_min": [0, 2],
                "a_mean": [0.5, 2.5],
                "a_max": [1, 3],
                "b_max": [6, 8],
                "a_98": [0.98, 2.98],
            },
            index=pd.Index(["a", "b"], name="group"),
            columns=["b_min", "a_min", "a_mean", "a_max", "b_max", "a_98"],
        )
        tm.assert_frame_equal(result, expected)

    def test_agg_relabel_non_identifier(self):
        df = pd.DataFrame(
            {"group": ["a", "a", "b", "b"], "A": [0, 1, 2, 3], "B": [5, 6, 7, 8]}
        )

        result = df.groupby("group").agg(**{"my col": ("A", "max")})
        expected = pd.DataFrame(
            {"my col": [1, 3]}, index=pd.Index(["a", "b"], name="group")
        )
        tm.assert_frame_equal(result, expected)

    def test_duplicate_no_raises(self):
        # GH 28426, if use same input function on same column,
        # no error should raise
        df = pd.DataFrame({"A": [0, 0, 1, 1], "B": [1, 2, 3, 4]})

        grouped = df.groupby("A").agg(a=("B", "min"), b=("B", "min"))
        expected = pd.DataFrame(
            {"a": [1, 3], "b": [1, 3]}, index=pd.Index([0, 1], name="A")
        )
        tm.assert_frame_equal(grouped, expected)

        quant50 = functools.partial(np.percentile, q=50)
        quant70 = functools.partial(np.percentile, q=70)
        quant50.__name__ = "quant50"
        quant70.__name__ = "quant70"

        test = pd.DataFrame(
            {"col1": ["a", "a", "b", "b", "b"], "col2": [1, 2, 3, 4, 5]}
        )

        grouped = test.groupby("col1").agg(
            quantile_50=("col2", quant50), quantile_70=("col2", quant70)
        )
        expected = pd.DataFrame(
            {"quantile_50": [1.5, 4.0], "quantile_70": [1.7, 4.4]},
            index=pd.Index(["a", "b"], name="col1"),
        )
        tm.assert_frame_equal(grouped, expected)

    def test_agg_relabel_with_level(self):
        df = pd.DataFrame(
            {"A": [0, 0, 1, 1], "B": [1, 2, 3, 4]},
            index=pd.MultiIndex.from_product([["A", "B"], ["a", "b"]]),
        )
        result = df.groupby(level=0).agg(
            aa=("A", "max"), bb=("A", "min"), cc=("B", "mean")
        )
        expected = pd.DataFrame(
            {"aa": [0, 1], "bb": [0, 1], "cc": [1.5, 3.5]}, index=["A", "B"]
        )
        tm.assert_frame_equal(result, expected)

    def test_agg_relabel_other_raises(self):
        df = pd.DataFrame({"A": [0, 0, 1], "B": [1, 2, 3]})
        grouped = df.groupby("A")
        match = "Must provide"
        with pytest.raises(TypeError, match=match):
            grouped.agg(foo=1)

        with pytest.raises(TypeError, match=match):
            grouped.agg()

        with pytest.raises(TypeError, match=match):
            grouped.agg(a=("B", "max"), b=(1, 2, 3))

    def test_missing_raises(self):
        df = pd.DataFrame({"A": [0, 1], "B": [1, 2]})
        with pytest.raises(KeyError, match="Column 'C' does not exist"):
            df.groupby("A").agg(c=("C", "sum"))

    def test_agg_namedtuple(self):
        df = pd.DataFrame({"A": [0, 1], "B": [1, 2]})
        result = df.groupby("A").agg(
            b=pd.NamedAgg("B", "sum"), c=pd.NamedAgg(column="B", aggfunc="count")
        )
        expected = df.groupby("A").agg(b=("B", "sum"), c=("B", "count"))
        tm.assert_frame_equal(result, expected)

    def test_mangled(self):
        df = pd.DataFrame({"A": [0, 1], "B": [1, 2], "C": [3, 4]})
        result = df.groupby("A").agg(b=("B", lambda x: 0), c=("C", lambda x: 1))
        expected = pd.DataFrame(
            {"b": [0, 0], "c": [1, 1]}, index=pd.Index([0, 1], name="A")
        )
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "agg_col1, agg_col2, agg_col3, agg_result1, agg_result2, agg_result3",
    [
        (
            (("y", "A"), "max"),
            (("y", "A"), np.min),
            (("y", "B"), "mean"),
            [1, 3],
            [0, 2],
            [5.5, 7.5],
        ),
        (
            (("y", "A"), lambda x: max(x)),
            (("y", "A"), lambda x: 1),
            (("y", "B"), "mean"),
            [1, 3],
            [1, 1],
            [5.5, 7.5],
        ),
        (
            pd.NamedAgg(("y", "A"), "max"),
            pd.NamedAgg(("y", "B"), np.mean),
            pd.NamedAgg(("y", "A"), lambda x: 1),
            [1, 3],
            [5.5, 7.5],
            [1, 1],
        ),
    ],
)
def test_agg_relabel_multiindex_column(
    agg_col1, agg_col2, agg_col3, agg_result1, agg_result2, agg_result3
):
    # GH 29422, add tests for multiindex column cases
    df = DataFrame(
        {"group": ["a", "a", "b", "b"], "A": [0, 1, 2, 3], "B": [5, 6, 7, 8]}
    )
    df.columns = pd.MultiIndex.from_tuples([("x", "group"), ("y", "A"), ("y", "B")])
    idx = pd.Index(["a", "b"], name=("x", "group"))

    result = df.groupby(("x", "group")).agg(a_max=(("y", "A"), "max"))
    expected = DataFrame({"a_max": [1, 3]}, index=idx)
    tm.assert_frame_equal(result, expected)

    result = df.groupby(("x", "group")).agg(
        col_1=agg_col1, col_2=agg_col2, col_3=agg_col3
    )
    expected = DataFrame(
        {"col_1": agg_result1, "col_2": agg_result2, "col_3": agg_result3}, index=idx
    )
    tm.assert_frame_equal(result, expected)


def test_agg_relabel_multiindex_raises_not_exist():
    # GH 29422, add test for raises senario when aggregate column does not exist
    df = DataFrame(
        {"group": ["a", "a", "b", "b"], "A": [0, 1, 2, 3], "B": [5, 6, 7, 8]}
    )
    df.columns = pd.MultiIndex.from_tuples([("x", "group"), ("y", "A"), ("y", "B")])

    with pytest.raises(KeyError, match="does not exist"):
        df.groupby(("x", "group")).agg(a=(("Y", "a"), "max"))


def test_agg_relabel_multiindex_duplicates():
    # GH29422, add test for raises senario when getting duplicates
    # GH28426, after this change, duplicates should also work if the relabelling is
    # different
    df = DataFrame(
        {"group": ["a", "a", "b", "b"], "A": [0, 1, 2, 3], "B": [5, 6, 7, 8]}
    )
    df.columns = pd.MultiIndex.from_tuples([("x", "group"), ("y", "A"), ("y", "B")])

    result = df.groupby(("x", "group")).agg(
        a=(("y", "A"), "min"), b=(("y", "A"), "min")
    )
    idx = pd.Index(["a", "b"], name=("x", "group"))
    expected = DataFrame({"a": [0, 2], "b": [0, 2]}, index=idx)
    tm.assert_frame_equal(result, expected)


def myfunc(s):
    return np.percentile(s, q=0.90)


@pytest.mark.parametrize("func", [lambda s: np.percentile(s, q=0.90), myfunc])
def test_lambda_named_agg(func):
    # see gh-28467
    animals = DataFrame(
        {
            "kind": ["cat", "dog", "cat", "dog"],
            "height": [9.1, 6.0, 9.5, 34.0],
            "weight": [7.9, 7.5, 9.9, 198.0],
        }
    )

    result = animals.groupby("kind").agg(
        mean_height=("height", "mean"), perc90=("height", func)
    )
    expected = DataFrame(
        [[9.3, 9.1036], [20.0, 6.252]],
        columns=["mean_height", "perc90"],
        index=Index(["cat", "dog"], name="kind"),
    )

    tm.assert_frame_equal(result, expected)


def test_aggregate_mixed_types():
    # GH 16916
    df = pd.DataFrame(
        data=np.array([0] * 9).reshape(3, 3), columns=list("XYZ"), index=list("abc")
    )
    df["grouping"] = ["group 1", "group 1", 2]
    result = df.groupby("grouping").aggregate(lambda x: x.tolist())
    expected_data = [[[0], [0], [0]], [[0, 0], [0, 0], [0, 0]]]
    expected = pd.DataFrame(
        expected_data,
        index=Index([2, "group 1"], dtype="object", name="grouping"),
        columns=Index(["X", "Y", "Z"], dtype="object"),
    )
    tm.assert_frame_equal(result, expected)


class TestLambdaMangling:
    def test_basic(self):
        df = pd.DataFrame({"A": [0, 0, 1, 1], "B": [1, 2, 3, 4]})
        result = df.groupby("A").agg({"B": [lambda x: 0, lambda x: 1]})

        expected = pd.DataFrame(
            {("B", "<lambda_0>"): [0, 0], ("B", "<lambda_1>"): [1, 1]},
            index=pd.Index([0, 1], name="A"),
        )
        tm.assert_frame_equal(result, expected)

    def test_mangle_series_groupby(self):
        gr = pd.Series([1, 2, 3, 4]).groupby([0, 0, 1, 1])
        result = gr.agg([lambda x: 0, lambda x: 1])
        expected = pd.DataFrame({"<lambda_0>": [0, 0], "<lambda_1>": [1, 1]})
        tm.assert_frame_equal(result, expected)

    @pytest.mark.xfail(reason="GH-26611. kwargs for multi-agg.")
    def test_with_kwargs(self):
        f1 = lambda x, y, b=1: x.sum() + y + b
        f2 = lambda x, y, b=2: x.sum() + y * b
        result = pd.Series([1, 2]).groupby([0, 0]).agg([f1, f2], 0)
        expected = pd.DataFrame({"<lambda_0>": [4], "<lambda_1>": [6]})
        tm.assert_frame_equal(result, expected)

        result = pd.Series([1, 2]).groupby([0, 0]).agg([f1, f2], 0, b=10)
        expected = pd.DataFrame({"<lambda_0>": [13], "<lambda_1>": [30]})
        tm.assert_frame_equal(result, expected)

    def test_agg_with_one_lambda(self):
        # GH 25719, write tests for DataFrameGroupby.agg with only one lambda
        df = pd.DataFrame(
            {
                "kind": ["cat", "dog", "cat", "dog"],
                "height": [9.1, 6.0, 9.5, 34.0],
                "weight": [7.9, 7.5, 9.9, 198.0],
            }
        )

        columns = ["height_sqr_min", "height_max", "weight_max"]
        expected = pd.DataFrame(
            {
                "height_sqr_min": [82.81, 36.00],
                "height_max": [9.5, 34.0],
                "weight_max": [9.9, 198.0],
            },
            index=pd.Index(["cat", "dog"], name="kind"),
            columns=columns,
        )

        # check pd.NameAgg case
        result1 = df.groupby(by="kind").agg(
            height_sqr_min=pd.NamedAgg(
                column="height", aggfunc=lambda x: np.min(x ** 2)
            ),
            height_max=pd.NamedAgg(column="height", aggfunc="max"),
            weight_max=pd.NamedAgg(column="weight", aggfunc="max"),
        )
        tm.assert_frame_equal(result1, expected)

        # check agg(key=(col, aggfunc)) case
        result2 = df.groupby(by="kind").agg(
            height_sqr_min=("height", lambda x: np.min(x ** 2)),
            height_max=("height", "max"),
            weight_max=("weight", "max"),
        )
        tm.assert_frame_equal(result2, expected)

    def test_agg_multiple_lambda(self):
        # GH25719, test for DataFrameGroupby.agg with multiple lambdas
        # with mixed aggfunc
        df = pd.DataFrame(
            {
                "kind": ["cat", "dog", "cat", "dog"],
                "height": [9.1, 6.0, 9.5, 34.0],
                "weight": [7.9, 7.5, 9.9, 198.0],
            }
        )
        columns = [
            "height_sqr_min",
            "height_max",
            "weight_max",
            "height_max_2",
            "weight_min",
        ]
        expected = pd.DataFrame(
            {
                "height_sqr_min": [82.81, 36.00],
                "height_max": [9.5, 34.0],
                "weight_max": [9.9, 198.0],
                "height_max_2": [9.5, 34.0],
                "weight_min": [7.9, 7.5],
            },
            index=pd.Index(["cat", "dog"], name="kind"),
            columns=columns,
        )

        # check agg(key=(col, aggfunc)) case
        result1 = df.groupby(by="kind").agg(
            height_sqr_min=("height", lambda x: np.min(x ** 2)),
            height_max=("height", "max"),
            weight_max=("weight", "max"),
            height_max_2=("height", lambda x: np.max(x)),
            weight_min=("weight", lambda x: np.min(x)),
        )
        tm.assert_frame_equal(result1, expected)

        # check pd.NamedAgg case
        result2 = df.groupby(by="kind").agg(
            height_sqr_min=pd.NamedAgg(
                column="height", aggfunc=lambda x: np.min(x ** 2)
            ),
            height_max=pd.NamedAgg(column="height", aggfunc="max"),
            weight_max=pd.NamedAgg(column="weight", aggfunc="max"),
            height_max_2=pd.NamedAgg(column="height", aggfunc=lambda x: np.max(x)),
            weight_min=pd.NamedAgg(column="weight", aggfunc=lambda x: np.min(x)),
        )
        tm.assert_frame_equal(result2, expected)
