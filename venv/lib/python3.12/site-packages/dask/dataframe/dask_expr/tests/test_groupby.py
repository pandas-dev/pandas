from __future__ import annotations

import re
from collections import OrderedDict
from functools import partial

import numpy as np
import pytest

import dask
from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr._groupby import Aggregation, GroupByUDFBlockwise
from dask.dataframe.dask_expr._reductions import TreeReduce
from dask.dataframe.dask_expr._shuffle import Shuffle, TaskShuffle, divisions_lru
from dask.dataframe.dask_expr.io import FromPandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq, xfail_gpu

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": list(range(10)) * 10, "y": range(100), "z": 1})
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=4)


def test_groupby_unsupported_by(pdf, df):
    assert_eq(df.groupby(df.x).sum(), pdf.groupby(pdf.x).sum())


@pytest.mark.parametrize("split_every", [None, 5])
@pytest.mark.parametrize(
    "api",
    ["sum", "mean", "min", "max", "prod", "first", "last", "var", "std", "idxmin"],
)
@pytest.mark.parametrize(
    "numeric_only",
    [
        pytest.param(True, marks=xfail_gpu("numeric_only not supported by cudf")),
        False,
    ],
)
def test_groupby_numeric(pdf, df, api, numeric_only, split_every):
    if not numeric_only and api in {"var", "std"}:
        pytest.xfail("not implemented")
    g = df.groupby("x")
    agg = getattr(g, api)(numeric_only=numeric_only, split_every=split_every)

    expect = getattr(pdf.groupby("x"), api)(numeric_only=numeric_only)
    assert_eq(agg, expect)

    g = df.y.groupby(df.x)
    agg = getattr(g, api)()

    expect = getattr(pdf.y.groupby(pdf.x), api)()
    assert_eq(agg, expect)

    g = df.groupby("x")
    agg = getattr(g, api)(numeric_only=numeric_only)["y"]

    expect = getattr(pdf.groupby("x"), api)(numeric_only=numeric_only)["y"]
    assert_eq(agg, expect)

    g = df.groupby([df.x])
    agg = getattr(g, api)(numeric_only=numeric_only)["y"]

    expect = getattr(pdf.groupby([pdf.x]), api)(numeric_only=numeric_only)["y"]
    assert_eq(agg, expect)

    g = df.groupby([df.x, df.z])
    agg = getattr(g, api)()
    expect = getattr(pdf.groupby([pdf.x, pdf.z]), api)(numeric_only=numeric_only)
    assert_eq(agg, expect)

    g = df.groupby([df.x, "z"])
    agg = getattr(g, api)()
    expect = getattr(pdf.groupby([pdf.x, "z"]), api)(numeric_only=numeric_only)
    assert_eq(agg, expect)

    pdf = pdf.set_index("x")
    df = from_pandas(pdf, npartitions=10, sort=False)
    g = df.groupby("x")
    agg = getattr(g, api)()
    expect = getattr(pdf.groupby("x"), api)(numeric_only=numeric_only)
    assert_eq(agg, expect)

    g = df.groupby(["x", "z"])
    agg = getattr(g, api)()
    expect = getattr(pdf.groupby(["x", "z"]), api)(numeric_only=numeric_only)
    assert_eq(agg, expect)


def test_groupby_reduction_optimize(pdf, df):
    df = df.replace(1, 5)
    agg = df.groupby(df.x).y.sum()
    expected_query = df[["x", "y"]]
    expected_query = expected_query.groupby(expected_query.x).y.sum()
    assert agg.optimize()._name == expected_query.optimize()._name
    expect = pdf.replace(1, 5).groupby(["x"]).y.sum()
    assert_eq(agg, expect)

    df2 = df[["y"]]
    agg = df2.groupby(df.x).y.sum()
    ops = [
        op for op in agg.expr.optimize(fuse=False).walk() if isinstance(op, FromPandas)
    ]
    assert len(ops) == 1
    assert ops[0].columns == ["x", "y"]

    df2 = df[["y"]]
    with pytest.warns(UserWarning, match="inferred from partial data"):
        agg = df2.groupby(df.x).y.apply(lambda x: x)
    ops = [
        op for op in agg.expr.optimize(fuse=False).walk() if isinstance(op, FromPandas)
    ]
    assert len(ops) == 1
    assert ops[0].columns == ["x", "y"]
    assert_eq(agg, pdf.replace(1, 5).groupby(pdf.replace(1, 5).x).y.apply(lambda x: x))


def test_std_columns_int():
    df = pd.DataFrame({0: [5], 1: [5]})
    ddf = from_pandas(df, npartitions=2)
    assert_eq(ddf.groupby(ddf[0]).std(), df.groupby(df[0]).std())


@pytest.mark.parametrize(
    "func",
    [
        "count",
        pytest.param(
            "value_counts", marks=xfail_gpu("value_counts not supported by cudf")
        ),
        "size",
        "head",
        "tail",
    ],
)
def test_groupby_no_numeric_only(pdf, func):
    pdf = pdf.drop(columns="z")
    df = from_pandas(pdf, npartitions=10)
    g = df.groupby("x")
    if func != "value_counts":
        agg = getattr(g, func)()

        expect = getattr(pdf.groupby("x"), func)()
        assert_eq(agg, expect)

    g = df.y.groupby(df.x)
    agg = getattr(g, func)()

    expect = getattr(pdf.y.groupby(pdf.x), func)()
    assert_eq(agg, expect)


def test_value_counts_split_out(pdf):
    df = from_pandas(pdf, npartitions=10)
    result = df.groupby("x").y.value_counts(split_out=True)
    expected = pdf.groupby("x").y.value_counts()
    assert_eq(result, expected)


def test_unique(df, pdf):
    result = df.groupby("x")["y"].unique()
    expected = pdf.groupby("x")["y"].unique()

    # Use explode because each DataFrame row is a list; equality fails
    assert_eq(result.explode(), expected.explode())

    result = df.y.groupby(df.x).unique()
    expected = pdf.y.groupby(pdf.x).unique()

    # Use explode because each DataFrame row is a list; equality fails
    assert_eq(result.explode(), expected.explode())


def test_groupby_mean_slice(pdf, df):
    g = df.groupby("x")
    agg = g.y.mean()

    expect = pdf.groupby("x").y.mean()
    assert_eq(agg, expect)


def test_groupby_agg_grouper_single(pdf):
    pdf = pdf[["x"]]
    df = from_pandas(pdf, npartitions=2)

    result = df.groupby("x")["x"].agg(["min", "max"])
    expected = pdf.groupby("x")["x"].agg(["min", "max"])
    assert_eq(result, expected)

    result = df.groupby("x")[["x"]].agg(["min", "max"])
    expected = pdf.groupby("x")[["x"]].agg(["min", "max"])
    assert_eq(result, expected)


def test_groupby_slice_agg_reduces(df, pdf):
    result = df.groupby("x")["y"].agg(["min", "max"])
    expected = pdf.groupby("x")["y"].agg(["min", "max"])
    assert_eq(result, expected)


def test_groupby_nunique(df, pdf):
    with pytest.raises(AttributeError):
        df.groupby("x").nunique()

    assert_eq(df.groupby("x").y.nunique(split_out=1), pdf.groupby("x").y.nunique())
    assert_eq(df.groupby("x").y.nunique(split_out=True), pdf.groupby("x").y.nunique())
    assert df.groupby("x").y.nunique().npartitions == df.npartitions
    assert_eq(df.y.groupby(df.x).nunique(split_out=1), pdf.y.groupby(pdf.x).nunique())

    pdf = pdf.add_prefix("x")
    df = from_pandas(pdf, npartitions=10)
    assert_eq(df.groupby("xx").xy.nunique(), pdf.groupby("xx").xy.nunique())
    assert_eq(df.xx.groupby(df.xy).nunique(), pdf.xx.groupby(pdf.xy).nunique())


def test_groupby_series(pdf, df):
    pdf_result = pdf.groupby(pdf.x).sum()
    result = df.groupby(df.x).sum()
    assert_eq(result, pdf_result)
    result = df.groupby("x").sum()
    assert_eq(result, pdf_result)

    df2 = from_pandas(pd.DataFrame({"a": [1, 2, 3]}), npartitions=2)

    with pytest.raises(NotImplementedError, match="DataFrames columns"):
        df.groupby(df2.a)


@pytest.mark.parametrize("group_keys", [True, False, None])
def test_groupby_group_keys(group_keys, pdf):
    pdf = pdf.set_index("x")
    df = from_pandas(pdf, npartitions=10)

    func = lambda g: g.copy()
    expected = pdf.groupby("x").apply(func)
    assert_eq(expected, df.groupby("x").apply(func, meta=expected))

    expected = pdf.groupby("x", group_keys=group_keys).apply(func)
    assert_eq(
        expected, df.groupby("x", group_keys=group_keys).apply(func, meta=expected)
    )


def test_dataframe_aggregations_multilevel(df, pdf):
    grouper = lambda df: [df["x"] > 2, df["y"] > 1]

    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        assert_eq(
            pdf.groupby(grouper(pdf)).sum(),
            df.groupby(grouper(df)).sum(split_out=2),
        )


@pytest.mark.parametrize(
    "spec",
    [
        {"x": "count"},
        {"x": ["count"]},
        {"x": ["count"], "y": "mean"},
        {"x": ["sum", "mean"]},
        ["min", "mean"],
        "sum",
        "median",
    ],
)
def test_groupby_agg(pdf, df, spec):
    g = df.groupby("x")
    agg = g.agg(spec)

    expect = pdf.groupby("x").agg(spec)
    assert_eq(agg, expect)

    g = df.groupby(["x", "y"])
    agg = g.agg(spec)

    expect = pdf.groupby(["x", "y"]).agg(spec)
    assert_eq(agg, expect)


@pytest.mark.parametrize("numeric_only", [True, False])
@pytest.mark.parametrize("api", ["cov", "corr"])
def test_groupby_cov(api, df, pdf, numeric_only):
    g = df.groupby("x")
    agg = getattr(g, api)(numeric_only=numeric_only)

    expect = getattr(pdf.groupby("x"), api)(numeric_only=numeric_only)
    assert_eq(agg, expect)


@pytest.mark.parametrize(
    "spec",
    [
        "sum",
        ["sum"],
        ["sum", "mean"],
    ],
)
def test_series_groupby_agg(pdf, df, spec):
    g = df.y.groupby(df.x)
    agg = g.agg(spec)

    expect = pdf.y.groupby(pdf.x).agg(spec)
    assert_eq(agg, expect)


def test_groupby_getitem_agg(pdf, df):
    assert_eq(df.groupby("x").y.sum(), pdf.groupby("x").y.sum())
    assert_eq(df.groupby("x")[["y"]].sum(), pdf.groupby("x")[["y"]].sum())


def test_groupby_agg_column_projection(pdf, df):
    g = df.groupby("x")
    agg = g.agg({"x": "count"}).simplify()

    assert list(agg.frame.columns) == ["x"]
    expect = pdf.groupby("x").agg({"x": "count"})
    assert_eq(agg, expect)


def test_groupby_split_every(pdf):
    df = from_pandas(pdf, npartitions=16)
    query = df.groupby("x").sum(split_out=1)
    tree_reduce_node = list(query.optimize(fuse=False).find_operations(TreeReduce))
    assert len(tree_reduce_node) == 1
    assert tree_reduce_node[0].split_every == 8

    query = df.groupby("x").aggregate({"y": "sum"}, split_out=1)
    tree_reduce_node = list(query.optimize(fuse=False).find_operations(TreeReduce))
    assert len(tree_reduce_node) == 1
    assert tree_reduce_node[0].split_every == 8


def test_groupby_index(pdf):
    pdf = pdf.set_index("x")
    df = from_pandas(pdf, npartitions=10)
    result = df.groupby(df.index).sum()
    expected = pdf.groupby(pdf.index).sum()
    assert_eq(result, expected)
    assert_eq(result["y"], expected["y"])

    result = df.groupby(df.index).var()
    expected = pdf.groupby(pdf.index).var()
    assert_eq(result, expected)
    assert_eq(result["y"], expected["y"])

    result = df.groupby(df.index).agg({"y": "sum"})
    expected = pdf.groupby(pdf.index).agg({"y": "sum"})
    assert_eq(result, expected)


def test_split_out_automatically():
    pdf = pd.DataFrame({"a": [1, 2, 3] * 1_000, "b": 1, "c": 1, "d": 1})
    df = from_pandas(pdf, npartitions=500)
    q = df.groupby("a").sum()
    assert q.optimize().npartitions == 34
    expected = pdf.groupby("a").sum()
    assert_eq(q, expected)

    q = df.groupby(["a", "b"]).sum()
    assert q.optimize().npartitions == 50
    expected = pdf.groupby(["a", "b"]).sum()
    assert_eq(q, expected)

    q = df.groupby(["a", "b", "c"]).sum()
    assert q.optimize().npartitions == 100
    expected = pdf.groupby(["a", "b", "c"]).sum()
    assert_eq(q, expected)


def test_split_out_sort_values_compute(pdf, df):
    divisions_lru.data = OrderedDict()
    result = df.groupby("x").sum(split_out=2).sort_values(by="y").compute()
    assert len(divisions_lru.data) == 0
    expected = pdf.groupby("x").sum().sort_values(by="y")
    assert_eq(result, expected)


def test_groupby_repartition_to_one(pdf, df):
    df = from_pandas(pdf, npartitions=25)
    result = (
        df.groupby("x", sort=True).sum(split_out=2).repartition(npartitions=1).compute()
    )
    expected = pdf.groupby("x").sum()
    assert_eq(result, expected)


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated ")
def test_groupby_apply(df, pdf):
    def test(x):
        x["new"] = x.sum().sum()
        return x

    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(df.groupby(df.x).apply(test), pdf.groupby(pdf.x).apply(test))
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(
            df.groupby(df.x, group_keys=False).apply(test),
            pdf.groupby(pdf.x, group_keys=False).apply(test),
        )
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(df.groupby("x").apply(test), pdf.groupby("x").apply(test))
    assert_eq(
        df.groupby("x").apply(test, meta=pdf.groupby("x").apply(test).head(0)),
        pdf.groupby("x").apply(test),
    )
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(
            df.groupby(["x", "y"]).apply(test), pdf.groupby(["x", "y"]).apply(test)
        )

    with pytest.warns(UserWarning, match="inferred from partial data"):
        query = df.groupby("x").apply(test).optimize(fuse=False)
    assert query.expr.find_operations(Shuffle)
    assert query.expr.find_operations(GroupByUDFBlockwise)

    with pytest.warns(UserWarning, match="inferred from partial data"):
        query = df.groupby("x")[["y"]].apply(test).simplify()
    with pytest.warns(UserWarning, match="inferred from partial data"):
        expected = df[["x", "y"]].groupby("x")[["y"]].apply(test).simplify()
    assert query._name == expected._name
    assert_eq(query, pdf.groupby("x")[["y"]].apply(test))


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated ")
def test_groupby_transform(df, pdf):
    def test(x):
        return x

    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(df.groupby(df.x).transform(test), pdf.groupby(pdf.x).transform(test))
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(df.groupby("x").transform(test), pdf.groupby("x").transform(test))
    assert_eq(
        df.groupby("x").transform(test, meta=pdf.groupby("x").transform(test).head(0)),
        pdf.groupby("x").transform(test),
    )

    with pytest.warns(UserWarning, match="inferred from partial data"):
        query = df.groupby("x").transform(test).optimize(fuse=False)
    assert query.expr.find_operations(Shuffle)
    assert query.expr.find_operations(GroupByUDFBlockwise)

    with pytest.warns(UserWarning, match="inferred from partial data"):
        query = df.groupby("x")[["y"]].transform(test).simplify()
    with pytest.warns(UserWarning, match="inferred from partial data"):
        expected = df[["x", "y"]].groupby("x")[["y"]].transform(test).simplify()
    assert query._name == expected._name
    assert_eq(query, pdf.groupby("x")[["y"]].transform(test))


def test_groupby_shift(df, pdf):
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(
            df.groupby(df.x).shift(periods=1), pdf.groupby(pdf.x).shift(periods=1)
        )
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(df.groupby("x").shift(periods=1), pdf.groupby("x").shift(periods=1))
    assert_eq(
        df.groupby("x").shift(
            periods=1, meta=pdf.groupby("x").shift(periods=1).head(0)
        ),
        pdf.groupby("x").shift(periods=1),
    )

    with pytest.warns(UserWarning, match="inferred from partial data"):
        query = df.groupby("x").shift(periods=1).optimize(fuse=False)
    assert query.expr.find_operations(Shuffle)
    assert query.expr.find_operations(GroupByUDFBlockwise)

    with pytest.warns(UserWarning, match="inferred from partial data"):
        query = df.groupby("x")[["y"]].shift(periods=1).simplify()
    with pytest.warns(UserWarning, match="inferred from partial data"):
        expected = df[["x", "y"]].groupby("x")[["y"]].shift(periods=1).simplify()
    assert query._name == expected._name
    assert_eq(query, pdf.groupby("x")[["y"]].shift(periods=1))


def test_size(pdf, df):
    assert_eq(df.groupby("x").agg("size"), pdf.groupby("x").agg("size"))


def test_groupby_numeric_only_lambda_caller(df, pdf):
    assert_eq(
        df.groupby(lambda x: x // 2).mean(numeric_only=False),
        pdf.groupby(lambda x: x // 2).mean(numeric_only=False),
    )


@pytest.mark.parametrize(
    "api", ["sum", "mean", "min", "max", "prod", "var", "std", "size"]
)
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("split_out", [1, 2])
def test_groupby_single_agg_split_out(pdf, df, api, sort, split_out):
    g = df.groupby("x", sort=sort)
    agg = getattr(g, api)(split_out=split_out)

    expect = getattr(pdf.groupby("x", sort=sort), api)()
    assert_eq(agg, expect, sort_results=not sort)

    g = df.y.groupby(df.x, sort=sort)
    agg = getattr(g, api)(split_out=split_out)
    expect = getattr(pdf.y.groupby(pdf.x, sort=sort), api)()
    assert_eq(agg, expect, sort_results=not sort)

    g = df.y.groupby([df.x, df.z], sort=sort)
    agg = getattr(g, api)(split_out=split_out)
    expect = getattr(pdf.y.groupby([pdf.x, pdf.z], sort=sort), api)()
    assert_eq(agg, expect, sort_results=not sort)


@pytest.mark.parametrize("cow", [True, False])
@pytest.mark.parametrize(
    "func",
    [
        lambda grouped: grouped.apply(lambda x: x.sum()),
        lambda grouped: grouped.transform(lambda x: x.sum()),
    ],
)
def test_apply_or_transform_shuffle_multilevel(pdf, df, func, cow):
    with pd.option_context("mode.copy_on_write", cow):
        grouper = lambda df: [df["x"] + 1, df["y"] + 1]

        with pytest.warns(UserWarning):
            # DataFrameGroupBy
            assert_eq(func(df.groupby(grouper(df))), func(pdf.groupby(grouper(pdf))))

            # SeriesGroupBy
            assert_eq(
                func(df.groupby(grouper(df))["z"]), func(pdf.groupby(grouper(pdf))["z"])
            )

            # DataFrameGroupBy with column slice
            assert_eq(
                func(df.groupby(grouper(df))[["z"]]),
                func(pdf.groupby(grouper(pdf))[["z"]]),
            )


@pytest.mark.parametrize(
    "spec",
    [
        {"x": "count"},
        {"x": ["count"]},
        {"x": ["count"], "y": "mean"},
        {"x": ["sum", "mean"]},
        ["min", "mean"],
        "sum",
    ],
)
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("split_out", [1, 2])
def test_groupby_agg_split_out(pdf, df, spec, sort, split_out):
    g = df.groupby("x", sort=sort)
    agg = g.agg(spec, split_out=split_out)

    expect = pdf.groupby("x", sort=sort).agg(spec)
    assert_eq(agg, expect, sort_results=not sort)


def test_groupby_reduction_shuffle(df, pdf):
    q = df.groupby("x").sum(split_out=True)
    assert q.optimize().npartitions == df.npartitions
    expected = pdf.groupby("x").sum()
    assert_eq(q, expected)


def test_groupby_projection_split_out(df, pdf):
    pdf_result = pdf.groupby("x")["y"].sum()
    result = df.groupby("x")["y"].sum(split_out=2)
    assert_eq(result, pdf_result)

    pdf_result = pdf.groupby("y")["x"].sum()
    df = from_pandas(pdf, npartitions=50)
    result = df.groupby("y")["x"].sum(split_out=2)
    assert_eq(result, pdf_result)


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated ")
def test_numeric_column_names():
    df = pd.DataFrame({0: [0, 1, 0, 1], 1: [1, 2, 3, 4], 2: [0, 1, 0, 1]})
    ddf = from_pandas(df, npartitions=2)
    assert_eq(ddf.groupby(0).sum(), df.groupby(0).sum())
    assert_eq(ddf.groupby([0, 2]).sum(), df.groupby([0, 2]).sum())
    expected = df.groupby(0).apply(lambda x: x)
    assert_eq(
        ddf.groupby(0).apply(lambda x: x, meta=expected),
        expected,
    )


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated ")
def test_apply_divisions(pdf):
    pdf = pdf.set_index("x")
    df = from_pandas(pdf, npartitions=10)
    with pytest.warns(UserWarning, match="inferred from partial data"):
        result = df.groupby(["x", "y"]).apply(lambda x: x)
    assert df.divisions == result.divisions
    assert_eq(result, pdf.groupby(["x", "y"]).apply(lambda x: x))


def test_groupby_co_aligned_grouper(df, pdf):
    assert_eq(
        df[["y"]].groupby(df["x"]).sum(),
        pdf[["y"]].groupby(pdf["x"]).sum(),
    )


@pytest.mark.parametrize("func", ["var", "std"])
@pytest.mark.parametrize("observed", [True, False])
@pytest.mark.parametrize("dropna", [True, False])
def test_groupby_var_dropna_observed(dropna, observed, func):
    df = pd.DataFrame(
        {
            "a": [11, 12, 31, 1, 2, 3, 4, 5, 6, 10],
            "b": pd.Categorical(values=[1] * 9 + [np.nan], categories=[1, 2]),
        }
    )
    ddf = from_pandas(df, npartitions=3)
    dd_result = getattr(ddf.groupby("b", observed=observed, dropna=dropna), func)()
    pdf_result = getattr(df.groupby("b", observed=observed, dropna=dropna), func)()
    assert_eq(dd_result, pdf_result)


def test_groupby_median(df, pdf):
    assert_eq(df.groupby("x").median(), pdf.groupby("x").median())
    q = df.groupby("x").median(split_out=2)
    assert q.optimize().npartitions == 2
    assert_eq(q, pdf.groupby("x").median())
    assert_eq(df.groupby("x")["y"].median(), pdf.groupby("x")["y"].median())
    assert_eq(df.groupby("x").median()["y"], pdf.groupby("x").median()["y"])
    assert_eq(
        df.groupby("x").median(split_every=2)["y"], pdf.groupby("x").median()["y"]
    )
    assert df.groupby("x").median(split_every=2).npartitions == df.npartitions // 2
    assert (
        df.groupby("x").median(split_every=2).optimize().npartitions
        == df.npartitions // 2
    )


def test_groupby_apply_args(df, pdf):
    with pytest.warns(UserWarning, match="inferred from partial data"):
        assert_eq(
            df.groupby("x").apply(lambda x, y: x + y, 1),
            pdf.groupby("x").apply(lambda x, y: x + y, 1),
        )


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated ")
def test_groupby_ffill_bfill(pdf):
    pdf["y"] = pdf["y"].astype("float64")

    pdf.iloc[np.arange(0, len(pdf) - 1, 3), 1] = np.nan
    df = from_pandas(pdf, npartitions=10)
    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        assert_eq(df.groupby("x").ffill(), pdf.groupby("x").ffill())
        assert_eq(df.groupby("x").bfill(), pdf.groupby("x").bfill())

        actual = df.groupby("x")["y"].ffill()
        expect = df[["x", "y"]].groupby("x")["y"].ffill()
        assert actual.optimize()._name == expect.optimize()._name
        assert_eq(actual, pdf.groupby("x")["y"].ffill())

        actual = df.groupby("x").ffill()["y"]
        expect = df[["x", "y"]].groupby("x").ffill()["y"]
        assert actual.optimize()._name == expect.optimize()._name
        assert_eq(actual, pdf.groupby("x")["y"].ffill())


@pytest.mark.parametrize(
    "by",
    [
        lambda df: "x",
        lambda df: df.x,
        lambda df: df.x + 1,
    ],
)
def test_get_group(df, pdf, by):
    ddgrouped = df.groupby(by(df))
    pdgrouped = pdf.groupby(by(pdf))

    # DataFrame
    assert_eq(ddgrouped.get_group(2), pdgrouped.get_group(2))
    assert_eq(ddgrouped.get_group(3), pdgrouped.get_group(3))
    # Series
    assert_eq(ddgrouped.y.get_group(3), pdgrouped.y.get_group(3))
    assert_eq(ddgrouped.y.get_group(2), pdgrouped.y.get_group(2))


def test_groupby_rolling():
    df = pd.DataFrame(
        {
            "column1": range(600),
            "group1": 5 * ["g" + str(i) for i in range(120)],
        },
        index=pd.date_range("20190101", periods=60).repeat(10),
    )

    ddf = from_pandas(df, npartitions=8)

    expected = df.groupby("group1").rolling("1D").sum()
    actual = ddf.groupby("group1").rolling("1D").sum()

    assert_eq(expected, actual, check_divisions=False)

    expected = df.groupby("group1").column1.rolling("1D").mean()
    actual = ddf.groupby("group1").column1.rolling("1D").mean()

    assert_eq(expected, actual, check_divisions=False)

    # Integer window w/ DateTimeIndex
    expected = df.groupby("group1").rolling(1).sum()
    actual = ddf.groupby("group1").rolling(1).sum()
    assert_eq(expected, actual, check_divisions=False)

    # Integer window w/o DateTimeIndex
    expected = df.reset_index(drop=True).groupby("group1").rolling(1).sum()
    actual = (
        from_pandas(df.reset_index(drop=True), npartitions=10)
        .groupby("group1")
        .rolling(1)
        .sum()
    )
    assert_eq(expected, actual, check_divisions=False)

    # Integer window fails w/ datetime in groupby
    with pytest.raises(pd.errors.DataError, match="Cannot aggregate non-numeric type"):
        df.reset_index().groupby("group1").rolling(1).sum()
    with pytest.raises(pd.errors.DataError, match="Cannot aggregate non-numeric type"):
        from_pandas(df.reset_index(), npartitions=10).groupby("group1").rolling(1).sum()


def test_rolling_groupby_projection():
    df = pd.DataFrame(
        {
            "column1": range(600),
            "a": 1,
            "group1": 5 * ["g" + str(i) for i in range(120)],
        },
        index=pd.date_range("20190101", periods=60).repeat(10),
    )

    ddf = from_pandas(df, npartitions=8)

    actual = ddf.groupby("group1").rolling("1D").sum()["column1"]
    expected = df.groupby("group1").rolling("1D").sum()["column1"]

    assert_eq(expected, actual, check_divisions=False)

    optimal = (
        ddf[["column1", "group1"]].groupby("group1").rolling("1D").sum()["column1"]
    )

    assert actual.optimize()._name == (optimal.optimize()._name)


def test_std_var_slice(pdf, df):
    assert_eq(df.groupby("x").y.std(), pdf.groupby("x").y.std())
    assert_eq(df.groupby("x").y.var(), pdf.groupby("x").y.var())


def test_groupby_error(df):
    with pytest.raises(KeyError):
        df.groupby("A")

    dp = df.groupby("y")

    msg = "Column not found: "
    with pytest.raises(KeyError, match=msg):
        dp["A"]

    msg = re.escape(
        "DataFrameGroupBy does not allow compute method."
        "Please chain it with an aggregation method (like ``.mean()``) or get a "
        "specific group using ``.get_group()`` before calling ``compute()``"
    )

    with pytest.raises(NotImplementedError, match=msg):
        dp.compute()


def test_groupby_dir(df):
    assert "y" in dir(df.groupby("x"))


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated ")
def test_groupby_udf_user_warning(df, pdf):
    def func(df):
        return df + 1

    expected = pdf.groupby("x").apply(func)
    with pytest.warns(UserWarning, match="`meta` is not specified"):
        assert_eq(expected, df.groupby("x").apply(func))

    expected = pdf.groupby("x").transform(func)
    with pytest.warns(UserWarning, match="`meta` is not specified"):
        assert_eq(expected, df.groupby("x").transform(func))


def test_groupby_index_array(pdf):
    pdf.index = pd.date_range(start="2020-12-31", freq="D", periods=len(pdf))
    df = from_pandas(pdf, npartitions=10)
    # pandas loses the freq for mac but keeps it on ubuntu
    expected = pdf.x.groupby(pdf.index).nunique()
    expected.index.freq = None

    assert_eq(
        df.x.groupby(df.index).nunique(),
        expected,
        check_names=False,
    )
    # pandas loses the freq for mac but keeps it on ubuntu
    expected = pdf.groupby(pdf.index).x.nunique()
    expected.index.freq = None
    assert_eq(
        df.groupby(df.index).x.nunique(),
        expected,
        check_names=False,
    )


def test_groupby_median_numeric_only():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1, "c": "a"})
    df = from_pandas(pdf, npartitions=2)

    assert_eq(
        df.groupby("b").median(numeric_only=True),
        pdf.groupby("b").median(numeric_only=True),
    )


def test_groupby_median_series():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    df = from_pandas(pdf, npartitions=2)
    assert_eq(
        df.a.groupby(df.b).median(),
        pdf.a.groupby(pdf.b).median(),
    )


@pytest.mark.parametrize("min_count", [0, 1, 2, 3])
@pytest.mark.parametrize("op", ("prod", "sum"))
def test_with_min_count(min_count, op):
    dfs = [
        pd.DataFrame(
            {
                "group": ["A", "A", "B"],
                "val1": [np.nan, 2, 3],
                "val2": [np.nan, 5, 6],
                "val3": [5, 4, 9],
            }
        ),
        pd.DataFrame(
            {
                "group": ["A", "A", "B"],
                "val1": [2, np.nan, np.nan],
                "val2": [np.nan, 5, 6],
                "val3": [5, 4, 9],
            }
        ),
    ]
    ddfs = [from_pandas(df, npartitions=4) for df in dfs]

    for df, ddf in zip(dfs, ddfs):
        assert_eq(
            getattr(df.groupby("group"), op)(min_count=min_count),
            getattr(ddf.groupby("group"), op)(min_count=min_count),
        )


@pytest.mark.parametrize("sel", ["a", "c", "d", ["a", "b"], ["c", "d"]])
@pytest.mark.parametrize("key", ["a", ["a", "b"]])
@pytest.mark.parametrize("func", ["cumsum", "cumprod", "cumcount"])
def test_cumulative(func, key, sel):
    df = pd.DataFrame(
        {
            "a": [1, 2, 6, 4, 4, 6, 4, 3, 7] * 6,
            "b": [4, 2, 7, 3, 3, 1, 1, 1, 2] * 6,
            "c": np.random.randn(54),
            "d": np.random.randn(54),
        },
        columns=["a", "b", "c", "d"],
    )
    df.iloc[[-18, -12, -6], -1] = np.nan
    ddf = from_pandas(df, npartitions=10)

    g, dg = (d.groupby(key)[sel] for d in (df, ddf))
    assert_eq(getattr(g, func)(), getattr(dg, func)())


@pytest.mark.parametrize("by", ["key1", ["key1", "key2"]])
@pytest.mark.parametrize(
    "slice_key",
    [
        3,
        "value",
        ["value"],
        ("value",),
        pd.Index(["value"]),
        pd.Series(["value"]),
    ],
)
def test_groupby_slice_getitem(by, slice_key):
    pdf = pd.DataFrame(
        {
            "key1": ["a", "b", "a"],
            "key2": ["c", "c", "c"],
            "value": [1, 2, 3],
            3: [1, 2, 3],
        }
    )

    ddf = from_pandas(pdf, npartitions=3)
    expect = pdf.groupby(by)[slice_key].count()
    got = ddf.groupby(by)[slice_key].count()
    assert_eq(expect, got)


@pytest.mark.parametrize("npartitions", (1, 2))
@pytest.mark.parametrize("func", ("cumsum", "cumprod", "cumcount"))
def test_series_groupby_cumfunc_with_named_index(npartitions, func):
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5, 6, 7], "y": [8, 9, 6, 2, 3, 5, 6]}
    ).set_index("x")
    ddf = from_pandas(df, npartitions)
    assert ddf.npartitions == npartitions
    expected = getattr(df["y"].groupby("x"), func)()
    result = getattr(ddf["y"].groupby("x"), func)()
    assert_eq(result, expected)


@pytest.mark.parametrize("attr", ("std", "var"))
def test_series_groupby_not_supported(df, attr):
    # See: https://github.com/dask-contrib/dask-expr/issues/840
    g = df.groupby("x")
    with pytest.raises(NotImplementedError, match="Please use `aggregate`"):
        getattr(g.x, attr)()


def test_dataframe_groupby_agg_custom_sum():
    pandas_spec = {"b": "sum"}
    dask_spec = {"b": Aggregation("sum", lambda s: s.sum(), lambda s0: s0.sum())}
    df = pd.DataFrame({"g": [0, 0, 1] * 3, "b": [1, 2, 3] * 3})
    ddf = from_pandas(df, npartitions=2)
    expected = df.groupby("g").aggregate(pandas_spec)
    result = ddf.groupby("g").aggregate(dask_spec)
    assert_eq(result, expected)


def test_groupby_size_drop_columns(df, pdf):
    result = df.groupby("x").size()
    assert result.simplify()._name == df[["x"]].groupby("x").size().simplify()._name
    assert_eq(result, pdf.groupby("x").size())

    result = df.groupby("x")[["y", "z"]].size()
    assert result.simplify()._name == df[["x"]].groupby("x").size().simplify()._name
    assert_eq(result, pdf.groupby("x")[["y", "z"]].size())

    assert_eq(df.groupby("x").y.size(), pdf.groupby("x").y.size())
    assert_eq(df.x.groupby(df.x).size(), pdf.x.groupby(pdf.x).size())

    result = df[["y", "z"]].groupby(df.x).size()
    assert result.simplify()._name == df[[]].groupby(df.x).size().simplify()._name
    assert_eq(result, pdf[["y", "z"]].groupby(pdf.x).size())

    pdf = pdf.set_index("x")
    df = from_pandas(pdf, npartitions=10)
    result = df.groupby("x").size()
    assert_eq(result, pdf.groupby("x").size())
    assert result.simplify()._name == df[[]].groupby("x").size().simplify()._name


@pytest.mark.filterwarnings("ignore:`meta` is not specified|DataFrameGroupBy.apply")
def test_groupby_respect_shuffle_context(df, pdf):
    def _check_task_shuffle(q):
        result = q.optimize(fuse=False)
        assert len([x for x in result.walk() if isinstance(x, TaskShuffle)]) > 0

    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        q = df.groupby("x").apply(lambda x: x)
    _check_task_shuffle(q)

    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        q = df.groupby("x").transform(lambda x: x)
    _check_task_shuffle(q)

    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        q = df.groupby("x").sum(split_out=True)
    _check_task_shuffle(q)

    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        q = df.groupby("x").y.nunique(split_out=True)
    _check_task_shuffle(q)

    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        q = df.groupby("x").median(split_out=True)
    _check_task_shuffle(q)


def test_groupby_agg_meta_error(df, pdf):
    result = df.groupby(["x"]).agg({"y": ["median", "std"]})
    expected = pdf.groupby(["x"]).agg({"y": ["median", "std"]})
    assert_eq(result, expected)


def test_groupby_aggregate_series_split_out(df, pdf):
    result = df.groupby("x").y.agg("sum", split_out=2)
    expected = pdf.groupby("x").y.agg("sum")
    assert_eq(result, expected)


def test_groupby_agg_rename_columns(df, pdf):
    result = df.groupby("x").y.agg(a="sum")
    expected = pdf.groupby("x").y.agg(a="sum")
    assert_eq(result, expected)

    result = df.groupby("x").agg(
        a=pd.NamedAgg("y", aggfunc="sum"),
        b=pd.NamedAgg("z", aggfunc=partial(np.std, ddof=1)),
    )
    expected = pdf.groupby("x").agg(
        a=pd.NamedAgg("y", aggfunc="sum"),
        b=pd.NamedAgg("z", aggfunc=partial(np.std, ddof=1)),
    )
    assert_eq(result, expected)


def test_get_group_multiple_keys():
    pdf = pd.DataFrame({"x": [1, 2, 3], "y": 1, "z": 0})
    df = from_pandas(pdf, npartitions=2)
    result = df.groupby(["x", "y"]).get_group((1, 1))
    expected = pdf.groupby(["x", "y"]).get_group((1, 1))
    assert_eq(result, expected)


def test_groupby_index_modified_divisions():
    date_range = pd.date_range(start="2023-01-01", end="2023-01-02", freq="1min")
    data = {
        "timestamp": date_range,
        "upper_bound_enter": [1] * len(date_range),
        "vwap": [2] * len(date_range),
    }
    pdf = pd.DataFrame(data, index=pd.Index(date_range, name="timestamp"))

    df = from_pandas(pdf, npartitions=8).repartition(freq="1D")
    assert_eq(
        df.groupby(df.index.dt.date).count(),
        pdf.groupby(pdf.index.date).count(),
    )


def test_groupby_getitem_apply_group_keys():
    pdf = pd.DataFrame(
        {
            "A": [0, 1] * 4,
            "B": [1] * 8,
        }
    )
    df = from_pandas(pdf, npartitions=4)
    result = df.groupby("A", group_keys=False).B.apply(lambda x: x, meta=("B", int))
    expected = pdf.groupby("A", group_keys=False).B.apply(lambda x: x)
    assert_eq(result, expected)
