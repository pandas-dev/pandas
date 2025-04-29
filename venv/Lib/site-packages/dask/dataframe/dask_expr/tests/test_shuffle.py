from __future__ import annotations

from collections import OrderedDict

import numpy as np
import pytest

import dask
from dask.dataframe._compat import PANDAS_GE_220
from dask.dataframe.dask_expr import from_pandas, new_collection
from dask.dataframe.dask_expr._expr import Assign, Blockwise, Filter
from dask.dataframe.dask_expr._reductions import NFirst, NLast
from dask.dataframe.dask_expr._repartition import RepartitionToFewer
from dask.dataframe.dask_expr._shuffle import (
    BaseSetIndexSortValues,
    TaskShuffle,
    divisions_lru,
)
from dask.dataframe.dask_expr.io import FromPandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq, xfail_gpu
from dask.dataframe.utils import pyarrow_strings_enabled

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    return pd.DataFrame({"x": list(range(20)) * 5, "y": range(100)})


@pytest.fixture
def df(pdf):
    return from_pandas(pdf, npartitions=10)


@pytest.mark.parametrize("ignore_index", [True, False])
@pytest.mark.parametrize("npartitions", [3, 6])
def test_disk_shuffle(ignore_index, npartitions, df):
    df2 = df.shuffle(
        "x",
        shuffle_method="disk",
        npartitions=npartitions,
        ignore_index=ignore_index,
    )

    # Check that the output partition count is correct
    assert df2.npartitions == (npartitions or df.npartitions)

    # Check the computed (re-ordered) result
    assert_eq(df, df2, check_index=not ignore_index, check_divisions=False)

    # Check that df was really partitioned by "x".
    # If any values of "x" can be found in multiple
    # partitions, this will fail
    df3 = df2["x"].map_partitions(lambda x: x.drop_duplicates())
    assert sorted(df3.compute().values) == list(range(20))

    # Check `partitions` after shuffle
    a = df2.partitions[1]
    b = df.shuffle(
        "y",
        shuffle_method="disk",
        npartitions=npartitions,
        ignore_index=ignore_index,
    ).partitions[1]
    assert set(a["x"].compute().values.tolist()).issubset(
        b["y"].compute().values.tolist()
    )

    # Check for culling
    assert len(a.optimize().dask) < len(df2.optimize().dask)


@pytest.mark.parametrize("ignore_index", [True, False])
@pytest.mark.parametrize("npartitions", [8, 12])
@pytest.mark.parametrize("max_branch", [32, 6])
def test_task_shuffle(ignore_index, npartitions, max_branch, df):
    df2 = df.shuffle(
        "x",
        shuffle_method="tasks",
        npartitions=npartitions,
        ignore_index=ignore_index,
        max_branch=max_branch,
    )

    # Check that the output partition count is correct
    assert df2.npartitions == (npartitions or df.npartitions)

    # Check the computed (re-ordered) result
    assert_eq(df, df2, check_index=not ignore_index, check_divisions=False)

    # Check that df was really partitioned by "x".
    # If any values of "x" can be found in multiple
    # partitions, this will fail
    df3 = df2["x"].map_partitions(lambda x: x.drop_duplicates())
    assert sorted(df3.compute().values) == list(range(20))

    # Check `partitions` after shuffle
    a = df2.partitions[1]
    b = df.shuffle(
        "y",
        shuffle_method="tasks",
        npartitions=npartitions,
        ignore_index=ignore_index,
    ).partitions[1]
    assert set(a["x"].compute().values.tolist()).issubset(
        b["y"].compute().values.tolist()
    )

    # Check for culling
    assert len(a.optimize().dask) < len(df2.optimize().dask)


@pytest.mark.parametrize("npartitions", [3, 12])
@pytest.mark.parametrize("max_branch", [32, 8])
def test_task_shuffle_index(npartitions, max_branch, pdf):
    pdf = pdf.set_index("x")
    df = from_pandas(pdf, 10)

    df2 = df.shuffle(
        "x",
        shuffle_method="tasks",
        npartitions=npartitions,
        max_branch=max_branch,
    )

    # Check that the output partition count is correct
    assert df2.npartitions == (npartitions or df.npartitions)

    # Check the computed (re-ordered) result
    assert_eq(df, df2, check_divisions=False)

    # Check that df was really partitioned by "x".
    # If any values of "x" can be found in multiple
    # partitions, this will fail
    df3 = df2.index.map_partitions(lambda x: x.drop_duplicates())
    assert sorted(df3.compute().values) == list(range(20))


def test_shuffle_str_column_not_in_dataframe(df):
    with pytest.raises(
        KeyError,
        match="Cannot shuffle on",
    ) as execinfo:
        df.shuffle(on="z")
    assert "z" in str(execinfo.value)


def test_shuffle_mixed_list_column_not_in_dataframe(df):
    with pytest.raises(
        KeyError,
        match="Cannot shuffle on",
    ) as execinfo:
        df.shuffle(["x", "z"])
    assert "z" in str(execinfo.value)
    assert "x" not in str(execinfo.value)


def test_shuffle_list_column_not_in_dataframe(df):
    with pytest.raises(KeyError, match=r"Cannot shuffle on") as excinfo:
        df.shuffle(["zz", "z"])
    assert "z" in str(excinfo.value)
    assert "zz" in str(excinfo.value)


def test_shuffle_column_columns(df):
    df.shuffle(df.columns[-1])


def test_shuffle_column_projection(df):
    df2 = df.shuffle("x")[["x"]].simplify()

    assert "y" not in df2.expr.operands[0].columns


def test_shuffle_reductions(df):
    assert df.shuffle("x").sum().simplify()._name == df.sum()._name


@pytest.mark.xfail(reason="Shuffle can't see the reduction through the Projection")
def test_shuffle_reductions_after_projection(df):
    assert df.shuffle("x").y.sum().simplify()._name == df.y.sum()._name


@pytest.mark.parametrize("ascending", [True, False])
@pytest.mark.parametrize("by", ["a", "b", ["a", "b"]])
@pytest.mark.parametrize("nelem", [10, 500])
def test_sort_values_(nelem, by, ascending):
    np.random.seed(0)
    df = pd.DataFrame()
    df["a"] = np.ascontiguousarray(np.arange(nelem)[::-1])
    df["b"] = np.arange(100, nelem + 100)
    ddf = from_pandas(df, npartitions=10)

    # run on single-threaded scheduler for debugging purposes
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(by=by, ascending=ascending)
    expect = df.sort_values(by=by, ascending=ascending)
    assert_eq(got, expect, check_index=False, sort_results=False)


@pytest.mark.parametrize("ascending", [True, False, [False, True], [True, False]])
@pytest.mark.parametrize("by", [["a", "b"], ["b", "a"]])
@pytest.mark.parametrize("nelem", [10, 500])
def test_sort_values_single_partition(nelem, by, ascending):
    np.random.seed(0)
    df = pd.DataFrame()
    df["a"] = np.ascontiguousarray(np.arange(nelem)[::-1])
    df["b"] = np.arange(100, nelem + 100)
    ddf = from_pandas(df, npartitions=1)

    # run on single-threaded scheduler for debugging purposes
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(by=by, ascending=ascending)
    expect = df.sort_values(by=by, ascending=ascending)
    assert_eq(got, expect, check_index=False)


@pytest.mark.parametrize("partition_size", [128e6, 128e5])
@pytest.mark.parametrize("upsample", [1.0, 2.0])
def test_set_index(df, pdf, upsample, partition_size):
    assert_eq(
        df.set_index("x", upsample=upsample, partition_size=partition_size),
        pdf.set_index("x"),
    )
    assert_eq(
        df.set_index(df.x, upsample=upsample, partition_size=partition_size),
        pdf.set_index(pdf.x),
    )

    with pytest.raises(NotImplementedError, match="does not yet support"):
        df.set_index(df)


def test_set_index_sorted(pdf):
    pdf = pdf.sort_values(by="y", ignore_index=True)
    pdf["z"] = pdf["x"]
    df = from_pandas(pdf, npartitions=10)
    q = df.set_index("y", sorted=True)
    assert_eq(q, pdf.set_index("y"))

    q = df.set_index("y", sorted=True)["x"]
    assert_eq(q, pdf.set_index("y")["x"])

    with pytest.raises(NotImplementedError, match="not yet support multi-indexes"):
        df.set_index(["y", "z"], sorted=True)

    with pytest.raises(TypeError, match="not supported by set_index"):
        df.set_index([df["y"], df["x"]], sorted=True)

    with pytest.raises(NotImplementedError, match="does not yet support"):
        df.set_index(["y", "z"], sorted=False)


def test_set_index_pre_sorted(pdf):
    pdf = pdf.sort_values(by="y", ignore_index=True)
    pdf["z"] = pdf["x"]
    df = from_pandas(pdf, npartitions=10)
    q = df.set_index("y")
    assert_eq(q, pdf.set_index("y"))
    result = q.optimize(fuse=False).expr
    assert all(isinstance(expr, (Blockwise, FromPandas)) for expr in result.walk())
    q = df.set_index(df.y)
    assert_eq(q, pdf.set_index(pdf.y))
    result = q.optimize(fuse=False).expr
    assert all(isinstance(expr, (Blockwise, FromPandas)) for expr in result.walk())

    q = df.set_index("y")["x"].optimize(fuse=False)
    expected = df[["x", "y"]].set_index("y")["x"].optimize(fuse=False)
    assert q._name == expected._name


@pytest.mark.parametrize("drop", (True, False))
@pytest.mark.parametrize("append", (True, False))
def test_set_index_no_sort(drop, append):
    df = pd.DataFrame({"col1": [2, 4, 1, 3, 5], "col2": [1, 2, 3, 4, 5]})
    ddf = from_pandas(df, npartitions=2)

    assert ddf.npartitions > 1

    # Default is sort=True
    # Index in ddf will be same values, but sorted
    df_result = df.set_index("col1")
    ddf_result = ddf.set_index("col1")
    assert ddf_result.known_divisions
    assert_eq(ddf_result, df_result.sort_index(), sort_results=False)

    # Unknown divisions and index remains unsorted when sort is False
    # and thus equal to pandas set_index, adding extra kwargs also supported by
    # pandas set_index to ensure they're forwarded.
    df_result = df.set_index("col1", drop=drop, append=append)
    ddf_result = ddf.set_index("col1", sort=False, drop=drop, append=append)
    assert not ddf_result.known_divisions
    assert_eq(ddf_result, df_result, sort_results=False)


def test_set_index_repartition(df, pdf):
    result = df.set_index("x", npartitions=2)
    assert result.npartitions == 2
    assert result.optimize(fuse=False).npartitions == 2
    assert_eq(result, pdf.set_index("x"))


def test_set_index_simplify(df, pdf):
    q = df.set_index("x")["y"].optimize(fuse=False)
    expected = df[["x", "y"]].set_index("x")["y"].optimize(fuse=False)
    assert q._name == expected._name

    q = df.set_index(df.x)["y"].optimize(fuse=False)
    expected = df[["y"]].set_index(df.x)["y"].optimize(fuse=False)
    assert q._name == expected._name


def test_set_index_numeric_columns():
    pdf = pd.DataFrame(
        {
            0: list("ABAABBABAA"),
            1: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            2: [1, 2, 3, 2, 1, 3, 2, 4, 2, 3],
        }
    )
    ddf = from_pandas(pdf, 3)
    assert_eq(ddf.set_index(0), pdf.set_index(0))


def test_set_index_without_sort(df, pdf):
    result = df.set_index("y", sort=False)
    assert_eq(result, pdf.set_index("y"))

    result = result.optimize(fuse=False)
    assert all(isinstance(ex, (FromPandas, Blockwise)) for ex in result.walk())

    result = df.set_index(df.y, sort=False)
    assert_eq(result, pdf.set_index(pdf.y))

    with pytest.raises(ValueError, match="Specifying npartitions with sort=False"):
        df.set_index("y", sort=False, npartitions=20)

    with pytest.raises(
        ValueError, match="Specifying npartitions with sort=False or sorted=True"
    ):
        df.set_index("y", sorted=True, npartitions=20)


@pytest.mark.parametrize("shuffle", [None, "tasks"])
def test_sort_values(df, pdf, shuffle):
    assert_eq(df.sort_values("x", shuffle_method=shuffle), pdf.sort_values("x"))
    assert_eq(
        df.sort_values("x", shuffle_method=shuffle, npartitions=2), pdf.sort_values("x")
    )
    pdf.iloc[5, 0] = -10
    df = from_pandas(pdf, npartitions=10)
    assert_eq(
        df.sort_values("x", shuffle_method=shuffle, upsample=2.0), pdf.sort_values("x")
    )

    assert_eq(
        df.sort_values(by=["x", "y"], shuffle_method=shuffle, ascending=[True, True]),
        pdf.sort_values(by=["x", "y"], ascending=[True, True]),
    )
    with pytest.raises(NotImplementedError, match="sorting by named columns"):
        df.sort_values(by=1, shuffle_method=shuffle)

    with pytest.raises(ValueError, match="must be either 'first' or 'last'"):
        df.sort_values(by="x", shuffle_method=shuffle, na_position="bla")


@pytest.mark.parametrize("shuffle", [None, "tasks"])
def test_sort_values_temporary_column_dropped(shuffle):
    pdf = pd.DataFrame(
        {"x": range(10), "y": [1, 2, 3, 4, 5] * 2, "z": ["cat", "dog"] * 5}
    )
    df = from_pandas(pdf, npartitions=2)
    _sorted = df.sort_values(["z"], shuffle_method=shuffle)
    result = _sorted.compute()
    assert "_partitions" not in result.columns


def test_sort_values_optimize(df, pdf):
    q = df.sort_values("x")["y"].optimize(fuse=False)
    expected = df[["x", "y"]].sort_values("x")["y"].optimize(fuse=False)
    assert q._name == expected._name

    q = df.sort_values("x")["x"].optimize(fuse=False)
    expected = df[["x"]].sort_values("x")["x"].optimize(fuse=False)
    assert q._name == expected._name


def test_set_index_single_partition(pdf):
    df = from_pandas(pdf, npartitions=1)
    assert_eq(df.set_index("x"), pdf.set_index("x"))


def test_set_index_list(df, pdf):
    assert_eq(df.set_index(["x"]), pdf.set_index(["x"]))


def test_sort_values_descending(df, pdf):
    assert_eq(
        df.sort_values(by="y", ascending=False),
        pdf.sort_values(by="y", ascending=False),
        sort_results=False,
    )


def test_sort_head_nlargest(df, pdf):
    a = df.sort_values("x", ascending=False).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x"], ascending=False)
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values("x", ascending=True).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x"], ascending=True)
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values("x", ascending=[False]).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x"], ascending=[False])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values("x", ascending=[True]).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x"], ascending=[True])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x"], ascending=[False]).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x"], ascending=[False])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x"], ascending=[True]).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x"], ascending=[True])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x", "y"], ascending=[False, False]).head(10, compute=False)
    b = NFirst(df, 10, _columns=["x", "y"], ascending=[False, False])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x", "y"], ascending=[True, True]).head(10, compute=False).expr
    b = NFirst(df, 10, _columns=["x", "y"], ascending=[True, True])
    assert a.optimize()._name == b.optimize()._name


def test_sort_tail_nsmallest(df, pdf):
    a = df.sort_values("x", ascending=False).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x"], ascending=False)
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values("x", ascending=True).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x"], ascending=True)
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values("x", ascending=[False]).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x"], ascending=[False])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values("x", ascending=[True]).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x"], ascending=[True])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x"], ascending=[False]).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x"], ascending=[False])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x"], ascending=[True]).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x"], ascending=[True])
    assert a.optimize()._name == b.optimize()._name

    with pytest.raises(ValueError, match="Length of ascending"):
        df.sort_values(["x", "y"], ascending=[False]).tail(10, compute=False)

    a = df.sort_values(["x", "y"], ascending=[True, True]).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x", "y"], ascending=[True, True])
    assert a.optimize()._name == b.optimize()._name

    a = df.sort_values(["x", "y"], ascending=[False, False]).tail(10, compute=False)
    b = NLast(df, 10, _columns=["x", "y"], ascending=[False, False])
    assert a.optimize()._name == b.optimize()._name


@pytest.mark.parametrize(
    "ascending",
    [
        pytest.param([True, False], id="[True, False]"),
        pytest.param([False, True], id="[False, True]"),
    ],
)
@pytest.mark.parametrize("npartitions", [1, 3])
def test_sort_values_conflicting_ascending_head_tail(pdf, ascending, npartitions):
    divisions_lru.data = OrderedDict()

    df = from_pandas(pdf, npartitions=npartitions)

    a = df.sort_values(by=["x", "y"], ascending=ascending).head(10, compute=False)
    b = new_collection(NFirst(df, _columns=["x", "y"], n=10, ascending=ascending))
    assert a.expr.optimize()._name == b.expr.optimize()._name
    assert len(divisions_lru) == 0
    assert_eq(
        a.compute(),
        pdf.sort_values(by=["x", "y"], ascending=ascending).head(10),
    )

    a = df.sort_values(by=["x", "y"], ascending=ascending).tail(10, compute=False)
    b = new_collection(NLast(df, _columns=["x", "y"], n=10, ascending=ascending))
    assert a.expr.optimize()._name == b.expr.optimize()._name
    assert len(divisions_lru) == 0
    assert_eq(
        a.compute(),
        pdf.sort_values(by=["x", "y"], ascending=ascending).tail(10),
    )


@xfail_gpu("cudf udf support")
def test_sort_head_nlargest_string(pdf):
    pdf["z"] = "a" + pdf.x.map(str)
    df = from_pandas(pdf, npartitions=5)
    a = df.sort_values("z", ascending=False).head(10, compute=False)
    assert_eq(a, pdf.sort_values("z", ascending=False).head(10))

    a = df.sort_values("z", ascending=True).head(10, compute=False)
    assert_eq(a, pdf.sort_values("z", ascending=True).head(10))

    a = df.sort_values("z", ascending=False).tail(10, compute=False)
    assert_eq(a, pdf.sort_values("z", ascending=False).tail(10))

    a = df.sort_values("z", ascending=True).tail(10, compute=False)
    assert_eq(a, pdf.sort_values("z", ascending=True).tail(10))


def test_set_index_head_nlargest(df, pdf):
    a = df.set_index("x").head(10, compute=False)
    b = new_collection(NFirst(df, 10, _columns="x", ascending=True)).set_index("x")
    assert a.optimize()._name == b.optimize()._name

    a = df.set_index("x").tail(10, compute=False)
    b = new_collection(NLast(df, 10, _columns="x", ascending=True)).set_index("x")
    assert a.optimize()._name == b.optimize()._name

    # These still work, even if we haven't optimized them yet
    df.set_index(df.x).head(3)
    # df.set_index([df.x, df.y]).head(3)


@pytest.mark.skipif(
    not pyarrow_strings_enabled() or not PANDAS_GE_220,
    reason="doesn't work without arrow",
)
@xfail_gpu("cudf udf support")
def test_set_index_head_nlargest_string(pdf):
    pdf["z"] = "a" + pdf.x.map(str)
    df = from_pandas(pdf, npartitions=5)
    print(df.dtypes)

    a = df.set_index("z").head(10, compute=False)
    assert_eq(a, pdf.set_index("z").sort_index().head(10))

    a = df.set_index("z").tail(10, compute=False)
    assert_eq(a, pdf.set_index("z").sort_index().tail(10))


def test_filter_sort(df):
    a = df.sort_values("x")
    a = a[a.y > 40]

    b = df[df.y > 40]
    b = b.sort_values("x")

    assert a.optimize()._name == b.optimize()._name


def test_sort_values_add():
    pdf = pd.DataFrame({"x": [1, 2, 3, 0, 1, 2, 4, 5], "y": 1})
    df = from_pandas(pdf, npartitions=2, sort=False)
    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        df = df.sort_values("x")
        df["z"] = df.x + df.y
        pdf = pdf.sort_values("x")
        pdf["z"] = pdf.x + pdf.y
        assert_eq(df, pdf, sort_results=False)


@pytest.mark.parametrize("null_value", [None, pd.NaT, pd.NA])
def test_index_nulls(null_value):
    "Setting the index with some non-numeric null raises error"
    df = pd.DataFrame(
        {"numeric": [1, 2, 3, 4], "non_numeric": ["foo", "bar", "foo", "bar"]}
    )
    ddf = from_pandas(df, npartitions=2)
    with pytest.raises(NotImplementedError, match="presence of nulls"):
        with pytest.warns(UserWarning):
            ddf.set_index(
                ddf["non_numeric"].map({"foo": "foo", "bar": null_value})
            ).compute()


@pytest.mark.parametrize("freq", ["16h", "-16h"])
def test_set_index_with_dask_dt_index(freq):
    values = {
        "x": [1, 2, 3, 4] * 3,
        "y": [10, 20, 30] * 4,
        "name": ["Alice", "Bob"] * 6,
    }
    date_index = pd.date_range(
        start="2022-02-22", freq=freq, periods=12
    ) - pd.Timedelta(seconds=30)
    df = pd.DataFrame(values, index=date_index)
    ddf = from_pandas(df, npartitions=3, sort=False)
    # specify a different date index entirely
    day_index = ddf.index.dt.floor("D")
    result = ddf.set_index(day_index)
    assert_eq(result, pd.DataFrame(values, index=date_index.floor("D")))


def test_set_index_sort_values_shuffle_options(df, pdf):
    q = df.set_index("x", shuffle_method="tasks", max_branch=10)
    shuffle = list(q.optimize().find_operations(TaskShuffle))[0]
    assert shuffle.options == {"max_branch": 10}
    assert_eq(q, pdf.set_index("x"))

    q = df.sort_values("x", shuffle_method="tasks", max_branch=10)
    sorter = list(q.optimize().find_operations(TaskShuffle))[0]
    assert sorter.options == {"max_branch": 10}
    assert_eq(q, pdf)


def test_set_index_predicate_pushdown(df, pdf):
    pdf = pdf.set_index("x")
    query = df.set_index("x")
    result = query[query.y > 5]
    expected = pdf[pdf.y > 5]
    assert_eq(result, expected)
    expected_query = df[df.y > 5].set_index("x").optimize()
    assert expected_query._name == result.optimize()._name

    result = query[query.index.to_series() > 5]
    assert_eq(result, pdf[pdf.index > 5])

    result = query[(query.index.to_series() > 5) & (query.y > -1)]
    assert_eq(result, pdf[(pdf.index > 5) & (pdf.y > -1)])


def test_set_index_with_explicit_divisions():
    pdf = pd.DataFrame({"x": [4, 1, 2, 5]}, index=[10, 20, 30, 40])

    df = from_pandas(pdf, npartitions=2)

    result = df.set_index("x", divisions=[1, 3, 5])
    assert result.divisions == (1, 3, 5)
    assert_eq(result, pdf.set_index("x"))


def test_set_index_npartitions_changes(pdf):
    df = from_pandas(pdf, npartitions=30)
    result = df.set_index("x", shuffle_method="disk")
    assert result.npartitions == result.optimize().npartitions
    assert_eq(result, pdf.set_index("x"))


def test_set_index_sorted_divisions(df):
    with pytest.raises(ValueError, match="must be the same length"):
        df.set_index("x", divisions=(1, 2, 3), sorted=True)


def test_set_index_sort_values_one_partition(pdf):
    divisions_lru.data = OrderedDict()
    df = from_pandas(pdf, sort=False)
    query = df.sort_values("x").optimize(fuse=False)
    assert query.divisions == (0, 99)
    assert_eq(pdf.sort_values("x"), query, sort_results=False)
    assert len(divisions_lru) == 0

    df = from_pandas(pdf, sort=False)
    query = df.set_index("x").optimize(fuse=False)
    assert query.divisions == (None, None)
    assert_eq(pdf.set_index("x"), query)
    assert len(divisions_lru) == 0

    df = from_pandas(pdf, sort=False, npartitions=2)
    query = df.set_index("x", npartitions=1).optimize(fuse=False)
    assert query.divisions == (None, None)
    assert_eq(pdf.set_index("x"), query)
    assert len(divisions_lru) == 0
    assert len(list(query.expr.find_operations(RepartitionToFewer))) > 0


def test_set_index_triggers_calc_when_accessing_divisions(pdf, df):
    divisions_lru.data = OrderedDict()
    query = df.set_index("x")
    assert len(divisions_lru.data) == 0
    divisions = query.divisions  # noqa: F841
    assert len(divisions_lru.data) == 1


def test_shuffle(df, pdf):
    result = df.shuffle(df.x)
    assert result.npartitions == df.npartitions
    assert_eq(result, pdf)

    result = df.shuffle(df[["x"]])
    assert result.npartitions == df.npartitions
    assert_eq(result, pdf)

    result = df[["y"]].shuffle(df[["x"]])
    assert result.npartitions == df.npartitions
    assert_eq(result, pdf[["y"]])

    with pytest.raises(TypeError, match="index must be aligned"):
        df.shuffle(df.x.repartition(npartitions=2))

    result = df.shuffle(df.x, npartitions=2)
    assert result.npartitions == 2
    assert_eq(result, pdf)


def test_empty_partitions():
    # See https://github.com/dask/dask/issues/2408
    df = pd.DataFrame({"a": list(range(10))})
    df["b"] = df["a"] % 3
    df["c"] = df["b"].astype(str)

    ddf = from_pandas(df, npartitions=3)
    ddf = ddf.set_index("b")
    ddf = ddf.repartition(npartitions=3)
    ddf.get_partition(0).compute()
    assert_eq(ddf, df.set_index("b"))

    ddf = ddf.set_index("c")
    assert_eq(ddf, df.set_index("b").set_index("c"))


def test_shuffle_no_assign(df, pdf):
    result = df.shuffle(df.x)
    q = result.optimize(fuse=False)
    assert len([x for x in q.walk() if isinstance(x, Assign)]) == 0


@pytest.mark.parametrize("meth", ["shuffle", "sort_values"])
def test_shuffle_filter_pushdown(pdf, meth):
    pdf["z"] = 1
    df = from_pandas(pdf, npartitions=10)
    result = getattr(df, meth)("x")
    result = result[result.x > 5.0]
    expected = getattr(df[df.x > 5.0], meth)("x")
    assert result.simplify()._name == expected._name

    result = getattr(df, meth)("x")
    result = result[result.x > 5.0][["x", "y"]]
    expected = df[["x", "y"]]
    expected = getattr(expected[expected.x > 5.0], meth)("x")
    assert result.simplify()._name == expected.simplify()._name

    result = getattr(df, meth)("x")[["x", "y"]]
    result = result[result.x > 5.0]
    expected = df[["x", "y"]]
    expected = getattr(expected[expected.x > 5.0], meth)("x")
    assert result.simplify()._name == expected.simplify()._name


@pytest.mark.parametrize("meth", ["set_index", "sort_values"])
def test_sort_values_avoid_overeager_filter_pushdown(meth):
    pdf1 = pd.DataFrame({"a": [4, 2, 3], "b": [1, 2, 3]})
    df = from_pandas(pdf1, npartitions=2)
    df = getattr(df, meth)("a")
    df = df[df.b > 2] + df.b.sum()
    result = df.simplify()
    assert isinstance(result.expr.left, Filter)
    assert isinstance(result.expr.left.frame, BaseSetIndexSortValues)


def test_set_index_filter_pushdown():
    pdf = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6, 7, 8] * 10, "y": 1, "z": 2})
    df = from_pandas(pdf, npartitions=10)
    result = df.set_index("x")
    result = result[result.y == 1]
    expected = df[df.y == 1].set_index("x")
    assert result.simplify()._name == expected._name

    result = df.set_index("x")
    result = result[result.y == 1][["y"]]
    expected = df[["x", "y"]]
    expected = expected[expected.y == 1].set_index("x")
    assert result.simplify()._name == expected.simplify()._name

    result = df.set_index("x")[["y"]]
    result = result[result.y == 1]
    expected = df[["x", "y"]]
    expected = expected[expected.y == 1].set_index("x")
    assert result.simplify()._name == expected.simplify()._name


def test_shuffle_index_shuffle(df):
    with pytest.raises(TypeError, match="Must shuffle on either "):
        df.shuffle()
    with pytest.raises(TypeError, match="Cannot shuffle on both"):
        df.shuffle("x", on_index=True)


def test_set_index_user_divisions_one_partition(pdf):
    df = from_pandas(pdf, npartitions=1)
    result = df.set_index("x", divisions=[0, 10, 20, 21])
    assert_eq(result, pdf.set_index("x"))


def test_set_index_divisions_npartitions(pdf):
    df = from_pandas(pdf, npartitions=1, sort=False)
    result = df.set_index("x", sort=True, npartitions=2)
    assert result.known_divisions
    assert result.optimize().known_divisions
    assert_eq(result, pdf.set_index("x"))


def test_sort_values_unordered_categorical():
    df = pd.DataFrame()
    df["a"] = np.arange(10)[::-1]
    df["b"] = df["a"].astype("category")
    ddf = from_pandas(df, npartitions=2)
    result = ddf.sort_values(by="b")
    expected = df.sort_values(by="b")
    assert_eq(result, expected, check_index=False)


def test_set_index_before_assign(df, pdf):
    result = df.set_index("x")
    result["z"] = result.y + 1
    expected = pdf.set_index("x")
    expected["z"] = expected.y + 1
    assert_eq(result["z"], expected["z"])


def test_set_index_shuffle_afterwards(pdf):
    ddf = from_pandas(pdf, npartitions=1)
    ddf = ddf.set_index("y", sort=True, divisions=[0, 10, 20, 100], shuffle="tasks")
    result = ddf.reset_index().y.unique()
    expected = pd.Series(pdf.y.unique(), name="y")
    assert_eq(result, expected, check_index=False)
