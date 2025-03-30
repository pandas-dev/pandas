from __future__ import annotations

import pytest

import dask.dataframe as dd
from dask.dataframe.dask_expr import from_pandas, optimize
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(100)})
    pdf["y"] = pdf.x * 10.0
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=10)


def test_simple(df):
    out = (df["x"] + df["y"]) - 1
    unfused = optimize(out, fuse=False)
    fused = optimize(out, fuse=True)

    # Should only get one task per partition
    # from_pandas is not fused together
    assert len(fused.dask) == df.npartitions + 10
    assert_eq(fused, unfused)


def test_with_non_fusable_on_top(df):
    out = (df["x"] + df["y"] - 1).sum()
    unfused = optimize(out, fuse=False)
    fused = optimize(out, fuse=True)

    assert len(fused.dask) < len(unfused.dask)
    assert_eq(fused, unfused)

    # Check that we still get fusion
    # after a non-blockwise operation as well
    fused_2 = optimize((out + 10) - 5, fuse=True)
    assert len(fused_2.dask) == len(fused.dask) + 1  # only one more task


def test_optimize_fusion_many():
    # Test that many `Blockwise`` operations,
    # originating from various IO operations,
    # can all be fused together
    a = from_pandas(pd.DataFrame({"x": range(100), "y": range(100)}), 10)
    b = from_pandas(pd.DataFrame({"a": range(100)}), 10)

    # some generic elemwise operations
    aa = a[["x"]] + 1
    aa["a"] = a["y"] + a["x"]
    aa["b"] = aa["x"] + 2
    series_a = aa[a["x"] > 1]["b"]

    bb = b[["a"]] + 1
    bb["b"] = b["a"] + b["a"]
    series_b = bb[b["a"] > 1]["b"]

    result = (series_a + series_b) + 1
    fused = optimize(result, fuse=True)
    unfused = optimize(result, fuse=False)
    assert fused.npartitions == a.npartitions
    # from_pandas is not fused together
    assert len(fused.dask) == fused.npartitions + 20
    assert_eq(fused, unfused)


def test_optimize_fusion_repeat(df):
    # Test that we can optimize a collection
    # more than once, and fusion still works

    original = df.copy()

    # some generic elemwise operations
    df["x"] += 1
    df["z"] = df.y
    df += 2

    # repeatedly call optimize after doing new fusable things
    fused = optimize(optimize(optimize(df) + 2).x)

    # from_pandas is not fused together
    assert len(fused.dask) - 10 == fused.npartitions == original.npartitions
    assert_eq(fused, df.x + 2)


def test_optimize_fusion_broadcast(df):
    # Check fusion with broadcated reduction
    result = ((df["x"] + 1) + df["y"].sum()) + 1
    fused = optimize(result)

    assert_eq(fused, result)
    assert len(fused.dask) < len(result.dask)


def test_persist_with_fusion(df):
    # Check that fusion works after persisting
    df = (df + 2).persist()
    out = (df.y + 1).sum()
    fused = optimize(out)

    assert_eq(out, fused)
    assert len(fused.dask) < len(out.dask)


def test_fuse_broadcast_deps():
    pdf = pd.DataFrame({"a": [1, 2, 3]})
    pdf2 = pd.DataFrame({"a": [2, 3, 4]})
    pdf3 = pd.DataFrame({"a": [3, 4, 5]})
    df = from_pandas(pdf, npartitions=1)
    df2 = from_pandas(pdf2, npartitions=1)
    df3 = from_pandas(pdf3, npartitions=2)

    query = df.merge(df2).merge(df3)
    # from_pandas is not fused together
    assert len(query.optimize().__dask_graph__()) == 2 + 4
    assert_eq(query, pdf.merge(pdf2).merge(pdf3))


def test_name(df):
    out = (df["x"] + df["y"]) - 1
    fused = optimize(out, fuse=True)
    assert "getitem" in str(fused.expr)
    assert "sub" in str(fused.expr)
    assert str(fused.expr) == str(fused.expr).lower()


def test_fusion_executes_only_once():
    times_called = []
    import pandas as pd

    def test(i):
        times_called.append(i)
        return pd.DataFrame({"a": [1, 2, 3], "b": 1})

    df = dd.from_map(test, [1], meta=[("a", "i8"), ("b", "i8")])
    df = df[df.a > 1]
    df.sum().compute()
    assert len(times_called) == 1
