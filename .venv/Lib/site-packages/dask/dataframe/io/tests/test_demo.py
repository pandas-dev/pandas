from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import dask
import dask.dataframe as dd
from dask.dataframe._compat import PANDAS_GE_220, tm
from dask.dataframe.utils import assert_eq, get_string_dtype

ME = "ME" if PANDAS_GE_220 else "M"


def test_make_timeseries():
    df = dd.demo.make_timeseries(
        "2000",
        "2015",
        {"A": float, "B": int, "C": str},
        freq="2D",
        partition_freq=f"6{ME}",
    )

    assert df.divisions[0] == pd.Timestamp("2000-01-31")
    assert df.divisions[-1] == pd.Timestamp("2014-07-31")
    tm.assert_index_equal(df.columns, pd.Index(["A", "B", "C"]))
    assert df["A"].head().dtype == float
    assert np.issubdtype(df["B"].head(), np.integer)
    assert df["C"].head().dtype == get_string_dtype()
    assert df.index.name == "timestamp"
    assert df.head().index.name == df.index.name
    assert df.divisions == tuple(pd.date_range(start="2000", end="2015", freq=f"6{ME}"))

    tm.assert_frame_equal(df.head(), df.head())

    a = dd.demo.make_timeseries(
        "2000",
        "2015",
        {"A": float, "B": int, "C": str},
        freq="2D",
        partition_freq=f"6{ME}",
        seed=123,
    )
    b = dd.demo.make_timeseries(
        "2000",
        "2015",
        {"A": float, "B": int, "C": str},
        freq="2D",
        partition_freq=f"6{ME}",
        seed=123,
    )
    c = dd.demo.make_timeseries(
        "2000",
        "2015",
        {"A": float, "B": int, "C": str},
        freq="2D",
        partition_freq=f"6{ME}",
        seed=456,
    )
    d = dd.demo.make_timeseries(
        "2000",
        "2015",
        {"A": float, "B": int, "C": str},
        freq="2D",
        partition_freq=f"3{ME}",
        seed=123,
    )
    e = dd.demo.make_timeseries(
        "2000",
        "2015",
        {"A": float, "B": int, "C": str},
        freq="1D",
        partition_freq=f"6{ME}",
        seed=123,
    )
    tm.assert_frame_equal(a.head(), b.head())
    assert not (a.head(10) == c.head(10)).all().all()
    assert a._name == b._name
    assert a._name != c._name
    assert a._name != d._name
    assert a._name != e._name


def test_make_timeseries_no_args():
    df = dd.demo.make_timeseries()
    assert 1 < df.npartitions < 1000
    assert len(df.columns) > 1
    assert len(set(df.dtypes)) > 1


def test_no_overlaps():
    df = dd.demo.make_timeseries(
        "2000", "2001", {"A": float}, freq="3h", partition_freq=f"3{ME}"
    )

    assert all(
        df.get_partition(i).index.max().compute()
        < df.get_partition(i + 1).index.min().compute()
        for i in range(df.npartitions - 2)
    )


def test_make_timeseries_keywords():
    df = dd.demo.make_timeseries(
        "2000",
        "2001",
        {"A": int, "B": int, "C": str},
        freq="1D",
        partition_freq=f"6{ME}",
        A_lam=1000000,
        B_lam=2,
    )
    a_cardinality = df.A.nunique()
    b_cardinality = df.B.nunique()

    aa, bb = dask.compute(a_cardinality, b_cardinality, scheduler="single-threaded")

    assert 100 < aa <= 10000000
    assert 1 < bb <= 100


def test_make_timeseries_fancy_keywords():
    df = dd.demo.make_timeseries(
        "2000",
        "2001",
        {"A_B": int, "B_": int, "C": str},
        freq="1D",
        partition_freq=f"6{ME}",
        A_B_lam=1000000,
        B__lam=2,
    )
    a_cardinality = df.A_B.nunique()
    b_cardinality = df.B_.nunique()

    aa, bb = dask.compute(a_cardinality, b_cardinality, scheduler="single-threaded")

    assert 100 < aa <= 10000000
    assert 1 < bb <= 100


def test_make_timeseries_getitem_compute():
    # See https://github.com/dask/dask/issues/7692

    df = dd.demo.make_timeseries()
    df2 = df[df.y > 0]
    df3 = df2.compute()
    assert df3["y"].min() > 0
    assert list(df.columns) == list(df3.columns)


def test_make_timeseries_column_projection():
    ddf = dd.demo.make_timeseries(
        "2001", "2002", freq="1D", partition_freq=f"3{ME}", seed=42
    )

    assert_eq(ddf[["x"]].compute(), ddf.compute()[["x"]])
    assert_eq(
        ddf.groupby("name").aggregate({"x": "sum", "y": "max"}).compute(),
        ddf.compute().groupby("name").aggregate({"x": "sum", "y": "max"}),
    )


@pytest.mark.parametrize("seed", [None, 42])
def test_with_spec(seed):
    """Make a dataset with default random columns"""
    from dask.dataframe.io.demo import DatasetSpec, with_spec

    spec = DatasetSpec(nrecords=10, npartitions=2)
    ddf = with_spec(spec, seed=seed)
    assert isinstance(ddf, dd.DataFrame)
    assert ddf.npartitions == 2
    assert ddf.columns.tolist() == ["i1", "f1", "c1", "s1"]
    assert ddf["i1"].dtype == "int64"
    assert ddf["f1"].dtype == float
    assert ddf["c1"].dtype.name == "category"
    assert ddf["s1"].dtype == get_string_dtype()
    res = ddf.compute()
    assert len(res) == 10


@pytest.mark.parametrize("seed", [None, 42])
def test_with_spec_non_default(seed):
    from dask.dataframe.io.demo import (
        ColumnSpec,
        DatasetSpec,
        RangeIndexSpec,
        with_spec,
    )

    spec = DatasetSpec(
        npartitions=3,
        nrecords=10,
        index_spec=RangeIndexSpec(dtype="int32", step=2),
        column_specs=[
            ColumnSpec(prefix="i", dtype="int32", low=1, high=100, random=True),
            ColumnSpec(prefix="f", dtype="float32", random=True),
            ColumnSpec(prefix="c", dtype="category", choices=["apple", "banana"]),
            ColumnSpec(prefix="s", dtype=str, length=15, random=True),
        ],
    )
    ddf = with_spec(spec, seed=seed)
    assert isinstance(ddf, dd.DataFrame)
    assert ddf.columns.tolist() == ["i1", "f1", "c1", "s1"]
    assert ddf.index.dtype == "int32"
    assert ddf["i1"].dtype == "int32"
    assert ddf["f1"].dtype == "float32"
    assert ddf["c1"].dtype.name == "category"
    assert ddf["s1"].dtype == get_string_dtype()
    res = ddf.compute().sort_index()
    assert len(res) == 10
    assert set(res.c1.cat.categories) == {"apple", "banana"}
    assert res.i1.min() >= 1
    assert res.i1.max() <= 100
    assert all(len(s) == 15 for s in res.s1.tolist())
    assert len(res.s1.unique()) <= 10


def test_with_spec_pyarrow():
    pytest.importorskip("pyarrow", "1.0.0", reason="pyarrow is required")
    from dask.dataframe.io.demo import ColumnSpec, DatasetSpec, with_spec

    spec = DatasetSpec(
        npartitions=1,
        nrecords=10,
        column_specs=[
            ColumnSpec(dtype="string[pyarrow]", length=10, random=True),
        ],
    )
    ddf = with_spec(spec, seed=42)
    assert isinstance(ddf, dd.DataFrame)
    assert ddf.columns.tolist() == ["string_pyarrow1"]
    assert ddf["string_pyarrow1"].dtype == "string[pyarrow]"
    res = ddf.compute()
    assert res["string_pyarrow1"].dtype == "string[pyarrow]"
    assert all(len(s) == 10 for s in res["string_pyarrow1"].tolist())


@pytest.mark.parametrize("seed", [None, 42])
def test_same_prefix_col_numbering(seed):
    from dask.dataframe.io.demo import ColumnSpec, DatasetSpec, with_spec

    spec = DatasetSpec(
        npartitions=1,
        nrecords=5,
        column_specs=[
            ColumnSpec(dtype=int),
            ColumnSpec(dtype=int),
            ColumnSpec(dtype=int),
            ColumnSpec(dtype=int),
        ],
    )
    ddf = with_spec(spec, seed=seed)
    assert ddf.columns.tolist() == ["int1", "int2", "int3", "int4"]


def test_with_spec_category_nunique():
    from dask.dataframe.io.demo import ColumnSpec, DatasetSpec, with_spec

    spec = DatasetSpec(
        npartitions=1,
        nrecords=20,
        column_specs=[
            ColumnSpec(dtype="category", nunique=10),
        ],
    )
    ddf = with_spec(spec, seed=42)
    res = ddf.compute()
    assert res.category1.cat.categories.tolist() == [
        "01",
        "02",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
    ]


@pytest.mark.parametrize("seed", [None, 42])
def test_with_spec_default_integer(seed):
    from dask.dataframe.io.demo import ColumnSpec, DatasetSpec, with_spec

    spec = DatasetSpec(
        npartitions=1,
        nrecords=5,
        column_specs=[
            ColumnSpec(dtype=int),
            ColumnSpec(dtype=int),
            ColumnSpec(dtype=int),
            ColumnSpec(dtype=int),
        ],
    )
    ddf = with_spec(spec, seed=seed)
    res = ddf.compute()
    for col in res.columns:
        assert 500 < res[col].min() < 1500
        assert 500 < res[col].max() < 1500


def test_with_spec_integer_method():
    from dask.dataframe.io.demo import ColumnSpec, DatasetSpec, with_spec

    spec = DatasetSpec(
        npartitions=1,
        nrecords=5,
        column_specs=[
            ColumnSpec(prefix="pois", dtype=int, method="poisson"),
            ColumnSpec(prefix="norm", dtype=int, method="normal"),
            ColumnSpec(prefix="unif", dtype=int, method="uniform"),
            ColumnSpec(prefix="binom", dtype=int, method="binomial", args=(100, 0.4)),
            ColumnSpec(prefix="choice", dtype=int, method="choice", args=(10,)),
            ColumnSpec(prefix="rand", dtype=int, random=True, low=0, high=10),
            ColumnSpec(prefix="rand", dtype=int, random=True),
        ],
    )
    ddf = with_spec(spec, seed=42)
    res = ddf.compute()
    assert res["pois1"].tolist() == [1002, 985, 947, 1003, 1017]
    assert res["norm1"].tolist() == [-1097, -276, 853, 272, 784]
    assert res["unif1"].tolist() == [772, 972, 798, 393, 656]
    assert res["binom1"].tolist() == [34, 46, 38, 37, 43]
    assert res["choice1"].tolist() == [0, 3, 1, 6, 6]
    assert res["rand1"].tolist() == [4, 6, 9, 4, 5]
    assert res["rand2"].tolist() == [883, 104, 192, 648, 256]


def test_with_spec_datetime_index():
    from dask.dataframe.io.demo import (
        ColumnSpec,
        DatasetSpec,
        DatetimeIndexSpec,
        with_spec,
    )

    spec = DatasetSpec(
        nrecords=10,
        index_spec=DatetimeIndexSpec(
            dtype="datetime64[ns]", freq="1h", start="2023-01-02", partition_freq="1D"
        ),
        column_specs=[ColumnSpec(dtype=int)],
    )
    ddf = with_spec(spec, seed=42)
    assert ddf.index.dtype == "datetime64[ns]"
    res = ddf.compute()
    assert len(res) == 10
