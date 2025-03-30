from __future__ import annotations

import multiprocessing as mp
import pickle
import random
import string
from concurrent.futures import ProcessPoolExecutor
from datetime import date, time
from decimal import Decimal
from functools import partial

import numpy as np
import pandas as pd
import pytest

import dask
import dask.dataframe as dd
from dask.base import compute_as_if_collection
from dask.dataframe._compat import PANDAS_GE_220, assert_categorical_equal, tm
from dask.dataframe.shuffle import maybe_buffered_partd, partitioning_index
from dask.dataframe.utils import assert_eq, make_meta

try:
    import pyarrow as pa
except ImportError:
    pa = None

dsk = {
    ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [1, 4, 7]}, index=[0, 1, 3]),
    ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [2, 5, 8]}, index=[5, 6, 8]),
    ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [3, 6, 9]}, index=[9, 9, 9]),
}
meta = make_meta(
    {"a": "i8", "b": "i8"}, index=pd.Index([], "i8"), parent_meta=pd.DataFrame()
)
d = dd.repartition(pd.concat(dsk.values()), [0, 4, 9, 9])
shuffle_func = lambda df, *args, **kwargs: df.shuffle(*args, **kwargs)
shuffle = lambda df, *args, **kwargs: df.shuffle(*args, **kwargs)
full = d.compute()


def test_shuffle(shuffle_method):
    s = shuffle_func(d, d.b, shuffle_method=shuffle_method)
    assert isinstance(s, dd.DataFrame)
    assert s.npartitions == d.npartitions
    assert set(s.dask).issuperset(d.dask)

    assert shuffle_func(d, d.b)._name == shuffle_func(d, d.b)._name


def test_default_partitions():
    assert shuffle(d, d.b).npartitions == d.npartitions


def test_shuffle_npartitions(shuffle_method):
    df = pd.DataFrame({"x": np.random.random(100)})
    ddf = dd.from_pandas(df, npartitions=10)
    s = shuffle(ddf, ddf.x, shuffle_method=shuffle_method, npartitions=17, max_branch=4)
    sc = s.compute()
    assert s.npartitions == 17
    assert set(s.dask).issuperset(set(ddf.dask))

    assert len(sc) == len(df)
    assert list(s.columns) == list(df.columns)
    assert set(map(tuple, sc.values.tolist())) == set(map(tuple, df.values.tolist()))


def test_shuffle_npartitions_lt_input_partitions(shuffle_method):
    df = pd.DataFrame({"x": np.random.random(100)})
    ddf = dd.from_pandas(df, npartitions=20)
    s = shuffle(ddf, ddf.x, shuffle_method=shuffle_method, npartitions=5, max_branch=2)
    sc = s.compute()
    assert s.npartitions == 5
    assert set(s.dask).issuperset(set(ddf.dask))

    assert len(sc) == len(df)
    assert list(s.columns) == list(df.columns)
    assert set(map(tuple, sc.values.tolist())) == set(map(tuple, df.values.tolist()))


def test_index_with_non_series(shuffle_method):
    from dask.dataframe.tests.test_multi import list_eq

    list_eq(
        shuffle(d, d.b, shuffle_method=shuffle_method),
        shuffle(d, "b", shuffle_method=shuffle_method),
    )


def test_index_with_dataframe(shuffle_method):
    res1 = shuffle(d, d[["b"]], shuffle_method=shuffle_method).compute()
    res2 = shuffle(d, ["b"], shuffle_method=shuffle_method).compute()
    res3 = shuffle(d, "b", shuffle_method=shuffle_method).compute()

    assert sorted(res1.values.tolist()) == sorted(res2.values.tolist())
    assert sorted(res1.values.tolist()) == sorted(res3.values.tolist())


def test_shuffle_from_one_partition_to_one_other(shuffle_method):
    df = pd.DataFrame({"x": [1, 2, 3]})
    a = dd.from_pandas(df, 1)

    for i in [1, 2]:
        b = shuffle(a, "x", npartitions=i, shuffle_method=shuffle_method)
        assert len(a.compute(scheduler="sync")) == len(b.compute(scheduler="sync"))


def test_shuffle_empty_partitions(shuffle_method):
    df = pd.DataFrame({"x": [1, 2, 3] * 10})
    ddf = dd.from_pandas(df, npartitions=3)
    s = shuffle(ddf, ddf.x, npartitions=6, shuffle_method=shuffle_method)
    parts = compute_as_if_collection(dd.DataFrame, s.dask, s.__dask_keys__())
    for p in parts:
        assert s.columns == p.columns


df2 = pd.DataFrame(
    {
        "i32": np.array([1, 2, 3] * 3, dtype="int32"),
        "f32": np.array([None, 2.5, 3.5] * 3, dtype="float32"),
        "cat": pd.Series(["a", "b", "c"] * 3).astype("category"),
        "obj": pd.Series(["d", "e", "f"] * 3),
        "bool": np.array([True, False, True] * 3),
        "dt": pd.Series(pd.date_range("20130101", periods=9)),
        "dt_tz": pd.Series(pd.date_range("20130101", periods=9, tz="US/Eastern")),
        "td": pd.Series(pd.timedelta_range("2000", periods=9)),
    }
)


def test_partitioning_index():
    res = partitioning_index(df2.i32, 3)
    assert ((res < 3) & (res >= 0)).all()
    assert len(np.unique(res)) > 1

    assert (partitioning_index(df2.i32, 3) == partitioning_index(df2.i32, 3)).all()

    res = partitioning_index(df2[["i32"]], 3)
    assert ((res < 3) & (res >= 0)).all()
    assert len(np.unique(res)) > 1

    res = partitioning_index(df2[["cat", "bool", "f32"]], 2)
    assert ((0 <= res) & (res < 2)).all()

    res = partitioning_index(df2.index, 4)
    assert ((res < 4) & (res >= 0)).all()
    assert len(np.unique(res)) > 1


def test_partitioning_index_categorical_on_values():
    df = pd.DataFrame({"a": list(string.ascii_letters), "b": [1, 2, 3, 4] * 13})
    df.a = df.a.astype("category")
    df2 = df.copy()
    df2.a = df2.a.cat.set_categories(list(reversed(df2.a.cat.categories)))

    res = partitioning_index(df.a, 5)
    res2 = partitioning_index(df2.a, 5)
    assert (res == res2).all()

    res = partitioning_index(df, 5)
    res2 = partitioning_index(df2, 5)
    assert (res == res2).all()


@pytest.mark.parametrize(
    "npartitions", [1, 4, 7, pytest.param(23, marks=pytest.mark.slow)]
)
def test_set_index_general(npartitions, shuffle_method):
    names = ["alice", "bob", "ricky"]
    df = pd.DataFrame(
        {
            "x": np.random.random(100),
            "y": np.random.random(100) // 0.2,
            "z": np.random.choice(names, 100),
        },
        index=np.random.random(100),
    )
    # Ensure extension dtypes work
    df = df.astype({"x": "Float64", "z": "string"})

    ddf = dd.from_pandas(df, npartitions=npartitions)

    assert_eq(df.set_index("x"), ddf.set_index("x", shuffle_method=shuffle_method))
    assert_eq(df.set_index("y"), ddf.set_index("y", shuffle_method=shuffle_method))
    assert_eq(df.set_index("z"), ddf.set_index("z", shuffle_method=shuffle_method))
    assert_eq(df.set_index(df.x), ddf.set_index(ddf.x, shuffle_method=shuffle_method))
    assert_eq(
        df.set_index(df.x + df.y),
        ddf.set_index(ddf.x + ddf.y, shuffle_method=shuffle_method),
    )
    assert_eq(
        df.set_index(df.x + 1), ddf.set_index(ddf.x + 1, shuffle_method=shuffle_method)
    )


@pytest.mark.parametrize(
    "string_dtype", ["string[python]", "string[pyarrow]", "object"]
)
def test_set_index_string(shuffle_method, string_dtype):
    if string_dtype == "string[pyarrow]":
        pytest.importorskip("pyarrow")
    names = ["alice", "bob", "ricky"]
    df = pd.DataFrame(
        {
            "x": np.random.random(100),
            "y": np.random.choice(names, 100),
        },
        index=np.random.random(100),
    )
    df = df.astype({"y": string_dtype})
    ddf = dd.from_pandas(df, npartitions=10)
    assert_eq(df.set_index("y"), ddf.set_index("y", shuffle_method=shuffle_method))


def test_set_index_self_index(shuffle_method):
    df = pd.DataFrame(
        {"x": np.random.random(100), "y": np.random.random(100) // 0.2},
        index=np.random.random(100),
    )

    a = dd.from_pandas(df, npartitions=4)
    b = a.set_index(a.index, shuffle_method=shuffle_method)
    assert a is b

    assert_eq(b, df.set_index(df.index))


def test_set_index_names(shuffle_method):
    if shuffle_method == "disk":
        pytest.xfail("dsk names in disk shuffle are not deterministic")

    df = pd.DataFrame(
        {"x": np.random.random(100), "y": np.random.random(100) // 0.2},
        index=np.random.random(100),
    )

    ddf = dd.from_pandas(df, npartitions=4)

    assert set(ddf.set_index("x", shuffle_method=shuffle_method).dask) == set(
        ddf.set_index("x", shuffle_method=shuffle_method).dask
    )
    assert set(ddf.set_index("x", shuffle_method=shuffle_method).dask) != set(
        ddf.set_index("y", shuffle_method=shuffle_method).dask
    )
    assert set(
        ddf.set_index("x", max_branch=4, shuffle_method=shuffle_method).dask
    ) != set(ddf.set_index("x", max_branch=3, shuffle_method=shuffle_method).dask)
    assert set(
        ddf.set_index("x", drop=True, shuffle_method=shuffle_method).dask
    ) != set(ddf.set_index("x", drop=False, shuffle_method=shuffle_method).dask)


ME = "ME" if PANDAS_GE_220 else "M"


def test_set_index_2(shuffle_method):
    df = dd.demo.make_timeseries(
        "2000",
        "2004",
        {"value": float, "name": str, "id": int},
        freq="2h",
        partition_freq=f"1{ME}",
        seed=1,
    )

    df2 = df.set_index("name", shuffle_method=shuffle_method)
    df2.value.sum().compute(scheduler="sync")


def test_set_index_3(shuffle_method):
    df = pd.DataFrame(np.random.random((10, 2)), columns=["x", "y"])
    ddf = dd.from_pandas(df, npartitions=5)

    ddf2 = ddf.set_index(
        "x", shuffle_method=shuffle_method, max_branch=2, npartitions=ddf.npartitions
    )
    df2 = df.set_index("x")
    assert_eq(df2, ddf2)
    assert ddf2.npartitions == ddf.npartitions


@pytest.mark.parametrize("drop", (True, False))
@pytest.mark.parametrize("append", (True, False))
def test_set_index_no_sort(drop, append):
    """
    GH10333 - Allow setting index on existing partitions without
    computing new divisions and repartitioning.
    """
    df = pd.DataFrame({"col1": [2, 4, 1, 3, 5], "col2": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

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


def test_shuffle_sort(shuffle_method):
    df = pd.DataFrame({"x": [1, 2, 3, 2, 1], "y": [9, 8, 7, 1, 5]})
    ddf = dd.from_pandas(df, npartitions=3)

    df2 = df.set_index("x").sort_index()
    ddf2 = ddf.set_index("x", shuffle_method=shuffle_method)

    assert_eq(ddf2, df2, sort_results=False)


def test_maybe_buffered_partd(tmp_path):
    import partd

    f = maybe_buffered_partd()
    p1 = f()
    assert isinstance(p1.partd, partd.Buffer)
    f2 = pickle.loads(pickle.dumps(f))
    assert not f2.buffer
    p2 = f2()
    assert isinstance(p2.partd, partd.File)

    f3 = maybe_buffered_partd(tempdir=tmp_path)
    p3 = f3()
    assert isinstance(p3.partd, partd.Buffer)
    contents = list(tmp_path.iterdir())
    assert len(contents) == 1
    assert contents[0].suffix == ".partd"
    assert contents[0].parent == tmp_path
    f4 = pickle.loads(pickle.dumps(f3))
    assert not f4.buffer
    assert f4.tempdir == tmp_path


def test_set_index_with_explicit_divisions():
    df = pd.DataFrame({"x": [4, 1, 2, 5]}, index=[10, 20, 30, 40])

    ddf = dd.from_pandas(df, npartitions=2)

    def throw(*args, **kwargs):
        raise Exception()

    with dask.config.set(scheduler=throw):
        ddf2 = ddf.set_index("x", divisions=[1, 3, 5])
    assert ddf2.divisions == (1, 3, 5)

    df2 = df.set_index("x")
    assert_eq(ddf2, df2)

    # Divisions must be sorted
    with pytest.raises(ValueError):
        ddf.set_index("x", divisions=[3, 1, 5])


def test_set_index_with_empty_divisions():
    df = pd.DataFrame({"x": [1, 2, 3, 4]})

    ddf = dd.from_pandas(df, npartitions=2)

    # Divisions must not be empty
    with pytest.raises(ValueError):
        ddf.set_index("x", divisions=[])


def test_set_index_divisions_2():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6], "y": list("abdabd")})
    ddf = dd.from_pandas(df, 2)

    result = ddf.set_index("y", divisions=["a", "c", "d"])
    assert result.divisions == ("a", "c", "d")

    assert list(result.compute(scheduler="sync").index[-2:]) == ["d", "d"]


@pytest.mark.slow
def test_set_index_consistent_divisions():
    # See https://github.com/dask/dask/issues/3867
    df = pd.DataFrame(
        {"x": np.random.random(100), "y": np.random.random(100) // 0.2},
        index=np.random.random(100),
    )
    ddf = dd.from_pandas(df, npartitions=4)
    ddf = ddf.clear_divisions()

    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(8, ctx) as pool:
        func = partial(_set_index, df=ddf, idx="x")
        divisions_set = set(pool.map(func, range(100)))
    assert len(divisions_set) == 1


def _set_index(i, df, idx):
    return df.set_index(idx).divisions


def make_part(n):
    return pd.DataFrame({"x": np.random.random(n), "y": np.random.random(n)})


def test_set_index_detects_sorted_data(shuffle_method):
    df = pd.DataFrame({"x": range(100), "y": range(100)})
    ddf = dd.from_pandas(df, npartitions=10, sort=False)
    ddf2 = ddf.set_index("x", shuffle_method=shuffle_method)
    assert len(ddf2.dask) < ddf.npartitions * 4


def test_set_index_sorts():
    # https://github.com/dask/dask/issues/2288
    vals = np.array(
        [
            1348550149000000000,
            1348550149000000000,
            1348558142000000000,
            1348558142000000000,
            1348585928000000000,
            1348585928000000000,
            1348600739000000000,
            1348601706000000000,
            1348600739000000000,
            1348601706000000000,
            1348614789000000000,
            1348614789000000000,
            1348621037000000000,
            1348621038000000000,
            1348621040000000000,
            1348621037000000000,
            1348621038000000000,
            1348621040000000000,
            1348637628000000000,
            1348638159000000000,
            1348638160000000000,
            1348638159000000000,
            1348638160000000000,
            1348637628000000000,
            1348646354000000000,
            1348646354000000000,
            1348659107000000000,
            1348657111000000000,
            1348659107000000000,
            1348657111000000000,
            1348672876000000000,
            1348672876000000000,
            1348682787000000000,
            1348681985000000000,
            1348682787000000000,
            1348681985000000000,
            1348728167000000000,
            1348728167000000000,
            1348730745000000000,
            1348730745000000000,
            1348750198000000000,
            1348750198000000000,
            1348750198000000000,
            1348753539000000000,
            1348753539000000000,
            1348753539000000000,
            1348754449000000000,
            1348754449000000000,
            1348761333000000000,
            1348761554000000000,
            1348761610000000000,
            1348761333000000000,
            1348761554000000000,
            1348761610000000000,
            1348782624000000000,
            1348782624000000000,
            1348782624000000000,
            1348782624000000000,
        ]
    )
    vals = pd.to_datetime(vals, unit="ns")
    breaks = [10, 36, 58]
    dfs = []

    for i in range(len(breaks)):
        lo = sum(breaks[:i])
        hi = sum(breaks[i : i + 1])

        dfs.append(pd.DataFrame({"timestamp": vals[lo:hi]}, index=range(lo, hi)))

    ddf = dd.concat(dfs).clear_divisions()
    assert ddf.set_index("timestamp").index.compute().is_monotonic_increasing is True


@pytest.mark.parametrize(
    "engine", ["pandas", pytest.param("cudf", marks=pytest.mark.gpu)]
)
def test_set_index(engine):
    if engine == "cudf":
        # NOTE: engine == "cudf" requires cudf/dask_cudf,
        # will be skipped by non-GPU CI.

        pytest.importorskip("dask_cudf")

    dsk = {
        ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [4, 2, 6]}, index=[0, 1, 3]),
        ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [3, 5, 8]}, index=[5, 6, 8]),
        ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [9, 1, 8]}, index=[9, 9, 9]),
    }

    d = dd.repartition(pd.concat(dsk.values()), [0, 4, 9, 9])
    full = d.compute()

    d2 = d.set_index("b", npartitions=3)
    assert d2.npartitions == 3
    assert d2.index.name == "b"
    assert_eq(d2, full.set_index("b"))

    d3 = d.set_index(d.b, npartitions=3)
    assert d3.npartitions == 3
    assert d3.index.name == "b"
    assert_eq(d3, full.set_index(full.b))

    d4 = d.set_index("b")
    assert d4.index.name == "b"
    assert_eq(d4, full.set_index("b"))

    d5 = d.set_index(["b"])
    assert d5.index.name == "b"
    assert_eq(d5, full.set_index(["b"]))


@pytest.mark.parametrize(
    "engine", ["pandas", pytest.param("cudf", marks=pytest.mark.gpu)]
)
def test_set_index_interpolate(engine):
    if engine == "cudf":
        # NOTE: engine == "cudf" requires cudf/dask_cudf,
        # will be skipped by non-GPU CI.

        cudf = pytest.importorskip("cudf")
        dask_cudf = pytest.importorskip("dask_cudf")

    df = pd.DataFrame({"x": [4, 1, 1, 3, 3], "y": [1.0, 1, 1, 1, 2]})

    if engine == "cudf":
        gdf = cudf.from_pandas(df)
        d = dask_cudf.from_cudf(gdf, npartitions=3)
    else:
        d = dd.from_pandas(df, 2)

    d1 = d.set_index("x", npartitions=3)
    assert d1.npartitions <= 3
    assert set(d1.divisions) == {1, 2, 4}

    d2 = d.set_index("y", npartitions=3)
    assert d2.divisions[0] == 1.0
    assert 1.0 < d2.divisions[1] < d2.divisions[2] < 2.0
    assert d2.divisions[3] == 2.0


@pytest.mark.parametrize(
    "engine", ["pandas", pytest.param("cudf", marks=pytest.mark.gpu)]
)
def test_set_index_interpolate_int(engine):
    if engine == "cudf":
        # NOTE: engine == "cudf" requires cudf/dask_cudf,
        # will be skipped by non-GPU CI.

        cudf = pytest.importorskip("cudf")
        dask_cudf = pytest.importorskip("dask_cudf")

    L = sorted(list(range(0, 200, 10)) * 2)
    df = pd.DataFrame({"x": 2 * L})

    if engine == "cudf":
        gdf = cudf.from_pandas(df)
        d = dask_cudf.from_cudf(gdf, npartitions=2)
    else:
        d = dd.from_pandas(df, 2)

    d1 = d.set_index("x", npartitions=10)
    assert all(np.issubdtype(type(x), np.integer) for x in d1.divisions)


def test_set_index_npartitions():
    # https://github.com/dask/dask/issues/6974
    data = pd.DataFrame(
        index=pd.Index(
            ["A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "C"]
        )
    )
    data = dd.from_pandas(data, npartitions=2)
    output = data.reset_index().set_index("index", npartitions=1)
    assert output.npartitions == 1


@pytest.mark.parametrize("unit", ["ns", "us"])
def test_set_index_datetime_precision(unit):
    # https://github.com/dask/dask/issues/6864

    df = pd.DataFrame(
        [
            [1567703791155681, 1],
            [1567703792155681, 2],
            [1567703790155681, 0],
            [1567703793155681, 3],
        ],
        columns=["ts", "rank"],
    )
    df.ts = pd.to_datetime(df.ts, unit=unit)
    ddf = dd.from_pandas(df, npartitions=2)
    ddf = ddf.set_index("ts")

    assert_eq(ddf, df.set_index("ts"))


@pytest.mark.parametrize("drop", [True, False])
def test_set_index_drop(drop):
    pdf = pd.DataFrame(
        {
            "A": list("ABAABBABAA"),
            "B": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "C": [1, 2, 3, 2, 1, 3, 2, 4, 2, 3],
        }
    )
    ddf = dd.from_pandas(pdf, 3)

    assert_eq(ddf.set_index("A", drop=drop), pdf.set_index("A", drop=drop))
    assert_eq(ddf.set_index("B", drop=drop), pdf.set_index("B", drop=drop))
    assert_eq(ddf.set_index("C", drop=drop), pdf.set_index("C", drop=drop))
    assert_eq(ddf.set_index(ddf.A, drop=drop), pdf.set_index(pdf.A, drop=drop))
    assert_eq(ddf.set_index(ddf.B, drop=drop), pdf.set_index(pdf.B, drop=drop))
    assert_eq(ddf.set_index(ddf.C, drop=drop), pdf.set_index(pdf.C, drop=drop))

    # numeric columns
    pdf = pd.DataFrame(
        {
            0: list("ABAABBABAA"),
            1: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            2: [1, 2, 3, 2, 1, 3, 2, 4, 2, 3],
        }
    )
    ddf = dd.from_pandas(pdf, 3)
    assert_eq(ddf.set_index(0, drop=drop), pdf.set_index(0, drop=drop))
    assert_eq(ddf.set_index(2, drop=drop), pdf.set_index(2, drop=drop))


def test_set_index_raises_error_on_bad_input():
    df = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7], "b": [7, 6, 5, 4, 3, 2, 1]})
    ddf = dd.from_pandas(df, 2)

    msg = r"Dask dataframe does not yet support multi-indexes"
    with pytest.raises(NotImplementedError) as err:
        ddf.set_index(["a", "b"])
    assert msg in str(err.value)

    with pytest.raises(NotImplementedError) as err:
        ddf.set_index([["a", "b"]])
    assert msg in str(err.value)

    with pytest.raises(NotImplementedError) as err:
        ddf.set_index([["a"]])
    assert msg in str(err.value)


def test_set_index_sorted_true():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 20, 40], "z": [4, 3, 2, 1]})
    a = dd.from_pandas(df, 2, sort=False).clear_divisions()
    assert not a.known_divisions

    b = a.set_index("x", sorted=True)
    assert b.known_divisions
    assert set(a.dask).issubset(set(b.dask))

    for drop in [True, False]:
        assert_eq(a.set_index("x", drop=drop), df.set_index("x", drop=drop))
        assert_eq(
            a.set_index(a.x, sorted=True, drop=drop), df.set_index(df.x, drop=drop)
        )
        assert_eq(
            a.set_index(a.x + 1, sorted=True, drop=drop),
            df.set_index(df.x + 1, drop=drop),
        )


def test_set_index_sorted_single_partition():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [1, 0, 1, 0]})
    ddf = dd.from_pandas(df, npartitions=1)
    assert_eq(ddf.set_index("x", sorted=True), df.set_index("x"))


def test_set_index_sorted_min_max_same():
    a = pd.DataFrame({"x": [1, 2, 3], "y": [0, 0, 0]})
    b = pd.DataFrame({"x": [1, 2, 3], "y": [1, 1, 1]})

    aa = dask.delayed(a)
    bb = dask.delayed(b)

    df = dd.from_delayed([aa, bb], meta=a)
    assert not df.known_divisions

    df2 = df.set_index("y", sorted=True)
    assert df2.divisions == (0, 1, 1)


def test_set_index_empty_partition():
    test_vals = [1, 2, 3]
    converters = [int, float, str, lambda x: pd.to_datetime(x, unit="ns")]

    for conv in converters:
        df = pd.DataFrame(
            [{"x": conv(i), "y": i} for i in test_vals], columns=["x", "y"]
        )
        ddf = dd.concat(
            [
                dd.from_pandas(df, npartitions=1),
                dd.from_pandas(df[df.y > df.y.max()], npartitions=1),
            ]
        )

        assert any(ddf.get_partition(p).compute().empty for p in range(ddf.npartitions))
        assert assert_eq(ddf.set_index("x"), df.set_index("x"))


@pytest.mark.parametrize(
    "converter", [int, float, str, lambda x: pd.to_datetime(x, unit="ns")]
)
def test_set_index_on_empty(converter):
    test_vals = [1, 2, 3, 4]

    df = pd.DataFrame([{"x": converter(x), "y": x} for x in test_vals])
    ddf = dd.from_pandas(df, npartitions=4)

    assert ddf.npartitions > 1

    actual = ddf[ddf.y > df.y.max()].set_index("x")
    expected = df[df.y > df.y.max()].set_index("x")

    assert assert_eq(actual, expected, check_freq=False)
    assert actual.npartitions == 1
    assert all(pd.isnull(d) for d in actual.divisions)


def test_set_index_categorical():
    # https://github.com/dask/dask/issues/5671
    order = list(reversed(string.ascii_letters))
    values = list(string.ascii_letters)
    random.shuffle(values)
    dtype = pd.api.types.CategoricalDtype(order, ordered=True)
    df = pd.DataFrame({"A": pd.Categorical(values, dtype=dtype), "B": 1})

    result = dd.from_pandas(df, npartitions=2).set_index("A")
    assert len(result) == len(df)

    # sorted with the metric defined by the Categorical
    divisions = pd.Categorical(result.divisions, dtype=dtype)
    assert_categorical_equal(divisions, divisions.sort_values())


def test_set_index_with_empty_and_overlap():
    # https://github.com/dask/dask/issues/8735
    df = pd.DataFrame(
        index=list(range(8)),
        data={
            "a": [1, 2, 2, 3, 3, 3, 4, 5],
            "b": [1, 1, 0, 0, 0, 1, 0, 0],
        },
    )
    ddf = dd.from_pandas(df, 4)
    result = ddf[ddf.b == 1].set_index("a", sorted=True)
    expected = df[df.b == 1].set_index("a")

    assert result.divisions == (1.0, 3.0, 3.0)
    assert_eq(result, expected)


def test_compute_divisions():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4], "y": [10, 20, 20, 40], "z": [4, 3, 2, 1]},
        index=[1, 3, 10, 20],
    )
    a = dd.from_pandas(df, 2, sort=False).clear_divisions()
    assert not a.known_divisions

    b = a.compute_current_divisions(set_divisions=True)
    assert_eq(a, b, check_divisions=False)
    assert b.known_divisions


def test_empty_partitions():
    # See https://github.com/dask/dask/issues/2408
    df = pd.DataFrame({"a": list(range(10))})
    df["b"] = df["a"] % 3
    df["c"] = df["b"].astype(str)

    ddf = dd.from_pandas(df, npartitions=3)
    ddf = ddf.set_index("b")
    ddf = ddf.repartition(npartitions=3)
    ddf.get_partition(0).compute()
    assert_eq(ddf, df.set_index("b"))

    ddf = ddf.set_index("c")
    assert_eq(ddf, df.set_index("b").set_index("c"))


@pytest.mark.slow
def test_gh_2730():
    large = pd.DataFrame({"KEY": np.arange(0, 50000)})
    small = pd.DataFrame({"KEY": np.arange(25, 500)})

    dd_left = dd.from_pandas(small, npartitions=3)
    dd_right = dd.from_pandas(large, npartitions=257)

    with dask.config.set({"dataframe.shuffle.method": "tasks", "scheduler": "sync"}):
        dd_merged = dd_left.merge(dd_right, how="inner", on="KEY")
        result = dd_merged.compute()

    expected = large.merge(small, how="inner", on="KEY")

    tm.assert_frame_equal(result.sort_values("KEY").reset_index(drop=True), expected)


def test_set_index_timestamp():
    df = pd.DataFrame({"A": pd.date_range("2000", periods=12, tz="US/Central"), "B": 1})
    ddf = dd.from_pandas(df, 2)
    divisions = (
        pd.Timestamp("2000-01-01 00:00:00-0600", tz="US/Central"),
        pd.Timestamp("2000-01-12 00:00:00-0600", tz="US/Central"),
    )

    # Note: `freq` is lost during round trip
    df2 = df.set_index("A")
    ddf_new_div = ddf.set_index("A", divisions=divisions)
    for ts1, ts2 in zip(divisions, ddf_new_div.divisions):
        assert ts1.timetuple() == ts2.timetuple()
        assert ts1.tz == ts2.tz

    assert_eq(df2, ddf_new_div, check_freq=False)
    assert_eq(df2, ddf.set_index("A"), check_freq=False)


def test_set_index_ea_dtype():
    pdf = pd.DataFrame({"a": 1, "b": pd.Series([1, 2], dtype="Int64")})
    ddf = dd.from_pandas(pdf, npartitions=2)
    pdf_result = pdf.set_index("b")
    ddf_result = ddf.set_index("b")
    assert_eq(ddf_result, pdf_result)


@pytest.mark.parametrize("compression", [None, "ZLib"])
def test_disk_shuffle_with_compression_option(compression):
    # test if dataframe shuffle works both with and without compression
    with dask.config.set({"dataframe.shuffle.compression": compression}):
        test_shuffle("disk")


def test_disk_shuffle_with_unknown_compression():
    # test if dask raises an error in case of fault config string
    with dask.config.set({"dataframe.shuffle.compression": "UNKOWN_COMPRESSION_ALGO"}):
        with pytest.raises(
            ImportError,
            match=(
                "Not able to import and load {} as compression algorithm."
                "Please check if the library is installed and supported by Partd.".format(
                    "UNKOWN_COMPRESSION_ALGO"
                )
            ),
        ):
            shuffle_func(d, d.b, shuffle_method="disk").compute()


def test_disk_shuffle_check_actual_compression():
    # test if the compression switch is really respected by testing the size of the actual partd-data on disk
    def generate_raw_partd_file(compression):
        # generate and write a dummy dataframe to disk and return the raw data bytes
        df1 = pd.DataFrame({"a": list(range(10000))})
        df1["b"] = (df1["a"] * 123).astype(str)
        with dask.config.set({"dataframe.shuffle.compression": compression}):
            p1 = maybe_buffered_partd(buffer=False, tempdir=None)()
            p1.append({"x": df1})
            # get underlying filename from partd - depending on nested structure of partd object
            filename = (
                p1.partd.partd.filename("x") if compression else p1.partd.filename("x")
            )
            with open(filename, "rb") as f:
                return f.read()

    # get compressed and uncompressed raw data
    uncompressed_data = generate_raw_partd_file(compression=None)
    compressed_data = generate_raw_partd_file(compression="BZ2")

    assert len(uncompressed_data) > len(compressed_data)


@pytest.mark.parametrize("ignore_index", [None, True, False])
@pytest.mark.parametrize(
    "on", ["id", "name", ["id", "name"], pd.Series(["id", "name"])]
)
@pytest.mark.parametrize("max_branch", [None, 4])
def test_dataframe_shuffle_on_arg(on, ignore_index, max_branch, shuffle_method):
    # Make sure DataFrame.shuffle API returns the same result
    # whether the ``on`` argument is a list of column names,
    # or a separate DataFrame with equivalent values...
    df_in = dask.datasets.timeseries(
        "2000",
        "2001",
        types={"value": float, "name": str, "id": int},
        freq="2h",
        partition_freq=f"1{ME}",
        seed=1,
    )
    if isinstance(on, str):
        ext_on = df_in[[on]].copy()
    else:
        ext_on = df_in[on].copy()
    df_out_1 = df_in.shuffle(
        on,
        shuffle_method=shuffle_method,
        ignore_index=ignore_index,
        max_branch=max_branch,
    )
    df_out_2 = df_in.shuffle(
        ext_on, shuffle_method=shuffle_method, ignore_index=ignore_index
    )

    assert_eq(df_out_1, df_out_2, check_index=(not ignore_index))

    # disk shuffling doesn't support ignore_index
    if ignore_index and shuffle_method == "tasks":
        assert df_out_1.index.dtype != df_in.index.dtype
    else:
        assert df_out_1.index.dtype == df_in.index.dtype


def test_set_index_overlap():
    A = pd.DataFrame({"key": [1, 2, 3, 4, 4, 5, 6, 7], "value": list("abcd" * 2)})
    a = dd.from_pandas(A, npartitions=2)
    a = a.set_index("key", sorted=True)
    b = a.repartition(divisions=a.divisions)
    assert_eq(a, b)


def test_set_index_overlap_2():
    df = pd.DataFrame(
        index=pd.Index(
            ["A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "C"],
            name="index",
        )
    )
    ddf = dd.from_pandas(df, npartitions=2)
    result = (
        ddf.reset_index().repartition(npartitions=8).set_index("index", sorted=True)
    )
    expected = df.reset_index().set_index("index")
    assert_eq(result, expected)
    assert result.npartitions == 8


def test_set_index_overlap_does_not_drop_rows_when_divisions_overlap():
    # https://github.com/dask/dask/issues/9339
    df = pd.DataFrame({"ts": [1, 1, 2, 2, 3, 3, 3, 3], "value": "abc"})
    ddf = dd.from_pandas(df, npartitions=3)

    expected = df.set_index("ts")
    actual = ddf.set_index("ts", sorted=True)

    assert_eq(expected, actual)


def test_compute_current_divisions_nan_partition():
    # Compute divisions 1 null partition
    a = d[d.a > 3].sort_values("a")
    divisions = a.compute_current_divisions("a")
    assert divisions == (4, 5, 8, 9)
    assert_eq(a, a, check_divisions=False)

    # Compute divisions with 0 null partitions
    a = d[d.a > 1].sort_values("a")
    divisions = a.compute_current_divisions("a")
    assert divisions == (2, 4, 7, 9)
    assert_eq(a, a, check_divisions=False)


def test_compute_current_divisions_overlap():
    A = pd.DataFrame({"key": [1, 2, 3, 4, 4, 5, 6, 7], "value": list("abcd" * 2)})
    a = dd.from_pandas(A, npartitions=2)
    with pytest.warns(UserWarning, match="Partitions have overlapping values"):
        divisions = a.compute_current_divisions("key")
        b = a.set_index("key", divisions=divisions)
        assert b.divisions == (1, 4, 7)
        assert [len(p) for p in b.partitions] == [3, 5]


def test_compute_current_divisions_overlap_2():
    data = pd.DataFrame(
        index=pd.Index(
            ["A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "C"],
            name="index",
        )
    )
    ddf1 = dd.from_pandas(data, npartitions=2)
    ddf2 = ddf1.clear_divisions().repartition(npartitions=8)
    with pytest.warns(UserWarning, match="Partitions have overlapping values"):
        ddf2.compute_current_divisions()


def test_set_index_nan_partition():
    d[d.a > 3].set_index("a")  # Set index with 1 null partition
    d[d.a > 1].set_index("a", sorted=True)  # Set sorted index with 0 null partitions
    a = d[d.a > 3].set_index("a", sorted=True)  # Set sorted index with 1 null partition
    assert_eq(a, a)


def test_set_index_with_dask_dt_index():
    values = {
        "x": [1, 2, 3, 4] * 3,
        "y": [10, 20, 30] * 4,
        "name": ["Alice", "Bob"] * 6,
    }
    date_index = pd.date_range(
        start="2022-02-22", freq="16h", periods=12
    ) - pd.Timedelta(seconds=30)
    df = pd.DataFrame(values, index=date_index)
    ddf = dd.from_pandas(df, npartitions=3)

    # specify a different date index entirely
    day_index = ddf.index.dt.floor("D")
    day_df = ddf.set_index(day_index)
    expected = dd.from_pandas(
        pd.DataFrame(values, index=date_index.floor("D")), npartitions=3
    )
    assert_eq(day_df, expected)

    # specify an index with shifted dates
    one_day = pd.Timedelta(days=1)
    next_day_df = ddf.set_index(ddf.index + one_day)
    expected = dd.from_pandas(
        pd.DataFrame(values, index=date_index + one_day), npartitions=3
    )
    assert_eq(next_day_df, expected)

    # try a different index type
    no_dates = dd.from_pandas(pd.DataFrame(values), npartitions=3)
    range_df = ddf.set_index(no_dates.index)
    expected = dd.from_pandas(pd.DataFrame(values), npartitions=3)
    assert_eq(range_df, expected)


def test_set_index_with_series_uses_fastpath():
    dates = pd.date_range(start="2022-02-22", freq="16h", periods=12) - pd.Timedelta(
        seconds=30
    )
    one_day = pd.Timedelta(days=1)
    df = pd.DataFrame(
        {
            "x": [1, 2, 3, 4] * 3,
            "y": [10, 20, 30] * 4,
            "name": ["Alice", "Bob"] * 6,
            "d1": dates + one_day,
            "d2": dates + one_day * 5,
        },
        index=dates,
    )
    ddf = dd.from_pandas(df, npartitions=3)

    res = ddf.set_index(ddf.d2 + one_day)
    expected = df.set_index(df.d2 + one_day)
    assert_eq(res, expected)


@pytest.mark.parametrize("ascending", [True, False])
@pytest.mark.parametrize("by", ["a", "b", ["a", "b"]])
@pytest.mark.parametrize("nelem", [10, 500])
def test_sort_values(nelem, by, ascending):
    np.random.seed(0)
    df = pd.DataFrame()
    df["a"] = np.ascontiguousarray(np.arange(nelem)[::-1])
    df["b"] = np.arange(100, nelem + 100)
    ddf = dd.from_pandas(df, npartitions=10)

    # run on single-threaded scheduler for debugging purposes
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(by=by, ascending=ascending)
    expect = df.sort_values(by=by, ascending=ascending)
    dd.assert_eq(got, expect, check_index=False, sort_results=False)


@pytest.mark.parametrize(
    "backend", ["pandas", pytest.param("cudf", marks=pytest.mark.gpu)]
)
@pytest.mark.parametrize("by", ["x", "z", ["x", "z"], ["z", "x"]])
@pytest.mark.parametrize("ascending", [True, False])
def test_sort_values_tasks_backend(backend, by, ascending):
    if backend == "cudf":
        pytest.importorskip("dask_cudf")
    pdf = pd.DataFrame(
        {"x": range(10), "y": [1, 2, 3, 4, 5] * 2, "z": ["cat", "dog"] * 5}
    )
    ddf = dd.from_pandas(pdf, npartitions=10)
    if backend == "cudf":
        ddf = ddf.to_backend(backend)

    expect = pdf.sort_values(by=by, ascending=ascending)
    got = dd.DataFrame.sort_values(
        ddf, by=by, ascending=ascending, shuffle_method="tasks"
    )
    dd.assert_eq(got, expect, sort_results=False)


@pytest.mark.parametrize("ascending", [True, False, [False, True], [True, False]])
@pytest.mark.parametrize("by", [["a", "b"], ["b", "a"]])
@pytest.mark.parametrize("nelem", [10, 500])
def test_sort_values_single_partition(nelem, by, ascending):
    np.random.seed(0)
    df = pd.DataFrame()
    df["a"] = np.ascontiguousarray(np.arange(nelem)[::-1])
    df["b"] = np.arange(100, nelem + 100)
    ddf = dd.from_pandas(df, npartitions=1)

    # run on single-threaded scheduler for debugging purposes
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(by=by, ascending=ascending)
    expect = df.sort_values(by=by, ascending=ascending)
    dd.assert_eq(got, expect, check_index=False)


@pytest.mark.parametrize("na_position", ["first", "last"])
@pytest.mark.parametrize("ascending", [True, False])
@pytest.mark.parametrize("by", ["a", "b", ["a", "b"]])
@pytest.mark.parametrize("nparts", [1, 5])
@pytest.mark.parametrize(
    "data",
    [
        {
            "a": list(range(50)) + [None] * 50 + list(range(50, 100)),
            "b": [None] * 100 + list(range(100, 150)),
        },
        {
            "a": list(range(15)) + [None] * 5,
            "b": list(reversed(range(20))),
        },
    ],
)
def test_sort_values_with_nulls(data, nparts, by, ascending, na_position):
    df = pd.DataFrame(data)
    ddf = dd.from_pandas(df, npartitions=nparts)

    # run on single-threaded scheduler for debugging purposes
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(by=by, ascending=ascending, na_position=na_position)
    expect = df.sort_values(by=by, ascending=ascending, na_position=na_position)
    dd.assert_eq(got, expect, check_index=False)


def test_shuffle_values_raises():
    df = pd.DataFrame({"a": [1, 3, 2]})
    ddf = dd.from_pandas(df, npartitions=3)
    with pytest.raises(
        ValueError, match="na_position must be either 'first' or 'last'"
    ):
        ddf.sort_values(by="a", na_position="invalid")


def test_shuffle_by_as_list():
    df = pd.DataFrame({"a": [1, 3, 2]})
    ddf = dd.from_pandas(df, npartitions=3)
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(by=["a"], ascending=True)
    expect = pd.DataFrame({"a": [1, 2, 3]})
    dd.assert_eq(got, expect, check_index=False)


@pytest.mark.parametrize("by", [["a", "b"], ["b", "a"]])
@pytest.mark.parametrize("nparts", [1, 10])
def test_sort_values_custom_function(by, nparts):
    df = pd.DataFrame({"a": [1, 2, 3] * 20, "b": [4, 5, 6, 7] * 15})
    ddf = dd.from_pandas(df, npartitions=nparts)

    def f(partition, by_columns, ascending, na_position, **kwargs):
        return partition.sort_values(
            by_columns, ascending=ascending, na_position=na_position
        )

    # run on single-threaded scheduler for debugging purposes
    with dask.config.set(scheduler="single-threaded"):
        got = ddf.sort_values(
            by=by[0], sort_function=f, sort_function_kwargs={"by_columns": by}
        )
    expect = df.sort_values(by=by)
    dd.assert_eq(got, expect, check_index=False)


def test_sort_values_bool_ascending():
    df = pd.DataFrame({"a": [1, 2, 3] * 20, "b": [4, 5, 6, 7] * 15})
    ddf = dd.from_pandas(df, npartitions=10)

    # attempt to sort with list of ascending booleans
    with pytest.raises(ValueError, match="length"):
        ddf.sort_values(by="a", ascending=[True, False])
    with pytest.raises(ValueError, match="length"):
        ddf.sort_values(by=["a", "b"], ascending=[True])
    assert_eq(
        ddf.sort_values(by="a", ascending=[True]),
        df.sort_values(by="a", ascending=[True]),
    )
    assert_eq(
        ddf.sort_values(by=["a", "b"], ascending=[True, False]),
        df.sort_values(by=["a", "b"], ascending=[True, False]),
    )


@pytest.mark.parametrize("npartitions", [1, 3])
def test_sort_values_timestamp(npartitions):
    # Regression test for https://github.com/dask/dask/issues/9641
    df = pd.DataFrame.from_records(
        [
            [pd.Timestamp("2002-01-11 21:00:01+0000", tz="UTC"), 4223, 54719.0],
            [pd.Timestamp("2002-01-14 21:00:01+0000", tz="UTC"), 6942, 19223.0],
            [pd.Timestamp("2002-01-15 21:00:01+0000", tz="UTC"), 12551, 72865.0],
            [pd.Timestamp("2002-01-23 21:00:01+0000", tz="UTC"), 6005, 57670.0],
            [pd.Timestamp("2002-01-29 21:00:01+0000", tz="UTC"), 2043, 58600.0],
            [pd.Timestamp("2002-02-01 21:00:01+0000", tz="UTC"), 6909, 8459.0],
            [pd.Timestamp("2002-01-14 21:00:01+0000", tz="UTC"), 5326, 77339.0],
            [pd.Timestamp("2002-01-14 21:00:01+0000", tz="UTC"), 4711, 54135.0],
            [pd.Timestamp("2002-01-22 21:00:01+0000", tz="UTC"), 103, 57627.0],
            [pd.Timestamp("2002-01-30 21:00:01+0000", tz="UTC"), 16862, 54458.0],
            [pd.Timestamp("2002-01-31 21:00:01+0000", tz="UTC"), 4143, 56280.0],
        ],
        columns=["time", "id1", "id2"],
    )
    ddf = dd.from_pandas(df, npartitions=npartitions)
    result = ddf.sort_values("time")
    expected = df.sort_values("time")
    assert_eq(result, expected)


@pytest.mark.skipif(pa is None, reason="Need pyarrow")
@pytest.mark.parametrize(
    "data, dtype",
    [
        (["a", "b"], "string[pyarrow]"),
        ([b"a", b"b"], "binary[pyarrow]"),
        # Should probably fix upstream, https://github.com/pandas-dev/pandas/issues/52590
        # (["a", "b"], pa.large_string()),
        # ([b"a", b"b"], pa.large_binary()),
        ([1, 2], "int64[pyarrow]"),
        ([1, 2], "float64[pyarrow]"),
        ([1, 2], "uint64[pyarrow]"),
        ([date(2022, 1, 1), date(1999, 12, 31)], "date32[pyarrow]"),
        (
            [pd.Timestamp("2022-01-01"), pd.Timestamp("2023-01-02")],
            "timestamp[ns][pyarrow]",
        ),
        ([Decimal("5"), Decimal("6.24")], "decimal128"),
        ([pd.Timedelta("1 day"), pd.Timedelta("20 days")], "duration[ns][pyarrow]"),
        ([time(12, 0), time(0, 12)], "time64[ns][pyarrow]"),
    ],
)
def test_set_index_pyarrow_dtype(data, dtype):
    if dtype == "decimal128":
        dtype = pd.ArrowDtype(pa.decimal128(10, 2))
    pdf = pd.DataFrame({"a": 1, "arrow_col": pd.Series(data, dtype=dtype)})
    ddf = dd.from_pandas(pdf, npartitions=2)
    pdf_result = pdf.set_index("arrow_col")
    ddf_result = ddf.set_index("arrow_col")
    assert_eq(ddf_result, pdf_result)


def test_shuffle_nulls_introduced():
    df1 = pd.DataFrame([[True], [False]] * 50, columns=["A"])
    df1["B"] = list(range(100))

    df2 = pd.DataFrame(
        [[2, 3], [109, 2], [345, 3], [50, 7], [95, 1]], columns=["B", "C"]
    )

    ddf1 = dd.from_pandas(df1, npartitions=10)
    ddf2 = dd.from_pandas(df2, npartitions=1)
    meta = pd.Series(dtype=int, index=pd.Index([], dtype=bool, name="A"), name="A")
    include_groups = {"include_groups": False} if PANDAS_GE_220 else {}
    result = (
        dd.merge(ddf1, ddf2, how="outer", on="B")
        .groupby("A")
        .apply(lambda df: len(df), meta=meta, **include_groups)
    )
    expected = (
        pd.merge(df1, df2, how="outer", on="B")
        .groupby("A")
        .apply(lambda df: len(df), **include_groups)
    )
    assert_eq(result, expected, check_names=False)
