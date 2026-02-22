from __future__ import annotations

from datetime import datetime

import numpy as np
import pandas as pd
import pytest

import dask
import dask.array as da
import dask.dataframe as dd
from dask import config
from dask.dataframe._compat import tm
from dask.dataframe.io.io import _meta_from_array, sorted_division_locations
from dask.dataframe.utils import assert_eq, get_string_dtype
from dask.delayed import Delayed, delayed

##########
# Arrays #
##########


def test_meta_from_array():
    x = np.array([[1, 2], [3, 4]], dtype=np.int64)
    res = _meta_from_array(x)
    assert isinstance(res, pd.DataFrame)
    assert res[0].dtype == np.int64
    assert res[1].dtype == np.int64
    tm.assert_index_equal(res.columns, pd.Index([0, 1]))

    x = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    res = _meta_from_array(x, columns=["a", "b"])
    assert isinstance(res, pd.DataFrame)
    assert res["a"].dtype == np.float64
    assert res["b"].dtype == np.float64
    tm.assert_index_equal(res.columns, pd.Index(["a", "b"]))

    with pytest.raises(ValueError):
        _meta_from_array(x, columns=["a", "b", "c"])

    np.random.seed(42)
    x = np.random.rand(201, 2)
    x = dd.from_array(x, chunksize=50, columns=["a", "b"])
    assert len(x.divisions) == 6  # Should be 5 partitions and the end


def test_meta_from_1darray():
    x = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    res = _meta_from_array(x)
    assert isinstance(res, pd.Series)
    assert res.dtype == np.float64

    x = np.array([1, 2, 3], dtype=np.object_)
    res = _meta_from_array(x, columns="x")
    assert isinstance(res, pd.Series)
    assert res.name == "x"
    assert res.dtype == np.object_

    x = np.array([1, 2, 3], dtype=np.object_)
    res = _meta_from_array(x, columns=["x"])
    assert isinstance(res, pd.DataFrame)
    assert res["x"].dtype == np.object_
    tm.assert_index_equal(res.columns, pd.Index(["x"]))

    with pytest.raises(ValueError):
        _meta_from_array(x, columns=["a", "b"])


def test_meta_from_recarray():
    x = np.array(
        [(i, i * 10) for i in range(10)], dtype=[("a", np.float64), ("b", np.int64)]
    )
    res = _meta_from_array(x)
    assert isinstance(res, pd.DataFrame)
    assert res["a"].dtype == np.float64
    assert res["b"].dtype == np.int64
    tm.assert_index_equal(res.columns, pd.Index(["a", "b"]))

    res = _meta_from_array(x, columns=["b", "a"])
    assert isinstance(res, pd.DataFrame)
    assert res["a"].dtype == np.float64
    assert res["b"].dtype == np.int64
    tm.assert_index_equal(res.columns, pd.Index(["b", "a"]))

    with pytest.raises(ValueError):
        _meta_from_array(x, columns=["a", "b", "c"])


def test_from_array():
    x = np.arange(10 * 3).reshape(10, 3)
    d = dd.from_array(x, chunksize=4)
    assert isinstance(d, dd.DataFrame)
    tm.assert_index_equal(d.columns, pd.Index([0, 1, 2]))
    assert d.divisions == (0, 4, 8, 9)
    assert (d.compute().values == x).all()

    d = dd.from_array(x, chunksize=4, columns=list("abc"))
    assert isinstance(d, dd.DataFrame)
    tm.assert_index_equal(d.columns, pd.Index(["a", "b", "c"]))
    assert d.divisions == (0, 4, 8, 9)
    assert (d.compute().values == x).all()

    with pytest.raises(ValueError):
        dd.from_array(np.ones(shape=(10, 10, 10)))


def test_from_array_with_record_dtype():
    x = np.array([(i, i * 10) for i in range(10)], dtype=[("a", "i4"), ("b", "i4")])
    d = dd.from_array(x, chunksize=4)
    assert isinstance(d, dd.DataFrame)
    assert list(d.columns) == ["a", "b"]
    assert d.divisions == (0, 4, 8, 9)

    assert (d.compute().to_records(index=False) == x).all()


def test_from_pandas_dataframe():
    a = list("aaaaaaabbbbbbbbccccccc")
    df = pd.DataFrame(
        dict(a=a, b=np.random.randn(len(a))),
        index=pd.date_range(start="20120101", periods=len(a)),
    )
    ddf = dd.from_pandas(df, 3)
    expected_layers = 3
    assert len(ddf.dask) == expected_layers
    assert len(ddf.divisions) == 4
    assert isinstance(ddf.divisions[0], type(df.index[0]))
    assert_eq(df, ddf)

    ddf = dd.from_pandas(df, chunksize=8)
    msg = "Exactly one of npartitions and chunksize must be specified."
    with pytest.raises(ValueError) as err:
        dd.from_pandas(df, npartitions=2, chunksize=2)
    assert msg in str(err.value)
    assert len(ddf.dask) == expected_layers
    assert len(ddf.divisions) == 4
    assert isinstance(ddf.divisions[0], type(df.index[0]))
    assert_eq(df, ddf)


def test_from_pandas_small():
    df = pd.DataFrame({"x": [1, 2, 3]})
    for i in [1, 2, 30]:
        a = dd.from_pandas(df, i)
        assert len(a.compute()) == 3
        assert a.divisions[0] == 0
        assert a.divisions[-1] == 2

        a = dd.from_pandas(df, chunksize=i)
        assert len(a.compute()) == 3
        assert a.divisions[0] == 0
        assert a.divisions[-1] == 2

    for sort in [True, False]:
        for i in [0, 2]:
            df = pd.DataFrame({"x": [0] * i})
            ddf = dd.from_pandas(df, npartitions=5, sort=sort)
            assert_eq(df, ddf)

            s = pd.Series([0] * i, name="x", dtype=int)
            ds = dd.from_pandas(s, npartitions=5, sort=sort)
            assert_eq(s, ds)


@pytest.mark.parametrize("n", [1, 2, 4, 5])
def test_from_pandas_npartitions_is_accurate(n):
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5, 6], "y": list("abdabd")}, index=[10, 20, 30, 40, 50, 60]
    )
    assert dd.from_pandas(df, npartitions=n).npartitions <= n


def test_from_pandas_series():
    n = 20
    s = pd.Series(np.random.randn(n), index=pd.date_range(start="20120101", periods=n))
    ds = dd.from_pandas(s, 3)
    assert len(ds.dask) == 3
    assert len(ds.divisions) == len(ds.dask) + 1
    assert isinstance(ds.divisions[0], type(s.index[0]))
    tm.assert_series_equal(s, ds.compute())

    ds = dd.from_pandas(s, chunksize=8)
    assert len(ds.dask) == 3
    assert len(ds.divisions) == len(ds.dask) + 1
    assert isinstance(ds.divisions[0], type(s.index[0]))
    tm.assert_series_equal(s, ds.compute())


def test_from_pandas_non_sorted():
    df = pd.DataFrame({"x": [1, 2, 3]}, index=[3, 1, 2])
    ddf = dd.from_pandas(df, npartitions=2, sort=False)
    assert not ddf.known_divisions
    assert_eq(df, ddf)

    ddf = dd.from_pandas(df, chunksize=2, sort=False)
    assert not ddf.known_divisions
    assert_eq(df, ddf)


def test_from_pandas_single_row():
    df = pd.DataFrame({"x": [1]}, index=[1])
    ddf = dd.from_pandas(df, npartitions=1)
    assert ddf.divisions == (1, 1)
    assert_eq(ddf, df)


def test_from_pandas_with_datetime_index():
    df = pd.DataFrame(
        {
            "Date": [
                "2015-08-28",
                "2015-08-27",
                "2015-08-26",
                "2015-08-25",
                "2015-08-24",
                "2015-08-21",
                "2015-08-20",
                "2015-08-19",
                "2015-08-18",
            ],
            "Val": list(range(9)),
        }
    )
    df.Date = df.Date.astype("datetime64[ns]")
    ddf = dd.from_pandas(df, 2)
    assert_eq(df, ddf)
    ddf = dd.from_pandas(df, chunksize=2)
    assert_eq(df, ddf)


@pytest.mark.parametrize("null_value", [None, pd.NaT, pd.NA])
def test_from_pandas_with_index_nulls(null_value):
    df = pd.DataFrame({"x": [1, 2, 3]}, index=["C", null_value, "A"])
    with pytest.raises(NotImplementedError, match="is non-numeric and contains nulls"):
        dd.from_pandas(df, npartitions=2, sort=False)


def test_from_pandas_with_wrong_args():
    df = pd.DataFrame({"x": [1, 2, 3]}, index=[3, 2, 1])
    with pytest.raises(TypeError, match="must be a pandas DataFrame or Series"):
        dd.from_pandas("foo")
    with pytest.raises(TypeError, match="provide npartitions as an int"):
        dd.from_pandas(df, npartitions=5.2, sort=False)
    with pytest.raises(TypeError, match="provide chunksize as an int"):
        dd.from_pandas(df, chunksize=18.27)


def test_from_pandas_chunksize_one():
    # See: https://github.com/dask/dask/issues/9218
    df = pd.DataFrame(np.random.randint(0, 10, size=(10, 4)), columns=list("ABCD"))
    ddf = dd.from_pandas(df, chunksize=1)
    num_rows = list(ddf.map_partitions(len).compute())
    # chunksize=1 with range index should
    # always have unit-length partitions
    assert num_rows == [1] * 10


@pytest.mark.parametrize(
    "index",
    [
        ["A", "B", "C", "C", "C", "C", "C", "C"],
        ["A", "B", "B", "B", "B", "B", "B", "C"],
        ["A", "A", "A", "A", "A", "B", "B", "C"],
    ],
)
def test_from_pandas_npartitions_duplicates(index):
    df = pd.DataFrame({"a": range(8), "index": index}).set_index("index")
    ddf = dd.from_pandas(df, npartitions=3)
    assert ddf.divisions == ("A", "B", "C", "C")


def test_from_pandas_convert_string_config():
    pytest.importorskip("pyarrow", reason="Requires pyarrow strings")
    # With `dataframe.convert-string=False`, strings should remain objects
    with dask.config.set({"dataframe.convert-string": False}):
        s = pd.Series(["foo", "bar", "ricky", "bobby"], index=["a", "b", "c", "d"])
        df = pd.DataFrame(
            {
                "x": [1, 2, 3, 4],
                "y": [5.0, 6.0, 7.0, 8.0],
                "z": ["foo", "bar", "ricky", "bobby"],
            },
            index=["a", "b", "c", "d"],
        )

        ds = dd.from_pandas(s, npartitions=2)
        ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(s, ds)
    assert_eq(df, ddf)

    # When `dataframe.convert-string = True`, dask should automatically
    # cast `object`s to pyarrow strings
    with dask.config.set({"dataframe.convert-string": True}):
        ds = dd.from_pandas(s, npartitions=2)
        ddf = dd.from_pandas(df, npartitions=2)

    s_pyarrow = s.astype("string[pyarrow]")
    s_pyarrow.index = s_pyarrow.index.astype("string[pyarrow]")
    df_pyarrow = df.astype({"z": "string[pyarrow]"})
    df_pyarrow.index = df_pyarrow.index.astype("string[pyarrow]")
    assert_eq(s_pyarrow, ds)
    assert_eq(df_pyarrow, ddf)


@pytest.mark.parametrize("index", [[1, 2, 3], [3, 2, 1]])
@pytest.mark.parametrize("sort", [True, False])
def test_from_pandas_immutable(sort, index):
    pdf = pd.DataFrame({"a": [1, 2, 3]}, index=index)
    expected = pdf.copy()
    df = dd.from_pandas(pdf, npartitions=2, sort=sort)
    pdf.iloc[0, 0] = 100
    assert_eq(df, expected)


@pytest.mark.gpu
def test_gpu_from_pandas_npartitions_duplicates():
    cudf = pytest.importorskip("cudf")

    index = ["A", "A", "A", "A", "A", "B", "B", "C"]
    df = cudf.DataFrame({"a": range(8), "index": index}).set_index("index")
    ddf = dd.from_pandas(df, npartitions=3)
    assert ddf.divisions == ("A", "B", "C", "C")


def test_DataFrame_from_dask_array():
    x = da.ones((10, 3), chunks=(4, 2))
    pdf = pd.DataFrame(np.ones((10, 3)), columns=["a", "b", "c"])
    df = dd.from_dask_array(x, ["a", "b", "c"])
    assert_eq(df, pdf)
    # dd.from_array should re-route to from_dask_array
    df2 = dd.from_array(x, columns=["a", "b", "c"])
    assert_eq(df, df2)


def test_DataFrame_from_dask_array_with_blockwise_ops():
    x = da.ones((10, 3), chunks=(4, 2))
    x *= 2
    pdf = pd.DataFrame(np.ones((10, 3)) * 2, columns=["a", "b", "c"])
    df = dd.from_dask_array(x, ["a", "b", "c"])
    assert_eq(df, pdf)


def test_Series_from_dask_array():
    x = da.ones(10, chunks=4)
    pser = pd.Series(np.ones(10), name="a")

    ser = dd.from_dask_array(x, "a")
    assert_eq(ser, pser)

    # Not passing a name should result in the name == None
    pser = pd.Series(np.ones(10))
    ser = dd.from_dask_array(x)
    assert_eq(ser, pser)

    # dd.from_array should re-route to from_dask_array
    ser2 = dd.from_array(x)
    assert_eq(ser, ser2)


@pytest.mark.parametrize("as_frame", [True, False])
def test_from_dask_array_index(as_frame):
    s = dd.from_pandas(pd.Series(range(10), index=list("abcdefghij")), npartitions=3)
    if as_frame:
        s = s.to_frame()
    result = dd.from_dask_array(s.values, index=s.index)
    assert_eq(s, result)


def test_from_dask_array_index_raises():
    x = da.random.uniform(size=(10,), chunks=(5,))
    with pytest.raises(ValueError, match="must be an instance"):
        dd.from_dask_array(x, index=pd.Index(np.arange(10)))

    a = dd.from_pandas(pd.Series(range(12)), npartitions=2)
    b = dd.from_pandas(pd.Series(range(12)), npartitions=4)
    with pytest.raises(ValueError, match=".*index.*numbers of blocks.*4 != 2"):
        dd.from_dask_array(a.values, index=b.index)


def test_from_array_raises_more_than_2D():
    x = da.ones((3, 3, 3), chunks=2)
    y = np.ones((3, 3, 3))

    with pytest.raises(ValueError, match="more than 2D array"):
        dd.from_dask_array(x)  # dask

    with pytest.raises(ValueError, match="more than 2D array"):
        dd.from_array(y)  # numpy


def test_from_dask_array_compat_numpy_array():
    x = da.ones((10, 3), chunks=(3, 3))
    y = np.ones((10, 3))
    d1 = dd.from_dask_array(x)  # dask
    p1 = pd.DataFrame(y)
    assert_eq(d1, p1)

    d2 = dd.from_array(y)  # numpy
    assert_eq(d2, d1)


def test_from_array_wrong_column_shape_error():
    x = da.ones((10, 3), chunks=(3, 3))
    with pytest.raises(ValueError, match="names must match width"):
        dd.from_dask_array(x, columns=["a"])  # dask

    y = np.ones((10, 3))
    with pytest.raises(ValueError, match="names must match width"):
        dd.from_array(y, columns=["a"])  # numpy


def test_from_array_with_column_names():
    x = da.ones((10, 3), chunks=(3, 3))
    y = np.ones((10, 3))
    d1 = dd.from_dask_array(x, columns=["a", "b", "c"])  # dask
    p1 = pd.DataFrame(y, columns=["a", "b", "c"])
    assert_eq(d1, p1)

    d2 = dd.from_array(y, columns=["a", "b", "c"])  # numpy
    assert_eq(d1, d2)


def test_from_dask_array_compat_numpy_array_1d():
    x = da.ones(10, chunks=3)
    y = np.ones(10)
    d1 = dd.from_dask_array(x)  # dask
    p1 = pd.Series(y)
    assert_eq(d1, p1)

    d2 = dd.from_array(y)  # numpy
    assert_eq(d2, d1)


def test_from_array_1d_with_column_names():
    x = da.ones(10, chunks=3)
    y = np.ones(10)
    d1 = dd.from_dask_array(x, columns="name")  # dask
    p1 = pd.Series(y, name="name")
    assert_eq(d1, p1)

    d2 = dd.from_array(x.compute(), columns="name")  # numpy
    assert_eq(d2, d1)


def test_from_array_1d_list_of_columns_gives_dataframe():
    x = da.ones(10, chunks=3)
    y = np.ones(10)
    # passing list via columns results in DataFrame
    d1 = dd.from_dask_array(x, columns=["name"])  # dask
    p1 = pd.DataFrame(y, columns=["name"])
    assert_eq(d1, p1)

    d2 = dd.from_array(y, columns=["name"])  # numpy
    assert_eq(d2, d1)


def test_from_dask_array_struct_dtype():
    x = np.array([(1, "a"), (2, "b")], dtype=[("a", "i4"), ("b", "object")])
    y = da.from_array(x, chunks=(1,))
    df = dd.from_dask_array(y)
    tm.assert_index_equal(df.columns, pd.Index(["a", "b"]))
    assert_eq(df, pd.DataFrame(x))

    assert_eq(
        dd.from_dask_array(y, columns=["b", "a"]), pd.DataFrame(x, columns=["b", "a"])
    )


def test_from_dask_array_unknown_chunks():
    # Series
    dx = da.Array(
        {("x", 0): np.arange(5), ("x", 1): np.arange(5, 11)},
        "x",
        ((np.nan, np.nan),),
        np.arange(1).dtype,
    )
    df = dd.from_dask_array(dx)
    assert isinstance(df, dd.Series)
    assert not df.known_divisions
    assert_eq(df, pd.Series(np.arange(11)), check_index=False)

    # DataFrame
    dsk = {("x", 0, 0): np.random.random((2, 3)), ("x", 1, 0): np.random.random((5, 3))}
    dx = da.Array(dsk, "x", ((np.nan, np.nan), (3,)), np.float64)
    df = dd.from_dask_array(dx)
    assert isinstance(df, dd.DataFrame)
    assert not df.known_divisions
    assert_eq(df, pd.DataFrame(dx.compute()), check_index=False)


@pytest.mark.parametrize(
    "chunksizes, expected_divisions",
    [
        pytest.param((1, 2, 3, 0), (0, 1, 3, 5, 5)),
        pytest.param((0, 1, 2, 3), (0, 0, 1, 3, 5)),
        pytest.param((1, 0, 2, 3), (0, 1, 1, 3, 5)),
    ],
)
def test_from_dask_array_empty_chunks(chunksizes, expected_divisions):
    monotonic_index = da.from_array(np.arange(6), chunks=chunksizes)
    df = dd.from_dask_array(monotonic_index)
    assert df.divisions == expected_divisions


def test_from_dask_array_unknown_width_error():
    dsk = {("x", 0, 0): np.random.random((2, 3)), ("x", 1, 0): np.random.random((5, 3))}
    dx = da.Array(dsk, "x", ((np.nan, np.nan), (np.nan,)), np.float64)
    with pytest.raises(ValueError, match="Shape along axis 1 must be known"):
        dd.from_dask_array(dx)


@pytest.mark.gpu
@pytest.mark.parametrize(
    "array_backend, df_backend",
    [("cupy", "cudf"), ("numpy", "pandas")],
)
def test_from_array_dispatching(array_backend, df_backend):
    # Check array -> dataframe dispatching
    array_lib = pytest.importorskip(array_backend)
    df_lib = pytest.importorskip(df_backend)

    with config.set({"array.backend": array_backend}):
        darr = da.ones(10)
    assert isinstance(darr._meta, array_lib.ndarray)

    ddf1 = dd.from_array(darr)  # Invokes `from_dask_array`
    ddf2 = dd.from_array(darr.compute())

    assert isinstance(ddf1._meta, df_lib.Series)
    assert isinstance(ddf2._meta, df_lib.Series)
    assert_eq(ddf1, ddf2)


def test_to_bag():
    a = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
    )
    ddf = dd.from_pandas(a, 2)

    assert ddf.to_bag().compute() == list(a.itertuples(False))
    assert ddf.to_bag(True).compute() == list(a.itertuples(True))
    assert ddf.to_bag(format="dict").compute() == [
        {"x": "a", "y": 2},
        {"x": "b", "y": 3},
        {"x": "c", "y": 4},
        {"x": "d", "y": 5},
    ]
    assert ddf.to_bag(True, format="dict").compute() == [
        {"index": 1.0, "x": "a", "y": 2},
        {"index": 2.0, "x": "b", "y": 3},
        {"index": 3.0, "x": "c", "y": 4},
        {"index": 4.0, "x": "d", "y": 5},
    ]
    assert ddf.x.to_bag(True).compute() == list(a.x.items())
    assert ddf.x.to_bag().compute() == list(a.x)

    assert ddf.x.to_bag(True, format="dict").compute() == [
        {"x": "a"},
        {"x": "b"},
        {"x": "c"},
        {"x": "d"},
    ]
    assert ddf.x.to_bag(format="dict").compute() == [
        {"x": "a"},
        {"x": "b"},
        {"x": "c"},
        {"x": "d"},
    ]


def test_to_bag_frame():
    from dask import get
    from dask.bag import Bag

    ddf = dd.from_pandas(
        pd.DataFrame(
            {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
            index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
        ),
        npartitions=2,
    )

    # Convert to bag, and check that
    # collection type has changed, but
    # partition data has not
    bagdf = ddf.to_bag(format="frame")
    assert isinstance(bagdf, Bag)
    assert_eq(get(bagdf.dask, (bagdf.name, 0)), ddf.partitions[0])
    assert_eq(get(bagdf.dask, (bagdf.name, 1)), ddf.partitions[1])


def test_to_records():
    pytest.importorskip("dask.array")
    from dask.array.utils import assert_eq

    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
    )
    ddf = dd.from_pandas(df, 2)

    assert_eq(
        df.to_records(), ddf.to_records(), check_type=False
    )  # TODO: make check_type pass


@pytest.mark.parametrize("lengths", [[2, 2], True])
def test_to_records_with_lengths(lengths):
    pytest.importorskip("dask.array")
    from dask.array.utils import assert_eq

    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
    )
    ddf = dd.from_pandas(df, 2)

    result = ddf.to_records(lengths=lengths)
    assert_eq(df.to_records(), result, check_type=False)  # TODO: make check_type pass

    assert isinstance(result, da.Array)

    expected_chunks = ((2, 2),)

    assert result.chunks == expected_chunks


def test_to_records_raises():
    pytest.importorskip("dask.array")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
    )
    ddf = dd.from_pandas(df, 2)

    with pytest.raises(ValueError):
        ddf.to_records(lengths=[2, 2, 2])
        pytest.fail("3 != 2")

    with pytest.raises(ValueError):
        ddf.to_records(lengths=5)
        pytest.fail("Unexpected value")


def test_from_delayed():
    df = pd.DataFrame(data=np.random.normal(size=(10, 4)), columns=list("abcd"))
    parts = [df.iloc[:1], df.iloc[1:3], df.iloc[3:6], df.iloc[6:10]]
    dfs = [delayed(parts.__getitem__)(i) for i in range(4)]
    meta = dfs[0].compute()

    my_len = lambda x: pd.Series([len(x)])

    for divisions in [None, [0, 1, 3, 6, 10]]:
        ddf = dd.from_delayed(dfs, meta=meta, divisions=divisions)
        assert_eq(ddf, df)
        assert list(ddf.map_partitions(my_len).compute()) == [1, 2, 3, 4]
        assert ddf.known_divisions == (divisions is not None)

        s = dd.from_delayed([d.a for d in dfs], meta=meta.a, divisions=divisions)
        assert_eq(s, df.a)
        assert list(s.map_partitions(my_len).compute()) == [1, 2, 3, 4]
        assert ddf.known_divisions == (divisions is not None)

    meta2 = [(c, "f8") for c in df.columns]
    assert_eq(dd.from_delayed(dfs, meta=meta2), df)
    assert_eq(dd.from_delayed([d.a for d in dfs], meta=("a", "f8")), df.a)

    with pytest.raises(ValueError):
        dd.from_delayed(dfs, meta=meta, divisions=[0, 1, 3, 6])

    with pytest.raises(ValueError) as e:
        dd.from_delayed(dfs, meta=meta.a).compute()
    assert str(e.value).startswith("Metadata mismatch found in `from_delayed`")


def test_from_delayed_to_dask_array():
    # Check that `from_delayed`` can be followed
    # by `to_dask_array` without breaking
    # optimization behavior
    # See: https://github.com/dask-contrib/dask-sql/issues/497
    from dask.blockwise import optimize_blockwise

    dfs = [delayed(pd.DataFrame)(np.ones((3, 2))) for i in range(3)]
    ddf = dd.from_delayed(dfs)
    arr = ddf.to_dask_array()

    # If we optimize this graph without calling
    # `fuse_roots`, the underlying `BlockwiseDep`
    # `mapping` keys will be 1-D (e.g. `(4,)`),
    # while the collection keys will be 2-D
    # (e.g. `(4, 0)`)
    keys = [k[0] for k in arr.__dask_keys__()]
    dsk = optimize_blockwise(arr.dask, keys=keys)
    dsk.cull(keys)

    result = arr.compute()
    assert result.shape == (9, 2)


def test_from_delayed_misordered_meta():
    df = pd.DataFrame(
        columns=["(1)", "(2)", "date", "ent", "val"],
        data=[range(i * 5, i * 5 + 5) for i in range(3)],
        index=range(3),
    )

    # meta with different order for columns
    misordered_meta = pd.DataFrame(
        columns=["date", "ent", "val", "(1)", "(2)"], data=[range(5)]
    )

    ddf = dd.from_delayed([delayed(lambda: df)()], meta=misordered_meta)

    with pytest.raises(ValueError) as info:
        # produces dataframe which does not match meta
        ddf.reset_index().compute(scheduler="sync")
    msg = (
        "The columns in the computed data do not match the columns in the"
        " provided metadata"
    )
    assert msg in str(info.value)


def test_to_delayed():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    # Frame
    a, b = ddf.to_delayed()
    assert isinstance(a, Delayed)
    assert isinstance(b, Delayed)
    assert_eq(a.compute(), df.iloc[:2])

    # Scalar
    x = ddf.x.sum()
    dx = x.to_delayed()
    assert isinstance(dx, Delayed)
    assert_eq(dx.compute(), x)


def test_to_delayed_optimize_graph():
    df = pd.DataFrame({"x": list(range(20))})
    ddf = dd.from_pandas(df, npartitions=20)
    ddf2 = (ddf + 1).loc[:2]

    # Frame
    d = ddf2.to_delayed()[0]
    assert len(d.dask) < 20
    d2 = ddf2.to_delayed(optimize_graph=False)[0]
    assert_eq(ddf2.get_partition(0), d.compute())
    assert_eq(ddf2.get_partition(0), d2.compute())

    # Scalar
    x = ddf2.x.sum()
    dx = x.to_delayed()
    dx2 = x.to_delayed(optimize_graph=False)
    assert_eq(dx.compute(), dx2.compute())


def test_from_dask_array_index_dtype():
    x = da.ones((10,), chunks=(5,))

    df = pd.DataFrame(
        {
            "date": pd.date_range("2019-01-01", periods=10, freq="1min"),
            "val1": list(range(10)),
        }
    )
    ddf = dd.from_pandas(df, npartitions=2).set_index("date")

    ddf2 = dd.from_dask_array(x, index=ddf.index, columns="val2")

    assert ddf.index.dtype == ddf2.index.dtype
    assert ddf.index.name == ddf2.index.name

    df = pd.DataFrame({"idx": np.arange(0, 1, 0.1), "val1": list(range(10))})
    ddf = dd.from_pandas(df, npartitions=2).set_index("idx")

    ddf2 = dd.from_dask_array(x, index=ddf.index, columns="val2")

    assert ddf.index.dtype == ddf2.index.dtype
    assert ddf.index.name == ddf2.index.name


@pytest.mark.parametrize(
    "vals",
    [
        ("A", "B"),
        (3, 4),
        (datetime(2020, 10, 1), datetime(2022, 12, 31)),
    ],
)
def test_from_map_simple(vals):
    # Simple test to ensure required inputs (func & iterable)
    # and basic kwargs work as expected for `from_map`

    def func(input, size=0):
        # Simple function to create Series with a
        # repeated value and index
        value, index = input
        return pd.Series([value] * size, index=[index] * size)

    iterable = [(vals[0], 1), (vals[1], 2)]
    ser = dd.from_map(func, iterable, size=2)
    expect = pd.Series(
        [vals[0], vals[0], vals[1], vals[1]],
        index=[1, 1, 2, 2],
    )

    # Check that result and partition count make sense
    assert ser.npartitions == len(iterable)
    assert_eq(ser, expect)


def test_from_map_multi():
    # Test that `iterables` can contain multiple Iterables

    func = lambda x, y: pd.DataFrame({"add": x + y})
    iterables = (
        [np.arange(2, dtype="int64"), np.arange(2, dtype="int64")],
        [np.array([2, 2], dtype="int64"), np.array([2, 2], dtype="int64")],
    )
    index = np.array([0, 1, 0, 1], dtype="int64")
    expect = pd.DataFrame({"add": np.array([2, 3, 2, 3], dtype="int64")}, index=index)

    ddf = dd.from_map(func, *iterables)
    assert_eq(ddf, expect)


def test_from_map_args():
    # Test that the optional `args` argument works as expected

    func = lambda x, y, z: pd.DataFrame({"add": x + y + z})
    iterable = [np.arange(2, dtype="int64"), np.arange(2, dtype="int64")]
    index = np.array([0, 1, 0, 1], dtype="int64")
    expect = pd.DataFrame({"add": np.array([5, 6, 5, 6], dtype="int64")}, index=index)

    ddf = dd.from_map(func, iterable, args=[2, 3])
    assert_eq(ddf, expect)


def test_from_map_divisions():
    # Test that `divisions` argument works as expected for `from_map`

    func = lambda x: pd.Series([x[0]] * 2, index=range(x[1], x[1] + 2))
    iterable = [("B", 0), ("C", 2)]
    divisions = (0, 2, 4)
    ser = dd.from_map(func, iterable, divisions=divisions)
    expect = pd.Series(
        ["B", "B", "C", "C"],
        index=[0, 1, 2, 3],
    )

    assert ser.divisions == divisions
    assert_eq(ser, expect)


def test_from_map_meta():
    # Test that `meta` can be specified to `from_map`,
    # and that `enforce_metadata` works as expected
    string_dtype = get_string_dtype()

    def func(x, s=0):
        df = pd.DataFrame({"x": [x] * s})
        return df

    iterable = ["A", "B"]

    expect = pd.DataFrame({"x": ["A", "A", "B", "B"]}, index=[0, 1, 0, 1])

    # First Check - Pass in valid metadata
    meta = pd.DataFrame({"x": pd.Series(["A"], dtype=string_dtype)}).iloc[:0]
    ddf = dd.from_map(func, iterable, meta=meta, s=2)
    assert_eq(ddf._meta, meta)
    assert_eq(ddf, expect)

    # Second Check - Pass in invalid metadata
    meta = pd.DataFrame({"a": pd.Series(["A"], dtype=string_dtype)}).iloc[:0]
    ddf = dd.from_map(func, iterable, meta=meta, s=2)
    assert_eq(ddf._meta, meta)

    # Third Check - Pass in invalid metadata,
    # but use `enforce_metadata=False`
    ddf = dd.from_map(func, iterable, meta=meta, enforce_metadata=False, s=2)
    assert_eq(ddf._meta, meta)
    assert_eq(ddf.compute(), expect)


def _generator():
    # Simple generator for test_from_map_other_iterables
    yield from enumerate(["A", "B", "C"])


@pytest.mark.parametrize(
    "iterable",
    [
        enumerate(["A", "B", "C"]),
        ((0, "A"), (1, "B"), (2, "C")),
        _generator(),
    ],
)
def test_from_map_other_iterables(iterable):
    # Test that iterable arguments to `from_map`
    # can be enumerate and generator
    # See: https://github.com/dask/dask/issues/9064

    def func(t):
        size = t[0] + 1
        x = t[1]
        return pd.Series([x] * size)

    ddf = dd.from_map(func, iterable)
    expect = pd.Series(
        ["A", "B", "B", "C", "C", "C"],
        index=[0, 0, 1, 0, 1, 2],
    )
    assert_eq(ddf.compute(), expect)


class MyFunc:
    projected: list[str] = []

    def __init__(self, columns=None):
        self.columns = columns

    def project_columns(self, columns):
        return MyFunc(columns)

    def __call__(self, t, columns=None):
        cols = self.columns or columns
        size = t[0] + 1
        x = t[1]
        df = pd.DataFrame({"A": [x] * size, "B": [10] * size})
        if cols is None:
            return df
        MyFunc.projected.extend(cols)
        return df[cols]


def test_from_map_column_projection():
    # Test that column projection works
    # as expected with from_map when
    # enforce_metadata=True

    ddf = dd.from_map(
        MyFunc(),
        enumerate([0, 1, 2]),
        label="myfunc",
        enforce_metadata=True,
    )
    expect = pd.DataFrame(
        {
            "A": [0, 1, 1, 2, 2, 2],
            "B": [10] * 6,
        },
        index=[0, 0, 1, 0, 1, 2],
    )
    assert_eq(ddf["A"], expect["A"])
    assert set(MyFunc.projected) == {"A"}
    assert_eq(ddf, expect)


@pytest.mark.gpu
@pytest.mark.parametrize("backend", ["pandas", "cudf"])
def test_from_dict_backends(backend):
    _lib = pytest.importorskip(backend)
    with config.set({"dataframe.backend": backend}):
        data = {"a": [1, 2, 3, 4], "B": [10, 11, 12, 13]}
        expected = _lib.DataFrame(data)

        # Check dd.from_dict API
        got = dd.from_dict(data, npartitions=2)
        assert_eq(expected, got)

        # Check from_dict classmethod
        got_classmethod = got.from_dict(data, npartitions=2)
        assert_eq(expected, got_classmethod)


@pytest.mark.parametrize(
    "backend", ["pandas", pytest.param("cudf", marks=pytest.mark.gpu)]
)
def test_sorted_division_locations_duplicates(backend):
    _lib = pytest.importorskip(backend)
    seq = _lib.Series([0, 0, 1, 2])
    divisions, locations = sorted_division_locations(seq, npartitions=2)
    assert divisions == [0, 1, 2]
    assert locations == [0, 2, 4]
