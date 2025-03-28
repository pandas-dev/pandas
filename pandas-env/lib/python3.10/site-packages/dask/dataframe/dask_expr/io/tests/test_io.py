from __future__ import annotations

import glob
import os

import numpy as np
import pyarrow.parquet as pq
import pytest

import dask.array as da
import dask.dataframe as dd
from dask import config
from dask.array.utils import assert_eq as array_assert_eq
from dask.dataframe import read_csv
from dask.dataframe.dask_expr import (
    DataFrame,
    from_array,
    from_dask_array,
    from_dict,
    from_map,
    from_pandas,
    optimize,
    read_parquet,
)
from dask.dataframe.dask_expr._expr import Replace
from dask.dataframe.dask_expr.io import FromArray, FromMap, ReadParquet, parquet
from dask.dataframe.dask_expr.tests._util import _backend_library
from dask.dataframe.utils import assert_eq, make_meta, pyarrow_strings_enabled

# Set DataFrame backend for this module
pd = _backend_library()


def _make_file(dir, format="parquet", df=None):
    fn = os.path.join(str(dir), f"myfile.{format}")
    if df is None:
        df = pd.DataFrame({c: range(10) for c in "abcde"})
    if format == "csv":
        df.to_csv(fn)
    elif format == "parquet":
        df.to_parquet(fn)
    else:
        ValueError(f"{format} not a supported format")
    return fn


def df(fn):
    return read_parquet(fn, columns=["a", "b", "c"])


def df_bc(fn):
    return read_parquet(fn, columns=["b", "c"])


@pytest.mark.parametrize(
    "input,expected",
    [
        (
            # Add -> Mul
            lambda fn: df(fn) + df(fn),
            lambda fn: df(fn) + df(fn),
        ),
        (
            # Column projection
            lambda fn: df(fn)[["b", "c"]],
            lambda fn: read_parquet(fn, columns=["b", "c"]),
        ),
        (
            # Compound
            lambda fn: 3 * (df(fn) + df(fn))[["b", "c"]],
            lambda fn: 3 * (df_bc(fn) + df_bc(fn)),
        ),
        (
            # Traverse Sum
            lambda fn: df(fn).sum()[["b", "c"]],
            lambda fn: df_bc(fn).sum(),
        ),
        (
            # Respect Sum keywords
            lambda fn: df(fn).sum(numeric_only=True)[["b", "c"]],
            lambda fn: df_bc(fn).sum(numeric_only=True),
        ),
    ],
)
def test_simplify(tmpdir, input, expected):
    fn = _make_file(tmpdir, format="parquet")
    result = input(fn).simplify()
    assert str(result.expr) == str(expected(fn).expr)


@pytest.mark.parametrize("fmt", ["parquet", "csv"])
def test_io_fusion(tmpdir, fmt):
    fn = _make_file(tmpdir, format=fmt)
    if fmt == "parquet":
        df = read_parquet(fn)
    else:
        df = read_csv(fn)
    df2 = optimize(df[["a", "b"]] + 1, fuse=True)

    # All tasks should be fused for each partition
    assert len(df2.dask) == df2.npartitions
    assert_eq(df2, df[["a", "b"]] + 1)
    assert all(
        tsk.data_producer
        for k, tsk in df.dask.items()
        if not k[0].startswith("_to_string")
    )


def test_csv_integer_usecols(tmpdir):
    fn = _make_file(tmpdir, format="csv")
    df = read_csv(fn, usecols=[0, 1], names=["x", "y"])
    result = df["x"]
    expected = pd.read_csv(fn, usecols=[0, 1], names=["x", "y"])["x"]
    assert_eq(result, expected)


def test_read_csv_keywords(tmpdir):
    fn = _make_file(tmpdir, format="csv")
    df = read_csv(
        fn,
        sep=",",
        names=["u", "v", "w", "x", "y", "z"],
        engine="pyarrow",
        header=None,
        dtype_backend="pyarrow",
    )
    expected = pd.read_csv(
        fn,
        sep=",",
        names=["u", "v", "w", "x", "y", "z"],
        engine="pyarrow",
        header=None,
        dtype_backend="pyarrow",
    )
    assert_eq(df, expected)


def test_io_fusion_blockwise(tmpdir):
    pdf = pd.DataFrame({c: range(10) for c in "abcdefghijklmn"})
    dd.from_pandas(pdf, 3).to_parquet(tmpdir)
    read_parq = read_parquet(tmpdir)
    df = read_parq["a"].fillna(10).optimize()
    assert df.npartitions == 2
    assert len(df.__dask_graph__()) == 2
    graph = (
        read_parquet(tmpdir)["a"]
        .repartition(npartitions=4)
        .optimize(fuse=False)
        .__dask_graph__()
    )
    assert any(
        f"{read_parq._expr._name.split('-')[0]}-fused" in key[0] for key in graph.keys()
    )


def test_repartition_io_fusion_blockwise(tmpdir):
    pdf = pd.DataFrame({c: range(10) for c in "ab"})
    dd.from_pandas(pdf, 10).to_parquet(tmpdir)
    df = read_parquet(tmpdir)["a"]
    df = df.repartition(npartitions=lambda x: max(x // 2, 1)).optimize()
    assert df.npartitions == 2


def test_io_fusion_merge(tmpdir):
    pdf = pd.DataFrame({c: range(10) for c in "ab"})
    pdf2 = pd.DataFrame({c: range(10) for c in "uvwxyz"})
    dd.from_pandas(pdf, 10).to_parquet(tmpdir)
    dd.from_pandas(pdf2, 10).to_parquet(tmpdir + "x")
    df = read_parquet(tmpdir)
    df2 = read_parquet(tmpdir + "x")
    result = df.merge(df2, left_on="a", right_on="w")[["a", "b", "u"]]
    assert_eq(
        result,
        pdf.merge(pdf2, left_on="a", right_on="w")[["a", "b", "u"]],
        check_index=False,
    )


def test_io_fusion_zero(tmpdir):
    pdf = pd.DataFrame({c: range(10) for c in "ab"})
    dd.from_pandas(pdf, 10).to_parquet(tmpdir)
    result = read_parquet(tmpdir, columns=[])
    assert_eq(result, pdf[[]])


@pytest.mark.parametrize("fmt", ["parquet", "csv", "pandas"])
def test_io_culling(tmpdir, fmt):
    pdf = pd.DataFrame({c: range(10) for c in "abcde"})
    if fmt == "parquet":
        dd.from_pandas(pdf, 2).to_parquet(tmpdir)
        df = read_parquet(tmpdir)
    elif fmt == "csv":
        dd.from_pandas(pdf, 2).to_csv(tmpdir)
        df = read_csv(tmpdir + "/*")
    else:
        df = from_pandas(pdf, 2)

    # Check that original pdf type is conserved
    with config.set({"dataframe.backend": "pandas"}):
        assert type(df.head()) == type(pdf)

    df = (df[["a", "b"]] + 1).partitions[1]
    df2 = optimize(df)

    # All tasks should be fused for the single output partition
    assert df2.npartitions == 1
    # from pandas is not fused together
    assert len(df2.dask) == df2.npartitions + (1 if fmt == "pandas" else 0)
    expected = pdf.iloc[5:][["a", "b"]] + 1
    assert_eq(df2, expected, check_index=False)

    def _check_culling(expr, partitions):
        """CHeck that _partitions is set to the expected value"""
        for dep in expr.dependencies():
            _check_culling(dep, partitions)
        if "_partitions" in expr._parameters:
            assert expr._partitions == partitions

    # Check that we still get culling without fusion
    df3 = optimize(df, fuse=False)
    _check_culling(df3.expr, [1])
    assert_eq(df3, expected, check_index=False)


def test_parquet_complex_filters(tmpdir):
    df = read_parquet(_make_file(tmpdir))
    pdf = df.compute()
    got = df["a"][df["b"] > df["b"].mean()]
    expect = pdf["a"][pdf["b"] > pdf["b"].mean()]

    assert_eq(got, expect)
    assert_eq(got.optimize(), expect)
    assert all(tsk.data_producer for tsk in df.dask.values())


@pytest.mark.parametrize("optimize", [True, False])
def test_to_dask_array(optimize):
    pdf = pd.DataFrame({"x": [1, 4, 3, 2, 0, 5]})
    df = from_pandas(pdf, npartitions=2)
    darr = df.to_dask_array(optimize=optimize)
    assert isinstance(darr, da.Array)
    array_assert_eq(darr, pdf.values)


@pytest.mark.parametrize(
    "fmt,read_func,read_cls",
    [("parquet", read_parquet, ReadParquet), ("csv", read_csv, FromMap)],
)
def test_combine_similar(tmpdir, fmt, read_func, read_cls):
    pdf = pd.DataFrame(
        {"x": [0, 1, 2, 3] * 4, "y": range(16), "z": [None, 1, 2, 3] * 4}
    )
    fn = _make_file(tmpdir, format=fmt, df=pdf)
    df = read_func(fn)
    df = df.replace(1, 100)
    df["xx"] = df.x != 0
    df["yy"] = df.y != 0
    got = df[["xx", "yy", "x"]].sum()

    pdf = pdf.replace(1, 100)
    pdf["xx"] = pdf.x != 0
    pdf["yy"] = pdf.y != 0
    expect = pdf[["xx", "yy", "x"]].sum()

    # Check correctness
    assert_eq(got, expect)
    assert_eq(got.optimize(fuse=False), expect)
    assert_eq(got.optimize(fuse=True), expect)

    # We should only have one ReadParquet/ReadCSV node,
    # and it should not include "z" in the column projection
    read_nodes = list(got.optimize(fuse=False).find_operations(read_cls))
    assert len(read_nodes) == 1
    assert set(read_nodes[0].columns) == {"x", "y"}

    # All Replace operations should also be the same
    replace_nodes = list(got.optimize(fuse=False).find_operations(Replace))
    assert len(replace_nodes) == 1


def test_combine_similar_no_projection_on_one_branch(tmpdir):
    pdf = pd.DataFrame(
        {"x": [0, 1, 2, 3] * 4, "y": range(16), "z": [None, 1, 2, 3] * 4}
    )
    fn = _make_file(tmpdir, format="parquet", df=pdf)
    df = read_parquet(fn)
    df["xx"] = df.x != 0

    pdf["xx"] = pdf.x != 0
    assert_eq(df, pdf)


@pytest.mark.parametrize("meta", [True, False])
@pytest.mark.parametrize("label", [None, "foo"])
@pytest.mark.parametrize("allow_projection", [True, False])
@pytest.mark.parametrize("enforce_metadata", [True, False])
def test_from_map(tmpdir, meta, label, allow_projection, enforce_metadata):
    pdf = pd.DataFrame({c: range(10) for c in "abcdefghijklmn"})
    dd.from_pandas(pdf, 3).to_parquet(tmpdir, write_index=False)
    files = sorted(glob.glob(str(tmpdir) + "/*.parquet"))
    if allow_projection:
        func = pd.read_parquet
    else:
        func = lambda *args, **kwargs: pd.read_parquet(*args, **kwargs)
    options = {
        "enforce_metadata": enforce_metadata,
        "label": label,
    }
    if meta:
        options["meta"] = pdf.iloc[:0]

    df = from_map(func, files, **options)
    assert_eq(df, pdf, check_index=False)
    assert_eq(df["a"], pdf["a"], check_index=False)
    assert_eq(df[["a"]], pdf[["a"]], check_index=False)
    assert_eq(df[["a", "b"]], pdf[["a", "b"]], check_index=False)
    assert all(
        tsk.data_producer
        for key, tsk in df.__dask_graph__().items()
        if not key[0].startswith("_to_string")
    )
    if allow_projection:
        assert all(
            tsk.data_producer
            for key, tsk in df[["a"]].optimize(fuse=False).__dask_graph__().items()
            if not key[0].startswith("_to_string")
        )

    if label:
        if pyarrow_strings_enabled():
            assert df.expr.frame._name.startswith(label)
        else:
            assert df.expr._name.startswith(label)

    if allow_projection:
        got = df[["a", "b"]].optimize(fuse=False)
        if pyarrow_strings_enabled():
            assert isinstance(got.expr.frame, FromMap)
            assert got.expr.frame.operand("columns") == ["a", "b"]
        else:
            assert isinstance(got.expr, FromMap)
            assert got.expr.operand("columns") == ["a", "b"]

    # Check that we can always pass columns up front
    if meta:
        options["meta"] = options["meta"][["a", "b"]]
    result = from_map(func, files, columns=["a", "b"], **options)
    assert_eq(result, pdf[["a", "b"]], check_index=False)
    if meta:
        options["meta"] = options["meta"][["a"]]
    result = from_map(func, files, columns="a", **options)
    assert_eq(result, pdf[["a"]], check_index=False)

    # Check the case that func returns a Series
    if meta:
        options["meta"] = options["meta"]["a"]
    result = from_map(lambda x: pd.read_parquet(x)["a"], files, **options)
    assert_eq(result, pdf["a"], check_index=False)


def func(path, columns):
    raise NotImplementedError("This shouldn't ever be called")


def func2(path, columns):
    return pd.DataFrame({"a": range(10), "b": range(10)})


def test_from_map_columns_required():
    with pytest.raises(TypeError, match=r"columns.*optional"):
        from_map(func, ["foo"])
    meta = make_meta(func2("foo", None))
    ddf = from_map(func2, ["foo"], meta=meta)
    ddf_expected = from_pandas(func2("foo", None), npartitions=1)

    assert_eq(ddf, ddf_expected, check_divisions=False)

    actual = from_map(func2, ["foo"], meta=meta)[["a"]].optimize()
    expected = from_map(func2, ["foo"], meta=meta, columns=["a"]).optimize()

    assert actual._name == expected._name


def func_object(path):
    return pd.DataFrame({"a": "x", "b": range(10)})


def test_from_map_string_conversion():
    result = from_map(func_object, ["foo"])
    dtype = "string" if pyarrow_strings_enabled() else "object"
    assert result.a.dtype == dtype
    assert result.a.compute().dtype == dtype


def test_from_array():
    arr = np.random.randint(1, 100, (100,))
    assert_eq(from_array(arr, chunksize=5), pd.Series(arr))
    assert_eq(from_array(arr, chunksize=5, columns="a"), pd.Series(arr, name="a"))
    columns = ["a", "b", "c", "d"]
    arr = np.random.randint(1, 100, (100, 4))
    assert_eq(
        from_array(arr, chunksize=5, columns=columns),
        pd.DataFrame(arr, columns=columns),
    )
    assert_eq(
        from_array(arr, chunksize=5, columns=columns)["a"],
        pd.DataFrame(arr, columns=columns)["a"],
    )
    assert_eq(
        from_array(arr, chunksize=5, columns=columns)[["a"]],
        pd.DataFrame(arr, columns=columns)[["a"]],
    )
    q = from_array(arr, chunksize=5, columns=columns)[["a"]].simplify()
    for expr in q.walk():
        if isinstance(expr, FromArray):
            assert expr.columns == ["a"]


@pytest.mark.parametrize("dtype", [object, str])
def test_from_array_string_conersion(dtype):
    arr = np.array(["a", "b", "c", "d"], dtype=dtype)
    result = from_array(arr, chunksize=2)
    dtype = "string" if pyarrow_strings_enabled() else "object"
    assert result.dtype == dtype
    assert result.compute().dtype == dtype


def test_from_dask_array():
    import dask.array as da

    arr = da.ones((20, 4), chunks=(2, 2))
    df = from_dask_array(arr, columns=["a", "b", "c", "d"])
    pdf = pd.DataFrame(arr.compute(), columns=["a", "b", "c", "d"])
    assert_eq(df, pdf)


@pytest.mark.parametrize("columns", [1, "aa"])
def test_from_dask_array_scalar_columns(columns):
    import dask.array as da

    arr = da.ones((20,), chunks=(2,))
    df = from_dask_array(arr, columns=columns)
    pdf = pd.Series(arr.compute(), name=columns)
    assert_eq(df, pdf)


def test_from_dask_array_projection():
    rng = np.random.default_rng()
    arr_np = rng.random((100, 10))
    arr = da.from_array(arr_np, chunks=(50, 10))
    pdf = pd.DataFrame(arr_np)
    df = from_dask_array(arr)
    # Check getitem[np.int64(0)]
    dd.assert_eq(pdf[pdf.columns[0]], df[df.columns[0]])
    # Check loc[:, np.int64(0)]
    dd.assert_eq(pdf.loc[:, pdf.columns[0]], df.loc[:, df.columns[0]])


def test_from_dict():
    data = {"a": [1, 2, 3, 4], "B": [10, 11, 12, 13]}
    result = from_dict(data, npartitions=2)
    expected = pd.DataFrame(data)
    assert_eq(result, expected)
    assert_eq(DataFrame.from_dict(data), expected)


@pytest.mark.parametrize("normalizer", ("filemetadata", "schema"))
def test_normalize_token_parquet_filemetadata_and_schema(tmpdir, normalizer):
    df = pd.DataFrame({c: range(10) for c in "abcde"})
    path = tmpdir / "df.parquet"
    df.to_parquet(path)

    meta = pq.read_metadata(path)

    if normalizer == "filemetadata":
        normalizer = parquet.normalize_pq_filemetadata
        obj = meta
    else:
        normalizer = parquet.normalize_pq_schema
        obj = meta.schema

    assert normalizer(obj) == normalizer(obj)


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


def test_from_array_partition_pruning():
    arr = np.random.random(size=(200, 2))
    df = dd.from_array(arr, chunksize=100)
    result = df.partitions[1]
    expected = arr[100:]
    assert_eq(result, pd.DataFrame(expected, index=list(range(100, 200))))
