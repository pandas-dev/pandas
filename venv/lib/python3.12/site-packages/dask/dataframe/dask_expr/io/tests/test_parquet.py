from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd
import pytest
from pyarrow import fs

import dask
from dask._compatibility import WINDOWS
from dask.dataframe._compat import PYARROW_GE_1500
from dask.dataframe.dask_expr import from_array, from_graph, from_pandas, read_parquet
from dask.dataframe.dask_expr._expr import Filter, Lengths, Literal
from dask.dataframe.dask_expr._reductions import Len
from dask.dataframe.dask_expr.io import FusedParquetIO, ReadParquet
from dask.dataframe.dask_expr.io.parquet import (
    _aggregate_statistics_to_file,
    _combine_stats,
    _extract_stats,
)
from dask.dataframe.utils import assert_eq
from dask.utils import key_split


def _make_file(dir, df=None, filename="myfile.parquet", **kwargs):
    fn = os.path.join(str(dir), filename)
    if df is None:
        df = pd.DataFrame({c: range(10) for c in "abcde"})
    df.to_parquet(fn, **kwargs)
    return fn


@pytest.fixture
def parquet_file(tmpdir):
    return _make_file(tmpdir)


@pytest.fixture(
    params=[
        pytest.param(
            "arrow",
            marks=pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0"),
        ),
        "fsspec",
    ]
)
def filesystem(request):
    return request.param


def test_parquet_len(tmpdir, filesystem):
    df = read_parquet(_make_file(tmpdir), filesystem=filesystem)
    pdf = df.compute()

    assert len(df[df.a > 5]) == len(pdf[pdf.a > 5])

    s = (df["b"] + 1).astype("Int32")
    assert len(s) == len(pdf)

    assert isinstance(Len(s.expr).optimize(), Literal)
    assert isinstance(Lengths(s.expr).optimize(), Literal)


def test_parquet_missing_stats(tmpdir, filesystem):
    _make_file(tmpdir)
    _make_file(tmpdir, write_statistics=["a", "b"], filename="bla.parquet")

    result = read_parquet(tmpdir, filesystem=filesystem)
    expected = pd.concat(
        [
            pd.DataFrame({c: range(10) for c in "abcde"}),
            pd.DataFrame({c: range(10) for c in "abcde"}),
        ]
    )
    assert_eq(result, expected, check_index=False)


@pytest.mark.parametrize("val", [np.nan, 1])
def test_parquet_all_na_column(tmpdir, filesystem, val):
    pdf = pd.DataFrame({"a": [np.nan] * 299 + [val], "b": [1, 2, 3] * 100})
    _make_file(tmpdir, df=pdf, filename="bla.parquet", row_group_size=100)
    result = read_parquet(tmpdir, filesystem=filesystem)
    assert_eq(result, pdf)


def test_parquet_len_filter(tmpdir, filesystem):
    df = read_parquet(_make_file(tmpdir), filesystem=filesystem)
    expr = Len(df[df.c > 0].expr)
    result = expr.simplify()
    for rp in result.find_operations(ReadParquet):
        assert rp.operand("columns") == ["c"] or rp.operand("columns") == []


def test_parquet_len_empty_dir():
    assert len(read_parquet(["empty_dir/*.parquet"])) == 0


@pytest.mark.skipif(WINDOWS, reason="directory empty on windows")
@pytest.mark.parametrize("write_metadata_file", [True, False])
def test_to_parquet(tmpdir, write_metadata_file):
    pdf = pd.DataFrame({"x": [1, 4, 3, 2, 0, 5]})
    df = from_pandas(pdf, npartitions=2)

    # Check basic parquet round trip
    df.to_parquet(tmpdir, write_metadata_file=write_metadata_file)
    df2 = read_parquet(tmpdir, calculate_divisions=True)
    assert_eq(df, df2)

    # Check overwrite behavior
    df["new"] = df["x"] + 1
    df.to_parquet(tmpdir, overwrite=True, write_metadata_file=write_metadata_file)
    df2 = read_parquet(tmpdir, calculate_divisions=True)
    assert_eq(df, df2)

    # Check that we cannot overwrite a path we are
    # reading from in the same graph
    with pytest.raises(ValueError, match="Cannot overwrite"):
        df2.to_parquet(tmpdir, overwrite=True)


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_pyarrow_filesystem(parquet_file):
    filesystem = fs.LocalFileSystem()

    df_pa = read_parquet(parquet_file, filesystem=filesystem)
    df = read_parquet(parquet_file)
    assert assert_eq(df, df_pa)
    assert all(tsk.data_producer for tsk in df.dask.values())
    assert all(tsk.data_producer for tsk in df.optimize().dask.values())


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
@pytest.mark.parametrize("dtype_backend", ["pyarrow", "numpy_nullable", None])
def test_pyarrow_filesystem_dtype_backend(parquet_file, dtype_backend):
    filesystem = fs.LocalFileSystem()

    df_pa = read_parquet(
        parquet_file, filesystem=filesystem, dtype_backend=dtype_backend
    )
    df = read_parquet(parquet_file, dtype_backend=dtype_backend)
    assert assert_eq(df, df_pa)


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
@pytest.mark.parametrize("types_mapper", [None, lambda x: None])
def test_pyarrow_filesystem_types_mapper(parquet_file, types_mapper):
    # This test isn't doing much other than ensuring the stuff is not raising
    # anywhere
    filesystem = fs.LocalFileSystem()

    df_pa = read_parquet(
        parquet_file,
        filesystem=filesystem,
        arrow_to_pandas={"types_mapper": types_mapper},
    )
    df = read_parquet(parquet_file, arrow_to_pandas={"types_mapper": types_mapper})
    assert assert_eq(df, df_pa)


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_pyarrow_filesystem_serialize(parquet_file):
    filesystem = fs.LocalFileSystem()

    df_pa = read_parquet(parquet_file, filesystem=filesystem)

    roundtripped = pickle.loads(pickle.dumps(df_pa.optimize().dask))
    roundtripped_df = from_graph(
        roundtripped,
        df_pa._meta,
        df_pa.divisions,
        df_pa.__dask_keys__(),
        key_split(df_pa._name),
    )
    assert assert_eq(df_pa, roundtripped_df)


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_pyarrow_filesystem_filters(parquet_file):
    filesystem = fs.LocalFileSystem()

    df_pa = read_parquet(parquet_file, filesystem=filesystem)
    df_pa = df_pa[df_pa.c == 1]
    expected = read_parquet(
        parquet_file, filesystem=filesystem, filters=[[("c", "==", 1)]]
    )
    assert df_pa.optimize()._name == expected.optimize()._name
    assert len(df_pa.compute()) == 1


second_parquet_file = parquet_file


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_pyarrow_filesystem_list_of_files(parquet_file, second_parquet_file):
    filesystem = fs.LocalFileSystem()

    result = read_parquet([parquet_file, second_parquet_file], filesystem=filesystem)
    expected = pd.read_parquet([parquet_file, second_parquet_file])
    assert_eq(result, expected, check_index=False)


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_partition_pruning(tmpdir):
    with dask.config.set({"dataframe.parquet.minimum-partition-size": 1}):
        filesystem = fs.LocalFileSystem()
        df = from_pandas(
            pd.DataFrame(
                {
                    "a": [1, 2, 3, 4, 5] * 10,
                    "b": range(50),
                }
            ),
            npartitions=2,
        )
        df.to_parquet(tmpdir, partition_on=["a"])
        ddf = read_parquet(tmpdir, filesystem=filesystem)
        ddf_filtered = read_parquet(
            tmpdir, filters=[[("a", "==", 1)]], filesystem=filesystem
        )
        assert ddf_filtered.npartitions == ddf.npartitions // 5

        ddf_optimize = read_parquet(tmpdir, filesystem=filesystem)
        ddf_optimize = ddf_optimize[ddf_optimize.a == 1].optimize()
        assert ddf_optimize.npartitions == ddf.npartitions // 5
        assert_eq(
            ddf_filtered,
            ddf_optimize,
            # FIXME ?
            check_names=False,
        )


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_predicate_pushdown(tmpdir):
    original = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5] * 10,
            "b": [0, 1, 2, 3, 4] * 10,
            "c": range(50),
            "d": [6, 7] * 25,
            "e": [8, 9] * 25,
        }
    )
    fn = _make_file(tmpdir, df=original)
    df = read_parquet(fn, filesystem="arrow")
    assert_eq(df, original)
    x = df[df.a == 5][df.c > 20]["b"]
    y = x.optimize(fuse=False)
    assert isinstance(y.expr.frame, FusedParquetIO)
    assert ("a", "==", 5) in y.expr.frame.operands[0].operand("filters")[0]
    assert ("c", ">", 20) in y.expr.frame.operands[0].operand("filters")[0]
    assert list(y.columns) == ["b"]

    # Check computed result
    y_result = y.compute()
    assert y_result.name == "b"
    assert len(y_result) == 6
    assert (y_result == 4).all()

    # Don't push down if replace is in there
    x = df[df.replace(5, 50).a == 5]["b"]
    y = x.optimize(fuse=False)
    assert isinstance(y.expr, Filter)
    assert len(y.compute()) == 0


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_predicate_pushdown_compound(tmpdir):
    pdf = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5] * 10,
            "b": [0, 1, 2, 3, 4] * 10,
            "c": range(50),
            "d": [6, 7] * 25,
            "e": [8, 9] * 25,
        }
    )
    fn = _make_file(tmpdir, df=pdf)
    df = read_parquet(fn, filesystem="arrow")

    # Test AND
    x = df[(df.a == 5) & (df.c > 20)]["b"]
    y = x.optimize(fuse=False)
    assert isinstance(y.expr.frame, FusedParquetIO)
    assert {("c", ">", 20), ("a", "==", 5)} == set(y.expr.frame.operands[0].filters[0])
    assert_eq(
        y,
        pdf[(pdf.a == 5) & (pdf.c > 20)]["b"],
        check_index=False,
    )

    # Test OR
    x = df[(df.a == 5) | (df.c > 20)]
    x = x[x.b != 0]["b"]
    y = x.optimize(fuse=False)
    assert isinstance(y.expr.frame, FusedParquetIO)
    filters = [
        set(y.expr.frame.operands[0].filters[0]),
        set(y.expr.frame.operands[0].filters[1]),
    ]
    assert {("c", ">", 20), ("b", "!=", 0)} in filters
    assert {("a", "==", 5), ("b", "!=", 0)} in filters
    expect = pdf[(pdf.a == 5) | (pdf.c > 20)]
    expect = expect[expect.b != 0]["b"]
    assert_eq(
        y,
        expect,
        check_index=False,
    )

    # Test OR and AND
    x = df[((df.a == 5) | (df.c > 20)) & (df.b != 0)]["b"]
    z = x.optimize(fuse=False)
    assert isinstance(z.expr.frame, FusedParquetIO)
    filters = [
        set(z.expr.frame.operands[0].filters[0]),
        set(z.expr.frame.operands[0].filters[1]),
    ]
    assert {("c", ">", 20), ("b", "!=", 0)} in filters
    assert {("a", "==", 5), ("b", "!=", 0)} in filters
    assert_eq(y, z)


@pytest.mark.xfail(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_aggregate_rg_stats_to_file(tmpdir):
    filesystem = fs.LocalFileSystem()
    fn = str(tmpdir)
    ddf = from_pandas(pd.DataFrame({"a": range(10)}), npartitions=1)
    ddf.to_parquet(fn)
    ddf = read_parquet(fn, filesystem=filesystem)
    frag = ddf._expr.fragments[0]
    # Make sure this doesn't raise. We'll test the actual aggregation below
    _aggregate_statistics_to_file([frag.metadata.to_dict()])
    # In reality, we'll strip the metadata down
    assert (
        len(_aggregate_statistics_to_file([_extract_stats(frag.metadata.to_dict())]))
        > 0
    )


def test_aggregate_statistics_to_file():
    file_in = {
        "top-level-file-stat": "not-interested",
        "row_groups": [
            # RG 1
            {
                "num_rows": 100,
                "total_byte_size": 1000,
                "columns": [
                    {
                        "total_compressed_size": 10,
                        "total_uncompressed_size": 15,
                        "statistics": {
                            "min": 0,
                            "max": 10,
                            "null_count": 5,
                            "distinct_count": None,
                            "new_powerful_yet_uknown_stat": 42,
                        },
                        "path_in_schema": "a",
                    },
                    {
                        "total_compressed_size": 7,
                        "total_uncompressed_size": 23,
                        "statistics": {
                            "min": 11,
                            "max": 20,
                            "null_count": 1,
                            "distinct_count": 12,
                            "new_powerful_yet_uknown_stat": 42,
                        },
                        "path_in_schema": "b",
                    },
                ],
            },
            # RG 2
            {
                "num_rows": 50,
                "total_byte_size": 500,
                "columns": [
                    {
                        "total_compressed_size": 5,
                        "total_uncompressed_size": 7,
                        "statistics": {
                            "min": 40,
                            "max": 50,
                            "null_count": 0,
                            "distinct_count": None,
                            "new_powerful_yet_uknown_stat": 42,
                        },
                        "path_in_schema": "a",
                    },
                    {
                        "total_compressed_size": 7,
                        "total_uncompressed_size": 23,
                        "statistics": {
                            "min": 0,
                            "max": 20,
                            "null_count": 0,
                            "distinct_count": None,
                            "new_powerful_yet_uknown_stat": 42,
                        },
                        "path_in_schema": "b",
                    },
                ],
            },
        ],
    }

    expected = {
        "top-level-file-stat": "not-interested",
        "num_rows": 150,
        "total_byte_size": 1500,
        "columns": [
            {
                "path_in_schema": "a",
                "total_compressed_size": 15,
                "total_uncompressed_size": 22,
                "statistics": {
                    "min": 0,
                    "max": 50,
                },
            },
            {
                "path_in_schema": "b",
                "total_compressed_size": 14,
                "total_uncompressed_size": 46,
                "statistics": {
                    "min": 0,
                    "max": 20,
                },
            },
        ],
    }
    actual = _aggregate_statistics_to_file([file_in])
    assert len(actual) == 1
    assert actual[0] == expected


def test_combine_statistics():
    file_in = [
        # File 1
        {
            "top-level-file-stat": "not-interested",
            "num_row_groups": 2,
            "row_groups": [
                # RG 1
                {
                    "num_rows": 200,
                    "total_byte_size": 1000,
                    "columns": [
                        {
                            "total_compressed_size": 10,
                            "total_uncompressed_size": 15,
                            "statistics": {
                                "min": 0,
                                "max": 10,
                                "null_count": 5,
                                "distinct_count": None,
                                "new_powerful_yet_uknown_stat": 42,
                            },
                            "path_in_schema": "a",
                        },
                        {
                            "total_compressed_size": 7,
                            "total_uncompressed_size": 23,
                            "statistics": {
                                "min": 11,
                                "max": 20,
                                "null_count": 1,
                                "distinct_count": 12,
                                "new_powerful_yet_uknown_stat": 42,
                            },
                            "path_in_schema": "b",
                        },
                    ],
                },
                # RG 2
                {
                    "num_rows": 50,
                    "total_byte_size": 500,
                    "columns": [
                        {
                            "total_compressed_size": 5,
                            "total_uncompressed_size": 7,
                            "statistics": {
                                "min": 40,
                                "max": 50,
                                "null_count": 0,
                                "distinct_count": None,
                                "new_powerful_yet_uknown_stat": 42,
                            },
                            "path_in_schema": "a",
                        },
                        {
                            "total_compressed_size": 7,
                            "total_uncompressed_size": 23,
                            "statistics": {
                                "min": 0,
                                "max": 20,
                                "null_count": 0,
                                "distinct_count": None,
                                "new_powerful_yet_uknown_stat": 42,
                            },
                            "path_in_schema": "b",
                        },
                    ],
                },
            ],
        },
        # File 2
        {
            "top-level-file-stat": "not-interested",
            "num_row_groups": 1,
            "row_groups": [
                # RG 1
                {
                    "num_rows": 100,
                    "total_byte_size": 2000,
                    "columns": [
                        {
                            "total_compressed_size": 10,
                            "total_uncompressed_size": 15,
                            "statistics": {
                                "min": 0,
                                "max": 10,
                                "null_count": 5,
                                "distinct_count": None,
                                "new_powerful_yet_uknown_stat": 42,
                            },
                            "path_in_schema": "a",
                        },
                        {
                            "total_compressed_size": 7,
                            "total_uncompressed_size": 23,
                            "statistics": {
                                "min": 11,
                                "max": 20,
                                "null_count": 1,
                                "distinct_count": 12,
                                "new_powerful_yet_uknown_stat": 42,
                            },
                            "path_in_schema": "b",
                        },
                    ],
                },
            ],
        },
    ]
    actual = _combine_stats(file_in)
    expected = {
        "num_rows": ((200 + 50) + 100) // 2,
        "total_byte_size": ((1000 + 500) + 2000) // 2,
        "num_row_groups": 1.5,
        "columns": [
            {
                "total_compressed_size": ((10 + 5) + 10) / 2,
                "total_uncompressed_size": ((15 + 7) + 15) / 2,
                "path_in_schema": "a",
            },
            {
                "total_compressed_size": ((7 + 7) + 7) / 2,
                "total_uncompressed_size": ((23 + 23) + 23) / 2,
                "path_in_schema": "b",
            },
        ],
    }
    assert actual == expected


def test_index_only_from_parquet(tmpdir):
    pdf = pd.DataFrame({"foo": range(5)}, index=range(50, 55))
    pdf.to_parquet(tmpdir + "/test.parquet")
    result = read_parquet(tmpdir + "/test.parquet").index
    assert_eq(result, pdf.index)


def test_timestamp_divisions(tmpdir):
    pdf = pd.DataFrame.from_dict(
        {
            "Date": ["11/26/2017", "11/26/2017"],
            "Time": ["17:00:00.067", "17:00:00.102"],
            "Volume": [403, 3],
        }
    )
    pdf["Timestamp"] = pd.to_datetime(pdf.Date) + pd.to_timedelta(pdf.Time)
    pdf.to_parquet(tmpdir + "/test.parquet")
    df = read_parquet(
        tmpdir + "/test.parquet", index="Timestamp", calculate_divisions=True
    )["Volume"]
    assert df.optimize().divisions == (
        pd.Timestamp("2017-11-26 17:00:00.067000"),
        pd.Timestamp("2017-11-26 17:00:00.102000"),
    )


def test_read_parquet_index_projection(tmpdir):
    df = from_array(
        np.zeros((201 * 10, 8), dtype=np.int64),
        columns=["val1", "val2", "val3", "val4", "val5", "val6", "val7", "tsprv"],
    )
    df = df.repartition(npartitions=201)
    df.to_parquet(tmpdir + "/index", write_index=False)

    df = read_parquet(tmpdir + "/index", calculate_divisions=True)
    result = df.assign(dts=df.index - df.tsprv)
    expected = df.compute()
    expected = expected.assign(dts=expected.index - expected.tsprv)

    assert_eq(result.dts.min(), expected.dts.min())


@pytest.mark.filterwarnings(
    # The stats collection encounters the async client and falls back to a
    # sync scheduler. Irrelevant for this test.
    "ignore:.*Client.*asynchronous:UserWarning",
)
def test_ensure_plan_computed_during_optimization(tmpdir, filesystem):

    pytest.importorskip("distributed")
    from distributed.utils_test import gen_cluster

    @gen_cluster(
        client=True,
        clean_kwargs={"threads": False},
        config={"admin.async-client-fallback": "sync"},
    )
    async def _test(c, s, *workers):
        # https://github.com/dask/dask/issues/11932

        df = read_parquet(_make_file(tmpdir), filesystem=filesystem)
        df = df[df.a == 1]
        res = await c.compute(df.size)
        assert res == 5

    _test()
