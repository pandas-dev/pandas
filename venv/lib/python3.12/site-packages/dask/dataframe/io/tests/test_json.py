from __future__ import annotations

import json
import os

import fsspec
import pandas as pd
import pytest
from packaging.version import Version

import dask
import dask.dataframe as dd
from dask.dataframe.utils import assert_eq
from dask.utils import tmpdir, tmpfile

df = pd.DataFrame({"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]})
ddf = dd.from_pandas(df, npartitions=2)


@pytest.mark.parametrize("orient", ["split", "records", "index", "columns", "values"])
def test_read_json_with_path_column(orient):
    with tmpfile("json") as f:
        df.to_json(f, orient=orient, lines=False)
        actual = dd.read_json(f, orient=orient, lines=False, include_path_column=True)
        actual_pd = pd.read_json(f, orient=orient, lines=False)
        # The default column name when include_path_colum is True is "path"
        # The paths on Windows are converted to forward slash somewhere in the file
        # reading chain in Dask, so we have to do the same here.
        actual_pd["path"] = pd.Series(
            (f.replace(os.sep, "/"),) * len(actual_pd), dtype="category"
        )
        assert actual.path.dtype == "category"
        assert_eq(actual, actual_pd)


def test_read_json_path_column_with_duplicate_name_is_error():
    with tmpfile("json") as f:
        df.to_json(f, orient="records", lines=False)
        with pytest.raises(ValueError, match="Files already contain"):
            dd.read_json(f, orient="records", lines=False, include_path_column="x")


def test_read_json_with_path_converter():
    path_column_name = "filenames"

    def path_converter(x):
        return "asdf.json"

    with tmpfile("json") as f:
        df.to_json(f, orient="records", lines=False)
        actual = dd.read_json(
            f,
            orient="records",
            lines=False,
            include_path_column=path_column_name,
            path_converter=path_converter,
        )
        actual_pd = pd.read_json(f, orient="records", lines=False)
        actual_pd[path_column_name] = pd.Series(
            (path_converter(f),) * len(actual_pd), dtype="category"
        )
        assert_eq(actual, actual_pd)


def test_read_orient_not_records_and_lines():
    with pytest.raises(ValueError, match="Line-delimited JSON"):
        dd.read_json("nofile.json", orient="split", lines=True)


def test_write_orient_not_records_and_lines():
    with tmpfile("json") as f:
        with pytest.raises(ValueError, match="Line-delimited JSON"):
            dd.to_json(ddf, f, orient="split", lines=True)


@pytest.mark.parametrize("blocksize", [5, 15, 33, 200, 90000])
def test_read_json_multiple_files_with_path_column(blocksize, tmpdir):
    fil1 = str(tmpdir.join("fil1.json")).replace(os.sep, "/")
    fil2 = str(tmpdir.join("fil2.json")).replace(os.sep, "/")
    df = pd.DataFrame({"x": range(5), "y": ["a", "b", "c", "d", "e"]})
    df2 = df.assign(x=df.x + 0.5)
    orient = "records"
    lines = True
    df.to_json(fil1, orient=orient, lines=lines)
    df2.to_json(fil2, orient=orient, lines=lines)
    path_dtype = pd.CategoricalDtype((fil1, fil2))
    df["path"] = pd.Series((fil1,) * len(df), dtype=path_dtype)
    df2["path"] = pd.Series((fil2,) * len(df2), dtype=path_dtype)
    sol = pd.concat([df, df2])
    res = dd.read_json(
        str(tmpdir.join("fil*.json")),
        orient=orient,
        lines=lines,
        include_path_column=True,
        blocksize=blocksize,
    )
    assert_eq(res, sol, check_index=False)


@pytest.mark.parametrize("orient", ["split", "records", "index", "columns", "values"])
def test_read_json_basic(orient):
    with tmpfile("json") as f:
        df.to_json(f, orient=orient, lines=False)
        actual = dd.read_json(f, orient=orient, lines=False)
        actual_pd = pd.read_json(f, orient=orient, lines=False)

        assert_eq(actual, actual_pd)
        if orient == "values":
            actual.columns = list(df.columns)
        assert_eq(actual, df)


@pytest.mark.parametrize("fkeyword", ["pandas", "json"])
def test_read_json_fkeyword(fkeyword):
    def _my_json_reader(*args, **kwargs):
        if fkeyword == "json":
            return pd.DataFrame.from_dict(json.load(*args))
        return pd.read_json(*args)

    with tmpfile("json") as f:
        df.to_json(f, orient="records", lines=False)
        actual = dd.read_json(f, orient="records", lines=False, engine=_my_json_reader)
        actual_pd = pd.read_json(f, orient="records", lines=False)
        assert_eq(actual, actual_pd)


@pytest.mark.parametrize("engine", ["ujson", pd.read_json])
def test_read_json_engine_str(engine):
    with tmpfile("json") as f:
        df.to_json(f, lines=False)
        got = dd.read_json(f, engine=engine, lines=False)
        assert_eq(got, df)


@pytest.mark.parametrize("orient", ["split", "records", "index", "columns", "values"])
def test_read_json_meta(orient, tmpdir):
    df = pd.DataFrame({"x": range(5), "y": ["a", "b", "c", "d", "e"]})
    df2 = df.assign(x=df.x + 0.5)
    lines = orient == "records"
    df.to_json(str(tmpdir.join("fil1.json")), orient=orient, lines=lines)
    df2.to_json(str(tmpdir.join("fil2.json")), orient=orient, lines=lines)
    sol = pd.concat([df, df2])
    meta = df2.iloc[:0]

    if orient == "values":
        # orient=values loses column names
        sol.columns = meta.columns = [0, 1]

    res = dd.read_json(
        str(tmpdir.join("fil*.json")), orient=orient, meta=meta, lines=lines
    )
    assert_eq(res, sol)

    if orient == "records":
        # Also check chunked version
        res = dd.read_json(
            str(tmpdir.join("fil*.json")),
            orient=orient,
            meta=meta,
            lines=True,
            blocksize=50,
        )
        assert_eq(res, sol, check_index=False)


@pytest.mark.parametrize("orient", ["split", "records", "index", "columns", "values"])
def test_write_json_basic(orient):
    with tmpdir() as path:
        fn = os.path.join(path, "1.json")
        df.to_json(fn, orient=orient, lines=False)
        actual = dd.read_json(fn, orient=orient, lines=False)
        if orient == "values":
            actual.columns = list(df.columns)
        assert_eq(actual, df)


def test_to_json_with_get():
    from dask.multiprocessing import get as mp_get

    flag = [False]

    def my_get(*args, **kwargs):
        flag[0] = True
        return mp_get(*args, **kwargs)

    df = pd.DataFrame({"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    with tmpdir() as dn:
        ddf.to_json(dn, compute_kwargs={"scheduler": my_get})
        assert flag[0]
        result = dd.read_json(os.path.join(dn, "*"))
        assert_eq(result, df, check_index=False)


def test_read_json_error():
    with tmpfile("json") as f:
        with pytest.raises(ValueError):
            df.to_json(f, orient="split", lines=True)
        df.to_json(f, orient="split", lines=False)
        with pytest.raises(ValueError):
            dd.read_json(f, orient="split", blocksize=1)


@pytest.mark.parametrize("block", [5, 15, 33, 200, 90000])
def test_read_chunked(block):
    with tmpdir() as path:
        fn = os.path.join(path, "1.json")
        df.to_json(fn, orient="records", lines=True)
        d = dd.read_json(fn, blocksize=block, sample=10)
        assert (d.npartitions > 1) or (block > 30)
        assert_eq(d, df, check_index=False)


@pytest.mark.parametrize("compression", [None, "gzip", "xz"])
def test_json_compressed(compression):
    with tmpdir() as path:
        dd.to_json(ddf, path, compression=compression)
        actual = dd.read_json(os.path.join(path, "*"), compression=compression)
        assert_eq(df, actual, check_index=False)


def test_read_json_inferred_compression():
    with tmpdir() as path:
        fn = os.path.join(path, "*.json.gz")
        dd.to_json(ddf, fn, compression="gzip")
        actual = dd.read_json(fn)
        assert_eq(df, actual, check_index=False)


@pytest.mark.skipif(
    Version(fsspec.__version__) == Version("2023.9.1"),
    reason="https://github.com/dask/dask/issues/10515",
)
def test_to_json_results():
    with tmpfile("json") as f:
        paths = ddf.to_json(f)
        assert paths == [os.path.join(f, f"{n}.part") for n in range(ddf.npartitions)]

    with tmpfile("json") as f:
        list_of_delayed = ddf.to_json(f, compute=False)
        paths = dask.compute(*list_of_delayed)
        # this is a tuple rather than a list since it's the output of dask.compute
        assert paths == tuple(
            os.path.join(f, f"{n}.part") for n in range(ddf.npartitions)
        )
