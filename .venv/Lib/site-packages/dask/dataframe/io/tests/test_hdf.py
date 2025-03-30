from __future__ import annotations

import os
import pathlib
from functools import lru_cache
from time import sleep

import numpy as np
import pandas as pd
import pytest
from packaging.version import Version

import dask
import dask.dataframe as dd
from dask._compatibility import PY_VERSION
from dask.core import get_deps
from dask.dataframe._compat import tm
from dask.dataframe.utils import assert_eq
from dask.utils import tmpdir, tmpfile

# there's no support in upstream for writing HDF with extension dtypes yet.
# see https://github.com/pandas-dev/pandas/issues/31199
pytestmark = pytest.mark.skip_with_pyarrow_strings  # no support for hdf yet


def test_to_hdf():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    a = dd.from_pandas(df, 2)

    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data")
        out = pd.read_hdf(fn, "/data")
        tm.assert_frame_equal(df, out[:])

    with tmpfile("h5") as fn:
        a.x.to_hdf(fn, "/data")
        out = pd.read_hdf(fn, "/data")
        tm.assert_series_equal(df.x, out[:])

    a = dd.from_pandas(df, 1)
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data")
        out = pd.read_hdf(fn, "/data")
        tm.assert_frame_equal(df, out[:])

    # test compute = False
    with tmpfile("h5") as fn:
        r = a.to_hdf(fn, "/data", compute=False)
        r.compute()
        out = pd.read_hdf(fn, "/data")
        tm.assert_frame_equal(df, out[:])


@pytest.mark.skipif(
    PY_VERSION >= Version("3.11"),
    reason="segfaults due to https://github.com/PyTables/PyTables/issues/977",
)
def test_to_hdf_multiple_nodes():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    a = dd.from_pandas(df, 2)
    df16 = pd.DataFrame(
        {
            "x": [
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
                "h",
                "i",
                "j",
                "k",
                "l",
                "m",
                "n",
                "o",
                "p",
            ],
            "y": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        },
        index=[
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ],
    )
    b = dd.from_pandas(df16, 16)

    # saving to multiple nodes
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data*")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df, out)

    # saving to multiple nodes making sure order is kept
    with tmpfile("h5") as fn:
        b.to_hdf(fn, "/data*")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df16, out)

    # saving to multiple datasets with custom name_function
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data_*", name_function=lambda i: "a" * (i + 1))
        out = dd.read_hdf(fn, "/data_*")
        assert_eq(df, out)

        out = pd.read_hdf(fn, "/data_a")
        tm.assert_frame_equal(out, df.iloc[:2])
        out = pd.read_hdf(fn, "/data_aa")
        tm.assert_frame_equal(out, df.iloc[2:])

    # test multiple nodes with hdf object
    with tmpfile("h5") as fn:
        with pd.HDFStore(fn) as hdf:
            b.to_hdf(hdf, "/data*")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df16, out)


def test_to_hdf_multiple_files():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    a = dd.from_pandas(df, 2)
    df16 = pd.DataFrame(
        {
            "x": [
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
                "h",
                "i",
                "j",
                "k",
                "l",
                "m",
                "n",
                "o",
                "p",
            ],
            "y": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        },
        index=[
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ],
    )
    b = dd.from_pandas(df16, 16)

    # saving to multiple files
    with tmpdir() as dn:
        fn = os.path.join(dn, "data_*.h5")
        a.to_hdf(fn, "/data")
        out = dd.read_hdf(fn, "/data")
        assert_eq(df, out)

    # saving to multiple files making sure order is kept
    with tmpdir() as dn:
        fn = os.path.join(dn, "data_*.h5")
        b.to_hdf(fn, "/data")
        out = dd.read_hdf(fn, "/data")
        assert_eq(df16, out)

    # saving to multiple files where first file is longer
    # https://github.com/dask/dask/issues/8023
    with tmpdir() as dn:
        fn1 = os.path.join(dn, "data_1.h5")
        fn2 = os.path.join(dn, "data_2.h5")
        b.to_hdf(fn1, "/data")
        a.to_hdf(fn2, "/data")
        out = dd.read_hdf([fn1, fn2], "/data")
        assert_eq(pd.concat([df16, df]), out)

    # saving to multiple files with custom name_function
    with tmpdir() as dn:
        fn = os.path.join(dn, "data_*.h5")
        a.to_hdf(fn, "/data", name_function=lambda i: "a" * (i + 1))
        out = dd.read_hdf(fn, "/data")
        assert_eq(df, out)

        out = pd.read_hdf(os.path.join(dn, "data_a.h5"), "/data")
        tm.assert_frame_equal(out, df.iloc[:2])
        out = pd.read_hdf(os.path.join(dn, "data_aa.h5"), "/data")
        tm.assert_frame_equal(out, df.iloc[2:])

    # test hdf object
    with tmpfile("h5") as fn:
        with pd.HDFStore(fn) as hdf:
            a.to_hdf(hdf, "/data*")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df, out)


def test_to_hdf_modes_multiple_nodes():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )

    # appending a single partition to existing data
    a = dd.from_pandas(df, 1)
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data2")
        a.to_hdf(fn, "/data*", mode="a")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(dd.concat([df, df]), out)

    # overwriting a file with a single partition
    a = dd.from_pandas(df, 1)
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data2")
        a.to_hdf(fn, "/data*", mode="w")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df, out)

    # appending two partitions to existing data
    a = dd.from_pandas(df, 2)
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data2")
        a.to_hdf(fn, "/data*", mode="a")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(dd.concat([df, df]), out)

    # overwriting a file with two partitions
    a = dd.from_pandas(df, 2)
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data2")
        a.to_hdf(fn, "/data*", mode="w")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df, out)

    # overwriting a single partition, keeping other partitions
    a = dd.from_pandas(df, 2)
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data1")
        a.to_hdf(fn, "/data2")
        a.to_hdf(fn, "/data*", mode="a", append=False)
        out = dd.read_hdf(fn, "/data*")
        assert_eq(dd.concat([df, df]), out)


def test_to_hdf_modes_multiple_files():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )

    # appending a single partition to existing data
    a = dd.from_pandas(df, 1)
    with tmpdir() as dn:
        fn = os.path.join(dn, "data*")
        a.to_hdf(os.path.join(dn, "data2"), "/data")
        a.to_hdf(fn, "/data", mode="a")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(dd.concat([df, df]), out)

    # appending two partitions to existing data
    a = dd.from_pandas(df, 2)
    with tmpdir() as dn:
        fn = os.path.join(dn, "data*")
        a.to_hdf(os.path.join(dn, "data2"), "/data")
        a.to_hdf(fn, "/data", mode="a")
        out = dd.read_hdf(fn, "/data")
        assert_eq(dd.concat([df, df]), out)

    # overwriting a file with two partitions
    a = dd.from_pandas(df, 2)
    with tmpdir() as dn:
        fn = os.path.join(dn, "data*")
        a.to_hdf(os.path.join(dn, "data1"), "/data")
        a.to_hdf(fn, "/data", mode="w")
        out = dd.read_hdf(fn, "/data")
        assert_eq(df, out)

    # overwriting a single partition, keeping other partitions
    a = dd.from_pandas(df, 2)
    with tmpdir() as dn:
        fn = os.path.join(dn, "data*")
        a.to_hdf(os.path.join(dn, "data1"), "/data")
        a.to_hdf(fn, "/data", mode="a", append=False)
        out = dd.read_hdf(fn, "/data")
        assert_eq(dd.concat([df, df]), out)


def dependency_depth(dsk):
    deps, _ = get_deps(dsk)

    @lru_cache(maxsize=None)
    def max_depth_by_deps(key):
        if not deps[key]:
            return 1

        d = 1 + max(max_depth_by_deps(dep_key) for dep_key in deps[key])
        return d

    return max(max_depth_by_deps(dep_key) for dep_key in deps.keys())


@pytest.mark.skipif(
    PY_VERSION >= Version("3.11"),
    reason="segfaults due to https://github.com/PyTables/PyTables/issues/977",
)
@pytest.mark.slow
def test_to_hdf_lock_delays():
    pytest.importorskip("tables")
    df16 = pd.DataFrame(
        {
            "x": [
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
                "h",
                "i",
                "j",
                "k",
                "l",
                "m",
                "n",
                "o",
                "p",
            ],
            "y": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        },
        index=[
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ],
    )
    a = dd.from_pandas(df16, 16)

    # adding artificial delays to make sure last tasks finish first
    # that's a way to simulate last tasks finishing last
    def delayed_nop(i):
        if i.iloc[1] < 10:
            sleep(0.1 * (10 - i.iloc[1]))
        return i

    # saving to multiple hdf nodes
    with tmpfile() as fn:
        a = a.apply(delayed_nop, axis=1, meta=a)
        a.to_hdf(fn, "/data*")
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df16, out)

    # saving to multiple hdf files
    # adding artificial delays to make sure last tasks finish first
    with tmpdir() as dn:
        fn = os.path.join(dn, "data*")
        a = a.apply(delayed_nop, axis=1, meta=a)
        a.to_hdf(fn, "/data")
        out = dd.read_hdf(fn, "/data")
        assert_eq(df16, out)


def test_to_hdf_exceptions():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    a = dd.from_pandas(df, 1)

    # triggering too many asterisks error
    with tmpdir() as dn:
        with pytest.raises(ValueError):
            fn = os.path.join(dn, "data_*.h5")
            a.to_hdf(fn, "/data_*")

    # triggering too many asterisks error
    with tmpfile() as fn:
        with pd.HDFStore(fn) as hdf:
            with pytest.raises(ValueError):
                a.to_hdf(hdf, "/data_*_*")


@pytest.mark.skipif(
    PY_VERSION >= Version("3.11"),
    reason="segfaults due to https://github.com/PyTables/PyTables/issues/977",
)
@pytest.mark.parametrize("scheduler", ["sync", "threads", "processes"])
@pytest.mark.parametrize("npartitions", [1, 4, 10])
def test_to_hdf_schedulers(scheduler, npartitions):
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {
            "x": [
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
                "h",
                "i",
                "j",
                "k",
                "l",
                "m",
                "n",
                "o",
                "p",
            ],
            "y": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        },
        index=[
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ],
    )
    a = dd.from_pandas(df, npartitions=npartitions)

    # test single file single node
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data", scheduler=scheduler)
        out = pd.read_hdf(fn, "/data")
        assert_eq(df, out)

    # test multiple files single node
    with tmpdir() as dn:
        fn = os.path.join(dn, "data_*.h5")
        a.to_hdf(fn, "/data", scheduler=scheduler)
        out = dd.read_hdf(fn, "/data")
        assert_eq(df, out)

    # test single file multiple nodes
    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data*", scheduler=scheduler)
        out = dd.read_hdf(fn, "/data*")
        assert_eq(df, out)


def test_to_hdf_kwargs():
    pytest.importorskip("tables")
    df = pd.DataFrame({"A": ["a", "aaaa"]})
    ddf = dd.from_pandas(df, npartitions=2)
    with tmpfile("h5") as fn:
        ddf.to_hdf(fn, "foo4", format="table", min_itemsize=4)
        df2 = pd.read_hdf(fn, "foo4")
        tm.assert_frame_equal(df, df2)

    # test shorthand 't' for table
    with tmpfile("h5") as fn:
        ddf.to_hdf(fn, "foo4", format="t", min_itemsize=4)
        df2 = pd.read_hdf(fn, "foo4")
        tm.assert_frame_equal(df, df2)


def test_to_fmt_warns():
    pytest.importorskip("tables")
    df16 = pd.DataFrame(
        {
            "x": [
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
                "h",
                "i",
                "j",
                "k",
                "l",
                "m",
                "n",
                "o",
                "p",
            ],
            "y": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        },
        index=[
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ],
    )
    a = dd.from_pandas(df16, 16)

    # testing warning when breaking order
    with tmpfile("h5") as fn:
        with pytest.warns(
            UserWarning, match="To preserve order between partitions name_function"
        ):
            a.to_hdf(fn, "/data*", name_function=str)

    with tmpdir() as dn:
        fn = os.path.join(dn, "data_*.csv")
        a.to_csv(fn, name_function=str)


@pytest.mark.parametrize(
    "data, compare",
    [
        (
            pd.DataFrame(
                {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]},
                index=[1.0, 2.0, 3.0, 4.0],
            ),
            tm.assert_frame_equal,
        ),
        (pd.Series([1, 2, 3, 4], name="a"), tm.assert_series_equal),
    ],
)
def test_read_hdf(data, compare):
    pytest.importorskip("tables")
    with tmpfile("h5") as fn:
        data.to_hdf(fn, key="/data")
        try:
            dd.read_hdf(fn, "data", chunksize=2, mode="r")
            assert False
        except TypeError as e:
            assert "format='table'" in str(e)

    with tmpfile("h5") as fn:
        data.to_hdf(fn, key="/data", format="table")
        a = dd.read_hdf(fn, "/data", chunksize=2, mode="r")
        assert a.npartitions == 2

        compare(a.compute(), data)

        compare(
            dd.read_hdf(fn, "/data", chunksize=2, start=1, stop=3, mode="r").compute(),
            pd.read_hdf(fn, "/data", start=1, stop=3),
        )

    with tmpfile("h5") as fn:
        sorted_data = data.sort_index()
        sorted_data.to_hdf(fn, key="/data", format="table")
        a = dd.read_hdf(fn, "/data", chunksize=2, sorted_index=True, mode="r")
        assert a.npartitions == 2

        compare(a.compute(), sorted_data)


def test_read_hdf_multiply_open():
    """Test that we can read from a file that's already opened elsewhere in
    read-only mode."""
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    with tmpfile("h5") as fn:
        df.to_hdf(fn, key="/data", format="table")
        with pd.HDFStore(fn, mode="r"):
            dd.read_hdf(fn, "/data", chunksize=2, mode="r")


@pytest.mark.skipif(
    PY_VERSION >= Version("3.11"),
    reason="segfaults due to https://github.com/PyTables/PyTables/issues/977",
)
def test_read_hdf_multiple():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {
            "x": [
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
                "h",
                "i",
                "j",
                "k",
                "l",
                "m",
                "n",
                "o",
                "p",
            ],
            "y": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        },
        index=[
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ],
    )
    a = dd.from_pandas(df, 16)

    with tmpfile("h5") as fn:
        a.to_hdf(fn, "/data*")
        r = dd.read_hdf(fn, "/data*", sorted_index=True)
        assert a.npartitions == r.npartitions
        assert a.divisions == r.divisions
        assert_eq(a, r)


def test_read_hdf_start_stop_values():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    with tmpfile("h5") as fn:
        df.to_hdf(fn, key="/data", format="table")

        with pytest.raises(ValueError, match="number of rows"):
            dd.read_hdf(fn, "/data", stop=10)

        with pytest.raises(ValueError, match="is above or equal to"):
            dd.read_hdf(fn, "/data", start=10)

        with pytest.raises(ValueError, match="positive integer"):
            dd.read_hdf(fn, "/data", chunksize=-1)


def test_hdf_globbing():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )

    with tmpdir() as tdir:
        df.to_hdf(os.path.join(tdir, "one.h5"), key="/foo/data", format="table")
        df.to_hdf(os.path.join(tdir, "two.h5"), key="/bar/data", format="table")
        df.to_hdf(os.path.join(tdir, "two.h5"), key="/foo/data", format="table")

        with dask.config.set(scheduler="sync"):
            res = dd.read_hdf(os.path.join(tdir, "one.h5"), "/*/data", chunksize=2)
            assert res.npartitions == 2
            tm.assert_frame_equal(res.compute(), df)

            res = dd.read_hdf(
                os.path.join(tdir, "one.h5"), "/*/data", chunksize=2, start=1, stop=3
            )
            expected = pd.read_hdf(
                os.path.join(tdir, "one.h5"), "/foo/data", start=1, stop=3
            )
            tm.assert_frame_equal(res.compute(), expected)

            res = dd.read_hdf(os.path.join(tdir, "two.h5"), "/*/data", chunksize=2)
            assert res.npartitions == 2 + 2
            tm.assert_frame_equal(res.compute(), pd.concat([df] * 2))

            res = dd.read_hdf(os.path.join(tdir, "*.h5"), "/foo/data", chunksize=2)
            assert res.npartitions == 2 + 2
            tm.assert_frame_equal(res.compute(), pd.concat([df] * 2))

            res = dd.read_hdf(os.path.join(tdir, "*.h5"), "/*/data", chunksize=2)
            assert res.npartitions == 2 + 2 + 2
            tm.assert_frame_equal(res.compute(), pd.concat([df] * 3))


def test_hdf_file_list():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )

    with tmpdir() as tdir:
        df.iloc[:2].to_hdf(
            os.path.join(tdir, "one.h5"), key="dataframe", format="table"
        )
        df.iloc[2:].to_hdf(
            os.path.join(tdir, "two.h5"), key="dataframe", format="table"
        )

        with dask.config.set(scheduler="sync"):
            input_files = [os.path.join(tdir, "one.h5"), os.path.join(tdir, "two.h5")]
            res = dd.read_hdf(input_files, "dataframe")
            tm.assert_frame_equal(res.compute(), df)


def test_read_hdf_pattern_pathlike():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )

    with tmpfile("h5") as fn:
        path = pathlib.Path(fn)
        df.to_hdf(path, key="dataframe", format="table")
        res = dd.read_hdf(path, "dataframe")
        assert_eq(res, df)


def test_to_hdf_path_pathlike():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    ddf = dd.from_pandas(df, npartitions=3)

    with tmpfile("h5") as fn:
        path = pathlib.Path(fn)
        ddf.to_hdf(path, "/data")
        res = pd.read_hdf(path, "/data")
        assert_eq(res, ddf)


def test_read_hdf_doesnt_segfault():
    pytest.importorskip("tables")
    with tmpfile("h5") as fn:
        N = 40
        df = pd.DataFrame(np.random.randn(N, 3))
        with pd.HDFStore(fn, mode="w") as store:
            store.append("/x", df)

        ddf = dd.read_hdf(fn, "/x", chunksize=2)
        assert len(ddf) == N


def test_hdf_filenames():
    pytest.importorskip("tables")
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [1, 2, 3, 4]}, index=[1.0, 2.0, 3.0, 4.0]
    )
    ddf = dd.from_pandas(df, npartitions=2)
    assert ddf.to_hdf("foo*.hdf5", "key") == ["foo0.hdf5", "foo1.hdf5"]
    os.remove("foo0.hdf5")
    os.remove("foo1.hdf5")


def test_hdf_path_exceptions():
    # single file doesn't exist
    with pytest.raises(IOError):
        dd.read_hdf("nonexistant_store_X34HJK", "/tmp")

    # a file from a list of files doesn't exist
    with pytest.raises(IOError):
        dd.read_hdf(["nonexistant_store_X34HJK", "nonexistant_store_UY56YH"], "/tmp")

    # list of files is empty
    with pytest.raises(ValueError):
        dd.read_hdf([], "/tmp")


def test_hdf_nonpandas_keys():
    # https://github.com/dask/dask/issues/5934
    # TODO: maybe remove this if/when pandas copes with all keys

    tables = pytest.importorskip("tables")
    import tables

    class Table1(tables.IsDescription):
        value1 = tables.Float32Col()

    class Table2(tables.IsDescription):
        value2 = tables.Float32Col()

    class Table3(tables.IsDescription):
        value3 = tables.Float32Col()

    with tmpfile("h5") as path:
        with tables.open_file(path, mode="a") as h5file:
            group = h5file.create_group("/", "group")
            t = h5file.create_table(group, "table1", Table1, "Table 1")
            row = t.row
            row["value1"] = 1
            row.append()
            t = h5file.create_table(group, "table2", Table2, "Table 2")
            row = t.row
            row["value2"] = 1
            row.append()
            t = h5file.create_table(group, "table3", Table3, "Table 3")
            row = t.row
            row["value3"] = 1
            row.append()

        # pandas keys should still work
        bar = pd.DataFrame(np.random.randn(10, 4))
        bar.to_hdf(path, key="/bar", format="table", mode="a")

        dd.read_hdf(path, "/group/table1")
        dd.read_hdf(path, "/group/table2")
        dd.read_hdf(path, "/group/table3")
        dd.read_hdf(path, "/bar")


def test_hdf_empty_dataframe(tmp_path):
    pytest.importorskip("tables")
    # https://github.com/dask/dask/issues/8707
    from dask.dataframe.io.hdf import dont_use_fixed_error_message

    df = pd.DataFrame({"A": [], "B": []}, index=[])
    df.to_hdf(tmp_path / "data.h5", format="fixed", key="df", mode="w")
    with pytest.raises(TypeError, match=dont_use_fixed_error_message):
        dd.read_hdf(tmp_path / "data.h5", "df")
