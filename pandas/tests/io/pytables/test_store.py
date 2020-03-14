import datetime
from datetime import timedelta
from distutils.version import LooseVersion
from io import BytesIO
import os
from pathlib import Path
import re
from warnings import catch_warnings, simplefilter

import numpy as np
import pytest

from pandas.compat import is_platform_little_endian, is_platform_windows
import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    Categorical,
    CategoricalIndex,
    DataFrame,
    DatetimeIndex,
    Index,
    Int64Index,
    MultiIndex,
    RangeIndex,
    Series,
    Timestamp,
    bdate_range,
    concat,
    date_range,
    isna,
    timedelta_range,
)
import pandas._testing as tm
from pandas.tests.io.pytables.common import (
    _maybe_remove,
    create_tempfile,
    ensure_clean_path,
    ensure_clean_store,
    safe_close,
    safe_remove,
    tables,
)

from pandas.io.pytables import (
    ClosedFileError,
    HDFStore,
    PossibleDataLossError,
    Term,
    read_hdf,
)

from pandas.io import pytables as pytables  # noqa: E402 isort:skip
from pandas.io.pytables import TableIterator  # noqa: E402 isort:skip


_default_compressor = "blosc"
ignore_natural_naming_warning = pytest.mark.filterwarnings(
    "ignore:object name:tables.exceptions.NaturalNameWarning"
)


@pytest.mark.single
class TestHDFStore:
    def test_format_type(self, setup_path):
        df = pd.DataFrame({"A": [1, 2]})
        with ensure_clean_path(setup_path) as path:
            with HDFStore(path) as store:
                store.put("a", df, format="fixed")
                store.put("b", df, format="table")

                assert store.get_storer("a").format_type == "fixed"
                assert store.get_storer("b").format_type == "table"

    def test_format_kwarg_in_constructor(self, setup_path):
        # GH 13291

        msg = "format is not a defined argument for HDFStore"

        with ensure_clean_path(setup_path) as path:
            with pytest.raises(ValueError, match=msg):
                HDFStore(path, format="table")

    def test_context(self, setup_path):
        path = create_tempfile(setup_path)
        try:
            with HDFStore(path) as tbl:
                raise ValueError("blah")
        except ValueError:
            pass
        finally:
            safe_remove(path)

        try:
            with HDFStore(path) as tbl:
                tbl["a"] = tm.makeDataFrame()

            with HDFStore(path) as tbl:
                assert len(tbl) == 1
                assert type(tbl["a"]) == DataFrame
        finally:
            safe_remove(path)

    def test_conv_read_write(self, setup_path):
        path = create_tempfile(setup_path)
        try:

            def roundtrip(key, obj, **kwargs):
                obj.to_hdf(path, key, **kwargs)
                return read_hdf(path, key)

            o = tm.makeTimeSeries()
            tm.assert_series_equal(o, roundtrip("series", o))

            o = tm.makeStringSeries()
            tm.assert_series_equal(o, roundtrip("string_series", o))

            o = tm.makeDataFrame()
            tm.assert_frame_equal(o, roundtrip("frame", o))

            # table
            df = DataFrame(dict(A=range(5), B=range(5)))
            df.to_hdf(path, "table", append=True)
            result = read_hdf(path, "table", where=["index>2"])
            tm.assert_frame_equal(df[df.index > 2], result)

        finally:
            safe_remove(path)

    def test_long_strings(self, setup_path):

        # GH6166
        df = DataFrame(
            {"a": tm.rands_array(100, size=10)}, index=tm.rands_array(100, size=10)
        )

        with ensure_clean_store(setup_path) as store:
            store.append("df", df, data_columns=["a"])

            result = store.select("df")
            tm.assert_frame_equal(df, result)

    def test_api(self, setup_path):

        # GH4584
        # API issue when to_hdf doesn't accept append AND format args
        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()
            df.iloc[:10].to_hdf(path, "df", append=True, format="table")
            df.iloc[10:].to_hdf(path, "df", append=True, format="table")
            tm.assert_frame_equal(read_hdf(path, "df"), df)

            # append to False
            df.iloc[:10].to_hdf(path, "df", append=False, format="table")
            df.iloc[10:].to_hdf(path, "df", append=True, format="table")
            tm.assert_frame_equal(read_hdf(path, "df"), df)

        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()
            df.iloc[:10].to_hdf(path, "df", append=True)
            df.iloc[10:].to_hdf(path, "df", append=True, format="table")
            tm.assert_frame_equal(read_hdf(path, "df"), df)

            # append to False
            df.iloc[:10].to_hdf(path, "df", append=False, format="table")
            df.iloc[10:].to_hdf(path, "df", append=True)
            tm.assert_frame_equal(read_hdf(path, "df"), df)

        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()
            df.to_hdf(path, "df", append=False, format="fixed")
            tm.assert_frame_equal(read_hdf(path, "df"), df)

            df.to_hdf(path, "df", append=False, format="f")
            tm.assert_frame_equal(read_hdf(path, "df"), df)

            df.to_hdf(path, "df", append=False)
            tm.assert_frame_equal(read_hdf(path, "df"), df)

            df.to_hdf(path, "df")
            tm.assert_frame_equal(read_hdf(path, "df"), df)

        with ensure_clean_store(setup_path) as store:

            path = store._path
            df = tm.makeDataFrame()

            _maybe_remove(store, "df")
            store.append("df", df.iloc[:10], append=True, format="table")
            store.append("df", df.iloc[10:], append=True, format="table")
            tm.assert_frame_equal(store.select("df"), df)

            # append to False
            _maybe_remove(store, "df")
            store.append("df", df.iloc[:10], append=False, format="table")
            store.append("df", df.iloc[10:], append=True, format="table")
            tm.assert_frame_equal(store.select("df"), df)

            # formats
            _maybe_remove(store, "df")
            store.append("df", df.iloc[:10], append=False, format="table")
            store.append("df", df.iloc[10:], append=True, format="table")
            tm.assert_frame_equal(store.select("df"), df)

            _maybe_remove(store, "df")
            store.append("df", df.iloc[:10], append=False, format="table")
            store.append("df", df.iloc[10:], append=True, format=None)
            tm.assert_frame_equal(store.select("df"), df)

        with ensure_clean_path(setup_path) as path:
            # Invalid.
            df = tm.makeDataFrame()

            msg = "Can only append to Tables"

            with pytest.raises(ValueError, match=msg):
                df.to_hdf(path, "df", append=True, format="f")

            with pytest.raises(ValueError, match=msg):
                df.to_hdf(path, "df", append=True, format="fixed")

            msg = r"invalid HDFStore format specified \[foo\]"

            with pytest.raises(TypeError, match=msg):
                df.to_hdf(path, "df", append=True, format="foo")

            with pytest.raises(TypeError, match=msg):
                df.to_hdf(path, "df", append=False, format="foo")

        # File path doesn't exist
        path = ""
        msg = f"File {path} does not exist"

        with pytest.raises(FileNotFoundError, match=msg):
            read_hdf(path, "df")

    def test_api_default_format(self, setup_path):

        # default_format option
        with ensure_clean_store(setup_path) as store:
            df = tm.makeDataFrame()

            pd.set_option("io.hdf.default_format", "fixed")
            _maybe_remove(store, "df")
            store.put("df", df)
            assert not store.get_storer("df").is_table

            msg = "Can only append to Tables"

            with pytest.raises(ValueError, match=msg):
                store.append("df2", df)

            pd.set_option("io.hdf.default_format", "table")
            _maybe_remove(store, "df")
            store.put("df", df)
            assert store.get_storer("df").is_table
            _maybe_remove(store, "df2")
            store.append("df2", df)
            assert store.get_storer("df").is_table

            pd.set_option("io.hdf.default_format", None)

        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()

            pd.set_option("io.hdf.default_format", "fixed")
            df.to_hdf(path, "df")
            with HDFStore(path) as store:
                assert not store.get_storer("df").is_table
            with pytest.raises(ValueError, match=msg):
                df.to_hdf(path, "df2", append=True)

            pd.set_option("io.hdf.default_format", "table")
            df.to_hdf(path, "df3")
            with HDFStore(path) as store:
                assert store.get_storer("df3").is_table
            df.to_hdf(path, "df4", append=True)
            with HDFStore(path) as store:
                assert store.get_storer("df4").is_table

            pd.set_option("io.hdf.default_format", None)

    def test_keys(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            store["a"] = tm.makeTimeSeries()
            store["b"] = tm.makeStringSeries()
            store["c"] = tm.makeDataFrame()

            assert len(store) == 3
            expected = {"/a", "/b", "/c"}
            assert set(store.keys()) == expected
            assert set(store) == expected

    def test_keys_ignore_hdf_softlink(self, setup_path):

        # GH 20523
        # Puts a softlink into HDF file and rereads

        with ensure_clean_store(setup_path) as store:

            df = DataFrame(dict(A=range(5), B=range(5)))
            store.put("df", df)

            assert store.keys() == ["/df"]

            store._handle.create_soft_link(store._handle.root, "symlink", "df")

            # Should ignore the softlink
            assert store.keys() == ["/df"]

    def test_iter_empty(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            # GH 12221
            assert list(store) == []

    def test_repr(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            repr(store)
            store.info()
            store["a"] = tm.makeTimeSeries()
            store["b"] = tm.makeStringSeries()
            store["c"] = tm.makeDataFrame()

            df = tm.makeDataFrame()
            df["obj1"] = "foo"
            df["obj2"] = "bar"
            df["bool1"] = df["A"] > 0
            df["bool2"] = df["B"] > 0
            df["bool3"] = True
            df["int1"] = 1
            df["int2"] = 2
            df["timestamp1"] = Timestamp("20010102")
            df["timestamp2"] = Timestamp("20010103")
            df["datetime1"] = datetime.datetime(2001, 1, 2, 0, 0)
            df["datetime2"] = datetime.datetime(2001, 1, 3, 0, 0)
            df.loc[df.index[3:6], ["obj1"]] = np.nan
            df = df._consolidate()._convert(datetime=True)

            with catch_warnings(record=True):
                simplefilter("ignore", pd.errors.PerformanceWarning)
                store["df"] = df

            # make a random group in hdf space
            store._handle.create_group(store._handle.root, "bah")

            assert store.filename in repr(store)
            assert store.filename in str(store)
            store.info()

        # storers
        with ensure_clean_store(setup_path) as store:

            df = tm.makeDataFrame()
            store.append("df", df)

            s = store.get_storer("df")
            repr(s)
            str(s)

    @ignore_natural_naming_warning
    def test_contains(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            store["a"] = tm.makeTimeSeries()
            store["b"] = tm.makeDataFrame()
            store["foo/bar"] = tm.makeDataFrame()
            assert "a" in store
            assert "b" in store
            assert "c" not in store
            assert "foo/bar" in store
            assert "/foo/bar" in store
            assert "/foo/b" not in store
            assert "bar" not in store

            # gh-2694: tables.NaturalNameWarning
            with catch_warnings(record=True):
                store["node())"] = tm.makeDataFrame()
            assert "node())" in store

    def test_versioning(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            store["a"] = tm.makeTimeSeries()
            store["b"] = tm.makeDataFrame()
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, "df1")
            store.append("df1", df[:10])
            store.append("df1", df[10:])
            assert store.root.a._v_attrs.pandas_version == "0.15.2"
            assert store.root.b._v_attrs.pandas_version == "0.15.2"
            assert store.root.df1._v_attrs.pandas_version == "0.15.2"

            # write a file and wipe its versioning
            _maybe_remove(store, "df2")
            store.append("df2", df)

            # this is an error because its table_type is appendable, but no
            # version info
            store.get_node("df2")._v_attrs.pandas_version = None

            msg = "'NoneType' object has no attribute 'startswith'"

            with pytest.raises(Exception, match=msg):
                store.select("df2")

    def test_mode(self, setup_path):

        df = tm.makeTimeDataFrame()

        def check(mode):

            with ensure_clean_path(setup_path) as path:

                # constructor
                if mode in ["r", "r+"]:
                    with pytest.raises(IOError):
                        HDFStore(path, mode=mode)

                else:
                    store = HDFStore(path, mode=mode)
                    assert store._handle.mode == mode
                    store.close()

            with ensure_clean_path(setup_path) as path:

                # context
                if mode in ["r", "r+"]:
                    with pytest.raises(IOError):
                        with HDFStore(path, mode=mode) as store:  # noqa
                            pass
                else:
                    with HDFStore(path, mode=mode) as store:
                        assert store._handle.mode == mode

            with ensure_clean_path(setup_path) as path:

                # conv write
                if mode in ["r", "r+"]:
                    with pytest.raises(IOError):
                        df.to_hdf(path, "df", mode=mode)
                    df.to_hdf(path, "df", mode="w")
                else:
                    df.to_hdf(path, "df", mode=mode)

                # conv read
                if mode in ["w"]:
                    msg = (
                        "mode w is not allowed while performing a read. "
                        r"Allowed modes are r, r\+ and a."
                    )
                    with pytest.raises(ValueError, match=msg):
                        read_hdf(path, "df", mode=mode)
                else:
                    result = read_hdf(path, "df", mode=mode)
                    tm.assert_frame_equal(result, df)

        def check_default_mode():

            # read_hdf uses default mode
            with ensure_clean_path(setup_path) as path:
                df.to_hdf(path, "df", mode="w")
                result = read_hdf(path, "df")
                tm.assert_frame_equal(result, df)

        check("r")
        check("r+")
        check("a")
        check("w")
        check_default_mode()

    def test_reopen_handle(self, setup_path):

        with ensure_clean_path(setup_path) as path:

            store = HDFStore(path, mode="a")
            store["a"] = tm.makeTimeSeries()

            # invalid mode change
            with pytest.raises(PossibleDataLossError):
                store.open("w")

            store.close()
            assert not store.is_open

            # truncation ok here
            store.open("w")
            assert store.is_open
            assert len(store) == 0
            store.close()
            assert not store.is_open

            store = HDFStore(path, mode="a")
            store["a"] = tm.makeTimeSeries()

            # reopen as read
            store.open("r")
            assert store.is_open
            assert len(store) == 1
            assert store._mode == "r"
            store.close()
            assert not store.is_open

            # reopen as append
            store.open("a")
            assert store.is_open
            assert len(store) == 1
            assert store._mode == "a"
            store.close()
            assert not store.is_open

            # reopen as append (again)
            store.open("a")
            assert store.is_open
            assert len(store) == 1
            assert store._mode == "a"
            store.close()
            assert not store.is_open

    def test_open_args(self, setup_path):

        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()

            # create an in memory store
            store = HDFStore(
                path, mode="a", driver="H5FD_CORE", driver_core_backing_store=0
            )
            store["df"] = df
            store.append("df2", df)

            tm.assert_frame_equal(store["df"], df)
            tm.assert_frame_equal(store["df2"], df)

            store.close()

            # the file should not have actually been written
            assert not os.path.exists(path)

    def test_flush(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            store["a"] = tm.makeTimeSeries()
            store.flush()
            store.flush(fsync=True)

    def test_get(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            store["a"] = tm.makeTimeSeries()
            left = store.get("a")
            right = store["a"]
            tm.assert_series_equal(left, right)

            left = store.get("/a")
            right = store["/a"]
            tm.assert_series_equal(left, right)

            with pytest.raises(KeyError, match="'No object named b in the file'"):
                store.get("b")

    @pytest.mark.parametrize(
        "where, expected",
        [
            (
                "/",
                {
                    "": ({"first_group", "second_group"}, set()),
                    "/first_group": (set(), {"df1", "df2"}),
                    "/second_group": ({"third_group"}, {"df3", "s1"}),
                    "/second_group/third_group": (set(), {"df4"}),
                },
            ),
            (
                "/second_group",
                {
                    "/second_group": ({"third_group"}, {"df3", "s1"}),
                    "/second_group/third_group": (set(), {"df4"}),
                },
            ),
        ],
    )
    def test_walk(self, where, expected, setup_path):
        # GH10143
        objs = {
            "df1": pd.DataFrame([1, 2, 3]),
            "df2": pd.DataFrame([4, 5, 6]),
            "df3": pd.DataFrame([6, 7, 8]),
            "df4": pd.DataFrame([9, 10, 11]),
            "s1": pd.Series([10, 9, 8]),
            # Next 3 items aren't pandas objects and should be ignored
            "a1": np.array([[1, 2, 3], [4, 5, 6]]),
            "tb1": np.array([(1, 2, 3), (4, 5, 6)], dtype="i,i,i"),
            "tb2": np.array([(7, 8, 9), (10, 11, 12)], dtype="i,i,i"),
        }

        with ensure_clean_store("walk_groups.hdf", mode="w") as store:
            store.put("/first_group/df1", objs["df1"])
            store.put("/first_group/df2", objs["df2"])
            store.put("/second_group/df3", objs["df3"])
            store.put("/second_group/s1", objs["s1"])
            store.put("/second_group/third_group/df4", objs["df4"])
            # Create non-pandas objects
            store._handle.create_array("/first_group", "a1", objs["a1"])
            store._handle.create_table("/first_group", "tb1", obj=objs["tb1"])
            store._handle.create_table("/second_group", "tb2", obj=objs["tb2"])

            assert len(list(store.walk(where=where))) == len(expected)
            for path, groups, leaves in store.walk(where=where):
                assert path in expected
                expected_groups, expected_frames = expected[path]
                assert expected_groups == set(groups)
                assert expected_frames == set(leaves)
                for leaf in leaves:
                    frame_path = "/".join([path, leaf])
                    obj = store.get(frame_path)
                    if "df" in leaf:
                        tm.assert_frame_equal(obj, objs[leaf])
                    else:
                        tm.assert_series_equal(obj, objs[leaf])

    def test_getattr(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            s = tm.makeTimeSeries()
            store["a"] = s

            # test attribute access
            result = store.a
            tm.assert_series_equal(result, s)
            result = getattr(store, "a")
            tm.assert_series_equal(result, s)

            df = tm.makeTimeDataFrame()
            store["df"] = df
            result = store.df
            tm.assert_frame_equal(result, df)

            # errors
            for x in ["d", "mode", "path", "handle", "complib"]:
                with pytest.raises(AttributeError):
                    getattr(store, x)

            # not stores
            for x in ["mode", "path", "handle", "complib"]:
                getattr(store, f"_{x}")

    def test_put(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            ts = tm.makeTimeSeries()
            df = tm.makeTimeDataFrame()
            store["a"] = ts
            store["b"] = df[:10]
            store["foo/bar/bah"] = df[:10]
            store["foo"] = df[:10]
            store["/foo"] = df[:10]
            store.put("c", df[:10], format="table")

            # not OK, not a table
            with pytest.raises(ValueError):
                store.put("b", df[10:], append=True)

            # node does not currently exist, test _is_table_type returns False
            # in this case
            _maybe_remove(store, "f")
            with pytest.raises(ValueError):
                store.put("f", df[10:], append=True)

            # can't put to a table (use append instead)
            with pytest.raises(ValueError):
                store.put("c", df[10:], append=True)

            # overwrite table
            store.put("c", df[:10], format="table", append=False)
            tm.assert_frame_equal(df[:10], store["c"])

    def test_put_string_index(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            index = Index([f"I am a very long string index: {i}" for i in range(20)])
            s = Series(np.arange(20), index=index)
            df = DataFrame({"A": s, "B": s})

            store["a"] = s
            tm.assert_series_equal(store["a"], s)

            store["b"] = df
            tm.assert_frame_equal(store["b"], df)

            # mixed length
            index = Index(
                ["abcdefghijklmnopqrstuvwxyz1234567890"]
                + [f"I am a very long string index: {i}" for i in range(20)]
            )
            s = Series(np.arange(21), index=index)
            df = DataFrame({"A": s, "B": s})
            store["a"] = s
            tm.assert_series_equal(store["a"], s)

            store["b"] = df
            tm.assert_frame_equal(store["b"], df)

    def test_put_compression(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            df = tm.makeTimeDataFrame()

            store.put("c", df, format="table", complib="zlib")
            tm.assert_frame_equal(store["c"], df)

            # can't compress if format='fixed'
            with pytest.raises(ValueError):
                store.put("b", df, format="fixed", complib="zlib")

    @td.skip_if_windows_python_3
    def test_put_compression_blosc(self, setup_path):
        df = tm.makeTimeDataFrame()

        with ensure_clean_store(setup_path) as store:

            # can't compress if format='fixed'
            with pytest.raises(ValueError):
                store.put("b", df, format="fixed", complib="blosc")

            store.put("c", df, format="table", complib="blosc")
            tm.assert_frame_equal(store["c"], df)

    def test_complibs_default_settings(self, setup_path):
        # GH15943
        df = tm.makeDataFrame()

        # Set complevel and check if complib is automatically set to
        # default value
        with ensure_clean_path(setup_path) as tmpfile:
            df.to_hdf(tmpfile, "df", complevel=9)
            result = pd.read_hdf(tmpfile, "df")
            tm.assert_frame_equal(result, df)

            with tables.open_file(tmpfile, mode="r") as h5file:
                for node in h5file.walk_nodes(where="/df", classname="Leaf"):
                    assert node.filters.complevel == 9
                    assert node.filters.complib == "zlib"

        # Set complib and check to see if compression is disabled
        with ensure_clean_path(setup_path) as tmpfile:
            df.to_hdf(tmpfile, "df", complib="zlib")
            result = pd.read_hdf(tmpfile, "df")
            tm.assert_frame_equal(result, df)

            with tables.open_file(tmpfile, mode="r") as h5file:
                for node in h5file.walk_nodes(where="/df", classname="Leaf"):
                    assert node.filters.complevel == 0
                    assert node.filters.complib is None

        # Check if not setting complib or complevel results in no compression
        with ensure_clean_path(setup_path) as tmpfile:
            df.to_hdf(tmpfile, "df")
            result = pd.read_hdf(tmpfile, "df")
            tm.assert_frame_equal(result, df)

            with tables.open_file(tmpfile, mode="r") as h5file:
                for node in h5file.walk_nodes(where="/df", classname="Leaf"):
                    assert node.filters.complevel == 0
                    assert node.filters.complib is None

        # Check if file-defaults can be overridden on a per table basis
        with ensure_clean_path(setup_path) as tmpfile:
            store = pd.HDFStore(tmpfile)
            store.append("dfc", df, complevel=9, complib="blosc")
            store.append("df", df)
            store.close()

            with tables.open_file(tmpfile, mode="r") as h5file:
                for node in h5file.walk_nodes(where="/df", classname="Leaf"):
                    assert node.filters.complevel == 0
                    assert node.filters.complib is None
                for node in h5file.walk_nodes(where="/dfc", classname="Leaf"):
                    assert node.filters.complevel == 9
                    assert node.filters.complib == "blosc"

    def test_complibs(self, setup_path):
        # GH14478
        df = tm.makeDataFrame()

        # Building list of all complibs and complevels tuples
        all_complibs = tables.filters.all_complibs
        # Remove lzo if its not available on this platform
        if not tables.which_lib_version("lzo"):
            all_complibs.remove("lzo")
        # Remove bzip2 if its not available on this platform
        if not tables.which_lib_version("bzip2"):
            all_complibs.remove("bzip2")

        all_levels = range(0, 10)
        all_tests = [(lib, lvl) for lib in all_complibs for lvl in all_levels]

        for (lib, lvl) in all_tests:
            with ensure_clean_path(setup_path) as tmpfile:
                gname = "foo"

                # Write and read file to see if data is consistent
                df.to_hdf(tmpfile, gname, complib=lib, complevel=lvl)
                result = pd.read_hdf(tmpfile, gname)
                tm.assert_frame_equal(result, df)

                # Open file and check metadata
                # for correct amount of compression
                h5table = tables.open_file(tmpfile, mode="r")
                for node in h5table.walk_nodes(where="/" + gname, classname="Leaf"):
                    assert node.filters.complevel == lvl
                    if lvl == 0:
                        assert node.filters.complib is None
                    else:
                        assert node.filters.complib == lib
                h5table.close()

    def test_put_integer(self, setup_path):
        # non-date, non-string index
        df = DataFrame(np.random.randn(50, 100))
        self._check_roundtrip(df, tm.assert_frame_equal, setup_path)

    @td.xfail_non_writeable
    def test_put_mixed_type(self, setup_path):
        df = tm.makeTimeDataFrame()
        df["obj1"] = "foo"
        df["obj2"] = "bar"
        df["bool1"] = df["A"] > 0
        df["bool2"] = df["B"] > 0
        df["bool3"] = True
        df["int1"] = 1
        df["int2"] = 2
        df["timestamp1"] = Timestamp("20010102")
        df["timestamp2"] = Timestamp("20010103")
        df["datetime1"] = datetime.datetime(2001, 1, 2, 0, 0)
        df["datetime2"] = datetime.datetime(2001, 1, 3, 0, 0)
        df.loc[df.index[3:6], ["obj1"]] = np.nan
        df = df._consolidate()._convert(datetime=True)

        with ensure_clean_store(setup_path) as store:
            _maybe_remove(store, "df")

            # PerformanceWarning
            with catch_warnings(record=True):
                simplefilter("ignore", pd.errors.PerformanceWarning)
                store.put("df", df)

            expected = store.get("df")
            tm.assert_frame_equal(expected, df)

    @pytest.mark.filterwarnings(
        "ignore:object name:tables.exceptions.NaturalNameWarning"
    )
    def test_append(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # this is allowed by almost always don't want to do it
            # tables.NaturalNameWarning):
            with catch_warnings(record=True):

                df = tm.makeTimeDataFrame()
                _maybe_remove(store, "df1")
                store.append("df1", df[:10])
                store.append("df1", df[10:])
                tm.assert_frame_equal(store["df1"], df)

                _maybe_remove(store, "df2")
                store.put("df2", df[:10], format="table")
                store.append("df2", df[10:])
                tm.assert_frame_equal(store["df2"], df)

                _maybe_remove(store, "df3")
                store.append("/df3", df[:10])
                store.append("/df3", df[10:])
                tm.assert_frame_equal(store["df3"], df)

                # this is allowed by almost always don't want to do it
                # tables.NaturalNameWarning
                _maybe_remove(store, "/df3 foo")
                store.append("/df3 foo", df[:10])
                store.append("/df3 foo", df[10:])
                tm.assert_frame_equal(store["df3 foo"], df)

                # dtype issues - mizxed type in a single object column
                df = DataFrame(data=[[1, 2], [0, 1], [1, 2], [0, 0]])
                df["mixed_column"] = "testing"
                df.loc[2, "mixed_column"] = np.nan
                _maybe_remove(store, "df")
                store.append("df", df)
                tm.assert_frame_equal(store["df"], df)

                # uints - test storage of uints
                uint_data = DataFrame(
                    {
                        "u08": Series(
                            np.random.randint(0, high=255, size=5), dtype=np.uint8
                        ),
                        "u16": Series(
                            np.random.randint(0, high=65535, size=5), dtype=np.uint16
                        ),
                        "u32": Series(
                            np.random.randint(0, high=2 ** 30, size=5), dtype=np.uint32
                        ),
                        "u64": Series(
                            [2 ** 58, 2 ** 59, 2 ** 60, 2 ** 61, 2 ** 62],
                            dtype=np.uint64,
                        ),
                    },
                    index=np.arange(5),
                )
                _maybe_remove(store, "uints")
                store.append("uints", uint_data)
                tm.assert_frame_equal(store["uints"], uint_data)

                # uints - test storage of uints in indexable columns
                _maybe_remove(store, "uints")
                # 64-bit indices not yet supported
                store.append("uints", uint_data, data_columns=["u08", "u16", "u32"])
                tm.assert_frame_equal(store["uints"], uint_data)

    def test_append_series(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # basic
            ss = tm.makeStringSeries()
            ts = tm.makeTimeSeries()
            ns = Series(np.arange(100))

            store.append("ss", ss)
            result = store["ss"]
            tm.assert_series_equal(result, ss)
            assert result.name is None

            store.append("ts", ts)
            result = store["ts"]
            tm.assert_series_equal(result, ts)
            assert result.name is None

            ns.name = "foo"
            store.append("ns", ns)
            result = store["ns"]
            tm.assert_series_equal(result, ns)
            assert result.name == ns.name

            # select on the values
            expected = ns[ns > 60]
            result = store.select("ns", "foo>60")
            tm.assert_series_equal(result, expected)

            # select on the index and values
            expected = ns[(ns > 70) & (ns.index < 90)]
            result = store.select("ns", "foo>70 and index<90")
            tm.assert_series_equal(result, expected)

            # multi-index
            mi = DataFrame(np.random.randn(5, 1), columns=["A"])
            mi["B"] = np.arange(len(mi))
            mi["C"] = "foo"
            mi.loc[3:5, "C"] = "bar"
            mi.set_index(["C", "B"], inplace=True)
            s = mi.stack()
            s.index = s.index.droplevel(2)
            store.append("mi", s)
            tm.assert_series_equal(store["mi"], s)

    def test_store_index_types(self, setup_path):
        # GH5386
        # test storing various index types

        with ensure_clean_store(setup_path) as store:

            def check(format, index):
                df = DataFrame(np.random.randn(10, 2), columns=list("AB"))
                df.index = index(len(df))

                _maybe_remove(store, "df")
                store.put("df", df, format=format)
                tm.assert_frame_equal(df, store["df"])

            for index in [
                tm.makeFloatIndex,
                tm.makeStringIndex,
                tm.makeIntIndex,
                tm.makeDateIndex,
            ]:

                check("table", index)
                check("fixed", index)

            # period index currently broken for table
            # seee GH7796 FIXME
            check("fixed", tm.makePeriodIndex)
            # check('table',tm.makePeriodIndex)

            # unicode
            index = tm.makeUnicodeIndex
            check("table", index)
            check("fixed", index)

    @pytest.mark.skipif(
        not is_platform_little_endian(), reason="reason platform is not little endian"
    )
    def test_encoding(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            df = DataFrame(dict(A="foo", B="bar"), index=range(5))
            df.loc[2, "A"] = np.nan
            df.loc[3, "B"] = np.nan
            _maybe_remove(store, "df")
            store.append("df", df, encoding="ascii")
            tm.assert_frame_equal(store["df"], df)

            expected = df.reindex(columns=["A"])
            result = store.select("df", Term("columns=A", encoding="ascii"))
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "val",
        [
            [b"E\xc9, 17", b"", b"a", b"b", b"c"],
            [b"E\xc9, 17", b"a", b"b", b"c"],
            [b"EE, 17", b"", b"a", b"b", b"c"],
            [b"E\xc9, 17", b"\xf8\xfc", b"a", b"b", b"c"],
            [b"", b"a", b"b", b"c"],
            [b"\xf8\xfc", b"a", b"b", b"c"],
            [b"A\xf8\xfc", b"", b"a", b"b", b"c"],
            [np.nan, b"", b"b", b"c"],
            [b"A\xf8\xfc", np.nan, b"", b"b", b"c"],
        ],
    )
    @pytest.mark.parametrize("dtype", ["category", object])
    def test_latin_encoding(self, setup_path, dtype, val):
        enc = "latin-1"
        nan_rep = ""
        key = "data"

        val = [x.decode(enc) if isinstance(x, bytes) else x for x in val]
        ser = pd.Series(val, dtype=dtype)

        with ensure_clean_path(setup_path) as store:
            ser.to_hdf(store, key, format="table", encoding=enc, nan_rep=nan_rep)
            retr = read_hdf(store, key)

        s_nan = ser.replace(nan_rep, np.nan)

        tm.assert_series_equal(s_nan, retr)

    def test_append_some_nans(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            df = DataFrame(
                {
                    "A": Series(np.random.randn(20)).astype("int32"),
                    "A1": np.random.randn(20),
                    "A2": np.random.randn(20),
                    "B": "foo",
                    "C": "bar",
                    "D": Timestamp("20010101"),
                    "E": datetime.datetime(2001, 1, 2, 0, 0),
                },
                index=np.arange(20),
            )
            # some nans
            _maybe_remove(store, "df1")
            df.loc[0:15, ["A1", "B", "D", "E"]] = np.nan
            store.append("df1", df[:10])
            store.append("df1", df[10:])
            tm.assert_frame_equal(store["df1"], df)

            # first column
            df1 = df.copy()
            df1.loc[:, "A1"] = np.nan
            _maybe_remove(store, "df1")
            store.append("df1", df1[:10])
            store.append("df1", df1[10:])
            tm.assert_frame_equal(store["df1"], df1)

            # 2nd column
            df2 = df.copy()
            df2.loc[:, "A2"] = np.nan
            _maybe_remove(store, "df2")
            store.append("df2", df2[:10])
            store.append("df2", df2[10:])
            tm.assert_frame_equal(store["df2"], df2)

            # datetimes
            df3 = df.copy()
            df3.loc[:, "E"] = np.nan
            _maybe_remove(store, "df3")
            store.append("df3", df3[:10])
            store.append("df3", df3[10:])
            tm.assert_frame_equal(store["df3"], df3)

    def test_append_all_nans(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            df = DataFrame(
                {"A1": np.random.randn(20), "A2": np.random.randn(20)},
                index=np.arange(20),
            )
            df.loc[0:15, :] = np.nan

            # nan some entire rows (dropna=True)
            _maybe_remove(store, "df")
            store.append("df", df[:10], dropna=True)
            store.append("df", df[10:], dropna=True)
            tm.assert_frame_equal(store["df"], df[-4:])

            # nan some entire rows (dropna=False)
            _maybe_remove(store, "df2")
            store.append("df2", df[:10], dropna=False)
            store.append("df2", df[10:], dropna=False)
            tm.assert_frame_equal(store["df2"], df)

            # tests the option io.hdf.dropna_table
            pd.set_option("io.hdf.dropna_table", False)
            _maybe_remove(store, "df3")
            store.append("df3", df[:10])
            store.append("df3", df[10:])
            tm.assert_frame_equal(store["df3"], df)

            pd.set_option("io.hdf.dropna_table", True)
            _maybe_remove(store, "df4")
            store.append("df4", df[:10])
            store.append("df4", df[10:])
            tm.assert_frame_equal(store["df4"], df[-4:])

            # nan some entire rows (string are still written!)
            df = DataFrame(
                {
                    "A1": np.random.randn(20),
                    "A2": np.random.randn(20),
                    "B": "foo",
                    "C": "bar",
                },
                index=np.arange(20),
            )

            df.loc[0:15, :] = np.nan

            _maybe_remove(store, "df")
            store.append("df", df[:10], dropna=True)
            store.append("df", df[10:], dropna=True)
            tm.assert_frame_equal(store["df"], df)

            _maybe_remove(store, "df2")
            store.append("df2", df[:10], dropna=False)
            store.append("df2", df[10:], dropna=False)
            tm.assert_frame_equal(store["df2"], df)

            # nan some entire rows (but since we have dates they are still
            # written!)
            df = DataFrame(
                {
                    "A1": np.random.randn(20),
                    "A2": np.random.randn(20),
                    "B": "foo",
                    "C": "bar",
                    "D": Timestamp("20010101"),
                    "E": datetime.datetime(2001, 1, 2, 0, 0),
                },
                index=np.arange(20),
            )

            df.loc[0:15, :] = np.nan

            _maybe_remove(store, "df")
            store.append("df", df[:10], dropna=True)
            store.append("df", df[10:], dropna=True)
            tm.assert_frame_equal(store["df"], df)

            _maybe_remove(store, "df2")
            store.append("df2", df[:10], dropna=False)
            store.append("df2", df[10:], dropna=False)
            tm.assert_frame_equal(store["df2"], df)

        # Test to make sure defaults are to not drop.
        # Corresponding to Issue 9382
        df_with_missing = DataFrame(
            {"col1": [0, np.nan, 2], "col2": [1, np.nan, np.nan]}
        )

        with ensure_clean_path(setup_path) as path:
            df_with_missing.to_hdf(path, "df_with_missing", format="table")
            reloaded = read_hdf(path, "df_with_missing")
            tm.assert_frame_equal(df_with_missing, reloaded)

    def test_read_missing_key_close_store(self, setup_path):
        # GH 25766
        with ensure_clean_path(setup_path) as path:
            df = pd.DataFrame({"a": range(2), "b": range(2)})
            df.to_hdf(path, "k1")

            with pytest.raises(KeyError, match="'No object named k2 in the file'"):
                pd.read_hdf(path, "k2")

            # smoke test to test that file is properly closed after
            # read with KeyError before another write
            df.to_hdf(path, "k2")

    def test_read_missing_key_opened_store(self, setup_path):
        # GH 28699
        with ensure_clean_path(setup_path) as path:
            df = pd.DataFrame({"a": range(2), "b": range(2)})
            df.to_hdf(path, "k1")

            store = pd.HDFStore(path, "r")

            with pytest.raises(KeyError, match="'No object named k2 in the file'"):
                pd.read_hdf(store, "k2")

            # Test that the file is still open after a KeyError and that we can
            # still read from it.
            pd.read_hdf(store, "k1")

    def test_append_frame_column_oriented(self, setup_path):
        with ensure_clean_store(setup_path) as store:

            # column oriented
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, "df1")
            store.append("df1", df.iloc[:, :2], axes=["columns"])
            store.append("df1", df.iloc[:, 2:])
            tm.assert_frame_equal(store["df1"], df)

            result = store.select("df1", "columns=A")
            expected = df.reindex(columns=["A"])
            tm.assert_frame_equal(expected, result)

            # selection on the non-indexable
            result = store.select("df1", ("columns=A", "index=df.index[0:4]"))
            expected = df.reindex(columns=["A"], index=df.index[0:4])
            tm.assert_frame_equal(expected, result)

            # this isn't supported
            with pytest.raises(TypeError):
                store.select("df1", "columns=A and index>df.index[4]")

    def test_append_with_different_block_ordering(self, setup_path):

        # GH 4096; using same frames, but different block orderings
        with ensure_clean_store(setup_path) as store:

            for i in range(10):

                df = DataFrame(np.random.randn(10, 2), columns=list("AB"))
                df["index"] = range(10)
                df["index"] += i * 10
                df["int64"] = Series([1] * len(df), dtype="int64")
                df["int16"] = Series([1] * len(df), dtype="int16")

                if i % 2 == 0:
                    del df["int64"]
                    df["int64"] = Series([1] * len(df), dtype="int64")
                if i % 3 == 0:
                    a = df.pop("A")
                    df["A"] = a

                df.set_index("index", inplace=True)

                store.append("df", df)

        # test a different ordering but with more fields (like invalid
        # combinate)
        with ensure_clean_store(setup_path) as store:

            df = DataFrame(np.random.randn(10, 2), columns=list("AB"), dtype="float64")
            df["int64"] = Series([1] * len(df), dtype="int64")
            df["int16"] = Series([1] * len(df), dtype="int16")
            store.append("df", df)

            # store additional fields in different blocks
            df["int16_2"] = Series([1] * len(df), dtype="int16")
            with pytest.raises(ValueError):
                store.append("df", df)

            # store multiple additional fields in different blocks
            df["float_3"] = Series([1.0] * len(df), dtype="float64")
            with pytest.raises(ValueError):
                store.append("df", df)

    def test_append_with_strings(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            with catch_warnings(record=True):

                def check_col(key, name, size):
                    assert (
                        getattr(store.get_storer(key).table.description, name).itemsize
                        == size
                    )

                # avoid truncation on elements
                df = DataFrame([[123, "asdqwerty"], [345, "dggnhebbsdfbdfb"]])
                store.append("df_big", df)
                tm.assert_frame_equal(store.select("df_big"), df)
                check_col("df_big", "values_block_1", 15)

                # appending smaller string ok
                df2 = DataFrame([[124, "asdqy"], [346, "dggnhefbdfb"]])
                store.append("df_big", df2)
                expected = concat([df, df2])
                tm.assert_frame_equal(store.select("df_big"), expected)
                check_col("df_big", "values_block_1", 15)

                # avoid truncation on elements
                df = DataFrame([[123, "asdqwerty"], [345, "dggnhebbsdfbdfb"]])
                store.append("df_big2", df, min_itemsize={"values": 50})
                tm.assert_frame_equal(store.select("df_big2"), df)
                check_col("df_big2", "values_block_1", 50)

                # bigger string on next append
                store.append("df_new", df)
                df_new = DataFrame(
                    [[124, "abcdefqhij"], [346, "abcdefghijklmnopqrtsuvwxyz"]]
                )
                with pytest.raises(ValueError):
                    store.append("df_new", df_new)

                # min_itemsize on Series index (GH 11412)
                df = tm.makeMixedDataFrame().set_index("C")
                store.append("ss", df["B"], min_itemsize={"index": 4})
                tm.assert_series_equal(store.select("ss"), df["B"])

                # same as above, with data_columns=True
                store.append(
                    "ss2", df["B"], data_columns=True, min_itemsize={"index": 4}
                )
                tm.assert_series_equal(store.select("ss2"), df["B"])

                # min_itemsize in index without appending (GH 10381)
                store.put("ss3", df, format="table", min_itemsize={"index": 6})
                # just make sure there is a longer string:
                df2 = df.copy().reset_index().assign(C="longer").set_index("C")
                store.append("ss3", df2)
                tm.assert_frame_equal(store.select("ss3"), pd.concat([df, df2]))

                # same as above, with a Series
                store.put("ss4", df["B"], format="table", min_itemsize={"index": 6})
                store.append("ss4", df2["B"])
                tm.assert_series_equal(
                    store.select("ss4"), pd.concat([df["B"], df2["B"]])
                )

                # with nans
                _maybe_remove(store, "df")
                df = tm.makeTimeDataFrame()
                df["string"] = "foo"
                df.loc[df.index[1:4], "string"] = np.nan
                df["string2"] = "bar"
                df.loc[df.index[4:8], "string2"] = np.nan
                df["string3"] = "bah"
                df.loc[df.index[1:], "string3"] = np.nan
                store.append("df", df)
                result = store.select("df")
                tm.assert_frame_equal(result, df)

        with ensure_clean_store(setup_path) as store:

            def check_col(key, name, size):
                assert getattr(
                    store.get_storer(key).table.description, name
                ).itemsize, size

            df = DataFrame(dict(A="foo", B="bar"), index=range(10))

            # a min_itemsize that creates a data_column
            _maybe_remove(store, "df")
            store.append("df", df, min_itemsize={"A": 200})
            check_col("df", "A", 200)
            assert store.get_storer("df").data_columns == ["A"]

            # a min_itemsize that creates a data_column2
            _maybe_remove(store, "df")
            store.append("df", df, data_columns=["B"], min_itemsize={"A": 200})
            check_col("df", "A", 200)
            assert store.get_storer("df").data_columns == ["B", "A"]

            # a min_itemsize that creates a data_column2
            _maybe_remove(store, "df")
            store.append("df", df, data_columns=["B"], min_itemsize={"values": 200})
            check_col("df", "B", 200)
            check_col("df", "values_block_0", 200)
            assert store.get_storer("df").data_columns == ["B"]

            # infer the .typ on subsequent appends
            _maybe_remove(store, "df")
            store.append("df", df[:5], min_itemsize=200)
            store.append("df", df[5:], min_itemsize=200)
            tm.assert_frame_equal(store["df"], df)

            # invalid min_itemsize keys
            df = DataFrame(["foo", "foo", "foo", "barh", "barh", "barh"], columns=["A"])
            _maybe_remove(store, "df")
            with pytest.raises(ValueError):
                store.append("df", df, min_itemsize={"foo": 20, "foobar": 20})

    def test_append_with_empty_string(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # with all empty strings (GH 12242)
            df = DataFrame({"x": ["a", "b", "c", "d", "e", "f", ""]})
            store.append("df", df[:-1], min_itemsize={"x": 1})
            store.append("df", df[-1:], min_itemsize={"x": 1})
            tm.assert_frame_equal(store.select("df"), df)

    def test_to_hdf_with_min_itemsize(self, setup_path):

        with ensure_clean_path(setup_path) as path:

            # min_itemsize in index with to_hdf (GH 10381)
            df = tm.makeMixedDataFrame().set_index("C")
            df.to_hdf(path, "ss3", format="table", min_itemsize={"index": 6})
            # just make sure there is a longer string:
            df2 = df.copy().reset_index().assign(C="longer").set_index("C")
            df2.to_hdf(path, "ss3", append=True, format="table")
            tm.assert_frame_equal(pd.read_hdf(path, "ss3"), pd.concat([df, df2]))

            # same as above, with a Series
            df["B"].to_hdf(path, "ss4", format="table", min_itemsize={"index": 6})
            df2["B"].to_hdf(path, "ss4", append=True, format="table")
            tm.assert_series_equal(
                pd.read_hdf(path, "ss4"), pd.concat([df["B"], df2["B"]])
            )

    @pytest.mark.parametrize(
        "format", [pytest.param("fixed", marks=td.xfail_non_writeable), "table"]
    )
    def test_to_hdf_errors(self, format, setup_path):

        data = ["\ud800foo"]
        ser = pd.Series(data, index=pd.Index(data))
        with ensure_clean_path(setup_path) as path:
            # GH 20835
            ser.to_hdf(path, "table", format=format, errors="surrogatepass")

            result = pd.read_hdf(path, "table", errors="surrogatepass")
            tm.assert_series_equal(result, ser)

    def test_append_with_data_columns(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            df = tm.makeTimeDataFrame()
            df.iloc[0, df.columns.get_loc("B")] = 1.0
            _maybe_remove(store, "df")
            store.append("df", df[:2], data_columns=["B"])
            store.append("df", df[2:])
            tm.assert_frame_equal(store["df"], df)

            # check that we have indices created
            assert store._handle.root.df.table.cols.index.is_indexed is True
            assert store._handle.root.df.table.cols.B.is_indexed is True

            # data column searching
            result = store.select("df", "B>0")
            expected = df[df.B > 0]
            tm.assert_frame_equal(result, expected)

            # data column searching (with an indexable and a data_columns)
            result = store.select("df", "B>0 and index>df.index[3]")
            df_new = df.reindex(index=df.index[4:])
            expected = df_new[df_new.B > 0]
            tm.assert_frame_equal(result, expected)

            # data column selection with a string data_column
            df_new = df.copy()
            df_new["string"] = "foo"
            df_new.loc[df_new.index[1:4], "string"] = np.nan
            df_new.loc[df_new.index[5:6], "string"] = "bar"
            _maybe_remove(store, "df")
            store.append("df", df_new, data_columns=["string"])
            result = store.select("df", "string='foo'")
            expected = df_new[df_new.string == "foo"]
            tm.assert_frame_equal(result, expected)

            # using min_itemsize and a data column
            def check_col(key, name, size):
                assert (
                    getattr(store.get_storer(key).table.description, name).itemsize
                    == size
                )

        with ensure_clean_store(setup_path) as store:
            _maybe_remove(store, "df")
            store.append(
                "df", df_new, data_columns=["string"], min_itemsize={"string": 30}
            )
            check_col("df", "string", 30)
            _maybe_remove(store, "df")
            store.append("df", df_new, data_columns=["string"], min_itemsize=30)
            check_col("df", "string", 30)
            _maybe_remove(store, "df")
            store.append(
                "df", df_new, data_columns=["string"], min_itemsize={"values": 30}
            )
            check_col("df", "string", 30)

        with ensure_clean_store(setup_path) as store:
            df_new["string2"] = "foobarbah"
            df_new["string_block1"] = "foobarbah1"
            df_new["string_block2"] = "foobarbah2"
            _maybe_remove(store, "df")
            store.append(
                "df",
                df_new,
                data_columns=["string", "string2"],
                min_itemsize={"string": 30, "string2": 40, "values": 50},
            )
            check_col("df", "string", 30)
            check_col("df", "string2", 40)
            check_col("df", "values_block_1", 50)

        with ensure_clean_store(setup_path) as store:
            # multiple data columns
            df_new = df.copy()
            df_new.iloc[0, df_new.columns.get_loc("A")] = 1.0
            df_new.iloc[0, df_new.columns.get_loc("B")] = -1.0
            df_new["string"] = "foo"

            sl = df_new.columns.get_loc("string")
            df_new.iloc[1:4, sl] = np.nan
            df_new.iloc[5:6, sl] = "bar"

            df_new["string2"] = "foo"
            sl = df_new.columns.get_loc("string2")
            df_new.iloc[2:5, sl] = np.nan
            df_new.iloc[7:8, sl] = "bar"
            _maybe_remove(store, "df")
            store.append("df", df_new, data_columns=["A", "B", "string", "string2"])
            result = store.select(
                "df", "string='foo' and string2='foo' and A>0 and B<0"
            )
            expected = df_new[
                (df_new.string == "foo")
                & (df_new.string2 == "foo")
                & (df_new.A > 0)
                & (df_new.B < 0)
            ]
            tm.assert_frame_equal(result, expected, check_index_type=False)

            # yield an empty frame
            result = store.select("df", "string='foo' and string2='cool'")
            expected = df_new[(df_new.string == "foo") & (df_new.string2 == "cool")]
            tm.assert_frame_equal(result, expected, check_index_type=False)

        with ensure_clean_store(setup_path) as store:
            # doc example
            df_dc = df.copy()
            df_dc["string"] = "foo"
            df_dc.loc[df_dc.index[4:6], "string"] = np.nan
            df_dc.loc[df_dc.index[7:9], "string"] = "bar"
            df_dc["string2"] = "cool"
            df_dc["datetime"] = Timestamp("20010102")
            df_dc = df_dc._convert(datetime=True)
            df_dc.loc[df_dc.index[3:5], ["A", "B", "datetime"]] = np.nan

            _maybe_remove(store, "df_dc")
            store.append(
                "df_dc", df_dc, data_columns=["B", "C", "string", "string2", "datetime"]
            )
            result = store.select("df_dc", "B>0")

            expected = df_dc[df_dc.B > 0]
            tm.assert_frame_equal(result, expected, check_index_type=False)

            result = store.select("df_dc", ["B > 0", "C > 0", "string == foo"])
            expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (df_dc.string == "foo")]
            tm.assert_frame_equal(result, expected, check_index_type=False)

        with ensure_clean_store(setup_path) as store:
            # doc example part 2
            np.random.seed(1234)
            index = date_range("1/1/2000", periods=8)
            df_dc = DataFrame(
                np.random.randn(8, 3), index=index, columns=["A", "B", "C"]
            )
            df_dc["string"] = "foo"
            df_dc.loc[df_dc.index[4:6], "string"] = np.nan
            df_dc.loc[df_dc.index[7:9], "string"] = "bar"
            df_dc.loc[:, ["B", "C"]] = df_dc.loc[:, ["B", "C"]].abs()
            df_dc["string2"] = "cool"

            # on-disk operations
            store.append("df_dc", df_dc, data_columns=["B", "C", "string", "string2"])

            result = store.select("df_dc", "B>0")
            expected = df_dc[df_dc.B > 0]
            tm.assert_frame_equal(result, expected)

            result = store.select("df_dc", ["B > 0", "C > 0", 'string == "foo"'])
            expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (df_dc.string == "foo")]
            tm.assert_frame_equal(result, expected)

    def test_create_table_index(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            with catch_warnings(record=True):

                def col(t, column):
                    return getattr(store.get_storer(t).table.cols, column)

                # data columns
                df = tm.makeTimeDataFrame()
                df["string"] = "foo"
                df["string2"] = "bar"
                store.append("f", df, data_columns=["string", "string2"])
                assert col("f", "index").is_indexed is True
                assert col("f", "string").is_indexed is True
                assert col("f", "string2").is_indexed is True

                # specify index=columns
                store.append(
                    "f2", df, index=["string"], data_columns=["string", "string2"]
                )
                assert col("f2", "index").is_indexed is False
                assert col("f2", "string").is_indexed is True
                assert col("f2", "string2").is_indexed is False

                # try to index a non-table
                _maybe_remove(store, "f2")
                store.put("f2", df)
                with pytest.raises(TypeError):
                    store.create_table_index("f2")

    def test_append_hierarchical(self, setup_path):
        index = MultiIndex(
            levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
            codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
            names=["foo", "bar"],
        )
        df = DataFrame(np.random.randn(10, 3), index=index, columns=["A", "B", "C"])

        with ensure_clean_store(setup_path) as store:
            store.append("mi", df)
            result = store.select("mi")
            tm.assert_frame_equal(result, df)

            # GH 3748
            result = store.select("mi", columns=["A", "B"])
            expected = df.reindex(columns=["A", "B"])
            tm.assert_frame_equal(result, expected)

        with ensure_clean_path("test.hdf") as path:
            df.to_hdf(path, "df", format="table")
            result = read_hdf(path, "df", columns=["A", "B"])
            expected = df.reindex(columns=["A", "B"])
            tm.assert_frame_equal(result, expected)

    def test_column_multiindex(self, setup_path):
        # GH 4710
        # recreate multi-indexes properly

        index = MultiIndex.from_tuples(
            [("A", "a"), ("A", "b"), ("B", "a"), ("B", "b")], names=["first", "second"]
        )
        df = DataFrame(np.arange(12).reshape(3, 4), columns=index)
        expected = df.copy()
        if isinstance(expected.index, RangeIndex):
            expected.index = Int64Index(expected.index)

        with ensure_clean_store(setup_path) as store:

            store.put("df", df)
            tm.assert_frame_equal(
                store["df"], expected, check_index_type=True, check_column_type=True
            )

            store.put("df1", df, format="table")
            tm.assert_frame_equal(
                store["df1"], expected, check_index_type=True, check_column_type=True
            )

            with pytest.raises(ValueError):
                store.put("df2", df, format="table", data_columns=["A"])
            with pytest.raises(ValueError):
                store.put("df3", df, format="table", data_columns=True)

        # appending multi-column on existing table (see GH 6167)
        with ensure_clean_store(setup_path) as store:
            store.append("df2", df)
            store.append("df2", df)

            tm.assert_frame_equal(store["df2"], concat((df, df)))

        # non_index_axes name
        df = DataFrame(
            np.arange(12).reshape(3, 4), columns=Index(list("ABCD"), name="foo")
        )
        expected = df.copy()
        if isinstance(expected.index, RangeIndex):
            expected.index = Int64Index(expected.index)

        with ensure_clean_store(setup_path) as store:

            store.put("df1", df, format="table")
            tm.assert_frame_equal(
                store["df1"], expected, check_index_type=True, check_column_type=True
            )

    def test_store_multiindex(self, setup_path):

        # validate multi-index names
        # GH 5527
        with ensure_clean_store(setup_path) as store:

            def make_index(names=None):
                return MultiIndex.from_tuples(
                    [
                        (datetime.datetime(2013, 12, d), s, t)
                        for d in range(1, 3)
                        for s in range(2)
                        for t in range(3)
                    ],
                    names=names,
                )

            # no names
            _maybe_remove(store, "df")
            df = DataFrame(np.zeros((12, 2)), columns=["a", "b"], index=make_index())
            store.append("df", df)
            tm.assert_frame_equal(store.select("df"), df)

            # partial names
            _maybe_remove(store, "df")
            df = DataFrame(
                np.zeros((12, 2)),
                columns=["a", "b"],
                index=make_index(["date", None, None]),
            )
            store.append("df", df)
            tm.assert_frame_equal(store.select("df"), df)

            # series
            _maybe_remove(store, "s")
            s = Series(np.zeros(12), index=make_index(["date", None, None]))
            store.append("s", s)
            xp = Series(np.zeros(12), index=make_index(["date", "level_1", "level_2"]))
            tm.assert_series_equal(store.select("s"), xp)

            # dup with column
            _maybe_remove(store, "df")
            df = DataFrame(
                np.zeros((12, 2)),
                columns=["a", "b"],
                index=make_index(["date", "a", "t"]),
            )
            with pytest.raises(ValueError):
                store.append("df", df)

            # dup within level
            _maybe_remove(store, "df")
            df = DataFrame(
                np.zeros((12, 2)),
                columns=["a", "b"],
                index=make_index(["date", "date", "date"]),
            )
            with pytest.raises(ValueError):
                store.append("df", df)

            # fully names
            _maybe_remove(store, "df")
            df = DataFrame(
                np.zeros((12, 2)),
                columns=["a", "b"],
                index=make_index(["date", "s", "t"]),
            )
            store.append("df", df)
            tm.assert_frame_equal(store.select("df"), df)

    def test_select_columns_in_where(self, setup_path):

        # GH 6169
        # recreate multi-indexes when columns is passed
        # in the `where` argument
        index = MultiIndex(
            levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
            codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
            names=["foo_name", "bar_name"],
        )

        # With a DataFrame
        df = DataFrame(np.random.randn(10, 3), index=index, columns=["A", "B", "C"])

        with ensure_clean_store(setup_path) as store:
            store.put("df", df, format="table")
            expected = df[["A"]]

            tm.assert_frame_equal(store.select("df", columns=["A"]), expected)

            tm.assert_frame_equal(store.select("df", where="columns=['A']"), expected)

        # With a Series
        s = Series(np.random.randn(10), index=index, name="A")
        with ensure_clean_store(setup_path) as store:
            store.put("s", s, format="table")
            tm.assert_series_equal(store.select("s", where="columns=['A']"), s)

    def test_mi_data_columns(self, setup_path):
        # GH 14435
        idx = pd.MultiIndex.from_arrays(
            [date_range("2000-01-01", periods=5), range(5)], names=["date", "id"]
        )
        df = pd.DataFrame({"a": [1.1, 1.2, 1.3, 1.4, 1.5]}, index=idx)

        with ensure_clean_store(setup_path) as store:
            store.append("df", df, data_columns=True)

            actual = store.select("df", where="id == 1")
            expected = df.iloc[[1], :]
            tm.assert_frame_equal(actual, expected)

    def test_pass_spec_to_storer(self, setup_path):

        df = tm.makeDataFrame()

        with ensure_clean_store(setup_path) as store:
            store.put("df", df)
            with pytest.raises(TypeError):
                store.select("df", columns=["A"])
            with pytest.raises(TypeError):
                store.select("df", where=[("columns=A")])

    @td.xfail_non_writeable
    def test_append_misc(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            df = tm.makeDataFrame()
            store.append("df", df, chunksize=1)
            result = store.select("df")
            tm.assert_frame_equal(result, df)

            store.append("df1", df, expectedrows=10)
            result = store.select("df1")
            tm.assert_frame_equal(result, df)

        # more chunksize in append tests
        def check(obj, comparator):
            for c in [10, 200, 1000]:
                with ensure_clean_store(setup_path, mode="w") as store:
                    store.append("obj", obj, chunksize=c)
                    result = store.select("obj")
                    comparator(result, obj)

        df = tm.makeDataFrame()
        df["string"] = "foo"
        df["float322"] = 1.0
        df["float322"] = df["float322"].astype("float32")
        df["bool"] = df["float322"] > 0
        df["time1"] = Timestamp("20130101")
        df["time2"] = Timestamp("20130102")
        check(df, tm.assert_frame_equal)

        # empty frame, GH4273
        with ensure_clean_store(setup_path) as store:

            # 0 len
            df_empty = DataFrame(columns=list("ABC"))
            store.append("df", df_empty)
            with pytest.raises(KeyError, match="'No object named df in the file'"):
                store.select("df")

            # repeated append of 0/non-zero frames
            df = DataFrame(np.random.rand(10, 3), columns=list("ABC"))
            store.append("df", df)
            tm.assert_frame_equal(store.select("df"), df)
            store.append("df", df_empty)
            tm.assert_frame_equal(store.select("df"), df)

            # store
            df = DataFrame(columns=list("ABC"))
            store.put("df2", df)
            tm.assert_frame_equal(store.select("df2"), df)

    def test_append_raise(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # test append with invalid input to get good error messages

            # list in column
            df = tm.makeDataFrame()
            df["invalid"] = [["a"]] * len(df)
            assert df.dtypes["invalid"] == np.object_
            with pytest.raises(TypeError):
                store.append("df", df)

            # multiple invalid columns
            df["invalid2"] = [["a"]] * len(df)
            df["invalid3"] = [["a"]] * len(df)
            with pytest.raises(TypeError):
                store.append("df", df)

            # datetime with embedded nans as object
            df = tm.makeDataFrame()
            s = Series(datetime.datetime(2001, 1, 2), index=df.index)
            s = s.astype(object)
            s[0:5] = np.nan
            df["invalid"] = s
            assert df.dtypes["invalid"] == np.object_
            with pytest.raises(TypeError):
                store.append("df", df)

            # directly ndarray
            with pytest.raises(TypeError):
                store.append("df", np.arange(10))

            # series directly
            with pytest.raises(TypeError):
                store.append("df", Series(np.arange(10)))

            # appending an incompatible table
            df = tm.makeDataFrame()
            store.append("df", df)

            df["foo"] = "foo"
            with pytest.raises(ValueError):
                store.append("df", df)

    def test_table_index_incompatible_dtypes(self, setup_path):
        df1 = DataFrame({"a": [1, 2, 3]})
        df2 = DataFrame({"a": [4, 5, 6]}, index=date_range("1/1/2000", periods=3))

        with ensure_clean_store(setup_path) as store:
            store.put("frame", df1, format="table")
            with pytest.raises(TypeError):
                store.put("frame", df2, format="table", append=True)

    def test_table_values_dtypes_roundtrip(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            df1 = DataFrame({"a": [1, 2, 3]}, dtype="f8")
            store.append("df_f8", df1)
            tm.assert_series_equal(df1.dtypes, store["df_f8"].dtypes)

            df2 = DataFrame({"a": [1, 2, 3]}, dtype="i8")
            store.append("df_i8", df2)
            tm.assert_series_equal(df2.dtypes, store["df_i8"].dtypes)

            # incompatible dtype
            with pytest.raises(ValueError):
                store.append("df_i8", df1)

            # check creation/storage/retrieval of float32 (a bit hacky to
            # actually create them thought)
            df1 = DataFrame(np.array([[1], [2], [3]], dtype="f4"), columns=["A"])
            store.append("df_f4", df1)
            tm.assert_series_equal(df1.dtypes, store["df_f4"].dtypes)
            assert df1.dtypes[0] == "float32"

            # check with mixed dtypes
            df1 = DataFrame(
                {
                    c: Series(np.random.randint(5), dtype=c)
                    for c in ["float32", "float64", "int32", "int64", "int16", "int8"]
                }
            )
            df1["string"] = "foo"
            df1["float322"] = 1.0
            df1["float322"] = df1["float322"].astype("float32")
            df1["bool"] = df1["float32"] > 0
            df1["time1"] = Timestamp("20130101")
            df1["time2"] = Timestamp("20130102")

            store.append("df_mixed_dtypes1", df1)
            result = store.select("df_mixed_dtypes1").dtypes.value_counts()
            result.index = [str(i) for i in result.index]
            expected = Series(
                {
                    "float32": 2,
                    "float64": 1,
                    "int32": 1,
                    "bool": 1,
                    "int16": 1,
                    "int8": 1,
                    "int64": 1,
                    "object": 1,
                    "datetime64[ns]": 2,
                }
            )
            result = result.sort_index()
            expected = expected.sort_index()
            tm.assert_series_equal(result, expected)

    def test_table_mixed_dtypes(self, setup_path):

        # frame
        df = tm.makeDataFrame()
        df["obj1"] = "foo"
        df["obj2"] = "bar"
        df["bool1"] = df["A"] > 0
        df["bool2"] = df["B"] > 0
        df["bool3"] = True
        df["int1"] = 1
        df["int2"] = 2
        df["timestamp1"] = Timestamp("20010102")
        df["timestamp2"] = Timestamp("20010103")
        df["datetime1"] = datetime.datetime(2001, 1, 2, 0, 0)
        df["datetime2"] = datetime.datetime(2001, 1, 3, 0, 0)
        df.loc[df.index[3:6], ["obj1"]] = np.nan
        df = df._consolidate()._convert(datetime=True)

        with ensure_clean_store(setup_path) as store:
            store.append("df1_mixed", df)
            tm.assert_frame_equal(store.select("df1_mixed"), df)

    def test_unimplemented_dtypes_table_columns(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            dtypes = [("date", datetime.date(2001, 1, 2))]

            # currently not supported dtypes ####
            for n, f in dtypes:
                df = tm.makeDataFrame()
                df[n] = f
                with pytest.raises(TypeError):
                    store.append(f"df1_{n}", df)

        # frame
        df = tm.makeDataFrame()
        df["obj1"] = "foo"
        df["obj2"] = "bar"
        df["datetime1"] = datetime.date(2001, 1, 2)
        df = df._consolidate()._convert(datetime=True)

        with ensure_clean_store(setup_path) as store:
            # this fails because we have a date in the object block......
            with pytest.raises(TypeError):
                store.append("df_unimplemented", df)

    @td.xfail_non_writeable
    @pytest.mark.skipif(
        LooseVersion(np.__version__) == LooseVersion("1.15.0"),
        reason=(
            "Skipping  pytables test when numpy version is "
            "exactly equal to 1.15.0: gh-22098"
        ),
    )
    def test_calendar_roundtrip_issue(self, setup_path):

        # 8591
        # doc example from tseries holiday section
        weekmask_egypt = "Sun Mon Tue Wed Thu"
        holidays = [
            "2012-05-01",
            datetime.datetime(2013, 5, 1),
            np.datetime64("2014-05-01"),
        ]
        bday_egypt = pd.offsets.CustomBusinessDay(
            holidays=holidays, weekmask=weekmask_egypt
        )
        dt = datetime.datetime(2013, 4, 30)
        dts = date_range(dt, periods=5, freq=bday_egypt)

        s = Series(dts.weekday, dts).map(Series("Mon Tue Wed Thu Fri Sat Sun".split()))

        with ensure_clean_store(setup_path) as store:

            store.put("fixed", s)
            result = store.select("fixed")
            tm.assert_series_equal(result, s)

            store.append("table", s)
            result = store.select("table")
            tm.assert_series_equal(result, s)

    def test_roundtrip_tz_aware_index(self, setup_path):
        # GH 17618
        time = pd.Timestamp("2000-01-01 01:00:00", tz="US/Eastern")
        df = pd.DataFrame(data=[0], index=[time])

        with ensure_clean_store(setup_path) as store:
            store.put("frame", df, format="fixed")
            recons = store["frame"]
            tm.assert_frame_equal(recons, df)
            assert recons.index[0].value == 946706400000000000

    def test_append_with_timedelta(self, setup_path):
        # GH 3577
        # append timedelta

        df = DataFrame(
            dict(
                A=Timestamp("20130101"),
                B=[
                    Timestamp("20130101") + timedelta(days=i, seconds=10)
                    for i in range(10)
                ],
            )
        )
        df["C"] = df["A"] - df["B"]
        df.loc[3:5, "C"] = np.nan

        with ensure_clean_store(setup_path) as store:

            # table
            _maybe_remove(store, "df")
            store.append("df", df, data_columns=True)
            result = store.select("df")
            tm.assert_frame_equal(result, df)

            result = store.select("df", where="C<100000")
            tm.assert_frame_equal(result, df)

            result = store.select("df", where="C<pd.Timedelta('-3D')")
            tm.assert_frame_equal(result, df.iloc[3:])

            result = store.select("df", "C<'-3D'")
            tm.assert_frame_equal(result, df.iloc[3:])

            # a bit hacky here as we don't really deal with the NaT properly

            result = store.select("df", "C<'-500000s'")
            result = result.dropna(subset=["C"])
            tm.assert_frame_equal(result, df.iloc[6:])

            result = store.select("df", "C<'-3.5D'")
            result = result.iloc[1:]
            tm.assert_frame_equal(result, df.iloc[4:])

            # fixed
            _maybe_remove(store, "df2")
            store.put("df2", df)
            result = store.select("df2")
            tm.assert_frame_equal(result, df)

    def test_remove(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            ts = tm.makeTimeSeries()
            df = tm.makeDataFrame()
            store["a"] = ts
            store["b"] = df
            _maybe_remove(store, "a")
            assert len(store) == 1
            tm.assert_frame_equal(df, store["b"])

            _maybe_remove(store, "b")
            assert len(store) == 0

            # nonexistence
            with pytest.raises(
                KeyError, match="'No object named a_nonexistent_store in the file'"
            ):
                store.remove("a_nonexistent_store")

            # pathing
            store["a"] = ts
            store["b/foo"] = df
            _maybe_remove(store, "foo")
            _maybe_remove(store, "b/foo")
            assert len(store) == 1

            store["a"] = ts
            store["b/foo"] = df
            _maybe_remove(store, "b")
            assert len(store) == 1

            # __delitem__
            store["a"] = ts
            store["b"] = df
            del store["a"]
            del store["b"]
            assert len(store) == 0

    def test_invalid_terms(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            with catch_warnings(record=True):

                df = tm.makeTimeDataFrame()
                df["string"] = "foo"
                df.loc[df.index[0:4], "string"] = "bar"

                store.put("df", df, format="table")

                # some invalid terms
                with pytest.raises(TypeError):
                    Term()

                # more invalid
                with pytest.raises(ValueError):
                    store.select("df", "df.index[3]")

                with pytest.raises(SyntaxError):
                    store.select("df", "index>")

        # from the docs
        with ensure_clean_path(setup_path) as path:
            dfq = DataFrame(
                np.random.randn(10, 4),
                columns=list("ABCD"),
                index=date_range("20130101", periods=10),
            )
            dfq.to_hdf(path, "dfq", format="table", data_columns=True)

            # check ok
            read_hdf(
                path, "dfq", where="index>Timestamp('20130104') & columns=['A', 'B']"
            )
            read_hdf(path, "dfq", where="A>0 or C>0")

        # catch the invalid reference
        with ensure_clean_path(setup_path) as path:
            dfq = DataFrame(
                np.random.randn(10, 4),
                columns=list("ABCD"),
                index=date_range("20130101", periods=10),
            )
            dfq.to_hdf(path, "dfq", format="table")

            with pytest.raises(ValueError):
                read_hdf(path, "dfq", where="A>0 or C>0")

    def test_same_name_scoping(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            import pandas as pd

            df = DataFrame(
                np.random.randn(20, 2), index=pd.date_range("20130101", periods=20)
            )
            store.put("df", df, format="table")
            expected = df[df.index > pd.Timestamp("20130105")]

            import datetime  # noqa

            result = store.select("df", "index>datetime.datetime(2013,1,5)")
            tm.assert_frame_equal(result, expected)

            from datetime import datetime  # noqa

            # technically an error, but allow it
            result = store.select("df", "index>datetime.datetime(2013,1,5)")
            tm.assert_frame_equal(result, expected)

            result = store.select("df", "index>datetime(2013,1,5)")
            tm.assert_frame_equal(result, expected)

    def test_series(self, setup_path):

        s = tm.makeStringSeries()
        self._check_roundtrip(s, tm.assert_series_equal, path=setup_path)

        ts = tm.makeTimeSeries()
        self._check_roundtrip(ts, tm.assert_series_equal, path=setup_path)

        ts2 = Series(ts.index, Index(ts.index, dtype=object))
        self._check_roundtrip(ts2, tm.assert_series_equal, path=setup_path)

        ts3 = Series(ts.values, Index(np.asarray(ts.index, dtype=object), dtype=object))
        self._check_roundtrip(
            ts3, tm.assert_series_equal, path=setup_path, check_index_type=False
        )

    def test_float_index(self, setup_path):

        # GH #454
        index = np.random.randn(10)
        s = Series(np.random.randn(10), index=index)
        self._check_roundtrip(s, tm.assert_series_equal, path=setup_path)

    @td.xfail_non_writeable
    def test_tuple_index(self, setup_path):

        # GH #492
        col = np.arange(10)
        idx = [(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)]
        data = np.random.randn(30).reshape((3, 10))
        DF = DataFrame(data, index=idx, columns=col)

        with catch_warnings(record=True):
            simplefilter("ignore", pd.errors.PerformanceWarning)
            self._check_roundtrip(DF, tm.assert_frame_equal, path=setup_path)

    @td.xfail_non_writeable
    @pytest.mark.filterwarnings("ignore::pandas.errors.PerformanceWarning")
    def test_index_types(self, setup_path):

        with catch_warnings(record=True):
            values = np.random.randn(2)

            func = lambda l, r: tm.assert_series_equal(l, r, check_index_type=True)

        with catch_warnings(record=True):
            ser = Series(values, [0, "y"])
            self._check_roundtrip(ser, func, path=setup_path)

        with catch_warnings(record=True):
            ser = Series(values, [datetime.datetime.today(), 0])
            self._check_roundtrip(ser, func, path=setup_path)

        with catch_warnings(record=True):
            ser = Series(values, ["y", 0])
            self._check_roundtrip(ser, func, path=setup_path)

        with catch_warnings(record=True):
            ser = Series(values, [datetime.date.today(), "a"])
            self._check_roundtrip(ser, func, path=setup_path)

        with catch_warnings(record=True):

            ser = Series(values, [0, "y"])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(values, [datetime.datetime.today(), 0])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(values, ["y", 0])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(values, [datetime.date.today(), "a"])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(values, [1.23, "b"])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(values, [1, 1.53])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(values, [1, 5])
            self._check_roundtrip(ser, func, path=setup_path)

            ser = Series(
                values, [datetime.datetime(2012, 1, 1), datetime.datetime(2012, 1, 2)]
            )
            self._check_roundtrip(ser, func, path=setup_path)

    def test_timeseries_preepoch(self, setup_path):

        dr = bdate_range("1/1/1940", "1/1/1960")
        ts = Series(np.random.randn(len(dr)), index=dr)
        try:
            self._check_roundtrip(ts, tm.assert_series_equal, path=setup_path)
        except OverflowError:
            pytest.skip("known failer on some windows platforms")

    @td.xfail_non_writeable
    @pytest.mark.parametrize(
        "compression", [False, pytest.param(True, marks=td.skip_if_windows_python_3)]
    )
    def test_frame(self, compression, setup_path):

        df = tm.makeDataFrame()

        # put in some random NAs
        df.values[0, 0] = np.nan
        df.values[5, 3] = np.nan

        self._check_roundtrip_table(
            df, tm.assert_frame_equal, path=setup_path, compression=compression
        )
        self._check_roundtrip(
            df, tm.assert_frame_equal, path=setup_path, compression=compression
        )

        tdf = tm.makeTimeDataFrame()
        self._check_roundtrip(
            tdf, tm.assert_frame_equal, path=setup_path, compression=compression
        )

        with ensure_clean_store(setup_path) as store:
            # not consolidated
            df["foo"] = np.random.randn(len(df))
            store["df"] = df
            recons = store["df"]
            assert recons._data.is_consolidated()

        # empty
        self._check_roundtrip(df[:0], tm.assert_frame_equal, path=setup_path)

    @td.xfail_non_writeable
    def test_empty_series_frame(self, setup_path):
        s0 = Series(dtype=object)
        s1 = Series(name="myseries", dtype=object)
        df0 = DataFrame()
        df1 = DataFrame(index=["a", "b", "c"])
        df2 = DataFrame(columns=["d", "e", "f"])

        self._check_roundtrip(s0, tm.assert_series_equal, path=setup_path)
        self._check_roundtrip(s1, tm.assert_series_equal, path=setup_path)
        self._check_roundtrip(df0, tm.assert_frame_equal, path=setup_path)
        self._check_roundtrip(df1, tm.assert_frame_equal, path=setup_path)
        self._check_roundtrip(df2, tm.assert_frame_equal, path=setup_path)

    @td.xfail_non_writeable
    @pytest.mark.parametrize(
        "dtype", [np.int64, np.float64, np.object, "m8[ns]", "M8[ns]"]
    )
    def test_empty_series(self, dtype, setup_path):
        s = Series(dtype=dtype)
        self._check_roundtrip(s, tm.assert_series_equal, path=setup_path)

    def test_can_serialize_dates(self, setup_path):

        rng = [x.date() for x in bdate_range("1/1/2000", "1/30/2000")]
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        self._check_roundtrip(frame, tm.assert_frame_equal, path=setup_path)

    def test_store_hierarchical(self, setup_path):
        index = MultiIndex(
            levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
            codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
            names=["foo", "bar"],
        )
        frame = DataFrame(np.random.randn(10, 3), index=index, columns=["A", "B", "C"])

        self._check_roundtrip(frame, tm.assert_frame_equal, path=setup_path)
        self._check_roundtrip(frame.T, tm.assert_frame_equal, path=setup_path)
        self._check_roundtrip(frame["A"], tm.assert_series_equal, path=setup_path)

        # check that the names are stored
        with ensure_clean_store(setup_path) as store:
            store["frame"] = frame
            recons = store["frame"]
            tm.assert_frame_equal(recons, frame)

    def test_store_index_name(self, setup_path):
        df = tm.makeDataFrame()
        df.index.name = "foo"

        with ensure_clean_store(setup_path) as store:
            store["frame"] = df
            recons = store["frame"]
            tm.assert_frame_equal(recons, df)

    def test_store_index_name_with_tz(self, setup_path):
        # GH 13884
        df = pd.DataFrame({"A": [1, 2]})
        df.index = pd.DatetimeIndex([1234567890123456787, 1234567890123456788])
        df.index = df.index.tz_localize("UTC")
        df.index.name = "foo"

        with ensure_clean_store(setup_path) as store:
            store.put("frame", df, format="table")
            recons = store["frame"]
            tm.assert_frame_equal(recons, df)

    @pytest.mark.parametrize("table_format", ["table", "fixed"])
    def test_store_index_name_numpy_str(self, table_format, setup_path):
        # GH #13492
        idx = pd.Index(
            pd.to_datetime([datetime.date(2000, 1, 1), datetime.date(2000, 1, 2)]),
            name="cols\u05d2",
        )
        idx1 = pd.Index(
            pd.to_datetime([datetime.date(2010, 1, 1), datetime.date(2010, 1, 2)]),
            name="rows\u05d0",
        )
        df = pd.DataFrame(np.arange(4).reshape(2, 2), columns=idx, index=idx1)

        # This used to fail, returning numpy strings instead of python strings.
        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", format=table_format)
            df2 = read_hdf(path, "df")

            tm.assert_frame_equal(df, df2, check_names=True)

            assert type(df2.index.name) == str
            assert type(df2.columns.name) == str

    def test_store_series_name(self, setup_path):
        df = tm.makeDataFrame()
        series = df["A"]

        with ensure_clean_store(setup_path) as store:
            store["series"] = series
            recons = store["series"]
            tm.assert_series_equal(recons, series)

    @td.xfail_non_writeable
    @pytest.mark.parametrize(
        "compression", [False, pytest.param(True, marks=td.skip_if_windows_python_3)]
    )
    def test_store_mixed(self, compression, setup_path):
        def _make_one():
            df = tm.makeDataFrame()
            df["obj1"] = "foo"
            df["obj2"] = "bar"
            df["bool1"] = df["A"] > 0
            df["bool2"] = df["B"] > 0
            df["int1"] = 1
            df["int2"] = 2
            return df._consolidate()

        df1 = _make_one()
        df2 = _make_one()

        self._check_roundtrip(df1, tm.assert_frame_equal, path=setup_path)
        self._check_roundtrip(df2, tm.assert_frame_equal, path=setup_path)

        with ensure_clean_store(setup_path) as store:
            store["obj"] = df1
            tm.assert_frame_equal(store["obj"], df1)
            store["obj"] = df2
            tm.assert_frame_equal(store["obj"], df2)

        # check that can store Series of all of these types
        self._check_roundtrip(
            df1["obj1"],
            tm.assert_series_equal,
            path=setup_path,
            compression=compression,
        )
        self._check_roundtrip(
            df1["bool1"],
            tm.assert_series_equal,
            path=setup_path,
            compression=compression,
        )
        self._check_roundtrip(
            df1["int1"],
            tm.assert_series_equal,
            path=setup_path,
            compression=compression,
        )

    @pytest.mark.filterwarnings(
        "ignore:\\nduplicate:pandas.io.pytables.DuplicateWarning"
    )
    def test_select_with_dups(self, setup_path):

        # single dtypes
        df = DataFrame(np.random.randn(10, 4), columns=["A", "A", "B", "B"])
        df.index = date_range("20130101 9:30", periods=10, freq="T")

        with ensure_clean_store(setup_path) as store:
            store.append("df", df)

            result = store.select("df")
            expected = df
            tm.assert_frame_equal(result, expected, by_blocks=True)

            result = store.select("df", columns=df.columns)
            expected = df
            tm.assert_frame_equal(result, expected, by_blocks=True)

            result = store.select("df", columns=["A"])
            expected = df.loc[:, ["A"]]
            tm.assert_frame_equal(result, expected)

        # dups across dtypes
        df = concat(
            [
                DataFrame(np.random.randn(10, 4), columns=["A", "A", "B", "B"]),
                DataFrame(
                    np.random.randint(0, 10, size=20).reshape(10, 2), columns=["A", "C"]
                ),
            ],
            axis=1,
        )
        df.index = date_range("20130101 9:30", periods=10, freq="T")

        with ensure_clean_store(setup_path) as store:
            store.append("df", df)

            result = store.select("df")
            expected = df
            tm.assert_frame_equal(result, expected, by_blocks=True)

            result = store.select("df", columns=df.columns)
            expected = df
            tm.assert_frame_equal(result, expected, by_blocks=True)

            expected = df.loc[:, ["A"]]
            result = store.select("df", columns=["A"])
            tm.assert_frame_equal(result, expected, by_blocks=True)

            expected = df.loc[:, ["B", "A"]]
            result = store.select("df", columns=["B", "A"])
            tm.assert_frame_equal(result, expected, by_blocks=True)

        # duplicates on both index and columns
        with ensure_clean_store(setup_path) as store:
            store.append("df", df)
            store.append("df", df)

            expected = df.loc[:, ["B", "A"]]
            expected = concat([expected, expected])
            result = store.select("df", columns=["B", "A"])
            tm.assert_frame_equal(result, expected, by_blocks=True)

    def test_overwrite_node(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            store["a"] = tm.makeTimeDataFrame()
            ts = tm.makeTimeSeries()
            store["a"] = ts

            tm.assert_series_equal(store["a"], ts)

    def test_select(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            with catch_warnings(record=True):

                # select with columns=
                df = tm.makeTimeDataFrame()
                _maybe_remove(store, "df")
                store.append("df", df)
                result = store.select("df", columns=["A", "B"])
                expected = df.reindex(columns=["A", "B"])
                tm.assert_frame_equal(expected, result)

                # equivalently
                result = store.select("df", [("columns=['A', 'B']")])
                expected = df.reindex(columns=["A", "B"])
                tm.assert_frame_equal(expected, result)

                # with a data column
                _maybe_remove(store, "df")
                store.append("df", df, data_columns=["A"])
                result = store.select("df", ["A > 0"], columns=["A", "B"])
                expected = df[df.A > 0].reindex(columns=["A", "B"])
                tm.assert_frame_equal(expected, result)

                # all a data columns
                _maybe_remove(store, "df")
                store.append("df", df, data_columns=True)
                result = store.select("df", ["A > 0"], columns=["A", "B"])
                expected = df[df.A > 0].reindex(columns=["A", "B"])
                tm.assert_frame_equal(expected, result)

                # with a data column, but different columns
                _maybe_remove(store, "df")
                store.append("df", df, data_columns=["A"])
                result = store.select("df", ["A > 0"], columns=["C", "D"])
                expected = df[df.A > 0].reindex(columns=["C", "D"])
                tm.assert_frame_equal(expected, result)

    def test_select_dtypes(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            # with a Timestamp data column (GH #2637)
            df = DataFrame(
                dict(ts=bdate_range("2012-01-01", periods=300), A=np.random.randn(300))
            )
            _maybe_remove(store, "df")
            store.append("df", df, data_columns=["ts", "A"])

            result = store.select("df", "ts>=Timestamp('2012-02-01')")
            expected = df[df.ts >= Timestamp("2012-02-01")]
            tm.assert_frame_equal(expected, result)

            # bool columns (GH #2849)
            df = DataFrame(np.random.randn(5, 2), columns=["A", "B"])
            df["object"] = "foo"
            df.loc[4:5, "object"] = "bar"
            df["boolv"] = df["A"] > 0
            _maybe_remove(store, "df")
            store.append("df", df, data_columns=True)

            expected = df[df.boolv == True].reindex(columns=["A", "boolv"])  # noqa
            for v in [True, "true", 1]:
                result = store.select("df", f"boolv == {v}", columns=["A", "boolv"])
                tm.assert_frame_equal(expected, result)

            expected = df[df.boolv == False].reindex(columns=["A", "boolv"])  # noqa
            for v in [False, "false", 0]:
                result = store.select("df", f"boolv == {v}", columns=["A", "boolv"])
                tm.assert_frame_equal(expected, result)

            # integer index
            df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
            _maybe_remove(store, "df_int")
            store.append("df_int", df)
            result = store.select("df_int", "index<10 and columns=['A']")
            expected = df.reindex(index=list(df.index)[0:10], columns=["A"])
            tm.assert_frame_equal(expected, result)

            # float index
            df = DataFrame(
                dict(
                    A=np.random.rand(20),
                    B=np.random.rand(20),
                    index=np.arange(20, dtype="f8"),
                )
            )
            _maybe_remove(store, "df_float")
            store.append("df_float", df)
            result = store.select("df_float", "index<10.0 and columns=['A']")
            expected = df.reindex(index=list(df.index)[0:10], columns=["A"])
            tm.assert_frame_equal(expected, result)

        with ensure_clean_store(setup_path) as store:

            # floats w/o NaN
            df = DataFrame(dict(cols=range(11), values=range(11)), dtype="float64")
            df["cols"] = (df["cols"] + 10).apply(str)

            store.append("df1", df, data_columns=True)
            result = store.select("df1", where="values>2.0")
            expected = df[df["values"] > 2.0]
            tm.assert_frame_equal(expected, result)

            # floats with NaN
            df.iloc[0] = np.nan
            expected = df[df["values"] > 2.0]

            store.append("df2", df, data_columns=True, index=False)
            result = store.select("df2", where="values>2.0")
            tm.assert_frame_equal(expected, result)

            # https://github.com/PyTables/PyTables/issues/282
            # bug in selection when 0th row has a np.nan and an index
            # store.append('df3',df,data_columns=True)
            # result = store.select(
            #    'df3', where='values>2.0')
            # tm.assert_frame_equal(expected, result)

            # not in first position float with NaN ok too
            df = DataFrame(dict(cols=range(11), values=range(11)), dtype="float64")
            df["cols"] = (df["cols"] + 10).apply(str)

            df.iloc[1] = np.nan
            expected = df[df["values"] > 2.0]

            store.append("df4", df, data_columns=True)
            result = store.select("df4", where="values>2.0")
            tm.assert_frame_equal(expected, result)

        # test selection with comparison against numpy scalar
        # GH 11283
        with ensure_clean_store(setup_path) as store:
            df = tm.makeDataFrame()

            expected = df[df["A"] > 0]

            store.append("df", df, data_columns=True)
            np_zero = np.float64(0)  # noqa
            result = store.select("df", where=["A>np_zero"])
            tm.assert_frame_equal(expected, result)

    def test_select_with_many_inputs(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            df = DataFrame(
                dict(
                    ts=bdate_range("2012-01-01", periods=300),
                    A=np.random.randn(300),
                    B=range(300),
                    users=["a"] * 50
                    + ["b"] * 50
                    + ["c"] * 100
                    + [f"a{i:03d}" for i in range(100)],
                )
            )
            _maybe_remove(store, "df")
            store.append("df", df, data_columns=["ts", "A", "B", "users"])

            # regular select
            result = store.select("df", "ts>=Timestamp('2012-02-01')")
            expected = df[df.ts >= Timestamp("2012-02-01")]
            tm.assert_frame_equal(expected, result)

            # small selector
            result = store.select(
                "df", "ts>=Timestamp('2012-02-01') & users=['a','b','c']"
            )
            expected = df[
                (df.ts >= Timestamp("2012-02-01")) & df.users.isin(["a", "b", "c"])
            ]
            tm.assert_frame_equal(expected, result)

            # big selector along the columns
            selector = ["a", "b", "c"] + [f"a{i:03d}" for i in range(60)]
            result = store.select(
                "df", "ts>=Timestamp('2012-02-01') and users=selector"
            )
            expected = df[(df.ts >= Timestamp("2012-02-01")) & df.users.isin(selector)]
            tm.assert_frame_equal(expected, result)

            selector = range(100, 200)
            result = store.select("df", "B=selector")
            expected = df[df.B.isin(selector)]
            tm.assert_frame_equal(expected, result)
            assert len(result) == 100

            # big selector along the index
            selector = Index(df.ts[0:100].values)
            result = store.select("df", "ts=selector")
            expected = df[df.ts.isin(selector.values)]
            tm.assert_frame_equal(expected, result)
            assert len(result) == 100

    def test_select_iterator(self, setup_path):

        # single table
        with ensure_clean_store(setup_path) as store:

            df = tm.makeTimeDataFrame(500)
            _maybe_remove(store, "df")
            store.append("df", df)

            expected = store.select("df")

            results = list(store.select("df", iterator=True))
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            results = list(store.select("df", chunksize=100))
            assert len(results) == 5
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            results = list(store.select("df", chunksize=150))
            result = concat(results)
            tm.assert_frame_equal(result, expected)

        with ensure_clean_path(setup_path) as path:

            df = tm.makeTimeDataFrame(500)
            df.to_hdf(path, "df_non_table")

            with pytest.raises(TypeError):
                read_hdf(path, "df_non_table", chunksize=100)

            with pytest.raises(TypeError):
                read_hdf(path, "df_non_table", iterator=True)

        with ensure_clean_path(setup_path) as path:

            df = tm.makeTimeDataFrame(500)
            df.to_hdf(path, "df", format="table")

            results = list(read_hdf(path, "df", chunksize=100))
            result = concat(results)

            assert len(results) == 5
            tm.assert_frame_equal(result, df)
            tm.assert_frame_equal(result, read_hdf(path, "df"))

        # multiple

        with ensure_clean_store(setup_path) as store:

            df1 = tm.makeTimeDataFrame(500)
            store.append("df1", df1, data_columns=True)
            df2 = tm.makeTimeDataFrame(500).rename(columns="{}_2".format)
            df2["foo"] = "bar"
            store.append("df2", df2)

            df = concat([df1, df2], axis=1)

            # full selection
            expected = store.select_as_multiple(["df1", "df2"], selector="df1")
            results = list(
                store.select_as_multiple(["df1", "df2"], selector="df1", chunksize=150)
            )
            result = concat(results)
            tm.assert_frame_equal(expected, result)

    def test_select_iterator_complete_8014(self, setup_path):

        # GH 8014
        # using iterator and where clause
        chunksize = 1e4

        # no iterator
        with ensure_clean_store(setup_path) as store:

            expected = tm.makeTimeDataFrame(100064, "S")
            _maybe_remove(store, "df")
            store.append("df", expected)

            beg_dt = expected.index[0]
            end_dt = expected.index[-1]

            # select w/o iteration and no where clause works
            result = store.select("df")
            tm.assert_frame_equal(expected, result)

            # select w/o iterator and where clause, single term, begin
            # of range, works
            where = f"index >= '{beg_dt}'"
            result = store.select("df", where=where)
            tm.assert_frame_equal(expected, result)

            # select w/o iterator and where clause, single term, end
            # of range, works
            where = f"index <= '{end_dt}'"
            result = store.select("df", where=where)
            tm.assert_frame_equal(expected, result)

            # select w/o iterator and where clause, inclusive range,
            # works
            where = f"index >= '{beg_dt}' & index <= '{end_dt}'"
            result = store.select("df", where=where)
            tm.assert_frame_equal(expected, result)

        # with iterator, full range
        with ensure_clean_store(setup_path) as store:

            expected = tm.makeTimeDataFrame(100064, "S")
            _maybe_remove(store, "df")
            store.append("df", expected)

            beg_dt = expected.index[0]
            end_dt = expected.index[-1]

            # select w/iterator and no where clause works
            results = list(store.select("df", chunksize=chunksize))
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # select w/iterator and where clause, single term, begin of range
            where = f"index >= '{beg_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # select w/iterator and where clause, single term, end of range
            where = f"index <= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # select w/iterator and where clause, inclusive range
            where = f"index >= '{beg_dt}' & index <= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            tm.assert_frame_equal(expected, result)

    def test_select_iterator_non_complete_8014(self, setup_path):

        # GH 8014
        # using iterator and where clause
        chunksize = 1e4

        # with iterator, non complete range
        with ensure_clean_store(setup_path) as store:

            expected = tm.makeTimeDataFrame(100064, "S")
            _maybe_remove(store, "df")
            store.append("df", expected)

            beg_dt = expected.index[1]
            end_dt = expected.index[-2]

            # select w/iterator and where clause, single term, begin of range
            where = f"index >= '{beg_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            rexpected = expected[expected.index >= beg_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, single term, end of range
            where = f"index <= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            rexpected = expected[expected.index <= end_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, inclusive range
            where = f"index >= '{beg_dt}' & index <= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            rexpected = expected[
                (expected.index >= beg_dt) & (expected.index <= end_dt)
            ]
            tm.assert_frame_equal(rexpected, result)

        # with iterator, empty where
        with ensure_clean_store(setup_path) as store:

            expected = tm.makeTimeDataFrame(100064, "S")
            _maybe_remove(store, "df")
            store.append("df", expected)

            end_dt = expected.index[-1]

            # select w/iterator and where clause, single term, begin of range
            where = f"index > '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            assert 0 == len(results)

    def test_select_iterator_many_empty_frames(self, setup_path):

        # GH 8014
        # using iterator and where clause can return many empty
        # frames.
        chunksize = int(1e4)

        # with iterator, range limited to the first chunk
        with ensure_clean_store(setup_path) as store:

            expected = tm.makeTimeDataFrame(100000, "S")
            _maybe_remove(store, "df")
            store.append("df", expected)

            beg_dt = expected.index[0]
            end_dt = expected.index[chunksize - 1]

            # select w/iterator and where clause, single term, begin of range
            where = f"index >= '{beg_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))
            result = concat(results)
            rexpected = expected[expected.index >= beg_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, single term, end of range
            where = f"index <= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))

            assert len(results) == 1
            result = concat(results)
            rexpected = expected[expected.index <= end_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, inclusive range
            where = f"index >= '{beg_dt}' & index <= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))

            # should be 1, is 10
            assert len(results) == 1
            result = concat(results)
            rexpected = expected[
                (expected.index >= beg_dt) & (expected.index <= end_dt)
            ]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause which selects
            # *nothing*.
            #
            # To be consistent with Python idiom I suggest this should
            # return [] e.g. `for e in []: print True` never prints
            # True.

            where = f"index <= '{beg_dt}' & index >= '{end_dt}'"
            results = list(store.select("df", where=where, chunksize=chunksize))

            # should be []
            assert len(results) == 0

    @pytest.mark.filterwarnings(
        "ignore:\\nthe :pandas.io.pytables.AttributeConflictWarning"
    )
    def test_retain_index_attributes(self, setup_path):

        # GH 3499, losing frequency info on index recreation
        df = DataFrame(
            dict(A=Series(range(3), index=date_range("2000-1-1", periods=3, freq="H")))
        )

        with ensure_clean_store(setup_path) as store:
            _maybe_remove(store, "data")
            store.put("data", df, format="table")

            result = store.get("data")
            tm.assert_frame_equal(df, result)

            for attr in ["freq", "tz", "name"]:
                for idx in ["index", "columns"]:
                    assert getattr(getattr(df, idx), attr, None) == getattr(
                        getattr(result, idx), attr, None
                    )

            # try to append a table with a different frequency
            with catch_warnings(record=True):
                df2 = DataFrame(
                    dict(
                        A=Series(
                            range(3), index=date_range("2002-1-1", periods=3, freq="D")
                        )
                    )
                )
                store.append("data", df2)

            assert store.get_storer("data").info["index"]["freq"] is None

            # this is ok
            _maybe_remove(store, "df2")
            df2 = DataFrame(
                dict(
                    A=Series(
                        range(3),
                        index=[
                            Timestamp("20010101"),
                            Timestamp("20010102"),
                            Timestamp("20020101"),
                        ],
                    )
                )
            )
            store.append("df2", df2)
            df3 = DataFrame(
                dict(
                    A=Series(
                        range(3), index=date_range("2002-1-1", periods=3, freq="D")
                    )
                )
            )
            store.append("df2", df3)

    @pytest.mark.filterwarnings(
        "ignore:\\nthe :pandas.io.pytables.AttributeConflictWarning"
    )
    def test_retain_index_attributes2(self, setup_path):
        with ensure_clean_path(setup_path) as path:

            with catch_warnings(record=True):

                df = DataFrame(
                    dict(
                        A=Series(
                            range(3), index=date_range("2000-1-1", periods=3, freq="H")
                        )
                    )
                )
                df.to_hdf(path, "data", mode="w", append=True)
                df2 = DataFrame(
                    dict(
                        A=Series(
                            range(3), index=date_range("2002-1-1", periods=3, freq="D")
                        )
                    )
                )
                df2.to_hdf(path, "data", append=True)

                idx = date_range("2000-1-1", periods=3, freq="H")
                idx.name = "foo"
                df = DataFrame(dict(A=Series(range(3), index=idx)))
                df.to_hdf(path, "data", mode="w", append=True)

            assert read_hdf(path, "data").index.name == "foo"

            with catch_warnings(record=True):

                idx2 = date_range("2001-1-1", periods=3, freq="H")
                idx2.name = "bar"
                df2 = DataFrame(dict(A=Series(range(3), index=idx2)))
                df2.to_hdf(path, "data", append=True)

            assert read_hdf(path, "data").index.name is None

    def test_frame_select(self, setup_path):

        df = tm.makeTimeDataFrame()

        with ensure_clean_store(setup_path) as store:
            store.put("frame", df, format="table")
            date = df.index[len(df) // 2]

            crit1 = Term("index>=date")
            assert crit1.env.scope["date"] == date

            crit2 = "columns=['A', 'D']"
            crit3 = "columns=A"

            result = store.select("frame", [crit1, crit2])
            expected = df.loc[date:, ["A", "D"]]
            tm.assert_frame_equal(result, expected)

            result = store.select("frame", [crit3])
            expected = df.loc[:, ["A"]]
            tm.assert_frame_equal(result, expected)

            # invalid terms
            df = tm.makeTimeDataFrame()
            store.append("df_time", df)
            with pytest.raises(ValueError):
                store.select("df_time", "index>0")

            # can't select if not written as table
            # store['frame'] = df
            # with pytest.raises(ValueError):
            #     store.select('frame', [crit1, crit2])

    def test_frame_select_complex(self, setup_path):
        # select via complex criteria

        df = tm.makeTimeDataFrame()
        df["string"] = "foo"
        df.loc[df.index[0:4], "string"] = "bar"

        with ensure_clean_store(setup_path) as store:
            store.put("df", df, format="table", data_columns=["string"])

            # empty
            result = store.select("df", 'index>df.index[3] & string="bar"')
            expected = df.loc[(df.index > df.index[3]) & (df.string == "bar")]
            tm.assert_frame_equal(result, expected)

            result = store.select("df", 'index>df.index[3] & string="foo"')
            expected = df.loc[(df.index > df.index[3]) & (df.string == "foo")]
            tm.assert_frame_equal(result, expected)

            # or
            result = store.select("df", 'index>df.index[3] | string="bar"')
            expected = df.loc[(df.index > df.index[3]) | (df.string == "bar")]
            tm.assert_frame_equal(result, expected)

            result = store.select(
                "df", '(index>df.index[3] & index<=df.index[6]) | string="bar"'
            )
            expected = df.loc[
                ((df.index > df.index[3]) & (df.index <= df.index[6]))
                | (df.string == "bar")
            ]
            tm.assert_frame_equal(result, expected)

            # invert
            result = store.select("df", 'string!="bar"')
            expected = df.loc[df.string != "bar"]
            tm.assert_frame_equal(result, expected)

            # invert not implemented in numexpr :(
            with pytest.raises(NotImplementedError):
                store.select("df", '~(string="bar")')

            # invert ok for filters
            result = store.select("df", "~(columns=['A','B'])")
            expected = df.loc[:, df.columns.difference(["A", "B"])]
            tm.assert_frame_equal(result, expected)

            # in
            result = store.select("df", "index>df.index[3] & columns in ['A','B']")
            expected = df.loc[df.index > df.index[3]].reindex(columns=["A", "B"])
            tm.assert_frame_equal(result, expected)

    def test_frame_select_complex2(self, setup_path):

        with ensure_clean_path(["parms.hdf", "hist.hdf"]) as paths:

            pp, hh = paths

            # use non-trivial selection criteria
            parms = DataFrame({"A": [1, 1, 2, 2, 3]})
            parms.to_hdf(pp, "df", mode="w", format="table", data_columns=["A"])

            selection = read_hdf(pp, "df", where="A=[2,3]")
            hist = DataFrame(
                np.random.randn(25, 1),
                columns=["data"],
                index=MultiIndex.from_tuples(
                    [(i, j) for i in range(5) for j in range(5)], names=["l1", "l2"]
                ),
            )

            hist.to_hdf(hh, "df", mode="w", format="table")

            expected = read_hdf(hh, "df", where="l1=[2, 3, 4]")

            # scope with list like
            l = selection.index.tolist()  # noqa
            store = HDFStore(hh)
            result = store.select("df", where="l1=l")
            tm.assert_frame_equal(result, expected)
            store.close()

            result = read_hdf(hh, "df", where="l1=l")
            tm.assert_frame_equal(result, expected)

            # index
            index = selection.index  # noqa
            result = read_hdf(hh, "df", where="l1=index")
            tm.assert_frame_equal(result, expected)

            result = read_hdf(hh, "df", where="l1=selection.index")
            tm.assert_frame_equal(result, expected)

            result = read_hdf(hh, "df", where="l1=selection.index.tolist()")
            tm.assert_frame_equal(result, expected)

            result = read_hdf(hh, "df", where="l1=list(selection.index)")
            tm.assert_frame_equal(result, expected)

            # scope with index
            store = HDFStore(hh)

            result = store.select("df", where="l1=index")
            tm.assert_frame_equal(result, expected)

            result = store.select("df", where="l1=selection.index")
            tm.assert_frame_equal(result, expected)

            result = store.select("df", where="l1=selection.index.tolist()")
            tm.assert_frame_equal(result, expected)

            result = store.select("df", where="l1=list(selection.index)")
            tm.assert_frame_equal(result, expected)

            store.close()

    def test_invalid_filtering(self, setup_path):

        # can't use more than one filter (atm)

        df = tm.makeTimeDataFrame()

        with ensure_clean_store(setup_path) as store:
            store.put("df", df, format="table")

            # not implemented
            with pytest.raises(NotImplementedError):
                store.select("df", "columns=['A'] | columns=['B']")

            # in theory we could deal with this
            with pytest.raises(NotImplementedError):
                store.select("df", "columns=['A','B'] & columns=['C']")

    def test_string_select(self, setup_path):
        # GH 2973
        with ensure_clean_store(setup_path) as store:

            df = tm.makeTimeDataFrame()

            # test string ==/!=
            df["x"] = "none"
            df.loc[df.index[2:7], "x"] = ""

            store.append("df", df, data_columns=["x"])

            result = store.select("df", "x=none")
            expected = df[df.x == "none"]
            tm.assert_frame_equal(result, expected)

            result = store.select("df", "x!=none")
            expected = df[df.x != "none"]
            tm.assert_frame_equal(result, expected)

            df2 = df.copy()
            df2.loc[df2.x == "", "x"] = np.nan

            store.append("df2", df2, data_columns=["x"])
            result = store.select("df2", "x!=none")
            expected = df2[isna(df2.x)]
            tm.assert_frame_equal(result, expected)

            # int ==/!=
            df["int"] = 1
            df.loc[df.index[2:7], "int"] = 2

            store.append("df3", df, data_columns=["int"])

            result = store.select("df3", "int=2")
            expected = df[df.int == 2]
            tm.assert_frame_equal(result, expected)

            result = store.select("df3", "int!=2")
            expected = df[df.int != 2]
            tm.assert_frame_equal(result, expected)

    def test_read_column(self, setup_path):

        df = tm.makeTimeDataFrame()

        with ensure_clean_store(setup_path) as store:
            _maybe_remove(store, "df")

            # GH 17912
            # HDFStore.select_column should raise a KeyError
            # exception if the key is not a valid store
            with pytest.raises(KeyError, match="No object named df in the file"):
                store.select_column("df", "index")

            store.append("df", df)
            # error
            with pytest.raises(
                KeyError, match=re.escape("'column [foo] not found in the table'")
            ):
                store.select_column("df", "foo")

            with pytest.raises(Exception):
                store.select_column("df", "index", where=["index>5"])

            # valid
            result = store.select_column("df", "index")
            tm.assert_almost_equal(result.values, Series(df.index).values)
            assert isinstance(result, Series)

            # not a data indexable column
            with pytest.raises(ValueError):
                store.select_column("df", "values_block_0")

            # a data column
            df2 = df.copy()
            df2["string"] = "foo"
            store.append("df2", df2, data_columns=["string"])
            result = store.select_column("df2", "string")
            tm.assert_almost_equal(result.values, df2["string"].values)

            # a data column with NaNs, result excludes the NaNs
            df3 = df.copy()
            df3["string"] = "foo"
            df3.loc[df3.index[4:6], "string"] = np.nan
            store.append("df3", df3, data_columns=["string"])
            result = store.select_column("df3", "string")
            tm.assert_almost_equal(result.values, df3["string"].values)

            # start/stop
            result = store.select_column("df3", "string", start=2)
            tm.assert_almost_equal(result.values, df3["string"].values[2:])

            result = store.select_column("df3", "string", start=-2)
            tm.assert_almost_equal(result.values, df3["string"].values[-2:])

            result = store.select_column("df3", "string", stop=2)
            tm.assert_almost_equal(result.values, df3["string"].values[:2])

            result = store.select_column("df3", "string", stop=-2)
            tm.assert_almost_equal(result.values, df3["string"].values[:-2])

            result = store.select_column("df3", "string", start=2, stop=-2)
            tm.assert_almost_equal(result.values, df3["string"].values[2:-2])

            result = store.select_column("df3", "string", start=-2, stop=2)
            tm.assert_almost_equal(result.values, df3["string"].values[-2:2])

            # GH 10392 - make sure column name is preserved
            df4 = DataFrame({"A": np.random.randn(10), "B": "foo"})
            store.append("df4", df4, data_columns=True)
            expected = df4["B"]
            result = store.select_column("df4", "B")
            tm.assert_series_equal(result, expected)

    def test_coordinates(self, setup_path):
        df = tm.makeTimeDataFrame()

        with ensure_clean_store(setup_path) as store:

            _maybe_remove(store, "df")
            store.append("df", df)

            # all
            c = store.select_as_coordinates("df")
            assert (c.values == np.arange(len(df.index))).all()

            # get coordinates back & test vs frame
            _maybe_remove(store, "df")

            df = DataFrame(dict(A=range(5), B=range(5)))
            store.append("df", df)
            c = store.select_as_coordinates("df", ["index<3"])
            assert (c.values == np.arange(3)).all()
            result = store.select("df", where=c)
            expected = df.loc[0:2, :]
            tm.assert_frame_equal(result, expected)

            c = store.select_as_coordinates("df", ["index>=3", "index<=4"])
            assert (c.values == np.arange(2) + 3).all()
            result = store.select("df", where=c)
            expected = df.loc[3:4, :]
            tm.assert_frame_equal(result, expected)
            assert isinstance(c, Index)

            # multiple tables
            _maybe_remove(store, "df1")
            _maybe_remove(store, "df2")
            df1 = tm.makeTimeDataFrame()
            df2 = tm.makeTimeDataFrame().rename(columns="{}_2".format)
            store.append("df1", df1, data_columns=["A", "B"])
            store.append("df2", df2)

            c = store.select_as_coordinates("df1", ["A>0", "B>0"])
            df1_result = store.select("df1", c)
            df2_result = store.select("df2", c)
            result = concat([df1_result, df2_result], axis=1)

            expected = concat([df1, df2], axis=1)
            expected = expected[(expected.A > 0) & (expected.B > 0)]
            tm.assert_frame_equal(result, expected)

        # pass array/mask as the coordinates
        with ensure_clean_store(setup_path) as store:

            df = DataFrame(
                np.random.randn(1000, 2), index=date_range("20000101", periods=1000)
            )
            store.append("df", df)
            c = store.select_column("df", "index")
            where = c[DatetimeIndex(c).month == 5].index
            expected = df.iloc[where]

            # locations
            result = store.select("df", where=where)
            tm.assert_frame_equal(result, expected)

            # boolean
            result = store.select("df", where=where)
            tm.assert_frame_equal(result, expected)

            # invalid
            with pytest.raises(ValueError):
                store.select("df", where=np.arange(len(df), dtype="float64"))

            with pytest.raises(ValueError):
                store.select("df", where=np.arange(len(df) + 1))

            with pytest.raises(ValueError):
                store.select("df", where=np.arange(len(df)), start=5)

            with pytest.raises(ValueError):
                store.select("df", where=np.arange(len(df)), start=5, stop=10)

            # selection with filter
            selection = date_range("20000101", periods=500)
            result = store.select("df", where="index in selection")
            expected = df[df.index.isin(selection)]
            tm.assert_frame_equal(result, expected)

            # list
            df = DataFrame(np.random.randn(10, 2))
            store.append("df2", df)
            result = store.select("df2", where=[0, 3, 5])
            expected = df.iloc[[0, 3, 5]]
            tm.assert_frame_equal(result, expected)

            # boolean
            where = [True] * 10
            where[-2] = False
            result = store.select("df2", where=where)
            expected = df.loc[where]
            tm.assert_frame_equal(result, expected)

            # start/stop
            result = store.select("df2", start=5, stop=10)
            expected = df[5:10]
            tm.assert_frame_equal(result, expected)

    def test_append_to_multiple(self, setup_path):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns="{}_2".format)
        df2["foo"] = "bar"
        df = concat([df1, df2], axis=1)

        with ensure_clean_store(setup_path) as store:

            # exceptions
            with pytest.raises(ValueError):
                store.append_to_multiple(
                    {"df1": ["A", "B"], "df2": None}, df, selector="df3"
                )

            with pytest.raises(ValueError):
                store.append_to_multiple({"df1": None, "df2": None}, df, selector="df3")

            with pytest.raises(ValueError):
                store.append_to_multiple("df1", df, "df1")

            # regular operation
            store.append_to_multiple(
                {"df1": ["A", "B"], "df2": None}, df, selector="df1"
            )
            result = store.select_as_multiple(
                ["df1", "df2"], where=["A>0", "B>0"], selector="df1"
            )
            expected = df[(df.A > 0) & (df.B > 0)]
            tm.assert_frame_equal(result, expected)

    def test_append_to_multiple_dropna(self, setup_path):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns="{}_2".format)
        df1.iloc[1, df1.columns.get_indexer(["A", "B"])] = np.nan
        df = concat([df1, df2], axis=1)

        with ensure_clean_store(setup_path) as store:

            # dropna=True should guarantee rows are synchronized
            store.append_to_multiple(
                {"df1": ["A", "B"], "df2": None}, df, selector="df1", dropna=True
            )
            result = store.select_as_multiple(["df1", "df2"])
            expected = df.dropna()
            tm.assert_frame_equal(result, expected)
            tm.assert_index_equal(store.select("df1").index, store.select("df2").index)

    @pytest.mark.xfail(
        run=False, reason="append_to_multiple_dropna_false is not raising as failed"
    )
    def test_append_to_multiple_dropna_false(self, setup_path):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns="{}_2".format)
        df1.iloc[1, df1.columns.get_indexer(["A", "B"])] = np.nan
        df = concat([df1, df2], axis=1)

        with ensure_clean_store(setup_path) as store:

            # dropna=False shouldn't synchronize row indexes
            store.append_to_multiple(
                {"df1a": ["A", "B"], "df2a": None}, df, selector="df1a", dropna=False
            )

            with pytest.raises(ValueError):
                store.select_as_multiple(["df1a", "df2a"])

            assert not store.select("df1a").index.equals(store.select("df2a").index)

    def test_select_as_multiple(self, setup_path):

        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns="{}_2".format)
        df2["foo"] = "bar"

        with ensure_clean_store(setup_path) as store:

            # no tables stored
            with pytest.raises(Exception):
                store.select_as_multiple(None, where=["A>0", "B>0"], selector="df1")

            store.append("df1", df1, data_columns=["A", "B"])
            store.append("df2", df2)

            # exceptions
            with pytest.raises(Exception):
                store.select_as_multiple(None, where=["A>0", "B>0"], selector="df1")

            with pytest.raises(Exception):
                store.select_as_multiple([None], where=["A>0", "B>0"], selector="df1")

            msg = "'No object named df3 in the file'"
            with pytest.raises(KeyError, match=msg):
                store.select_as_multiple(
                    ["df1", "df3"], where=["A>0", "B>0"], selector="df1"
                )

            with pytest.raises(KeyError, match=msg):
                store.select_as_multiple(["df3"], where=["A>0", "B>0"], selector="df1")

            with pytest.raises(KeyError, match="'No object named df4 in the file'"):
                store.select_as_multiple(
                    ["df1", "df2"], where=["A>0", "B>0"], selector="df4"
                )

            # default select
            result = store.select("df1", ["A>0", "B>0"])
            expected = store.select_as_multiple(
                ["df1"], where=["A>0", "B>0"], selector="df1"
            )
            tm.assert_frame_equal(result, expected)
            expected = store.select_as_multiple(
                "df1", where=["A>0", "B>0"], selector="df1"
            )
            tm.assert_frame_equal(result, expected)

            # multiple
            result = store.select_as_multiple(
                ["df1", "df2"], where=["A>0", "B>0"], selector="df1"
            )
            expected = concat([df1, df2], axis=1)
            expected = expected[(expected.A > 0) & (expected.B > 0)]
            tm.assert_frame_equal(result, expected)

            # multiple (diff selector)
            result = store.select_as_multiple(
                ["df1", "df2"], where="index>df2.index[4]", selector="df2"
            )
            expected = concat([df1, df2], axis=1)
            expected = expected[5:]
            tm.assert_frame_equal(result, expected)

            # test exception for diff rows
            store.append("df3", tm.makeTimeDataFrame(nper=50))
            with pytest.raises(ValueError):
                store.select_as_multiple(
                    ["df1", "df3"], where=["A>0", "B>0"], selector="df1"
                )

    @pytest.mark.skipif(
        LooseVersion(tables.__version__) < LooseVersion("3.1.0"),
        reason=("tables version does not support fix for nan selection bug: GH 4858"),
    )
    def test_nan_selection_bug_4858(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            df = DataFrame(dict(cols=range(6), values=range(6)), dtype="float64")
            df["cols"] = (df["cols"] + 10).apply(str)
            df.iloc[0] = np.nan

            expected = DataFrame(
                dict(cols=["13.0", "14.0", "15.0"], values=[3.0, 4.0, 5.0]),
                index=[3, 4, 5],
            )

            # write w/o the index on that particular column
            store.append("df", df, data_columns=True, index=["cols"])
            result = store.select("df", where="values>2.0")
            tm.assert_frame_equal(result, expected)

    def test_start_stop_table(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # table
            df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
            store.append("df", df)

            result = store.select("df", "columns=['A']", start=0, stop=5)
            expected = df.loc[0:4, ["A"]]
            tm.assert_frame_equal(result, expected)

            # out of range
            result = store.select("df", "columns=['A']", start=30, stop=40)
            assert len(result) == 0
            expected = df.loc[30:40, ["A"]]
            tm.assert_frame_equal(result, expected)

    def test_start_stop_multiple(self, setup_path):

        # GH 16209
        with ensure_clean_store(setup_path) as store:

            df = DataFrame({"foo": [1, 2], "bar": [1, 2]})

            store.append_to_multiple(
                {"selector": ["foo"], "data": None}, df, selector="selector"
            )
            result = store.select_as_multiple(
                ["selector", "data"], selector="selector", start=0, stop=1
            )
            expected = df.loc[[0], ["foo", "bar"]]
            tm.assert_frame_equal(result, expected)

    def test_start_stop_fixed(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # fixed, GH 8287
            df = DataFrame(
                dict(A=np.random.rand(20), B=np.random.rand(20)),
                index=pd.date_range("20130101", periods=20),
            )
            store.put("df", df)

            result = store.select("df", start=0, stop=5)
            expected = df.iloc[0:5, :]
            tm.assert_frame_equal(result, expected)

            result = store.select("df", start=5, stop=10)
            expected = df.iloc[5:10, :]
            tm.assert_frame_equal(result, expected)

            # out of range
            result = store.select("df", start=30, stop=40)
            expected = df.iloc[30:40, :]
            tm.assert_frame_equal(result, expected)

            # series
            s = df.A
            store.put("s", s)
            result = store.select("s", start=0, stop=5)
            expected = s.iloc[0:5]
            tm.assert_series_equal(result, expected)

            result = store.select("s", start=5, stop=10)
            expected = s.iloc[5:10]
            tm.assert_series_equal(result, expected)

            # sparse; not implemented
            df = tm.makeDataFrame()
            df.iloc[3:5, 1:3] = np.nan
            df.iloc[8:10, -2] = np.nan

    def test_select_filter_corner(self, setup_path):

        df = DataFrame(np.random.randn(50, 100))
        df.index = [f"{c:3d}" for c in df.index]
        df.columns = [f"{c:3d}" for c in df.columns]

        with ensure_clean_store(setup_path) as store:
            store.put("frame", df, format="table")

            crit = "columns=df.columns[:75]"
            result = store.select("frame", [crit])
            tm.assert_frame_equal(result, df.loc[:, df.columns[:75]])

            crit = "columns=df.columns[:75:2]"
            result = store.select("frame", [crit])
            tm.assert_frame_equal(result, df.loc[:, df.columns[:75:2]])

    def test_path_pathlib(self, setup_path):
        df = tm.makeDataFrame()

        result = tm.round_trip_pathlib(
            lambda p: df.to_hdf(p, "df"), lambda p: pd.read_hdf(p, "df")
        )
        tm.assert_frame_equal(df, result)

    @pytest.mark.parametrize("start, stop", [(0, 2), (1, 2), (None, None)])
    def test_contiguous_mixed_data_table(self, start, stop, setup_path):
        # GH 17021
        # ValueError when reading a contiguous mixed-data table ft. VLArray
        df = DataFrame(
            {
                "a": Series([20111010, 20111011, 20111012]),
                "b": Series(["ab", "cd", "ab"]),
            }
        )

        with ensure_clean_store(setup_path) as store:
            store.append("test_dataset", df)

            result = store.select("test_dataset", start=start, stop=stop)
            tm.assert_frame_equal(df[start:stop], result)

    def test_path_pathlib_hdfstore(self, setup_path):
        df = tm.makeDataFrame()

        def writer(path):
            with pd.HDFStore(path) as store:
                df.to_hdf(store, "df")

        def reader(path):
            with pd.HDFStore(path) as store:
                return pd.read_hdf(store, "df")

        result = tm.round_trip_pathlib(writer, reader)
        tm.assert_frame_equal(df, result)

    def test_pickle_path_localpath(self, setup_path):
        df = tm.makeDataFrame()
        result = tm.round_trip_pathlib(
            lambda p: df.to_hdf(p, "df"), lambda p: pd.read_hdf(p, "df")
        )
        tm.assert_frame_equal(df, result)

    def test_path_localpath_hdfstore(self, setup_path):
        df = tm.makeDataFrame()

        def writer(path):
            with pd.HDFStore(path) as store:
                df.to_hdf(store, "df")

        def reader(path):
            with pd.HDFStore(path) as store:
                return pd.read_hdf(store, "df")

        result = tm.round_trip_localpath(writer, reader)
        tm.assert_frame_equal(df, result)

    def _check_roundtrip(self, obj, comparator, path, compression=False, **kwargs):

        options = {}
        if compression:
            options["complib"] = _default_compressor

        with ensure_clean_store(path, "w", **options) as store:
            store["obj"] = obj
            retrieved = store["obj"]
            comparator(retrieved, obj, **kwargs)

    def _check_double_roundtrip(
        self, obj, comparator, path, compression=False, **kwargs
    ):
        options = {}
        if compression:
            options["complib"] = compression or _default_compressor

        with ensure_clean_store(path, "w", **options) as store:
            store["obj"] = obj
            retrieved = store["obj"]
            comparator(retrieved, obj, **kwargs)
            store["obj"] = retrieved
            again = store["obj"]
            comparator(again, obj, **kwargs)

    def _check_roundtrip_table(self, obj, comparator, path, compression=False):
        options = {}
        if compression:
            options["complib"] = _default_compressor

        with ensure_clean_store(path, "w", **options) as store:
            store.put("obj", obj, format="table")
            retrieved = store["obj"]

            comparator(retrieved, obj)

    def test_multiple_open_close(self, setup_path):
        # gh-4409: open & close multiple times

        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()
            df.to_hdf(path, "df", mode="w", format="table")

            # single
            store = HDFStore(path)
            assert "CLOSED" not in store.info()
            assert store.is_open

            store.close()
            assert "CLOSED" in store.info()
            assert not store.is_open

        with ensure_clean_path(setup_path) as path:

            if pytables._table_file_open_policy_is_strict:

                # multiples
                store1 = HDFStore(path)

                with pytest.raises(ValueError):
                    HDFStore(path)

                store1.close()
            else:

                # multiples
                store1 = HDFStore(path)
                store2 = HDFStore(path)

                assert "CLOSED" not in store1.info()
                assert "CLOSED" not in store2.info()
                assert store1.is_open
                assert store2.is_open

                store1.close()
                assert "CLOSED" in store1.info()
                assert not store1.is_open
                assert "CLOSED" not in store2.info()
                assert store2.is_open

                store2.close()
                assert "CLOSED" in store1.info()
                assert "CLOSED" in store2.info()
                assert not store1.is_open
                assert not store2.is_open

                # nested close
                store = HDFStore(path, mode="w")
                store.append("df", df)

                store2 = HDFStore(path)
                store2.append("df2", df)
                store2.close()
                assert "CLOSED" in store2.info()
                assert not store2.is_open

                store.close()
                assert "CLOSED" in store.info()
                assert not store.is_open

                # double closing
                store = HDFStore(path, mode="w")
                store.append("df", df)

                store2 = HDFStore(path)
                store.close()
                assert "CLOSED" in store.info()
                assert not store.is_open

                store2.close()
                assert "CLOSED" in store2.info()
                assert not store2.is_open

        # ops on a closed store
        with ensure_clean_path(setup_path) as path:

            df = tm.makeDataFrame()
            df.to_hdf(path, "df", mode="w", format="table")

            store = HDFStore(path)
            store.close()

            with pytest.raises(ClosedFileError):
                store.keys()

            with pytest.raises(ClosedFileError):
                "df" in store

            with pytest.raises(ClosedFileError):
                len(store)

            with pytest.raises(ClosedFileError):
                store["df"]

            with pytest.raises(AttributeError):
                store.df

            with pytest.raises(ClosedFileError):
                store.select("df")

            with pytest.raises(ClosedFileError):
                store.get("df")

            with pytest.raises(ClosedFileError):
                store.append("df2", df)

            with pytest.raises(ClosedFileError):
                store.put("df3", df)

            with pytest.raises(ClosedFileError):
                store.get_storer("df2")

            with pytest.raises(ClosedFileError):
                store.remove("df2")

            with pytest.raises(ClosedFileError, match="file is not open"):
                store.select("df")

    def test_pytables_native_read(self, datapath, setup_path):
        with ensure_clean_store(
            datapath("io", "data", "legacy_hdf/pytables_native.h5"), mode="r"
        ) as store:
            d2 = store["detector/readout"]
            assert isinstance(d2, DataFrame)

    @pytest.mark.skipif(
        is_platform_windows(), reason="native2 read fails oddly on windows"
    )
    def test_pytables_native2_read(self, datapath, setup_path):
        with ensure_clean_store(
            datapath("io", "data", "legacy_hdf", "pytables_native2.h5"), mode="r"
        ) as store:
            str(store)
            d1 = store["detector"]
            assert isinstance(d1, DataFrame)

    @td.xfail_non_writeable
    def test_legacy_table_fixed_format_read_py2(self, datapath, setup_path):
        # GH 24510
        # legacy table with fixed format written in Python 2
        with ensure_clean_store(
            datapath("io", "data", "legacy_hdf", "legacy_table_fixed_py2.h5"), mode="r"
        ) as store:
            result = store.select("df")
            expected = pd.DataFrame(
                [[1, 2, 3, "D"]],
                columns=["A", "B", "C", "D"],
                index=pd.Index(["ABC"], name="INDEX_NAME"),
            )
            tm.assert_frame_equal(expected, result)

    def test_legacy_table_fixed_format_read_datetime_py2(self, datapath, setup_path):
        # GH 31750
        # legacy table with fixed format and datetime64 column written in Python 2
        with ensure_clean_store(
            datapath("io", "data", "legacy_hdf", "legacy_table_fixed_datetime_py2.h5"),
            mode="r",
        ) as store:
            result = store.select("df")
            expected = pd.DataFrame(
                [[pd.Timestamp("2020-02-06T18:00")]],
                columns=["A"],
                index=pd.Index(["date"]),
            )
            tm.assert_frame_equal(expected, result)

    def test_legacy_table_read_py2(self, datapath, setup_path):
        # issue: 24925
        # legacy table written in Python 2
        with ensure_clean_store(
            datapath("io", "data", "legacy_hdf", "legacy_table_py2.h5"), mode="r"
        ) as store:
            result = store.select("table")

        expected = pd.DataFrame({"a": ["a", "b"], "b": [2, 3]})
        tm.assert_frame_equal(expected, result)

    def test_copy(self, setup_path):

        with catch_warnings(record=True):

            def do_copy(f, new_f=None, keys=None, propindexes=True, **kwargs):
                try:
                    store = HDFStore(f, "r")

                    if new_f is None:
                        import tempfile

                        fd, new_f = tempfile.mkstemp()

                    tstore = store.copy(
                        new_f, keys=keys, propindexes=propindexes, **kwargs
                    )

                    # check keys
                    if keys is None:
                        keys = store.keys()
                    assert set(keys) == set(tstore.keys())

                    # check indices & nrows
                    for k in tstore.keys():
                        if tstore.get_storer(k).is_table:
                            new_t = tstore.get_storer(k)
                            orig_t = store.get_storer(k)

                            assert orig_t.nrows == new_t.nrows

                            # check propindixes
                            if propindexes:
                                for a in orig_t.axes:
                                    if a.is_indexed:
                                        assert new_t[a.name].is_indexed

                finally:
                    safe_close(store)
                    safe_close(tstore)
                    try:
                        os.close(fd)
                    except (OSError, ValueError):
                        pass
                    safe_remove(new_f)

            # new table
            df = tm.makeDataFrame()

            try:
                path = create_tempfile(setup_path)
                st = HDFStore(path)
                st.append("df", df, data_columns=["A"])
                st.close()
                do_copy(f=path)
                do_copy(f=path, propindexes=False)
            finally:
                safe_remove(path)

    def test_store_datetime_fractional_secs(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            dt = datetime.datetime(2012, 1, 2, 3, 4, 5, 123456)
            series = Series([0], [dt])
            store["a"] = series
            assert store["a"].index[0] == dt

    def test_tseries_indices_series(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            idx = tm.makeDateIndex(10)
            ser = Series(np.random.randn(len(idx)), idx)
            store["a"] = ser
            result = store["a"]

            tm.assert_series_equal(result, ser)
            assert result.index.freq == ser.index.freq
            tm.assert_class_equal(result.index, ser.index, obj="series index")

            idx = tm.makePeriodIndex(10)
            ser = Series(np.random.randn(len(idx)), idx)
            store["a"] = ser
            result = store["a"]

            tm.assert_series_equal(result, ser)
            assert result.index.freq == ser.index.freq
            tm.assert_class_equal(result.index, ser.index, obj="series index")

    def test_tseries_indices_frame(self, setup_path):

        with ensure_clean_store(setup_path) as store:
            idx = tm.makeDateIndex(10)
            df = DataFrame(np.random.randn(len(idx), 3), index=idx)
            store["a"] = df
            result = store["a"]

            tm.assert_frame_equal(result, df)
            assert result.index.freq == df.index.freq
            tm.assert_class_equal(result.index, df.index, obj="dataframe index")

            idx = tm.makePeriodIndex(10)
            df = DataFrame(np.random.randn(len(idx), 3), idx)
            store["a"] = df
            result = store["a"]

            tm.assert_frame_equal(result, df)
            assert result.index.freq == df.index.freq
            tm.assert_class_equal(result.index, df.index, obj="dataframe index")

    def test_unicode_index(self, setup_path):

        unicode_values = ["\u03c3", "\u03c3\u03c3"]

        # PerformanceWarning
        with catch_warnings(record=True):
            simplefilter("ignore", pd.errors.PerformanceWarning)
            s = Series(np.random.randn(len(unicode_values)), unicode_values)
            self._check_roundtrip(s, tm.assert_series_equal, path=setup_path)

    def test_unicode_longer_encoded(self, setup_path):
        # GH 11234
        char = "\u0394"
        df = pd.DataFrame({"A": [char]})
        with ensure_clean_store(setup_path) as store:
            store.put("df", df, format="table", encoding="utf-8")
            result = store.get("df")
            tm.assert_frame_equal(result, df)

        df = pd.DataFrame({"A": ["a", char], "B": ["b", "b"]})
        with ensure_clean_store(setup_path) as store:
            store.put("df", df, format="table", encoding="utf-8")
            result = store.get("df")
            tm.assert_frame_equal(result, df)

    @td.xfail_non_writeable
    def test_store_datetime_mixed(self, setup_path):

        df = DataFrame({"a": [1, 2, 3], "b": [1.0, 2.0, 3.0], "c": ["a", "b", "c"]})
        ts = tm.makeTimeSeries()
        df["d"] = ts.index[:3]
        self._check_roundtrip(df, tm.assert_frame_equal, path=setup_path)

    # FIXME: don't leave commented-out code
    # def test_cant_write_multiindex_table(self):
    #     # for now, #1848
    #     df = DataFrame(np.random.randn(10, 4),
    #                    index=[np.arange(5).repeat(2),
    #                           np.tile(np.arange(2), 5)])
    #
    #     with pytest.raises(Exception):
    #         store.put('foo', df, format='table')

    def test_append_with_diff_col_name_types_raises_value_error(self, setup_path):
        df = DataFrame(np.random.randn(10, 1))
        df2 = DataFrame({"a": np.random.randn(10)})
        df3 = DataFrame({(1, 2): np.random.randn(10)})
        df4 = DataFrame({("1", 2): np.random.randn(10)})
        df5 = DataFrame({("1", 2, object): np.random.randn(10)})

        with ensure_clean_store(setup_path) as store:
            name = f"df_{tm.rands(10)}"
            store.append(name, df)

            for d in (df2, df3, df4, df5):
                with pytest.raises(ValueError):
                    store.append(name, d)

    def test_query_with_nested_special_character(self, setup_path):
        df = DataFrame(
            {
                "a": ["a", "a", "c", "b", "test & test", "c", "b", "e"],
                "b": [1, 2, 3, 4, 5, 6, 7, 8],
            }
        )
        expected = df[df.a == "test & test"]
        with ensure_clean_store(setup_path) as store:
            store.append("test", df, format="table", data_columns=True)
            result = store.select("test", 'a = "test & test"')
        tm.assert_frame_equal(expected, result)

    def test_categorical(self, setup_path):

        with ensure_clean_store(setup_path) as store:

            # Basic
            _maybe_remove(store, "s")
            s = Series(
                Categorical(
                    ["a", "b", "b", "a", "a", "c"],
                    categories=["a", "b", "c", "d"],
                    ordered=False,
                )
            )
            store.append("s", s, format="table")
            result = store.select("s")
            tm.assert_series_equal(s, result)

            _maybe_remove(store, "s_ordered")
            s = Series(
                Categorical(
                    ["a", "b", "b", "a", "a", "c"],
                    categories=["a", "b", "c", "d"],
                    ordered=True,
                )
            )
            store.append("s_ordered", s, format="table")
            result = store.select("s_ordered")
            tm.assert_series_equal(s, result)

            _maybe_remove(store, "df")
            df = DataFrame({"s": s, "vals": [1, 2, 3, 4, 5, 6]})
            store.append("df", df, format="table")
            result = store.select("df")
            tm.assert_frame_equal(result, df)

            # Dtypes
            _maybe_remove(store, "si")
            s = Series([1, 1, 2, 2, 3, 4, 5]).astype("category")
            store.append("si", s)
            result = store.select("si")
            tm.assert_series_equal(result, s)

            _maybe_remove(store, "si2")
            s = Series([1, 1, np.nan, 2, 3, 4, 5]).astype("category")
            store.append("si2", s)
            result = store.select("si2")
            tm.assert_series_equal(result, s)

            # Multiple
            _maybe_remove(store, "df2")
            df2 = df.copy()
            df2["s2"] = Series(list("abcdefg")).astype("category")
            store.append("df2", df2)
            result = store.select("df2")
            tm.assert_frame_equal(result, df2)

            # Make sure the metadata is OK
            info = store.info()
            assert "/df2   " in info
            # assert '/df2/meta/values_block_0/meta' in info
            assert "/df2/meta/values_block_1/meta" in info

            # unordered
            _maybe_remove(store, "s2")
            s = Series(
                Categorical(
                    ["a", "b", "b", "a", "a", "c"],
                    categories=["a", "b", "c", "d"],
                    ordered=False,
                )
            )
            store.append("s2", s, format="table")
            result = store.select("s2")
            tm.assert_series_equal(result, s)

            # Query
            _maybe_remove(store, "df3")
            store.append("df3", df, data_columns=["s"])
            expected = df[df.s.isin(["b", "c"])]
            result = store.select("df3", where=['s in ["b","c"]'])
            tm.assert_frame_equal(result, expected)

            expected = df[df.s.isin(["b", "c"])]
            result = store.select("df3", where=['s = ["b","c"]'])
            tm.assert_frame_equal(result, expected)

            expected = df[df.s.isin(["d"])]
            result = store.select("df3", where=['s in ["d"]'])
            tm.assert_frame_equal(result, expected)

            expected = df[df.s.isin(["f"])]
            result = store.select("df3", where=['s in ["f"]'])
            tm.assert_frame_equal(result, expected)

            # Appending with same categories is ok
            store.append("df3", df)

            df = concat([df, df])
            expected = df[df.s.isin(["b", "c"])]
            result = store.select("df3", where=['s in ["b","c"]'])
            tm.assert_frame_equal(result, expected)

            # Appending must have the same categories
            df3 = df.copy()
            df3["s"].cat.remove_unused_categories(inplace=True)

            with pytest.raises(ValueError):
                store.append("df3", df3)

            # Remove, and make sure meta data is removed (its a recursive
            # removal so should be).
            result = store.select("df3/meta/s/meta")
            assert result is not None
            store.remove("df3")

            with pytest.raises(
                KeyError, match="'No object named df3/meta/s/meta in the file'"
            ):
                store.select("df3/meta/s/meta")

    def test_categorical_conversion(self, setup_path):

        # GH13322
        # Check that read_hdf with categorical columns doesn't return rows if
        # where criteria isn't met.
        obsids = ["ESP_012345_6789", "ESP_987654_3210"]
        imgids = ["APF00006np", "APF0001imm"]
        data = [4.3, 9.8]

        # Test without categories
        df = DataFrame(dict(obsids=obsids, imgids=imgids, data=data))

        # We are expecting an empty DataFrame matching types of df
        expected = df.iloc[[], :]
        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", format="table", data_columns=True)
            result = read_hdf(path, "df", where="obsids=B")
            tm.assert_frame_equal(result, expected)

        # Test with categories
        df.obsids = df.obsids.astype("category")
        df.imgids = df.imgids.astype("category")

        # We are expecting an empty DataFrame matching types of df
        expected = df.iloc[[], :]
        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", format="table", data_columns=True)
            result = read_hdf(path, "df", where="obsids=B")
            tm.assert_frame_equal(result, expected)

    def test_categorical_nan_only_columns(self, setup_path):
        # GH18413
        # Check that read_hdf with categorical columns with NaN-only values can
        # be read back.
        df = pd.DataFrame(
            {
                "a": ["a", "b", "c", np.nan],
                "b": [np.nan, np.nan, np.nan, np.nan],
                "c": [1, 2, 3, 4],
                "d": pd.Series([None] * 4, dtype=object),
            }
        )
        df["a"] = df.a.astype("category")
        df["b"] = df.b.astype("category")
        df["d"] = df.b.astype("category")
        expected = df
        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", format="table", data_columns=True)
            result = read_hdf(path, "df")
            tm.assert_frame_equal(result, expected)

    def test_duplicate_column_name(self, setup_path):
        df = DataFrame(columns=["a", "a"], data=[[0, 0]])

        with ensure_clean_path(setup_path) as path:
            with pytest.raises(ValueError):
                df.to_hdf(path, "df", format="fixed")

            df.to_hdf(path, "df", format="table")
            other = read_hdf(path, "df")

            tm.assert_frame_equal(df, other)
            assert df.equals(other)
            assert other.equals(df)

    def test_round_trip_equals(self, setup_path):
        # GH 9330
        df = DataFrame({"B": [1, 2], "A": ["x", "y"]})

        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", format="table")
            other = read_hdf(path, "df")
            tm.assert_frame_equal(df, other)
            assert df.equals(other)
            assert other.equals(df)

    def test_preserve_timedeltaindex_type(self, setup_path):
        # GH9635
        # Storing TimedeltaIndexed DataFrames in fixed stores did not preserve
        # the type of the index.
        df = DataFrame(np.random.normal(size=(10, 5)))
        df.index = timedelta_range(start="0s", periods=10, freq="1s", name="example")

        with ensure_clean_store(setup_path) as store:

            store["df"] = df
            tm.assert_frame_equal(store["df"], df)

    def test_columns_multiindex_modified(self, setup_path):
        # BUG: 7212
        # read_hdf store.select modified the passed columns parameters
        # when multi-indexed.

        df = DataFrame(np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE"))
        df.index.name = "letters"
        df = df.set_index(keys="E", append=True)

        data_columns = df.index.names + df.columns.tolist()
        with ensure_clean_path(setup_path) as path:
            df.to_hdf(
                path,
                "df",
                mode="a",
                append=True,
                data_columns=data_columns,
                index=False,
            )
            cols2load = list("BCD")
            cols2load_original = list(cols2load)
            df_loaded = read_hdf(path, "df", columns=cols2load)  # noqa
            assert cols2load_original == cols2load

    @ignore_natural_naming_warning
    def test_to_hdf_with_object_column_names(self, setup_path):
        # GH9057
        # Writing HDF5 table format should only work for string-like
        # column types

        types_should_fail = [
            tm.makeIntIndex,
            tm.makeFloatIndex,
            tm.makeDateIndex,
            tm.makeTimedeltaIndex,
            tm.makePeriodIndex,
        ]
        types_should_run = [
            tm.makeStringIndex,
            tm.makeCategoricalIndex,
            tm.makeUnicodeIndex,
        ]

        for index in types_should_fail:
            df = DataFrame(np.random.randn(10, 2), columns=index(2))
            with ensure_clean_path(setup_path) as path:
                with catch_warnings(record=True):
                    msg = "cannot have non-object label DataIndexableCol"
                    with pytest.raises(ValueError, match=msg):
                        df.to_hdf(path, "df", format="table", data_columns=True)

        for index in types_should_run:
            df = DataFrame(np.random.randn(10, 2), columns=index(2))
            with ensure_clean_path(setup_path) as path:
                with catch_warnings(record=True):
                    df.to_hdf(path, "df", format="table", data_columns=True)
                    result = pd.read_hdf(path, "df", where=f"index = [{df.index[0]}]")
                    assert len(result)

    def test_read_hdf_open_store(self, setup_path):
        # GH10330
        # No check for non-string path_or-buf, and no test of open store
        df = DataFrame(np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE"))
        df.index.name = "letters"
        df = df.set_index(keys="E", append=True)

        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", mode="w")
            direct = read_hdf(path, "df")
            store = HDFStore(path, mode="r")
            indirect = read_hdf(store, "df")
            tm.assert_frame_equal(direct, indirect)
            assert store.is_open
            store.close()

    def test_read_hdf_iterator(self, setup_path):
        df = DataFrame(np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE"))
        df.index.name = "letters"
        df = df.set_index(keys="E", append=True)

        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", mode="w", format="t")
            direct = read_hdf(path, "df")
            iterator = read_hdf(path, "df", iterator=True)
            assert isinstance(iterator, TableIterator)
            indirect = next(iterator.__iter__())
            tm.assert_frame_equal(direct, indirect)
            iterator.store.close()

    def test_read_hdf_errors(self, setup_path):
        df = DataFrame(np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE"))

        with ensure_clean_path(setup_path) as path:
            with pytest.raises(IOError):
                read_hdf(path, "key")

            df.to_hdf(path, "df")
            store = HDFStore(path, mode="r")
            store.close()

            with pytest.raises(IOError):
                read_hdf(store, "df")

    def test_read_hdf_generic_buffer_errors(self):
        with pytest.raises(NotImplementedError):
            read_hdf(BytesIO(b""), "df")

    def test_invalid_complib(self, setup_path):
        df = DataFrame(np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE"))
        with ensure_clean_path(setup_path) as path:
            with pytest.raises(ValueError):
                df.to_hdf(path, "df", complib="foolib")

    # GH10443

    def test_read_nokey(self, setup_path):
        df = DataFrame(np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE"))

        # Categorical dtype not supported for "fixed" format. So no need
        # to test with that dtype in the dataframe here.
        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", mode="a")
            reread = read_hdf(path)
            tm.assert_frame_equal(df, reread)
            df.to_hdf(path, "df2", mode="a")

            with pytest.raises(ValueError):
                read_hdf(path)

    def test_read_nokey_table(self, setup_path):
        # GH13231
        df = DataFrame({"i": range(5), "c": Series(list("abacd"), dtype="category")})

        with ensure_clean_path(setup_path) as path:
            df.to_hdf(path, "df", mode="a", format="table")
            reread = read_hdf(path)
            tm.assert_frame_equal(df, reread)
            df.to_hdf(path, "df2", mode="a", format="table")

            with pytest.raises(ValueError):
                read_hdf(path)

    def test_read_nokey_empty(self, setup_path):
        with ensure_clean_path(setup_path) as path:
            store = HDFStore(path)
            store.close()

            with pytest.raises(ValueError):
                read_hdf(path)

    def test_read_from_pathlib_path(self, setup_path):

        # GH11773
        expected = DataFrame(
            np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE")
        )
        with ensure_clean_path(setup_path) as filename:
            path_obj = Path(filename)

            expected.to_hdf(path_obj, "df", mode="a")
            actual = read_hdf(path_obj, "df")

        tm.assert_frame_equal(expected, actual)

    @td.skip_if_no("py.path")
    def test_read_from_py_localpath(self, setup_path):

        # GH11773
        from py.path import local as LocalPath

        expected = DataFrame(
            np.random.rand(4, 5), index=list("abcd"), columns=list("ABCDE")
        )
        with ensure_clean_path(setup_path) as filename:
            path_obj = LocalPath(filename)

            expected.to_hdf(path_obj, "df", mode="a")
            actual = read_hdf(path_obj, "df")

        tm.assert_frame_equal(expected, actual)

    def test_query_long_float_literal(self, setup_path):
        # GH 14241
        df = pd.DataFrame({"A": [1000000000.0009, 1000000000.0011, 1000000000.0015]})

        with ensure_clean_store(setup_path) as store:
            store.append("test", df, format="table", data_columns=True)

            cutoff = 1000000000.0006
            result = store.select("test", f"A < {cutoff:.4f}")
            assert result.empty

            cutoff = 1000000000.0010
            result = store.select("test", f"A > {cutoff:.4f}")
            expected = df.loc[[1, 2], :]
            tm.assert_frame_equal(expected, result)

            exact = 1000000000.0011
            result = store.select("test", f"A == {exact:.4f}")
            expected = df.loc[[1], :]
            tm.assert_frame_equal(expected, result)

    def test_query_compare_column_type(self, setup_path):
        # GH 15492
        df = pd.DataFrame(
            {
                "date": ["2014-01-01", "2014-01-02"],
                "real_date": date_range("2014-01-01", periods=2),
                "float": [1.1, 1.2],
                "int": [1, 2],
            },
            columns=["date", "real_date", "float", "int"],
        )

        with ensure_clean_store(setup_path) as store:
            store.append("test", df, format="table", data_columns=True)

            ts = pd.Timestamp("2014-01-01")  # noqa
            result = store.select("test", where="real_date > ts")
            expected = df.loc[[1], :]
            tm.assert_frame_equal(expected, result)

            for op in ["<", ">", "=="]:
                # non strings to string column always fail
                for v in [2.1, True, pd.Timestamp("2014-01-01"), pd.Timedelta(1, "s")]:
                    query = f"date {op} v"
                    with pytest.raises(TypeError):
                        store.select("test", where=query)

                # strings to other columns must be convertible to type
                v = "a"
                for col in ["int", "float", "real_date"]:
                    query = f"{col} {op} v"
                    with pytest.raises(ValueError):
                        store.select("test", where=query)

                for v, col in zip(
                    ["1", "1.1", "2014-01-01"], ["int", "float", "real_date"]
                ):
                    query = f"{col} {op} v"
                    result = store.select("test", where=query)

                    if op == "==":
                        expected = df.loc[[0], :]
                    elif op == ">":
                        expected = df.loc[[1], :]
                    else:
                        expected = df.loc[[], :]
                    tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize("format", ["fixed", "table"])
    def test_read_hdf_series_mode_r(self, format, setup_path):
        # GH 16583
        # Tests that reading a Series saved to an HDF file
        # still works if a mode='r' argument is supplied
        series = tm.makeFloatSeries()
        with ensure_clean_path(setup_path) as path:
            series.to_hdf(path, key="data", format=format)
            result = pd.read_hdf(path, key="data", mode="r")
        tm.assert_series_equal(result, series)

    def test_fspath(self):
        with tm.ensure_clean("foo.h5") as path:
            with pd.HDFStore(path) as store:
                assert os.fspath(store) == str(path)

    def test_read_py2_hdf_file_in_py3(self, datapath):
        # GH 16781

        # tests reading a PeriodIndex DataFrame written in Python2 in Python3

        # the file was generated in Python 2.7 like so:
        #
        # df = pd.DataFrame([1.,2,3], index=pd.PeriodIndex(
        #              ['2015-01-01', '2015-01-02', '2015-01-05'], freq='B'))
        # df.to_hdf('periodindex_0.20.1_x86_64_darwin_2.7.13.h5', 'p')

        expected = pd.DataFrame(
            [1.0, 2, 3],
            index=pd.PeriodIndex(["2015-01-01", "2015-01-02", "2015-01-05"], freq="B"),
        )

        with ensure_clean_store(
            datapath(
                "io", "data", "legacy_hdf", "periodindex_0.20.1_x86_64_darwin_2.7.13.h5"
            ),
            mode="r",
        ) as store:
            result = store["p"]
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("where", ["", (), (None,), [], [None]])
    def test_select_empty_where(self, where):
        # GH26610

        # Using keyword `where` as '' or (), or [None], etc
        # while reading from HDF store raises
        # "SyntaxError: only a single expression is allowed"

        df = pd.DataFrame([1, 2, 3])
        with ensure_clean_path("empty_where.h5") as path:
            with pd.HDFStore(path) as store:
                store.put("df", df, "t")
                result = pd.read_hdf(store, "df", where=where)
                tm.assert_frame_equal(result, df)

    @pytest.mark.parametrize(
        "idx",
        [
            date_range("2019", freq="D", periods=3, tz="UTC"),
            CategoricalIndex(list("abc")),
        ],
    )
    def test_to_hdf_multiindex_extension_dtype(self, idx, setup_path):
        # GH 7775
        mi = MultiIndex.from_arrays([idx, idx])
        df = pd.DataFrame(0, index=mi, columns=["a"])
        with ensure_clean_path(setup_path) as path:
            with pytest.raises(NotImplementedError, match="Saving a MultiIndex"):
                df.to_hdf(path, "df")
