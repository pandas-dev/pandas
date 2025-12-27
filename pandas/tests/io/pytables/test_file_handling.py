import os
import uuid

import numpy as np
import pytest

from pandas.compat import (
    is_platform_linux,
    is_platform_little_endian,
    is_platform_mac,
)
from pandas.errors import (
    ClosedFileError,
    PossibleDataLossError,
)

from pandas import (
    DataFrame,
    HDFStore,
    Index,
    Series,
    _testing as tm,
    date_range,
    read_hdf,
)

from pandas.io import pytables
from pandas.io.pytables import Term

tables = pytest.importorskip("tables")
pytestmark = [pytest.mark.single_cpu]


@pytest.mark.parametrize("mode", ["r", "r+", "a", "w"])
def test_mode(temp_h5_path, mode, using_infer_string):
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=date_range("2000-01-01", periods=10, freq="B"),
    )
    msg = r"[\S]* does not exist"
    doesnt_exist = f"{uuid.uuid4()}.h5"

    # constructor
    if mode in ["r", "r+"]:
        with pytest.raises(OSError, match=msg):
            HDFStore(doesnt_exist, mode=mode)

    else:
        with HDFStore(temp_h5_path, mode=mode) as store:
            assert store._handle.mode == mode

    # context
    if mode in ["r", "r+"]:
        with pytest.raises(OSError, match=msg):
            with HDFStore(doesnt_exist, mode=mode) as store:
                pass
    else:
        with HDFStore(temp_h5_path, mode=mode) as store:
            assert store._handle.mode == mode

    # conv write
    if mode in ["r", "r+"]:
        with pytest.raises(OSError, match=msg):
            df.to_hdf(doesnt_exist, key="df", mode=mode)
        df.to_hdf(temp_h5_path, key="df", mode="w")
    else:
        df.to_hdf(temp_h5_path, key="df", mode=mode)

    # conv read
    if mode in ["w"]:
        msg = (
            "mode w is not allowed while performing a read. "
            r"Allowed modes are r, r\+ and a."
        )
        with pytest.raises(ValueError, match=msg):
            read_hdf(temp_h5_path, "df", mode=mode)
    else:
        result = read_hdf(temp_h5_path, "df", mode=mode)
        if using_infer_string:
            df.columns = df.columns.astype("str")
        tm.assert_frame_equal(result, df)


def test_default_mode(temp_h5_path, using_infer_string):
    # read_hdf uses default mode
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=date_range("2000-01-01", periods=10, freq="B"),
    )
    df.to_hdf(temp_h5_path, key="df", mode="w")
    result = read_hdf(temp_h5_path, "df")
    expected = df.copy()
    if using_infer_string:
        expected.columns = expected.columns.astype("str")
    tm.assert_frame_equal(result, expected)


def test_reopen_handle(temp_h5_path):
    store = HDFStore(temp_h5_path, mode="a")
    store["a"] = Series(
        np.arange(10, dtype=np.float64), index=date_range("2020-01-01", periods=10)
    )

    msg = (
        r"Re-opening the file \[[\S]*\] with mode \[a\] will delete the "
        "current file!"
    )
    # invalid mode change
    with pytest.raises(PossibleDataLossError, match=msg):
        store.open("w")

    store.close()
    assert not store.is_open

    # truncation ok here
    store.open("w")
    assert store.is_open
    assert len(store) == 0
    store.close()
    assert not store.is_open

    store = HDFStore(temp_h5_path, mode="a")
    store["a"] = Series(
        np.arange(10, dtype=np.float64), index=date_range("2020-01-01", periods=10)
    )

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


def test_open_args(using_infer_string):
    not_written = f"{uuid.uuid4()}.h5"
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )

    # create an in memory store
    store = HDFStore(
        not_written, mode="a", driver="H5FD_CORE", driver_core_backing_store=0
    )
    store["df"] = df
    store.append("df2", df)

    expected = df.copy()
    if using_infer_string:
        expected.index = expected.index.astype("str")
        expected.columns = expected.columns.astype("str")

    tm.assert_frame_equal(store["df"], expected)
    tm.assert_frame_equal(store["df2"], expected)

    store.close()

    # the file should not have actually been written
    assert not os.path.exists(not_written)


def test_flush(temp_h5_path):
    with HDFStore(temp_h5_path, mode="w") as store:
        store["a"] = Series(range(5))
        store.flush()
        store.flush(fsync=True)


def test_complibs_default_settings(temp_h5_path, using_infer_string):
    # GH15943
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )

    # Set complevel and check if complib is automatically set to
    # default value
    df.to_hdf(temp_h5_path, key="df", complevel=9)
    result = read_hdf(temp_h5_path, "df")
    expected = df.copy()
    if using_infer_string:
        expected.index = expected.index.astype("str")
        expected.columns = expected.columns.astype("str")
    tm.assert_frame_equal(result, expected)

    with tables.open_file(temp_h5_path, mode="r") as h5file:
        for node in h5file.walk_nodes(where="/df", classname="Leaf"):
            assert node.filters.complevel == 9
            assert node.filters.complib == "zlib"

    # Set complib and check to see if compression is disabled
    df.to_hdf(temp_h5_path, key="df", complib="zlib")
    result = read_hdf(temp_h5_path, "df")
    expected = df.copy()
    if using_infer_string:
        expected.index = expected.index.astype("str")
        expected.columns = expected.columns.astype("str")
    tm.assert_frame_equal(result, expected)

    with tables.open_file(temp_h5_path, mode="r") as h5file:
        for node in h5file.walk_nodes(where="/df", classname="Leaf"):
            assert node.filters.complevel == 0
            assert node.filters.complib is None

    # Check if not setting complib or complevel results in no compression
    df.to_hdf(temp_h5_path, key="df")
    result = read_hdf(temp_h5_path, "df")
    expected = df.copy()
    if using_infer_string:
        expected.index = expected.index.astype("str")
        expected.columns = expected.columns.astype("str")
    tm.assert_frame_equal(result, expected)

    with tables.open_file(temp_h5_path, mode="r") as h5file:
        for node in h5file.walk_nodes(where="/df", classname="Leaf"):
            assert node.filters.complevel == 0
            assert node.filters.complib is None


def test_complibs_default_settings_override(temp_h5_path):
    # Check if file-defaults can be overridden on a per table basis
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )
    store = HDFStore(temp_h5_path)
    store.append("dfc", df, complevel=9, complib="blosc")
    store.append("df", df)
    store.close()

    with tables.open_file(temp_h5_path, mode="r") as h5file:
        for node in h5file.walk_nodes(where="/df", classname="Leaf"):
            assert node.filters.complevel == 0
            assert node.filters.complib is None
        for node in h5file.walk_nodes(where="/dfc", classname="Leaf"):
            assert node.filters.complevel == 9
            assert node.filters.complib == "blosc"


@pytest.mark.parametrize("lvl", range(10))
@pytest.mark.parametrize("lib", tables.filters.all_complibs)
@pytest.mark.filterwarnings("ignore:object name is not a valid")
def test_complibs(tmp_path, lvl, lib, request):
    # GH14478
    if is_platform_linux() and lib == "blosc2" and lvl != 0:
        request.applymarker(pytest.mark.xfail(reason=f"Fails for {lib} on Linux"))
    df = DataFrame(
        np.ones((30, 4)), columns=list("ABCD"), index=np.arange(30).astype(np.str_)
    )

    # Remove lzo if its not available on this platform
    if not tables.which_lib_version("lzo"):
        pytest.skip("lzo not available")
    # Remove bzip2 if its not available on this platform
    if not tables.which_lib_version("bzip2"):
        pytest.skip("bzip2 not available")

    tmpfile = tmp_path / f"{lvl}_{lib}.h5"
    gname = f"{lvl}_{lib}"

    # Write and read file to see if data is consistent
    df.to_hdf(tmpfile, key=gname, complib=lib, complevel=lvl)
    result = read_hdf(tmpfile, gname)
    tm.assert_frame_equal(result, df)

    is_mac = is_platform_mac()

    # Open file and check metadata for correct amount of compression
    with tables.open_file(tmpfile, mode="r") as h5table:
        for node in h5table.walk_nodes(where="/" + gname, classname="Leaf"):
            assert node.filters.complevel == lvl
            if lvl == 0:
                assert node.filters.complib is None
            elif is_mac and lib == "blosc2":
                res = node.filters.complib
                assert res in [lib, "blosc2:blosclz"], res
            else:
                assert node.filters.complib == lib


@pytest.mark.skipif(
    not is_platform_little_endian(), reason="reason platform is not little endian"
)
def test_encoding(temp_hdfstore):
    df = DataFrame({"A": "foo", "B": "bar"}, index=range(5))
    df.loc[2, "A"] = np.nan
    df.loc[3, "B"] = np.nan
    temp_hdfstore.append("df", df, encoding="ascii")
    tm.assert_frame_equal(temp_hdfstore["df"], df)

    expected = df.reindex(columns=["A"])
    result = temp_hdfstore.select("df", Term("columns=A", encoding="ascii"))
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
@pytest.mark.parametrize("dtype", ["category", None])
def test_latin_encoding(temp_h5_path, dtype, val):
    enc = "latin-1"
    nan_rep = ""
    key = "data"

    val = [x.decode(enc) if isinstance(x, bytes) else x for x in val]
    ser = Series(val, dtype=dtype)

    ser.to_hdf(temp_h5_path, key=key, format="table", encoding=enc, nan_rep=nan_rep)
    retr = read_hdf(temp_h5_path, key)

    # TODO:(3.0): once Categorical replace deprecation is enforced,
    #  we may be able to re-simplify the construction of s_nan
    if dtype == "category":
        if nan_rep in ser.cat.categories:
            s_nan = ser.cat.remove_categories([nan_rep])
        else:
            s_nan = ser
    else:
        s_nan = ser.replace(nan_rep, np.nan)

    tm.assert_series_equal(s_nan, retr)


def test_multiple_open_close(temp_h5_path):
    # gh-4409: open & close multiple times

    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )
    df.to_hdf(temp_h5_path, key="df", mode="w", format="table")

    # single
    store = HDFStore(temp_h5_path)
    assert "CLOSED" not in store.info()
    assert store.is_open

    store.close()
    assert "CLOSED" in store.info()
    assert not store.is_open

    if pytables._table_file_open_policy_is_strict:
        # multiples
        store1 = HDFStore(temp_h5_path)
        msg = (
            r"The file [\S]* is already opened\.  Please close it before "
            r"reopening in write mode\."
        )
        with pytest.raises(ValueError, match=msg):
            HDFStore(temp_h5_path)

        store1.close()
    else:
        # multiples
        store1 = HDFStore(temp_h5_path)
        store2 = HDFStore(temp_h5_path)

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
        store = HDFStore(temp_h5_path, mode="w")
        store.append("df", df)

        store2 = HDFStore(temp_h5_path)
        store2.append("df2", df)
        store2.close()
        assert "CLOSED" in store2.info()
        assert not store2.is_open

        store.close()
        assert "CLOSED" in store.info()
        assert not store.is_open

        # double closing
        store = HDFStore(temp_h5_path, mode="w")
        store.append("df", df)

        store2 = HDFStore(temp_h5_path)
        store.close()
        assert "CLOSED" in store.info()
        assert not store.is_open

        store2.close()
        assert "CLOSED" in store2.info()
        assert not store2.is_open

    # ops on a closed store
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )
    df.to_hdf(temp_h5_path, key="df", mode="w", format="table")

    store = HDFStore(temp_h5_path)
    store.close()

    msg = r"[\S]* file is not open!"
    with pytest.raises(ClosedFileError, match=msg):
        store.keys()

    with pytest.raises(ClosedFileError, match=msg):
        "df" in store

    with pytest.raises(ClosedFileError, match=msg):
        len(store)

    with pytest.raises(ClosedFileError, match=msg):
        store["df"]

    with pytest.raises(ClosedFileError, match=msg):
        store.select("df")

    with pytest.raises(ClosedFileError, match=msg):
        store.get("df")

    with pytest.raises(ClosedFileError, match=msg):
        store.append("df2", df)

    with pytest.raises(ClosedFileError, match=msg):
        store.put("df3", df)

    with pytest.raises(ClosedFileError, match=msg):
        store.get_storer("df2")

    with pytest.raises(ClosedFileError, match=msg):
        store.remove("df2")

    with pytest.raises(ClosedFileError, match=msg):
        store.select("df")

    msg = "'HDFStore' object has no attribute 'df'"
    with pytest.raises(AttributeError, match=msg):
        store.df


def test_fspath(temp_h5_path):
    with HDFStore(temp_h5_path) as store:
        assert os.fspath(store) == str(temp_h5_path)
