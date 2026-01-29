from contextlib import closing
import re

import numpy as np
import pytest

from pandas.compat import is_platform_windows

import pandas as pd
from pandas import (
    DataFrame,
    HDFStore,
    Index,
    Series,
    _testing as tm,
    date_range,
    read_hdf,
)

from pandas.io.pytables import TableIterator

pytestmark = [pytest.mark.single_cpu]


def test_read_missing_key_close_store(temp_h5_path):
    # GH 25766
    df = DataFrame({"a": range(2), "b": range(2)})
    df.to_hdf(temp_h5_path, key="k1")

    with pytest.raises(KeyError, match="'No object named k2 in the file'"):
        read_hdf(temp_h5_path, "k2")

    # smoke test to test that file is properly closed after
    # read with KeyError before another write
    df.to_hdf(temp_h5_path, key="k2")


def test_read_index_error_close_store(temp_h5_path):
    # GH 25766
    df = DataFrame({"A": [], "B": []}, index=[])
    df.to_hdf(temp_h5_path, key="k1")

    with pytest.raises(IndexError, match=r"list index out of range"):
        read_hdf(temp_h5_path, "k1", stop=0)

    # smoke test to test that file is properly closed after
    # read with IndexError before another write
    df.to_hdf(temp_h5_path, key="k1")


def test_read_missing_key_opened_store(temp_h5_path):
    # GH 28699
    df = DataFrame({"a": range(2), "b": range(2)})
    df.to_hdf(temp_h5_path, key="k1")

    with HDFStore(temp_h5_path, "r") as store:
        with pytest.raises(KeyError, match="'No object named k2 in the file'"):
            read_hdf(store, "k2")

        # Test that the file is still open after a KeyError and that we can
        # still read from it.
        read_hdf(store, "k1")


def test_read_column(temp_hdfstore):
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD")),
        index=date_range("2000-01-01", periods=10, freq="B"),
    )

    # GH 17912
    # HDFStore.select_column should raise a KeyError
    # exception if the key is not a valid store
    with pytest.raises(KeyError, match="No object named df in the file"):
        temp_hdfstore.select_column("df", "index")

    temp_hdfstore.append("df", df)
    # error
    with pytest.raises(
        KeyError, match=re.escape("'column [foo] not found in the table'")
    ):
        temp_hdfstore.select_column("df", "foo")

    msg = re.escape("select_column() got an unexpected keyword argument 'where'")
    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.select_column("df", "index", where=["index>5"])

    # valid
    result = temp_hdfstore.select_column("df", "index")
    tm.assert_almost_equal(result.values, Series(df.index).values)
    assert isinstance(result, Series)

    # not a data indexable column
    msg = re.escape(
        "column [values_block_0] can not be extracted individually; "
        "it is not data indexable"
    )
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.select_column("df", "values_block_0")

    # a data column
    df2 = df.copy()
    df2["string"] = "foo"
    temp_hdfstore.append("df2", df2, data_columns=["string"])
    result = temp_hdfstore.select_column("df2", "string")
    tm.assert_almost_equal(result.values, df2["string"].values)

    # a data column with NaNs, result excludes the NaNs
    df3 = df.copy()
    df3["string"] = "foo"
    df3.loc[df3.index[4:6], "string"] = np.nan
    temp_hdfstore.append("df3", df3, data_columns=["string"])
    result = temp_hdfstore.select_column("df3", "string")
    tm.assert_almost_equal(result.values, df3["string"].values)

    # start/stop
    result = temp_hdfstore.select_column("df3", "string", start=2)
    tm.assert_almost_equal(result.values, df3["string"].values[2:])

    result = temp_hdfstore.select_column("df3", "string", start=-2)
    tm.assert_almost_equal(result.values, df3["string"].values[-2:])

    result = temp_hdfstore.select_column("df3", "string", stop=2)
    tm.assert_almost_equal(result.values, df3["string"].values[:2])

    result = temp_hdfstore.select_column("df3", "string", stop=-2)
    tm.assert_almost_equal(result.values, df3["string"].values[:-2])

    result = temp_hdfstore.select_column("df3", "string", start=2, stop=-2)
    tm.assert_almost_equal(result.values, df3["string"].values[2:-2])

    result = temp_hdfstore.select_column("df3", "string", start=-2, stop=2)
    tm.assert_almost_equal(result.values, df3["string"].values[-2:2])

    # GH 10392 - make sure column name is preserved
    df4 = DataFrame({"A": np.random.default_rng(2).standard_normal(10), "B": "foo"})
    temp_hdfstore.append("df4", df4, data_columns=True)
    expected = df4["B"]
    result = temp_hdfstore.select_column("df4", "B")
    tm.assert_series_equal(result, expected)


def test_pytables_native_read(datapath):
    with HDFStore(
        datapath("io", "data", "legacy_hdf/pytables_native.h5"), mode="r"
    ) as store:
        d2 = store["detector/readout"]
    assert isinstance(d2, DataFrame)


@pytest.mark.skipif(is_platform_windows(), reason="native2 read fails oddly on windows")
def test_pytables_native2_read(datapath):
    with HDFStore(
        datapath("io", "data", "legacy_hdf", "pytables_native2.h5"), mode="r"
    ) as store:
        str(store)
        d1 = store["detector"]
    assert isinstance(d1, DataFrame)


def test_read_hdf_open_store(temp_h5_path, using_infer_string):
    # GH10330
    # No check for non-string path_or-buf, and no test of open store
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=list("abcd"),
        columns=list("ABCDE"),
    )
    df.index.name = "letters"
    df = df.set_index(keys="E", append=True)

    df.to_hdf(temp_h5_path, key="df", mode="w")
    direct = read_hdf(temp_h5_path, "df")
    with HDFStore(temp_h5_path, mode="r") as store:
        indirect = read_hdf(store, "df")
        tm.assert_frame_equal(direct, indirect)
        assert store.is_open


def test_read_hdf_index_not_view(temp_h5_path):
    # GH 37441
    # Ensure that the index of the DataFrame is not a view
    # into the original recarray that pytables reads in
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=[0, 1, 2, 3],
        columns=list("ABCDE"),
    )

    df.to_hdf(temp_h5_path, key="df", mode="w", format="table")

    df2 = read_hdf(temp_h5_path, "df")
    assert df2.index._data.base is None
    tm.assert_frame_equal(df, df2)


def test_read_hdf_iterator(temp_h5_path):
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=list("abcd"),
        columns=list("ABCDE"),
    )
    df.index.name = "letters"
    df = df.set_index(keys="E", append=True)

    df.to_hdf(temp_h5_path, key="df", mode="w", format="t")
    direct = read_hdf(temp_h5_path, "df")
    iterator = read_hdf(temp_h5_path, "df", iterator=True)
    with closing(iterator.store):
        assert isinstance(iterator, TableIterator)
        indirect = next(iterator.__iter__())
    tm.assert_frame_equal(direct, indirect)


def test_read_nokey(temp_h5_path):
    # GH10443
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=list("abcd"),
        columns=list("ABCDE"),
    )

    # Categorical dtype not supported for "fixed" format. So no need
    # to test with that dtype in the dataframe here.
    df.to_hdf(temp_h5_path, key="df", mode="a")
    reread = read_hdf(temp_h5_path)
    tm.assert_frame_equal(df, reread)
    df.to_hdf(temp_h5_path, key="df2", mode="a")

    msg = "key must be provided when HDF5 file contains multiple datasets."
    with pytest.raises(ValueError, match=msg):
        read_hdf(temp_h5_path)


def test_read_nokey_table(temp_h5_path):
    # GH13231
    df = DataFrame({"i": range(5), "c": Series(list("abacd"), dtype="category")})

    df.to_hdf(temp_h5_path, key="df", mode="a", format="table")
    reread = read_hdf(temp_h5_path)
    tm.assert_frame_equal(df, reread)
    df.to_hdf(temp_h5_path, key="df2", mode="a", format="table")

    msg = "key must be provided when HDF5 file contains multiple datasets."
    with pytest.raises(ValueError, match=msg):
        read_hdf(temp_h5_path)


def test_read_nokey_empty(temp_h5_path):
    store = HDFStore(temp_h5_path)
    store.close()
    msg = re.escape(
        "Dataset(s) incompatible with Pandas data types, not table, or no "
        "datasets found in HDF5 file."
    )
    with pytest.raises(ValueError, match=msg):
        read_hdf(temp_h5_path)


def test_read_from_pathlib_path(temp_h5_path):
    # GH11773
    expected = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=list("abcd"),
        columns=list("ABCDE"),
    )

    expected.to_hdf(temp_h5_path, key="df", mode="a")
    actual = read_hdf(temp_h5_path, key="df")

    tm.assert_frame_equal(expected, actual)


@pytest.mark.parametrize("format", ["fixed", "table"])
def test_read_hdf_series_mode_r(temp_h5_path, format):
    # GH 16583
    # Tests that reading a Series saved to an HDF file
    # still works if a mode='r' argument is supplied
    series = Series(range(10), dtype=np.float64)
    series.to_hdf(temp_h5_path, key="data", format=format)
    result = read_hdf(temp_h5_path, key="data", mode="r")
    tm.assert_series_equal(result, series)


def test_read_infer_string(temp_h5_path):
    # GH#54431
    df = DataFrame({"a": ["a", "b", None]})
    df.to_hdf(temp_h5_path, key="data", format="table")
    with pd.option_context("future.infer_string", True):
        result = read_hdf(temp_h5_path, key="data", mode="r")
    expected = DataFrame(
        {"a": ["a", "b", None]},
        dtype=pd.StringDtype(na_value=np.nan),
        columns=Index(["a"], dtype=pd.StringDtype(na_value=np.nan)),
    )
    tm.assert_frame_equal(result, expected)


def test_hdfstore_read_datetime64_unit_s(temp_hdfstore):
    # GH 59004
    df_s = DataFrame(["2001-01-01", "2002-02-02"], dtype="datetime64[s]")
    temp_hdfstore.put("df_s", df_s)
    df_fromstore = temp_hdfstore.get("df_s")
    tm.assert_frame_equal(df_s, df_fromstore)
