import re

import numpy as np
import pytest

from pandas._libs.tslibs import Timestamp

import pandas as pd
from pandas import (
    DataFrame,
    HDFStore,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
    concat,
    date_range,
)
from pandas.util import _test_decorators as td

pytestmark = [pytest.mark.single_cpu]


def test_format_type(temp_hdfstore):
    df = DataFrame({"A": [1, 2]})
    temp_hdfstore.put("a", df, format="fixed")
    temp_hdfstore.put("b", df, format="table")

    assert temp_hdfstore.get_storer("a").format_type == "fixed"
    assert temp_hdfstore.get_storer("b").format_type == "table"


def test_format_kwarg_in_constructor(temp_h5_path):
    # GH 13291

    msg = "format is not a defined argument for HDFStore"

    with pytest.raises(ValueError, match=msg):
        HDFStore(temp_h5_path, format="table")


def test_api_default_format(temp_hdfstore):
    # default_format option
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD")),
        index=Index([f"i-{i}" for i in range(30)]),
    )

    with pd.option_context("io.hdf.default_format", "fixed"):
        temp_hdfstore.put("df", df)
        assert not temp_hdfstore.get_storer("df").is_table

        msg = "Can only append to Tables"
        with pytest.raises(ValueError, match=msg):
            temp_hdfstore.append("df2", df)

    with pd.option_context("io.hdf.default_format", "table"):
        temp_hdfstore.remove("df")
        temp_hdfstore.put("df", df)
        assert temp_hdfstore.get_storer("df").is_table

        temp_hdfstore.append("df2", df)
        assert temp_hdfstore.get_storer("df").is_table


def test_api_default_format_path(temp_h5_path):
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD")),
        index=Index([f"i-{i}" for i in range(30)]),
    )

    with pd.option_context("io.hdf.default_format", "fixed"):
        df.to_hdf(temp_h5_path, key="df")
        with HDFStore(temp_h5_path) as store:
            assert not store.get_storer("df").is_table
        msg = "Can only append to Tables"
        with pytest.raises(ValueError, match=msg):
            df.to_hdf(temp_h5_path, key="df2", append=True)

    with pd.option_context("io.hdf.default_format", "table"):
        df.to_hdf(temp_h5_path, key="df3")
        with HDFStore(temp_h5_path) as store:
            assert store.get_storer("df3").is_table
        df.to_hdf(temp_h5_path, key="df4", append=True)
        with HDFStore(temp_h5_path) as store:
            assert store.get_storer("df4").is_table


def test_put(temp_hdfstore):
    store = temp_hdfstore
    ts = Series(
        np.arange(10, dtype=np.float64), index=date_range("2020-01-01", periods=10)
    )
    df = DataFrame(
        np.random.default_rng(2).standard_normal((20, 4)),
        columns=Index(list("ABCD")),
        index=date_range("2000-01-01", periods=20, freq="B"),
    )
    store["a"] = ts
    store["b"] = df[:10]
    store["foo/bar/bah"] = df[:10]
    store["foo"] = df[:10]
    store["/foo"] = df[:10]
    store.put("c", df[:10], format="table")

    # not OK, not a table
    msg = "Can only append to Tables"
    with pytest.raises(ValueError, match=msg):
        store.put("b", df[10:], append=True)

    # node does not currently exist, test _is_table_type returns False
    # in this case
    with pytest.raises(ValueError, match=msg):
        store.put("f", df[10:], append=True)

    # can't put to a table (use append instead)
    with pytest.raises(ValueError, match=msg):
        store.put("c", df[10:], append=True)

    # overwrite table
    store.put("c", df[:10], format="table", append=False)
    tm.assert_frame_equal(df[:10], store["c"])


def test_put_string_index(temp_hdfstore):
    store = temp_hdfstore
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


def test_put_compression(temp_hdfstore):
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD")),
        index=date_range("2000-01-01", periods=10, freq="B"),
    )

    temp_hdfstore.put("c", df, format="table", complib="zlib")
    tm.assert_frame_equal(temp_hdfstore["c"], df)

    # can't compress if format='fixed'
    msg = "Compression not supported on Fixed format stores"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put("b", df, format="fixed", complib="zlib")


@td.skip_if_windows
def test_put_compression_blosc(temp_hdfstore):
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD")),
        index=date_range("2000-01-01", periods=10, freq="B"),
    )
    # can't compress if format='fixed'
    msg = "Compression not supported on Fixed format stores"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put("b", df, format="fixed", complib="blosc")

    temp_hdfstore.put("c", df, format="table", complib="blosc")
    tm.assert_frame_equal(temp_hdfstore["c"], df)


def test_put_datetime_ser(temp_hdfstore, performance_warning, using_infer_string):
    # https://github.com/pandas-dev/pandas/pull/60663
    ser = Series(3 * [Timestamp("20010102").as_unit("ns")])
    temp_hdfstore.put("ser", ser)
    expected = ser.copy()
    result = temp_hdfstore.get("ser")
    tm.assert_series_equal(result, expected)


def test_put_mixed_type(temp_hdfstore, performance_warning, using_infer_string):
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD")),
        index=date_range("2000-01-01", periods=10, freq="B"),
    )
    df["obj1"] = "foo"
    df["obj2"] = "bar"
    df["bool1"] = df["A"] > 0
    df["bool2"] = df["B"] > 0
    df["bool3"] = True
    df["int1"] = 1
    df["int2"] = 2
    df["timestamp1"] = Timestamp("20010102").as_unit("ns")
    df["timestamp2"] = Timestamp("20010103").as_unit("ns")
    df["datetime1"] = Timestamp("20010102").as_unit("ns")
    df["datetime2"] = Timestamp("20010103").as_unit("ns")
    df.loc[df.index[3:6], ["obj1"]] = np.nan
    df = df._consolidate()

    warning = None if using_infer_string else performance_warning
    with tm.assert_produces_warning(warning):
        temp_hdfstore.put("df", df)

    expected = temp_hdfstore.get("df")
    tm.assert_frame_equal(expected, df)


def test_put_str_frame(temp_hdfstore, performance_warning, string_dtype_arguments):
    # https://github.com/pandas-dev/pandas/pull/60663
    dtype = pd.StringDtype(*string_dtype_arguments)
    df = DataFrame({"a": pd.array(["x", pd.NA, "y"], dtype=dtype)})

    temp_hdfstore.put("df", df)
    expected_dtype = "str" if dtype.na_value is np.nan else "string"
    expected = df.astype(expected_dtype)
    result = temp_hdfstore.get("df")
    tm.assert_frame_equal(result, expected)


def test_put_str_series(temp_hdfstore, performance_warning, string_dtype_arguments):
    # https://github.com/pandas-dev/pandas/pull/60663
    dtype = pd.StringDtype(*string_dtype_arguments)
    ser = Series(["x", pd.NA, "y"], dtype=dtype)

    temp_hdfstore.put("ser", ser)
    expected_dtype = "str" if dtype.na_value is np.nan else "string"
    expected = ser.astype(expected_dtype)
    result = temp_hdfstore.get("ser")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("format", ["table", "fixed"])
@pytest.mark.parametrize(
    "index",
    [
        Index([str(i) for i in range(10)]),
        Index(np.arange(10, dtype=float)),
        Index(np.arange(10)),
        date_range("2020-01-01", periods=10),
        pd.period_range("2020-01-01", periods=10),
    ],
)
def test_store_index_types(temp_hdfstore, format, index):
    # GH5386
    # test storing various index types
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 2)),
        columns=list("AB"),
        index=index,
    )
    temp_hdfstore.put("df", df, format=format)
    tm.assert_frame_equal(df, temp_hdfstore["df"])


def test_column_multiindex(temp_hdfstore, using_infer_string):
    # GH 4710
    # recreate multi-indexes properly

    index = MultiIndex.from_tuples(
        [("A", "a"), ("A", "b"), ("B", "a"), ("B", "b")], names=["first", "second"]
    )
    df = DataFrame(np.arange(12).reshape(3, 4), columns=index)
    expected = df.set_axis(df.index.to_numpy())

    temp_hdfstore.put("df", df)
    tm.assert_frame_equal(
        temp_hdfstore["df"], expected, check_index_type=True, check_column_type=True
    )

    temp_hdfstore.put("df1", df, format="table")
    tm.assert_frame_equal(
        temp_hdfstore["df1"], expected, check_index_type=True, check_column_type=True
    )

    msg = re.escape("cannot use a multi-index on axis [1] with data_columns ['A']")
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put("df2", df, format="table", data_columns=["A"])
    msg = re.escape("cannot use a multi-index on axis [1] with data_columns True")
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put("df3", df, format="table", data_columns=True)


def test_column_multiindex_existing(temp_hdfstore, using_infer_string):
    # appending multi-column on existing table (see GH 6167)

    index = MultiIndex.from_tuples(
        [("A", "a"), ("A", "b"), ("B", "a"), ("B", "b")], names=["first", "second"]
    )
    df = DataFrame(np.arange(12).reshape(3, 4), columns=index)
    temp_hdfstore.append("df2", df)
    temp_hdfstore.append("df2", df)

    tm.assert_frame_equal(temp_hdfstore["df2"], concat((df, df)))


def test_column_multiindex_non_index_axes(temp_hdfstore, using_infer_string):
    df = DataFrame(np.arange(12).reshape(3, 4), columns=Index(list("ABCD"), name="foo"))
    expected = df.set_axis(df.index.to_numpy())

    temp_hdfstore.put("df1", df, format="table")
    tm.assert_frame_equal(
        temp_hdfstore["df1"], expected, check_index_type=True, check_column_type=True
    )


def test_store_multiindex(temp_hdfstore):
    # validate multi-index names
    # GH 5527

    def make_index(names=None):
        dti = date_range("2013-12-01", "2013-12-02")
        mi = MultiIndex.from_product([dti, range(2), range(3)], names=names)
        return mi

    # no names
    df = DataFrame(np.zeros((12, 2)), columns=["a", "b"], index=make_index())
    temp_hdfstore.append("df", df)
    tm.assert_frame_equal(temp_hdfstore.select("df"), df)

    # partial names
    temp_hdfstore.remove("df")
    df = DataFrame(
        np.zeros((12, 2)),
        columns=["a", "b"],
        index=make_index(["date", None, None]),
    )
    temp_hdfstore.append("df", df)
    tm.assert_frame_equal(temp_hdfstore.select("df"), df)

    # series
    ser = Series(np.zeros(12), index=make_index(["date", None, None]))
    temp_hdfstore.append("ser", ser)
    xp = Series(np.zeros(12), index=make_index(["date", "level_1", "level_2"]))
    tm.assert_series_equal(temp_hdfstore.select("ser"), xp)

    # dup with column
    temp_hdfstore.remove("df")
    df = DataFrame(
        np.zeros((12, 2)),
        columns=["a", "b"],
        index=make_index(["date", "a", "t"]),
    )
    msg = "duplicate names/columns in the multi-index when storing as a table"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df", df)

    # dup within level
    temp_hdfstore.remove("df")
    df = DataFrame(
        np.zeros((12, 2)),
        columns=["a", "b"],
        index=make_index(["date", "date", "date"]),
    )
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df", df)

    # fully names
    temp_hdfstore.remove("df")
    df = DataFrame(
        np.zeros((12, 2)),
        columns=["a", "b"],
        index=make_index(["date", "s", "t"]),
    )
    temp_hdfstore.append("df", df)
    tm.assert_frame_equal(temp_hdfstore.select("df"), df)


@pytest.mark.parametrize("format", ["fixed", "table"])
def test_store_periodindex(temp_h5_path, format):
    # GH 7796
    # test of PeriodIndex in HDFStore
    df = DataFrame(
        np.random.default_rng(2).standard_normal((5, 1)),
        index=pd.period_range("20220101", freq="M", periods=5),
    )

    df.to_hdf(temp_h5_path, key="df", mode="w", format=format)
    expected = pd.read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(df, expected)
