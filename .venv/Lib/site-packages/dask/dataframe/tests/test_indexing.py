from __future__ import annotations

import contextlib

import numpy as np
import pandas as pd
import pytest

import dask
import dask.dataframe as dd
from dask.base import tokenize
from dask.dataframe._compat import (
    PANDAS_GE_210,
    PANDAS_GE_220,
    PANDAS_GE_300,
    IndexingError,
    tm,
)
from dask.dataframe.indexing import _coerce_loc_index
from dask.dataframe.utils import assert_eq, make_meta, pyarrow_strings_enabled

dsk = {
    ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=[0, 1, 3]),
    ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [3, 2, 1]}, index=[5, 6, 8]),
    ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [0, 0, 0]}, index=[9, 9, 9]),
}
meta = make_meta(
    {"a": "i8", "b": "i8"}, index=pd.Index([], "i8"), parent_meta=pd.DataFrame()
)
d = dd.repartition(pd.concat(dsk.values()), divisions=[0, 5, 9, 9])
full = d.compute()


def test_loc():
    assert d.loc[3:8].divisions[0] == 3
    assert d.loc[3:8].divisions[-1] == 8

    assert d.loc[5].divisions == (5, 5)

    assert_eq(d.loc[5], full.loc[5:5])
    assert_eq(d.loc[3:8], full.loc[3:8])
    assert_eq(d.loc[:8], full.loc[:8])
    assert_eq(d.loc[3:], full.loc[3:])
    assert_eq(d.loc[[5]], full.loc[[5]])

    assert_eq(d.a.loc[5], full.a.loc[5:5])
    assert_eq(d.a.loc[3:8], full.a.loc[3:8])
    assert_eq(d.a.loc[:8], full.a.loc[:8])
    assert_eq(d.a.loc[3:], full.a.loc[3:])
    assert_eq(d.a.loc[[5]], full.a.loc[[5]])
    assert_eq(d.a.loc[[]], full.a.loc[[]])
    assert_eq(d.a.loc[np.array([])], full.a.loc[np.array([])])

    pytest.raises(KeyError, lambda: d.loc[1000])
    assert_eq(d.loc[1000:], full.loc[1000:])
    assert_eq(d.loc[1000:2000], full.loc[1000:2000])
    assert_eq(d.loc[:-1000], full.loc[:-1000])
    assert_eq(d.loc[-2000:-1000], full.loc[-2000:-1000])

    assert sorted(d.loc[5].dask) == sorted(d.loc[5].dask)
    assert sorted(d.loc[5].dask) != sorted(d.loc[6].dask)


def test_loc_non_informative_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4]}, index=[10, 20, 30, 40])
    ddf = dd.from_pandas(df, npartitions=2, sort=True).clear_divisions()
    assert not ddf.known_divisions

    ddf.loc[20:30].compute(scheduler="sync")

    assert_eq(ddf.loc[20:30], df.loc[20:30])

    df = pd.DataFrame({"x": [1, 2, 3, 4]}, index=[10, 20, 20, 40])
    ddf = dd.from_pandas(df, npartitions=2, sort=True)
    assert_eq(ddf.loc[20], df.loc[20:20])


def test_loc_with_text_dates():
    A = dd._compat.makeTimeSeries().iloc[:5]
    B = dd._compat.makeTimeSeries().iloc[5:]
    s = dd.repartition(
        pd.concat([A, B]), divisions=[A.index.min(), B.index.min(), B.index.max()]
    )

    assert s.loc["2000":"2010"].divisions == s.divisions
    assert_eq(s.loc["2000":"2010"], s)
    assert len(s.loc["2000-01-03":"2000-01-05"].compute()) == 3


def test_loc_with_series():
    assert_eq(d.loc[d.a % 2 == 0], full.loc[full.a % 2 == 0])

    assert sorted(d.loc[d.a % 2 == 0].dask) == sorted(d.loc[d.a % 2 == 0].dask)
    assert sorted(d.loc[d.a % 2 == 0].dask) != sorted(d.loc[d.a % 3 == 0].dask)


def test_loc_with_array():
    assert_eq(d.loc[(d.a % 2 == 0).values], full.loc[(full.a % 2 == 0).values])


def test_loc_with_function():
    assert_eq(d.loc[lambda df: df["a"] > 3, :], full.loc[lambda df: df["a"] > 3, :])

    def _col_loc_fun(_df):
        return _df.columns.str.contains("b")

    assert_eq(d.loc[:, _col_loc_fun], full.loc[:, _col_loc_fun])


def test_loc_with_array_different_partition():
    df = pd.DataFrame(
        np.random.randn(20, 5),
        index=list("abcdefghijklmnopqrst"),
        columns=list("ABCDE"),
    )
    ddf = dd.from_pandas(df, 3)

    assert_eq(ddf.loc[(ddf.A > 0).values], df.loc[(df.A > 0).values])
    with pytest.raises(ValueError):
        ddf.loc[(ddf.A > 0).repartition(["a", "g", "k", "o", "t"]).values]


def test_loc_with_series_different_partition():
    df = pd.DataFrame(
        np.random.randn(20, 5),
        index=list("abcdefghijklmnopqrst"),
        columns=list("ABCDE"),
    )
    ddf = dd.from_pandas(df, 3)

    assert_eq(ddf.loc[ddf.A > 0], df.loc[df.A > 0])
    assert_eq(
        ddf.loc[(ddf.A > 0).repartition(["a", "g", "k", "o", "t"])], df.loc[df.A > 0]
    )


def test_loc_with_non_boolean_series():
    df = pd.Series(
        np.random.randn(20),
        index=list("abcdefghijklmnopqrst"),
    )
    ddf = dd.from_pandas(df, 3)

    s = pd.Series(list("bdmnat"))
    ds = dd.from_pandas(s, npartitions=3)

    msg = (
        "Cannot index with non-boolean dask Series. Try passing computed values instead"
    )
    with pytest.raises(KeyError, match=msg):
        ddf.loc[ds]

    assert_eq(ddf.loc[s], df.loc[s])

    ctx = contextlib.nullcontext()
    if pyarrow_strings_enabled():
        ctx = pytest.warns(
            UserWarning, match="converting pandas extension dtypes to arrays"
        )
    with ctx:
        with pytest.raises(KeyError, match=msg):
            ddf.loc[ds.values]

    assert_eq(ddf.loc[s.values], df.loc[s])

    ddf = ddf.clear_divisions()
    with pytest.raises(KeyError, match=msg):
        ddf.loc[ds]

    with pytest.raises(
        KeyError, match="Cannot index with list against unknown division"
    ):
        ddf.loc[s]


def test_loc2d():
    # index indexer is always regarded as slice for duplicated values
    assert_eq(d.loc[5, "a"], full.loc[5:5, "a"])
    # assert_eq(d.loc[[5], 'a'], full.loc[[5], 'a'])
    assert_eq(d.loc[5, ["a"]], full.loc[5:5, ["a"]])
    # assert_eq(d.loc[[5], ['a']], full.loc[[5], ['a']])

    assert_eq(d.loc[3:8, "a"], full.loc[3:8, "a"])
    assert_eq(d.loc[:8, "a"], full.loc[:8, "a"])
    assert_eq(d.loc[3:, "a"], full.loc[3:, "a"])
    assert_eq(d.loc[[8], "a"], full.loc[[8], "a"])

    assert_eq(d.loc[3:8, ["a"]], full.loc[3:8, ["a"]])
    assert_eq(d.loc[:8, ["a"]], full.loc[:8, ["a"]])
    assert_eq(d.loc[3:, ["a"]], full.loc[3:, ["a"]])

    # 3d
    with pytest.raises(IndexingError):
        d.loc[3, 3, 3]

    # Series should raise
    with pytest.raises(IndexingError):
        d.a.loc[3, 3]

    with pytest.raises(IndexingError):
        d.a.loc[3:, 3]

    with pytest.raises(IndexingError):
        d.a.loc[d.a % 2 == 0, 3]


def test_loc2d_with_known_divisions():
    df = pd.DataFrame(
        np.random.randn(20, 5),
        index=list("abcdefghijklmnopqrst"),
        columns=list("ABCDE"),
    )
    ddf = dd.from_pandas(df, 3)

    assert_eq(ddf.loc["a", "A"], df.loc[["a"], "A"])
    assert_eq(ddf.loc["a", ["A"]], df.loc[["a"], ["A"]])
    assert_eq(ddf.loc["a":"o", "A"], df.loc["a":"o", "A"])
    assert_eq(ddf.loc["a":"o", ["A"]], df.loc["a":"o", ["A"]])
    assert_eq(ddf.loc[["n"], ["A"]], df.loc[["n"], ["A"]])
    assert_eq(ddf.loc[["a", "c", "n"], ["A"]], df.loc[["a", "c", "n"], ["A"]])
    assert_eq(ddf.loc[["t", "b"], ["A"]], df.loc[["t", "b"], ["A"]])
    assert_eq(
        ddf.loc[["r", "r", "c", "g", "h"], ["A"]],
        df.loc[["r", "r", "c", "g", "h"], ["A"]],
    )


def test_loc2d_with_unknown_divisions():
    df = pd.DataFrame(
        np.random.randn(20, 5),
        index=list("abcdefghijklmnopqrst"),
        columns=list("ABCDE"),
    )
    ddf = dd.from_pandas(df, 3).clear_divisions()

    assert ddf.known_divisions is False

    assert_eq(ddf.loc["a", "A"], df.loc[["a"], "A"])
    assert_eq(ddf.loc["a", ["A"]], df.loc[["a"], ["A"]])
    assert_eq(ddf.loc["a":"o", "A"], df.loc["a":"o", "A"])
    assert_eq(ddf.loc["a":"o", ["A"]], df.loc["a":"o", ["A"]])


def test_getitem():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8, 9],
            "B": [9, 8, 7, 6, 5, 4, 3, 2, 1],
            "C": [True, False, True] * 3,
        },
        columns=list("ABC"),
    )
    ddf = dd.from_pandas(df, 2)
    assert_eq(ddf["A"], df["A"])
    # check cache consistency
    tm.assert_series_equal(ddf["A"]._meta, ddf._meta["A"])

    assert_eq(ddf[["A", "B"]], df[["A", "B"]])
    tm.assert_frame_equal(ddf[["A", "B"]]._meta, ddf._meta[["A", "B"]])

    assert_eq(ddf[ddf.C], df[df.C])
    tm.assert_series_equal(ddf.C._meta, ddf._meta.C)

    assert_eq(ddf[ddf.C.repartition([0, 2, 5, 8])], df[df.C])

    pytest.raises(KeyError, lambda: df["X"])
    pytest.raises(KeyError, lambda: df[["A", "X"]])
    pytest.raises(AttributeError, lambda: df.X)

    # not str/unicode
    df = pd.DataFrame(np.random.randn(10, 5))
    ddf = dd.from_pandas(df, 2)
    assert_eq(ddf[0], df[0])
    assert_eq(ddf[[1, 2]], df[[1, 2]])

    pytest.raises(KeyError, lambda: df[8])
    pytest.raises(KeyError, lambda: df[[1, 8]])


def test_getitem_slice():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8, 9],
            "B": [9, 8, 7, 6, 5, 4, 3, 2, 1],
            "C": [True, False, True] * 3,
        },
        index=list("abcdefghi"),
    )
    ddf = dd.from_pandas(df, 3)
    assert_eq(ddf["a":"e"], df["a":"e"])
    assert_eq(ddf["a":"b"], df["a":"b"])
    assert_eq(ddf["f":], df["f":])


def test_getitem_integer_slice():
    df = pd.DataFrame({"A": range(6)})
    ddf = dd.from_pandas(df, 2)
    # integer slicing is iloc based
    with pytest.raises(NotImplementedError):
        ddf[1:3]

    df = pd.DataFrame({"A": range(6)}, index=[1.0, 2.0, 3.0, 5.0, 10.0, 11.0])
    ddf = dd.from_pandas(df, 2)
    # except for float dtype indexes
    ctx = contextlib.nullcontext()
    if PANDAS_GE_210 and not PANDAS_GE_300:
        ctx = pytest.warns(FutureWarning, match="float-dtype index")
    elif PANDAS_GE_300:
        ctx = pytest.raises(NotImplementedError)
    with ctx:
        assert_eq(ddf[2:8], df[2:8])
    with ctx:
        assert_eq(ddf[2:], df[2:])
    with ctx:
        assert_eq(ddf[:8], df[:8])


def test_loc_on_numpy_datetimes():
    df = pd.DataFrame(
        {"x": [1, 2, 3]}, index=list(map(np.datetime64, ["2014", "2015", "2016"]))
    )
    a = dd.from_pandas(df, 2)
    a = a.repartition(divisions=tuple(map(np.datetime64, a.divisions)))

    assert_eq(a.loc["2014":"2015"], a.loc["2014":"2015"])


def test_loc_on_pandas_datetimes():
    df = pd.DataFrame(
        {"x": [1, 2, 3]}, index=list(map(pd.Timestamp, ["2014", "2015", "2016"]))
    )
    a = dd.from_pandas(df, 2)
    a = a.repartition(divisions=tuple(map(pd.Timestamp, a.divisions)))

    assert_eq(a.loc["2014":"2015"], a.loc["2014":"2015"])


def test_loc_datetime_no_freq():
    # https://github.com/dask/dask/issues/2389

    datetime_index = pd.date_range("2016-01-01", "2016-01-31", freq="12h")
    datetime_index.freq = None  # FORGET FREQUENCY
    df = pd.DataFrame({"num": range(len(datetime_index))}, index=datetime_index)

    ddf = dd.from_pandas(df, npartitions=1)
    slice_ = slice("2016-01-03", "2016-01-05")
    result = ddf.loc[slice_, :]
    expected = df.loc[slice_, :]
    assert_eq(result, expected)


def test_coerce_loc_index():
    for t in [pd.Timestamp, np.datetime64]:
        assert isinstance(_coerce_loc_index([t("2014")], "2014"), t)


ME = "ME" if PANDAS_GE_220 else "M"


def test_loc_timestamp_str():
    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.date_range("2011-01-01", freq="h", periods=100),
    )
    ddf = dd.from_pandas(df, 10)

    # partial string slice
    assert_eq(df.loc["2011-01-02"], ddf.loc["2011-01-02"])
    assert_eq(df.loc["2011-01-02":"2011-01-10"], ddf.loc["2011-01-02":"2011-01-10"])
    # same reso, dask result is always DataFrame
    assert_eq(
        df.loc["2011-01-02 10:00"].to_frame().T,
        ddf.loc["2011-01-02 10:00"],
        check_freq=False,
    )

    # series
    assert_eq(df.A.loc["2011-01-02"], ddf.A.loc["2011-01-02"], check_freq=False)
    assert_eq(
        df.A.loc["2011-01-02":"2011-01-10"],
        ddf.A.loc["2011-01-02":"2011-01-10"],
        check_freq=False,
    )

    # slice with timestamp (dask result must be DataFrame)
    assert_eq(
        df.loc[pd.Timestamp("2011-01-02")].to_frame().T,
        ddf.loc[pd.Timestamp("2011-01-02")],
        check_freq=False,
    )
    assert_eq(
        df.loc[pd.Timestamp("2011-01-02") : pd.Timestamp("2011-01-10")],
        ddf.loc[pd.Timestamp("2011-01-02") : pd.Timestamp("2011-01-10")],
        check_freq=False,
    )
    assert_eq(
        df.loc[pd.Timestamp("2011-01-02 10:00")].to_frame().T,
        ddf.loc[pd.Timestamp("2011-01-02 10:00")],
        check_freq=False,
    )

    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.date_range("2011-01-01", freq=ME, periods=100),
    )
    ddf = dd.from_pandas(df, 50)
    assert_eq(df.loc["2011-01"], ddf.loc["2011-01"])
    assert_eq(df.loc["2011"], ddf.loc["2011"])

    assert_eq(df.loc["2011-01":"2012-05"], ddf.loc["2011-01":"2012-05"])
    assert_eq(df.loc["2011":"2015"], ddf.loc["2011":"2015"])

    # series
    assert_eq(df.B.loc["2011-01"], ddf.B.loc["2011-01"])
    assert_eq(df.B.loc["2011"], ddf.B.loc["2011"])

    assert_eq(df.B.loc["2011-01":"2012-05"], ddf.B.loc["2011-01":"2012-05"])
    assert_eq(df.B.loc["2011":"2015"], ddf.B.loc["2011":"2015"])


def test_getitem_timestamp_str():
    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.date_range("2011-01-01", freq="h", periods=100),
    )
    ddf = dd.from_pandas(df, 10)
    assert_eq(df["2011-01-02":"2011-01-10"], ddf["2011-01-02":"2011-01-10"])

    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.date_range("2011-01-01", freq="D", periods=100),
    )
    ddf = dd.from_pandas(df, 50)
    assert_eq(df.loc["2011-01"], ddf.loc["2011-01"])
    assert_eq(df.loc["2011"], ddf.loc["2011"])

    assert_eq(df["2011-01":"2012-05"], ddf["2011-01":"2012-05"])
    assert_eq(df["2011":"2015"], ddf["2011":"2015"])


def test_loc_period_str():
    # .loc with PeriodIndex doesn't support partial string indexing
    # https://github.com/pydata/pandas/issues/13429
    # -> this started working in pandas 1.1
    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.period_range("2011-01-01", freq="h", periods=100),
    )
    ddf = dd.from_pandas(df, 10)

    # partial string slice
    assert_eq(df.loc["2011-01-02"], ddf.loc["2011-01-02"])
    assert_eq(df.loc["2011-01-02":"2011-01-10"], ddf.loc["2011-01-02":"2011-01-10"])
    # same reso, dask result is always DataFrame

    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.period_range("2011-01-01", freq="D", periods=100),
    )
    ddf = dd.from_pandas(df, 50)
    assert_eq(df.loc["2011-01"], ddf.loc["2011-01"])
    assert_eq(df.loc["2011"], ddf.loc["2011"])

    assert_eq(df.loc["2011-01":"2012-05"], ddf.loc["2011-01":"2012-05"])
    assert_eq(df.loc["2011":"2015"], ddf.loc["2011":"2015"])


def test_getitem_period_str():
    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.period_range("2011-01-01", freq="h", periods=100),
    )
    ddf = dd.from_pandas(df, 10)
    assert_eq(df["2011-01-02":"2011-01-10"], ddf["2011-01-02":"2011-01-10"])
    # same reso, dask result is always DataFrame

    df = pd.DataFrame(
        {"A": np.random.randn(100), "B": np.random.randn(100)},
        index=pd.period_range("2011-01-01", freq="D", periods=100),
    )
    ddf = dd.from_pandas(df, 50)

    assert_eq(df["2011-01":"2012-05"], ddf["2011-01":"2012-05"])
    assert_eq(df["2011":"2015"], ddf["2011":"2015"])


@pytest.mark.parametrize(
    "index",
    [
        pd.date_range("2011-01-01", freq="h", periods=100),  # time index
        range(100),  # numerical index
    ],
)
def test_to_series(index):
    df = pd.DataFrame({"A": np.random.randn(100)}, index=index)
    ddf = dd.from_pandas(df, 10)

    expected = df.index.to_series()
    actual = ddf.index.to_series()

    assert actual.known_divisions
    assert_eq(expected, actual)


@pytest.mark.parametrize(
    "index",
    [
        pd.date_range("2011-01-01", freq="h", periods=100),  # time index
        range(100),  # numerical index
    ],
)
def test_to_frame(index):
    df = pd.DataFrame({"A": np.random.randn(100)}, index=index)
    ddf = dd.from_pandas(df, 10)

    expected = df.index.to_frame()
    actual = ddf.index.to_frame()

    assert actual.known_divisions
    assert_eq(expected, actual)

    # test name option
    assert_eq(df.index.to_frame(name="foo"), ddf.index.to_frame(name="foo"))


@pytest.mark.parametrize("indexer", [0, [0], [0, 1], [1, 0], [False, True, True]])
def test_iloc(indexer):
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    ddf = dd.from_pandas(df, 2)

    result = ddf.iloc[:, indexer]
    expected = df.iloc[:, indexer]

    assert_eq(result, expected)


def test_iloc_series():
    s = pd.Series([1, 2, 3])
    ds = dd.from_pandas(s, 2)
    with pytest.raises(AttributeError):
        ds.iloc[:]


def test_iloc_raises():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    ddf = dd.from_pandas(df, 2)

    with pytest.raises(NotImplementedError):
        ddf.iloc[[0, 1], :]

    with pytest.raises(NotImplementedError):
        ddf.iloc[[0, 1], [0, 1]]

    with pytest.raises(ValueError):
        ddf.iloc[[0, 1], [0, 1], [1, 2]]

    with pytest.raises(IndexError):
        ddf.iloc[:, [5, 6]]


@pytest.mark.xfail(reason="duplicated columns")
def test_iloc_duplicate_columns():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    ddf = dd.from_pandas(df, 2)
    df.columns = ["A", "A", "C"]
    ddf.columns = ["A", "A", "C"]

    selection = ddf.iloc[:, 2]
    # Check that `iloc` is called instead of getitem
    assert any([key.startswith("iloc") for key in selection.dask.layers.keys()])

    select_first = ddf.iloc[:, 1]
    assert_eq(select_first, df.iloc[:, 1])

    select_zeroth = ddf.iloc[:, 0]
    assert_eq(select_zeroth, df.iloc[:, 0])

    select_list_cols = ddf.iloc[:, [0, 2]]
    assert_eq(select_list_cols, df.iloc[:, [0, 2]])

    select_negative = ddf.iloc[:, -1:-3:-1]
    assert_eq(select_negative, df.iloc[:, -1:-3:-1])


def test_iloc_dispatch_to_getitem():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    ddf = dd.from_pandas(df, 2)

    select_first = ddf.iloc[:, 1]
    assert_eq(select_first, df.iloc[:, 1])

    select_zeroth = ddf.iloc[:, 0]
    assert_eq(select_zeroth, df.iloc[:, 0])

    select_list_cols = ddf.iloc[:, [0, 2]]
    assert_eq(select_list_cols, df.iloc[:, [0, 2]])

    select_negative = ddf.iloc[:, -1:-3:-1]
    assert_eq(select_negative, df.iloc[:, -1:-3:-1])


def test_iloc_out_of_order_selection():
    df = pd.DataFrame({"A": [1] * 100, "B": [2] * 100, "C": [3] * 100, "D": [4] * 100})
    ddf = dd.from_pandas(df, 2)
    ddf = ddf[["C", "A", "B"]]
    a = ddf.iloc[:, 0]
    b = ddf.iloc[:, 1]
    c = ddf.iloc[:, 2]

    assert a.name == "C"
    assert b.name == "A"
    assert c.name == "B"

    a1, b1, c1 = dask.compute(a, b, c)

    assert a1.name == "C"
    assert b1.name == "A"
    assert c1.name == "B"


def test_pandas_nullable_boolean_data_type():
    s1 = pd.Series([0, 1, 2])
    s2 = pd.Series([True, False, pd.NA], dtype="boolean")

    ddf1 = dd.from_pandas(s1, npartitions=1)
    ddf2 = dd.from_pandas(s2, npartitions=1)

    assert_eq(ddf1[ddf2], s1[s2])
    assert_eq(ddf1.loc[ddf2], s1.loc[s2])


def test_deterministic_hashing_series():
    obj = pd.Series([0, 1, 2])

    dask_df = dd.from_pandas(obj, npartitions=1)

    ddf1 = dask_df.loc[0:1]
    ddf2 = dask_df.loc[0:1]

    assert tokenize(ddf1) == tokenize(ddf2)

    ddf2 = dask_df.loc[0:2]
    assert ddf1._name != ddf2._name


@pytest.mark.gpu
def test_gpu_loc():
    cudf = pytest.importorskip("cudf")
    cupy = pytest.importorskip("cupy")

    index = [1, 5, 10, 11, 12, 100, 200, 300]
    df = cudf.DataFrame({"a": range(8), "index": index}).set_index("index")
    ddf = dd.from_pandas(df, npartitions=3)
    cdf_index = cudf.Series([1, 100, 300])
    cupy_index = cupy.array([1, 100, 300])

    assert_eq(ddf.loc[cdf_index], df.loc[cupy_index])
    assert_eq(ddf.loc[cupy_index], df.loc[cupy_index])
