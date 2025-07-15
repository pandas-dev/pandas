from __future__ import annotations

import re
import warnings
from collections.abc import Iterable

import numpy as np
import pandas as pd
import pytest
from packaging.version import Version

import dask
import dask.dataframe as dd
from dask.dataframe._compat import PANDAS_GE_300, tm
from dask.dataframe.core import apply_and_enforce
from dask.dataframe.utils import (
    UNKNOWN_CATEGORIES,
    assert_eq,
    check_matching_columns,
    check_meta,
    is_dataframe_like,
    is_index_like,
    is_series_like,
    make_meta,
    meta_frame_constructor,
    meta_nonempty,
    meta_series_constructor,
    pyarrow_strings_enabled,
    raise_on_meta_error,
    valid_divisions,
)
from dask.local import get_sync


def test_make_meta():
    df = pd.DataFrame(
        {"a": [1, 2, 3], "b": list("abc"), "c": [1.0, 2.0, 3.0]}, index=[10, 20, 30]
    )

    # Pandas dataframe
    meta = make_meta(df)
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes).all()
    assert isinstance(meta.index, type(df.index))
    # - ensure no references to original data arrays are kept
    for col in "abc":
        meta_pointer = meta[col].values.__array_interface__["data"][0]
        df_pointer = df[col].values.__array_interface__["data"][0]
        assert meta_pointer != df_pointer
    meta_pointer = meta.index.values.__array_interface__["data"][0]
    df_pointer = df.index.values.__array_interface__["data"][0]
    assert meta_pointer != df_pointer

    # Pandas series
    meta = make_meta(df.a)
    assert len(meta) == 0
    assert meta.dtype == df.a.dtype
    assert isinstance(meta.index, type(df.index))
    # - ensure no references to original data arrays are kept
    meta_pointer = meta.values.__array_interface__["data"][0]
    df_pointer = df.a.values.__array_interface__["data"][0]
    assert meta_pointer != df_pointer
    meta_pointer = meta.index.values.__array_interface__["data"][0]
    df_pointer = df.index.values.__array_interface__["data"][0]
    assert meta_pointer != df_pointer

    # Pandas index
    meta = make_meta(df.index)
    assert isinstance(meta, type(df.index))
    assert len(meta) == 0
    # - ensure no references to original data arrays are kept
    meta_pointer = meta.values.__array_interface__["data"][0]
    df_pointer = df.index.values.__array_interface__["data"][0]
    assert meta_pointer != df_pointer

    # Dask object
    ddf = dd.from_pandas(df, npartitions=2)
    assert make_meta(ddf) is ddf._meta

    # Dict
    meta = make_meta({"a": "i8", "b": "O", "c": "f8"})
    assert isinstance(meta, pd.DataFrame)
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes).all()
    assert isinstance(meta.index, pd.RangeIndex)

    # List
    meta = make_meta([("a", "i8"), ("c", "f8"), ("b", "O")])
    assert (meta.columns == ["a", "c", "b"]).all()
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes[meta.dtypes.index]).all()
    assert isinstance(meta.index, pd.RangeIndex)

    # Tuple
    meta = make_meta(("a", "i8"))
    assert isinstance(meta, pd.Series)
    assert len(meta) == 0
    assert meta.dtype == "i8"
    assert meta.name == "a"

    # Iterable
    class CustomMetadata(Iterable):
        """Custom class iterator returning pandas types."""

        def __init__(self, max=0):
            self.types = [("a", "i8"), ("c", "f8"), ("b", "O")]

        def __iter__(self):
            self.n = 0
            return self

        def __next__(self):
            if self.n < len(self.types):
                ret = self.types[self.n]
                self.n += 1
                return ret
            else:
                raise StopIteration

    meta = make_meta(CustomMetadata())
    assert (meta.columns == ["a", "c", "b"]).all()
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes[meta.dtypes.index]).all()
    assert isinstance(meta.index, pd.RangeIndex)

    # With index
    idx = pd.Index([1, 2], name="foo")
    meta = make_meta(
        {"a": "i8", "b": "i4"},
        index=idx,
    )
    assert type(meta.index) is type(idx)
    assert meta.index.dtype == "int64"
    assert len(meta.index) == 0

    meta = make_meta(("a", "i8"), index=idx)
    assert type(meta.index) is type(idx)
    assert meta.index.dtype == "int64"
    assert len(meta.index) == 0

    # Categoricals
    meta = make_meta({"a": "category"}, parent_meta=df)
    assert len(meta.a.cat.categories) == 1
    assert meta.a.cat.categories[0] == UNKNOWN_CATEGORIES
    meta = make_meta(("a", "category"), parent_meta=df)
    assert len(meta.cat.categories) == 1
    assert meta.cat.categories[0] == UNKNOWN_CATEGORIES

    # Categorials with Index
    meta = make_meta({"a": "category", "b": "int64"}, index=idx)
    assert len(meta.a.cat.categories) == 1
    assert meta.index.dtype == "int64"
    assert meta.index.empty

    # Numpy scalar
    meta = make_meta(np.float64(1.0), parent_meta=df)
    assert isinstance(meta, np.float64)

    # Python scalar
    meta = make_meta(1.0, parent_meta=df)
    assert isinstance(meta, np.float64)

    # Timestamp
    x = pd.Timestamp(2000, 1, 1)
    meta = make_meta(x, parent_meta=df)
    assert meta is x

    # DatetimeTZDtype
    x = pd.DatetimeTZDtype(tz="UTC")
    meta = make_meta(x)
    assert meta == pd.Timestamp(1, tz=x.tz, unit=x.unit)

    # Dtype expressions
    meta = make_meta("i8", parent_meta=df)
    assert isinstance(meta, np.int64)
    meta = make_meta(float, parent_meta=df)
    assert isinstance(meta, np.dtype(float).type)
    meta = make_meta(np.dtype("bool"), parent_meta=df)
    assert isinstance(meta, np.bool_)
    assert pytest.raises(TypeError, lambda: make_meta(None))


def test_meta_nonempty():
    df1 = pd.DataFrame(
        {
            "A": pd.Categorical(["Alice", "Bob", "Carol"]),
            "B": list("abc"),
            "C": "bar",
            "D": np.float32(1),
            "E": np.int32(1),
            "F": pd.Timestamp("2016-01-01"),
            "G": pd.date_range("2016-01-01", periods=3, tz="America/New_York"),
            "H": pd.Timedelta("1 hours"),
            "I": np.void(b" "),
            "J": pd.Categorical([UNKNOWN_CATEGORIES] * 3),
            "K": pd.Categorical([None, None, None]),
        },
        columns=list("DCBAHGFEIJK"),
    )
    df2 = df1.iloc[0:0]
    df3 = meta_nonempty(df2)
    assert (df3.dtypes == df2.dtypes).all()
    assert df3["A"][0] == "Alice"
    assert df3["B"][0] == "foo"
    assert df3["C"][0] == "foo"
    assert df3["D"][0] == np.float32(1)
    assert df3["D"][0].dtype == "f4"
    assert df3["E"][0] == np.int32(1)
    assert df3["E"][0].dtype == "i4"
    assert df3["F"][0] == pd.Timestamp("1970-01-01 00:00:00")
    assert df3["G"][0] == pd.Timestamp("1970-01-01 00:00:00", tz="America/New_York")
    assert df3["H"][0] == pd.Timedelta("1")
    assert df3["I"][0] == "foo"
    assert df3["J"][0] == UNKNOWN_CATEGORIES
    assert len(df3["K"].cat.categories) == 0

    s = meta_nonempty(df2["A"])
    assert s.dtype == df2["A"].dtype
    assert (df3["A"] == s).all()


def test_meta_duplicated():
    df = pd.DataFrame(columns=["A", "A", "B"])
    res = meta_nonempty(df)

    exp = pd.DataFrame(
        [["foo", "foo", "foo"], ["foo", "foo", "foo"]],
        index=meta_nonempty(df.index),
        columns=["A", "A", "B"],
    )
    tm.assert_frame_equal(res, exp)


def test_meta_nonempty_empty_categories():
    for dtype in ["O", "f8", "M8[ns]"]:
        # Index
        idx = pd.CategoricalIndex(
            [], pd.Index([], dtype=dtype), ordered=True, name="foo"
        )
        res = meta_nonempty(idx)
        assert type(res) is pd.CategoricalIndex
        assert type(res.categories) is type(idx.categories)
        assert res.ordered == idx.ordered
        assert res.name == idx.name
        # Series
        s = idx.to_series()
        res = meta_nonempty(s)
        assert res.dtype == "category"
        assert s.dtype == "category"
        assert type(res.cat.categories) is type(s.cat.categories)
        assert res.cat.ordered == s.cat.ordered
        assert res.name == s.name


def test_meta_nonempty_index():
    idx = pd.RangeIndex(1, name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.RangeIndex
    assert res.name == idx.name

    idx = pd.Index([1], name="foo", dtype="int")
    res = meta_nonempty(idx)
    assert type(res) is type(idx)
    assert res.dtype == np.int_
    assert res.name == idx.name

    idx = pd.Index(["a"], name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.Index
    assert res.name == idx.name

    idx = pd.DatetimeIndex(["1970-01-01"], freq="D", tz="America/New_York", name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.DatetimeIndex
    assert res.tz == idx.tz
    assert res.freq == idx.freq
    assert res.name == idx.name

    idx = pd.PeriodIndex(["1970-01-01"], freq="D", name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.PeriodIndex
    assert res.freq == idx.freq
    assert res.name == idx.name

    idx = pd.TimedeltaIndex([pd.Timedelta(1, "D")], freq="D", name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.TimedeltaIndex
    assert res.freq == idx.freq
    assert res.name == idx.name

    idx = pd.CategoricalIndex(["xyx"], ["xyx", "zzz"], ordered=True, name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.CategoricalIndex
    assert (res.categories == idx.categories).all()
    assert res.ordered == idx.ordered
    assert res.name == idx.name

    idx = pd.CategoricalIndex([], [UNKNOWN_CATEGORIES], ordered=True, name="foo")
    res = meta_nonempty(idx)
    assert type(res) is pd.CategoricalIndex
    assert res.ordered == idx.ordered
    assert res.name == idx.name

    levels = [pd.Index([1], name="a"), pd.Index([1.0], name="b")]
    codes = [[0], [0]]
    idx = pd.MultiIndex(levels=levels, names=["a", "b"], codes=codes)
    res = meta_nonempty(idx)
    assert type(res) is pd.MultiIndex
    for idx1, idx2 in zip(idx.levels, res.levels):
        assert type(idx1) is type(idx2)
        assert idx1.name == idx2.name
    assert res.names == idx.names

    levels = [
        pd.Index([1], name="a"),
        pd.CategoricalIndex(data=["xyx"], categories=["xyx"], name="b"),
        pd.TimedeltaIndex([np.timedelta64(1, "D")], name="timedelta"),
    ]

    codes = [[0], [0], [0]]

    idx = pd.MultiIndex(levels=levels, names=["a", "b", "timedelta"], codes=codes)
    res = meta_nonempty(idx)
    assert type(res) is pd.MultiIndex
    for idx1, idx2 in zip(idx.levels, res.levels):
        assert type(idx1) is type(idx2)
        assert idx1.name == idx2.name
    assert res.names == idx.names


def test_meta_nonempty_uint64index():
    idx = pd.Index([1], name="foo", dtype="uint64")
    res = meta_nonempty(idx)
    assert type(res) is type(idx)
    assert res.dtype == "uint64"
    assert res.name == idx.name


def test_meta_nonempty_scalar():
    meta = meta_nonempty(np.float64(1.0))
    assert isinstance(meta, np.float64)

    x = pd.Timestamp(2000, 1, 1)
    meta = meta_nonempty(x)
    assert meta is x

    # DatetimeTZDtype
    x = pd.DatetimeTZDtype(tz="UTC")
    meta = meta_nonempty(x)
    assert meta == pd.Timestamp(1, tz=x.tz, unit=x.unit)


def test_raise_on_meta_error():
    try:
        with raise_on_meta_error():
            raise RuntimeError("Bad stuff")
    except Exception as e:
        assert e.args[0].startswith("Metadata inference failed.\n")
        assert "RuntimeError" in e.args[0]
    else:
        assert False, "should have errored"

    try:
        with raise_on_meta_error("myfunc"):
            raise RuntimeError("Bad stuff")
    except Exception as e:
        assert e.args[0].startswith("Metadata inference failed in `myfunc`.\n")
        assert "RuntimeError" in e.args[0]
    else:
        assert False, "should have errored"


def test_check_meta():
    df = pd.DataFrame(
        {
            "a": ["x", "y", "z"],
            "b": [True, False, True],
            "c": [1, 2.5, 3.5],
            "d": [1, 2, 3],
            "e": pd.Categorical(["x", "y", "z"]),
            "f": pd.Series([1, 2, 3], dtype=np.uint64),
        }
    )
    meta = df.iloc[:0]

    # DataFrame metadata passthrough if correct
    assert check_meta(df, meta) is df
    # Series metadata passthrough if correct
    e = df.e
    assert check_meta(e, meta.e) is e
    # numeric_equal means floats and ints are equivalent
    d = df.d
    f = df.f
    assert check_meta(d, meta.d.astype("f8"), numeric_equal=True) is d
    assert check_meta(f, meta.f.astype("f8"), numeric_equal=True) is f
    assert check_meta(f, meta.f.astype("i8"), numeric_equal=True) is f

    # Series metadata error
    with pytest.raises(ValueError) as err:
        check_meta(d, meta.d.astype("f8"), numeric_equal=False)
    series = "pandas.core.series.Series" if not PANDAS_GE_300 else "pandas.Series"
    assert str(err.value) == (
        "Metadata mismatch found.\n"
        "\n"
        f"Partition type: `{series}`\n"
        "+----------+---------+\n"
        "|          | dtype   |\n"
        "+----------+---------+\n"
        "| Found    | int64   |\n"
        "| Expected | float64 |\n"
        "+----------+---------+"
    )

    # DataFrame metadata error
    meta2 = meta.astype({"a": "category", "d": "f8"})[["a", "b", "c", "d"]]
    df2 = df[["a", "b", "d", "e"]]
    with pytest.raises(ValueError) as err:
        check_meta(df2, meta2, funcname="from_delayed")
    frame = "pandas.core.frame.DataFrame" if not PANDAS_GE_300 else "pandas.DataFrame"

    exp = (
        "Metadata mismatch found in `from_delayed`.\n"
        "\n"
        f"Partition type: `{frame}`\n"
        "+--------+----------+----------+\n"
        "| Column | Found    | Expected |\n"
        "+--------+----------+----------+\n"
        "| 'a'    | object   | category |\n"
        "| 'c'    | -        | float64  |\n"
        "| 'e'    | category | -        |\n"
        "+--------+----------+----------+"
    )
    assert str(err.value) == exp

    # pandas dtype metadata error
    with pytest.raises(ValueError) as err:
        check_meta(df.a, pd.Series([], dtype="string"), numeric_equal=False)
    assert str(err.value) == (
        "Metadata mismatch found.\n"
        "\n"
        f"Partition type: `{series}`\n"
        "+----------+--------+\n"
        "|          | dtype  |\n"
        "+----------+--------+\n"
        "| Found    | object |\n"
        "| Expected | string |\n"
        "+----------+--------+"
    )


def test_check_matching_columns_raises_appropriate_errors():
    df = pd.DataFrame(columns=["a", "b", "c"])

    meta = pd.DataFrame(columns=["b", "a", "c"])
    with pytest.raises(
        ValueError,
        match="Order of columns does not match."
        "\nActual:   \\['a', 'b', 'c'\\]"
        "\nExpected: \\['b', 'a', 'c'\\]",
    ):
        assert check_matching_columns(meta, df)

    meta = pd.DataFrame(columns=["a", "b", "c", "d"])
    with pytest.raises(ValueError, match="Missing: \\['d'\\]"):
        assert check_matching_columns(meta, df)

    meta = pd.DataFrame(columns=["a", "b"])
    with pytest.raises(ValueError, match="Extra:   \\['c'\\]"):
        assert check_matching_columns(meta, df)


def test_check_meta_typename():
    df = pd.DataFrame({"x": []})
    ddf = dd.from_pandas(df, npartitions=1)
    check_meta(df, df)
    with pytest.raises(Exception) as info:
        check_meta(ddf, df)

    assert "dask" in str(info.value)
    assert "pandas" in str(info.value)


@pytest.mark.parametrize("frame_value_counts", [True, False])
def test_is_dataframe_like(monkeypatch, frame_value_counts):
    # When we drop support for pandas 1.0, this compat check can
    # be dropped
    if frame_value_counts:
        monkeypatch.setattr(pd.DataFrame, "value_counts", lambda x: None, raising=False)

    df = pd.DataFrame({"x": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=1)

    assert is_dataframe_like(df)
    assert is_dataframe_like(ddf)
    assert not is_dataframe_like(df.x)
    assert not is_dataframe_like(ddf.x)
    assert not is_dataframe_like(df.index)
    assert not is_dataframe_like(ddf.index)
    assert not is_dataframe_like(pd.DataFrame)

    assert not is_series_like(df)
    assert not is_series_like(ddf)
    assert is_series_like(df.x)
    assert is_series_like(ddf.x)
    assert not is_series_like(df.index)
    assert not is_series_like(ddf.index)
    assert not is_series_like(pd.Series)

    assert not is_index_like(df)
    assert not is_index_like(ddf)
    assert not is_index_like(df.x)
    assert not is_index_like(ddf.x)
    assert is_index_like(df.index)
    assert is_index_like(ddf.index)
    assert not is_index_like(pd.Index)

    # The following checks support of class wrappers, which
    # requires the comparisons of `x.__class__` instead of `type(x)`
    class DataFrameWrapper:
        __class__ = pd.DataFrame

    wrap = DataFrameWrapper()
    wrap.dtypes = None
    wrap.columns = None
    assert is_dataframe_like(wrap)

    class SeriesWrapper:
        __class__ = pd.Series

    wrap = SeriesWrapper()
    wrap.dtype = None
    wrap.name = None
    assert is_series_like(wrap)

    class IndexWrapper:
        __class__ = pd.Index

    wrap = IndexWrapper()
    wrap.dtype = None
    wrap.name = None
    assert is_index_like(wrap)


def test_apply_and_enforce_message():
    def func():
        return pd.DataFrame(columns=["A", "B", "C"], index=[0])

    meta = pd.DataFrame(columns=["A", "D"], index=[0])
    with pytest.raises(ValueError, match="Extra: *['B', 'C']"):
        apply_and_enforce(_func=func, _meta=meta)

    with pytest.raises(ValueError, match=re.escape("Missing: ['D']")):
        apply_and_enforce(_func=func, _meta=meta)


def test_nonempty_series_sparse():
    ser = pd.Series(pd.array([0, 1], dtype="Sparse"))
    with warnings.catch_warnings(record=True) as record:
        meta_nonempty(ser)
    assert not record


def test_nonempty_series_nullable_float():
    ser = pd.Series([], dtype="Float64")
    non_empty = meta_nonempty(ser)
    assert non_empty.dtype == "Float64"


def test_assert_eq_sorts():
    df = pd.DataFrame({"A": np.linspace(0, 1, 10), "B": np.random.random(10)})
    df_s = df.sort_values("B")
    assert_eq(df, df_s)
    with pytest.raises(AssertionError):
        assert_eq(df, df_s, sort_results=False)

    df_sr = df_s.reset_index(drop=True)
    assert_eq(df, df_sr, check_index=False)
    with pytest.raises(AssertionError):
        assert_eq(df, df_sr)
    with pytest.raises(AssertionError):
        assert_eq(df, df_sr, check_index=False, sort_results=False)

    ddf = dd.from_pandas(df, npartitions=2)
    ddf_s = ddf.sort_values(["B"])
    assert_eq(df, ddf_s)
    with pytest.raises(AssertionError):
        assert_eq(df, ddf_s, sort_results=False)

    ddf_sr = ddf_s.reset_index(drop=True)
    assert_eq(df, ddf_sr, check_index=False)
    with pytest.raises(AssertionError):
        assert_eq(df, ddf_sr, check_index=False, sort_results=False)


def test_assert_eq_scheduler():
    using_custom_scheduler = False

    def custom_scheduler(*args, **kwargs):
        nonlocal using_custom_scheduler
        try:
            using_custom_scheduler = True
            return get_sync(*args, **kwargs)
        finally:
            using_custom_scheduler = False

    def check_custom_scheduler(part: pd.DataFrame) -> pd.DataFrame:
        assert using_custom_scheduler, "not using custom scheduler"
        return part + 1

    df = pd.DataFrame({"x": [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = ddf.map_partitions(check_custom_scheduler, meta=ddf)

    with pytest.raises(AssertionError, match="not using custom scheduler"):
        # NOTE: we compare `ddf2` to itself in order to test both sides of the `assert_eq` logic.
        assert_eq(ddf2, ddf2)

    assert_eq(ddf2, ddf2, scheduler=custom_scheduler)
    with dask.config.set(scheduler=custom_scheduler):
        assert_eq(ddf2, ddf2, scheduler=None)


@pytest.mark.parametrize(
    "data",
    [
        pd.DataFrame([0]),
        pd.Series([0]),
        pd.Index([0]),
        dd.from_dict({"x": [0]}, npartitions=1),
        dd.from_dict({"x": [0]}, npartitions=1).x,
        dd.from_dict({"x": [0]}, npartitions=1).index,
    ],
)
def test_meta_constructor_utilities(data):
    assert meta_series_constructor(data) is pd.Series
    assert meta_frame_constructor(data) is pd.DataFrame


@pytest.mark.parametrize(
    "data",
    [
        dd.from_dict({"x": [0]}, npartitions=1).x.values,
        np.array([0]),
    ],
)
def test_meta_constructor_utilities_raise(data):
    with pytest.raises(TypeError, match="not supported by meta_series"):
        meta_series_constructor(data)
    with pytest.raises(TypeError, match="not supported by meta_frame"):
        meta_frame_constructor(data)


@pytest.mark.parametrize(
    "divisions, valid",
    [
        ([1, 2, 3], True),
        ([3, 2, 1], False),
        ([1, 1, 1], False),
        ([0, 1, 1], True),
        ((1, 2, 3), True),
        (123, False),
        ([0, float("nan"), 1], False),
    ],
)
def test_valid_divisions(divisions, valid):
    assert valid_divisions(divisions) == valid


def test_pyarrow_strings_enabled():
    try:
        import pyarrow as pa
    except ImportError:
        pa = None

    # If `pyarrow>=12` are installed, then default to using pyarrow strings
    if (
        dask.config.get("dataframe.convert-string") in (True, None)
        and pa is not None
        and Version(pa.__version__) >= Version("12.0.0")
    ):
        assert pyarrow_strings_enabled() is True
    else:
        assert pyarrow_strings_enabled() is False

    # Regardless of dependencies that are installed, always obey
    # the `dataframe.convert-string` config value if it's specified
    with dask.config.set({"dataframe.convert-string": False}):
        assert pyarrow_strings_enabled() is False

    with dask.config.set({"dataframe.convert-string": True}):
        assert pyarrow_strings_enabled() is True
