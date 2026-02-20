from __future__ import annotations

import contextlib

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
import dask.dataframe as dd
from dask.dataframe._compat import PANDAS_GE_210, PANDAS_GE_300
from dask.dataframe._pyarrow import to_pyarrow_string
from dask.dataframe.utils import assert_eq, pyarrow_strings_enabled


@contextlib.contextmanager
def ensure_removed(obj, attr):
    """Ensure that an attribute added to 'obj' during the test is
    removed when we're done"""
    try:
        yield
    finally:
        try:
            delattr(obj, attr)
        except AttributeError:
            pass
        obj._accessors.discard(attr)


class MyAccessor:
    def __init__(self, obj):
        self.obj = obj
        self.item = "item"

    @property
    def prop(self):
        return self.item

    def method(self):
        return self.item


@pytest.mark.parametrize(
    "obj, registrar",
    [
        (dd.Series, dd.extensions.register_series_accessor),
        (dd.DataFrame, dd.extensions.register_dataframe_accessor),
        (dd.Index, dd.extensions.register_index_accessor),
    ],
)
def test_register(obj, registrar):
    if obj is dd.Index:
        pytest.skip("from_pandas doesn't support Index")
    with ensure_removed(obj, "mine"):
        before = set(dir(obj))
        registrar("mine")(MyAccessor)
        instance = dd.from_pandas(obj._partition_type([], dtype=float), 2)
        assert instance.mine.prop == "item"
        after = set(dir(obj))
        assert (before ^ after) == {"mine"}
        assert "mine" in obj._accessors


def test_accessor_works():
    with ensure_removed(dd.Series, "mine"):
        dd.extensions.register_series_accessor("mine")(MyAccessor)

        a = pd.Series([1, 2])
        b = dd.from_pandas(a, 2)
        assert b.mine.obj is b

        assert b.mine.prop == "item"
        assert b.mine.method() == "item"


@pytest.fixture
def df_ddf():
    import numpy as np

    df = pd.DataFrame(
        {
            "str_col": ["abc", "bcd", "cdef", "DEFG"],
            "int_col": [1, 2, 3, 4],
            "dt_col": np.array(
                [int(1e9), int(1.1e9), int(1.2e9), None], dtype="M8[ns]"
            ),
        },
        index=["E", "f", "g", "h"],
    )

    df["string_col"] = df["str_col"].astype("string")
    df.loc["E", "string_col"] = pd.NA

    ddf = dd.from_pandas(df, 2)

    return df, ddf


@pytest.mark.filterwarnings("ignore:The behavior of DatetimeProperties")
@pytest.mark.xfail(PANDAS_GE_300, reason="divisions are incorrect")
def test_dt_accessor(df_ddf):
    df, ddf = df_ddf

    assert "date" in dir(ddf.dt_col.dt)

    # pandas loses Series.name via datetime accessor
    # see https://github.com/pydata/pandas/issues/10712
    assert_eq(ddf.dt_col.dt.date, df.dt_col.dt.date, check_names=False)
    if PANDAS_GE_210:
        warning_ctx = pytest.warns(FutureWarning, match="will return a Series")
    else:
        warning_ctx = contextlib.nullcontext()
    # to_pydatetime returns a numpy array in pandas, but a Series in dask
    # pandas will start returning a Series with 3.0 as well
    with warning_ctx:
        ddf_result = ddf.dt_col.dt.to_pydatetime()
    with warning_ctx:
        pd_result = pd.Series(
            df.dt_col.dt.to_pydatetime(), index=df.index, dtype=object
        )
    assert_eq(ddf_result, pd_result)

    assert set(ddf.dt_col.dt.date.dask) == set(ddf.dt_col.dt.date.dask)
    assert set(ddf.dt_col.dt.to_pydatetime().dask) == set(
        ddf.dt_col.dt.to_pydatetime().dask
    )


def test_dt_accessor_not_available(df_ddf):
    df, ddf = df_ddf

    # Not available on invalid dtypes
    with pytest.raises(AttributeError) as exc:
        ddf.str_col.dt
    assert ".dt accessor" in str(exc.value)


def test_str_accessor(df_ddf):
    df, ddf = df_ddf

    # implemented methods are present in tab completion
    assert "upper" in dir(ddf.str_col.str)
    assert "upper" in dir(ddf.string_col.str)
    assert "upper" in dir(ddf.index.str)

    # not implemented methods don't show up
    assert "get_dummies" not in dir(ddf.str_col.str)
    assert not hasattr(ddf.str_col.str, "get_dummies")

    # Test simple method on both series and index
    assert_eq(ddf.str_col.str.upper(), df.str_col.str.upper())
    assert set(ddf.str_col.str.upper().dask) == set(ddf.str_col.str.upper().dask)

    assert_eq(ddf.string_col.str.upper(), df.string_col.str.upper())
    assert set(ddf.string_col.str.upper().dask) == set(ddf.string_col.str.upper().dask)

    assert_eq(ddf.index.str.upper(), df.index.str.upper())
    assert set(ddf.index.str.upper().dask) == set(ddf.index.str.upper().dask)

    # make sure to pass through args & kwargs
    # NOTE: when using pyarrow strings, `.str.contains(...)` will return a result
    # with `boolean` dtype, while using object strings returns a `bool`. We cast
    # the pandas DataFrame here to ensure pandas and Dask return the same dtype.
    ctx = contextlib.nullcontext()
    if pyarrow_strings_enabled():
        df.str_col = to_pyarrow_string(df.str_col)
        if not PANDAS_GE_210:
            ctx = pytest.warns(
                pd.errors.PerformanceWarning, match="Falling back on a non-pyarrow"
            )
    assert_eq(
        ddf.str_col.str.contains("a"),
        df.str_col.str.contains("a"),
    )
    assert_eq(ddf.string_col.str.contains("a"), df.string_col.str.contains("a"))
    assert set(ddf.str_col.str.contains("a").dask) == set(
        ddf.str_col.str.contains("a").dask
    )

    with ctx:
        expected = df.str_col.str.contains("d", case=False)
        assert_eq(
            ddf.str_col.str.contains("d", case=False),
            expected,
        )
        assert set(ddf.str_col.str.contains("d", case=False).dask) == set(
            ddf.str_col.str.contains("d", case=False).dask
        )

    for na in [True, False]:
        assert_eq(
            ddf.str_col.str.contains("a", na=na),
            df.str_col.str.contains("a", na=na),
        )
        assert set(ddf.str_col.str.contains("a", na=na).dask) == set(
            ddf.str_col.str.contains("a", na=na).dask
        )

    for regex in [True, False]:
        assert_eq(
            ddf.str_col.str.contains("a", regex=regex),
            df.str_col.str.contains("a", regex=regex),
        )
        assert set(ddf.str_col.str.contains("a", regex=regex).dask) == set(
            ddf.str_col.str.contains("a", regex=regex).dask
        )


def test_str_accessor_not_available(df_ddf):
    df, ddf = df_ddf

    # Not available on invalid dtypes
    with pytest.raises(AttributeError) as exc:
        ddf.int_col.str
    assert ".str accessor" in str(exc.value)

    assert "str" not in dir(ddf.int_col)


def test_str_accessor_getitem(df_ddf):
    df, ddf = df_ddf
    assert_eq(ddf.str_col.str[:2], df.str_col.str[:2])
    assert_eq(ddf.str_col.str[1], df.str_col.str[1])


def test_str_accessor_extractall(df_ddf):
    df, ddf = df_ddf
    assert_eq(
        ddf.str_col.str.extractall("(.*)b(.*)"), df.str_col.str.extractall("(.*)b(.*)")
    )


@pytest.mark.parametrize("method", ["removeprefix", "removesuffix"])
def test_str_accessor_removeprefix_removesuffix(df_ddf, method):
    df, ddf = df_ddf
    prefix = df.str_col.iloc[0][:2]
    suffix = df.str_col.iloc[0][-2:]
    missing = "definitely a missing prefix/suffix"

    def call(df, arg):
        return getattr(df.str_col.str, method)(arg)

    assert_eq(call(ddf, prefix), call(df, prefix))
    assert_eq(call(ddf, suffix), call(df, suffix))
    assert_eq(call(ddf, missing), call(df, missing))


def test_str_accessor_cat(df_ddf):
    df, ddf = df_ddf
    sol = df.str_col.str.cat(df.str_col.str.upper(), sep=":")
    assert_eq(ddf.str_col.str.cat(ddf.str_col.str.upper(), sep=":"), sol)
    assert_eq(ddf.str_col.str.cat(df.str_col.str.upper(), sep=":"), sol)
    assert_eq(
        ddf.str_col.str.cat([ddf.str_col.str.upper(), df.str_col.str.lower()], sep=":"),
        df.str_col.str.cat([df.str_col.str.upper(), df.str_col.str.lower()], sep=":"),
    )
    assert_eq(ddf.str_col.str.cat(sep=":"), df.str_col.str.cat(sep=":"))

    for o in ["foo", ["foo"]]:
        with pytest.raises(TypeError):
            ddf.str_col.str.cat(o)


def test_str_accessor_cat_none():
    s = pd.Series(["a", "a", "b", "b", "c", np.nan], name="foo")
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(ds.str.cat(), s.str.cat())
    assert_eq(ds.str.cat(na_rep="-"), s.str.cat(na_rep="-"))
    assert_eq(ds.str.cat(sep="_", na_rep="-"), s.str.cat(sep="_", na_rep="-"))


@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_str_accessor_split_noexpand(method):
    def call(obj, *args, **kwargs):
        return getattr(obj.str, method)(*args, **kwargs)

    s = pd.Series(["a b c d", "aa bb cc dd", "aaa bbb ccc dddd"], name="foo")
    ds = dd.from_pandas(s, npartitions=2)

    for n in [1, 2, 3]:
        assert_eq(call(s, n=n, expand=False), call(ds, n=n, expand=False))

    assert call(ds, n=1, expand=False).name == "foo"


@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_str_accessor_split_expand(method):
    def call(obj, *args, **kwargs):
        return getattr(obj.str, method)(*args, **kwargs)

    s = pd.Series(
        ["a b c d", "aa bb cc dd", "aaa bbb ccc dddd"], index=["row1", "row2", "row3"]
    )
    ds = dd.from_pandas(s, npartitions=2)

    for n in [1, 2, 3]:
        assert_eq(call(s, n=n, expand=True), call(ds, n=n, expand=True))

    with pytest.raises(NotImplementedError) as info:
        call(ds, expand=True)

    assert "n=" in str(info.value)

    s = pd.Series(["a,bcd,zz,f", "aabb,ccdd,z,kk", "aaabbb,cccdddd,l,pp"])
    ds = dd.from_pandas(s, npartitions=2)

    for n in [1, 2, 3]:
        assert_eq(
            call(s, pat=",", n=n, expand=True), call(ds, pat=",", n=n, expand=True)
        )


@pytest.mark.xfail(reason="Need to pad columns")
def test_str_accessor_split_expand_more_columns():
    s = pd.Series(["a b c d", "aa", "aaa bbb ccc dddd"])
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(s.str.split(n=3, expand=True), ds.str.split(n=3, expand=True))

    s = pd.Series(["a b c", "aa bb cc", "aaa bbb ccc"])
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(
        ds.str.split(n=10, expand=True),
        s.str.split(n=10, expand=True),
    )


@pytest.mark.parametrize("index", [None, [0]], ids=["range_index", "other index"])
def test_str_split_no_warning(index):
    df = pd.DataFrame({"a": ["a\nb"]}, index=index)
    ddf = dd.from_pandas(df, npartitions=1)

    pd_a = df["a"].str.split("\n", n=1, expand=True)
    dd_a = ddf["a"].str.split("\n", n=1, expand=True)

    assert_eq(dd_a, pd_a)


def test_string_nullable_types(df_ddf):
    df, ddf = df_ddf
    assert_eq(ddf.string_col.str.count("A"), df.string_col.str.count("A"))
    assert_eq(ddf.string_col.str.isalpha(), df.string_col.str.isalpha())
