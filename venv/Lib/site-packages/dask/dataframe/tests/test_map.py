from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import dask.dataframe as dd
from dask.dataframe._compat import (
    PANDAS_GE_210,
    PANDAS_GE_220,
    PANDAS_GE_300,
    PANDAS_GE_310,
)
from dask.dataframe.utils import assert_eq


@pytest.mark.parametrize("base_npart", [1, 4])
@pytest.mark.parametrize("map_npart", [1, 3])
@pytest.mark.parametrize("sorted_index", [False, True])
@pytest.mark.parametrize("sorted_map_index", [False, True])
def test_series_map(base_npart, map_npart, sorted_index, sorted_map_index):
    base = pd.Series(np.arange(100, dtype="i2"))
    if not sorted_index:
        index = np.arange(100)
        np.random.shuffle(index)
        base.index = index
    map_index = np.arange(0, 1000, 10, dtype="i2")
    mapper = pd.Series(np.random.randint(50, size=len(map_index)), index=map_index)
    if not sorted_map_index:
        map_index = np.array(map_index)
        np.random.shuffle(map_index)
        mapper.index = map_index
    expected = base.map(mapper)
    dask_base = dd.from_pandas(base, npartitions=base_npart, sort=False)
    dask_map = dd.from_pandas(mapper, npartitions=map_npart, sort=False)
    with pytest.warns(UserWarning, match="meta"):
        result = dask_base.map(dask_map)
    assert_eq(expected, result)


def test_series_map2():
    df = pd.DataFrame(
        {"a": range(9), "b": [4, 5, 6, 1, 2, 3, 0, 0, 0]},
        index=pd.Index([0, 1, 3, 5, 6, 8, 9, 9, 9], name="myindex"),
    )
    ddf = dd.from_pandas(df, npartitions=3)
    with pytest.warns(UserWarning, match="meta"):
        assert_eq(ddf.a.map(lambda x: x + 1), df.a.map(lambda x: x + 1))
        lk = {v: v + 1 for v in df.a.values}
        assert_eq(ddf.a.map(lk), df.a.map(lk))
        assert_eq(ddf.b.map(lk), df.b.map(lk))
        lk = pd.Series(lk)
        assert_eq(ddf.a.map(lk), df.a.map(lk))
        assert_eq(ddf.b.map(lk), df.b.map(lk))
    assert_eq(ddf.b.map(lk, meta=ddf.b), df.b.map(lk))
    assert_eq(ddf.b.map(lk, meta=("b", "i8")), df.b.map(lk))


def test_series_map_meta():
    ser = pd.Series(np.arange(100), dtype="i2")
    mapper = pd.Series(np.random.randint(50, size=len(ser), dtype="u1"))
    expected = ser.map(mapper)
    dask_base = dd.from_pandas(ser, npartitions=5)
    dask_map = dd.from_pandas(mapper, npartitions=5)
    with pytest.warns(UserWarning):
        result = dask_base.map(dask_map)
    assert_eq(expected, result)


def test_series_map_meta_2():
    pdf = pd.DataFrame({"x": range(100)})
    df = dd.from_pandas(pdf, npartitions=10)

    expected = pdf.x.map(lambda x: x + 1)
    result = df.x.map(lambda x: x + 1, meta=expected.iloc[:0])
    assert_eq(result, expected)

    result = df.x.map(df.x + 1, meta=expected.iloc[:0])
    assert_eq(result, expected)

    result = df.x.map(df.x + 1, meta=("x", "int64"))
    assert_eq(result, expected)

    result = df.x.map((df.x + 1).compute(), meta=("x", "int64"))
    assert_eq(result, expected)

    pdf.index.name = "a"
    df = dd.from_pandas(pdf, npartitions=10)
    expected = pdf.x.map(lambda x: x + 1)
    result = df.x.map(lambda x: x + 1, meta=("x", "int64"))
    assert_eq(result, expected)


def test_index_map():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert ddf.known_divisions is True

    with pytest.warns(UserWarning):
        cleared = ddf.index.map(lambda x: x * 10)
    assert cleared.known_divisions is False

    with pytest.warns(UserWarning):
        applied = ddf.index.map(lambda x: x * 10, is_monotonic=True)
    assert applied.known_divisions is True
    assert applied.divisions == tuple(x * 10 for x in ddf.divisions)


@pytest.mark.skipif(not PANDAS_GE_210, reason="Not available before")
@pytest.mark.parametrize("na_action", [None, "ignore"])
def test_dataframe_map(na_action):
    df = pd.DataFrame({"x": [1, 2, 3, np.nan], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)
    with pytest.warns(UserWarning, match="meta"):
        assert_eq(
            ddf.map(lambda x: x + 1, na_action=na_action),
            df.map(lambda x: x + 1, na_action=na_action),
        )
        assert_eq(ddf.map(lambda x: (x, x)), df.map(lambda x: (x, x)))


def test_series_meta_raises():
    # Raise when we use a user defined function
    s = pd.Series(["abcd", "abcd"])
    ds = dd.from_pandas(s, npartitions=2)
    try:
        ds.map(lambda x: x[3])
    except ValueError as e:
        assert "meta=" in str(e)


@pytest.mark.skipif(PANDAS_GE_210, reason="Added in pandas 2.1.0")
def test_dataframe_map_raises():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)
    with pytest.raises(NotImplementedError, match="DataFrame.map requires pandas"):
        ddf.map(lambda x: x + 1)


def test_map_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert ddf.known_divisions is True
    with pytest.warns(UserWarning, match="meta"):
        cleared = ddf.index.map(lambda x: x * 10)
    assert cleared.known_divisions is False

    with pytest.warns(UserWarning, match="meta"):
        applied = ddf.index.map(lambda x: x * 10, is_monotonic=True)
    assert applied.known_divisions is True
    assert applied.divisions == tuple(x * 10 for x in ddf.divisions)


@pytest.mark.xfail(
    reason="Co-aligned Series.map(Series) lookup semantics: each partition "
    "only sees its mapper slice, not the full lookup table",
)
def test_series_map_co_aligned_lookup():
    base = pd.Series([1, 2])
    dbase = dd.from_pandas(base, npartitions=2)
    expected = base.map(base)
    result = dbase.map(dbase, meta=expected)
    assert_eq(result, expected)


def assert_series_map_dtype(base: pd.Series | pd.Index, func: pd.Series) -> None:
    """Test that Series.map() behaves identically in Pandas and Dask
    when `func` is a Series, and that Dask's behaviour is coherent in the edge cases
    where the output dtype changes dynamically.
    """
    if isinstance(base, pd.Index):
        dbase = dd.from_pandas(pd.Series(base, index=base), npartitions=2).index
        assert_eq(dbase, base)
    else:
        dbase = dd.from_pandas(base, npartitions=2)
    dfunc = dd.from_pandas(func, npartitions=2)

    expect = base.map(func)

    # Note: expect.dtype != func.dtype in edge cases
    result = dbase.map(dfunc, meta=expect)

    # Note that this will concatenate partitions which may have
    # mismatched dtypes
    assert_eq(result, expect)


series_map_dtypes = pytest.mark.parametrize(
    "values,dtype",
    [
        # NumPy dtypes
        ([10, 20], np.int16),
        ([True, False], np.bool_),
        ([10.0, 20.0], np.float32),
        # Pandas nullable dtypes
        ([10, 20], pd.Int16Dtype()),
        ([True, False], pd.BooleanDtype()),
        ([10.0, 20.0], pd.Float32Dtype()),
        # PyArrow nullable dtypes
        ([10, 20], "int16[pyarrow]"),
        ([True, False], "bool[pyarrow]"),
        ([10.0, 20.0], "float32[pyarrow]"),
        # Strings
        (["a", "b"], object),
        (["foo", "bar"], "str"),
    ],
)


@series_map_dtypes
@pytest.mark.parametrize("base_cls", [pd.Series, pd.Index])
def test_series_map_dtype_all_mapped(base_cls, values, dtype):
    """When all keys mapped, the result dtype is equal to the mapper.dtype
    (basic use case).
    """
    base = base_cls([1, 2])
    mapper = pd.Series(values, index=[2, 1], dtype=dtype)
    assert_series_map_dtype(base, mapper)


@series_map_dtypes
@pytest.mark.parametrize("base_cls", [pd.Series, pd.Index])
def test_series_map_dtype_unmapped(xfail, base_cls, values, dtype):
    """When some keys fail to map, the output dtype can change dynamically;
    in Dask, this can't be reflected in the meta dtype.

    Note: in this test, the mismatch happens in only one partition, which means that
    finalize() will need to deal with mismatched dtypes in the final concatenation
    step.
    """
    if base_cls is pd.Index and dtype is np.bool_ and not PANDAS_GE_300:
        xfail("Index.map dtype mismatch")

    base = base_cls([1, 2])
    mapper = pd.Series(values, index=[2, 3], dtype=dtype)
    assert_series_map_dtype(base, mapper)


@series_map_dtypes
@pytest.mark.parametrize("base_dtype", [np.float32, pd.Int16Dtype(), "int16[pyarrow]"])
@pytest.mark.parametrize("base_cls", [pd.Series, pd.Index])
def test_series_map_dtype_null_base(xfail, base_cls, base_dtype, values, dtype):
    """A NULL in the base triggers dynamic dtype changes
    when the mapper's dtype is not nullable.

    Note: Pandas' behaviour changed in version 3.1
    """
    if base_cls is pd.Index and base_dtype in ("Int16", "int16[pyarrow]"):
        xfail("Index with nullable base dtype: NA in divisions")
    elif (
        base_cls is pd.Index
        and base_dtype is np.float32
        and dtype in (object, "str")
        # Test passes with Pandas 2.1 but not 2.0 or 2.2+ for some reason
        and (PANDAS_GE_220 or not PANDAS_GE_210)
    ):
        xfail("Index.map dtype mismatch: float32 NaN base with object/str mapper")
    elif (
        base_dtype in ("Int16", "int16[pyarrow]") and dtype is object and PANDAS_GE_310
    ):
        xfail("https://github.com/pandas-dev/pandas/issues/65735")

    base = base_cls([1, None], dtype=base_dtype)
    mapper_idx = pd.Index([2, 3], dtype=base_dtype)
    mapper = pd.Series(values, index=mapper_idx, dtype=dtype)
    assert_series_map_dtype(base, mapper)
