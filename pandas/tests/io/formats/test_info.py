from io import StringIO
import re
import sys
import textwrap

import numpy as np
import pytest

from pandas.compat import PYPY

import pandas as pd


def test_info_categorical_column():

    # make sure it works
    n = 2500
    df = pd.DataFrame({"int64": np.random.randint(100, size=n)})
    df["category"] = pd.Series(
        np.array(list("abcdefghij")).take(np.random.randint(0, 10, size=n))
    ).astype("category")
    df.isna()
    buf = StringIO()
    df.info(buf=buf)

    df2 = df[df["category"] == "d"]
    buf = StringIO()
    df2.info(buf=buf)


def test_info(float_frame, datetime_frame):
    io = StringIO()
    float_frame.info(buf=io)
    datetime_frame.info(buf=io)

    frame = pd.DataFrame(np.random.randn(5, 3))

    frame.info()
    frame.info(verbose=False)


def test_info_verbose():
    buf = StringIO()
    size = 1001
    start = 5
    frame = pd.DataFrame(np.random.randn(3, size))
    frame.info(verbose=True, buf=buf)

    res = buf.getvalue()
    header = " #    Column  Dtype  \n---   ------  -----  "
    assert header in res

    frame.info(verbose=True, buf=buf)
    buf.seek(0)
    lines = buf.readlines()
    assert len(lines) > 0

    for i, line in enumerate(lines):
        if i >= start and i < start + size:
            line_nr = f" {i - start} "
            assert line.startswith(line_nr)


def test_info_memory():
    # https://github.com/pandas-dev/pandas/issues/21056
    df = pd.DataFrame({"a": pd.Series([1, 2], dtype="i8")})
    buf = StringIO()
    df.info(buf=buf)
    result = buf.getvalue()
    bytes = float(df.memory_usage().sum())
    expected = textwrap.dedent(
        f"""\
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 2 entries, 0 to 1
    Data columns (total 1 columns):
     #   Column  Non-Null Count  Dtype
    ---  ------  --------------  -----
     0   a       2 non-null      int64
    dtypes: int64(1)
    memory usage: {bytes} bytes
    """
    )
    assert result == expected


def test_info_wide():
    from pandas import set_option, reset_option

    io = StringIO()
    df = pd.DataFrame(np.random.randn(5, 101))
    df.info(buf=io)

    io = StringIO()
    df.info(buf=io, max_cols=101)
    rs = io.getvalue()
    assert len(rs.splitlines()) > 100
    xp = rs

    set_option("display.max_info_columns", 101)
    io = StringIO()
    df.info(buf=io)
    assert rs == xp
    reset_option("display.max_info_columns")


def test_info_duplicate_columns():
    io = StringIO()

    # it works!
    frame = pd.DataFrame(np.random.randn(1500, 4), columns=["a", "a", "b", "b"])
    frame.info(buf=io)


def test_info_duplicate_columns_shows_correct_dtypes():
    # GH11761
    io = StringIO()

    frame = pd.DataFrame([[1, 2.0]], columns=["a", "a"])
    frame.info(buf=io)
    io.seek(0)
    lines = io.readlines()
    assert " 0   a       1 non-null      int64  \n" == lines[5]
    assert " 1   a       1 non-null      float64\n" == lines[6]


def test_info_shows_column_dtypes():
    dtypes = [
        "int64",
        "float64",
        "datetime64[ns]",
        "timedelta64[ns]",
        "complex128",
        "object",
        "bool",
    ]
    data = {}
    n = 10
    for i, dtype in enumerate(dtypes):
        data[i] = np.random.randint(2, size=n).astype(dtype)
    df = pd.DataFrame(data)
    buf = StringIO()
    df.info(buf=buf)
    res = buf.getvalue()
    header = (
        " #   Column  Non-Null Count  Dtype          \n"
        "---  ------  --------------  -----          "
    )
    assert header in res
    for i, dtype in enumerate(dtypes):
        name = f" {i:d}   {i:d}       {n:d} non-null     {dtype}"
        assert name in res


def test_info_max_cols():
    df = pd.DataFrame(np.random.randn(10, 5))
    for len_, verbose in [(5, None), (5, False), (12, True)]:
        # For verbose always      ^ setting  ^ summarize ^ full output
        with pd.option_context("max_info_columns", 4):
            buf = StringIO()
            df.info(buf=buf, verbose=verbose)
            res = buf.getvalue()
            assert len(res.strip().split("\n")) == len_

    for len_, verbose in [(12, None), (5, False), (12, True)]:

        # max_cols not exceeded
        with pd.option_context("max_info_columns", 5):
            buf = StringIO()
            df.info(buf=buf, verbose=verbose)
            res = buf.getvalue()
            assert len(res.strip().split("\n")) == len_

    for len_, max_cols in [(12, 5), (5, 4)]:
        # setting truncates
        with pd.option_context("max_info_columns", 4):
            buf = StringIO()
            df.info(buf=buf, max_cols=max_cols)
            res = buf.getvalue()
            assert len(res.strip().split("\n")) == len_

        # setting wouldn't truncate
        with pd.option_context("max_info_columns", 5):
            buf = StringIO()
            df.info(buf=buf, max_cols=max_cols)
            res = buf.getvalue()
            assert len(res.strip().split("\n")) == len_


def test_info_memory_usage():
    # Ensure memory usage is displayed, when asserted, on the last line
    dtypes = [
        "int64",
        "float64",
        "datetime64[ns]",
        "timedelta64[ns]",
        "complex128",
        "object",
        "bool",
    ]
    data = {}
    n = 10
    for i, dtype in enumerate(dtypes):
        data[i] = np.random.randint(2, size=n).astype(dtype)
    df = pd.DataFrame(data)
    buf = StringIO()

    # display memory usage case
    df.info(buf=buf, memory_usage=True)
    res = buf.getvalue().splitlines()
    assert "memory usage: " in res[-1]

    # do not display memory usage case
    df.info(buf=buf, memory_usage=False)
    res = buf.getvalue().splitlines()
    assert "memory usage: " not in res[-1]

    df.info(buf=buf, memory_usage=True)
    res = buf.getvalue().splitlines()

    # memory usage is a lower bound, so print it as XYZ+ MB
    assert re.match(r"memory usage: [^+]+\+", res[-1])

    df.iloc[:, :5].info(buf=buf, memory_usage=True)
    res = buf.getvalue().splitlines()

    # excluded column with object dtype, so estimate is accurate
    assert not re.match(r"memory usage: [^+]+\+", res[-1])

    # Test a DataFrame with duplicate columns
    dtypes = ["int64", "int64", "int64", "float64"]
    data = {}
    n = 100
    for i, dtype in enumerate(dtypes):
        data[i] = np.random.randint(2, size=n).astype(dtype)
    df = pd.DataFrame(data)
    df.columns = dtypes

    df_with_object_index = pd.DataFrame({"a": [1]}, index=["foo"])
    df_with_object_index.info(buf=buf, memory_usage=True)
    res = buf.getvalue().splitlines()
    assert re.match(r"memory usage: [^+]+\+", res[-1])

    df_with_object_index.info(buf=buf, memory_usage="deep")
    res = buf.getvalue().splitlines()
    assert re.match(r"memory usage: [^+]+$", res[-1])

    # Ensure df size is as expected
    # (cols * rows * bytes) + index size
    df_size = df.memory_usage().sum()
    exp_size = len(dtypes) * n * 8 + df.index.nbytes
    assert df_size == exp_size

    # Ensure number of cols in memory_usage is the same as df
    size_df = np.size(df.columns.values) + 1  # index=True; default
    assert size_df == np.size(df.memory_usage())

    # assert deep works only on object
    assert df.memory_usage().sum() == df.memory_usage(deep=True).sum()

    # test for validity
    pd.DataFrame(1, index=["a"], columns=["A"]).memory_usage(index=True)
    pd.DataFrame(1, index=["a"], columns=["A"]).index.nbytes
    df = pd.DataFrame(
        data=1, index=pd.MultiIndex.from_product([["a"], range(1000)]), columns=["A"],
    )
    df.index.nbytes
    df.memory_usage(index=True)
    df.index.values.nbytes

    mem = df.memory_usage(deep=True).sum()
    assert mem > 0


@pytest.mark.skipif(PYPY, reason="on PyPy deep=True doesn't change result")
def test_info_memory_usage_deep_not_pypy():
    df_with_object_index = pd.DataFrame({"a": [1]}, index=["foo"])
    assert (
        df_with_object_index.memory_usage(index=True, deep=True).sum()
        > df_with_object_index.memory_usage(index=True).sum()
    )

    df_object = pd.DataFrame({"a": ["a"]})
    assert df_object.memory_usage(deep=True).sum() > df_object.memory_usage().sum()


@pytest.mark.skipif(not PYPY, reason="on PyPy deep=True does not change result")
def test_info_memory_usage_deep_pypy():
    df_with_object_index = pd.DataFrame({"a": [1]}, index=["foo"])
    assert (
        df_with_object_index.memory_usage(index=True, deep=True).sum()
        == df_with_object_index.memory_usage(index=True).sum()
    )

    df_object = pd.DataFrame({"a": ["a"]})
    assert df_object.memory_usage(deep=True).sum() == df_object.memory_usage().sum()


@pytest.mark.skipif(PYPY, reason="PyPy getsizeof() fails by design")
def test_usage_via_getsizeof():
    df = pd.DataFrame(
        data=1, index=pd.MultiIndex.from_product([["a"], range(1000)]), columns=["A"],
    )
    mem = df.memory_usage(deep=True).sum()
    # sys.getsizeof will call the .memory_usage with
    # deep=True, and add on some GC overhead
    diff = mem - sys.getsizeof(df)
    assert abs(diff) < 100


def test_info_memory_usage_qualified():

    buf = StringIO()
    df = pd.DataFrame(1, columns=list("ab"), index=[1, 2, 3])
    df.info(buf=buf)
    assert "+" not in buf.getvalue()

    buf = StringIO()
    df = pd.DataFrame(1, columns=list("ab"), index=list("ABC"))
    df.info(buf=buf)
    assert "+" in buf.getvalue()

    buf = StringIO()
    df = pd.DataFrame(
        1, columns=list("ab"), index=pd.MultiIndex.from_product([range(3), range(3)]),
    )
    df.info(buf=buf)
    assert "+" not in buf.getvalue()

    buf = StringIO()
    df = pd.DataFrame(
        1,
        columns=list("ab"),
        index=pd.MultiIndex.from_product([range(3), ["foo", "bar"]]),
    )
    df.info(buf=buf)
    assert "+" in buf.getvalue()


def test_info_memory_usage_bug_on_multiindex():
    # GH 14308
    # memory usage introspection should not materialize .values

    from string import ascii_uppercase as uppercase

    def memory_usage(f):
        return f.memory_usage(deep=True).sum()

    N = 100
    M = len(uppercase)
    index = pd.MultiIndex.from_product(
        [list(uppercase), pd.date_range("20160101", periods=N)], names=["id", "date"],
    )
    df = pd.DataFrame({"value": np.random.randn(N * M)}, index=index)

    unstacked = df.unstack("id")
    assert df.values.nbytes == unstacked.values.nbytes
    assert memory_usage(df) > memory_usage(unstacked)

    # high upper bound
    assert memory_usage(unstacked) - memory_usage(df) < 2000


def test_info_categorical():
    # GH14298
    idx = pd.CategoricalIndex(["a", "b"])
    df = pd.DataFrame(np.zeros((2, 2)), index=idx, columns=idx)

    buf = StringIO()
    df.info(buf=buf)
