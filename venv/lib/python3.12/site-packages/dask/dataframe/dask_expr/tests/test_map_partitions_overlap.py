from __future__ import annotations

import datetime

import numpy as np
import pytest

from dask.dataframe.dask_expr import from_pandas, map_overlap, map_partitions
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(100)})
    pdf["y"] = pdf.x // 7  # Not unique; duplicates span different partitions
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=10)


def test_map_partitions(df):
    def combine_x_y(x, y, foo=None):
        assert foo == "bar"
        return x + y

    df2 = df.map_partitions(combine_x_y, df + 1, foo="bar")
    assert_eq(df2, df + (df + 1))

    df2 = map_partitions(combine_x_y, df, df + 1, foo="bar")
    assert_eq(df2, df + (df + 1))


def test_map_partitions_broadcast(df):
    def combine_x_y(x, y, val, foo=None):
        assert foo == "bar"
        return x + y + val

    df2 = df.map_partitions(combine_x_y, df["x"].sum(), 123, foo="bar")
    assert_eq(df2, df + df["x"].sum() + 123)
    assert_eq(df2.optimize(), df + df["x"].sum() + 123)


@pytest.mark.parametrize("opt", [True, False])
def test_map_partitions_merge(opt):
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)})
    df1 = from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame({"x": range(0, 20, 2), "z": range(10)})
    df2 = from_pandas(pdf2, 1)

    # Partition-wise merge with map_partitions
    df3 = df1.map_partitions(
        lambda l, r: l.merge(r, on="x"),
        df2,
        enforce_metadata=False,
        clear_divisions=True,
    )

    # Check result with/without fusion
    expect = pdf1.merge(pdf2, on="x")
    df3 = (df3.optimize() if opt else df3)[list(expect.columns)]
    assert_eq(df3, expect, check_index=False)


def test_map_overlap():
    def func(x):
        x = x + x.sum()
        return x

    idx = pd.date_range("2020-01-01", periods=5, freq="D")
    pdf = pd.DataFrame(1, index=idx, columns=["a"])
    df = from_pandas(pdf, npartitions=2)

    result = df.map_overlap(func, before=0, after="2D")
    expected = pd.DataFrame([5, 5, 5, 3, 3], index=idx, columns=["a"])
    assert_eq(result, expected)
    result = df.map_overlap(func, before=0, after=1)
    assert_eq(result, expected)

    # Bug in dask/dask
    # result = df.map_overlap(func, before=0, after="1D")
    # expected = pd.DataFrame([4, 4, 4, 3, 3], index=idx, columns=["a"])
    # assert_eq(result, expected)

    result = df.map_overlap(func, before="2D", after=0)
    expected = pd.DataFrame(4, index=idx, columns=["a"])
    assert_eq(result, expected, check_index=False)

    result = df.map_overlap(func, before=1, after=0)
    assert_eq(result, expected, check_index=False)


def test_map_overlap_raises():
    def func(x):
        x = x + x.sum()
        return x

    idx = pd.date_range("2020-01-01", periods=5, freq="D")
    pdf = pd.DataFrame(1, index=idx, columns=["a"])
    df = from_pandas(pdf, npartitions=2)

    with pytest.raises(NotImplementedError, match="is less than"):
        df.map_overlap(func, before=5, after=0).compute()

    with pytest.raises(NotImplementedError, match="is less than"):
        df.map_overlap(func, before=0, after=5).compute()

    with pytest.raises(NotImplementedError, match="is less than"):
        df.map_overlap(func, before=0, after="5D").compute()

    with pytest.raises(ValueError, match="positive"):
        df.map_overlap(func, before=-1, after=5).compute()

    with pytest.raises(ValueError, match="positive"):
        df.map_overlap(func, before=1, after=-5).compute()


@pytest.mark.parametrize("npartitions", [1, 4])
def test_map_overlap(npartitions, pdf, df):
    def shifted_sum(df, before, after, c=0):
        a = df.shift(before)
        b = df.shift(-after)
        return df + a + b + c

    for before, after in [(0, 3), (3, 0), (3, 3), (0, 0)]:
        # DataFrame
        res = df.map_overlap(shifted_sum, before, after, before, after, c=2)
        sol = shifted_sum(pdf, before, after, c=2)
        assert_eq(res, sol)

        # Series
        res = df.x.map_overlap(shifted_sum, before, after, before, after, c=2)
        sol = shifted_sum(pdf.x, before, after, c=2)
        assert_eq(res, sol)


def test_map_overlap_divisions(df, pdf):
    result = df.shift(2)
    assert result.divisions == result.optimize().divisions
    result = df.ffill()
    assert result.divisions == result.optimize().divisions
    result = df.bfill()
    assert result.divisions == result.optimize().divisions

    pdf.index = pd.date_range("2019-12-31", freq="s", periods=len(pdf))
    df = from_pandas(pdf, npartitions=10)
    result = df.shift(freq="2s")
    assert result.known_divisions
    assert not result.optimize().known_divisions


def test_map_partitions_partition_info(df):
    partitions = {i: df.get_partition(i).compute() for i in range(df.npartitions)}

    def f(x, partition_info=None):
        assert partition_info is not None
        assert "number" in partition_info
        assert "division" in partition_info
        assert partitions[partition_info["number"]].equals(x)
        assert partition_info["division"] == x.index.min()
        return x

    df = df.map_partitions(f, meta=df._meta)
    result = df.compute(scheduler="single-threaded")
    assert type(result) == pd.DataFrame


def test_map_overlap_provide_meta():
    df = pd.DataFrame(
        {"x": [1, 2, 4, 7, 11], "y": [1.0, 2.0, 3.0, 4.0, 5.0]}
    ).rename_axis("myindex")
    ddf = from_pandas(df, npartitions=2)

    # Provide meta spec, but not full metadata
    res = ddf.map_overlap(
        lambda df: df.rolling(2).sum(), 2, 0, meta={"x": "i8", "y": "i8"}
    )
    sol = df.rolling(2).sum()
    assert_eq(res, sol)

    res = map_overlap(
        lambda df: df.rolling(2).sum(), ddf, 2, 0, meta={"x": "i8", "y": "i8"}
    )
    sol = df.rolling(2).sum()
    assert_eq(res, sol)


def test_map_overlap_errors(df):
    # Non-integer
    func = lambda x, *args, **kwargs: x
    with pytest.raises(ValueError):
        df.map_overlap(func, 0.5, 3, 0, 2, c=2)

    # Negative
    with pytest.raises(ValueError):
        df.map_overlap(func, 0, -5, 0, 2, c=2)

    # Partition size < window size
    with pytest.raises(NotImplementedError):
        df.map_overlap(func, 0, 100, 0, 100, c=2).compute()

    # Timedelta offset with non-datetime
    with pytest.raises(TypeError):
        df.map_overlap(func, pd.Timedelta("1s"), pd.Timedelta("1s"), 0, 2, c=2)

    # String timedelta offset with non-datetime
    with pytest.raises(TypeError):
        df.map_overlap(func, "1s", "1s", 0, 2, c=2)


@pytest.mark.parametrize("clear_divisions", [True, False])
def test_align_dataframes(clear_divisions):
    df1 = pd.DataFrame({"A": [1, 2, 3, 3, 2, 3], "B": [1, 2, 3, 4, 5, 6]})
    df2 = pd.DataFrame({"A": [3, 1, 2], "C": [1, 2, 3]})
    ddf1 = from_pandas(df1, npartitions=2)
    ddf2 = from_pandas(df2, npartitions=1)

    actual = ddf1.map_partitions(
        pd.merge, df2, align_dataframes=False, left_on="A", right_on="A", how="left"
    )
    expected = pd.merge(df1, df2, left_on="A", right_on="A", how="left")
    assert_eq(actual, expected, check_index=False, check_divisions=False)

    actual = ddf2.map_partitions(
        pd.merge,
        ddf1,
        align_dataframes=False,
        clear_divisions=clear_divisions,
        left_on="A",
        right_on="A",
        how="right",
    )
    expected = pd.merge(df2, df1, left_on="A", right_on="A", how="right")
    assert_eq(actual, expected, check_index=False, check_divisions=False)


def shifted_sum(df, before, after, c=0):
    a = df.shift(before)
    b = df.shift(-after)
    return df + a + b + c


@pytest.mark.parametrize("use_dask_input", [True, False])
@pytest.mark.parametrize("npartitions", [1, 4])
@pytest.mark.parametrize("enforce_metadata", [True, False])
@pytest.mark.parametrize("transform_divisions", [True, False])
@pytest.mark.parametrize("align_dataframes", [True, False])
@pytest.mark.parametrize(
    "overlap_setup",
    [
        (0, 3),
        (3, 0),
        (3, 3),
        (0, 0),
    ],
)
def test_map_overlap_multiple_dataframes(
    use_dask_input,
    npartitions,
    enforce_metadata,
    transform_divisions,
    align_dataframes,
    overlap_setup,
    pdf,
):
    before, after = overlap_setup

    ddf = pdf
    ddf2 = pdf * 2
    if use_dask_input:
        ddf = from_pandas(ddf, npartitions)
        ddf2 = from_pandas(ddf2, 2 if align_dataframes else npartitions)

    def get_shifted_sum_arg(overlap):
        return (
            overlap.seconds - 1 if isinstance(overlap, datetime.timedelta) else overlap
        )

    before_shifted_sum, after_shifted_sum = get_shifted_sum_arg(
        before
    ), get_shifted_sum_arg(after)

    # DataFrame
    res = map_overlap(
        shifted_sum,
        ddf,
        before,
        after,
        before_shifted_sum,
        after_shifted_sum,
        ddf2,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    sol = shifted_sum(pdf, before_shifted_sum, after_shifted_sum, pdf * 2)
    assert_eq(res, sol)

    # Series
    res = map_overlap(
        shifted_sum,
        ddf.x,
        before,
        after,
        before_shifted_sum,
        after_shifted_sum,
        ddf2.x,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    sol = shifted_sum(pdf.x, before_shifted_sum, after_shifted_sum, pdf.x * 2)
    assert_eq(res, sol)


def test_map_partitions_propagates_index_metadata():
    index = pd.Series(list("abcde"), name="myindex")
    df = pd.DataFrame(
        {"A": np.arange(5, dtype=np.int32), "B": np.arange(10, 15, dtype=np.int32)},
        index=index,
    )
    ddf = from_pandas(df, npartitions=2)
    res = ddf.map_partitions(
        lambda df: df.assign(C=df.A + df.B),
        meta=[("A", "i4"), ("B", "i4"), ("C", "i4")],
    )
    sol = df.assign(C=df.A + df.B)
    assert_eq(res, sol)

    res = ddf.map_partitions(lambda df: df.rename_axis("newindex"))
    sol = df.rename_axis("newindex")
    assert_eq(res, sol)


def test_token_given(df, pdf):
    def my_func(df):
        return df + 1

    result = df.map_partitions(my_func, token="token_given")
    assert result._name.split("-")[0] == "token_given"
    assert_eq(result, pdf + 1)

    result = df.map_partitions(my_func)
    assert result._name.split("-")[0] == "my_func"
    assert_eq(result, pdf + 1)
