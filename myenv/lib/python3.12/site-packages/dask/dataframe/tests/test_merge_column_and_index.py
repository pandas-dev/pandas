from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import dask.dataframe as dd
from dask.dataframe._compat import PANDAS_GE_210
from dask.dataframe.utils import assert_eq


# Fixtures
# ========
@pytest.fixture
def df_left():
    # Create frame with 10 partitions
    # Frame has 11 distinct idx values
    partition_sizes = np.array([3, 4, 2, 5, 3, 2, 5, 9, 4, 7, 4])
    idx = [i for i, s in enumerate(partition_sizes) for _ in range(s)]
    k = [i for s in partition_sizes for i in range(s)]
    vi = range(len(k))

    return pd.DataFrame(dict(idx=idx, k=k, v1=vi)).set_index(["idx"])


@pytest.fixture
def df_right():
    # Create frame with 10 partitions
    # Frame has 11 distinct idx values
    partition_sizes = np.array([4, 2, 5, 3, 2, 5, 9, 4, 7, 4, 8])
    idx = [i for i, s in enumerate(partition_sizes) for _ in range(s)]
    k = [i for s in partition_sizes for i in range(s)]
    vi = range(len(k))

    return pd.DataFrame(dict(idx=idx, k=k, v1=vi)).set_index(["idx"])


@pytest.fixture
def ddf_left(df_left):
    # Create frame with 10 partitions
    # Skip division on 2 so there is one mismatch with ddf_right
    return dd.repartition(df_left, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11])


@pytest.fixture
def ddf_left_unknown(ddf_left):
    return ddf_left.clear_divisions()


@pytest.fixture
def ddf_left_single(df_left):
    return dd.from_pandas(df_left, npartitions=1, sort=False)


@pytest.fixture
def ddf_right(df_right):
    # Create frame with 10 partitions
    # Skip division on 3 so there is one mismatch with ddf_left
    return dd.repartition(df_right, [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11])


@pytest.fixture
def ddf_right_unknown(ddf_right):
    return ddf_right.clear_divisions()


@pytest.fixture
def ddf_right_single(df_right):
    return dd.from_pandas(df_right, npartitions=1, sort=False)


@pytest.fixture
def ddf_right_double(df_right):
    return dd.from_pandas(df_right, npartitions=2, sort=False)


@pytest.fixture
def ddf_left_double(df_left):
    return dd.from_pandas(df_left, npartitions=2, sort=False)


@pytest.fixture(params=["inner", "left", "right", "outer"])
def how(request):
    return request.param


@pytest.fixture(params=["idx", ["idx"], ["idx", "k"], ["k", "idx"]])
def on(request):
    return request.param


# Tests
# =====
def test_merge_known_to_known(
    df_left, df_right, ddf_left, ddf_right, on, how, shuffle_method
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(ddf_right, on=on, how=how, shuffle_method=shuffle_method)

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(range(12)))
    assert len(result.__dask_graph__()) < 80


@pytest.mark.parametrize("how", ["inner", "left"])
def test_merge_known_to_single(
    df_left, df_right, ddf_left, ddf_right_single, on, how, shuffle_method
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(
        ddf_right_single, on=on, how=how, shuffle_method=shuffle_method
    )

    # Assertions
    assert_eq(result, expected)
    assert result.divisions == ddf_left.divisions
    assert len(result.__dask_graph__()) < 30


@pytest.mark.parametrize("how", ["inner", "right"])
def test_merge_single_to_known(
    df_left, df_right, ddf_left_single, ddf_right, on, how, shuffle_method
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left_single.merge(
        ddf_right, on=on, how=how, shuffle_method=shuffle_method
    )

    # Assertions
    assert_eq(result, expected)
    assert result.divisions == ddf_right.divisions
    assert len(result.__dask_graph__()) < 30


def test_merge_known_to_unknown(
    df_left, df_right, ddf_left, ddf_right_unknown, on, how, shuffle_method
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(
        ddf_right_unknown, on=on, how=how, shuffle_method=shuffle_method
    )

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(None for _ in range(11)))


def test_merge_unknown_to_known(
    df_left, df_right, ddf_left_unknown, ddf_right, on, how, shuffle_method
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left_unknown.merge(
        ddf_right, on=on, how=how, shuffle_method=shuffle_method
    )

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(None for _ in range(11)))


def test_merge_unknown_to_unknown(
    df_left,
    df_right,
    ddf_left_unknown,
    ddf_right_unknown,
    on,
    how,
    shuffle_method,
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Merge unknown to unknown
    result = ddf_left_unknown.merge(
        ddf_right_unknown, on=on, how=how, shuffle_method=shuffle_method
    )

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(None for _ in range(11)))


@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
@pytest.mark.parametrize("shuffle_method", ["tasks", "disk"])
def test_join(how, shuffle_method):
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)})
    df1 = dd.from_pandas(pdf1, 4)
    pdf2 = pd.DataFrame({"z": range(10)}, index=pd.Index(range(10), name="a"))
    df2 = dd.from_pandas(pdf2, 2)

    # Partition-wise merge with map_partitions
    df3 = df1.join(
        df2.clear_divisions(), on="x", how=how, shuffle_method=shuffle_method
    )

    # Check result with/without fusion
    expect = pdf1.join(pdf2, on="x", how=how)
    assert_eq(df3, expect, check_index=False, check_dtype=False)


@pytest.mark.parametrize("how", ["inner", "left"])
def test_merge_known_to_double_bcast_right(
    df_left, df_right, ddf_left, ddf_right_double, on, how, shuffle_method
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(
        ddf_right_double, on=on, how=how, shuffle_method=shuffle_method, broadcast=True
    )

    # Assertions
    assert_eq(result, expected)
    # Hash join used in disk-shuffling doesn't preserve divisions.
    if shuffle_method == "task":
        assert_eq(result.divisions, ddf_left.divisions)


@pytest.mark.parametrize("how", ["inner", "right"])
@pytest.mark.parametrize("broadcast", [True, 0.75])
def test_merge_known_to_double_bcast_left(
    df_left, df_right, ddf_left_double, ddf_right, on, shuffle_method, how, broadcast
):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left_double.merge(
        ddf_right, on=on, how=how, broadcast=broadcast, shuffle_method=shuffle_method
    )

    # Assertions
    assert_eq(result, expected)
    # Hash join used in disk-shuffling doesn't preserve divisions.
    if shuffle_method == "task":
        assert_eq(result.divisions, ddf_right.divisions)

    # Check that culling works
    result.head(1)


@pytest.mark.skipif(PANDAS_GE_210, reason="breaks with pandas=2.1.0+")
@pytest.mark.parametrize("repartition", [None, 4])
def test_merge_column_with_nulls(repartition):
    # See: https://github.com/dask/dask/issues/7558

    df1 = pd.DataFrame({"a": ["0", "0", None, None, None, None, "5", "7", "15", "33"]})
    df2 = pd.DataFrame({"c": ["1", "2", "3", "4"], "b": ["0", "5", "7", "15"]})
    df1_d = dd.from_pandas(df1, npartitions=4)
    df2_d = dd.from_pandas(df2, npartitions=3).set_index("b")
    if repartition:
        df2_d = df2_d.repartition(npartitions=repartition)

    pandas_result = df1.merge(
        df2.set_index("b"), how="left", left_on="a", right_index=True
    )
    dask_result = df1_d.merge(df2_d, how="left", left_on="a", right_index=True)

    assert_eq(dask_result, pandas_result)
