import datetime as dt
from datetime import datetime
from warnings import catch_warnings

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    MultiIndex,
    Series,
    concat,
)
import pandas._testing as tm
from pandas.core.arrays import SparseArray
from pandas.core.construction import create_series_with_explicit_dtype


@pytest.mark.parametrize("pdt", [pd.Series, pd.DataFrame])
@pytest.mark.parametrize("dt", np.sctypes["float"])
def test_concat_no_unnecessary_upcast(dt, pdt):
    # GH 13247
    dims = pdt(dtype=object).ndim

    dfs = [
        pdt(np.array([1], dtype=dt, ndmin=dims)),
        pdt(np.array([np.nan], dtype=dt, ndmin=dims)),
        pdt(np.array([5], dtype=dt, ndmin=dims)),
    ]
    x = pd.concat(dfs)
    assert x.values.dtype == dt


@pytest.mark.parametrize("pdt", [create_series_with_explicit_dtype, pd.DataFrame])
@pytest.mark.parametrize("dt", np.sctypes["int"])
def test_concat_will_upcast(dt, pdt):
    with catch_warnings(record=True):
        dims = pdt().ndim
        dfs = [
            pdt(np.array([1], dtype=dt, ndmin=dims)),
            pdt(np.array([np.nan], ndmin=dims)),
            pdt(np.array([5], dtype=dt, ndmin=dims)),
        ]
        x = pd.concat(dfs)
        assert x.values.dtype == "float64"


def test_concat_categorical_tz():
    # GH-23816
    a = pd.Series(pd.date_range("2017-01-01", periods=2, tz="US/Pacific"))
    b = pd.Series(["a", "b"], dtype="category")
    result = pd.concat([a, b], ignore_index=True)
    expected = pd.Series(
        [
            pd.Timestamp("2017-01-01", tz="US/Pacific"),
            pd.Timestamp("2017-01-02", tz="US/Pacific"),
            "a",
            "b",
        ]
    )
    tm.assert_series_equal(result, expected)


def test_concat_categorical_unchanged():
    # GH-12007
    # test fix for when concat on categorical and float
    # coerces dtype categorical -> float
    df = pd.DataFrame(pd.Series(["a", "b", "c"], dtype="category", name="A"))
    ser = pd.Series([0, 1, 2], index=[0, 1, 3], name="B")
    result = pd.concat([df, ser], axis=1)
    expected = pd.DataFrame(
        {
            "A": pd.Series(["a", "b", "c", np.nan], dtype="category"),
            "B": pd.Series([0, 1, np.nan, 2], dtype="float"),
        }
    )
    tm.assert_equal(result, expected)


def test_concat_datetimeindex_freq():
    # GH 3232
    # Monotonic index result
    dr = pd.date_range("01-Jan-2013", periods=100, freq="50L", tz="UTC")
    data = list(range(100))
    expected = pd.DataFrame(data, index=dr)
    result = pd.concat([expected[:50], expected[50:]])
    tm.assert_frame_equal(result, expected)

    # Non-monotonic index result
    result = pd.concat([expected[50:], expected[:50]])
    expected = pd.DataFrame(data[50:] + data[:50], index=dr[50:].append(dr[:50]))
    expected.index._data.freq = None
    tm.assert_frame_equal(result, expected)


def test_concat_sparse():
    # GH 23557
    a = pd.Series(SparseArray([0, 1, 2]))
    expected = pd.DataFrame(data=[[0, 0], [1, 1], [2, 2]]).astype(
        pd.SparseDtype(np.int64, 0)
    )
    result = pd.concat([a, a], axis=1)
    tm.assert_frame_equal(result, expected)


def test_concat_dense_sparse():
    # GH 30668
    a = pd.Series(pd.arrays.SparseArray([1, None]), dtype=float)
    b = pd.Series([1], dtype=float)
    expected = pd.Series(data=[1, None, 1], index=[0, 1, 0]).astype(
        pd.SparseDtype(np.float64, None)
    )
    result = pd.concat([a, b], axis=0)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("test_series", [True, False])
def test_concat_copy_index(test_series, axis):
    # GH 29879
    if test_series:
        ser = Series([1, 2])
        comb = concat([ser, ser], axis=axis, copy=True)
        assert comb.index is not ser.index
    else:
        df = DataFrame([[1, 2], [3, 4]], columns=["a", "b"])
        comb = concat([df, df], axis=axis, copy=True)
        assert comb.index is not df.index
        assert comb.columns is not df.columns


def test_concat_multiindex_datetime_object_index():
    # https://github.com/pandas-dev/pandas/issues/11058
    s = Series(
        ["a", "b"],
        index=MultiIndex.from_arrays(
            [[1, 2], Index([dt.date(2013, 1, 1), dt.date(2014, 1, 1)], dtype="object")],
            names=["first", "second"],
        ),
    )
    s2 = Series(
        ["a", "b"],
        index=MultiIndex.from_arrays(
            [[1, 2], Index([dt.date(2013, 1, 1), dt.date(2015, 1, 1)], dtype="object")],
            names=["first", "second"],
        ),
    )
    expected = DataFrame(
        [["a", "a"], ["b", np.nan], [np.nan, "b"]],
        index=MultiIndex.from_arrays(
            [
                [1, 2, 2],
                DatetimeIndex(
                    ["2013-01-01", "2014-01-01", "2015-01-01"],
                    dtype="datetime64[ns]",
                    freq=None,
                ),
            ],
            names=["first", "second"],
        ),
    )
    result = concat([s, s2], axis=1)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("keys", [["e", "f", "f"], ["f", "e", "f"]])
def test_duplicate_keys(keys):
    # GH 33654
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    s1 = Series([7, 8, 9], name="c")
    s2 = Series([10, 11, 12], name="d")
    result = concat([df, s1, s2], axis=1, keys=keys)
    expected_values = [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]
    expected_columns = pd.MultiIndex.from_tuples(
        [(keys[0], "a"), (keys[0], "b"), (keys[1], "c"), (keys[2], "d")]
    )
    expected = DataFrame(expected_values, columns=expected_columns)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "obj",
    [
        tm.SubclassedDataFrame({"A": np.arange(0, 10)}),
        tm.SubclassedSeries(np.arange(0, 10), name="A"),
    ],
)
def test_concat_preserves_subclass(obj):
    # GH28330 -- preserve subclass

    result = concat([obj, obj])
    assert isinstance(result, type(obj))


def test_concat_preserves_extension_int64_dtype():
    # GH 24768
    df_a = pd.DataFrame({"a": [-1]}, dtype="Int64")
    df_b = pd.DataFrame({"b": [1]}, dtype="Int64")
    result = pd.concat([df_a, df_b], ignore_index=True)
    expected = pd.DataFrame({"a": [-1, None], "b": [None, 1]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)
