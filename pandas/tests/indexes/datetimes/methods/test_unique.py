from datetime import (
    datetime,
    timedelta,
)

import pandas as pd
import pandas._testing as tm


def test_unique(tz_naive_fixture):
    idx = pd.DatetimeIndex(["2017"] * 2, tz=tz_naive_fixture)
    expected = idx[:1]

    result = idx.unique()
    tm.assert_index_equal(result, expected)
    # GH#21737
    # Ensure the underlying data is consistent
    assert result[0] == expected[0]


def test_index_unique(rand_series_with_duplicate_datetimeindex):
    dups = rand_series_with_duplicate_datetimeindex
    index = dups.index

    uniques = index.unique()
    expected = pd.DatetimeIndex(
        [
            datetime(2000, 1, 2),
            datetime(2000, 1, 3),
            datetime(2000, 1, 4),
            datetime(2000, 1, 5),
        ],
        dtype=index.dtype,
    )
    assert uniques.dtype == index.dtype  # sanity
    tm.assert_index_equal(uniques, expected)
    assert index.nunique() == 4

    # GH#2563
    assert isinstance(uniques, pd.DatetimeIndex)

    dups_local = index.tz_localize("US/Eastern")
    dups_local.name = "foo"
    result = dups_local.unique()
    expected = pd.DatetimeIndex(expected, name="foo")
    expected = expected.tz_localize("US/Eastern")
    assert result.tz is not None
    assert result.name == "foo"
    tm.assert_index_equal(result, expected)


def test_index_unique2():
    # NaT, note this is excluded
    arr = [1370745748 + t for t in range(20)] + [pd.NaT._value]
    idx = pd.DatetimeIndex(arr * 3)
    tm.assert_index_equal(idx.unique(), pd.DatetimeIndex(arr))
    assert idx.nunique() == 20
    assert idx.nunique(dropna=False) == 21


def test_index_unique3():
    arr = [
        pd.Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t) for t in range(20)
    ] + [pd.NaT]
    idx = pd.DatetimeIndex(arr * 3)
    tm.assert_index_equal(idx.unique(), pd.DatetimeIndex(arr))
    assert idx.nunique() == 20
    assert idx.nunique(dropna=False) == 21


def test_is_unique_monotonic(rand_series_with_duplicate_datetimeindex):
    index = rand_series_with_duplicate_datetimeindex.index
    assert not index.is_unique
