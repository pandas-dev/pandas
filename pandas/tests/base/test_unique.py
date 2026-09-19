import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.base.common import allow_na_ops


def test_unique(index_or_series_obj):
    obj = index_or_series_obj
    obj = np.repeat(obj, range(1, len(obj) + 1))
    result = obj.unique()

    # dict.fromkeys preserves the order
    unique_values = list(dict.fromkeys(obj._values))
    if isinstance(obj, pd.MultiIndex):
        expected = pd.MultiIndex.from_tuples(unique_values)
        expected.names = obj.names
        tm.assert_index_equal(result, expected, exact=True, check_freq=False)
    elif isinstance(obj, pd.Index):
        expected = pd.Index(unique_values, dtype=obj.dtype)
        if isinstance(obj.dtype, pd.DatetimeTZDtype):
            expected = expected.normalize()
        tm.assert_index_equal(result, expected, exact=True, check_freq=False)
    else:
        expected = np.array(unique_values)
        tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("null_obj", [np.nan, None])
def test_unique_null(null_obj, index_or_series_obj, using_nan_is_na):
    obj = index_or_series_obj

    if not allow_na_ops(obj):
        pytest.skip("type doesn't allow for NA operations")
    elif len(obj) < 1:
        pytest.skip("Test doesn't make sense on empty data")
    elif isinstance(obj, pd.MultiIndex):
        pytest.skip(f"MultiIndex can't hold '{null_obj}'")
    elif (
        null_obj is not None
        and not using_nan_is_na
        and obj.dtype in ["Int64", "UInt16", "Float32"]
    ):
        pytest.skip("NaN is not a valid NA for this dtype.")

    obj = obj.copy(deep=True)
    values = obj._values
    values[0:2] = null_obj

    klass = type(obj)
    repeated_values = np.repeat(values, range(1, len(values) + 1))
    obj = klass(repeated_values, dtype=obj.dtype)
    result = obj.unique()

    unique_values_raw = dict.fromkeys(obj._values)
    # because np.nan == np.nan is False, but None == None is True
    # np.nan would be duplicated, whereas None wouldn't
    unique_values_not_null = [val for val in unique_values_raw if not pd.isnull(val)]
    unique_values = [null_obj, *unique_values_not_null]

    if isinstance(obj, pd.Index):
        expected = pd.Index(unique_values, dtype=obj.dtype)
        if isinstance(obj.dtype, pd.DatetimeTZDtype):
            result = result.normalize()
            expected = expected.normalize()
        tm.assert_index_equal(result, expected, exact=True)
    else:
        expected = np.array(unique_values, dtype=obj.dtype)
        tm.assert_numpy_array_equal(result, expected)


def test_nunique(index_or_series_obj):
    obj = index_or_series_obj
    obj = np.repeat(obj, range(1, len(obj) + 1))
    expected = len(obj.unique())
    assert obj.nunique(dropna=False) == expected


@pytest.mark.parametrize("null_obj", [np.nan, None])
def test_nunique_null(null_obj, index_or_series_obj):
    obj = index_or_series_obj

    if not allow_na_ops(obj):
        pytest.skip("type doesn't allow for NA operations")
    elif isinstance(obj, pd.MultiIndex):
        pytest.skip(f"MultiIndex can't hold '{null_obj}'")

    obj = obj.copy(deep=True)
    values = obj._values
    values[0:2] = null_obj

    klass = type(obj)
    repeated_values = np.repeat(values, range(1, len(values) + 1))
    obj = klass(repeated_values, dtype=obj.dtype)

    if isinstance(obj, pd.CategoricalIndex):
        assert obj.nunique() == len(obj.categories)
        assert obj.nunique(dropna=False) == len(obj.categories) + 1
    else:
        num_unique_values = len(obj.unique())
        assert obj.nunique() == max(0, num_unique_values - 1)
        assert obj.nunique(dropna=False) == max(0, num_unique_values)


@pytest.mark.single_cpu
def test_unique_bad_unicode(index_or_series):
    # regression test for #34550
    uval = "\ud83d"  # smiley emoji

    obj = index_or_series([uval] * 2, dtype=object)
    result = obj.unique()

    if isinstance(obj, pd.Index):
        expected = pd.Index(["\ud83d"], dtype=object)
        tm.assert_index_equal(result, expected, exact=True)
    else:
        expected = np.array(["\ud83d"], dtype=object)
        tm.assert_numpy_array_equal(result, expected)


@pytest.mark.single_cpu
def test_unique_distinct_bad_unicode(index_or_series):
    # GH#34550 distinct non-utf8-encodable values must not collapse together.
    #  Repeated values can pass even when the underlying buffers dangle.
    uvals = [chr(0xD800 + i) for i in range(50)]
    arr = np.empty(len(uvals), dtype=object)
    arr[:] = uvals

    obj = index_or_series(arr, dtype=object)
    result = obj.unique()

    assert len(result) == len(uvals)
    assert set(result) == set(uvals)

    codes, uniques = pd.factorize(arr)
    tm.assert_numpy_array_equal(codes, np.arange(len(uvals), dtype=np.intp))
    assert set(uniques) == set(uvals)


@pytest.mark.parametrize(
    "uvals",
    [
        ["", "\x00"],
        ["x\x00y", "x\x00z"],
        ["a", "a\x00", "a\x00\x00", "\x00a"],
    ],
)
def test_unique_embedded_null(index_or_series, uvals):
    # GH#34551 object-dtype strings route through StringHashTable, which used
    #  to key on a NUL-terminated C string and so collapsed distinct values
    #  sharing a prefix up to their first NUL.
    #  The values are repeated deliberately: Index.unique short-circuits on
    #  is_unique for an all-distinct input, so only a duplicated input reaches
    #  StringHashTable at all.
    vals = uvals * 2
    arr = np.empty(len(vals), dtype=object)
    arr[:] = vals

    result = index_or_series(arr, dtype=object).unique()
    assert list(result) == uvals

    codes, uniques = pd.factorize(arr)
    expected = np.tile(np.arange(len(uvals), dtype=np.intp), 2)
    tm.assert_numpy_array_equal(codes, expected)
    assert list(uniques) == uvals

    ser = pd.Series(arr, dtype=object)
    assert ser.nunique() == len(uvals)

    # groupby factorizes through the same table
    grouped = pd.DataFrame({"key": ser}).groupby("key")
    assert grouped.ngroups == len(uvals)


def test_nunique_dropna(dropna):
    # GH37566
    ser = pd.Series(["yes", "yes", pd.NA, np.nan, None, pd.NaT])
    res = ser.nunique(dropna)
    assert res == 1 if dropna else 5
