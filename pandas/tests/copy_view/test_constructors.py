import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

# -----------------------------------------------------------------------------
# Copy/view behaviour for Series / DataFrame constructors


def test_series_from_series(using_copy_on_write):
    # Case: constructing a Series from another Series object follows CoW rules:
    # a new object is returned and thus mutations are not propagated
    ser = Series([1, 2, 3], name="name")

    # default is copy=False -> new Series is a shallow copy / view of original
    result = Series(ser)

    # the shallow copy still shares memory
    assert np.shares_memory(ser.values, result.values)

    if using_copy_on_write:
        assert result._mgr.refs is not None

    if using_copy_on_write:
        # mutating new series copy doesn't mutate original
        result.iloc[0] = 0
        assert ser.iloc[0] == 1
        # mutating triggered a copy-on-write -> no longer shares memory
        assert not np.shares_memory(ser.values, result.values)
    else:
        # mutating shallow copy does mutate original
        result.iloc[0] = 0
        assert ser.iloc[0] == 0
        # and still shares memory
        assert np.shares_memory(ser.values, result.values)

    # the same when modifying the parent
    result = Series(ser)

    if using_copy_on_write:
        # mutating original doesn't mutate new series
        ser.iloc[0] = 0
        assert result.iloc[0] == 1
    else:
        # mutating original does mutate shallow copy
        ser.iloc[0] = 0
        assert result.iloc[0] == 0


def test_series_from_series_with_reindex(using_copy_on_write):
    # Case: constructing a Series from another Series with specifying an index
    # that potentially requires a reindex of the values
    ser = Series([1, 2, 3], name="name")

    # passing an index that doesn't actually require a reindex of the values
    # -> without CoW we get an actual mutating view
    for index in [
        ser.index,
        ser.index.copy(),
        list(ser.index),
        ser.index.rename("idx"),
    ]:
        result = Series(ser, index=index)
        assert np.shares_memory(ser.values, result.values)
        result.iloc[0] = 0
        if using_copy_on_write:
            assert ser.iloc[0] == 1
        else:
            assert ser.iloc[0] == 0

    # ensure that if an actual reindex is needed, we don't have any refs
    # (mutating the result wouldn't trigger CoW)
    result = Series(ser, index=[0, 1, 2, 3])
    assert not np.shares_memory(ser.values, result.values)
    if using_copy_on_write:
        assert result._mgr.refs is None or result._mgr.refs[0] is None


@pytest.mark.parametrize("dtype", [None, "int64"])
@pytest.mark.parametrize("columns", [None, ["a", "b"], ["a", "b", "c"]])
def test_dataframe_from_dict_of_series(using_copy_on_write, columns, dtype):
    # Case: constructing a DataFrame from Series objects with copy=False
    # has to do a lazy following CoW rules
    # (the default for DataFrame(dict) is still to copy to ensure consolidation)
    s1 = Series([1, 2, 3])
    s2 = Series([4, 5, 6])
    s1_orig = s1.copy()
    expected = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, dtype=dtype, columns=columns)

    result = DataFrame({"a": s1, "b": s2}, columns=columns, dtype=dtype, copy=False)

    # the shallow copy still shares memory
    assert np.shares_memory(get_array(result, "a"), s1.values)

    # mutating the new dataframe doesn't mutate original
    result.iloc[0, 0] = 10
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), s1.values)
        tm.assert_series_equal(s1, s1_orig)
    else:
        assert s1.iloc[0] == 10

    # the same when modifying the parent series
    result = DataFrame({"a": s1, "b": s2}, columns=columns, dtype=dtype, copy=False)
    s1.iloc[0] = 10
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), s1.values)
        tm.assert_frame_equal(result, expected)
    else:
        assert result.iloc[0, 0] == 10
