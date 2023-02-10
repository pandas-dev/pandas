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


@pytest.mark.parametrize("dtype", [None, "int64"])
def test_series_from_series(dtype, using_copy_on_write):
    # Case: constructing a Series from another Series object follows CoW rules:
    # a new object is returned and thus mutations are not propagated
    ser = Series([1, 2, 3], name="name")

    # default is copy=False -> new Series is a shallow copy / view of original
    result = Series(ser, dtype=dtype)

    # the shallow copy still shares memory
    assert np.shares_memory(ser.values, result.values)

    if using_copy_on_write:
        assert result._mgr.blocks[0].refs.has_reference()

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
    result = Series(ser, dtype=dtype)

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
        assert not result._mgr.blocks[0].refs.has_reference()


@pytest.mark.parametrize("func", [lambda x: x, lambda x: x._mgr])
@pytest.mark.parametrize("columns", [None, ["a"]])
def test_dataframe_constructor_mgr_or_df(using_copy_on_write, columns, func):
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()

    new_df = DataFrame(func(df))

    assert np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
    new_df.iloc[0] = 100

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
        tm.assert_frame_equal(df, df_orig)
    else:
        assert np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
        tm.assert_frame_equal(df, new_df)
