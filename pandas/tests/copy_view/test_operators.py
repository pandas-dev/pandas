import operator

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

arith_operators = [
    [operator.add, operator.iadd],
    [operator.sub, operator.isub],
    [operator.mul, operator.imul],
    [operator.truediv, operator.itruediv],
    [operator.floordiv, operator.ifloordiv],
    [operator.pow, operator.ipow],
]


@pytest.mark.parametrize("op,inplace_op", arith_operators)
def test_inplace_arith_series_scalar(op, inplace_op):
    ser = Series([2, 4, 6], dtype="float64")
    data = get_array(ser)
    expected = op(ser, 2)
    inplace_op(ser, 2)
    assert np.shares_memory(get_array(ser), data)
    tm.assert_numpy_array_equal(data, get_array(ser))
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize("op,inplace_op", arith_operators)
def test_inplace_arith_series_scalar_ref(using_copy_on_write, op, inplace_op):
    ser = Series([2, 4, 6], dtype="float64")
    ser1 = ser
    ser_orig = ser.copy()
    view = ser[:]
    expected = op(ser, 2)
    inplace_op(ser, 2)
    if using_copy_on_write:
        assert not np.shares_memory(get_array(ser), get_array(view))
        tm.assert_series_equal(ser_orig, view)
    else:
        assert np.shares_memory(get_array(ser), get_array(view))
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert ser is ser1
    assert ser._mgr is ser1._mgr
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize("op,inplace_op", arith_operators)
def test_inplace_arith_series_series(op, inplace_op):
    ser = Series([2, 4, 6], dtype="float64")
    other = Series([4, 6, 8], dtype="float64")
    data = get_array(ser)
    expected = op(ser, other)
    inplace_op(ser, other)
    assert np.shares_memory(get_array(ser), data)
    tm.assert_numpy_array_equal(data, get_array(ser))
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize("op,inplace_op", arith_operators)
def test_inplace_arith_series_series_ref(using_copy_on_write, op, inplace_op):
    ser = Series([2, 4, 6], dtype="float64")
    other = Series([4, 6, 8], dtype="float64")
    ser1 = ser
    ser_orig = ser.copy()
    view = ser[:]
    expected = op(ser, other)
    inplace_op(ser, other)
    if using_copy_on_write:
        assert not np.shares_memory(get_array(ser), get_array(view))
        tm.assert_series_equal(ser_orig, view)
    else:
        assert np.shares_memory(get_array(ser), get_array(view))
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert ser is ser1
    assert ser._mgr is ser1._mgr
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize(
    "op,inplace_op",
    # TODO: Add support for operating inplace to the rest of the arithmetic operators
    [[operator.add, operator.iadd]],
)
def test_inplace_arith_df_scalar(using_copy_on_write, op, inplace_op):
    df = DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    # Values will all be 1 block
    values = df._mgr.blocks[0].values
    df1 = df
    expected = op(df, 2)
    inplace_op(df, 2)
    for col in df.columns:
        assert np.shares_memory(get_array(df, col), values)
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert df is df1
    assert df._mgr is df1._mgr
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize(
    "op,inplace_op",
    # TODO: Add support for operating inplace to the rest of the arithmetic operators
    [[operator.add, operator.iadd]],
)
def test_inplace_arith_df_scalar_ref(using_copy_on_write, op, inplace_op):
    df = DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    df_orig = df.copy()
    # Values will all be 1 block
    values = df._mgr.blocks[0].values
    df1 = df
    view = df[:]
    expected = op(df, 2)
    inplace_op(df, 2)
    for col in df.columns:
        if using_copy_on_write:
            assert not np.shares_memory(get_array(df, col), values)
        else:
            assert np.shares_memory(get_array(df, col), values)
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert df is df1
    assert df._mgr is df1._mgr
    tm.assert_frame_equal(df, expected)
    if using_copy_on_write:
        tm.assert_frame_equal(view, df_orig)


@pytest.mark.parametrize(
    "op,inplace_op",
    # TODO: Add support for operating inplace to the rest of the arithmetic operators
    [[operator.add, operator.iadd]],
)
def test_inplace_arith_df_series(using_copy_on_write, op, inplace_op):
    df = DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    # Values will all be 1 block
    values = df._mgr.blocks[0].values
    df1 = df
    ser = Series([1, 2, 3], index=["a", "b", "c"], dtype="float64")
    expected = op(df, ser)
    inplace_op(df, ser)
    for col in df.columns:
        assert np.shares_memory(get_array(df, col), values)
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert df is df1
    assert df._mgr is df1._mgr
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize(
    "op,inplace_op",
    # TODO: Add support for operating inplace to the rest of the arithmetic operators
    [[operator.add, operator.iadd]],
)
def test_inplace_arith_df_series_ref(using_copy_on_write, op, inplace_op):
    df = DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    df_orig = df.copy()
    view = df[:]
    # Values will all be 1 block
    values = df._mgr.blocks[0].values
    df1 = df
    ser = Series([1, 2, 3], index=["a", "b", "c"], dtype="float64")
    expected = op(df, ser)
    inplace_op(df, ser)
    for col in df.columns:
        if using_copy_on_write:
            assert not np.shares_memory(get_array(df, col), values)
        else:
            assert np.shares_memory(get_array(df, col), values)
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert df is df1
    assert df._mgr is df1._mgr
    tm.assert_frame_equal(df, expected)
    if using_copy_on_write:
        tm.assert_frame_equal(view, df_orig)


@pytest.mark.parametrize(
    "op,inplace_op",
    # TODO: Add support for operating inplace to the rest of the arithmetic operators
    [[operator.add, operator.iadd]],
)
def test_inplace_arith_df_df(using_copy_on_write, op, inplace_op):
    df = DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    other = df.copy()
    # Values will all be 1 block
    values = df._mgr.blocks[0].values
    df1 = df
    expected = op(df, other)
    inplace_op(df, other)
    for col in df.columns:
        assert np.shares_memory(get_array(df, col), values)
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert df is df1
    assert df._mgr is df1._mgr
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize(
    "op,inplace_op",
    # TODO: Add support for operating inplace to the rest of the arithmetic operators
    [[operator.add, operator.iadd]],
)
def test_inplace_arith_df_df_ref(using_copy_on_write, op, inplace_op):
    df = DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    other = DataFrame(
        np.array([[2, 3, 4], [5, 6, 7], [8, 9, 10]], dtype="float64"),
        columns=["a", "b", "c"],
    )
    df_orig = df.copy()
    view = df[:]
    # Values will all be 1 block
    values = df._mgr.blocks[0].values
    df1 = df
    expected = op(df, other)
    inplace_op(df, other)
    for col in df.columns:
        if using_copy_on_write:
            assert not np.shares_memory(get_array(df, col), values)
        else:
            assert np.shares_memory(get_array(df, col), values)
    # Identity checks for the inplace op (check that an inplace op returns itself)
    assert df is df1
    assert df._mgr is df1._mgr
    tm.assert_frame_equal(df, expected)
    if using_copy_on_write:
        tm.assert_frame_equal(view, df_orig)
