import operator

import numpy as np
import pytest

from pandas import Series
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
def test_inplace_arithmetic_series(op, inplace_op):
    ser = Series([2, 4, 6], dtype="float64")
    data = get_array(ser)
    expected = op(ser, 2)
    inplace_op(ser, 2)
    assert np.shares_memory(get_array(ser), data)
    tm.assert_numpy_array_equal(data, get_array(ser))
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize("op,inplace_op", arith_operators)
def test_inplace_arithmetic_series_with_reference(using_copy_on_write, op, inplace_op):
    ser = Series([2, 4, 6], dtype="float64")
    ser_orig = ser.copy()
    view = ser[:]
    expected = op(ser, 2)
    inplace_op(ser, 2)
    if using_copy_on_write:
        assert not np.shares_memory(get_array(ser), get_array(view))
        tm.assert_series_equal(ser_orig, view)
    else:
        assert np.shares_memory(get_array(ser), get_array(view))
    tm.assert_series_equal(ser, expected)


# TODO:
#  Series <-> Series test
#  DF <-> DF
#  DF <-> Series
#  Alignment?
