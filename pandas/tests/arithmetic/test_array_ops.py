import operator

import numpy as np
import pytest

from pandas.core.dtypes.missing import isna

import pandas._testing as tm
from pandas.core.ops.array_ops import (
    comparison_op,
    na_logical_op,
)


def test_na_logical_op_2d():
    left = np.arange(8).reshape(4, 2)
    right = left.astype(object)
    right[0, 0] = np.nan

    # Check that we fall back to the vec_binop branch
    with pytest.raises(TypeError, match="unsupported operand type"):
        operator.or_(left, right)

    result = na_logical_op(left, right, operator.or_)
    expected = right
    tm.assert_numpy_array_equal(result, expected)


def test_object_comparison_2d():
    left = np.arange(9).reshape(3, 3).astype(object)
    right = left.T

    result = comparison_op(left, right, operator.eq)
    expected = np.eye(3).astype(bool)
    tm.assert_numpy_array_equal(result, expected)

    # Ensure that cython doesn't raise on non-writeable arg, which
    #  we can get from np.broadcast_to
    right.flags.writeable = False
    result = comparison_op(left, right, operator.ne)
    tm.assert_numpy_array_equal(result, ~expected)


@pytest.mark.parametrize("rvalues", [1, [1, 1, 1], np.nan, None])
@pytest.mark.parametrize(
    "op", [operator.eq, operator.ne, operator.lt, operator.le, operator.gt, operator.ge]
)
def test_comparison_for_subclasses(rvalues, op):
    # GH#63205 Ensure subclasses of ndarray are correctly handled in comparison_op
    # Define a custom ndarray subclass
    class TestArray(np.ndarray):
        def __new__(cls, input_array):
            return np.asarray(input_array).view(cls)

        def __array_finalize__(self, obj) -> None:
            self._is_test_array = True

    def expected_with_na_handling(lvalues, rvalues, op):
        # Similar to comparison_op, handle zerodim arrays with na value separately
        if (rvalues.ndim == 0) and isna(rvalues.item()):
            # numpy does not like comparisons vs None
            if op is operator.ne:
                return np.ones(lvalues.shape, dtype=bool)
            else:
                return np.zeros(lvalues.shape, dtype=bool)
        return op(lvalues, rvalues)

    # Define test data
    lvalues = [1, 2, 3]

    # Test with both ndarray and TestArray
    result = comparison_op(np.array(lvalues), np.array(rvalues), op)
    expected = expected_with_na_handling(np.array(lvalues), np.array(rvalues), op)
    tm.assert_numpy_array_equal(result, expected)

    result = comparison_op(TestArray(lvalues), TestArray(rvalues), op)
    expected = expected_with_na_handling(TestArray(lvalues), TestArray(rvalues), op)
    tm.assert_numpy_array_equal(result, expected)
