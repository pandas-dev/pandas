import operator

import numpy as np
import pytest

import pandas._testing as tm

from pandas import DataFrame
from pandas.core.ops.array_ops import (
    comparison_op,
    na_logical_op,
)
from pandas.core.computation.expressions import (
    _evaluate_numexpr,
    _evaluate_standard,
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


def test_str_comparison():
    a = np.array(range(10 ** 5))
    right = DataFrame(a, dtype=np.int64)
    left = "    "
    
    result = right == left
    expected = DataFrame(np.zeros(right.shape, dtype=bool))
    tm.assert_frame_equal(result, expected)
