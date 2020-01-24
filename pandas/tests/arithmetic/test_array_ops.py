import operator

import numpy as np
import pytest

import pandas._testing as tm
from pandas.core.ops.array_ops import na_logical_op


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
