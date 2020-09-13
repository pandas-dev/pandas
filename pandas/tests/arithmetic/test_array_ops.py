import operator

import numpy as np
import pytest

from pandas import Index, Series, array
import pandas._testing as tm
from pandas.core.ops.array_ops import comparison_op, na_logical_op


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


@pytest.fixture(params=[Index, Series, tm.to_array])
def box_pandas_1d_array(request):
    """
    Fixture to test behavior for Index, Series and tm.to_array classes
    """
    return request.param


@pytest.mark.parametrize(
    "data, expected", [([0, 1, 2], [0, 2, 4])],
)
def test_extension_array_add(box_pandas_1d_array, box_1d_array, data, expected):
    # GH22606 Verify operators with Extension Array and list-likes

    arr = array(data, dtype="Int64")
    container = box_pandas_1d_array(arr)
    left = container + box_1d_array(data)
    right = box_1d_array(data) + container

    assert left.tolist() == expected
    assert right.tolist() == expected
