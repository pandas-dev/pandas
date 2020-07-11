import numpy as np
import pytest

from pandas import Series


@pytest.mark.parametrize(
    "input_data", [[np.array([1, 2]), np.array([1, 2])]],
)
def test_equals_list_array(input_data):
    # GH20676 Verify equals operator for list of Numpy arrays
    ser = Series(input_data)
    assert ser.equals(ser)
