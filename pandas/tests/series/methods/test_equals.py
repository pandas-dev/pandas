import numpy as np

from pandas import Series


def test_equals_list_array():
    # GH20676 Verify equals operator for list of Numpy arrays
    arr = np.array([1, 2])
    ser = Series([arr, arr])
    assert ser.equals(ser)
