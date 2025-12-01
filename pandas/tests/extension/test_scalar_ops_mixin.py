import numpy as np

import pandas as pd
from pandas.core.arrays.base import ExtensionScalarOpsMixin
import pandas.testing as tm


class CustomAddArray(pd.api.extensions.ExtensionArray, ExtensionScalarOpsMixin):
    def __init__(self, values):
        self._data = np.array(values)

    @property
    def dtype(self):
        return self._data.dtype

    def __getitem__(self, idx):
        return self._data[idx]

    def __len__(self):
        return len(self._data)

    def __array__(self):
        return self._data

    def _from_sequence(self, scalars, dtype=None):
        return type(self)(scalars)

    def _cast_pointwise_result(self, arr):
        return type(self)(arr)

    # Custom __add__ implementation
    def __add__(self, other):
        for i in range(len(self._data)):
            self._data[i] += other * 2
        return self


# Test fallback logic for arithmetic ops
def test_add_arithmetic_ops_custom():
    array_add = CustomAddArray([1, 2, 3])
    expected_add = CustomAddArray([7, 8, 9])

    # Add mixin ops
    CustomAddArray._add_arithmetic_ops()
    array_add += 3

    # Should use custom add (elementwise equality)
    tm.assert_numpy_array_equal(array_add._data, expected_add._data)

    # Check that another op (e.g., __sub__) is present and works
    assert hasattr(CustomAddArray, "__sub__")
    array_sub = CustomAddArray([1, 2, 3])
    expected_sub = CustomAddArray([0, 1, 2])

    array_sub -= 1

    assert isinstance(array_sub, CustomAddArray)
    tm.assert_numpy_array_equal(array_sub._data, expected_sub._data)


# Test fallback logic for comparison ops
class CustomEqArray(pd.api.extensions.ExtensionArray, ExtensionScalarOpsMixin):
    def __init__(self, values):
        self._data = np.array(values)

    @property
    def dtype(self):
        return self._data.dtype

    def __getitem__(self, idx):
        return self._data[idx]

    def __len__(self):
        return len(self._data)

    def __array__(self):
        return self._data

    def _from_sequence(self, scalars, dtype=None):
        return type(self)(scalars)

    def _cast_pointwise_result(self, arr):
        return type(self)(arr)

    # Custom __eq__ implementation
    def __eq__(self, other):
        # Dummy implementation
        for i in range(len(self._data)):
            if self._data[i] != other:
                return False
        return True


def test_add_comparison_ops_custom():
    arr_true = CustomEqArray([1, 1, 1])
    arr_false = CustomEqArray([1, 2, 3])
    CustomEqArray._add_comparison_ops()

    # Test custom __eq__ implementation
    result_true = arr_true == 1
    result_false = arr_false == 1
    assert result_true
    assert not result_false

    assert hasattr(CustomEqArray, "__ne__")
