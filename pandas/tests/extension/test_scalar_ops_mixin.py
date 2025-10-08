import numpy as np

import pandas as pd
from pandas.core.arrays.base import ExtensionScalarOpsMixin


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
        return "custom_add"


# Test fallback logic for arithmetic ops
def test_add_arithmetic_ops_custom():
    arr = CustomAddArray([1, 2, 3])
    # Remove __add__ if present, then add custom
    CustomAddArray.__add__ = lambda self, other: "custom_add"
    # Add mixin ops
    CustomAddArray._add_arithmetic_ops()
    result = arr + 1

    # Should use custom
    assert result == "custom_add"

    # Check that another op (e.g., __sub__) is present and works
    assert hasattr(CustomAddArray, "__sub__")
    sub_result = arr - 1
    assert isinstance(sub_result, CustomAddArray)


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


def test_add_comparison_ops_custom():
    arr = CustomEqArray([1, 2, 3])
    CustomEqArray.__eq__ = lambda self, other: self != other
    CustomEqArray._add_comparison_ops()
    result = arr == 1

    assert not result
    # Check that another op (e.g., __ne__) is present and works
    assert hasattr(CustomEqArray, "__ne__")
    ne_result = arr != 1
    assert isinstance(ne_result, np.ndarray) or isinstance(ne_result, CustomEqArray)
