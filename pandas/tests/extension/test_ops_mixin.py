import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.base import ExtensionScalarOpsMixin


class _BaseArray(pd.api.extensions.ExtensionArray, ExtensionScalarOpsMixin):
    def __init__(self, data) -> None:
        self._data = np.asarray(data)

    @property
    def dtype(self):
        return self._data.dtype

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __len__(self) -> int:
        return len(self._data)

    def isna(self):
        return np.isnan(self._data.astype(float))

    def take(self, indices, allow_fill=False, fill_value=None):
        if allow_fill:
            result = np.full(len(indices), fill_value, dtype=self._data.dtype)
            mask = np.asarray(indices) >= 0
            result[mask] = self._data[np.asarray(indices)[mask]]
        else:
            result = self._data[np.asarray(indices)]
        return type(self)(result)

    def copy(self):
        return type(self)(self._data.copy())

    @classmethod
    def _from_sequence(cls, scalars, *, dtype=None, copy=False):
        return cls(scalars)

    def _cast_pointwise_result(self, arr):
        return type(self)(arr)


class CustomAddArray(_BaseArray):
    def __add__(self, other):
        return type(self)(self._data + other * 2)


CustomAddArray._add_arithmetic_ops()
CustomAddArray._add_comparison_ops()


class TestCustomArithmeticOp:
    def test_custom_add_is_preserved(self):
        arr = CustomAddArray([1, 2, 3])
        result = arr + 3
        tm.assert_numpy_array_equal(result._data, np.array([7, 8, 9]))

    def test_mixin_fallback_sub_is_registered(self):
        assert hasattr(CustomAddArray, "__sub__")

    def test_mixin_fallback_sub_works(self):
        arr = CustomAddArray([4, 5, 6])
        result = arr - 1
        assert isinstance(result, CustomAddArray)
        tm.assert_numpy_array_equal(result._data, np.array([3, 4, 5]))


class CustomEqArray(_BaseArray):
    def __eq__(self, other):
        return bool(np.all(self._data == other))


CustomEqArray._add_comparison_ops()


class TestCustomComparisonOp:
    def test_custom_eq_is_preserved(self):
        assert CustomEqArray([1, 1, 1]) == 1
        assert not (CustomEqArray([1, 2, 3]) == 1)

    def test_mixin_fallback_ne_is_registered(self):
        assert hasattr(CustomEqArray, "__ne__")

    def test_mixin_fallback_ne_works(self):
        arr = CustomEqArray([1, 1, 1])
        result = arr != 1
        assert isinstance(result, np.ndarray) or isinstance(result, CustomEqArray)

    def test_mixin_fallback_lt_is_registered(self):
        assert hasattr(CustomEqArray, "__lt__")


class PlainArray(_BaseArray):
    pass


PlainArray._add_arithmetic_ops()
PlainArray._add_comparison_ops()

_ALL_ARITHMETIC = [
    "__add__",
    "__radd__",
    "__sub__",
    "__rsub__",
    "__mul__",
    "__rmul__",
    "__pow__",
    "__rpow__",
    "__mod__",
    "__rmod__",
    "__floordiv__",
    "__rfloordiv__",
    "__truediv__",
    "__rtruediv__",
    "__divmod__",
    "__rdivmod__",
]
_ALL_COMPARISON = ["__eq__", "__ne__", "__lt__", "__gt__", "__le__", "__ge__"]


@pytest.mark.parametrize("op_name", _ALL_ARITHMETIC)
def test_all_arithmetic_ops_registered_when_none_defined(op_name):
    assert hasattr(PlainArray, op_name)


@pytest.mark.parametrize("op_name", _ALL_COMPARISON)
def test_all_comparison_ops_registered_when_none_defined(op_name):
    assert hasattr(PlainArray, op_name)
