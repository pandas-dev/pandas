from collections import abc
import numbers

import numpy as np
import pytest

from pandas.core.dtypes import dtypes
from pandas.core.dtypes.common import is_extension_array_dtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import ExtensionArray


class DummyClass:
    pass


class DummyDtype(dtypes.ExtensionDtype):

    type = DummyClass
    name = "dummy"

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError(f"Cannot construct a '{cls}' from '{string}'")

    @classmethod
    def construct_array_type(cls):
        return DummyArray


class DummyArray(ExtensionArray):

    _dtype = DummyDtype

    def __init__(self, data):
        self.data = np.array(data)

    def __array__(self, dtype=None):
        return self.data

    @property
    def dtype(self):
        return DummyArray._dtype()

    def astype(self, dtype, copy=True):
        # we don't support anything but a single dtype
        if isinstance(dtype, self._dtype):
            if copy:
                return type(self)(self.data)
            return self

        return np.array(self, dtype=dtype, copy=copy)

    def __len__(self):
        return len(self.data)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        if isinstance(scalars, cls._dtype.type):
            scalars = [scalars]
        return DummyArray(scalars)

    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take

        data = self.astype(object)

        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value

        result = take(data, indices, fill_value=fill_value, allow_fill=allow_fill)
        return self._from_sequence(result, dtype=self.dtype)

    def isna(self):
        return np.array([x is self.dtype.na_value for x in self.data], dtype="bool")

    def __getitem__(self, idx):
        if isinstance(idx, numbers.Integral):
            return self.data[idx]
        elif isinstance(idx, (abc.Iterable, slice)):
            return DummyArray(self.data[idx])
        else:
            raise TypeError("Index type not supported", idx)


class TestExtensionArrayDtype:
    @pytest.mark.parametrize(
        "values",
        [
            pd.Categorical([]),
            pd.Categorical([]).dtype,
            pd.Series(pd.Categorical([])),
            DummyDtype(),
            DummyArray(np.array([1, 2])),
        ],
    )
    def test_is_extension_array_dtype(self, values):
        assert is_extension_array_dtype(values)

    @pytest.mark.parametrize("values", [np.array([]), pd.Series(np.array([]))])
    def test_is_not_extension_array_dtype(self, values):
        assert not is_extension_array_dtype(values)


def test_astype():

    arr = DummyArray(np.array([1, 2, 3]))
    expected = np.array([1, 2, 3], dtype=object)

    result = arr.astype(object)
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype("object")
    tm.assert_numpy_array_equal(result, expected)


def test_astype_no_copy():
    arr = DummyArray(np.array([1, 2, 3], dtype=np.int64))
    result = arr.astype(arr.dtype, copy=False)

    assert arr is result

    result = arr.astype(arr.dtype)
    assert arr is not result


@pytest.mark.parametrize("dtype", [dtypes.CategoricalDtype(), dtypes.IntervalDtype()])
def test_is_extension_array_dtype(dtype):
    assert isinstance(dtype, dtypes.ExtensionDtype)
    assert is_extension_array_dtype(dtype)


@pytest.mark.parametrize("na_value", [np.nan, pd.NA, None])
def test_empty_series_construction(monkeypatch, na_value):
    monkeypatch.setattr(DummyDtype, "na_value", na_value)
    result = pd.Series(index=[1, 2, 3], dtype=DummyDtype())
    expected = pd.Series([na_value] * 3, index=[1, 2, 3], dtype=DummyDtype())
    tm.assert_series_equal(result, expected)
