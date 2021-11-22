from __future__ import annotations

from typing import Callable

import numpy as np
import pytest

from pandas._typing import FloatFormatType

from pandas.core.dtypes import dtypes
from pandas.core.dtypes.common import is_extension_array_dtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import ExtensionArray


class DummyDtype(dtypes.ExtensionDtype):
    type = object
    name = "dummy"


class DummyArray(ExtensionArray):
    def __init__(self, data):
        self.data = data

    def __array__(self, dtype):
        return self.data

    @property
    def dtype(self):
        return DummyDtype()

    def astype(self, dtype, copy=True):
        # we don't support anything but a single dtype
        if isinstance(dtype, DummyDtype):
            if copy:
                return type(self)(self.data)
            return self

        return np.array(self, dtype=dtype, copy=copy)

    def __len__(self) -> int:
        return len(self.data)


class DummyArrayNoAsarray(DummyArray):
    def __array__(self, dtype=None):
        raise ValueError("Cannot be converted to an array!")

    def _format_array(
        self,
        formatter: Callable | None,
        *,
        float_format: FloatFormatType,
        na_rep: str = "NaN",
        digits: int,
        space: str | int,
        justify: str = "right",
        decimal: str = ".",
        leading_space: bool | None = True,
        quoting: int | None = None,
    ):
        return ["<MyEA>" for _ in self.data]


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


def test_repr_no_conversion():
    # https://github.com/pandas-dev/pandas/issues/26837#issuecomment-967268492
    # Validates
    s = pd.Series(DummyArrayNoAsarray([1]))
    repr(s)  # OK!

    df = pd.DataFrame({"A": DummyArrayNoAsarray([1])}, copy=False)
    repr(df)  # OK!
