import numpy as np
import pytest

from pandas.core.dtypes import dtypes
from pandas.core.dtypes.common import is_extension_array_dtype

import pandas as pd
from pandas.core.arrays import ExtensionArray
import pandas.util.testing as tm


class DummyDtype(dtypes.ExtensionDtype):
    pass


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


@pytest.mark.parametrize(
    "array",
    [
        pd.Series([1, None], dtype="Int64"),
        pd.Series(["2019", "2020"], dtype="datetime64[ns, UTC]"),
        pd.Series([0, 0], dtype="timedelta64[ns]"),
    ],
)
def test_compare_unequal_to_string(array):
    # GH 28930
    result = array == "a"
    expected = pd.Series([False, False])

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "array",
    [
        pd.Series([1, None], dtype="Int64"),
        pd.Series(["2019", "2020"], dtype="datetime64[ns, UTC]"),
        pd.Series([0, 0], dtype="timedelta64[ns]"),
    ],
)
@pytest.mark.parametrize("op", ["__lt__", "__le__", "__gt__", "__ge__"])
def test_compare_to_string_invalid(array, op):
    # GH 28930
    method = getattr(array, op)

    with pytest.raises(TypeError):
        method("a")
