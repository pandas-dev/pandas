import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas.core.dtypes import dtypes
from pandas.core.dtypes.common import is_extension_array_dtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import ExtensionArray


class DummyDtype(dtypes.ExtensionDtype):
    pass


class DummyArray(ExtensionArray):
    def __init__(self, data) -> None:
        self.data = data

    def __array__(self, dtype=None, copy=None):
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
        elif not copy:
            return np.asarray(self, dtype=dtype)
        else:
            return np.array(self, dtype=dtype, copy=copy)


class NoDunderDtype(dtypes.ExtensionDtype):
    name = "nodunder"
    type = object


class NoDunderArray(ExtensionArray):
    # relies on the default ExtensionArray.__array__
    def __init__(self, data) -> None:
        self.data = data

    def __getitem__(self, item):
        return self.data[item]

    def __len__(self) -> int:
        return len(self.data)

    def copy(self):
        return type(self)(self.data.copy())

    @property
    def dtype(self):
        return NoDunderDtype()


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


class CapturingStringArray(pd.arrays.StringArray):
    """Extend StringArray to capture arguments to __getitem__"""

    def __getitem__(self, item):
        self.last_item_arg = item
        return super().__getitem__(item)


def test_ellipsis_index():
    # GH#42430 1D slices over extension types turn into N-dimensional slices
    #  over ExtensionArrays
    dtype = pd.StringDtype()
    df = pd.DataFrame(
        {
            "col1": CapturingStringArray(
                np.array(["hello", "world"], dtype=object), dtype=dtype
            )
        }
    )
    _ = df.iloc[:1]

    # String comparison because there's no native way to compare slices.
    # Before the fix for GH#42430, last_item_arg would get set to the 2D slice
    # (Ellipsis, slice(None, 1, None))
    out = df["col1"]._values.last_item_arg
    assert str(out) == "slice(None, 1, None)"


@pytest.mark.parametrize(
    "func, expected",
    [
        (np.asarray, np.array([1, 2, 3])),
        (np.array, np.array([1, 2, 3])),
        (lambda arr: np.asarray(arr, dtype="float64"), np.array([1.0, 2.0, 3.0])),
        (lambda arr: arr.to_numpy(), np.array([1, 2, 3])),
        (lambda arr: pd.Series(arr).to_numpy(), np.array([1, 2, 3])),
    ],
)
def test_default_array_deprecated(func, expected):
    # GH#45093
    arr = NoDunderArray(np.array([1, 2, 3]))
    msg = "NoDunderArray does not implement __array__"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = func(arr)
    tm.assert_numpy_array_equal(result, expected)


def test_default_array_copy_false_raises():
    # GH#45093
    arr = NoDunderArray(np.array([1, 2, 3]))
    msg = "NoDunderArray does not implement __array__"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        with pytest.raises(ValueError, match="Unable to avoid copy"):
            np.asarray(arr, copy=False)


@pytest.mark.parametrize(
    "elements", [[[1, 2], [3, 4]], [[1], [2, 3]], [(1, 2), (3, 4)]]
)
def test_default_array_list_like_elements(elements):
    # GH#45093 list-like elements stay elements of a 1-D array
    data = np.empty(len(elements), dtype=object)
    data[:] = elements
    arr = NoDunderArray(data)
    msg = "NoDunderArray does not implement __array__"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = np.asarray(arr)
    tm.assert_numpy_array_equal(result, data)

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = pd.Index(arr, dtype=object)
    assert result.tolist() == elements

    ser = pd.Series([0, "x"], dtype=object)
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        ser[:] = arr
    assert ser.tolist() == elements


def test_default_array_list_like_elements_numeric_dtype():
    # GH#45093 a numeric dtype keeps numpy's 2-D result
    data = np.empty(2, dtype=object)
    data[:] = [(1, 2), (3, 4)]
    arr = NoDunderArray(data)
    msg = "NoDunderArray does not implement __array__"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = np.asarray(arr, dtype="float64")
    tm.assert_numpy_array_equal(result, np.array([[1.0, 2.0], [3.0, 4.0]]))


def test_default_array_defers_to_later_base():
    # GH#45093
    class ArrayMixin:
        def __array__(self, dtype=None, copy=None):
            return np.array([f"copy={copy}"] * len(self))

    class WithMixin(NoDunderArray, ArrayMixin):
        pass

    arr = WithMixin(np.array([1, 2]))
    with tm.assert_produces_warning(None):
        result = np.asarray(arr)
        copied = np.array(arr, copy=True)
    tm.assert_numpy_array_equal(result, np.array(["copy=None"] * 2))
    tm.assert_numpy_array_equal(copied, np.array(["copy=True"] * 2))


def test_default_array_2d():
    # GH#45093
    class TwoDimArray(NoDunderArray):
        @property
        def shape(self):
            return self.data.shape

        @property
        def ndim(self):
            return self.data.ndim

    data = np.arange(6.0).reshape(2, 3)
    arr = TwoDimArray(data)
    msg = "TwoDimArray does not implement __array__"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = np.asarray(arr)
    tm.assert_numpy_array_equal(result, data)
