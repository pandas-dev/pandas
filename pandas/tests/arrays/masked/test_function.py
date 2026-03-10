import numpy as np
import pytest

from pandas.core.dtypes.common import is_integer_dtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import BaseMaskedArray

arrays = [pd.array([1, 2, 3, None], dtype=dtype) for dtype in tm.ALL_INT_EA_DTYPES]
arrays += [
    pd.array([0.141, -0.268, 5.895, None], dtype=dtype) for dtype in tm.FLOAT_EA_DTYPES
]


@pytest.fixture(params=arrays, ids=[a.dtype.name for a in arrays])
def data(request):
    """
    Fixture returning parametrized 'data' array with different integer and
    floating point types
    """
    return request.param


@pytest.fixture
def numpy_dtype(data):
    """
    Fixture returning numpy dtype from 'data' input array.
    """
    # For integer dtype, the numpy conversion must be done to float
    if is_integer_dtype(data):
        numpy_dtype = float
    else:
        numpy_dtype = data.dtype.type
    return numpy_dtype


def test_round(data, numpy_dtype):
    # No arguments
    result = data.round()
    np_result = np.round(data.to_numpy(dtype=numpy_dtype, na_value=None))
    exp_np = np_result.astype(object)
    exp_np[data.isna()] = pd.NA
    expected = pd.array(exp_np, dtype=data.dtype)
    tm.assert_extension_array_equal(result, expected)

    # Decimals argument
    result = data.round(decimals=2)
    np_result = np.round(data.to_numpy(dtype=numpy_dtype, na_value=None), decimals=2)
    exp_np = np_result.astype(object)
    exp_np[data.isna()] = pd.NA
    expected = pd.array(exp_np, dtype=data.dtype)
    tm.assert_extension_array_equal(result, expected)


def test_tolist(data):
    result = data.tolist()
    expected = list(data)
    tm.assert_equal(result, expected)


def test_to_numpy():
    # GH#56991

    class MyStringArray(BaseMaskedArray):
        dtype = pd.StringDtype()
        _dtype_cls = pd.StringDtype
        _internal_fill_value = pd.NA

    arr = MyStringArray(
        values=np.array(["a", "b", "c"]), mask=np.array([False, True, False])
    )
    result = arr.to_numpy()
    expected = np.array(["a", pd.NA, "c"])
    tm.assert_numpy_array_equal(result, expected)


class TestAnyAll2D:
    """Tests for any/all on 2D masked arrays with axis parameter."""

    @pytest.fixture
    def arr2d(self):
        # [[True,  False],
        #  [True,  <NA> ]]
        from pandas.core.arrays import BooleanArray

        data = np.array([[True, False], [True, True]], dtype=bool)
        mask = np.array([[False, False], [False, True]], dtype=bool)
        return BooleanArray._simple_new(data, mask)

    @pytest.mark.parametrize("method", ["any", "all"])
    def test_axis_none(self, arr2d, method):
        # axis=None reduces to scalar; should match flattened 1D result
        result = getattr(arr2d, method)(axis=None)
        expected = getattr(arr2d.ravel(), method)()
        assert result == expected

    def test_any_axis0(self, arr2d):
        result = arr2d.any(axis=0)
        expected = pd.array([True, False], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_any_axis1(self, arr2d):
        result = arr2d.any(axis=1)
        expected = pd.array([True, True], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_all_axis0(self, arr2d):
        result = arr2d.all(axis=0)
        expected = pd.array([True, False], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_all_axis1(self, arr2d):
        result = arr2d.all(axis=1)
        expected = pd.array([False, True], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_any_skipna_false_axis0(self, arr2d):
        # col0: True, True -> True; col1: False, NA -> NA (Kleene: False|NA=NA)
        result = arr2d.any(axis=0, skipna=False)
        expected = pd.array([True, pd.NA], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_any_skipna_false_axis1(self, arr2d):
        # row0: True, False -> True; row1: True, NA -> True (Kleene: True|NA=True)
        result = arr2d.any(axis=1, skipna=False)
        expected = pd.array([True, True], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_all_skipna_false_axis0(self, arr2d):
        # col0: True, True -> True (no NAs); col1: False, NA -> False (Kleene)
        result = arr2d.all(axis=0, skipna=False)
        expected = pd.array([True, False], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    def test_all_skipna_false_axis1(self, arr2d):
        # row0: True, False -> False; row1: True, NA -> NA (Kleene: True&NA=NA)
        result = arr2d.all(axis=1, skipna=False)
        expected = pd.array([False, pd.NA], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

    @pytest.mark.parametrize("method", ["any", "all"])
    def test_all_na_column(self, method):
        # Column that is all-NA should return falsey/truthy for skipna=True
        from pandas.core.arrays import IntegerArray

        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)

        result = getattr(arr, method)(axis=0)
        if method == "any":
            # col0: 1,3 -> True; col1: all-NA -> False
            expected = pd.array([True, False], dtype="boolean")
        else:
            # col0: 1,3 -> True; col1: all-NA -> True
            expected = pd.array([True, True], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)
