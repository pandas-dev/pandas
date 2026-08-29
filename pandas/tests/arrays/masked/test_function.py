import numpy as np
import pytest

from pandas.core.dtypes.common import is_integer_dtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import (
    BaseMaskedArray,
    BooleanArray,
    FloatingArray,
    IntegerArray,
)

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
        # GH#64510
        # [[True,  False],
        #  [True,  <NA> ]]
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


class TestReductions2D:
    """Tests for reductions on 2D masked arrays with axis parameter."""

    @pytest.fixture
    def int_arr2d(self):
        # [[1,  NA],
        #  [3,  4 ]]
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, False]], dtype=bool)
        return IntegerArray._simple_new(data, mask)

    @pytest.fixture
    def float_arr2d(self):
        # [[1.0,  NA ],
        #  [3.0,  4.0]]
        data = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        mask = np.array([[False, True], [False, False]], dtype=bool)
        return FloatingArray._simple_new(data, mask)

    @pytest.mark.parametrize("method", ["sum", "prod", "mean", "var", "std"])
    def test_axis_none(self, int_arr2d, method):
        # axis=None reduces to scalar; should match the flattened 1D result
        result = getattr(int_arr2d, method)(axis=None)
        expected = getattr(int_arr2d.ravel(), method)()
        assert result == expected or (pd.isna(result) and pd.isna(expected))

    @pytest.mark.parametrize("method", ["min", "max"])
    def test_minmax_axis_none(self, int_arr2d, method):
        result = getattr(int_arr2d, method)(axis=None)
        expected = getattr(int_arr2d.ravel(), method)()
        assert result == expected

    # --- sum ---

    def test_sum_axis0(self, int_arr2d):
        result = int_arr2d.sum(axis=0)
        # col0: 1+3=4, col1: 0+4=4 (NA skipped, identity=0, min_count=0)
        expected = pd.array([4, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_sum_axis1(self, int_arr2d):
        result = int_arr2d.sum(axis=1)
        # row0: 1+0=1 (NA skipped), row1: 3+4=7
        expected = pd.array([1, 7], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_sum_axis0_min_count(self, int_arr2d):
        # min_count=1: col1 has 1 valid value -> not NA
        result = int_arr2d.sum(axis=0, min_count=1)
        expected = pd.array([4, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_sum_axis0_min_count_all_na(self):
        # Column that is all-NA with min_count=1 -> NA
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.sum(axis=0, min_count=1)
        expected = pd.array([4, pd.NA], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    # --- prod ---

    def test_prod_axis0(self, int_arr2d):
        result = int_arr2d.prod(axis=0)
        # col0: 1*3=3, col1: 1*4=4 (NA skipped, identity=1, min_count=0)
        expected = pd.array([3, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_prod_axis1(self, int_arr2d):
        result = int_arr2d.prod(axis=1)
        # row0: 1*1=1 (NA skipped), row1: 3*4=12
        expected = pd.array([1, 12], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    # --- mean ---

    def test_mean_axis0(self, int_arr2d):
        result = int_arr2d.mean(axis=0)
        # col0: (1+3)/2=2.0, col1: 4/1=4.0
        expected = pd.array([2.0, 4.0], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

    def test_mean_axis1(self, int_arr2d):
        result = int_arr2d.mean(axis=1)
        # row0: 1/1=1.0, row1: (3+4)/2=3.5
        expected = pd.array([1.0, 3.5], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

    def test_mean_axis0_all_na(self):
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.mean(axis=0)
        expected = pd.array([2.0, pd.NA], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

    # --- var / std ---

    def test_var_axis0(self, float_arr2d):
        result = float_arr2d.var(axis=0, ddof=0)
        # col0: var([1,3])=1.0, col1: var([4])=0.0
        expected = pd.array([1.0, 0.0], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

    def test_std_axis0(self, float_arr2d):
        result = float_arr2d.std(axis=0, ddof=0)
        expected = pd.array([1.0, 0.0], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

    # --- min / max ---

    def test_min_axis0(self, int_arr2d):
        result = int_arr2d.min(axis=0)
        # col0: min(1,3)=1, col1: min(4)=4 (NA skipped)
        expected = pd.array([1, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_max_axis0(self, int_arr2d):
        result = int_arr2d.max(axis=0)
        # col0: max(1,3)=3, col1: max(4)=4
        expected = pd.array([3, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_min_axis1(self, int_arr2d):
        result = int_arr2d.min(axis=1)
        # row0: min(1)=1, row1: min(3,4)=3
        expected = pd.array([1, 3], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_max_axis1(self, int_arr2d):
        result = int_arr2d.max(axis=1)
        # row0: max(1)=1, row1: max(3,4)=4
        expected = pd.array([1, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_min_axis0_all_na(self):
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.min(axis=0)
        expected = pd.array([1, pd.NA], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    def test_max_axis0_all_na(self):
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.max(axis=0)
        expected = pd.array([3, pd.NA], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

    # --- skipna=False ---

    @pytest.mark.parametrize("method", ["sum", "prod", "mean", "var", "std"])
    def test_skipna_false_axis0(self, int_arr2d, method):
        # col0: no NAs -> valid result, col1: has NA -> NA
        kwargs = {}
        if method in ("var", "std"):
            kwargs["ddof"] = 0
        result = getattr(int_arr2d, method)(axis=0, skipna=False, **kwargs)
        assert isinstance(result, BaseMaskedArray)
        assert not pd.isna(result[0])
        assert pd.isna(result[1])

    @pytest.mark.parametrize("method", ["min", "max"])
    def test_minmax_skipna_false_axis0(self, int_arr2d, method):
        result = getattr(int_arr2d, method)(axis=0, skipna=False)
        assert isinstance(result, BaseMaskedArray)
        assert not pd.isna(result[0])
        assert pd.isna(result[1])
