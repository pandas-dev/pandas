"""
Tests for subclasses of NDArrayBackedExtensionArray
"""
import numpy as np
import pytest

from pandas.compat.numpy import IS_NEP18_ACTIVE

from pandas import (
    CategoricalIndex,
    date_range,
)
import pandas._testing as tm
from pandas.core.arrays import (
    Categorical,
    DatetimeArray,
    PandasArray,
    TimedeltaArray,
)


class TestEmpty:
    def test_empty_categorical(self):
        ci = CategoricalIndex(["a", "b", "c"], ordered=True)
        dtype = ci.dtype

        # case with int8 codes
        shape = (4,)
        result = Categorical._empty(shape, dtype=dtype)
        assert isinstance(result, Categorical)
        assert result.shape == shape
        assert result._ndarray.dtype == np.int8

        # case where repr would segfault if we didn't override base implementation
        result = Categorical._empty((4096,), dtype=dtype)
        assert isinstance(result, Categorical)
        assert result.shape == (4096,)
        assert result._ndarray.dtype == np.int8
        repr(result)

        # case with int16 codes
        ci = CategoricalIndex(list(range(512)) * 4, ordered=False)
        dtype = ci.dtype
        result = Categorical._empty(shape, dtype=dtype)
        assert isinstance(result, Categorical)
        assert result.shape == shape
        assert result._ndarray.dtype == np.int16

    def test_empty_dt64tz(self):
        dti = date_range("2016-01-01", periods=2, tz="Asia/Tokyo")
        dtype = dti.dtype

        shape = (0,)
        result = DatetimeArray._empty(shape, dtype=dtype)
        assert result.dtype == dtype
        assert isinstance(result, DatetimeArray)
        assert result.shape == shape

    def test_empty_dt64(self):
        shape = (3, 9)
        result = DatetimeArray._empty(shape, dtype="datetime64[ns]")
        assert isinstance(result, DatetimeArray)
        assert result.shape == shape

    def test_empty_td64(self):
        shape = (3, 9)
        result = TimedeltaArray._empty(shape, dtype="m8[ns]")
        assert isinstance(result, TimedeltaArray)
        assert result.shape == shape

    def test_empty_pandas_array(self):
        arr = PandasArray(np.array([1, 2]))
        dtype = arr.dtype

        shape = (3, 9)
        result = PandasArray._empty(shape, dtype=dtype)
        assert isinstance(result, PandasArray)
        assert result.dtype == dtype
        assert result.shape == shape


pytestmark = pytest.mark.skipif(
    not IS_NEP18_ACTIVE,
    reason="__array_function__ is not enabled by default until numpy 1.17",
)


class ArrayFunctionTests:
    # Tests for subclasses that do not explicitly support 2D yet.
    def test_delete_no_axis(self, array):
        # with no axis, operates on flattened version
        result = np.delete(array, 1)

        backing = np.delete(array._ndarray.ravel(), 1)
        expected = array._from_backing_data(backing)
        tm.assert_equal(result, expected)

    def test_repeat(self, array):
        result = np.repeat(array, 2)

        backing = np.repeat(array._ndarray.ravel(), 2)
        expected = array._from_backing_data(backing)
        tm.assert_equal(result, expected)

    def test_may_share_memory(self, array):
        assert np.may_share_memory(array, array)
        assert np.may_share_memory(array, array[:1])

        # partially-sharing memory
        assert np.may_share_memory(array[:-1], array[1:])

        assert not np.may_share_memory(array[:1], array[-1:])

        assert np.may_share_memory(array, array.ravel())
        assert np.may_share_memory(array, array.T)

        assert not np.may_share_memory(array, np.arange(0))

        assert np.may_share_memory(array, array._ndarray)
        # assert np.may_share_memory(array, Series(array.ravel("K")))  # nope!

    def test_ndim(self, array):
        assert np.ndim(array) == array.ndim

    def test_shape(self, array):
        assert np.shape(array) == array.shape


class ArrayFunctionTests2D(ArrayFunctionTests):
    @pytest.mark.parametrize("axis", [0, 1])
    def test_delete_axis(self, array, axis):
        result = np.delete(array, 1, axis=axis)
        if axis == 0:
            assert result.shape == (array.shape[0] - 1, array.shape[1])
        else:
            assert result.shape == (array.shape[0], array.shape[1] - 1)

        backing = np.delete(array._ndarray, 1, axis=axis)
        expected = array._from_backing_data(backing)
        tm.assert_equal(result, expected)

        # axis as an arg instead of as a kwarg
        result = np.delete(array, 1, axis)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("axis", [0, 1])
    def test_repeat_axis(self, array, axis):
        result = np.repeat(array, 2, axis=axis)

        backing = np.repeat(array._ndarray, 2, axis=axis)
        expected = array._from_backing_data(backing)
        tm.assert_equal(result, expected)

        # axis as an arg instead of a kwarg
        result = np.repeat(array, 2, axis)
        tm.assert_equal(result, expected)

    def test_atleast_2d(self, array):
        result = np.atleast_2d(array)

        assert result.ndim >= 2

        if array.ndim == 1:
            assert result.shape == (1, array.size)
        else:
            assert result.shape == array.shape

    def test_broadcast_to(self, array):
        array = array.reshape(-1, 1)
        shape = (array.shape[0], 3)

        result = np.broadcast_to(array, shape)

        backing = np.broadcast_to(array._ndarray, shape)
        expected = array._from_backing_data(backing)
        tm.assert_equal(result, expected)

        # TODO check that a copy was not made?

    def test_transpose(self, array):
        result = np.transpose(array, (0, 1))
        expected = array.transpose((0, 1))
        tm.assert_equal(result, expected)

        # no shape passed
        result = np.transpose(array)
        backing = np.transpose(array._ndarray)
        expected = array._from_backing_data(backing)
        tm.assert_equal(result, expected)

    def test_swapaxes(self, array):
        result = np.swapaxes(array, 1, 0)
        expected = array.swapaxes(1, 0)
        tm.assert_equal(result, expected)


class TestDatetimeArray(ArrayFunctionTests2D):
    @pytest.fixture(params=[2])
    def array(self, request):
        ndim = request.param
        dti = date_range("1994-05-12", periods=12, tz="US/Pacific")
        dta = dti._data
        if ndim == 2:
            dta = dta.reshape(3, 4)
        return dta


class TestTimedeltaArray(ArrayFunctionTests2D):
    @pytest.fixture(params=[2])
    def array(self, request):
        ndim = request.param
        dti = date_range("1994-05-12", periods=12, tz="US/Pacific")
        tdi = dti - dti[0]
        tda = tdi._data
        if ndim == 2:
            tda = tda.reshape(3, 4)
        return tda


class TestPeriodArray(ArrayFunctionTests2D):
    @pytest.fixture(params=[2])
    def array(self, request):
        ndim = request.param
        dti = date_range("1994-05-12", periods=12)
        pa = dti._data.to_period("D")
        if ndim == 2:
            pa = pa.reshape(3, 4)
        return pa


class TestPandasArray(ArrayFunctionTests):
    @pytest.fixture
    def array(self):
        return PandasArray(np.arange(12))


class TestCategorical(ArrayFunctionTests):
    @pytest.fixture
    def array(self):
        dti = date_range("1994-05-12", periods=12, tz="US/Pacific")
        return Categorical(dti)
