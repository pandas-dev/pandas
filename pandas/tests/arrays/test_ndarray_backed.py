"""
Tests for subclasses of NDArrayBackedExtensionArray
"""
import numpy as np
import pytest

from pandas import (
    CategoricalIndex,
    date_range,
)
from pandas.core.arrays import (
    Categorical,
    DatetimeArray,
    PandasArray,
    PeriodArray,
    TimedeltaArray,
)


@pytest.fixture(
    params=[Categorical, DatetimeArray, TimedeltaArray, PeriodArray, PandasArray]
)
def ea_subclass(request):
    """
    Fixture for subclasses of NDArrayBackedExtensionArray.
    """
    return request.param


class TestEmpty:
    # def test_empty(self, ea_subclass):

    def test_empty_categorical(self):
        ci = CategoricalIndex(["a", "b", "c"], ordered=True)
        dtype = ci.dtype

        # case with int8 codes
        shape = (4,)
        result = Categorical.empty(shape, dtype=dtype)
        assert isinstance(result, Categorical)
        assert result.shape == shape
        assert result._ndarray.dtype == np.int8

        # case where repr would segfault if we didn't override base implementation
        result = Categorical.empty((4096,), dtype=dtype)
        assert isinstance(result, Categorical)
        assert result.shape == (4096,)
        assert result._ndarray.dtype == np.int8
        repr(result)

        # case with int16 codes
        ci = CategoricalIndex(list(range(512)) * 4, ordered=False)
        dtype = ci.dtype
        result = Categorical.empty(shape, dtype=dtype)
        assert isinstance(result, Categorical)
        assert result.shape == shape
        assert result._ndarray.dtype == np.int16

    def test_empty_dt64tz(self):
        dti = date_range("2016-01-01", periods=2, tz="Asia/Tokyo")
        dtype = dti.dtype

        shape = (0,)
        result = DatetimeArray.empty(shape, dtype=dtype)
        assert result.dtype == dtype
        assert isinstance(result, DatetimeArray)
        assert result.shape == shape

    def test_empty_dt64(self):
        shape = (3, 9)
        result = DatetimeArray.empty(shape, dtype="datetime64[ns]")
        assert isinstance(result, DatetimeArray)
        assert result.shape == shape

    def test_empty_td64(self):
        shape = (3, 9)
        result = TimedeltaArray.empty(shape, dtype="m8[ns]")
        assert isinstance(result, TimedeltaArray)
        assert result.shape == shape

    def test_empty_pandas_array(self):
        arr = PandasArray(np.array([1, 2]))
        dtype = arr.dtype

        shape = (3, 9)
        result = PandasArray.empty(shape, dtype=dtype)
        assert isinstance(result, PandasArray)
        assert result.dtype == dtype
        assert result.shape == shape
