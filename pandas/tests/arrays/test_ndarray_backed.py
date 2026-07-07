"""
Tests for subclasses of NDArrayBackedExtensionArray
"""

import pickle

import numpy as np
import pytest

from pandas._libs.arrays import NDArrayBacked

from pandas import (
    CategoricalIndex,
    MultiIndex,
    date_range,
    timedelta_range,
)
import pandas._testing as tm
from pandas.core.arrays import (
    Categorical,
    DatetimeArray,
    NumpyExtensionArray,
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
        arr = NumpyExtensionArray(np.array([1, 2]))
        dtype = arr.dtype

        shape = (3, 9)
        result = NumpyExtensionArray._empty(shape, dtype=dtype)
        assert isinstance(result, NumpyExtensionArray)
        assert result.dtype == dtype
        assert result.shape == shape


@pytest.mark.parametrize(
    "arr",
    [
        Categorical(["a", "b", "a"]),
        date_range("2011-01-01", periods=3)._data,
        timedelta_range("1 Day", periods=3)._data,
    ],
)
def test_setstate_2tuple_without_attrs_dict(arr):
    # GH#63078, GH#62820: Cython 3.2's auto-pickle hands __setstate__ a bare
    #  (dtype, ndarray) 2-tuple with no trailing attrs dict, which previously
    #  raised NotImplementedError.
    result = type(arr).__new__(type(arr))
    NDArrayBacked.__setstate__(result, (arr.dtype, arr._ndarray))
    tm.assert_extension_array_equal(result, arr)


def test_pickle_datetime_multiindex_level():
    # GH#63078: pickling a datetime level of a MultiIndex raised
    #  NotImplementedError in NDArrayBacked.__setstate__
    mi = MultiIndex.from_product(
        [date_range("2011-01-01", periods=3), [1, 2]], names=["date", "id"]
    )
    lev = mi.get_level_values("date")
    tm.assert_index_equal(pickle.loads(pickle.dumps(lev)), lev)
