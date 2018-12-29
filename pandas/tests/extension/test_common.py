import numpy as np
import pytest

from pandas.compat import StringIO

from pandas.core.dtypes import dtypes
from pandas.core.dtypes.common import is_extension_array_dtype

import pandas as pd
from pandas.core.arrays import ExtensionArray
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype, UInt8Dtype, UInt16Dtype,
    UInt32Dtype, UInt64Dtype, integer_array)
import pandas.util.testing as tm


def make_data():
    return (list(range(1, 9)) + [np.nan] + list(range(10, 98))
            + [np.nan] + [99, 100])


@pytest.fixture(params=[Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype,
                        UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype])
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    return integer_array(make_data(), dtype=dtype)


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


class TestExtensionArrayDtype(object):

    @pytest.mark.parametrize('values', [
        pd.Categorical([]),
        pd.Categorical([]).dtype,
        pd.Series(pd.Categorical([])),
        DummyDtype(),
        DummyArray(np.array([1, 2])),
    ])
    def test_is_extension_array_dtype(self, values):
        assert is_extension_array_dtype(values)

    @pytest.mark.parametrize('values', [
        np.array([]),
        pd.Series(np.array([])),
    ])
    def test_is_not_extension_array_dtype(self, values):
        assert not is_extension_array_dtype(values)


def test_astype():

    arr = DummyArray(np.array([1, 2, 3]))
    expected = np.array([1, 2, 3], dtype=object)

    result = arr.astype(object)
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype('object')
    tm.assert_numpy_array_equal(result, expected)


def test_astype_no_copy():
    arr = DummyArray(np.array([1, 2, 3], dtype=np.int64))
    result = arr.astype(arr.dtype, copy=False)

    assert arr is result

    result = arr.astype(arr.dtype)
    assert arr is not result


@pytest.mark.parametrize('dtype', [
    dtypes.DatetimeTZDtype('ns', 'US/Central'),
])
def test_is_not_extension_array_dtype(dtype):
    assert not isinstance(dtype, dtypes.ExtensionDtype)
    assert not is_extension_array_dtype(dtype)


@pytest.mark.parametrize('dtype', [
    dtypes.CategoricalDtype(),
    dtypes.IntervalDtype(),
])
def test_is_extension_array_dtype(dtype):
    assert isinstance(dtype, dtypes.ExtensionDtype)
    assert is_extension_array_dtype(dtype)


@pytest.mark.parametrize('engine', ['c', 'python'])
def test_EA_types(engine):
    df = pd.DataFrame({'Int': pd.Series([1, 2, 3], dtype='Int64'),
                       'A': [1, 2, 1]})
    data = df.to_csv(index=False)
    result = pd.read_csv(StringIO(data), dtype={'Int': Int64Dtype},
                         engine=engine)
    assert result is not None
    tm.assert_frame_equal(result, df)

    df = pd.DataFrame({'Int': pd.Series([1, 2, 3], dtype='Int8'),
                       'A': [1, 2, 1]})
    data = df.to_csv(index=False)
    result = pd.read_csv(StringIO(data), dtype={'Int': 'Int8'},
                         engine=engine)
    assert result is not None
    tm.assert_frame_equal(result, df)
