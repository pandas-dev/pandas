import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import integer_array
from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
)


@pytest.fixture(
    params=[
        Int8Dtype,
        Int16Dtype,
        Int32Dtype,
        Int64Dtype,
        UInt8Dtype,
        UInt16Dtype,
        UInt32Dtype,
        UInt64Dtype,
    ]
)
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    return integer_array(
        list(range(8)) + [np.nan] + list(range(10, 98)) + [np.nan] + [99, 100],
        dtype=dtype,
    )


def test_uses_pandas_na():
    a = pd.array([1, None], dtype=pd.Int64Dtype())
    assert a[1] is pd.NA


def test_from_dtype_from_float(data):
    # construct from our dtype & string dtype
    dtype = data.dtype

    # from float
    expected = pd.Series(data)
    result = pd.Series(data.to_numpy(na_value=np.nan, dtype="float"), dtype=str(dtype))
    tm.assert_series_equal(result, expected)

    # from int / list
    expected = pd.Series(data)
    result = pd.Series(np.array(data).tolist(), dtype=str(dtype))
    tm.assert_series_equal(result, expected)

    # from int / array
    expected = pd.Series(data).dropna().reset_index(drop=True)
    dropped = np.array(data.dropna()).astype(np.dtype((dtype.type)))
    result = pd.Series(dropped, dtype=str(dtype))
    tm.assert_series_equal(result, expected)
