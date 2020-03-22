import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas.core.dtypes.generic import ABCIndexClass

import pandas as pd
import pandas._testing as tm
from pandas.api.types import is_float, is_float_dtype, is_integer, is_scalar
from pandas.core.arrays import IntegerArray, integer_array
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
from pandas.tests.extension.base import BaseOpsUtil


def make_data():
    return list(range(8)) + [np.nan] + list(range(10, 98)) + [np.nan] + [99, 100]


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
    return integer_array(make_data(), dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    return integer_array([np.nan, 1], dtype=dtype)


@pytest.fixture(params=["data", "data_missing"])
def all_data(request, data, data_missing):
    """Parametrized fixture giving 'data' and 'data_missing'"""
    if request.param == "data":
        return data
    elif request.param == "data_missing":
        return data_missing
