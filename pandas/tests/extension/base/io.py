import pytest
import pandas as pd
import numpy as np
from pandas.compat import StringIO
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype, UInt8Dtype, UInt16Dtype,
    UInt32Dtype, UInt64Dtype, integer_array,
)
from .base import BaseExtensionTests


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


class ExtensionParsingTests(BaseExtensionTests):
    def test_EA_types(self):
        df = pd.DataFrame({'Int': pd.Series([1, 2, 3], dtype='Int64'),
                           'A': [1, 2, 1]})
        data = df.to_csv(index=False)
        result = pd.read_csv(StringIO(data), dtype={'Int': Int64Dtype})
        assert result is not None

        df = pd.DataFrame({'Int': pd.Series([1, 2, 3], dtype='Int8'),
                           'A': [1, 2, 1]})
        data = df.to_csv(index=False)
        result = pd.read_csv(StringIO(data), dtype={'Int': 'Int8'})
        assert result is not None
