import numpy as np
import pytest

from pandas.compat import StringIO

import pandas as pd
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype, UInt8Dtype, UInt16Dtype,
    UInt32Dtype, UInt64Dtype, integer_array)

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

    @pytest.mark.parametrize('engine', ['c', 'python'])
    def test_EA_types(self, engine, data):
        df = pd.DataFrame({'Int': pd.Series(data, dtype=str(data.dtype)),
                           'A': data})
        data = df.to_csv(index=False)
        result = pd.read_csv(StringIO(data), dtype={'Int': str(data.dtype)},
                             engine=engine)
        assert result is not None
