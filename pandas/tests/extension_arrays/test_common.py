import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.common import is_extension_array_dtype
from pandas.core.dtypes.dtypes import ExtensionDtype


class DummyDtype(ExtensionDtype):
    pass


class DummyArray(ExtensionArray):

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype):
        return self.data


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


def test_astype_raises():
    arr = DummyArray(np.array([1, 2, 3]))

    # type  int for py2
    # class int for py3
    xpr = ("DummyArray can only be coerced to 'object' dtype, not "
           "'<.* 'int'>'")

    with tm.assert_raises_regex(ValueError, xpr):
        arr.astype(int)
