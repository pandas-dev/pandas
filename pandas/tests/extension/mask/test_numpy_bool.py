"""
This file contains a minimal set of tests for compliance with the extension
array interface test suite, and should contain no other tests.
The test suite for the full functionality of the array is located in
`pandas/tests/arrays/`.

The tests in this file are inherited from the BaseExtensionTests, and only
minimal tweaks should be applied to get the tests passing (by overwriting a
parent method).

Additional tests should either be added to one of the BaseExtensionTests
classes (if they are relevant for the extension interface for all dtypes), or
be added to the array-specific tests in `pandas/tests/arrays/`.

"""
import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays.mask._numpy import NumpyMaskArray, NumpyMaskDtype
from pandas.tests.extension import base
import pandas.util.testing as tm


@pytest.fixture
def dtype():
    return NumpyMaskDtype()


@pytest.fixture
def data():
    return NumpyMaskArray.from_scalars(
        np.random.randint(0, 2, size=100, dtype=bool))


@pytest.fixture
def data_missing():
    pytest.skip("not supported in NumpyMaskArray")


class BaseNumpyTests(object):
    pass


class TestDtype(BaseNumpyTests, base.BaseDtypeTests):
    pass


class TestInterface(BaseNumpyTests, base.BaseInterfaceTests):
    pass


class TestConstructors(BaseNumpyTests, base.BaseConstructorsTests):
    def test_from_dtype(self, data):
        pytest.skip("GH-22666")


class TestReduceBoolean(base.BaseBooleanReduceTests):

    @pytest.mark.parametrize('skipna', [True, False])
    def test_reduce_series(
            self, data, only_numeric_reductions, skipna):
        op_name = only_numeric_reductions
        s = pd.Series(data)
        with pytest.raises(TypeError):
            getattr(s, op_name)(skipna=skipna)

    @pytest.mark.parametrize('skipna', [True, False])
    def test_reduce_series_non_numeric(
            self, data, only_non_numeric_reductions, skipna):
        op_name = only_non_numeric_reductions
        s = pd.Series(data)
        if op_name == 'sum':
            self.check_reduce(s, op_name, skipna)
        else:
            self.check_reduce_bool(s, op_name, skipna)


def test_is_bool_dtype(data):
    assert pd.api.types.is_bool_dtype(data)
    assert pd.core.common.is_bool_indexer(data)
    s = pd.Series(range(len(data)))
    result = s[data]
    expected = s[np.asarray(data)]
    tm.assert_series_equal(result, expected)
