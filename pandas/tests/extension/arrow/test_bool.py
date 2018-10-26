import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.tests.extension import base

from .bool import ArrowBoolArray, ArrowBoolDtype

pytest.importorskip('pyarrow', minversion="0.10.0")



@pytest.fixture
def dtype():
    return ArrowBoolDtype()


@pytest.fixture
def data():
    return ArrowBoolArray.from_scalars(np.random.randint(0, 2, size=100,
                                       dtype=bool))


@pytest.fixture
def data_missing():
    return ArrowBoolArray.from_scalars([None, True])


class BaseArrowTests(object):
    pass


class TestDtype(BaseArrowTests, base.BaseDtypeTests):
    def test_array_type_with_arg(self, data, dtype):
        pytest.skip("GH-22666")


class TestInterface(BaseArrowTests, base.BaseInterfaceTests):
    def test_repr(self, data):
        raise pytest.skip("TODO")


class TestConstructors(BaseArrowTests, base.BaseConstructorsTests):
    def test_from_dtype(self, data):
        pytest.skip("GH-22666")


class TestReduce(base.BaseNoReduceTests):
    def test_reduce_series_boolean(self):
        pass


class TestReduceBoolean(base.BaseBooleanReduceTests):
    pass


def test_is_bool_dtype(data):
    assert pd.api.types.is_bool_dtype(data)
    assert pd.core.common.is_bool_indexer(data)
    s = pd.Series(range(len(data)))
    result = s[data]
    expected = s[np.asarray(data)]
    tm.assert_series_equal(result, expected)
