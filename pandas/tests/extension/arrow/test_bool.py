import numpy as np
import pytest

import pandas as pd
from pandas.tests.extension import base
import pandas.util.testing as tm

pytest.importorskip("pyarrow", minversion="0.10.0")

from .bool import ArrowBoolArray, ArrowBoolDtype  # isort:skip


@pytest.fixture
def dtype():
    return ArrowBoolDtype()


@pytest.fixture
def data():
    values = np.random.randint(0, 2, size=100, dtype=bool)
    values[1] = ~values[0]
    return ArrowBoolArray.from_scalars(values)


@pytest.fixture
def data_missing():
    return ArrowBoolArray.from_scalars([None, True])


class BaseArrowTests:
    pass


class TestDtype(BaseArrowTests, base.BaseDtypeTests):
    def test_array_type_with_arg(self, data, dtype):
        pytest.skip("GH-22666")


class TestInterface(BaseArrowTests, base.BaseInterfaceTests):
    def test_copy(self, data):
        # __setitem__ does not work, so we only have a smoke-test
        data.copy()


class TestConstructors(BaseArrowTests, base.BaseConstructorsTests):
    def test_from_dtype(self, data):
        pytest.skip("GH-22666")

    # seems like some bug in isna on empty BoolArray returning floats.
    @pytest.mark.xfail(reason="bad is-na for empty data")
    def test_from_sequence_from_cls(self, data):
        super().test_from_sequence_from_cls(data)


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
